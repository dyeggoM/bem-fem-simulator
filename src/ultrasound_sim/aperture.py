"""Fuente por apertura plana sobre un baffle rígido infinito (Fase 4B).

Modelo (``docs/plan.md``, Fase 4B)
----------------------------------
La apertura vive en el plano ``x = 0``, se extiende a lo largo de ``y`` y radia hacia ``x > 0``.
Sobre el baffle infinito rígido la presión en el fluido es la integral de Rayleigh 2D con la imagen
del baffle (factor 2):

    p(x) = 2·i·ω·ρ · ∫_apertura v_n(y) · G(|x − y|) dy ,      G(r) = (i/4)·H₀⁽¹⁾(k·r)

con ``ω = 2π·f`` y ``k`` de :meth:`materials.Material.wavenumber` (complejo si hay atenuación, y se
evalúa tal cual: ``H₀⁽¹⁾`` acepta argumento complejo). Unidades: ``[ω·ρ·v_n] = Pa/m`` y ``dy`` en
metros ⇒ Pa.

Nada de esto es una aproximación dentro del modelo: no hay factor de oblicuidad añadido, porque el
factor 2 ya contiene la imagen del baffle (es la integral de Rayleigh-Sommerfeld exacta para baffle
plano rígido). Lo único aproximado es la cuadratura.

Convención (``docs/decisions.md`` D4): tiempo ``exp(−i·ω·t)``, ondas salientes ``∝ exp(+i·k·r)``,
y el fluido está en ``x > 0``. Para una apertura puntual el prefactor compuesto es
``2i·ω·ρ·(i/4) = −ω·ρ/2``, o sea ``p = −(ω·ρ/2)·v_n·H₀⁽¹⁾(k·r)``; con ``v_n`` real la presión sobre
el eje es la de una fuente 2D saliente con ese signo.

Cuadratura
----------
Gauss-Legendre por apertura. Cerca del campo cercano el integrando tiene una casi-singularidad
logarítmica (``H₀`` diverge como ``ln r``), así que:

* si la distancia del punto al segmento es menor que ``refinement`` veces el subsegmento base
  (``width/quad_order``), la apertura se subdivide por bisección;
* la bisección se detiene cuando el subsegmento mide menos que ``distancia/refinement`` (es decir,
  el integrando ya es suave sobre él), con un tope de profundidad.

El test de convergencia mide contra una cuadratura 10× más fina (orden ×10 y ``refinement`` ×10).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import hankel1

from .materials import TWO_PI, Material

__all__ = [
    "Aperture",
    "DEFAULT_QUAD_ORDER",
    "DEFAULT_REFINEMENT",
    "far_field_pattern",
    "pressure",
]

#: Orden de Gauss-Legendre por subsegmento (16 nodos integran un subsegmento suave al redondeo).
DEFAULT_QUAD_ORDER = 16

#: Criterio de cercanía / graduación, en unidades del subsegmento base ``width/quad_order``.
DEFAULT_REFINEMENT = 2.0

#: Tope de bisección (2^40 subdivisiones relativas; nunca se alcanza con ``x > 0`` razonable).
_MAX_DEPTH = 40

#: Puntos por bloque en la evaluación vectorizada (acota la memoria de ``(m, order)``).
_BLOCK = 256


@dataclass(frozen=True)
class Aperture:
    """Una apertura (un transductor) sobre el plano ``x = 0``, centrada en ``y = center_y``.

    Atributos
    ---------
    center_y : centro de la apertura sobre el eje ``y`` [m].
    width    : ancho ``D`` [m] (la apertura cubre ``center_y ± width/2``).
    amplitude: velocidad normal compleja ``v_n`` [m/s]; su fase es el delay del elemento.
               ``|amplitude|`` es la velocidad normal prescrita (uniforme sobre la apertura).
    """

    center_y: float
    width: float
    amplitude: complex = 1.0

    def __post_init__(self) -> None:
        if not (float(self.width) > 0.0):
            raise ValueError(f"Aperture: width debe ser > 0; recibido {self.width!r}")
        if not np.isfinite(float(self.center_y)):
            raise ValueError(f"Aperture: center_y debe ser finito; recibido {self.center_y!r}")

    @property
    def edges(self) -> tuple[float, float]:
        """``(y_min, y_max)`` de la apertura [m]."""
        half = 0.5 * float(self.width)
        return (float(self.center_y) - half, float(self.center_y) + half)

    def __repr__(self) -> str:  # pragma: no cover - cosmético
        y0, y1 = self.edges
        return (f"Aperture(y=[{y0 * 1e3:.3f}, {y1 * 1e3:.3f}] mm, "
                f"width={self.width * 1e3:.3f} mm, amplitude={self.amplitude:.3g})")


# --------------------------------------------------------------------------- #
# cuadratura
# --------------------------------------------------------------------------- #


def _unit_gauss(order: int) -> tuple[np.ndarray, np.ndarray]:
    """Nodos y pesos de Gauss-Legendre en ``[0, 1]`` (``Σw = 1``)."""
    if order < 2:
        raise ValueError(f"quad_order debe ser >= 2; recibido {order!r}")
    x, w = leggauss(int(order))
    return 0.5 * (x + 1.0), 0.5 * w


def _graded_nodes(y0: float, y1: float, px: float, py: float,
                  unit_nodes: np.ndarray, unit_weights: np.ndarray,
                  refinement: float) -> tuple[np.ndarray, np.ndarray]:
    """Nodos/pesos sobre ``[y0, y1]`` graduados hacia el punto ``(px, py)``.

    Bisección recursiva: un subsegmento se parte mientras su longitud exceda
    ``distancia(punto, subsegmento)/refinement``. Los subsegmentos lejanos no se tocan (su
    integrando ya es suave), así que el costo crece sólo logarítmicamente al acercarse.
    """
    nodes: list[np.ndarray] = []
    weights: list[np.ndarray] = []
    stack = [(float(y0), float(y1), 0)]
    while stack:
        a, b, depth = stack.pop()
        length = b - a
        closest = min(max(py, a), b)
        dist = math.hypot(px, py - closest)
        if depth >= _MAX_DEPTH or length <= dist / refinement:
            y = a + length * unit_nodes
            nodes.append(y)
            weights.append(length * unit_weights)
        else:
            mid = 0.5 * (a + b)
            stack.append((a, mid, depth + 1))
            stack.append((mid, b, depth + 1))
    return np.concatenate(nodes), np.concatenate(weights)


def _distance_to_segment(px: np.ndarray, py: np.ndarray, y0: float, y1: float) -> np.ndarray:
    """Distancia de cada punto ``(px, py)`` al segmento ``[y0, y1]`` sobre ``x = 0``."""
    closest = np.clip(py, y0, y1)
    return np.hypot(px, py - closest)


def _aperture_integral(pts: np.ndarray, aperture: Aperture, k: complex,
                       quad_order: int, refinement: float) -> np.ndarray:
    """``∫_apertura G(|x − y|) dy`` en cada punto de ``pts`` ``(M, 2)`` ⇒ ``(M,)`` complejo."""
    y0, y1 = aperture.edges
    width = y1 - y0
    unit_nodes, unit_weights = _unit_gauss(quad_order)
    nodes = y0 + width * unit_nodes
    weights = width * unit_weights

    out = np.empty(pts.shape[0], dtype=complex)
    h = width / quad_order                      # subsegmento base
    for start in range(0, pts.shape[0], _BLOCK):
        block = pts[start:start + _BLOCK]
        r = np.hypot(block[:, :1], nodes[None, :] - block[:, 1:2])
        out[start:start + block.shape[0]] = (0.25j * hankel1(0, k * r)) @ weights

        # pares cercanos: la cuadratura fija pierde dígitos, se gradúa hacia el punto
        d = _distance_to_segment(block[:, 0], block[:, 1], y0, y1)
        for i in np.nonzero(d < refinement * h)[0]:
            gn, gw = _graded_nodes(y0, y1, block[i, 0], block[i, 1],
                                   unit_nodes, unit_weights, refinement)
            r = np.hypot(block[i, 0], gn - block[i, 1])
            out[start + i] = np.sum(0.25j * hankel1(0, k * r) * gw)
    return out


# --------------------------------------------------------------------------- #
# API
# --------------------------------------------------------------------------- #


def _as_points(points) -> np.ndarray:
    pts = np.asarray(points, dtype=float)
    if pts.ndim == 1:
        if pts.shape != (2,):
            raise ValueError(f"un punto suelto debe tener forma (2,); recibido {pts.shape}")
        pts = pts[None, :]
    if pts.ndim != 2 or pts.shape[1] != 2:
        raise ValueError(f"points debe tener forma (M, 2); recibido {pts.shape}")
    if pts.shape[0] == 0:
        raise ValueError("points está vacío")
    if not np.all(np.isfinite(pts)):
        raise ValueError("points contiene valores no finitos")
    return pts


def pressure(points, frequency: float, material: Material, apertures,
             quad_order: int = DEFAULT_QUAD_ORDER,
             refinement: float = DEFAULT_REFINEMENT) -> np.ndarray:
    """Presión en ``points`` [Pa] radiada por ``apertures`` sobre baffle rígido.

    ``p(x) = 2·i·ω·ρ · Σ_j A_j ∫_{apertura j} G(|x − y|) dy``,
    ``G = (i/4)H₀⁽¹⁾(k·r)``, ``A_j`` la amplitud compleja (velocidad normal) de cada apertura.

    Parámetros
    ----------
    points : ``(M, 2)`` [m]; todos deben cumplir ``x > 0`` (el modelo no está definido sobre el
        baffle ni detrás). Se acepta ``(2,)`` para un punto suelto.
    frequency : [Hz].
    material : :class:`materials.Material`; de ahí salen ``c``, ``ρ`` y ``α(f)``.
    apertures : :class:`Aperture` o secuencia de ellas.
    quad_order : nodos de Gauss-Legendre por subsegmento base.
    refinement : criterio de cercanía/graduación (ver docstring del módulo).

    Devuelve
    --------
    ``(M,)`` complejo [Pa]. La presión es lineal en las amplitudes: duplicar ``amplitude`` duplica
    ``p`` exactamente (en punto flotante, hasta redondeo).
    """
    if not isinstance(material, Material):
        raise TypeError("material debe ser un materials.Material")
    if not (float(frequency) > 0.0):
        raise ValueError(f"frequency debe ser > 0; recibido {frequency!r}")
    if not (float(refinement) > 0.0):
        raise ValueError(f"refinement debe ser > 0; recibido {refinement!r}")

    pts = _as_points(points)
    if np.any(pts[:, 0] <= 0.0):
        bad = int(np.count_nonzero(pts[:, 0] <= 0.0))
        raise ValueError(
            "points debe tener x > 0 (el fluido está delante del baffle, x = 0 es la apertura); "
            f"{bad} punto(s) con x <= 0, x mínimo = {pts[:, 0].min():.6g} m"
        )

    if isinstance(apertures, Aperture):
        apertures = [apertures]
    aps = list(apertures)
    if not aps:
        raise ValueError("apertures está vacío")
    for a in aps:
        if not isinstance(a, Aperture):
            raise TypeError(f"apertures debe contener Aperture; recibido {type(a).__name__}")

    k = material.wavenumber(frequency)
    omega = TWO_PI * float(frequency)
    prefactor = 2.0j * omega * material.rho

    total = np.zeros(pts.shape[0], dtype=complex)
    for a in aps:
        total += complex(a.amplitude) * _aperture_integral(pts, a, k, quad_order, refinement)
    return prefactor * total


def far_field_pattern(angles, frequency: float, material: Material, apertures,
                      radius: float = 1.0,
                      quad_order: int = DEFAULT_QUAD_ORDER,
                      refinement: float = DEFAULT_REFINEMENT) -> np.ndarray:
    """Campo complejo a distancia ``radius`` en las direcciones ``angles`` (patrón de radiación).

    ``angles`` [rad] se mide desde el eje ``x`` (``0`` = broadside, ``+`` hacia ``+y``); se admiten
    sólo ``|θ| < π/2`` porque los puntos deben caer en ``x > 0``. Los puntos de evaluación son
    ``radius · (cos θ, sin θ)``.

    Para ``r`` grande el campo tiende a ``p ≈ C·(e^{ikr}/√r)·D(θ)`` (2D: decae como ``1/√r``), así
    que ``|D(θ)|`` se recupera normalizando ``|p|`` por su máximo. No se aplica ninguna corrección
    asintótica: se evalúa la integral exacta, de modo que el usuario elige cuán lejos muestrea.

    Devuelve
    --------
    ``(M,)`` complejo [Pa]. Para comparar con el factor de arreglo conviene
    ``np.abs(p) / np.abs(p).max()``.
    """
    if not (float(radius) > 0.0):
        raise ValueError(f"radius debe ser > 0; recibido {radius!r}")
    theta = np.atleast_1d(np.asarray(angles, dtype=float))
    if np.any(np.abs(theta) >= 0.5 * np.pi):
        raise ValueError(
            "angles debe cumplir |θ| < π/2: con el baffle en x = 0, más allá el punto no está en "
            f"el fluido (x = radius·cosθ > 0). Máximo |θ| recibido = {np.abs(theta).max():.6g} rad"
        )
    pts = float(radius) * np.column_stack([np.cos(theta), np.sin(theta)])
    return pressure(pts, frequency, material, apertures, quad_order=quad_order,
                    refinement=refinement)
