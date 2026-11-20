"""Apertura plana 3D sobre baffle rígido infinito: integral de Rayleigh (Fase 6a).

Modelo
------
La apertura vive en el plano ``z = 0`` y radia hacia ``z > 0``. Sobre un baffle infinito rígido la
presión es la integral de Rayleigh-Sommerfeld con la imagen del baffle (factor 2)::

    p(x) = 2·i·ω·ρ · ∫_S v_n(y) · G(|x − y|) dS ,      G(r) = e^{i·k·r} / (4π·r)

que es el análogo 3D exacto de la integral 2D de `aperture.py` (allí ``G = (i/4)·H₀⁽¹⁾``). Con
``ω/k = c`` el prefactor se puede reescribir como ``(i·ω·ρ/2π) ∫ v_n·e^{ikr}/r dS``, que es la forma
clásica de Rayleigh-Sommerfeld. No hay factor de oblicuidad que agregar: el factor 2 ya contiene la
imagen del baffle. Lo único aproximado es la cuadratura.

Unidades (SI): ``[ω·ρ·v_n] = Pa/m²`` y ``dS`` en m² ⇒ Pa.

Convención (D4): tiempo ``e^{−iωt}``, ondas salientes ``∝ e^{+i·k·r}``, y el fluido en ``z > 0``.

Verificación analítica en el eje
--------------------------------
Para un pistón circular uniforme (radio ``a``, velocidad normal ``v₀``) la integral tiene forma
cerrada en el eje, que es lo que valida el módulo contra la solución de Rayleigh-Sommerfeld::

    p(0,0,z) = ρ·c·v₀·( e^{i·k·R} − e^{i·k·z} ),      R = √(a² + z²)

(se obtiene integrando en anillos: ``∫ ρ'/√(ρ'²+z²)·e^{ik√(ρ'²+z²)} dρ'`` es exacta). Su módulo es
``2·ρ·c·v₀·|sin(k(R−z)/2)|``, con los máximos de interferencia en ``R − z = (2m+1)·λ/2``.
:func:`rayleigh_sommerfeld_axis` implementa esa fórmula.

Foco
----
Una apertura enfocada lleva la ley de fase ``v_n(y) = v₀·exp(−i·k·|F − y|)`` (D13: la ley de foco
es ``−k|F−y|``; el signo ``+`` es para apuntar, no para enfocar). Con eso la integral de Rayleigh
sigue siendo exacta dentro del modelo: el foco es una prescripción de fase sobre la apertura, no una
aproximación.

Cuadratura
----------
- **Círculo**: Gauss-Legendre en ``ρ ∈ [0, a]`` (con el jacobiano ``ρ``) y trapecio uniforme en
  ``θ ∈ [0, 2π)``, que para un integrando periódico es espectralmente exacto.
- **Rectángulo**: Gauss-Legendre tensorial en ``x`` e ``y``.

El integrando tiene un pico de anchura ``z`` cuando el punto de evaluación está cerca de la
apertura: la precisión se degrada para ``z ≲ a/quad_order`` (o ``z`` menor que el espaciado de
nodos). No es un problema para campo lejano, axial o el plano focal, que es el uso previsto.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.polynomial.legendre import leggauss

from .materials import TWO_PI, Material

__all__ = [
    "Aperture3D",
    "DEFAULT_QUAD_ORDER",
    "pressure",
    "axis_pressure",
    "axial_plane",
    "rayleigh_sommerfeld_axis",
]

#: Nodos de Gauss-Legendre por eje de la cuadratura (el error decae espectralmente).
DEFAULT_QUAD_ORDER = 64

#: Puntos por bloque en la evaluación vectorizada (acota la memoria de ``(bloque, nodos)``).
_BLOCK = 256


@dataclass(frozen=True)
class Aperture3D:
    """Apertura plana en ``z = 0``: pistón circular o rectangular.

    Atributos
    ---------
    kind : ``"circle"`` o ``"rectangle"``.
    radius : radio [m] (sólo círculo).
    width, height : lados en ``x`` e ``y`` [m] (sólo rectángulo).
    center : ``(x₀, y₀)`` centro de la apertura [m].
    amplitude : velocidad normal compleja ``v_n`` [m/s] uniforme; ``|amplitude|`` es la amplitud y
        su fase el delay global.
    focus : punto focal ``(Fx, Fy, Fz)`` [m] o ``None``. Si se da, la fase sobre la apertura es
        ``exp(−i·k·|F − y|)`` (foco real, no paraxial).

    Construir con :meth:`circle` o :meth:`rectangle`.
    """

    kind: str
    radius: float | None = None
    width: float | None = None
    height: float | None = None
    center: tuple[float, float] = (0.0, 0.0)
    amplitude: complex = 1.0
    focus: tuple[float, float, float] | None = None

    def __post_init__(self) -> None:
        if self.kind not in ("circle", "rectangle"):
            raise ValueError(f"kind debe ser 'circle' o 'rectangle'; recibido {self.kind!r}")
        if self.kind == "circle":
            if self.radius is None or not (float(self.radius) > 0.0):
                raise ValueError(f"un círculo necesita radius > 0 [m]; recibido {self.radius!r}")
            if self.width is not None or self.height is not None:
                raise ValueError("un círculo no lleva width/height")
        else:
            for name, value in (("width", self.width), ("height", self.height)):
                if value is None or not (float(value) > 0.0):
                    raise ValueError(f"un rectángulo necesita {name} > 0 [m]; recibido {value!r}")
            if self.radius is not None:
                raise ValueError("un rectángulo no lleva radius")
        c = tuple(float(v) for v in self.center)
        if len(c) != 2 or not np.all(np.isfinite(c)):
            raise ValueError(f"center debe ser (x₀, y₀) finito; recibido {self.center!r}")
        if not np.isfinite(np.abs(complex(self.amplitude))):
            raise ValueError(f"amplitude debe ser finita; recibido {self.amplitude!r}")
        if self.focus is not None:
            f = tuple(float(v) for v in self.focus)
            if len(f) != 3 or not np.all(np.isfinite(f)):
                raise ValueError(f"focus debe ser (Fx, Fy, Fz) finito; recibido {self.focus!r}")
            if not f[2] > 0.0:
                raise ValueError(
                    f"focus debe estar delante del baffle (Fz > 0); recibido Fz = {f[2]:.6g} m"
                )
            object.__setattr__(self, "focus", f)
        object.__setattr__(self, "center", c)

    @classmethod
    def circle(cls, radius: float, center=(0.0, 0.0), amplitude: complex = 1.0,
               focus=None) -> "Aperture3D":
        """Pistón circular de radio ``radius`` [m]."""
        return cls(kind="circle", radius=float(radius), center=center, amplitude=amplitude,
                   focus=focus)

    @classmethod
    def rectangle(cls, width: float, height: float, center=(0.0, 0.0), amplitude: complex = 1.0,
                  focus=None) -> "Aperture3D":
        """Pistón rectangular de lados ``width`` × ``height`` [m]."""
        return cls(kind="rectangle", width=float(width), height=float(height), center=center,
                   amplitude=amplitude, focus=focus)

    @property
    def area(self) -> float:
        """Área de la apertura [m²]."""
        radius = 0.0 if self.radius is None else float(self.radius)
        width = 0.0 if self.width is None else float(self.width)
        height = 0.0 if self.height is None else float(self.height)
        if self.kind == "circle":
            return float(np.pi * radius**2)
        return float(width * height)

    def bounds(self) -> tuple[float, float, float, float]:
        """``(x_min, x_max, y_min, y_max)`` de la apertura [m]."""
        x0, y0 = self.center
        if self.kind == "circle":
            r = 0.0 if self.radius is None else float(self.radius)
            return (x0 - r, x0 + r, y0 - r, y0 + r)
        hw = 0.5 * (0.0 if self.width is None else float(self.width))
        hh = 0.5 * (0.0 if self.height is None else float(self.height))
        return (x0 - hw, x0 + hw, y0 - hh, y0 + hh)


# --------------------------------------------------------------------------- #
# Cuadratura
# --------------------------------------------------------------------------- #


def _gauss_01(order: int) -> tuple[np.ndarray, np.ndarray]:
    """Nodos y pesos de Gauss-Legendre en ``[0, 1]`` (``Σw = 1``)."""
    if int(order) < 2:
        raise ValueError(f"quad_order debe ser >= 2; recibido {order!r}")
    x, w = leggauss(int(order))
    return 0.5 * (x + 1.0), 0.5 * w


def _aperture_quadrature(aperture: Aperture3D, order: int) -> tuple[np.ndarray, np.ndarray]:
    """Nodos ``(Q, 3)`` [m] y pesos ``(Q,)`` [m²] de la cuadratura sobre la apertura."""
    u, wu = _gauss_01(order)
    x0, y0 = aperture.center
    if aperture.kind == "circle":
        a = 0.0 if aperture.radius is None else float(aperture.radius)
        rho = a * u
        w_rho = a * wu * rho  # jacobiano ρ dρ
        theta = TWO_PI * np.arange(order) / order
        w_theta = np.full(order, TWO_PI / order)  # trapecio periódico: Σ = 2π = ∫dθ
        rr, tt = np.meshgrid(rho, theta, indexing="ij")
        ww = np.outer(w_rho, w_theta)
        nodes = np.column_stack([
            (x0 + rr * np.cos(tt)).ravel(),
            (y0 + rr * np.sin(tt)).ravel(),
            np.zeros(rr.size),
        ])
        return nodes, ww.ravel()
    hw = 0.5 * (0.0 if aperture.width is None else float(aperture.width))
    hh = 0.5 * (0.0 if aperture.height is None else float(aperture.height))
    xs = x0 - hw + 2.0 * hw * u
    ys = y0 - hh + 2.0 * hh * u
    gx, gy = np.meshgrid(xs, ys, indexing="ij")
    ww = np.outer(2.0 * hw * wu, 2.0 * hh * wu)
    nodes = np.column_stack([gx.ravel(), gy.ravel(), np.zeros(gx.size)])
    return nodes, ww.ravel()


def _as_points(points) -> np.ndarray:
    """Valida y normaliza ``points`` a ``(M, 3)`` [m]."""
    pts = np.asarray(points, dtype=float)
    if pts.ndim == 1:
        if pts.shape != (3,):
            raise ValueError(f"un punto suelto debe tener forma (3,); recibido {pts.shape}")
        pts = pts[None, :]
    if pts.ndim != 2 or pts.shape[1] != 3:
        raise ValueError(f"points debe tener forma (M, 3); recibido {pts.shape}")
    if pts.shape[0] == 0:
        raise ValueError("points está vacío")
    if not np.all(np.isfinite(pts)):
        raise ValueError("points contiene valores no finitos (NaN o inf)")
    if np.any(pts[:, 2] <= 0.0):
        bad = int(np.count_nonzero(pts[:, 2] <= 0.0))
        raise ValueError(
            "points debe tener z > 0: el fluido está delante del baffle (z = 0 es el plano de la "
            f"apertura, donde 1/r es singular); {bad} punto(s) con z <= 0, z mínimo = "
            f"{pts[:, 2].min():.6g} m"
        )
    return pts


# --------------------------------------------------------------------------- #
# API
# --------------------------------------------------------------------------- #


def pressure(points, frequency: float, material: Material, aperture: Aperture3D,
             quad_order: int = DEFAULT_QUAD_ORDER) -> np.ndarray:
    """Presión compleja [Pa] radiada por la apertura sobre baffle rígido.

    ``p(x) = 2·i·ω·ρ · ∫_S v_n(y)·e^{i·k·|x−y|}/(4π·|x−y|) dS``, con la fase de foco aplicada a
    ``v_n`` si ``aperture.focus`` no es ``None``.

    Parameters
    ----------
    points : ``(M, 3)`` posiciones de evaluación [m], todas con ``z > 0``.
    frequency : [Hz].
    material : :class:`materials.Material`; de ahí salen ``c``, ``ρ`` y ``α(f)``.
    aperture : :class:`Aperture3D`.
    quad_order : nodos de Gauss por eje de la cuadratura.

    Returns
    -------
    ``(M,)`` complejo [Pa]. Es lineal en ``aperture.amplitude``.
    """
    if not isinstance(material, Material):
        raise TypeError("material debe ser un materials.Material")
    if not isinstance(aperture, Aperture3D):
        raise TypeError(f"aperture debe ser Aperture3D; recibido {type(aperture).__name__}")
    f = float(frequency)
    if not np.isfinite(f) or f <= 0.0:
        raise ValueError(f"frequency debe ser finita y > 0 [Hz]; recibido {frequency!r}")
    if int(quad_order) != quad_order or int(quad_order) < 2:
        raise ValueError(f"quad_order debe ser un entero >= 2; recibido {quad_order!r}")

    pts = _as_points(points)
    nodes, weights = _aperture_quadrature(aperture, int(quad_order))

    k = material.wavenumber(f)
    omega = TWO_PI * f
    prefactor = 2.0j * omega * material.rho

    # velocidad normal compleja en cada nodo de la apertura (fase de foco incluida)
    v_n = np.full(nodes.shape[0], complex(aperture.amplitude), dtype=complex)
    if aperture.focus is not None:
        focus = np.asarray(aperture.focus, dtype=float)
        v_n = v_n * np.exp(-1j * k * np.linalg.norm(nodes - focus[None, :], axis=1))
    source = weights * v_n

    out = np.empty(pts.shape[0], dtype=complex)
    for start in range(0, pts.shape[0], _BLOCK):
        block = pts[start:start + _BLOCK]
        r = np.linalg.norm(block[:, None, :] - nodes[None, :, :], axis=2)
        green = np.exp(1j * k * r) / (4.0 * np.pi * r)
        out[start:start + block.shape[0]] = green @ source
    return prefactor * out


def axis_pressure(z, frequency: float, material: Material, aperture: Aperture3D,
                  quad_order: int = DEFAULT_QUAD_ORDER) -> np.ndarray:
    """Presión en el eje de la apertura, ``(x₀, y₀, z)``, para cada ``z`` [m] de ``z``.

    Devuelve un array de la misma forma que ``z`` (escalar incluido).
    """
    zs = np.asarray(z, dtype=float)
    flat = np.atleast_1d(zs).reshape(-1)
    if np.any(~np.isfinite(flat)) or np.any(flat <= 0.0):
        raise ValueError("z debe ser finito y > 0 [m] (el fluido está en z > 0)")
    x0, y0 = aperture.center
    pts = np.column_stack([np.full(flat.size, x0), np.full(flat.size, y0), flat])
    values = pressure(pts, frequency, material, aperture, quad_order=quad_order)
    if np.ndim(z) == 0:
        return np.asarray(values[0])
    return values.reshape(zs.shape)


def axial_plane(frequency: float, material: Material, aperture: Aperture3D,
                x_extent: tuple[float, float], z_extent: tuple[float, float],
                nx: int, nz: int, quad_order: int = DEFAULT_QUAD_ORDER
                ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Campo en el plano axial ``y = y₀`` (el plano que contiene el eje de la apertura).

    Parameters
    ----------
    x_extent : ``(x_min, x_max)`` [m].
    z_extent : ``(z_min, z_max)`` [m], con ``z_min > 0``.
    nx, nz : número de puntos en ``x`` y en ``z``.

    Returns
    -------
    ``(X, Z, P)`` : mallas ``(nz, nx)`` [m] y campo complejo ``(nz, nx)`` [Pa].
    """
    x0, x1 = (float(v) for v in x_extent)
    z0, z1 = (float(v) for v in z_extent)
    if not (x1 > x0):
        raise ValueError(f"x_extent debe ser creciente; recibido {x_extent!r}")
    if not (z1 > z0 > 0.0):
        raise ValueError(f"z_extent debe cumplir z1 > z0 > 0 [m]; recibido {z_extent!r}")
    if int(nx) < 2 or int(nz) < 2:
        raise ValueError(f"nx y nz deben ser >= 2; recibido {nx!r}, {nz!r}")

    xs = np.linspace(x0, x1, int(nx))
    zs = np.linspace(z0, z1, int(nz))
    gx, gz = np.meshgrid(xs, zs, indexing="xy")  # (nz, nx)
    pts = np.column_stack([gx.ravel(), np.full(gx.size, aperture.center[1]), gz.ravel()])
    values = pressure(pts, frequency, material, aperture, quad_order=quad_order)
    return gx, gz, values.reshape(gx.shape)


def rayleigh_sommerfeld_axis(z, frequency: float, material: Material, radius: float,
                             amplitude: complex = 1.0) -> np.ndarray:
    """Solución exacta en el eje de un pistón circular uniforme: ``ρ·c·v₀·(e^{ikR} − e^{ikz})``.

    Es la referencia analítica del módulo (ver el docstring de arriba). ``R = √(a² + z²)``,
    ``a = radius``, ``v₀ = amplitude``.

    Parameters
    ----------
    z : distancia al plano de la apertura [m], ``z > 0``.
    radius : radio del pistón [m].
    amplitude : velocidad normal compleja ``v₀`` [m/s].

    Returns
    -------
    Presión compleja [Pa], con la forma de ``z``.
    """
    if not isinstance(material, Material):
        raise TypeError("material debe ser un materials.Material")
    f = float(frequency)
    if not np.isfinite(f) or f <= 0.0:
        raise ValueError(f"frequency debe ser finita y > 0 [Hz]; recibido {frequency!r}")
    a = float(radius)
    if not np.isfinite(a) or a <= 0.0:
        raise ValueError(f"radius debe ser finito y > 0 [m]; recibido {radius!r}")
    zs = np.asarray(z, dtype=float)
    flat = np.atleast_1d(zs).reshape(-1)
    if np.any(~np.isfinite(flat)) or np.any(flat <= 0.0):
        raise ValueError("z debe ser finito y > 0 [m]")

    k = material.wavenumber(f)
    big_r = np.sqrt(a**2 + flat**2)
    values = material.rho * material.c * complex(amplitude) * (
        np.exp(1j * k * big_r) - np.exp(1j * k * flat)
    )
    if np.ndim(z) == 0:
        return np.asarray(values[0])
    return values.reshape(zs.shape)
