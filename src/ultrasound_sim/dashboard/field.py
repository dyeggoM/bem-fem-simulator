"""Motor del dashboard: campo de un conjunto arbitrario de fuentes sobre una malla rectangular.

La idea de diseño es que *una fuente sea un dato, no un parámetro*: cada elemento se describe con
``{"y": posición [m], "width": ancho [m], "amp": amplitud, "phase": fase [rad]}``. Mover una fuente es
cambiar ``y``; agregar una es agregar un diccionario. El arreglo enfocado no es más que una forma de
generar esa lista (``array_sources``), así que la UI puede empezar con un arreglo y después editar
las fuentes una por una sin que el motor cambie.

El cálculo es el mismo de :mod:`ultrasound_sim.aperture` (integral de Rayleigh sobre baffle rígido,
``p = 2·i·ω·ρ·Σ_j A_j ∫ G dy``, ``G = (i/4)H₀⁽¹⁾``), pero vectorizado sobre la malla: los nodos de
Gauss de todas las fuentes se concatenan y se evaluán de una sola vez por fila. Como el integrando es
suave en la malla (no hay puntos singulares sobre la apertura), Gauss de 8 nodos da error de máquina;
verificado contra ``aperture.pressure`` en ``tests/test_dashboard.py``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import hankel1

from .. import focus as focus_mod
from ..materials import Material

#: Nodos de Gauss por fuente. El integrando es analítico en la malla (x > 0), así que 8 sobran:
#: contra ``aperture.pressure`` el error medido es 1e-15.
DEFAULT_QUAD = 8


@dataclass(frozen=True)
class Source:
    """Una apertura plana sobre el baffle ``x = 0``."""

    y: float          # centro [m]
    width: float      # ancho [m]
    amp: float = 1.0
    phase: float = 0.0

    @property
    def amplitude(self) -> complex:
        return self.amp * complex(math.cos(self.phase), math.sin(self.phase))

    @classmethod
    def from_dict(cls, d: dict) -> "Source":
        return cls(
            y=float(d["y"]),
            width=float(d.get("width", 1e-3)),
            amp=float(d.get("amp", 1.0)),
            phase=float(d.get("phase", 0.0)),
        )


def array_sources(
    n_elements: int,
    aperture: float,
    focus_x: float,
    focus_y: float,
    frequency: float,
    material: Material,
) -> list[Source]:
    """Fuentes de un arreglo lineal uniforme enfocado en ``(focus_x, focus_y)``.

    ``pitch = aperture / n_elements`` y ``φ_j = −k_real·|F − y_j|`` (la misma ley que
    :meth:`arrays.LinearArray.focused`, para que el dashboard y el resto del paquete no discrepen).
    """
    if n_elements < 1:
        raise ValueError(f"n_elements debe ser >= 1; recibido {n_elements!r}")
    if aperture <= 0.0:
        raise ValueError(f"aperture debe ser > 0; recibido {aperture!r}")
    pitch = aperture / n_elements
    k = material.wavenumber(frequency).real
    fuentes = []
    for j in range(n_elements):
        y = -aperture / 2.0 + pitch * (j + 0.5)
        fase = -k * math.hypot(focus_x, focus_y - y)
        fuentes.append(Source(y=y, width=pitch, amp=1.0, phase=fase))
    return fuentes


def evaluate(
    xs: np.ndarray,
    ys: np.ndarray,
    sources: list[Source],
    frequency: float,
    material: Material,
    quad: int = DEFAULT_QUAD,
) -> np.ndarray:
    """Presión compleja ``(len(ys), len(xs))`` en la malla ``xs × ys``."""
    if not sources:
        raise ValueError("sources está vacío")
    if np.any(np.asarray(xs) <= 0.0):
        raise ValueError("la malla debe tener x > 0 (el fluido está delante del baffle)")

    k = complex(material.wavenumber(frequency))
    omega = 2.0 * math.pi * frequency

    xi, wi = leggauss(quad)
    ys_nodos, pesos, amps = [], [], []
    for s in sources:
        if s.width <= 0.0:
            raise ValueError(f"source con width <= 0: {s!r}")
        ys_nodos.append(s.y + 0.5 * s.width * xi)
        pesos.append(0.5 * s.width * wi)
        amps.append(np.full(quad, s.amplitude, dtype=complex))
    yq = np.concatenate(ys_nodos)
    wq = np.concatenate(pesos)
    aq = np.concatenate(amps)

    xs = np.asarray(xs, dtype=float)
    out = np.empty((len(ys), xs.size), dtype=complex)
    for i, y0 in enumerate(np.asarray(ys, dtype=float)):
        r = np.hypot(xs[:, None], y0 - yq[None, :])
        out[i] = (aq * wq * (0.25j * hankel1(0, k * r))).sum(axis=1)
    return out * (2j * omega * material.rho)


def metrics(
    xs: np.ndarray,
    ys: np.ndarray,
    p: np.ndarray,
    *,
    pitch: float | None = None,
    frequency: float | None = None,
    material: Material | None = None,
) -> dict:
    """Métricas del campo: pico, FWHM lateral y axial, lóbulo lateral y ángulo de grataje.

    Convenciones (importan, porque equivocarlas da números plausibles y falsos):

    - ``p`` tiene forma ``(len(ys), len(xs))``: **filas = y, columnas = x**.
    - *axial* = a lo largo del haz, ``x`` variable a ``y`` fijo en el pico → ``intensidad[iy, :]``.
    - *lateral* = a través del haz, ``y`` variable a ``x`` fijo en el pico → ``intensidad[:, ix]``.
    - El FWHM se devuelve en **mm** (``lateral_fwhm`` trabaja en metros).
    - Si el lóbulo no cruza −6 dB dentro de la ventana, se devuelve ``None``: es preferible a un
      número que en realidad mide el recorte de la malla.
    - El lóbulo lateral es el máximo **después del primer nulo**, no el flanco del lóbulo principal.
    """
    intensidad = np.abs(p) ** 2
    tope = float(intensidad.max())
    if tope <= 0.0:
        raise ValueError("campo idénticamente nulo")
    iy, ix = (int(v) for v in np.unravel_index(int(intensidad.argmax()), intensidad.shape))

    def _fwhm(coord: np.ndarray, linea: np.ndarray) -> float | None:
        objetivo = 0.5 * float(linea.max())
        if linea[0] > objetivo or linea[-1] > objetivo:
            return None  # el lóbulo no se cierra dentro de la ventana: FWHM recortado
        try:
            return float(focus_mod.lateral_fwhm(coord, linea)) * 1e3
        except (ValueError, IndexError):
            return None

    axial = intensidad[iy, :]       # varía x
    lateral = intensidad[:, ix]     # varía y

    def _primer_nulo(v: np.ndarray) -> int | None:
        """Índice del primer mínimo local recorriendo desde el pico (``v[0]`` = pico)."""
        subidas = np.nonzero(np.diff(v) > 0.0)[0]
        return int(subidas[0]) if subidas.size else None

    niveles = []
    for rama in (lateral[iy:], lateral[: iy + 1][::-1]):
        k = _primer_nulo(rama)
        if k is not None and k + 1 < rama.size:
            niveles.append(float(rama[k + 1:].max()))
    fuera = max(niveles) if niveles else 0.0

    grin = None
    if pitch and frequency and material is not None:
        lam = material.wavelength(frequency)
        cociente = lam / pitch
        if abs(cociente) < 1.0:
            grin = float(math.degrees(math.asin(cociente)))

    return {
        "peak_x_mm": float(xs[ix]) * 1e3,
        "peak_y_mm": float(ys[iy]) * 1e3,
        "p_max": float(np.abs(p).max()),
        "p2_max": tope,
        "fwhm_lateral_mm": _fwhm(np.asarray(ys, float), lateral),
        "fwhm_axial_mm": _fwhm(np.asarray(xs, float), axial),
        "sidelobe_db": float(10.0 * math.log10(fuera / tope)) if fuera > 0.0 else None,
        "grating_deg": grin,
        "db_min": float(10.0 * math.log10(max(tope / 1e6, 1e-30))),
    }
