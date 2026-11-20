"""Arreglos lineales de aperturas con delays (Fase 4B).

Un **transductor** es una :class:`aperture.Aperture`; un **arreglo** son ``N`` aperturas contiguas
(ancho de elemento ``w = pitch``, sin huecos entre elementos) sobre el mismo plano ``x = 0``, cada
una con su amplitud compleja = velocidad normal con delay.

Fases (convención ``exp(−i·ω·t)``, para que todo llegue **en fase** al objetivo)
-------------------------------------------------------------------------------
Con ``k_real = Re(k)`` (la atenuación no cambia la fase, sólo la amplitud):

* enfocar en ``F``: ``φ_j = −k_real · |F − y_j|``
* apuntar a la dirección ``d = (cosθ, sinθ)`` (θ desde ``+x``): ``φ_j = +k_real · y_j · sinθ``

La amplitud del elemento es ``|a|·exp(i·φ_j)``, así que en el blanco todos los términos suman con la
misma fase. El frente de onda queda cuantizado por elemento (fase constante sobre cada pistón): eso
es exactamente lo que distingue un arreglo de una apertura continua, y es el origen del lóbulo de
grataje cuando ``pitch`` crece (ver ``tests/test_array_focus.py``).

**Signo del apuntamiento** (discrepancia con ``docs/plan.md``, resuelta por física y verificada
contra ``array_factor``): el plan escribe ``φ_j = −k_real·(y_j·d)``. Eso apunta al ángulo
*espejado* (``−θ``), porque en este modelo un punto a distancia ``r`` y ángulo ``θ`` recibe
``exp(+i·k·r)`` con ``r ≈ R − y_j·sinθ``, o sea fase ``−k·y_j·sinθ``: para compensarla hace falta
``+k·y_j·sinθ``. La fórmula del plan con el signo menos es consistente sólo si ``d`` se toma como la
dirección *blanco → apertura*. La misma cuenta sale del propio criterio de foco al alejarse:
``−k|R·d − y_j| ≈ −kR + k·y_j·sinθ``. Se usa ``+``, que es lo que exige "llegar en fase al
objetivo" y lo que hace coincidir el máximo medido con el pico de ``array_factor(steer_angle=θ)``
(``tests/test_array_focus.py``).

Factor de arreglo
-----------------
Para ``N`` fuentes puntuales equiespaciadas ``d``, con apuntamiento ``θ₀``:

    ψ = k·d·(sinθ − sinθ₀),      AF(θ) = |sin(Nψ/2) / (N·sin(ψ/2))|

normalizado a 1 en el máximo (``ψ → 0``). ``array_factor`` lo evalúa de forma estable numéricamente
(``sinc(N ψ/2)/sinc(ψ/2)``, sin cancelación catastrófica en ``ψ ≈ 0``).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .aperture import Aperture
from .materials import Material

__all__ = ["LinearArray", "array_factor", "grating_lobe_sine"]


@dataclass(frozen=True)
class LinearArray:
    """Arreglo lineal de ``n_elements`` aperturas contiguas de ancho ``pitch`` sobre ``x = 0``.

    Atributos
    ---------
    center_y  : centro geométrico del arreglo sobre ``y`` [m].
    n_elements: cantidad de elementos ``N``.
    pitch     : paso entre centros de elemento y ancho de cada elemento ``w = pitch`` [m]
                (contiguos, sin kerf).
    """

    center_y: float
    n_elements: int
    pitch: float

    def __post_init__(self) -> None:
        if int(self.n_elements) < 1:
            raise ValueError(f"n_elements debe ser >= 1; recibido {self.n_elements!r}")
        if not (float(self.pitch) > 0.0):
            raise ValueError(f"pitch debe ser > 0; recibido {self.pitch!r}")
        if not np.isfinite(float(self.center_y)):
            raise ValueError("center_y debe ser finito")

    # ---- geometría ---------------------------------------------------------

    @property
    def n(self) -> int:
        return int(self.n_elements)

    @property
    def width(self) -> float:
        """Ancho total del arreglo, ``N·pitch`` [m]."""
        return self.n * float(self.pitch)

    @property
    def centers(self) -> np.ndarray:
        """``(N,)`` centros ``y_j`` de cada elemento [m]."""
        offsets = (np.arange(self.n) - 0.5 * (self.n - 1)) * float(self.pitch)
        return float(self.center_y) + offsets

    def __repr__(self) -> str:  # pragma: no cover - cosmético
        return (f"LinearArray(center_y={self.center_y * 1e3:.3f} mm, N={self.n}, "
                f"pitch={self.pitch * 1e3:.3f} mm, width={self.width * 1e3:.3f} mm)")

    # ---- aperturas ---------------------------------------------------------

    def apertures(self, amplitudes=None) -> list[Aperture]:
        """Aperturas de ancho ``pitch``, una por elemento, con las amplitudes dadas.

        ``amplitudes``: escalar o ``(N,)`` complejo; ``None`` ⇒ todas 1 (arreglo uniforme en fase).
        """
        if amplitudes is None:
            amps = np.ones(self.n, dtype=complex)
        else:
            amps = np.asarray(amplitudes, dtype=complex)
            if amps.ndim == 0:
                amps = np.full(self.n, complex(amps))
            if amps.shape != (self.n,):
                raise ValueError(f"amplitudes debe ser escalar o forma ({self.n},); "
                                 f"recibido {amps.shape}")
        return [Aperture(float(y), float(self.pitch), complex(a))
                for y, a in zip(self.centers, amps)]

    def _delays(self, phases: np.ndarray) -> list[Aperture]:
        return self.apertures(np.exp(1.0j * phases))

    def focused(self, focus_point, frequency: float, material: Material) -> list[Aperture]:
        """Aperturas con la fase que enfoca en ``focus_point = (x, y)`` [m] (``x > 0``).

        ``φ_j = −k_real·|F − y_j|``; en ``F`` todas las contribuciones llegan en fase.
        """
        if not isinstance(material, Material):
            raise TypeError("material debe ser un materials.Material")
        f = np.asarray(focus_point, dtype=float).reshape(-1)
        if f.shape != (2,):
            raise ValueError(f"focus_point debe ser (x, y); recibido {focus_point!r}")
        if not (f[0] > 0.0):
            raise ValueError(f"focus_point debe tener x > 0 (frente al baffle); recibido x={f[0]!r}")
        k_real = material.wavenumber(float(frequency)).real
        r = np.hypot(f[0], f[1] - self.centers)
        return self._delays(-k_real * r)

    def steered(self, angle: float, frequency: float, material: Material) -> list[Aperture]:
        """Aperturas con la fase que apunta a la dirección ``(cos θ, sin θ)`` [rad desde +x].

        ``φ_j = +k_real·y_j·sinθ`` — el signo que hace llegar en fase al blanco (ver el docstring
        del módulo: el plan escribe ``−k·(y_j·d)``, que apunta al ángulo espejado).

        Verificación: ``far_field_pattern`` de este arreglo tiene su máximo en ``θ`` y coincide con
        ``array_factor(..., steer_angle=θ)``.
        """
        if not isinstance(material, Material):
            raise TypeError("material debe ser un materials.Material")
        if not (abs(math.sin(float(angle))) < 1.0):
            raise ValueError(f"angle debe cumplir |θ| < π/2; recibido {angle!r} rad")
        k_real = material.wavenumber(float(frequency)).real
        return self._delays(k_real * self.centers * math.sin(float(angle)))


def array_factor(angles, n_elements: int, pitch: float, frequency: float, material: Material,
                 steer_angle: float = 0.0) -> np.ndarray:
    """``|sin(Nψ/2)/(N·sin(ψ/2))|`` para ``N`` fuentes puntuales con paso ``pitch``.

    ``ψ = k_real·pitch·(sinθ − sinθ₀)``; ``angles`` y ``steer_angle`` en rad desde broadside
    (eje ``x``). Devuelve ``(M,)`` real en ``[0, 1]``, con el límite ``ψ→0 ⇒ 1``.

    Es el patrón de un arreglo de **fuentes puntuales**: no incluye el factor de elemento (el
    ``sinc`` del pistón de ancho ``pitch``). Para el campo real de las aperturas usar
    :func:`aperture.far_field_pattern` con ``LinearArray.apertures()``.
    """
    if not isinstance(material, Material):
        raise TypeError("material debe ser un materials.Material")
    n = int(n_elements)
    if n < 1:
        raise ValueError(f"n_elements debe ser >= 1; recibido {n_elements!r}")
    if not (float(pitch) > 0.0):
        raise ValueError(f"pitch debe ser > 0; recibido {pitch!r}")
    k_real = material.wavenumber(float(frequency)).real
    theta = np.asarray(angles, dtype=float)
    psi = k_real * float(pitch) * (np.sin(theta) - math.sin(float(steer_angle)))
    half = 0.5 * psi
    # sin(N·half)/(N·sin(half)) = sinc(N·half/π)/sinc(half/π): estable en half → 0
    af = np.sinc(n * half / np.pi) / np.sinc(half / np.pi)
    return np.abs(af)


def grating_lobe_sine(pitch: float, wavelength: float, steer_angle: float = 0.0):
    """``sinθ_m = sinθ₀ + m·λ/pitch`` para los órdenes visibles ``|sinθ_m| <= 1``.

    Devuelve la lista de ``(m, sinθ)`` visibles (m ≠ 0). Útil para reportar dónde *debería* aparecer
    un lóbulo de grataje y contrastarlo con el campo medido.
    """
    if not (float(pitch) > 0.0 and float(wavelength) > 0.0):
        raise ValueError("pitch y wavelength deben ser > 0")
    s0 = math.sin(float(steer_angle))
    out = []
    m_max = int(np.ceil((1.0 + abs(s0)) * float(pitch) / float(wavelength))) + 1
    for m in range(-m_max, m_max + 1):
        if m == 0:
            continue
        s = s0 + m * float(wavelength) / float(pitch)
        if abs(s) <= 1.0:
            out.append((m, s))
    return out
