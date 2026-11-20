"""Materiales con procedencia explícita.

Convención única del repo (docs/decisions.md D4):

- tiempo ``exp(-i·ω·t)``  →  ondas salientes ``∝ exp(+i·k·r)``
- Helmholtz: ``∇²p + k²p = 0``
- atenuación por ley de potencias: ``α(f) = α0 · (f/f_ref)^n``  [Np/m]
- número de onda complejo: ``k = ω/c + i·α(f)``  (parte imaginaria positiva ⇒ ``exp(i·k·x)`` decae)

Ningún número de esta librería es inventado: cada entrada dice de dónde sale en ``provenance``.
Lo que el plan no fija queda marcado ``VERIFY`` y no se usa en validaciones.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

TWO_PI = 2.0 * math.pi

#: Tolerancia con la que se considera "sin atenuación".
_NO_ATTENUATION = 1e-12


@dataclass(frozen=True)
class Material:
    """Fluido acústico lineal (o aproximación fluida de un sólido).

    Atributos
    ---------
    name : nombre legible.
    c    : velocidad longitudinal real [m/s].
    rho  : densidad [kg/m³].
    alpha0 : atenuación [Np/m] a ``f_ref``; ``None`` o 0 ⇒ medio sin pérdidas.
    alpha_power : exponente ``n`` de la ley de potencias (``n = 1`` ⇒ lineal en f).
    f_ref : frecuencia de referencia [Hz] para ``alpha0``.
    provenance : de dónde salió el número ("plan", "literatura: ...", "VERIFY ...").
    """

    name: str
    c: float
    rho: float
    alpha0: float | None = None
    alpha_power: float = 1.0
    f_ref: float = 1.0e6
    provenance: str = ""

    def __post_init__(self) -> None:
        if not (self.c > 0.0):
            raise ValueError(f"{self.name}: c debe ser > 0 (recibido {self.c!r})")
        if not (self.rho > 0.0):
            raise ValueError(f"{self.name}: rho debe ser > 0 (recibido {self.rho!r})")
        if self.alpha0 is not None and self.alpha0 < 0.0:
            raise ValueError(f"{self.name}: alpha0 no puede ser negativo")
        if not (self.f_ref > 0.0):
            raise ValueError(f"{self.name}: f_ref debe ser > 0")

    # ---- magnitudes derivadas -------------------------------------------------

    @property
    def impedance(self) -> float:
        """Impedancia acústica característica ``Z = ρ·c`` [Pa·s/m]."""
        return self.rho * self.c

    @property
    def lossless(self) -> bool:
        return self.alpha0 is None or abs(self.alpha0) <= _NO_ATTENUATION

    def attenuation(self, frequency: float) -> float:
        """Atenuación en amplitud [Np/m] a ``frequency`` [Hz]."""
        alpha0 = self.alpha0
        if alpha0 is None or abs(alpha0) <= _NO_ATTENUATION:
            return 0.0
        if frequency < 0.0:
            raise ValueError("la frecuencia no puede ser negativa")
        return alpha0 * (frequency / self.f_ref) ** self.alpha_power

    def attenuation_db_per_m(self, frequency: float) -> float:
        """Atenuación en dB/m (= 8.686 · α[Np/m])."""
        return 8.685889638065035 * self.attenuation(frequency)

    def wavenumber(self, frequency: float) -> complex:
        """``k = ω/c + i·α(f)`` [1/m] con la convención e^{-iωt}."""
        if frequency < 0.0:
            raise ValueError("la frecuencia no puede ser negativa")
        k_real = TWO_PI * frequency / self.c
        return complex(k_real, self.attenuation(frequency))

    def wavelength(self, frequency: float) -> float:
        """Longitud de onda en el medio [m] (parte real del número de onda)."""
        if frequency <= 0.0:
            raise ValueError("la frecuencia debe ser > 0")
        return TWO_PI / self.wavenumber(frequency).real


@dataclass(frozen=True)
class ElasticSolid:
    """Sólido elástico isotrópico, definido por velocidades (Fase 2).

    Se parametriza por ``c_p`` (longitudinal) y ``c_s`` (corte) porque son las magnitudes que
    aparecen en la literatura de hueso; los módulos se derivan:
    ``μ = ρ·c_s²``, ``λ = ρ·(c_p² − 2·c_s²)``.
    """

    name: str
    rho: float
    c_p: float
    c_s: float
    provenance: str = ""

    def __post_init__(self) -> None:
        if not (self.rho > 0.0):
            raise ValueError(f"{self.name}: rho debe ser > 0")
        if not (self.c_p > 0.0 and self.c_s > 0.0):
            raise ValueError(f"{self.name}: c_p y c_s deben ser > 0")
        if self.c_s >= self.c_p:
            raise ValueError(
                f"{self.name}: c_s ({self.c_s}) debe ser menor que c_p ({self.c_p})"
            )

    @property
    def mu(self) -> float:
        """Segundo módulo de Lamé [Pa]."""
        return self.rho * self.c_s**2

    @property
    def lam(self) -> float:
        """Primer módulo de Lamé [Pa]."""
        return self.rho * (self.c_p**2 - 2.0 * self.c_s**2)

    @property
    def bulk_modulus(self) -> float:
        """Módulo volumétrico ``K = λ + 2μ/3`` [Pa]."""
        return self.lam + 2.0 * self.mu / 3.0

    @property
    def impedance_longitudinal(self) -> float:
        return self.rho * self.c_p


# ---------------------------------------------------------------------------
# Librería
# ---------------------------------------------------------------------------
# Del plan de Diego (docs/plan.md):
#   fluidos cocleares ≈ agua: c ≈ 1500 m/s, ρ ≈ 1000 kg/m³
#   hueso: c_long ≈ 3000–4000 m/s, ρ ≈ 1900–2000 kg/m³  ("ballpark, verify")
# Los valores en el medio del rango son un punto medio aritmético, no un dato medido.

WATER = Material(
    name="water",
    c=1500.0,
    rho=1000.0,
    provenance="plan: c≈1500 m/s, ρ≈1000 kg/m³ (sin pérdidas)",
)

PERILYMPH = Material(
    name="perilymph",
    c=1500.0,
    rho=1000.0,
    provenance="plan: se comporta como agua (c≈1500, ρ≈1000); atenuación no especificada",
)

ENDOLYMPH = Material(
    name="endolymph",
    c=1500.0,
    rho=1000.0,
    provenance="plan: se comporta como agua (c≈1500, ρ≈1000); atenuación no especificada",
)

SOFT_TISSUE = Material(
    name="soft_tissue",
    c=1540.0,
    rho=1050.0,
    provenance=(
        "VERIFY — el plan dice 'fluido acústico con atenuación' pero no fija c, ρ ni α. "
        "1540 m/s y 1050 kg/m³ son valores convencionales de tejido blando en literatura de "
        "ultrasonido: confirmar fuente antes de usarlos en una validación"
    ),
)

BONE_CORTICAL = Material(
    name="bone_cortical",
    c=3500.0,
    rho=1950.0,
    provenance=(
        "plan: c 3000–4000 m/s, ρ 1900–2000 kg/m³ ('ballpark, verify contra literatura'). "
        "3500 / 1950 = punto medio del rango; atenuación no especificada. "
        "Consulta a Sim-Datos (2026-11-19): no hay conjunto experimental citable específico de la "
        "cápsula ótica humana, así que esto es un SUSTITUTO ISOTRÓPICO declarado, no propiedades "
        "medidas. Para referencia, hueso temporal humano: c_P = 1752–3285 m/s a 1.03 MHz "
        "(Ammi 2008) y cortical del temporal ρ ≈ 1868–1922 kg/m³ (Peterson & Dechow 2003): "
        "3500 m/s está en el extremo alto y conviene declararlo así"
    ),
)

#: Materiales con números que el plan no fija y que no deben entrar en validaciones sin fuente.
NEEDS_VERIFICATION = frozenset({SOFT_TISSUE.name, BONE_CORTICAL.name})


def fluid_from_impedance(name: str, impedance: float, rho: float) -> Material:
    """Construye un fluido a partir de ``Z = ρ·c`` (cómodo para benchmarks de R/T)."""
    return Material(name=name, c=impedance / rho, rho=rho, provenance="derivado de Z y rho")
