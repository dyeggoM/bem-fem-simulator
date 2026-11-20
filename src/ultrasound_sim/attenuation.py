"""Leyes de atenuación del hueso, con conversión correcta y procedencia verificada.

Convención (docs/decisions.md D4): ``k = ω/c + i·α`` con ``e^{−iωt}``, así que el campo decae como
``p ∝ exp(−α·z)`` y ``α`` es un coeficiente de **amplitud** en Np/m.

**La trampa, bien dicha esta vez** (corregida tras la consulta a Sim-Datos):

- El coeficiente en dB/cm es una pérdida *en dB*, y el dB se define con 20·log para amplitud y con
  10·log para intensidad. Resultado: **la pérdida en dB/cm es el MISMO número para amplitud y para
  intensidad** (13.3 dB/cm son 13.3 dB/cm de intensidad).
- Lo que se duplica al pasar de amplitud a intensidad es el **exponente en Np**: la intensidad cae
  como ``exp(−2·α·z)`` con el mismo ``α`` en Np/m.
- Por eso: ``α[Np/m] = α[dB/cm]·100/8.6859`` es el exponente de amplitud que va al solver, y la
  pérdida en dB a distancia ``d`` es ``8.686·α·d`` — no ``2·8.686·α·d``. Duplicar los dB es un error
  de 2× que no rompe ningún test numérico; sólo miente en el informe.

Procedencia (contrastada con la fuente, ver D19 y D20):
- Pinton 2012: 13.3 ± 0.97 dB/cm a 1 MHz, lineal en f; medido como atenuación de presión.
- White 2006 (tabla 3): puntos medidos en calvarias, 14 ± 17 / 53 ± 43 / 70 ± 28 Np/m a 272 / 548 /
  840 kHz. **Discrepa de Pichardo por un factor ~5 y no se promedian.**
- Pichardo 2011: 33.8 / 231 / 334 Np/m a 270 / 836 / 1402 kHz. Son coeficientes *efectivos*
  identificados por inversión, no mediciones locales, y las muestras estaban fijadas en formol.
- **Rechazado**: ``27.31·f^1.722``, que circula como "Pichardo", no reproduce sus datos.
"""


from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .materials import Material

__all__ = [
    "DB_PER_NEPER",
    "PowerLaw",
    "MeasuredAttenuation",
    "SKULL_PINTON_2012",
    "SKULL_WHITE2006",
    "SKULL_WHITE2006_POINTS",
    "SKULL_PICHARDO_2011_POINTS",
    "SKULL_PICHARDO_2011_FIT",
    "REJECTED_MODELS",
    "db_cm_to_np_m",
    "np_m_to_db_cm",
    "skull_material",
]

#: 20·log10(e): dB por neper de AMPLITUD.
DB_PER_NEPER = 8.685889638065035


def db_cm_to_np_m(alpha_db_cm) -> np.ndarray | float:
    """dB/cm (amplitud) → Np/m."""
    out = np.asarray(alpha_db_cm, dtype=float) * 100.0 / DB_PER_NEPER
    return float(out) if out.ndim == 0 else out


def np_m_to_db_cm(alpha_np_m) -> np.ndarray | float:
    """Np/m (amplitud) → dB/cm."""
    out = np.asarray(alpha_np_m, dtype=float) * DB_PER_NEPER / 100.0
    return float(out) if out.ndim == 0 else out


@dataclass(frozen=True)
class PowerLaw:
    """``α(f) = α_ref · (f / f_ref)^n``, en Np/m de amplitud."""

    name: str
    alpha_ref: float
    exponent: float
    f_ref: float
    provenance: str
    verified: bool
    source_url: str = ""

    def __post_init__(self) -> None:
        if not (self.alpha_ref >= 0.0):
            raise ValueError(f"{self.name}: alpha_ref no puede ser negativo")
        if not (self.exponent > 0.0):
            raise ValueError(f"{self.name}: el exponente debe ser > 0")
        if not (self.f_ref > 0.0):
            raise ValueError(f"{self.name}: f_ref debe ser > 0")

    @classmethod
    def from_db_cm_per_mhz(
        cls,
        name: str,
        alpha_db_cm_at_1mhz: float,
        exponent: float,
        provenance: str,
        verified: bool,
        source_url: str = "",
    ) -> "PowerLaw":
        """Desde la forma habitual ``α = a·(f[MHz])^n`` [dB/cm]."""
        return cls(
            name=name,
            alpha_ref=float(db_cm_to_np_m(alpha_db_cm_at_1mhz)),
            exponent=float(exponent),
            f_ref=1.0e6,
            provenance=provenance,
            verified=verified,
            source_url=source_url,
        )

    def alpha(self, frequency) -> np.ndarray | float:
        """Atenuación de amplitud [Np/m] a ``frequency`` [Hz]."""
        f = np.asarray(frequency, dtype=float)
        if np.any(f < 0.0):
            raise ValueError("la frecuencia no puede ser negativa")
        out = self.alpha_ref * (f / self.f_ref) ** self.exponent
        return float(out) if out.ndim == 0 else out

    def alpha_db_cm(self, frequency) -> np.ndarray | float:
        return np_m_to_db_cm(self.alpha(frequency))

    def loss_db(self, frequency, distance: float) -> float:
        """Pérdida en dB a ``distance`` [m]: ``8.686·α·d``.

        Es la misma cifra para amplitud y para intensidad (el dB ya está definido distinto para
        cada una); lo que se duplica en intensidad es el exponente en Np, no los dB.
        """
        return float(self.alpha(frequency) * float(distance) * DB_PER_NEPER)

    def intensity_loss_db(self, frequency, distance: float) -> float:
        """Alias histórico de :meth:`loss_db`; la pérdida en dB no cambia con la magnitud."""
        return self.loss_db(frequency, distance)

    def as_material(self, name: str, c: float, rho: float) -> Material:
        """Material equivalente para el solver: ``α = α0·(f/f_ref)^n``."""
        return Material(
            name=name,
            c=float(c),
            rho=float(rho),
            alpha0=float(self.alpha_ref),
            alpha_power=float(self.exponent),
            f_ref=float(self.f_ref),
            provenance=f"atenuacion: {self.name} — {self.provenance}",
        )


@dataclass(frozen=True)
class MeasuredAttenuation:
    """Atenuación medida en frecuencias puntuales, con interpolación log-log."""

    name: str
    frequencies_hz: np.ndarray
    alphas_np_m: np.ndarray
    provenance: str
    verified: bool
    source_url: str = ""

    def __post_init__(self) -> None:
        f = np.asarray(self.frequencies_hz, dtype=float)
        a = np.asarray(self.alphas_np_m, dtype=float)
        if f.ndim != 1 or a.shape != f.shape or f.size < 2:
            raise ValueError(f"{self.name}: hacen falta >= 2 pares (f, α) alineados")
        if np.any(f <= 0.0) or np.any(a < 0.0):
            raise ValueError(f"{self.name}: frecuencias > 0 y α >= 0")
        if np.any(np.diff(f) <= 0.0):
            raise ValueError(f"{self.name}: las frecuencias deben ser crecientes")
        object.__setattr__(self, "frequencies_hz", f)
        object.__setattr__(self, "alphas_np_m", a)

    def alpha(self, frequency) -> np.ndarray | float:
        """Interpolación log-log (ley de potencia por tramos); extrapola con los extremos."""
        f = np.asarray(frequency, dtype=float)
        if np.any(f <= 0.0):
            raise ValueError("la frecuencia debe ser > 0")
        out = np.exp(
            np.interp(
                np.log(f),
                np.log(self.frequencies_hz),
                np.log(np.maximum(self.alphas_np_m, 1e-30)),
            )
        )
        # por debajo/encima del rango, np.interp mantiene el extremo: se reescala con la ley local
        lo, hi = self.frequencies_hz[0], self.frequencies_hz[-1]
        n_lo = np.log(self.alphas_np_m[1] / self.alphas_np_m[0]) / np.log(
            self.frequencies_hz[1] / self.frequencies_hz[0])
        n_hi = np.log(self.alphas_np_m[-1] / self.alphas_np_m[-2]) / np.log(
            self.frequencies_hz[-1] / self.frequencies_hz[-2])
        out = np.where(f < lo, self.alphas_np_m[0] * (f / lo) ** n_lo, out)
        out = np.where(f > hi, self.alphas_np_m[-1] * (f / hi) ** n_hi, out)
        return float(out) if out.ndim == 0 else out

    def fit_power_law(self, f_ref: float = 1.0e6) -> PowerLaw:
        """Ajuste por mínimos cuadrados en log-log a los puntos medidos."""
        ln_f = np.log(self.frequencies_hz)
        ln_a = np.log(self.alphas_np_m)
        n, ln_c = np.polyfit(ln_f, ln_a, 1)
        alpha_ref = float(np.exp(ln_c) * f_ref**n)
        return PowerLaw(
            name=f"{self.name} (ajuste propio, n={n:.2f})",
            alpha_ref=alpha_ref,
            exponent=float(n),
            f_ref=float(f_ref),
            provenance=(
                f"ajuste log-log de este repo sobre los puntos medidos de {self.name}: "
                f"{np.round(self.frequencies_hz / 1e3, 0).tolist()} kHz → "
                f"{np.round(self.alphas_np_m, 1).tolist()} Np/m. {self.provenance}"
            ),
            verified=self.verified,
            source_url=self.source_url,
        )


# ---------------------------------------------------------------------------
# Librería
# ---------------------------------------------------------------------------

SKULL_PINTON_2012 = PowerLaw.from_db_cm_per_mhz(
    name="craneo_pinton_2012",
    alpha_db_cm_at_1mhz=13.3,
    exponent=1.0,
    provenance=(
        "Pinton et al. 2012, IEEE TUFFC 39(1):299-307: atenuación medida 13.3 ± 0.97 dB/cm "
        "(presión) en muestras de ~8 mm a 1 MHz; dependencia lineal usada en la literatura "
        "transcraneal. La absorción sola es 2.7 dB/cm: el resto es reflexión/dispersión/"
        "conversión de modo"
    ),
    verified=True,
    source_url="https://pubmed.ncbi.nlm.nih.gov/22225300/",
)

SKULL_WHITE2006 = PowerLaw.from_db_cm_per_mhz(
    name="craneo_white2006_treeby2014",
    alpha_db_cm_at_1mhz=8.83,
    exponent=1.43,
    provenance=(
        "White et al. 2006 (cráneo completo) en la forma de Treeby & Cox 2014, citada en "
        "Robertson et al. 2017 Phys Med Biol: α = 8.83·f^1.43 dB/cm (formalinizado)"
    ),
    verified=True,
    source_url="https://bug.medphys.ucl.ac.uk/papers/2017-Robertson-PMB.pdf",
)

SKULL_PICHARDO_2011_POINTS = MeasuredAttenuation(
    name="craneo_pichardo_2011",
    frequencies_hz=np.array([270.0e3, 836.0e3, 1402.0e3]),
    alphas_np_m=np.array([33.8, 231.0, 334.0]),
    provenance=(
        "Pichardo, Sin & Hynynen 2011, Phys Med Biol 56:219-250: promedios medidos en cráneos "
        "humanos frescos (33.8±7.6 / 231±15 / 334±42 Np/m). Los tres puntos son el dato; "
        "cualquier ley de potencia sobre ellos es un ajuste"
    ),
    verified=True,
    source_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC3166773/",
)

SKULL_WHITE2006_POINTS = MeasuredAttenuation(
    name="craneo_white2006_medidos",
    frequencies_hz=np.array([272.0e3, 548.0e3, 840.0e3]),
    alphas_np_m=np.array([14.0, 53.0, 70.0]),
    provenance=(
        "White, Clement & Hynynen 2006, Ultrasound Med Biol 32:1085, tabla 3 (calvarias humanas "
        "fijadas, agua a 21.9±1.0 °C, coeficientes ajustados de transmisión experimental). "
        "±17 / ±43 / ±28 Np/m. OJO: discrepa de Pichardo 2011 por 2.4× a 272 kHz y 3.3× a ~838 kHz; "
        "no son intercambiables ni promediables (muestras y método de inversión distintos)"
    ),
    verified=True,
    source_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC1560344/",
)

SKULL_PICHARDO_2011_FIT = SKULL_PICHARDO_2011_POINTS.fit_power_law()
# Aviso medido, no supuesto: los tres puntos NO son colineales en log-log (270→836 kHz da n≈1.68;
# 836→1402 kHz da n≈0.72). El ajuste global se desvía hasta −20 % en 836 kHz. Para trabajo
# cuantitativo hay que usar SKULL_PICHARDO_2011_POINTS (interpolación por tramos, exacta en los
# nodos), no esta ley.


#: Modelos que circularon pero NO reproducen la fuente. Se guardan para que nadie los reintroduzca.
REJECTED_MODELS = {
    "pichardo_27.31_f^1.722": (
        "no reproduce los datos de Pichardo 2011: a 1.4 MHz da 49 dB/cm frente a 29 medidos "
        "(coincide sólo en el extremo de 270 kHz). El ajuste de sus tres puntos da n ≈ 1.39"
    ),
    "ilustrativo_20_f": (
        "el propio documento que lo propone lo declara modelo ilustrativo, no una curva del "
        "cráneo cortical"
    ),
}


def skull_material(
    law: PowerLaw,
    name: str | None = None,
    c: float = 2800.0,
    rho: float = 1900.0,
) -> Material:
    """Material de cráneo homogéneo equivalente con la atenuación dada.

    ``c = 2800 m/s`` y ``rho = 1900 kg/m³`` son los valores de *cráneo homogéneo equivalente* usados
    en simulación transcraneal (distintos de la cortical pura, c ≈ 3500). Documentado acá porque
    mezclar los dos juegos de números es un error fácil de cometer.
    """
    m = law.as_material(name or f"craneo_{law.name}", c=c, rho=rho)
    return Material(
        name=m.name,
        c=m.c,
        rho=m.rho,
        alpha0=m.alpha0,
        alpha_power=m.alpha_power,
        f_ref=m.f_ref,
        provenance=(
            m.provenance
            + f". c={c} m/s y rho={rho} kg/m³ son cráneo homogéneo equivalente "
            "(no cortical pura)"
        ),
    )
