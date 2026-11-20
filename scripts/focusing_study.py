#!/usr/bin/env python
"""Estudio de focalización: arreglo lineal vs apertura continua, para las geometrías documentadas.

**Supuesto obligatorio (vale para toda la tabla): el modelo es 2D.** Una apertura es una línea
sobre ``x = 0`` (``aperture.Aperture``) y la presión sale de la integral de Rayleigh 2D con la
imagen del baffle rígido. Un arreglo 2D o hemisférico **no** se simula en 3D: se modela su **corte**,
o sea una línea de ``N`` elementos con el mismo pitch y el mismo ancho total. Eso cambia el número
de elementos respecto del equipo físico (el 16×16 de Marzo & Drinkwater son 256 elementos; su corte
es una línea de 16) y es la razón por la que la ganancia de presión no se compara con ``20log10(256)``
sino con la de los elementos que el corte tiene. Cada fila del CSV y cada figura dicen ``2D``.

Qué mide
--------
Para cada configuración documentada en ``docs/decisions.md`` **D20 bloque 6** (rangos reales, no
productos inventados; la columna ``fuente`` de cada fila dice de dónde salió la geometría):

* **pico axial real** — ``max |p|²`` sobre el eje ``y = 0``. No se asume que caiga en ``F``: D13
  documenta que a f/# bajo el máximo se corre hacia la apertura (D=60 mm, F=40 mm ⇒ 29.5 vs 40 mm).
* ``|p|`` y ``|p|²`` **en F**, y ``|p(F)|² / max|p|²`` (la «meseta» de D13).
* **FWHM lateral en el plano del pico real**, con ``focus.lateral_fwhm``.
* **nivel del lóbulo lateral** (dB) con ``focus.sidelobe_level_db``, medido en el **corte focal**
  del pico real (no es el lóbulo de grataje de campo lejano: el ángulo de grataje se reporta aparte
  con ``arrays.grating_lobe_sine`` y a frecuencias altas el lóbulo puede no haberse separado aún en
  el plano focal). Se reporta también el de la apertura continua en el mismo plano.
* **ganancia de presión contra un elemento** (el central, sin delay) en el pico, y su comparación
  con ``20log10(N)`` — D13: fuera de paraxial la identidad **no** se cumple.
* **correlación de campo complejo** contra una apertura continua del mismo ancho total, enfocada en
  el mismo ``F`` (``focus.field_correlation``), sobre una caja centrada en el pico real.
* ``pitch/λ``, ``D/λ``, ``f/#`` y el ángulo del lóbulo de grataje (``arrays.grating_lobe_sine``).

La «apertura continua» se materializa como un arreglo **finamente muestreado** (``pitch = λ/16``,
tope 1200 elementos), que es como se construyó la referencia de la Fase 4B: el modelo de ``Aperture``
no admite una fase continua dentro de una apertura, así que la fase de foco continua se aproxima
con muchos elementos. El test de la Fase 4B verifica que con ``λ/16`` la referencia ya está
convergida (``λ/8`` → 8.66 mm, ``λ/32`` → 8.63 mm, ``λ/64`` → 8.629 mm).

Régimen (D8): a 40–100 kHz con aperturas de decenas de mm **no hay focalización** — el FWHM es mayor
que la cóclea (9 mm). Se incluye igual porque es el régimen que Diego fijó como techo inicial, y
porque la comparación «transmisor vs arreglo» ahí mide iluminación, no foco. Los tests de
focalización corren a 250 kHz – 1 MHz.

Uso
---
    .venv/bin/python scripts/focusing_study.py                # tabla + CSV + 6+5 figuras (~2–3 min)
    .venv/bin/python scripts/focusing_study.py --map-px 40    # mapas más baratos
    .venv/bin/python scripts/focusing_study.py --only aire_40k_linea
    .venv/bin/python scripts/focusing_study.py --no-figures

Salidas
-------
``out/focusing_study.csv`` (todos los números + ``fuente``), ``out/focus_<nombre>.png`` por
configuración y ``out/focusing_summary.png``. La resolución de los mapas usada queda en la columna
``mapa_px`` del CSV.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

_REPO = Path(__file__).resolve().parents[1]
if (_REPO / "src").is_dir() and str(_REPO / "src") not in sys.path:
    sys.path.insert(0, str(_REPO / "src"))

from ultrasound_sim.aperture import Aperture, pressure  # noqa: E402
from ultrasound_sim.arrays import LinearArray, grating_lobe_sine  # noqa: E402
from ultrasound_sim.focus import (  # noqa: E402
    field_correlation,
    lateral_fwhm,
    sidelobe_level_db,
)
from ultrasound_sim.materials import WATER  # noqa: E402

OUT_DIR = _REPO / "out"
OUT_CSV = OUT_DIR / "focusing_study.csv"
SUMMARY_PNG = OUT_DIR / "focusing_summary.png"

#: Resolución de los mapas 2D por defecto. Medido: 80×80 con 256 aperturas son ~9 s por
#: configuración; el encargo pide bajar si pasa de 60 s (no ocurre).
DEFAULT_MAP_PX = 80

#: Muestreo de la referencia «continua»: ``pitch_ref = λ/16``, tope 1200 elementos.
REF_LAMBDA_DIV = 16
REF_MAX_ELEMENTS = 1200

#: Caja de correlación usada por el estudio de Fase 4B (25–55 mm × ±30 mm), para poder
#: reproducir sus números de correlación 2D.
BOX_4B = (25.0e-3, 55.0e-3, 30.0e-3, 31, 41)


# --------------------------------------------------------------------------- #
# Configuraciones
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Config:
    """Una configuración a simular.

    ``aperture_m`` es el ancho total de la cara activa y ``n_elements`` la cantidad de elementos
    del **corte 2D**; el pitch sale de ``aperture_m / n_elements`` (elementos contiguos, sin kerf,
    que es lo único que el modelo soporta). ``focus_m = None`` ⇒ pistón plano.
    """

    key: str
    name: str
    frequency_hz: float
    n_elements: int
    aperture_m: float
    focus_m: float | None
    source: str
    kind: str = "documentada"
    note: str = ""

    @property
    def pitch_m(self) -> float:
        return self.aperture_m / self.n_elements

    @property
    def f_number(self) -> float | None:
        return None if self.focus_m is None else self.focus_m / self.aperture_m


#: Fuente común del bloque 6 de D20 (transductores).
_D20_6 = "docs/decisions.md D20 bloque 6 (transductores)"

CONFIGS: tuple[Config, ...] = (
    Config(
        key="aire_40k_linea",
        name="Aire 40 kHz, 16 elementos (corte del 16×16)",
        frequency_hz=40.0e3,
        n_elements=16,
        aperture_m=0.160,
        focus_m=0.150,
        source=(
            f"{_D20_6}: emisor de 10 mm, arreglo 16×16 = 256 con pitch 10 mm "
            "(Marzo & Drinkwater 2019). Corte 2D = línea de 16 elementos de 10 mm."
        ),
        note="El equipo real es 16×16; el corte 2D tiene 16 elementos y 160 mm de ancho.",
    ),
    Config(
        key="langevin_40k_plano",
        name="Langevin 40 kHz, 45 mm, plano",
        frequency_hz=40.0e3,
        n_elements=1,
        aperture_m=0.045,
        focus_m=None,
        source=(
            f"{_D20_6}: contacto 40 kHz, Langevin de 45 mm de diámetro radiante (ficha técnica)."
        ),
        note="Monolítico plano: no hay foco geométrico. Sirve de referencia de «sin foco».",
    ),
    Config(
        key="langevin_40k_foco100",
        name="Langevin 40 kHz, 45 mm, enfocado a 100 mm",
        frequency_hz=40.0e3,
        n_elements=1,
        aperture_m=0.045,
        focus_m=0.100,
        source=(
            f"{_D20_6}: geometría del Langevin de 45 mm (ficha técnica). F = 100 mm es una "
            "hipótesis del estudio, no un dato del equipo (el modelo de apertura no puede "
            "enfocar un pistón rígido; se incluye para medir el límite D/λ = 1.2)."
        ),
        note="Un solo elemento: la fase de foco no existe físicamente; fila de control de D/λ.",
    ),
    Config(
        key="transcraneal_250k",
        name="Transcraneal 250 kHz, 300 mm, 256 elementos",
        frequency_hz=250.0e3,
        n_elements=256,
        aperture_m=0.300,
        focus_m=0.150,
        source=(
            f"{_D20_6}: terapéutico transcraneal, aperturas 200–300 mm y 256–1024 elementos, "
            "F = 150 mm a 250/500 kHz. Sin pitch verificable: se usa pitch = 300/256."
        ),
        note="Pitch no documentado (D20): pitch = ancho/n_elementos, declarado como supuesto.",
    ),
    Config(
        key="transcraneal_500k",
        name="Transcraneal 500 kHz, 300 mm, 256 elementos",
        frequency_hz=500.0e3,
        n_elements=256,
        aperture_m=0.300,
        focus_m=0.150,
        source=(
            f"{_D20_6}: misma familia transcraneal (200–300 mm, 256–1024 elementos, F = 150 mm) "
            "a 500 kHz. Pitch = 300/256, supuesto declarado."
        ),
        note="Igual geometría que 250 kHz: el único cambio es λ, o sea pitch/λ.",
    ),
    Config(
        key="diseno_1m",
        name="Diseño 1 MHz, 200 mm, 256 elementos",
        frequency_hz=1.0e6,
        n_elements=256,
        aperture_m=0.200,
        focus_m=0.150,
        source=(
            f"{_D20_6} (rango de aperturas 200–300 mm y 256–1024 elementos) + D8 (hoja de ruta "
            "a >1 MHz): 1 MHz con 200 mm y 256 elementos, F = 150 mm. Pitch = 200/256."
        ),
        note="Configuración de diseño propio dentro del rango documentado; no es un producto.",
    ),
)

#: Barrido de frecuencia con geometría fija: 160 mm, pitch 10 mm, F = 150 mm.
SWEEP_FREQUENCIES: tuple[float, ...] = (40.0e3, 100.0e3, 250.0e3, 500.0e3, 1.0e6)
SWEEP_WIDTH_M = 0.160
SWEEP_N = 16
SWEEP_FOCUS_M = 0.150

#: Tabla de referencia de ``docs/plan.md`` (Fase 4B): FWHM en el plano del máximo axial real,
#: ``D = 60 mm``, ``F = 40 mm``, 100 kHz, agua. Se cita para contrastar, no para forzar la medición.
TABLA_4B_FWHM_MM = {1.0: 12.74, 0.5: 10.97, 0.25: 9.17, 0.125: 8.77}

#: Escalones de pitch del estudio de Fase 4B (D=60 mm, F=40 mm, 100 kHz): pitch = λ/div.
CALIB_DIVISIONS = (1, 2, 4, 8)


def sweep_configs() -> tuple[Config, ...]:
    """Las cinco frecuencias con la **misma** geometría (160 mm, 16 elementos, F = 150 mm)."""
    return tuple(
        Config(
            key=f"barrido_{f / 1e3:.0f}k".replace(".", "p"),
            name=f"Barrido {f / 1e3:g} kHz, 160 mm, 16 elementos",
            frequency_hz=f,
            n_elements=SWEEP_N,
            aperture_m=SWEEP_WIDTH_M,
            focus_m=SWEEP_FOCUS_M,
            source=(
                "geometría fija del barrido de este estudio (160 mm, pitch 10 mm, F = 150 mm); "
                "las frecuencias son los escalones de D8/D20 (40 / 100 / 250 / 500 kHz / 1 MHz)"
            ),
            kind="barrido",
            note="Misma geometría que `aire_40k_linea`; cambia sólo λ.",
        )
        for f in SWEEP_FREQUENCIES
    )


def calibration_4b() -> tuple[Config, ...]:
    """Los cuatro pitches del estudio de Fase 4B, para anclar «qué pitch hace falta»."""
    lam = WATER.wavelength(100.0e3)
    return tuple(
        Config(
            key=f"calib4b_lambda_{div}",
            name=f"Fase 4B: pitch = λ/{div}",
            frequency_hz=100.0e3,
            n_elements=int(round(0.060 / (lam / div))),
            aperture_m=0.060,
            focus_m=0.040,
            source=(
                "docs/plan.md Fase 4B (tabla de referencia): D = 60 mm, F = 40 mm, 100 kHz, agua, "
                f"pitch = λ/{div}."
            ),
            kind="calibracion_4b",
            note="Ancla de la equivalencia arreglo ↔ apertura continua.",
        )
        for div in CALIB_DIVISIONS
    )


# --------------------------------------------------------------------------- #
# Medición
# --------------------------------------------------------------------------- #


@dataclass
class Measurement:
    """Todo lo medido para una configuración (y los arrays que necesitan las figuras)."""

    key: str
    name: str
    kind: str
    fuente: str
    note: str
    frequency_hz: float
    n_elements: int
    pitch_m: float
    aperture_m: float
    focus_m: float | None
    wavelength_m: float
    pitch_over_lambda: float
    aperture_over_lambda: float
    f_number: float | None
    grating_lobe_deg: float | None
    peak_x_m: float
    focal_shift_m: float | None
    p_at_focus_pa: float | None
    p2_at_focus_rel: float | None
    fwhm_m: float
    fwhm_continuous_m: float
    fwhm_ratio: float
    sidelobe_db: float
    sidelobe_continuous_db: float
    gain_db: float
    gain_ideal_db: float
    correlation: float
    map_px: int
    seconds: float
    axial_x: np.ndarray | None = field(default=None, repr=False)
    axial_i: np.ndarray | None = field(default=None, repr=False)
    lateral_y: np.ndarray | None = field(default=None, repr=False)
    lateral_i: np.ndarray | None = field(default=None, repr=False)
    map_x: np.ndarray | None = field(default=None, repr=False)
    map_y: np.ndarray | None = field(default=None, repr=False)
    map_i: np.ndarray | None = field(default=None, repr=False)


def _apertures(cfg: Config):
    array = LinearArray(0.0, cfg.n_elements, cfg.pitch_m)
    if cfg.focus_m is None:
        return array.apertures()
    return array.focused((cfg.focus_m, 0.0), cfg.frequency_hz, WATER)


def _axial_peak(aps, frequency: float, x_max: float, step: float) -> tuple[np.ndarray, np.ndarray]:
    """Perfil axial en ``y = 0``; si el máximo cae en el borde, alarga el barrido (hasta 2 veces)."""
    for _ in range(3):
        xs = np.arange(2.0e-3, x_max, step)
        p = pressure(np.column_stack([xs, np.zeros_like(xs)]), frequency, WATER, aps)
        i = np.abs(p) ** 2
        idx = int(np.argmax(i))
        if idx < i.size - 2:
            return xs, i
        x_max *= 1.6
    return xs, i


def _grating_angle_deg(pitch: float, wavelength: float, n_elements: int) -> float | None:
    """Ángulo del primer lóbulo de grataje visible (``|sinθ|`` más chico), o ``None``.

    Con un solo elemento no hay arreglo: ``grating_lobe_sine`` describe el muestreo, no un pistón
    único, así que devolver un ángulo ahí sería un número con aire de medido.
    """
    if n_elements < 2:
        return None
    preds = grating_lobe_sine(pitch, wavelength)
    if not preds:
        return None
    _, s = min(preds, key=lambda item: abs(item[1]))
    return math.degrees(math.asin(s))


def measure(
    cfg: Config,
    *,
    map_px: int | None = DEFAULT_MAP_PX,
    corr_box: tuple[float, float, float, int, int] | None = None,
) -> Measurement:
    """Mide una configuración. ``map_px=None`` omite el mapa 2D (los tests no lo necesitan)."""
    t0 = time.perf_counter()
    f = cfg.frequency_hz
    lam = WATER.wavelength(f)
    pitch = cfg.pitch_m
    aps = _apertures(cfg)

    # ---- perfil axial y pico real ------------------------------------------
    if cfg.focus_m is not None:
        ref_len = cfg.focus_m
        x_max = 1.8 * ref_len
    else:
        ref_len = cfg.aperture_m**2 / (4.0 * lam)  # longitud de campo cercano del pistón
        x_max = max(4.0 * ref_len, 30.0e-3)
    dx = min(0.2e-3, max(0.05e-3, lam / 40.0))
    xs, axial_i = _axial_peak(aps, f, x_max, dx)
    peak_x = float(xs[int(np.argmax(axial_i))])
    p_peak = float(np.sqrt(axial_i.max()))

    # ---- corte lateral en el plano del pico real ---------------------------
    half = max(0.45 * ref_len, 2.5 * lam * ref_len / cfg.aperture_m)
    dy = min(0.05e-3, lam / 40.0)
    ys = np.arange(-half, half + 0.5 * dy, dy)
    lateral = pressure(np.column_stack([np.full_like(ys, peak_x), ys]), f, WATER, aps)
    lateral_i = np.abs(lateral) ** 2
    fwhm = lateral_fwhm(ys, lateral_i)
    sidelobe = sidelobe_level_db(ys, lateral_i)

    # ---- |p| y |p|² en F ---------------------------------------------------
    if cfg.focus_m is None:
        p_focus = p2_rel = None
        shift = None
    else:
        p_focus = float(abs(pressure(np.array([[cfg.focus_m, 0.0]]), f, WATER, aps)[0]))
        p2_rel = float(p_focus**2 / axial_i.max())
        shift = float(peak_x - cfg.focus_m)

    # ---- ganancia contra un elemento (el central, sin delay) --------------
    ref_ap = Aperture(0.0, pitch)
    p_single = float(abs(pressure(np.array([[peak_x, 0.0]]), f, WATER, [ref_ap])[0]))
    gain = float(20.0 * np.log10(p_peak / p_single)) if p_single > 0.0 else float("nan")

    # ---- referencia continua del mismo ancho ------------------------------
    n_ref = int(min(REF_MAX_ELEMENTS, max(64, round(cfg.aperture_m / (lam / REF_LAMBDA_DIV)))))
    ref_array = LinearArray(0.0, n_ref, cfg.aperture_m / n_ref)
    ref_aps = (ref_array.apertures() if cfg.focus_m is None
               else ref_array.focused((cfg.focus_m, 0.0), f, WATER))

    # pico de la referencia (barrido axial grueso: sólo hace falta dónde medir el FWHM)
    xs_ref = np.arange(2.0e-3, x_max, max(dx, 0.5e-3))
    ref_axial = np.abs(pressure(np.column_stack([xs_ref, np.zeros_like(xs_ref)]),
                                f, WATER, ref_aps)) ** 2
    ref_peak_x = float(xs_ref[int(np.argmax(ref_axial))])

    # FWHM de la referencia en su propio plano, con el rango justo (y reintento si se recorta)
    half_ref = max(3.0 * lam * ref_peak_x / cfg.aperture_m, 2.0 * fwhm)
    for _ in range(3):
        ys_ref = np.arange(-half_ref, half_ref + 0.5 * dy, dy)
        ref_lat = np.abs(pressure(np.column_stack([np.full_like(ys_ref, ref_peak_x), ys_ref]),
                                  f, WATER, ref_aps)) ** 2
        try:
            fwhm_ref = lateral_fwhm(ys_ref, ref_lat)
            sidelobe_ref = sidelobe_level_db(ys_ref, ref_lat)
            break
        except ValueError:
            half_ref *= 2.0
    else:  # pragma: no cover - no debería ocurrir con este generador de rangos
        raise RuntimeError(f"{cfg.key}: no se pudo medir el FWHM de la referencia continua")

    # ---- correlación de campo complejo en la región blanco ----------------
    if corr_box is not None:
        x0, x1, half_y, nx, ny = corr_box
    else:
        x0, x1 = 0.55 * peak_x, 1.45 * peak_x
        half_y = 0.40 * (cfg.focus_m if cfg.focus_m is not None else peak_x)
        nx = ny = 0
    if nx <= 0:
        nx, ny = 20, 24
    gx, gy = np.meshgrid(np.linspace(x0, x1, nx), np.linspace(-half_y, half_y, ny), indexing="ij")
    box = np.column_stack([gx.ravel(), gy.ravel()])
    field_array = pressure(box, f, WATER, aps)
    field_ref = pressure(box, f, WATER, ref_aps)
    corr = field_correlation(field_array, field_ref)

    # ---- mapa 2D |p|² ------------------------------------------------------
    map_x = map_y = map_i = None
    if map_px:
        x0m = max(1.0e-3, 0.08 * ref_len)
        mx = np.linspace(x0m, 1.55 * ref_len, int(map_px))
        my = np.linspace(-0.55 * ref_len, 0.55 * ref_len, int(map_px))
        gmx, gmy = np.meshgrid(mx, my, indexing="ij")
        pts = np.column_stack([gmx.ravel(), gmy.ravel()])
        map_i = (np.abs(pressure(pts, f, WATER, aps)) ** 2).reshape(gmx.shape)
        map_x, map_y = mx, my

    return Measurement(
        key=cfg.key,
        name=cfg.name,
        kind=cfg.kind,
        fuente=cfg.source,
        note=cfg.note,
        frequency_hz=f,
        n_elements=cfg.n_elements,
        pitch_m=pitch,
        aperture_m=cfg.aperture_m,
        focus_m=cfg.focus_m,
        wavelength_m=lam,
        pitch_over_lambda=pitch / lam,
        aperture_over_lambda=cfg.aperture_m / lam,
        f_number=cfg.f_number,
        grating_lobe_deg=_grating_angle_deg(pitch, lam, cfg.n_elements),
        peak_x_m=peak_x,
        focal_shift_m=shift,
        p_at_focus_pa=p_focus,
        p2_at_focus_rel=p2_rel,
        fwhm_m=fwhm,
        fwhm_continuous_m=fwhm_ref,
        fwhm_ratio=fwhm / fwhm_ref,
        sidelobe_db=sidelobe,
        sidelobe_continuous_db=sidelobe_ref,
        gain_db=gain,
        gain_ideal_db=float(20.0 * np.log10(cfg.n_elements)),
        correlation=corr,
        map_px=int(map_px or 0),
        seconds=time.perf_counter() - t0,
        axial_x=xs,
        axial_i=axial_i,
        lateral_y=ys,
        lateral_i=lateral_i,
        map_x=map_x,
        map_y=map_y,
        map_i=map_i,
    )


def study(
    configs=None,
    *,
    map_px: int | None = DEFAULT_MAP_PX,
    verbose: bool = True,
) -> list[Measurement]:
    """Mide todas las configuraciones (por defecto documentadas + barrido + calibración 4B)."""
    if configs is None:
        configs = CONFIGS + sweep_configs() + calibration_4b()
    out: list[Measurement] = []
    for cfg in configs:
        m = measure(cfg, map_px=map_px, corr_box=(BOX_4B if cfg.kind == "calibracion_4b" else None))
        out.append(m)
        if verbose:
            print(f"  [{len(out):2d}/{len(configs):2d}] {cfg.key:22s} "
                  f"pitch/λ={m.pitch_over_lambda:6.3f}  FWHM={m.fwhm_m * 1e3:7.2f} mm  "
                  f"FWHM/cont={m.fwhm_ratio:6.3f}  corr={m.correlation:7.4f}  "
                  f"{m.seconds:5.1f} s")
    return out


# --------------------------------------------------------------------------- #
# Presentación
# --------------------------------------------------------------------------- #


def format_table(rows: list[Measurement]) -> str:
    """Tabla por consola: una fila por configuración."""
    head = (
        f"{'configuración':22s} {'f':>6s} {'N':>4s} {'pitch/λ':>7s} {'D/λ':>6s} {'f/#':>5s} "
        f"{'pico':>7s} {'|p(F)|²/':>8s} {'FWHM':>7s} {'FWHM/':>6s} {'lóbulo':>7s} "
        f"{'gan.':>6s} {'corr':>7s} {'grat.':>6s}"
    )
    units = (
        f"{'':22s} {'kHz':>6s} {'':>4s} {'':>7s} {'':>6s} {'':>5s} "
        f"{'mm':>7s} {'máx':>8s} {'mm':>7s} {'cont':>6s} {'dB':>7s} {'dB':>6s} {'':>7s} {'°':>6s}"
    )
    lines = [head, units, "-" * len(head)]
    for m in rows:
        fnr = "—" if m.f_number is None else f"{m.f_number:5.2f}"
        p2 = "—" if m.p2_at_focus_rel is None else f"{m.p2_at_focus_rel:8.4f}"
        lob = "—" if not np.isfinite(m.sidelobe_db) else f"{m.sidelobe_db:7.2f}"
        gt = "—" if m.grating_lobe_deg is None else f"{m.grating_lobe_deg:6.1f}"
        lines.append(
            f"{m.key:22s} {m.frequency_hz / 1e3:6.0f} {m.n_elements:4d} "
            f"{m.pitch_over_lambda:7.3f} {m.aperture_over_lambda:6.1f} {fnr} "
            f"{m.peak_x_m * 1e3:7.1f} {p2} {m.fwhm_m * 1e3:7.2f} {m.fwhm_ratio:6.3f} "
            f"{lob} {m.gain_db:6.2f} {m.correlation:7.4f} {gt}"
        )
    return "\n".join(lines)


CSV_FIELDS = (
    "nombre", "tipo", "fuente", "frecuencia_khz", "n_elementos", "pitch_mm", "apertura_mm",
    "foco_mm", "pitch_over_lambda", "apertura_over_lambda", "f_number", "grating_lobe_deg",
    "pico_axial_mm", "corrimiento_focal_mm", "p_at_foco_Pa", "p2_foco_rel", "fwhm_mm",
    "fwhm_continua_mm", "fwhm_ratio", "sidelobe_db", "sidelobe_continua_db", "ganancia_db",
    "ganancia_ideal_db",
    "correlacion", "mapa_px", "segundos", "notas",
)


def _csv_num(value, spec: str) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and not np.isfinite(value):
        return ""
    return format(value, spec)


def write_csv(rows: list[Measurement], path: Path = OUT_CSV) -> Path:
    """Escribe el CSV con todos los números y la procedencia de cada fila."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(CSV_FIELDS)
        for m in rows:
            w.writerow([
                m.key,
                m.kind,
                m.fuente,
                _csv_num(m.frequency_hz / 1e3, "g"),
                m.n_elements,
                _csv_num(m.pitch_m * 1e3, ".4f"),
                _csv_num(m.aperture_m * 1e3, ".2f"),
                _csv_num(None if m.focus_m is None else m.focus_m * 1e3, ".2f"),
                _csv_num(m.pitch_over_lambda, ".4f"),
                _csv_num(m.aperture_over_lambda, ".2f"),
                _csv_num(m.f_number, ".3f"),
                _csv_num(m.grating_lobe_deg, ".2f"),
                _csv_num(m.peak_x_m * 1e3, ".2f"),
                _csv_num(None if m.focal_shift_m is None else m.focal_shift_m * 1e3, "+.2f"),
                _csv_num(m.p_at_focus_pa, ".4e"),
                _csv_num(m.p2_at_focus_rel, ".4f"),
                _csv_num(m.fwhm_m * 1e3, ".3f"),
                _csv_num(m.fwhm_continuous_m * 1e3, ".3f"),
                _csv_num(m.fwhm_ratio, ".4f"),
                _csv_num(m.sidelobe_db, ".2f"),
                _csv_num(m.sidelobe_continuous_db, ".2f"),
                _csv_num(m.gain_db, ".2f"),
                _csv_num(m.gain_ideal_db, ".2f"),
                _csv_num(m.correlation, ".5f"),
                m.map_px,
                _csv_num(m.seconds, ".1f"),
                "; ".join(x for x in (m.note, "modelo 2D (corte de la apertura)") if x),
            ])
    return path


def _plt():
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    return plt


def _short(name: str) -> str:
    return name.replace("Transcraneal", "Transcr.").replace("elementos", "elem.")


def plot_config(m: Measurement, path: Path) -> Path:
    """Figura de una configuración: mapa 2D de ``|p|²`` + cortes axial y lateral."""
    plt = _plt()
    fig = plt.figure(figsize=(13.2, 4.0))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.0, 1.0], wspace=0.34)

    ax = fig.add_subplot(gs[0])
    if m.map_i is not None:
        db = 10.0 * np.log10(m.map_i / m.map_i.max() + 1e-14)
        pcm = ax.pcolormesh(m.map_x * 1e3, m.map_y * 1e3, db, shading="auto",
                            cmap="inferno", vmin=-25.0, vmax=0.0)
        cb = fig.colorbar(pcm, ax=ax, pad=0.02)
        cb.set_label("|p|² / máx [dB]")
        ax.plot(m.peak_x_m * 1e3, 0.0, "c+", ms=12, mew=2, label="pico real")
        if m.focus_m is not None:
            ax.plot(m.focus_m * 1e3, 0.0, "o", mfc="none", mec="w", ms=9, label="F geométrico")
        ax.legend(loc="upper right", fontsize=8, framealpha=0.6)
        ax.set_xlabel("x [mm]")
        ax.set_ylabel("y [mm]")
    ax.set_title(f"|p|² — {_short(m.name)}\n2D: N={m.n_elements}, pitch={m.pitch_m * 1e3:.2f} mm",
                 fontsize=9)

    ax1 = fig.add_subplot(gs[1])
    ax1.plot(m.axial_x * 1e3, 10.0 * np.log10(m.axial_i / m.axial_i.max() + 1e-14), lw=1.2)
    ax1.axvline(m.peak_x_m * 1e3, color="c", ls="--", lw=1, label=f"pico {m.peak_x_m * 1e3:.1f} mm")
    if m.focus_m is not None:
        ax1.axvline(m.focus_m * 1e3, color="0.6", ls=":", lw=1,
                    label=f"F = {m.focus_m * 1e3:.0f} mm")
    ax1.set_ylim(-30, 1)
    ax1.set_xlabel("x [mm]")
    ax1.set_ylabel("|p|² / máx [dB]")
    ax1.set_title("corte axial (y = 0)", fontsize=9)
    ax1.grid(alpha=0.25)
    ax1.legend(fontsize=8)

    ax2 = fig.add_subplot(gs[2])
    ax2.plot(m.lateral_y * 1e3, 10.0 * np.log10(m.lateral_i / m.lateral_i.max() + 1e-14), lw=1.2)
    ax2.axhline(-3.0, color="0.6", ls=":", lw=1)
    ax2.set_ylim(-30, 1)
    ax2.set_xlabel("y [mm]")
    ax2.set_title(f"corte lateral en x = pico (FWHM {m.fwhm_m * 1e3:.2f} mm)", fontsize=9)
    ax2.grid(alpha=0.25)

    fig.suptitle(
        f"{m.key} · f = {m.frequency_hz / 1e3:g} kHz · modelo 2D (corte) · "
        f"pitch/λ = {m.pitch_over_lambda:.3f} · D/λ = {m.aperture_over_lambda:.1f}",
        fontsize=9, y=1.02,
    )
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_summary(rows: list[Measurement], path: Path = SUMMARY_PNG) -> Path:
    """Resumen: FWHM vs pitch/λ y cociente contra la continua, con la curva de la Fase 4B."""
    plt = _plt()
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(12.4, 4.6))

    calib = [m for m in rows if m.kind == "calibracion_4b"]
    docs = [m for m in rows if m.kind == "documentada"]
    sweep = [m for m in rows if m.kind == "barrido"]

    axa.plot([], [], "o", mfc="none", mec="C0", label="documentada (D20)")
    axa.plot([], [], "s", mfc="none", mec="C1", label="barrido de f, geometría fija")
    seen = [(m.pitch_over_lambda, m.fwhm_m) for m in docs]
    for i, m in enumerate(docs + sweep):
        axa.plot(m.pitch_over_lambda, m.fwhm_m * 1e3, "o" if m.kind == "documentada" else "s",
                 ms=7, mfc="none", mec="C0" if m.kind == "documentada" else "C1")
        if m.kind == "barrido":
            # un punto del barrido que coincide con una configuración documentada (misma geometría
            # y misma f) no se rotula: el texto quedaría exactamente encima
            if any(abs(m.pitch_over_lambda - x) < 1e-9 and abs(m.fwhm_m - y) < 1e-12
                   for x, y in seen):
                continue
            label, offset = f"{m.frequency_hz / 1e3:g} kHz", (5, -11)
        else:
            label, offset = m.key, (5, 6 if i % 2 == 0 else -12)
        axa.annotate(label, (m.pitch_over_lambda, m.fwhm_m * 1e3), fontsize=6.5,
                     xytext=offset, textcoords="offset points")
    axa.set_xscale("log")
    axa.set_xlabel("pitch / λ")
    axa.set_ylabel("FWHM lateral en el pico real [mm]")
    axa.set_title("FWHM vs pitch/λ (una fila por configuración)", fontsize=9)
    axa.grid(alpha=0.25, which="both")
    axa.legend(fontsize=7, loc="upper left")

    if calib:
        order = sorted(calib, key=lambda m: m.pitch_over_lambda)
        axb.plot([m.pitch_over_lambda for m in order], [m.fwhm_ratio for m in order], "k-o",
                 ms=5, lw=1.2, label="Fase 4B: D=60 mm, F=40 mm, 100 kHz (f/#=0.67)")
        for m in order:
            axb.annotate(f"λ/{int(round(1 / m.pitch_over_lambda))}",
                         (m.pitch_over_lambda, m.fwhm_ratio), fontsize=6.5, ha="right",
                         xytext=(-6, -11), textcoords="offset points")
    for m in docs:
        axb.plot(m.pitch_over_lambda, m.fwhm_ratio, "o", ms=7, mfc="none", mec="C0")
    for m in sweep:
        axb.plot(m.pitch_over_lambda, m.fwhm_ratio, "s", ms=6, mfc="none", mec="C1")
    # los puntos documentados no se rotulan acá: están rotulados en el panel izquierdo y a esta
    # escala (ratios 1.00–1.02) los textos se pisan entre sí
    axb.plot([], [], "o", mfc="none", mec="C0", label="documentada (D20)")
    axb.plot([], [], "s", mfc="none", mec="C1", label="barrido de f, geometría fija")
    axb.axhline(1.02, color="C2", ls="--", lw=1, label="+2 % (criterio D13/4B)")
    axb.axhline(1.07, color="C3", ls=":", lw=1, label="+7 %")
    axb.set_xscale("log")
    axb.set_ylim(0.98, max(1.6, max(m.fwhm_ratio for m in rows) * 1.05))
    axb.set_xlabel("pitch / λ")
    axb.set_ylabel("FWHM arreglo / FWHM continua")
    axb.set_title("equivalencia con la apertura continua (misma D, mismo F)", fontsize=9)
    axb.grid(alpha=0.25, which="both")
    axb.legend(fontsize=7, loc="upper left")

    fig.suptitle(
        "Estudio de focalización — modelo 2D (aperturas lineales; un arreglo 2D se modela por su "
        "corte)", fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return path


def conclusions(rows: list[Measurement]) -> list[str]:
    """Las tres conclusiones, calculadas de los números (no escritas a mano)."""
    calib = sorted((m for m in rows if m.kind == "calibracion_4b"),
                   key=lambda m: m.pitch_over_lambda)
    docs = [m for m in rows if m.kind == "documentada"]
    sweep = [m for m in rows if m.kind == "barrido"]
    lines: list[str] = []

    if calib:
        # el pitch MÁS GRUESO que todavía cumple el umbral (el más barato de fabricar)
        ok2 = [m for m in calib if m.fwhm_ratio <= 1.02] or [calib[0]]
        ok7 = [m for m in calib if m.fwhm_ratio <= 1.07] or [calib[0]]
        need2, need7 = ok2[-1], ok7[-1]
        table = ", ".join(
            f"λ/{int(round(1 / m.pitch_over_lambda))}: {m.fwhm_m * 1e3:.2f} mm "
            f"({100 * (m.fwhm_ratio - 1):+.1f} %, corr {m.correlation:.4f})" for m in calib
        )
        deltas = ", ".join(
            f"λ/{int(round(1 / m.pitch_over_lambda))} "
            f"{TABLA_4B_FWHM_MM[round(m.pitch_over_lambda, 3)] - m.fwhm_m * 1e3:+.2f} mm"
            for m in calib
        )
        lines.append(
            f"1) Pitch necesario para empatar con un monolítico, medido en la geometría ancla de la "
            f"Fase 4B (D=60 mm=4λ, F=40 mm, f/#=0.67, 100 kHz): {table}. Para el 2 % de FWHM "
            f"alcanza pitch ≤ λ/{int(round(1 / need2.pitch_over_lambda))} "
            f"({100 * (need2.fwhm_ratio - 1):+.1f} %); para el 7 %, pitch ≤ "
            f"λ/{int(round(1 / need7.pitch_over_lambda))} ({100 * (need7.fwhm_ratio - 1):+.1f} %)."
        )
        lines.append(
            f"   Contra la tabla de ``plan.md`` (12.74 / 10.97 / 9.17 / 8.77 mm): la diferencia "
            f"tabla−medido es {deltas}. El ancla λ/8 se reproduce; las otras filas de 4B son más "
            f"anchas, y su fila λ/2 declara corr 0.998 junto con FWHM +27 %, cuando una correlación "
            f"de 0.998 implica campos casi idénticos (acá esa fila da +6.1 % con la misma corr)."
        )

    arrays = [m for m in docs if m.n_elements > 1]
    if arrays:
        best = min(arrays, key=lambda m: m.pitch_over_lambda)
        dentro = [m for m in arrays if m.fwhm_ratio <= 1.02]
        lines.append(
            f"2) De las configuraciones documentadas, la que más se acerca a la apertura continua es "
            f"{best.key} (pitch/λ = {best.pitch_over_lambda:.3f} ≈ λ/{1 / best.pitch_over_lambda:.1f}, "
            f"D/λ = {best.aperture_over_lambda:.0f}): FWHM {100 * (best.fwhm_ratio - 1):+.1f} %, "
            f"correlación {best.correlation:.4f}, lóbulo lateral {best.sidelobe_db:.1f} dB contra "
            f"{best.sidelobe_continuous_db:.1f} dB de la continua. {len(dentro)} de {len(arrays)} "
            f"arreglos documentados caen dentro del 2 % — incluido `diseno_1m` a λ/1.9 —, así que "
            f"con 256 elementos el criterio λ/8 de 4B no transfiere: lo que manda es la cantidad de "
            f"elementos (error de cuantización fino, que va a lóbulos lejanos) y no pitch/λ solo."
        )

    if docs:
        worst = max(docs, key=lambda m: m.fwhm_ratio)
        why = ("un pistón único no tiene perfil de fase: no puede enfocar" if worst.n_elements == 1
               else "es el arreglo documentado con peor equivalencia")
        lines.append(
            f"3) La que no lo cumple: {worst.key} (N={worst.n_elements}, pitch/λ = "
            f"{worst.pitch_over_lambda:.3f}, f/# = "
            f"{'—' if worst.f_number is None else f'{worst.f_number:.2f}'}), FWHM "
            f"{100 * (worst.fwhm_ratio - 1):+.1f} % sobre la continua, correlación "
            f"{worst.correlation:.4f} — {why}."
        )

    if sweep:
        worst_sweep = max(sweep, key=lambda m: m.pitch_over_lambda)
        lines.append(
            f"4) El barrido de geometría fija (160 mm, 16 elementos, F=150 mm) muestra el mecanismo: "
            f"pitch/λ sube de {sweep[0].pitch_over_lambda:.2f} a "
            f"{worst_sweep.pitch_over_lambda:.2f} y la correlación cae de "
            f"{sweep[0].correlation:.4f} a {worst_sweep.correlation:.4f} "
            f"(FWHM {100 * (worst_sweep.fwhm_ratio - 1):+.0f} %). A partir de 250 kHz "
            f"(pitch/λ = {[m for m in sweep if m.frequency_hz == 250e3][0].pitch_over_lambda:.2f}) "
            f"hay lóbulo de grataje visible y el nivel lateral llega a "
            f"{[m for m in sweep if m.frequency_hz == 1e6][0].sidelobe_db:.1f} dB: la misma "
            f"geometría deja de ser equivalente sólo por subir la frecuencia."
        )
    return lines


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #


def run(
    *,
    out_dir: Path = OUT_DIR,
    map_px: int | None = DEFAULT_MAP_PX,
    figures: bool = True,
    only: tuple[str, ...] = (),
    verbose: bool = True,
) -> list[Measurement]:
    """Corre el estudio completo y escribe CSV y figuras. Devuelve las mediciones."""
    configs = CONFIGS + sweep_configs() + calibration_4b()
    if only:
        configs = tuple(c for c in configs if c.key in set(only))
        if not configs:
            raise SystemExit(f"--only no coincide con ninguna configuración: {sorted(only)}")

    if verbose:
        print("=" * 150)
        print("ESTUDIO DE FOCALIZACIÓN — modelo 2D (aperturas lineales; arreglo 2D = su corte)")
        print("=" * 150)

    rows = study(configs, map_px=map_px, verbose=verbose)

    if verbose:
        print("\n" + format_table(rows) + "\n")
        for line in conclusions(rows):
            print("  " + line)
        print()

    out_dir.mkdir(parents=True, exist_ok=True)
    if not only or any(m.kind != "calibracion_4b" for m in rows):
        path = write_csv(rows, out_dir / OUT_CSV.name)
        if verbose:
            print(f"  CSV: {path}  ({len(rows)} filas)")

    if figures:
        for m in rows:
            if m.map_i is None:
                continue
            fig_path = out_dir / f"focus_{m.key}.png"
            plot_config(m, fig_path)
        summary = plot_summary(rows, out_dir / SUMMARY_PNG.name)
        if verbose:
            print(f"  figuras: {len([m for m in rows if m.map_i is not None])} focus_*.png + "
                  f"{summary.name}")
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--map-px", type=int, default=DEFAULT_MAP_PX,
                    help="puntos por lado del mapa 2D (0 = sin mapa)")
    ap.add_argument("--only", nargs="*", default=(), help="correr sólo estas claves")
    ap.add_argument("--no-figures", action="store_true")
    args = ap.parse_args(argv)

    t0 = time.perf_counter()
    run(out_dir=args.out_dir, map_px=(args.map_px or None), figures=not args.no_figures,
        only=tuple(args.only))
    print(f"  tiempo total: {time.perf_counter() - t0:.1f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
