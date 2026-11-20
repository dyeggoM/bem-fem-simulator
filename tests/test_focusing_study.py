"""Estudio de focalización (``scripts/focusing_study.py``): los cinco criterios del encargo.

1. **Regresión contra la Fase 4B**: con D = 60 mm, pitch = λ/8, F = 40 mm, 100 kHz y agua, el FWHM
   medido reproduce el del estudio 4B (8.77 mm) dentro de ±0.15 mm y el pico axial cae en la banda
   medida 25–33 mm (D13: el máximo no está en F).
2. **Lóbulo de grataje**: el ángulo medido coincide con ``grating_lobe_sine`` dentro de 0.5°.
   El caso literal (pitch = λ, steer = 0) tiene la predicción en ±90°, justo donde el factor de
   elemento de un pistón de ancho = pitch = λ tiene un **cero exacto**: no hay lóbulo medible. Eso se
   documenta con números y el criterio se verifica donde el lóbulo es interior (pitch = λ, steer
   −30° → +30°, amplitud igual al principal).
3. **Barrido de frecuencia con geometría fija** (160 mm, 16 elementos, F = 150 mm): ``pitch/λ``
   monótono creciente con f y correlación con la apertura continua monótona decreciente.
4. **Procedencia**: toda configuración del CSV tiene ``fuente`` no vacía y ``n_elementos`` > 0.
5. **Figuras**: se generan y pesan más de 10 kB.

Costos: los tests 1–4 corren en ~45 s (el barrido de frecuencia pesa ~28 s y los otros ~17 s; las
aperturas son de 16–32 elementos). El test 5 es ``slow`` porque corre el pipeline completo con mapas
reducidos (``--map-px 40``).
"""

from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

from ultrasound_sim.aperture import far_field_pattern
from ultrasound_sim.arrays import LinearArray, grating_lobe_sine
from ultrasound_sim.materials import WATER

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "focusing_study.py"

#: Ancla de la Fase 4B (``docs/plan.md``, tabla de referencia).
ANCHOR_FWHM_MM = 8.77
ANCHOR_FWHM_TOL_MM = 0.15
ANCHOR_PEAK_BAND_MM = (25.0, 33.0)

#: Valores de la tabla de 4B, para contrastar (no se exigen: la medición fina difiere, ver el
#: docstring de ``test_fase_4b_anchor``).
TABLA_4B_FWHM_MM = {1.0: 12.74, 0.5: 10.97, 0.25: 9.17, 0.125: 8.77}
TABLA_4B_CORR = {1.0: 0.904, 0.5: 0.998, 0.25: 0.9999, 0.125: 1.0}


def _load_script():
    """Carga ``scripts/focusing_study.py`` como módulo (``scripts/`` no es un paquete)."""
    name = "focusing_study"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, SCRIPT)
    if spec is None or spec.loader is None:
        raise ImportError(f"no se pudo cargar {SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module  # dataclass necesita el módulo registrado
    spec.loader.exec_module(module)
    return module


FS = _load_script()


def _cfg(key: str, f: float, n: int, width: float, focus: float | None):
    """Configuración de prueba con la procedencia de la Fase 4B."""
    return FS.Config(
        key=key,
        name=key,
        frequency_hz=f,
        n_elements=n,
        aperture_m=width,
        focus_m=focus,
        source="test: geometría de la Fase 4B (docs/plan.md)",
    )


# --------------------------------------------------------------------------- #
# (1) regresión contra la Fase 4B
# --------------------------------------------------------------------------- #


def test_fase_4b_anchor_fwhm_and_peak():
    """D = 60 mm, pitch = λ/8, F = 40 mm, 100 kHz: FWHM 8.77 ± 0.15 mm y pico en 25–33 mm.

    Medido acá: FWHM = 8.66 mm con el barrido axial fino (paso 0.05 mm) y pico en 29.6 mm. La
    diferencia de 0.11 mm con el ancla 8.77 no es física nueva: 8.77 corresponde a medir en el plano
    del pico que da la cuadratura del paquete (30.4–30.8 mm, D13) y no en el máximo fino (29.6 mm);
    a x = 30.4 mm esta misma apertura da 8.765 mm. Los dos valores están dentro de la banda pedida.

    La referencia continua de este estudio (mismo ancho, pitch λ/16) da 8.632 mm, o sea que el
    arreglo de λ/8 está a +0.3 % — consistente con la correlación ~1 que reporta 4B.
    """
    m = FS.measure(_cfg("ancla_4b", 100.0e3, 32, 0.060, 0.040), map_px=None)
    fwhm_mm = m.fwhm_m * 1e3
    peak_mm = m.peak_x_m * 1e3

    print("\n(1) ancla Fase 4B: D=60 mm, pitch=λ/8, F=40 mm, 100 kHz")
    print(f"    pico axial = {peak_mm:.3f} mm (banda 25–33 mm; corrimiento "
          f"{m.focal_shift_m * 1e3:+.2f} mm hacia la apertura)")
    print(f"    FWHM = {fwhm_mm:.3f} mm vs ancla {ANCHOR_FWHM_MM} ± {ANCHOR_FWHM_TOL_MM} mm "
          f"(Δ {fwhm_mm - ANCHOR_FWHM_MM:+.3f} mm)")
    print(f"    FWHM de la continua (λ/16) = {m.fwhm_continuous_m * 1e3:.3f} mm ⇒ "
          f"arreglo/continua = {m.fwhm_ratio:.4f} (+{100 * (m.fwhm_ratio - 1):.2f} %)")
    print(f"    |p(F)|²/máx = {m.p2_at_focus_rel:.4f}, correlación = {m.correlation:.6f}, "
          f"lóbulo lateral = {m.sidelobe_db:.2f} dB")

    assert abs(fwhm_mm - ANCHOR_FWHM_MM) <= ANCHOR_FWHM_TOL_MM, (
        f"FWHM {fwhm_mm:.3f} mm fuera de {ANCHOR_FWHM_MM} ± {ANCHOR_FWHM_TOL_MM} mm")
    assert ANCHOR_PEAK_BAND_MM[0] <= peak_mm <= ANCHOR_PEAK_BAND_MM[1], (
        f"pico axial en {peak_mm:.2f} mm, fuera de la banda medida 25–33 mm")
    assert m.focal_shift_m < 0.0, "el corrimiento focal tiene que ir hacia la apertura (D13)"
    assert m.fwhm_ratio < 1.02, "λ/8 tiene que empatar con la continua dentro del 2 %"


def test_tabla_4b_parcialmente_reproducida():
    """Qué filas de la tabla de 4B reproduce la medición fina y cuáles no (se reporta, no se fuerza).

    Medido con la caja de correlación de 4B (25–55 mm × ±30 mm) y la referencia λ/16:
    λ → 10.96 mm (+27.0 %, corr 0.956), λ/2 → 9.16 (+6.1 %, 0.998), λ/4 → 8.75 (+1.4 %, 0.9999),
    λ/8 → 8.66 (+0.3 %, 1.0000). La tabla de ``plan.md`` dice 12.74 / 10.97 / 9.17 / 8.77 y
    correlaciones 0.904 / 0.998 / 0.9999 / ~1.

    Coinciden el ancla λ/8 (0.11 mm) y las correlaciones de λ/2, λ/4 y λ/8. Las filas λ/2 y λ/4
    quedan 1.8 y 0.4 mm por debajo de lo medido acá, y la fila λ no se reproduce (10.96 vs 12.74;
    corr 0.956 vs 0.904). La correlación es casi independiente del muestreo de la referencia
    (λ/8…λ/64 dan 0.9562–0.9573), así que las diferencias no vienen de la referencia continua.
    Se deja escrito porque el criterio del encargo sólo fija la fila λ/8, y maquillar las otras
    sería reportar un número que no sale del código actual.
    """
    rows = {}
    for div in (1, 2, 4, 8):
        m = FS.measure(_cfg(f"calib_{div}", 100.0e3, int(round(0.060 / (0.015 / div))),
                            0.060, 0.040), map_px=None, corr_box=FS.BOX_4B)
        rows[1.0 / div] = m

    print("\n    pitch/λ   FWHM [mm]   tabla 4B   Δ [mm]   corr     corr 4B")
    for ratio, m in sorted(rows.items()):
        print(f"    {ratio:6.3f}   {m.fwhm_m * 1e3:8.3f}   {TABLA_4B_FWHM_MM[ratio]:8.2f}   "
              f"{m.fwhm_m * 1e3 - TABLA_4B_FWHM_MM[ratio]:+6.3f}   {m.correlation:.5f}  "
              f"{TABLA_4B_CORR[ratio]:.4f}")

    # lo que sí se exige: FWHM creciente con el pitch y correlación decreciente
    ordered = [rows[r] for r in sorted(rows)]
    assert all(b.fwhm_m > a.fwhm_m for a, b in zip(ordered, ordered[1:])), (
        "el FWHM tiene que crecer al engrosar el pitch")
    assert all(b.correlation <= a.correlation + 1e-9 for a, b in zip(ordered, ordered[1:])), (
        "la correlación con la continua no puede subir al engrosar el pitch")
    # el ancla se reproduce; las otras filas de 4B son más anchas (se deja dicho, no se relaja)
    anchor = rows[0.125]
    assert abs(anchor.fwhm_m * 1e3 - TABLA_4B_FWHM_MM[0.125]) < 0.2
    assert rows[0.5].fwhm_m * 1e3 < TABLA_4B_FWHM_MM[0.5] - 1.0, (
        "la fila λ/2 de la tabla 4B no se reproduce con la medición fina")


# --------------------------------------------------------------------------- #
# (2) lóbulo de grataje
# --------------------------------------------------------------------------- #


def test_grating_lobe_literal_case_is_degenerate():
    """``pitch = λ``, ``steer = 0``: la predicción es ±90° y ahí el factor de elemento es nulo.

    Con elementos contiguos de ancho ``w = pitch``, el patrón total es ``AF(θ)·sinc(w·sinθ/λ)``. En
    la dirección de grataje ``sinθ_g = λ/pitch`` el argumento del ``sinc`` vale exactamente 1 — un
    cero — así que **el lóbulo de grataje no existe como lóbulo**: lo que aparece fuera del principal
    es un lóbulo lateral normal. Medido (N=16, λ=15 mm): el máximo fuera del principal está en 5.1°
    con −13.2 dB, y a 89.5° el campo está > 40 dB abajo.

    Esto no es un fallo del modelo ni de ``grating_lobe_sine``: la fórmula da el ángulo correcto
    (90°); lo que se anula es la amplitud. El criterio de 0.5° se verifica en el caso interior
    (``test_grating_lobe_angle_matches_prediction``).
    """
    lam = WATER.wavelength(100.0e3)
    predictions = grating_lobe_sine(lam, lam)
    angles = [float(np.rad2deg(np.arcsin(s))) for _, s in predictions]
    print(f"\n(2a) pitch=λ={lam * 1e3:.1f} mm, steer=0: grating_lobe_sine → "
          f"{[(m, round(float(s), 4)) for m, s in predictions]} = {angles}°")

    n = 16
    theta = np.deg2rad(np.arange(-89.9, 89.901, 0.05))
    array = LinearArray(0.0, n, lam)
    field = far_field_pattern(theta, 100.0e3, WATER, array.steered(0.0, 100.0e3, WATER),
                             radius=20.0)
    intensity = np.abs(field) ** 2
    intensity /= intensity.max()
    outside = np.abs(np.rad2deg(theta)) > 5.0
    index = int(np.argmax(np.where(outside, intensity, -1.0)))
    at_grazing = float(intensity[np.argmin(np.abs(np.rad2deg(theta) - 89.0))])
    print(f"    máximo fuera del principal: θ={np.rad2deg(theta[index]):+.2f}°, "
          f"nivel={10 * np.log10(intensity[index]):+.2f} dB (lóbulo lateral normal)")
    print(f"    campo a 89° (dirección de grataje): {10 * np.log10(at_grazing):+.1f} dB "
          f"(cero del factor de elemento, sinc(1)=0)")

    assert all(abs(abs(a) - 90.0) < 1e-6 for a in angles), (
        f"la predicción tiene que ser ±90°: {angles}")
    assert 10 * np.log10(at_grazing) < -40.0, "en 90° el factor de elemento tiene que anular el campo"
    assert np.rad2deg(theta[index]) < 10.0, "el lóbulo más alto fuera del principal es un lateral"


def test_grating_lobe_angle_matches_prediction():
    """``pitch = λ``, steer −30°: el lóbulo de grataje medido coincide con la predicción en 0.5°.

    Con ``sinθ_g = sin(−30°) + λ/pitch = 0.5`` la predicción es +30.0°. Medido con N=32: +29.96°
    (error 0.04°). Con N=8 el mismo pico cae en 29.38° porque el factor de elemento (``sinc``) tira
    del máximo compuesto hacia broadside: el ancho del lóbulo del array factor en ``sinθ`` es ``2/N``
    y el ``sinc`` no es plano sobre él. Se usa N=32, que es el régimen donde la comparación
    ``AF`` ↔ campo es limpia.

    La amplitud relativa al principal es 0.0 dB, que es la identidad ``sinc(x+1) = −sinc(x)``: con
    ``pitch = λ`` y apuntando a −30°, los dos órdenes tienen el mismo factor de elemento.
    """
    lam = WATER.wavelength(100.0e3)
    steer = np.deg2rad(-30.0)
    predicted = grating_lobe_sine(lam, lam, steer)
    target = float(np.rad2deg(np.arcsin(predicted[0][1])))
    print(f"\n(2b) pitch=λ={lam * 1e3:.1f} mm, N=32, steer=−30°: predicción "
          f"sinθ_g={predicted[0][1]:.4f} ⇒ {target:+.2f}°")

    theta = np.deg2rad(np.arange(-85.0, 85.001, 0.01))
    array = LinearArray(0.0, 32, lam)
    field = far_field_pattern(theta, 100.0e3, WATER, array.steered(steer, 100.0e3, WATER),
                             radius=20.0)
    intensity = np.abs(field) ** 2
    intensity /= intensity.max()
    peaks = np.nonzero((intensity[1:-1] > intensity[:-2])
                       & (intensity[1:-1] >= intensity[2:]))[0] + 1
    main = peaks[int(np.argmin(np.abs(theta[peaks] - steer)))]
    outside = peaks[np.abs(np.rad2deg(theta[peaks]) - np.rad2deg(theta[main])) > 5.0]
    assert outside.size > 0, "no se encontró el lóbulo de grataje"
    top = outside[int(np.argmax(intensity[outside]))]
    measured = float(np.rad2deg(theta[top]))
    level = float(10 * np.log10(intensity[top]))
    print(f"    principal en {np.rad2deg(theta[main]):+.3f}°; lóbulo secundario en "
          f"{measured:+.3f}° a {level:+.2f} dB")
    print(f"    |medido − predicción| = {abs(measured - target):.3f}° (tolerancia 0.5°)")

    assert abs(measured - target) <= 0.5, (
        f"lóbulo medido en {measured:.3f}°, predicción {target:.3f}°")
    assert level > -1.0, "con pitch=λ los dos órdenes tienen la misma amplitud"


# --------------------------------------------------------------------------- #
# (3) barrido de frecuencia con geometría fija
# --------------------------------------------------------------------------- #


def test_frequency_sweep_monotone_pitch_and_decreasing_correlation():
    """160 mm / 16 elementos / F = 150 mm en 40, 100, 250, 500 y 1000 kHz.

    ``pitch/λ = 0.267 … 6.667`` (estrictamente creciente, porque ``pitch`` es fijo y ``λ`` baja) y
    la correlación con la apertura continua cae 0.9999 → 0.1844. Es el mecanismo de la equivalencia:
    la misma geometría deja de aproximar una apertura continua sólo por subir la frecuencia, y a
    partir de 250 kHz el lóbulo de grataje ya es visible en campo lejano.
    """
    ms = [FS.measure(c, map_px=None) for c in FS.sweep_configs()]

    print("\n(3) barrido 160 mm / 16 elementos / F=150 mm")
    print("      f [kHz]   pitch/λ   D/λ    FWHM [mm]   FWHM/cont   corr      grat. [°]")
    for m in ms:
        grat = "—" if m.grating_lobe_deg is None else f"{m.grating_lobe_deg:.2f}"
        print(f"      {m.frequency_hz / 1e3:7.0f}   {m.pitch_over_lambda:7.3f}   "
              f"{m.aperture_over_lambda:5.1f}   {m.fwhm_m * 1e3:8.2f}   {m.fwhm_ratio:8.3f}   "
              f"{m.correlation:.4f}   {grat:>8s}")

    ratios = [m.pitch_over_lambda for m in ms]
    corr = [m.correlation for m in ms]
    assert all(b > a for a, b in zip(ratios, ratios[1:])), f"pitch/λ no monótono: {ratios}"
    assert all(b < a for a, b in zip(corr, corr[1:])), f"correlación no decreciente: {corr}"
    assert corr[0] > 0.999 and corr[-1] < 0.30, (
        f"la caída de correlación tiene que ser grande: {corr[0]:.4f} → {corr[-1]:.4f}")
    # el lóbulo de grataje sólo existe cuando λ/pitch <= 1, o sea f >= 150 kHz
    assert all(m.grating_lobe_deg is None for m in ms if m.frequency_hz < 150.0e3)
    assert all(m.grating_lobe_deg is not None for m in ms if m.frequency_hz >= 250.0e3)
    # y la ganancia deja de seguir a 20log10(N) cuando el foco se vuelve no paraxial (D13)
    assert all(abs(m.gain_db - m.gain_ideal_db) > 0.5 for m in ms), (
        "la ganancia de N no debería ser 20log10(N) en f/# ≈ 0.94")


# --------------------------------------------------------------------------- #
# (4) procedencia en el CSV
# --------------------------------------------------------------------------- #

REQUIRED_COLUMNS = {
    "nombre", "tipo", "fuente", "n_elementos", "pitch_mm", "apertura_mm", "pitch_over_lambda",
    "fwhm_mm", "fwhm_continua_mm", "fwhm_ratio", "sidelobe_db", "correlacion", "mapa_px",
}


def test_every_configuration_declares_provenance():
    """Toda configuración del estudio (documentadas, barrido y calibración) tiene fuente y N > 0."""
    configs = FS.CONFIGS + FS.sweep_configs() + FS.calibration_4b()
    print(f"\n(4) {len(configs)} configuraciones declaradas")
    for cfg in configs:
        print(f"    {cfg.key:22s} N={cfg.n_elements:4d} D={cfg.aperture_m * 1e3:6.1f} mm  "
              f"f={cfg.frequency_hz / 1e3:7.0f} kHz  fuente={cfg.source[:60]}…")
        assert cfg.source.strip(), f"{cfg.key}: fuente vacía"
        assert cfg.n_elements > 0, f"{cfg.key}: n_elements = {cfg.n_elements}"
        assert cfg.pitch_m > 0.0
    # las seis configuraciones documentadas de D20 (bloque 6) están todas
    assert [c.key for c in FS.CONFIGS] == [
        "aire_40k_linea", "langevin_40k_plano", "langevin_40k_foco100",
        "transcraneal_250k", "transcraneal_500k", "diseno_1m",
    ]


def test_csv_rows_have_provenance(tmp_path):
    """El CSV generado: columnas del encargo, ``fuente`` no vacía y ``n_elementos`` > 0 por fila.

    Se regenera el CSV para dos configuraciones en ``tmp_path`` (no depende de que ``out/`` exista)
    y si además está el CSV completo del estudio, se valida también ese.
    """
    rows = FS.run(
        out_dir=tmp_path,
        map_px=None,
        figures=False,
        only=("aire_40k_linea", "calib4b_lambda_8"),
        verbose=False,
    )
    path = tmp_path / "focusing_study.csv"
    with path.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        assert REQUIRED_COLUMNS <= set(reader.fieldnames or ())
        written = list(reader)
    assert len(written) == len(rows) == 2
    for row in written:
        assert row["fuente"].strip(), f"{row['nombre']}: fuente vacía"
        assert int(row["n_elementos"]) > 0
        assert float(row["pitch_over_lambda"]) > 0.0
        assert float(row["fwhm_mm"]) > 0.0
        assert row["tipo"] in {"documentada", "barrido", "calibracion_4b"}
        assert "2D" in row["notas"]

    full = REPO / "out" / "focusing_study.csv"
    if full.is_file():
        with full.open(encoding="utf-8") as fh:
            allrows = list(csv.DictReader(fh))
        assert len(allrows) >= 15
        assert all(r["fuente"].strip() and int(r["n_elementos"]) > 0 for r in allrows)
        print(f"\n(4) CSV completo: {len(allrows)} filas, todas con fuente y n_elementos")
    else:
        print("\n(4) CSV completo ausente (correr scripts/focusing_study.py); se validó el parcial")


# --------------------------------------------------------------------------- #
# (5) figuras
# --------------------------------------------------------------------------- #


@pytest.mark.slow
def test_figures_are_generated_and_not_empty(tmp_path):
    """Un PNG por configuración + el resumen, todos > 10 kB (con mapas reducidos a 40×40)."""
    FS.run(
        out_dir=tmp_path,
        map_px=40,
        figures=True,
        only=("aire_40k_linea", "transcraneal_250k", "calib4b_lambda_8"),
        verbose=False,
    )
    figures = sorted(tmp_path.glob("focus_*.png"))
    summary = tmp_path / "focusing_summary.png"
    print(f"\n(5) {len(figures)} figuras por configuración + {summary.name}")
    for fig in figures + [summary]:
        size = fig.stat().st_size
        print(f"    {fig.name:34s} {size / 1e3:8.1f} kB")
        assert fig.is_file() and size > 10_000, f"{fig.name} vacía o demasiado chica ({size} B)"
    assert len(figures) == 3
