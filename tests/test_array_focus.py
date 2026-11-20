"""Fase 4B — arreglos lineales con delays: lóbulo lateral, foco, equivalencia con apertura continua.

Criterios (``docs/plan.md``, Fase 4B), y lo que se mide de verdad:

(a) lóbulo lateral de un arreglo uniforme ``N=8``, ``pitch=λ/2``, dentro de 2 dB de la fórmula de
    array factor ``|sin(Nψ/2)/(N sin(ψ/2))|``;
(b) foco de un arreglo enfocado en ``F=(40 mm, 0)``: pico dentro de 1 mm de ``F`` y ganancia
    dentro de 1 dB de ``N`` en presión;
(c) equivalencia arreglo vs apertura continua del mismo ancho total (``N·pitch``): correlación de
    campo complejo > 0.95 y FWHM dentro del 5 %, para ``pitch = λ/4`` y ``λ/2``;
(d) el estudio que importa: ``pitch = λ`` (lóbulo de grataje) — ángulo del lóbulo y caída de
    correlación.

Los tres primeros salen como salen; el punto (b) **no se cumple con los parámetros del plan** y el
test lo dice con el número medido y el diagnóstico (ver ``test_focus_peak_at_geometric_focus``).
"""

import numpy as np
import pytest

from ultrasound_sim.aperture import Aperture, far_field_pattern, pressure
from ultrasound_sim.arrays import LinearArray, array_factor, grating_lobe_sine
from ultrasound_sim.focus import (field_correlation, intensity_gain_db, lateral_fwhm,
                                  peak_location, sidelobe_level_db)
from ultrasound_sim.materials import WATER

FREQUENCY = 100.0e3
LAMBDA = WATER.wavelength(FREQUENCY)          # 15 mm
N_ELEMENTS = 8
PITCH = 0.5 * LAMBDA                           # 7.5 mm
FOCUS = np.array([40.0e-3, 0.0])
FAR_RADIUS = 5.0


# --------------------------------------------------------------------------- #
# métricas: se validan contra definiciones cerradas antes de usarlas para juzgar el campo
# --------------------------------------------------------------------------- #


def test_metric_definitions():
    """``lateral_fwhm``, ``sidelobe_level_db``, ``intensity_gain_db`` y ``field_correlation``.

    Si una métrica está mal, todas las comparaciones de arriba mienten sin fallar; por eso se
    contrastan con valores analíticos: FWHM de una gaussiana (``2√(ln2)·σ``), lóbulo lateral de un
    patrón sintético, ``6.02 dB`` al duplicar la presión y correlación invariante a fase global.
    """
    sigma = 2.0e-3
    y = np.arange(-12.0e-3, 12.0e-3, 1.0e-4)
    gaussian = np.exp(-(y / sigma) ** 2)
    fwhm = lateral_fwhm(y, gaussian)
    analytic = 2.0 * sigma * np.sqrt(np.log(2.0))
    print(f"\n  gaussiana σ={sigma * 1e3:.2f} mm: FWHM medido = {fwhm * 1e3:.6f} mm, "
          f"2√(ln2)σ = {analytic * 1e3:.6f} mm, error = {abs(fwhm / analytic - 1):.2e}")

    single = sidelobe_level_db(y, gaussian)
    synthetic = gaussian + 0.02 * ((y > 6.0e-3) & (y < 7.0e-3))
    sidelobe = sidelobe_level_db(y, synthetic)
    print(f"    lóbulo lateral: gaussiana sola = {single} dB;  con un lóbulo de 2 % = "
          f"{sidelobe:.3f} dB (10log10 0.02 = {10 * np.log10(0.02):.3f} dB)")

    field = np.array([1.0 + 1.0j, 2.0 - 0.5j, -0.3j])
    gain = intensity_gain_db(field, 2.0 * field)
    rotated = field * np.exp(1.0j * 0.7)
    correlation = field_correlation(field, rotated)
    anticorrelated = field_correlation(field, 1.0j * field.conj())
    print(f"    duplicar presión: ganancia = {gain:.4f} dB (10log10 4 = {10 * np.log10(4):.4f} dB);  "
          f"correlación con fase global e^(i0.7) = {correlation:.12f};  "
          f"con el conjugado = {anticorrelated:.6f}")

    grid = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 2.0]])
    peak = peak_location(grid, np.array([1.0, 5.0, 3.0]))
    print(f"    peak_location = {peak}")

    assert abs(fwhm / analytic - 1.0) < 1e-3
    assert single == -np.inf
    assert abs(sidelobe - 10.0 * np.log10(0.02)) < 0.6
    assert abs(gain - 10.0 * np.log10(4.0)) < 1e-12
    assert abs(correlation - 1.0) < 1e-12
    assert anticorrelated < 0.9
    assert np.array_equal(peak, grid[1])


def local_maxima(y: np.ndarray) -> np.ndarray:
    """Índices de máximos locales estrictos."""
    return np.nonzero((y[1:-1] > y[:-2]) & (y[1:-1] >= y[2:]))[0] + 1


def focused_region_field(array: LinearArray, n_x: int = 31, n_y: int = 41) -> tuple:
    """Campo de un arreglo enfocado en ``F`` sobre la región blanco (caja alrededor de ``F``)."""
    xs = np.linspace(25.0e-3, 55.0e-3, n_x)
    ys = np.linspace(-30.0e-3, 30.0e-3, n_y)
    grid_x, grid_y = np.meshgrid(xs, ys, indexing="ij")
    points = np.column_stack([grid_x.ravel(), grid_y.ravel()])
    field = pressure(points, FREQUENCY, WATER, array.focused(FOCUS, FREQUENCY, WATER))
    return points, field


# --------------------------------------------------------------------------- #
# (a) lóbulo lateral del arreglo uniforme vs array factor
# --------------------------------------------------------------------------- #


def test_uniform_array_sidelobe_matches_array_factor():
    """``N=8``, ``pitch=λ/2``: primer lóbulo lateral medido vs ``|sin(Nψ/2)/(N sin(ψ/2))|``."""
    array = LinearArray(0.0, N_ELEMENTS, PITCH)
    theta = np.deg2rad(np.arange(-80.0, 80.001, 0.01))
    field = far_field_pattern(theta, FREQUENCY, WATER, array.apertures(), radius=FAR_RADIUS)
    intensity = np.abs(field) ** 2
    intensity /= intensity.max()
    factor = array_factor(theta, N_ELEMENTS, PITCH, FREQUENCY, WATER)

    measured_db = sidelobe_level_db(np.rad2deg(theta), intensity)
    factor_db = sidelobe_level_db(np.rad2deg(theta), factor ** 2)
    theory_db = 20.0 * np.log10(1.0 / (N_ELEMENTS * np.sin(1.5 * np.pi / N_ELEMENTS)))

    print(f"\n  N={N_ELEMENTS}, pitch=λ/2={PITCH * 1e3:.2f} mm, ancho total="
          f"{array.width * 1e3:.1f} mm ({array.width / LAMBDA:.1f}λ), f={FREQUENCY / 1e3:.0f} kHz, "
          f"r={FAR_RADIUS:g} m")
    print(f"    lóbulo lateral: medido = {measured_db:.3f} dB   array factor (muestreado) = "
          f"{factor_db:.3f} dB   fórmula 20log10(1/(N sin(3π/2N))) = {theory_db:.3f} dB   "
          f"|medido - fórmula| = {abs(measured_db - theory_db):.3f} dB")

    peaks_measured = local_maxima(intensity)
    peaks_factor = local_maxima(factor)
    order = np.argsort(theta[peaks_measured])
    print("    lóbulos (θ, medido, AF, AF·elemento):")
    for index in peaks_measured[order]:
        element = abs(np.sinc(PITCH * np.sin(theta[index]) / LAMBDA))
        print(f"      θ={np.rad2deg(theta[index]):+8.3f}°  medido={10 * np.log10(intensity[index]):+8.3f} dB"
              f"  AF={10 * np.log10(factor[index] ** 2 + 1e-300):+8.3f} dB"
              f"  AF·elem={10 * np.log10(factor[index] ** 2 * element ** 2 + 1e-300):+8.3f} dB")

    # el campo es exactamente AF × factor de elemento (aperturas contiguas e idénticas en el lejano)
    element_taper = np.abs(np.sinc(PITCH * np.sin(theta) / LAMBDA))
    product = factor * element_taper
    product /= product.max()
    deviation = np.abs(np.abs(field) / np.abs(field).max() - product)
    print(f"    identidad AF·elemento: desvío máximo = {deviation.max():.2e} "
          f"(sólo se aparta donde el patrón cae > 40 dB)")
    print("    picos del AF (teoría, °): "
          f"{[round(float(np.rad2deg(t)), 2) for t in theta[peaks_factor[:6]]]}")
    assert abs(measured_db - theory_db) < 2.0
    assert deviation.max() < 0.02


def axial_profile(array: LinearArray, frequency: float = FREQUENCY,
                  x_min: float = 5.0e-3, x_max: float = 70.0e-3,
                  step: float = 0.05e-3):
    """Perfil axial ``|p|²`` del arreglo enfocado en ``F``, sobre el eje ``y = 0``."""
    xs = np.arange(x_min, x_max, step)
    points = np.column_stack([xs, np.zeros_like(xs)])
    field = pressure(points, frequency, WATER, array.focused(FOCUS, frequency, WATER))
    return xs, field


# --------------------------------------------------------------------------- #
# (b) foco: posición y ganancia
# --------------------------------------------------------------------------- #


def test_focus_peak_position_and_flatness():
    """Posición del máximo axial y cuánto vale ``|p(F)|²`` frente a ese máximo.

    Se reporta el número tal cual: el máximo axial **no** cae en ``F`` (ver el test siguiente), y a
    100 kHz el foco es una meseta de ~1 dB a lo largo de 20 mm, no un punto.
    """
    array = LinearArray(0.0, N_ELEMENTS, PITCH)
    xs, field = axial_profile(array)
    intensity = np.abs(field) ** 2
    peak_x = float(xs[int(np.argmax(intensity))])
    at_focus = float(intensity[int(np.argmin(np.abs(xs - FOCUS[0])))])
    print(f"\n  N={N_ELEMENTS}, pitch=λ/2, ancho={array.width * 1e3:.1f} mm "
          f"({array.width / LAMBDA:.1f}λ, f/# = {FOCUS[0] / array.width:.2f}), "
          f"F = {FOCUS[0] * 1e3:.0f} mm, f = {FREQUENCY / 1e3:.0f} kHz")
    print(f"    pico axial = {peak_x * 1e3:.3f} mm   (F = {FOCUS[0] * 1e3:.1f} mm)   "
          f"desviación = {(peak_x - FOCUS[0]) * 1e3:+.3f} mm")
    print(f"    |p(F)|²/|p|²max = {at_focus / intensity.max():.4f} "
          f"({10 * np.log10(at_focus / intensity.max()):+.2f} dB)")
    print(f"    nivel axial a ±5 mm de F: "
          f"{10 * np.log10(intensity[np.argmin(np.abs(xs - FOCUS[0] - 5e-3))] / intensity.max()):+.2f} dB / "
          f"{10 * np.log10(intensity[np.argmin(np.abs(xs - FOCUS[0] + 5e-3))] / intensity.max()):+.2f} dB")
    location = peak_location(np.column_stack([xs, np.zeros_like(xs)]), field)
    print(f"    peak_location() = ({location[0] * 1e3:.3f}, {location[1] * 1e3:.3f}) mm")
    assert location[1] == 0.0
    # la meseta: a 100 kHz F no está a más de 1 dB del máximo (y por eso "el pico" no la distingue)
    assert 10.0 * np.log10(at_focus / intensity.max()) > -1.0


def test_focal_shift_is_measured_and_points_toward_the_aperture():
    """El foco geométrico NO es el máximo axial a 100 kHz: se mide y se fija esa banda.

    Con D=60 mm (4λ), F=40 mm, f/# = 0.67 y semiángulo 37°, el máximo axial cae en ~29.5 mm, o sea
    10.5 mm *hacia la apertura*. Tres vías independientes coinciden: integral directa por trapecio
    20 001 puntos → 29.5 mm; cuadratura del paquete → 30.8 mm; la misma apertura física a 500 kHz
    enfoca en 39.4 mm (caso paraxial, ``test_focus_peak_in_focusing_regime``).

    Por eso el criterio "``|pico − F| < 1 mm`` a 100 kHz" no es físico y se reemplaza por: banda
    medida del corrimiento + dirección correcta + el requisito de que en F la intensidad esté a
    menos de 1 dB del máximo (``test_focus_is_a_plateau``), que es lo que hace útil el foco.
    """
    array = LinearArray(0.0, N_ELEMENTS, PITCH)
    xs, field = axial_profile(array)
    intensity = np.abs(field) ** 2
    peak_x = float(xs[int(np.argmax(intensity))])
    print(f"\n  pico axial = {peak_x * 1e3:.2f} mm  (F = {FOCUS[0] * 1e3:.0f} mm, "
          f"corrimiento {(peak_x - FOCUS[0]) * 1e3:+.2f} mm)")
    assert 25e-3 < peak_x < 33e-3, (
        f"pico axial en {peak_x * 1e3:.2f} mm, fuera de la banda medida 25–33 mm")
    assert peak_x < FOCUS[0], "el corrimiento focal va hacia la apertura, no al revés"


def test_focus_peak_in_focusing_regime():
    """Mismo criterio de 1 mm, con la misma apertura física (60 mm) en el régimen de focalización.

    ``500 kHz`` ⇒ ``D = 20λ``, ``D²/(λF) = 30``: el foco geométrico sí es el máximo axial, y el
    FWHM lateral coincide con ``λ·F/D`` (difracción). Sirve de control: el modelo de apertura
    enfoca donde debe cuando la geometría lo permite.
    """
    frequency = 500.0e3
    wavelength = WATER.wavelength(frequency)
    width = 60.0e-3
    n_elements = 64
    array = LinearArray(0.0, n_elements, width / n_elements)
    apertures = array.focused(FOCUS, frequency, WATER)

    xs = np.arange(5.0e-3, 90.0e-3, 0.05e-3)
    axial = np.abs(pressure(np.column_stack([xs, np.zeros_like(xs)]),
                            frequency, WATER, apertures)) ** 2
    peak_x = float(xs[int(np.argmax(axial))])

    ys = np.arange(-15.0e-3, 15.0e-3, 0.02e-3)
    lateral = np.abs(pressure(np.column_stack([np.full_like(ys, peak_x), ys]),
                              frequency, WATER, apertures)) ** 2
    fwhm = lateral_fwhm(ys, lateral)
    diffraction = wavelength * FOCUS[0] / width
    print(f"\n  f={frequency / 1e3:.0f} kHz, {n_elements} elementos de {array.pitch * 1e3:.3f} mm "
          f"(λ/{wavelength / array.pitch:.1f}), ancho {width * 1e3:.0f} mm = {width / wavelength:.0f}λ, "
          f"D²/(λF) = {width ** 2 / (wavelength * FOCUS[0]):.1f}")
    print(f"    pico axial = {peak_x * 1e3:.3f} mm (F = {FOCUS[0] * 1e3:.1f} mm), desviación = "
          f"{abs(peak_x - FOCUS[0]) * 1e3:.3f} mm;  FWHM lateral = {fwhm * 1e3:.3f} mm vs "
          f"λF/D = {diffraction * 1e3:.3f} mm")
    assert abs(peak_x - FOCUS[0]) < 1.0e-3
    assert abs(fwhm / diffraction - 1.0) < 0.15


def test_focus_gain_and_sidelobe():
    """Ganancia de presión en ``F`` respecto de un elemento, y reparto de la energía en el foco.

    ``|p_arreglo(F)| = Σ_j |p_j(F)|`` porque en ``F`` todas las contribuciones están en fase (es la
    definición de los delays). La ganancia depende de con qué elemento se compare:

    * contra la **media** de los elementos: exactamente ``N`` (``20·log10 N``);
    * contra el elemento **central**: ``17.28 dB`` vs ``18.06 dB`` ⇒ 0.78 dB;
    * contra el elemento del **borde** (más lejos de ``F``): ``19.10 dB``.
    """
    array = LinearArray(0.0, N_ELEMENTS, PITCH)
    apertures = array.focused(FOCUS, FREQUENCY, WATER)
    array_pressure = abs(pressure(FOCUS, FREQUENCY, WATER, apertures)[0])

    per_element = []
    for aperture in apertures:
        per_element.append(abs(pressure(FOCUS, FREQUENCY, WATER, [aperture])[0]))
    per_element = np.array(per_element)
    gains = 20.0 * np.log10(array_pressure / per_element)
    central = N_ELEMENTS // 2
    ideal = 20.0 * np.log10(N_ELEMENTS)

    print(f"\n  |p_arreglo(F)| = {array_pressure:.6e} Pa; 20log10(N) = {ideal:.3f} dB")
    for j, (gain, value) in enumerate(zip(gains, per_element)):
        print(f"    elemento {j} (y={apertures[j].center_y * 1e3:+7.2f} mm): |p_j(F)|={value:.4e} Pa  "
              f"ganancia={gain:6.3f} dB")
    mean_gain = 20.0 * np.log10(array_pressure / per_element.mean())
    print(f"    ganancia vs media de elementos = {mean_gain:.4f} dB "
          f"(20log10 N = {ideal:.4f} dB: identidad de la suma coherente)")
    print(f"    ganancia vs elemento central  = {gains[central]:.3f} dB  "
          f"(dentro de 1 dB de 20log10 N: {abs(gains[central] - ideal):.3f} dB)")
    print(f"    ganancia vs elemento de borde  = {gains[0]:.3f} dB  (dif {abs(gains[0] - ideal):.3f} dB)")
    print(f"    referencia 10log10(N) = {10 * np.log10(N_ELEMENTS):.3f} dB: es la ganancia de "
          f"INTENSIDAD con potencia total igual, no la de N elementos de igual amplitud")

    # la versión "potencia total igual": amplitud 1/√N por elemento ⇒ 10log10(N) de intensidad
    phases = np.exp(1.0j * np.angle([a.amplitude for a in apertures]))
    equal_power = array.apertures(phases / np.sqrt(N_ELEMENTS))
    ref_field = pressure(FOCUS, FREQUENCY, WATER, [Aperture(apertures[central].center_y, PITCH)])
    test_field = pressure(FOCUS, FREQUENCY, WATER, equal_power)
    gain_db = intensity_gain_db(ref_field, test_field)
    print(f"    amplitud 1/√N (potencia total igual): ganancia de intensidad = {gain_db:.3f} dB "
          f"vs 10log10 N = {10 * np.log10(N_ELEMENTS):.3f} dB "
          f"(dif {abs(gain_db - 10 * np.log10(N_ELEMENTS)):.3f} dB)")

    ys = np.arange(-25.0e-3, 25.0e-3, 0.02e-3)
    lateral = np.abs(pressure(np.column_stack([np.full_like(ys, FOCUS[0]), ys]),
                              FREQUENCY, WATER, apertures)) ** 2
    fwhm = lateral_fwhm(ys, lateral)
    sidelobe = sidelobe_level_db(ys, lateral)
    print(f"    FWHM lateral en x=F = {fwhm * 1e3:.3f} mm (λF/D = "
          f"{LAMBDA * FOCUS[0] / array.width * 1e3:.3f} mm);  lóbulo lateral = {sidelobe:.2f} dB")

    assert abs(gains[central] - ideal) < 1.0
    assert abs(gain_db - 10.0 * np.log10(N_ELEMENTS)) < 1.0
    # identidad de la suma coherente: |p_arr| = Σ|p_j| salvo la fase residual de cada apertura
    # extendida (medido 2e-5 relativo, o sea 1.6e-4 dB)
    assert abs(mean_gain - ideal) < 1e-3


# --------------------------------------------------------------------------- #
# (c) equivalencia arreglo ↔ apertura continua del mismo ancho
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("pitch_in_lambda", [0.25, 0.5])
def test_array_equals_continuous_aperture(pitch_in_lambda):
    """``N`` aperturas contiguas de ancho ``pitch`` vs una apertura de ancho ``N·pitch``.

    Con elementos contiguos y la misma amplitud, la distribución de velocidad es **idéntica** a la
    de la apertura continua (no hay kerf ni apodización): la correlación tiene que dar 1 y el FWHM
    tiene que coincidir. El test es entonces una verificación de la aritmética de elementos (centros,
    bordes, suma de aperturas) y de que la cuadratura no introduce costuras entre elementos.

    Se agrega, como comparación con contenido, el arreglo **enfocado** contra la misma apertura
    muestreada finamente (32 elementos): ahí sí hay una diferencia física, la fase cuantizada por
    elemento.
    """
    pitch = pitch_in_lambda * LAMBDA
    array = LinearArray(0.0, N_ELEMENTS, pitch)
    continuous = Aperture(0.0, array.width)

    ys = np.arange(-60.0e-3, 60.0e-3, 0.2e-3)
    line = np.column_stack([np.full_like(ys, FOCUS[0]), ys])
    array_field = pressure(line, FREQUENCY, WATER, array.apertures())
    continuous_field = pressure(line, FREQUENCY, WATER, [continuous])
    correlation = field_correlation(array_field, continuous_field)
    fwhm_array = lateral_fwhm(ys, np.abs(array_field) ** 2)
    fwhm_continuous = lateral_fwhm(ys, np.abs(continuous_field) ** 2)
    ratio = fwhm_array / fwhm_continuous

    fine = LinearArray(0.0, 32, array.width / 32)
    points, focused_coarse = focused_region_field(array)
    _, focused_fine = focused_region_field(fine)
    focused_correlation = field_correlation(focused_coarse, focused_fine)

    print(f"\n  pitch=λ/{1 / pitch_in_lambda:g} = {pitch * 1e3:.3f} mm, ancho total = "
          f"{array.width * 1e3:.1f} mm ({array.width / LAMBDA:.1f}λ), N={N_ELEMENTS}")
    print(f"    uniforme (x=40 mm, |y|<=60 mm): correlación = {correlation:.12f}   "
          f"FWHM arreglo = {fwhm_array * 1e3:.6f} mm, continuo = {fwhm_continuous * 1e3:.6f} mm, "
          f"cociente = {ratio:.9f}")
    print(f"    enfocado (caja 25..55 mm × ±30 mm) vs 32 elementos finos: correlación = "
          f"{focused_correlation:.6f}")
    assert correlation > 0.95
    assert abs(ratio - 1.0) < 0.05
    assert focused_correlation > 0.95


# --------------------------------------------------------------------------- #
# (d) pitch = λ: lóbulo de grataje
# --------------------------------------------------------------------------- #


def test_grating_lobe_at_pitch_lambda():
    """``pitch = λ``: dónde aparece el lóbulo de grataje y cuánto degrada el campo.

    Con ``d = λ`` el criterio ``sinθ_g = sinθ₀ + m·λ/d`` sólo deja el orden ``m = ∓1`` dentro del
    espacio visible si se apunta hacia atrás. Se apunta a ``θ₀ = −30°``: el lóbulo principal queda
    en −29.4° y el de grataje en +29.4°, **con la misma amplitud** (0.0 dB de diferencia) porque con
    ``d = λ`` los dos órdenes son igualmente fuertes. Con ``pitch = λ/2`` no hay ningún orden visible
    y el lóbulo más alto fuera del principal queda en −12.0 dB (a −8.0°), o sea un lateral normal.

    Correlación del campo enfocado en la región blanco (contra la misma apertura muestreada fina):
    ``pitch = λ/2`` → 0.998, ``pitch = λ`` → 0.90.
    """
    theta = np.deg2rad(np.arange(-85.0, 85.001, 0.01))
    steer = np.deg2rad(-30.0)
    predictions = grating_lobe_sine(LAMBDA, LAMBDA, steer)
    print(f"\n  pitch = λ = {LAMBDA * 1e3:.1f} mm, N={N_ELEMENTS}, apuntando a −30°")
    print(f"    predicción sinθ_g = sinθ₀ + m·λ/d: "
          f"{[(m, round(float(s), 4), round(float(np.rad2deg(np.arcsin(s))), 2)) for m, s in predictions]}")

    results = {}
    for label, pitch in (("λ", LAMBDA), ("λ/2", 0.5 * LAMBDA)):
        array = LinearArray(0.0, N_ELEMENTS, pitch)
        field = far_field_pattern(theta, FREQUENCY, WATER, array.steered(steer, FREQUENCY, WATER),
                                  radius=20.0)
        intensity = np.abs(field) ** 2
        intensity /= intensity.max()
        peaks = local_maxima(intensity)
        # lóbulo principal: el pico más cercano a la dirección apuntada
        main = peaks[int(np.argmin(np.abs(theta[peaks] - steer)))]
        outside = peaks[np.abs(theta[peaks] - steer) > np.deg2rad(5.0)]
        top = outside[int(np.argmax(intensity[outside]))] if outside.size else None
        results[label] = (float(np.rad2deg(theta[main])), float(np.rad2deg(theta[top]))
                          if top is not None else None,
                          float(10.0 * np.log10(intensity[top])) if top is not None else None)
        print(f"    pitch={label:4s}: lóbulo principal en {np.rad2deg(theta[main]):+7.3f}° "
              f"({10 * np.log10(intensity[main]):+.2f} dB); más alto fuera de ±5°: "
              f"θ={np.rad2deg(theta[top]):+7.3f}° (sinθ={np.sin(theta[top]):+.4f}), "
              f"nivel={10 * np.log10(intensity[top]):+.2f} dB relativo al principal")
        print(f"      AF(θ) en esa dirección = {array_factor(theta[top], N_ELEMENTS, pitch, FREQUENCY, WATER, steer):.3f}"
              f";  " + ("orden de grataje predicho" if grating_lobe_sine(pitch, LAMBDA, steer)
                        else "ningún orden de grataje visible"))

    main_angle, grating_angle, grating_level = results["λ"]
    assert abs(main_angle + 30.0) < 1.0, f"el haz no apunta a −30°: {main_angle:.2f}°"
    assert grating_angle is not None
    assert abs(grating_angle - 30.0) < 1.5, (
        f"lóbulo de grataje medido en {grating_angle:.2f}°, esperado sinθ_g = 0.5 ⇒ 30°")
    assert grating_level > -1.0, "con d = λ el orden m=±1 tiene la amplitud del principal"
    assert abs(results["λ/2"][0] + 30.0) < 1.0
    assert results["λ/2"][2] < -10.0, "con d = λ/2 no debe haber ningún lóbulo de grataje"

    # caída de correlación en la región blanco (campo enfocado)
    print("    correlación del campo enfocado en la región blanco (vs 32 elementos del mismo ancho):")
    for pitch_in_lambda in (0.25, 0.5, 1.0):
        array = LinearArray(0.0, N_ELEMENTS, pitch_in_lambda * LAMBDA)
        fine = LinearArray(0.0, 32, array.width / 32)
        points, coarse = focused_region_field(array)
        _, reference = focused_region_field(fine)
        print(f"      pitch=λ/{1 / pitch_in_lambda:g}: correlación = "
              f"{field_correlation(coarse, reference):.6f}   "
              f"(pico en x = {peak_location(points, coarse)[0] * 1e3:.2f} mm)")
    array_lambda = LinearArray(0.0, N_ELEMENTS, LAMBDA)
    fine_lambda = LinearArray(0.0, 32, array_lambda.width / 32)
    _, coarse_lambda = focused_region_field(array_lambda)
    _, reference_lambda = focused_region_field(fine_lambda)
    array_half = LinearArray(0.0, N_ELEMENTS, 0.5 * LAMBDA)
    fine_half = LinearArray(0.0, 32, array_half.width / 32)
    _, coarse_half = focused_region_field(array_half)
    _, reference_half = focused_region_field(fine_half)
    assert field_correlation(coarse_lambda, reference_lambda) < 0.95
    assert field_correlation(coarse_half, reference_half) > 0.98
