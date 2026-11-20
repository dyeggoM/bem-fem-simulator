"""Leyes de atenuación: conversión correcta, procedencia verificada y modelos rechazados.

Lo que se prueba acá es que la conversión dB/cm ↔ Np/m no esté al revés y que las tres fuentes
reproducen sus números publicados. El caso "Pichardo 27.31 f^1.722" queda documentado como rechazado.
"""

import numpy as np
import pytest

from ultrasound_sim import attenuation, materials
from ultrasound_sim.analytic import helmholtz1d


def test_db_cm_to_np_m_is_the_amplitude_conversion():
    # 20·log10(e) = 8.6859 dB por neper de amplitud
    assert attenuation.db_cm_to_np_m(13.3) == pytest.approx(153.1210, rel=1e-5)
    assert attenuation.np_m_to_db_cm(153.1210) == pytest.approx(13.3, rel=1e-5)
    # ida y vuelta
    for a in (0.5, 8.83, 27.3):
        assert attenuation.np_m_to_db_cm(attenuation.db_cm_to_np_m(a)) == pytest.approx(a)


def test_pinton_matches_the_published_value_and_the_hand_computed_one():
    law = attenuation.SKULL_PINTON_2012
    assert law.alpha_db_cm(1.0e6) == pytest.approx(13.3, rel=1e-6)
    assert law.alpha(1.0e6) == pytest.approx(153.121, rel=1e-5)
    # el número que aparecía en el documento compartido (76.6 Np/m a 500 kHz)
    assert law.alpha(500e3) == pytest.approx(76.56, rel=1e-3)
    # lineal en f
    assert law.alpha(2e6) == pytest.approx(2 * law.alpha(1e6), rel=1e-12)
    assert law.verified is True and "Pinton" in law.provenance


def test_db_loss_is_the_same_for_amplitude_and_intensity():
    """Corregido tras la consulta a Sim-Datos: el dB no se duplica, se duplica el Np.

    p = p0·exp(−αz) y I = I0·exp(−2αz) con el MISMO α en Np/m. Pero dB = 20·log(p) = 10·log(I),
    así que la pérdida en dB por cm es idéntica: 13.3 dB/cm de amplitud son 13.3 dB/cm de
    intensidad. Duplicar los dB es un error de 2× que ningún test numérico del solver detectaría.
    """
    law = attenuation.SKULL_PINTON_2012
    d = 6e-3
    amp_db = 20.0 * np.log10(np.exp(-law.alpha(1e6) * d))
    assert law.loss_db(1e6, d) == pytest.approx(-amp_db, rel=1e-9)
    # 13.3 dB/cm ⇒ 13.3 dB en 1 cm (no 26.6)
    assert law.loss_db(1e6, 0.01) == pytest.approx(13.3, rel=1e-4)
    # el exponente de intensidad en Np sí es el doble, y eso NO cambia los dB
    alpha_int = 2.0 * law.alpha(1e6)
    assert 10.0 * np.log10(np.exp(alpha_int * 0.01)) == pytest.approx(13.3, rel=1e-4)


def test_pichardo_measured_points_are_reproduced_exactly():
    pts = attenuation.SKULL_PICHARDO_2011_POINTS
    for f, a in zip(pts.frequencies_hz, pts.alphas_np_m):
        assert pts.alpha(f) == pytest.approx(float(a), rel=1e-12)
    # interpolación monótona entre puntos
    mid = pts.alpha(500e3)
    assert float(pts.alphas_np_m[0]) < mid < float(pts.alphas_np_m[1])


def test_pichardo_data_are_not_a_clean_power_law():
    """Hallazgo, no un fallo: sus tres puntos no son colineales en log-log.

    270→836 kHz exige n ≈ 1.68 y 836→1402 kHz n ≈ 0.72. Un ajuste global único se desvía ~20 % en el
    medio, así que la recomendación es usar los puntos medidos (interpolación por tramos), no la ley.
    """
    law = attenuation.SKULL_PICHARDO_2011_FIT
    pts = attenuation.SKULL_PICHARDO_2011_POINTS
    desvios = [(law.alpha(f) / float(a) - 1.0) for f, a in
               zip(pts.frequencies_hz, pts.alphas_np_m)]
    print(f"\n  ajuste global Pichardo: n = {law.exponent:.2f} "
          f"({law.alpha_db_cm(1e6):.1f} dB/cm a 1 MHz); desvíos en los nodos = "
          f"{[round(float(d) * 100, 1) for d in desvios]} %")
    # el peor desvío está entre 15 y 30 %: es lo que hay, y por eso no se usa como ley
    assert 0.15 < max(abs(d) for d in desvios) < 0.30
    # la interpolación por tramos, en cambio, es exacta en los nodos
    for f, a in zip(pts.frequencies_hz, pts.alphas_np_m):
        assert pts.alpha(f) == pytest.approx(float(a), rel=1e-12)


def test_the_rejected_model_really_does_not_fit_the_source():
    """27.31·f^1.722 daría ~49 dB/cm a 1.4 MHz; Pichardo mide 29. Queda fuera de la librería."""
    a, n = 27.31, 1.722
    predicho = a * 1.402**n
    medido = attenuation.np_m_to_db_cm(334.0)
    print(f"\n  27.31·f^1.722 → {predicho:.1f} dB/cm a 1.4 MHz; medido {medido:.1f} dB/cm")
    assert medido == pytest.approx(29.0, abs=0.5)
    assert predicho > 1.5 * medido, "el modelo rechazado no debe parecerse a los datos"
    assert "pichardo_27.31_f^1.722" in attenuation.REJECTED_MODELS


def test_skull_material_feeds_the_solver_consistently():
    m = attenuation.skull_material(attenuation.SKULL_PINTON_2012)
    assert m.alpha0 == pytest.approx(153.121, rel=1e-5)
    assert m.alpha_power == pytest.approx(1.0)
    assert m.f_ref == pytest.approx(1e6)
    # el solver usa k = ω/c + i·α(f): la parte imaginaria es la atenuación de amplitud
    assert m.wavenumber(1e6).imag == pytest.approx(153.121, rel=1e-5)
    # la procedencia tiene que dejar claro que es cráneo "equivalente", no cortical pura
    assert "cortical pura" in m.provenance
    assert "2800" in m.provenance and "1900" in m.provenance


def test_white_and_pichardo_disagree_and_are_not_averaged():
    """Dos estudios serios del mismo hueso difieren 2.4×–3.3×. Promediarlos sería inventar."""
    p = attenuation.SKULL_PICHARDO_2011_POINTS
    w = attenuation.SKULL_WHITE2006_POINTS
    bajo = float(p.alpha(272e3)) / float(w.alpha(272e3))
    alto = float(p.alpha(840e3)) / float(w.alpha(840e3))
    print(f"\n  Pichardo/White = {bajo:.2f}× a 272 kHz y {alto:.2f}× a 840 kHz")
    assert 2.0 < bajo < 3.0
    assert 3.0 < alto < 4.0
    # y cada uno carga su procedencia, para que nadie los mezcle por accidente
    assert "White" in w.provenance and "Pichardo" in p.provenance
    assert "no son intercambiables" in w.provenance


def test_no_dB_doubling_in_the_reported_losses():
    """Regresión del error de la consulta: la pérdida en dB es la misma para amplitud e intensidad."""
    law = attenuation.SKULL_PINTON_2012
    for d in (1e-3, 6e-3, 2e-2):
        assert law.loss_db(1e6, d) == pytest.approx(
            20.0 * np.log10(np.exp(law.alpha(1e6) * d)), rel=1e-12)


@pytest.mark.parametrize(
    "law", [attenuation.SKULL_PINTON_2012, attenuation.SKULL_WHITE2006]
)
def test_transmission_through_6mm_of_skull_with_and_without_losses(law):
    """El pago: cuánto cuesta ignorar la atenuación (placa de 6 mm, incidencia normal)."""
    water = materials.WATER
    thickness = 6e-3
    rows = []
    for f in (100e3, 500e3, 1e6):
        skull_lossless = materials.Material(
            name="skull_sin_perdidas", c=2800.0, rho=1900.0, provenance="control sin pérdidas")
        skull_lossy = attenuation.skull_material(law, c=2800.0, rho=1900.0)
        t_lossless = abs(helmholtz1d.stack_response(
            f, [water, skull_lossless, water], [thickness]).T)
        t_lossy = abs(helmholtz1d.stack_response(
            f, [water, skull_lossy, water], [thickness]).T)
        rows.append((f, t_lossless, t_lossy))
        assert 0.0 < t_lossy < t_lossless, "la atenuación sólo puede sacar energía"
        print(f"    {f/1e3:4.0f} kHz: |T| sin pérdidas {t_lossless:.4f} → con pérdidas "
              f"{t_lossy:.4f}  ({20*np.log10(t_lossy/t_lossless):+.2f} dB)")
    # a 1 MHz la atenuación tiene que morder de verdad: al menos 3 dB de caída en amplitud
    f_last, t_ll, t_l = rows[-1]
    assert 20 * np.log10(t_l / t_ll) < -3.0, rows[-1]


def test_decreasing_frequency_helps_but_does_not_remove_the_loss():
    law = attenuation.SKULL_PINTON_2012
    d = 6e-3
    losses = {f: law.intensity_loss_db(f, d) for f in (100e3, 250e3, 500e3, 1e6)}
    print(f"\n  pérdida de intensidad en 6 mm: "
          + ", ".join(f"{f/1e3:.0f} kHz → {v:.1f} dB" for f, v in losses.items()))
    assert losses[100e3] < losses[250e3] < losses[500e3] < losses[1e6]
    assert losses[100e3] > 0.5, "a 100 kHz, 6 mm de cráneo cuestan más de medio dB"
    assert losses[1e6] > 7.0, "a 1 MHz son más de 7 dB: ahí la atenuación manda"
