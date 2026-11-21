"""Barrido de apertura: FWHM lateral contra el paso entre elementos (D25/D26 medido en la UI).

La pregunta de diseño del proyecto es "¿qué pitch necesito para que el FWHM no se despegue de una
apertura continua?". El barrido la responde resolviendo un campo por punto.

Ojo con el caso que este test ataja: la primera versión **no barría**. Usaba las claves `n` y `d` en el
payload de `calcular`, que espera `n_elements` y `aperture_mm`; no fallaba, simplemente las ignoraba y
devolvía el mismo campo cinco veces (FWHM 9.12 mm y sidelobe −9.24 dB idénticos con 8 y con 64
elementos). Un barrido que no barre es peor que no tenerlo, porque parece un resultado.
"""

from __future__ import annotations

import pytest

from ultrasound_sim.dashboard import sweeps


@pytest.fixture(scope="module")
def barrido():
    return sweeps.barrido_apertura({
        "desde": 8, "hasta": 64, "puntos": 4, "apertura_mm": 60.0, "foco_mm": 90.0,
        "frecuencia_khz": 250.0, "material": "agua", "modelo": "baffle",
    })


def test_el_barrido_barre_de_verdad(barrido):
    """El FWHM tiene que CAMBIAR con el pitch. El caso del bug era un FWHM constante."""
    fwhm = [f["fwhm_lateral_mm"] for f in barrido["filas"]]
    print(f"\n  FWHM por pitch: {[round(v, 3) for v in fwhm]}")
    assert len(set(round(v, 6) for v in fwhm)) > 1, "todos los FWHM iguales: no está barriendo"
    # con pitch grueso se despega de la referencia continua
    assert fwhm[0] > fwhm[-1]


def test_el_pitch_grueso_pierde_contra_la_apertura_continua(barrido):
    """Con `pitch = 1.25 λ` el FWHM se desvía más del 2 %; con pitch fino empata."""
    filas = barrido["filas"]
    grueso, fino = filas[0], filas[-1]
    print(f"\n  pitch {grueso['pitch_over_lambda']:.2f} λ → desvío "
          f"{100 * grueso['desvio_vs_continuo']:.1f} % · pitch {fino['pitch_over_lambda']:.2f} λ → "
          f"{100 * fino['desvio_vs_continuo']:.1f} %")
    assert grueso["desvio_vs_continuo"] > 0.05
    assert fino["desvio_vs_continuo"] < 0.02


def test_el_pitch_que_cumple_el_2pct_queda_cerca_de_lambda_4(barrido):
    """D25/D26 dicen `pitch ≤ λ/4` para ~2 %. El barrido tiene que caer en ese entorno, no en otro
    orden de magnitud."""
    pitch = barrido["resumen"]["pitch_con_2pct"]
    print(f"\n  pitch que cumple el 2 %: {pitch:.3f} λ (λ/4 = 0.25)")
    assert pitch is not None
    assert 0.15 <= pitch <= 0.60


def test_la_referencia_continua_se_resuelve_y_es_mas_fina(barrido):
    ref = barrido["referencia_continua"]
    print(f"\n  referencia (n=400): FWHM {ref['fwhm_lateral_mm']:.3f} mm, pitch "
          f"{ref['pitch_over_lambda']:.4f} λ")
    assert ref["fwhm_lateral_mm"] > 0
    assert ref["fwhm_lateral_mm"] <= barrido["filas"][-1]["fwhm_lateral_mm"] + 1e-9


def test_los_parametros_fuera_de_rango_se_rechazan():
    with pytest.raises(ValueError):
        sweeps.barrido_apertura({"desde": 64, "hasta": 8})
    with pytest.raises(ValueError):
        sweeps.barrido_apertura({"puntos": 1})
    with pytest.raises(ValueError):
        sweeps.barrido_apertura({"apertura_mm": 1.0})
