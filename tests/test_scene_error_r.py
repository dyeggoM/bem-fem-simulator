"""La escena reporta el error de R como amplitud, y el fasor aparte con su desfase explicado (D52).

Historia: `error_R_pct` daba 100-200 % al lado de un sello VALIDADO con 0.19 % de error de amplitud.
La causa eran dos cosas: R se comparaba como fasor y T como módulo (descuido), y los fasores estaban
referidos a orígenes distintos (x=0 el numérico, la primera interfaz el analítico), con un giro de
propagación e^{2ik t1} en medio. Se comprobó contra |ΔR|/|R| = 2·sin(k·t1) en cuatro frecuencias.
"""

from __future__ import annotations

import numpy as np
import pytest

from ultrasound_sim.dashboard import scene

BASE = {
    "ppw": 8, "alto_m": 0.04,
    "capas": [{"material": "tejido", "espesor_m": 0.020}, {"material": "hueso", "espesor_m": 0.010}],
}


@pytest.fixture(scope="module")
def corridas():
    return {
        f: scene.resolver_escena({**BASE, "frecuencia_hz": f * 1e3})["metricas"]
        for f in (100.0, 250.0)
    }


def test_el_error_de_R_que_se_muestra_es_de_amplitud(corridas):
    """El número comparable con `error_T_pct` tiene que ser el del módulo, no el del fasor. Antes daba
    199 % donde la amplitud coincidía al 0.19 %."""
    for f, m in corridas.items():
        print(f"\n  {f:.0f} kHz: R módulo {m['error_R_pct']:.3f} % · fasor {m['error_R_fasor_pct']:.2f} %"
              f" · error_R_abs {m['error_R_abs']:.5f}")
        assert m["error_R_pct"] < 5.0, "el error de amplitud no puede ser de cientos por ciento"
        # y tiene que ser coherente con el error absoluto: mismo número en magnitudes parecidas
        assert abs(m["error_R_pct"] / 100.0 - m["error_R_abs"] / m["R_analitico"]) < 0.01


def test_el_desfase_de_referencia_se_reporta_y_coincide_con_la_propagacion(corridas):
    """El desfase que se sacó para alinear los fasores es `2·k·t₁`: si no coincidiera, la alineación
    estaría corrigiendo otra cosa (y tapando un error)."""
    material = scene.MATERIALES["tejido"]
    for f, m in corridas.items():
        k = abs(complex(material.wavenumber(f * 1e3)))
        esperado = np.degrees(2.0 * k * 0.020)
        print(f"\n  {f:.0f} kHz: desfase {m['desfase_referencia_R_deg']:.1f}°"
              f" · 2k·t₁ {esperado:.1f}°")
        assert m["desfase_referencia_R_deg"] == pytest.approx(esperado, rel=1e-6)


def test_el_fasor_corregido_es_del_orden_del_error_de_fase_de_la_malla(corridas):
    """El fasor alineado queda en el orden del error de FASE (unos grados con ppw 8), no en cientos
    por ciento. Se acota por arriba para que no vuelva a pasar desapercibido un 199 %."""
    for f, m in corridas.items():
        print(f"\n  {f:.0f} kHz: fasor {m['error_R_fasor_pct']:.2f} %")
        assert m["error_R_fasor_pct"] < 30.0
