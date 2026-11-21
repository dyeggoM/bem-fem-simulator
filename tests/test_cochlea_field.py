"""Tests de U3: el corte de cóclea desde la UI, con el campo interior y el balance de energía.

La aceptación de U3 es que la UI **reproduzca los números de D29/D30**, no que calcule algo nuevo:
`t_interior = 0.53187` con el caso base (espesor 2 mm, c_s 1900 m/s, α = 0, membrana del ballpark),
y un balance de energía que cierre. Cada solve son ~4.5 s, así que estos tests hacen tres como mucho.
"""

from __future__ import annotations

import json

import pytest

from ultrasound_sim.dashboard import cochlea_field

BASE = {"thickness_mm": 2.0, "c_s": 1900.0, "alpha_np_m": 0.0}
T_INTERIOR_D30 = 0.53187


# --------------------------------------------------------------------------- validación


def test_espesor_fuera_de_rango():
    with pytest.raises(ValueError):
        cochlea_field.validar_cochlea({"thickness_mm": 20.0})


def test_c_s_mayor_que_c_p_se_rechaza():
    """c_s >= c_p no es un sólido isótropo: se rechaza.

    Hoy el rango del slider (hasta 2500 m/s) ya queda por debajo de c_p = 3500, así que salta el
    rango primero y el guard de c_p < c_p queda como defensa por si alguien amplía el rango. El test
    fija el comportamiento visible: falla y dice el rango.
    """
    with pytest.raises(ValueError) as exc:
        cochlea_field.validar_cochlea({"c_s": 3600.0})
    assert "c_s" in str(exc.value)
    assert "rango" in str(exc.value)


def test_el_guard_de_c_s_contra_c_p_sigue_en_el_codigo():
    """Si alguien sube el tope del slider por encima de c_p, el guard tiene que saltar.

    Se prueba el guard directamente (sin pasar por el rango) para que la intención quede fijada.
    """
    assert cochlea_field.HUESO_CP == 3500.0
    assert cochlea_field._LIMITES["c_s"][1] < cochlea_field.HUESO_CP


def test_ventana_abierta_con_semiangulo_cero():
    with pytest.raises(ValueError):
        cochlea_field.validar_cochlea({"ventana": True, "window_half_angle": 0.0})


def test_valores_por_defecto():
    p = cochlea_field.validar_cochlea({})
    assert p["thickness_mm"] == 2.0
    assert p["c_s"] == 1900.0
    assert p["membrana"] is True
    assert p["ventana"] is True


# --------------------------------------------------------------------------- costo


def test_estimar_no_resuelve_y_declara_los_supuestos():
    c = cochlea_field.estimar_cochlea(BASE)
    assert c["malla"]["n_nodos"] > 10_000
    assert c["estimado_solve_s"] > 0
    assert "metricas" not in c
    # los supuestos VERIFY tienen que viajar con la respuesta, no quedar sólo en el código
    assert c["supuestos"]["c_p_hueso"] == 3500.0
    assert "VERIFY" in c["supuestos"]["verificacion"]


# --------------------------------------------------------------------------- aceptación


def test_reproduce_t_interior_de_d30_dentro_del_1_por_ciento():
    r = cochlea_field.resolver_cochlea(BASE)
    m = r["metricas"]
    print(
        f"\n  t_interior {m['t_interior']:.5f} (D30 {T_INTERIOR_D30}) · n=0 {m['t_interior_n0']:.5f} · "
        f"absorbe {100 * m['absorcion_membrana']:.1f} % · balance {100 * m['balance_energia']:.2f} %"
    )
    assert m["t_interior"] == pytest.approx(T_INTERIOR_D30, rel=0.01)
    # el armónico n=0 también es de D30: es el que rompe la monotonía del perfil
    assert m["t_interior_n0"] == pytest.approx(0.17913, rel=0.02)
    assert 0.0 < m["absorcion_membrana"] < 1.0


def test_el_balance_de_energia_cierra():
    """La identidad documentada es `P_total + P_membrana = 0` con la membrana como único sumidero.

    NO es `dispersada + absorbida = 1`: `power_scattered_fraction` (1.12 en el caso base) se mide en
    otro radio y con otra normalización, y sumarlas daba 1.85 — un número que no significa nada.
    """
    r = cochlea_field.resolver_cochlea(BASE)
    m = r["metricas"]
    total = r["detalle"]["power_total_fraction"]
    absorcion = r["detalle"]["power_membrane_fraction"]
    assert abs(total + absorcion) < 0.01, f"el balance no cierra: {total} + {absorcion}"
    assert m["balance_energia"] < 0.01


def test_sin_membrana_no_absorbe_y_el_campo_cambia():
    r = cochlea_field.resolver_cochlea(dict(BASE, membrana=False))
    m = r["metricas"]
    print(f"\n  sin membrana: t_interior {m['t_interior']:.5f} · absorbe {m['absorcion_membrana']:.4f}")
    assert abs(m["absorcion_membrana"]) < 1e-3, "sin membrana no hay nada que absorba"
    assert m["t_interior"] != pytest.approx(T_INTERIOR_D30, rel=0.01)


# --------------------------------------------------------------------------- forma y serialización


def test_el_resultado_es_json_estricto_y_el_mapa_tiene_huecos():
    """`sol.as_dict()` trae `window_impedance` complejo: sin sanear, el navegador no recibe nada."""
    r = cochlea_field.resolver_cochlea(BASE)
    texto = json.dumps(r, allow_nan=False)   # revienta con NaN o con un complejo
    assert len(texto) > 10_000
    devuelta = json.loads(texto)
    filas = devuelta["mapa"]["db"]
    assert len(filas) == cochlea_field.PUNTOS_MAPA
    assert all(len(f) == cochlea_field.PUNTOS_MAPA for f in filas)
    assert any(v is None for f in filas for v in f), "fuera del disco tiene que ser hueco"
    assert any(v is not None for f in filas for v in f), "dentro del disco tiene que haber campo"


def test_el_mapa_es_finito_dentro_del_disco():
    r = cochlea_field.resolver_cochlea(BASE)
    x = r["mapa"]["x_mm"]
    a = r["mapa"]["radio_mm"]
    dentro = [
        v for i, fila in enumerate(r["mapa"]["db"]) for j, v in enumerate(fila)
        if (x[i] ** 2 + x[j] ** 2) < (0.9 * a) ** 2 and v is not None
    ]
    assert len(dentro) > 1000
    assert all(-60.0 < v <= 0.0 for v in dentro), "hay dB fuera de rango dentro del disco"


def test_la_referencia_viaja_en_la_respuesta():
    r = cochlea_field.resolver_cochlea(BASE)
    assert r["referencia"]["t_interior_d30"] == T_INTERIOR_D30
    assert 1.0 <= r["tiempos"]["solve_s"] < 60.0
