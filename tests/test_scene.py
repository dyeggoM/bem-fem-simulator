"""Tests de la escena 2D (Fase U1): validación, malla y R/T contra el analítico.

El criterio de aceptación de U1 es el de la Fase 1: R y T de una capa contra la matriz de
transferencia **dentro del 2 % en amplitud**. Se cumple con ppw >= 10 en los tres casos medidos,
y el caso de contraste fuerte (agua/hueso/agua, R = 0.90) es el que de verdad ejerce la métrica:
con agua/tejido R vale 0.012 y el error absoluto que se mide (0.005 a ppw 16) es casi todo piso de
discretización, no error de física.
"""

from __future__ import annotations

import numpy as np
import pytest

from ultrasound_sim.dashboard import scene


def _capas(*pares):
    return [{"material": m, "espesor_m": e} for m, e in pares]


ESCENA_AGUA_TEJIDO = {
    "frecuencia_hz": 100.0e3,
    "ppw": 12,
    "alto_m": 20.0e-3,
    "capas": _capas(("agua", 20.0e-3), ("tejido", 15.0e-3), ("agua", 20.0e-3)),
}

ESCENA_AGUA_HUESO = {
    "frecuencia_hz": 100.0e3,
    "ppw": 12,
    "alto_m": 20.0e-3,
    "capas": _capas(("agua", 20.0e-3), ("hueso", 10.0e-3), ("agua", 20.0e-3)),
}


# --------------------------------------------------------------------------- validación


def test_material_desconocido_falla_con_las_opciones():
    malo = dict(ESCENA_AGUA_TEJIDO, capas=_capas(("agua", 10e-3), ("unobtainium", 5e-3)))
    with pytest.raises(ValueError) as exc:
        scene.validar_escena(malo)
    assert "unobtainium" in str(exc.value)
    assert "agua" in str(exc.value)


def test_incidencia_oblicua_se_rechaza_con_la_razon():
    """No es una limitación de alcance: arriba y abajo las paredes son rígidas y la onda se
    reflejaría ahí. Rechazarla con el motivo es mejor que devolver un campo que parece válido."""
    with pytest.raises(ValueError) as exc:
        scene.validar_escena(dict(ESCENA_AGUA_TEJIDO, angulo_deg=20.0))
    assert "incidencia normal" in str(exc.value)


def test_malla_gigante_se_rechaza_antes_de_colgar():
    grande = dict(ESCENA_AGUA_TEJIDO, ppw=40, alto_m=100e-3)
    # el tope vive en malla_de_escena: validar la escena no arma la malla, así que un ppw alto
    # pasa la validación y revienta al mallar. Es a propósito: el costo se mide, no se supone.
    with pytest.raises(ValueError) as exc:
        scene.malla_de_escena(scene.validar_escena(grande))
    mensaje = str(exc.value)
    assert "tope interactivo" in mensaje
    assert "ppw" in mensaje


def test_cantidad_de_capas_fuera_de_rango():
    with pytest.raises(ValueError):
        scene.validar_escena(dict(ESCENA_AGUA_TEJIDO, capas=_capas(("agua", 10e-3))))


def test_espesor_no_positivo():
    with pytest.raises(ValueError):
        scene.validar_escena(
            dict(ESCENA_AGUA_TEJIDO, capas=_capas(("agua", 10e-3), ("tejido", 0.0)))
        )


def test_nan_en_frecuencia():
    with pytest.raises(ValueError):
        scene.validar_escena(dict(ESCENA_AGUA_TEJIDO, frecuencia_hz=float("nan")))


# --------------------------------------------------------------------------- malla


def test_malla_es_cerrada_conserva_area_y_respeta_las_interfaces():
    escena = scene.validar_escena(ESCENA_AGUA_TEJIDO)
    xs, _ys, malla, capa = scene.malla_de_escena(escena)

    # área total = suma de las áreas elementales (sin solapes ni huecos)
    from ultrasound_sim.meshing import triangle_areas

    area_el = triangle_areas(malla.nodes, malla.triangles)
    assert np.isclose(area_el.sum(), xs[-1] * escena.alto_m, rtol=1e-12)

    # las interfaces caen exactamente en xs (la malla es conforme por construcción)
    interfaces = np.cumsum(escena.espesores())[:-1]
    for interfaz in interfaces:
        assert np.any(np.isclose(xs, interfaz, rtol=0.0, atol=1e-12))

    # número de capas distintas en la asignación por elemento
    assert set(capa.tolist()) == {0, 1, 2}
    assert capa.min() == 0 and capa.max() == 2


def test_malla_es_reproducible_bit_a_bit():
    """Guardar y cargar la escena tiene que devolver la MISMA malla: si cambia el conteo de
    elementos cambian los números y se pierde la comparación entre corridas (D33)."""
    a = scene.validar_escena(ESCENA_AGUA_TEJIDO)
    b = scene.validar_escena(dict(ESCENA_AGUA_TEJIDO))  # mismo dict, otro objeto
    c = scene.validar_escena(a.a_dict())  # ida y vuelta por JSON

    xs_a, ys_a, malla_a, capa_a = scene.malla_de_escena(a)
    xs_c, ys_c, malla_c, capa_c = scene.malla_de_escena(c)

    assert np.array_equal(xs_a, xs_c)
    assert np.array_equal(ys_a, ys_c)
    assert np.array_equal(malla_a.nodes, malla_c.nodes)
    assert np.array_equal(malla_a.triangles, malla_c.triangles)
    assert np.array_equal(capa_a, capa_c)
    assert np.array_equal(malla_a.nodes, scene.malla_de_escena(b)[2].nodes)


def test_ida_y_vuelta_por_json_conserva_la_escena():
    import json

    escena = scene.validar_escena(ESCENA_AGUA_TEJIDO)
    texto = json.dumps(escena.a_dict())
    vuelta = scene.validar_escena(json.loads(texto))
    assert vuelta.a_dict() == escena.a_dict()


# --------------------------------------------------------------------------- costo


def test_estimar_costo_no_resuelve_y_da_lo_mismo_que_la_malla():
    costo = scene.estimar_costo(ESCENA_AGUA_TEJIDO)
    _xs, _ys, malla, _ = scene.malla_de_escena(scene.validar_escena(ESCENA_AGUA_TEJIDO))
    assert costo["grilla"]["n_elementos"] == malla.triangles.shape[0]
    assert costo["grilla"]["n_nodos"] == malla.nodes.shape[0]
    assert costo["estimado_solve_s"] > 0.0
    assert len(costo["grilla"]["interfaces_mm"]) == 2


# --------------------------------------------------------------------------- R/T (aceptación U1)


@pytest.mark.parametrize("ppw", [10, 12, 16])
def test_rt_agua_hueso_agua_contraste_fuerte(ppw):
    """Caso que de verdad ejerce la métrica: R = 0.90, T = 0.43."""
    resultado = scene.resolver_escena(dict(ESCENA_AGUA_HUESO, ppw=ppw))
    m = resultado["metricas"]
    print(
        f"\n  ppw={ppw}: R {m['R']:.4f} vs {m['R_analitico']:.4f} (abs {m['error_R_abs']:.4f}) | "
        f"T {m['T']:.4f} vs {m['T_analitico']:.4f} (abs {m['error_T_abs']:.4f}) | "
        f"energía {m['energia']:.4f}"
    )
    assert m["error_R_abs"] < 0.02, f"R fuera del 2 % en amplitud: {m['error_R_abs']:.4f}"
    assert m["error_T_abs"] < 0.02, f"T fuera del 2 % en amplitud: {m['error_T_abs']:.4f}"
    assert m["energia"] == pytest.approx(1.0, abs=0.02)


@pytest.mark.parametrize("ppw", [12, 16])
def test_rt_agua_tejido_agua_dentro_del_2_por_ciento(ppw):
    """El caso del criterio U1, con el material real del proyecto."""
    resultado = scene.resolver_escena(dict(ESCENA_AGUA_TEJIDO, ppw=ppw))
    m = resultado["metricas"]
    print(
        f"\n  ppw={ppw}: R {m['R']:.4f} vs {m['R_analitico']:.4f} | "
        f"T {m['T']:.4f} vs {m['T_analitico']:.4f} | energía {m['energia']:.4f}"
    )
    assert m["error_R_abs"] < 0.02
    assert m["error_T_abs"] < 0.02


def test_rt_mejora_o_se_mantiene_al_refinar():
    """Refinar no puede empeorar el error de T: es el sanity de que la malla está sirviendo."""
    errores = []
    for ppw in (8, 12, 20):
        m = scene.resolver_escena(dict(ESCENA_AGUA_HUESO, ppw=ppw))["metricas"]
        errores.append(m["error_T_abs"])
    assert errores[-1] <= errores[0] + 1e-6, f"T empeora al refinar: {errores}"


# --------------------------------------------------------------------------- forma del resultado


def test_resultado_tiene_lo_que_la_ui_necesita():
    resultado = scene.resolver_escena(ESCENA_AGUA_TEJIDO)
    for clave in ("ok", "escena", "grilla", "campo", "cortes", "metricas", "tiempos"):
        assert clave in resultado
    n_y = len(resultado["grilla"]["y_mm"])
    n_x = len(resultado["grilla"]["x_mm"])
    assert np.asarray(resultado["campo"]["db"]).shape == (n_y, n_x)
    assert np.asarray(resultado["campo"]["amplitud"]).shape == (n_y, n_x)
    assert len(resultado["cortes"]["axial"]) == n_x
    assert len(resultado["cortes"]["lateral_en_x_pico"]) == n_y
    assert resultado["tiempos"]["solve_s"] > 0.0
    assert np.isfinite(np.asarray(resultado["campo"]["db"])).all()


def test_el_campo_no_es_nan_ni_inf_en_ninguna_escena_valida():
    for escena_datos in (ESCENA_AGUA_TEJIDO, ESCENA_AGUA_HUESO):
        resultado = scene.resolver_escena(escena_datos)
        campo = np.asarray(resultado["campo"]["amplitud"])
        assert np.isfinite(campo).all()
        assert campo.max() > 0.0
