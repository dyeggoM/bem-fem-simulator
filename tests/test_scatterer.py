"""Tests de U2: dispersor cilíndrico por BEM, FEM y la serie, con la discrepancia a la vista.

Criterio de aceptación (el de la Fase 4): los tres caminos tienen que coincidir **dentro del 2 %**
en el campo evaluado sobre un círculo de radio 1.5a. Los números de referencia medidos en D12
están en el docstring de cada test.

El test que más vale acá no es "FEM ≈ BEM" —eso podría pasar con los dos equivocados— sino
"cada uno contra la serie", que es la verdad de referencia.
"""

from __future__ import annotations

import numpy as np

import json

import pytest

from ultrasound_sim.dashboard import scatterer

CRITERIO = 2.0  # %


# --------------------------------------------------------------------------- validación


def test_ka_fuera_de_rango_se_rechaza():
    with pytest.raises(ValueError):
        scatterer.validar_dispersor({"ka": 20.0})


def test_ecuacion_desconocida_falla_con_las_opciones():
    with pytest.raises(ValueError) as exc:
        scatterer.validar_dispersor({"equation": "galerkin"})
    assert "standard" in str(exc.value)


def test_abc_order_invalido():
    with pytest.raises(ValueError):
        scatterer.validar_dispersor({"abc_order": 3})


def test_paneles_auto_se_normaliza_a_none():
    assert scatterer.validar_dispersor({"n_panels": "auto"})["n_panels"] is None
    assert scatterer.validar_dispersor({})["n_panels"] is None
    assert scatterer.validar_dispersor({"n_panels": 32})["n_panels"] == 32


def test_malla_gigante_se_rechaza_con_que_bajar():
    with pytest.raises(ValueError) as exc:
        scatterer.estimar_dispersor({"ka": 5.0, "ppw": 40, "r_out": 6.0})
    assert "tope interactivo" in str(exc.value)
    assert "ppw" in str(exc.value)


# --------------------------------------------------------------------------- costo


def test_estimar_no_resuelve_y_da_la_geometria():
    c = scatterer.estimar_dispersor({"ka": 1.0})
    assert c["malla"]["n_elementos"] > 0
    assert c["paneles_bem"] > 0
    # ka = k·a = 1 con k de agua a 100 kHz: a = λ/2π = 15/6.283 mm
    assert c["geo"]["radio_mm"] == pytest.approx(2.387, abs=0.01)
    assert c["geo"]["longitud_onda_mm"] == pytest.approx(15.0, abs=1e-6)
    assert "metricas" not in c


# --------------------------------------------------------------------------- aceptación


@pytest.mark.parametrize("ka, esperado", [(1.0, (0.97, 0.94, 0.05)), (2.0, (0.58, 0.60, 0.14))])
def test_los_tres_metodos_coinciden_dentro_del_2_por_ciento(ka, esperado):
    """Cifras de D12: ka=1 → 0.97/0.94/0.05; ka=2 → 0.59/0.60/0.14 (medido 0.58/0.60/0.14)."""
    r = scatterer.resolver_dispersor({"ka": ka, "ppw": 20, "r_out": 3.0})
    m = r["metricas"]
    print(
        f"\n  ka={ka}: FEM−BEM {m['fem_vs_bem_pct']:.2f} % | FEM−serie {m['fem_vs_serie_pct']:.2f} % | "
        f"BEM−serie {m['bem_vs_serie_pct']:.2f} %"
    )
    assert m["fem_vs_bem_pct"] < CRITERIO
    assert m["fem_vs_serie_pct"] < CRITERIO
    assert m["bem_vs_serie_pct"] < CRITERIO
    # y además cerca de lo ya publicado: si se mueve, quiero enterarme
    e_fb, e_fs, e_bs = esperado
    assert m["fem_vs_bem_pct"] == pytest.approx(e_fb, abs=0.15)
    assert m["fem_vs_serie_pct"] == pytest.approx(e_fs, abs=0.15)
    assert m["bem_vs_serie_pct"] == pytest.approx(e_bs, abs=0.15)


def test_burton_miller_y_chief_tambien_dan_el_mismo_campo():
    """Las tres ecuaciones resuelven la misma EDP: a ka=1 sin resonancia tienen que coincidir.

    BM con eta fijo empeora el condicionamiento (medido: 2.37 → 12.78) y por eso su error contra
    la serie sube de 0.05 % a 0.24 %, pero sigue muy dentro del criterio. Es el precio conocido de
    regularizar para poder cruzar resonancias en un barrido.
    """
    base = scatterer.resolver_dispersor({"ka": 1.0, "equation": "standard"})
    for ecuacion in ("burton_miller", "chief"):
        r = scatterer.resolver_dispersor({"ka": 1.0, "equation": ecuacion})
        m = r["metricas"]
        print(
            f"\n  {ecuacion}: FEM−BEM {m['fem_vs_bem_pct']:.2f} % | BEM−serie "
            f"{m['bem_vs_serie_pct']:.2f} % | cond {m['condicion']:.2f}"
        )
        assert m["bem_vs_serie_pct"] < CRITERIO
        assert m["fem_vs_bem_pct"] < CRITERIO
        assert r["metricas"]["max_bem"] == pytest.approx(base["metricas"]["max_bem"], rel=0.02)


# --------------------------------------------------------------------------- forma y serialización


def test_el_resultado_es_json_valido_sin_nan():
    """El mapa tiene huecos dentro del cilindro, donde no hay campo.

    Si se manda `NaN` (que es lo que hace `json.dumps` por defecto con un float NaN) el navegador
    rechaza el cuerpo ENTERO con "Unexpected token 'N'". Se manda `None` → `null`, y este test lo
    fija: `allow_nan=False` hace que cualquier NaN futuro reviente acá en vez de en el navegador.
    """
    r = scatterer.resolver_dispersor({"ka": 1.0})
    texto = json.dumps(r, allow_nan=False)  # revienta si quedó un NaN
    assert len(texto) > 1000
    devuelta = json.loads(texto)
    filas = devuelta["mapa"]["db"]
    planos = [v for fila in filas for v in fila if v is None]
    assert planos, "el centro del mapa tiene que ser hueco (dentro del cilindro)"
    assert all(v is None or isinstance(v, (int, float)) for fila in filas for v in fila)


def test_mapa_y_comparacion_tienen_la_forma_que_la_ui_necesita():
    r = scatterer.resolver_dispersor({"ka": 1.0})
    n = scatterer.PUNTOS_MAPA
    assert len(r["mapa"]["x_mm"]) == n
    assert len(r["mapa"]["db"]) == n
    assert all(len(fila) == n for fila in r["mapa"]["db"])
    c = r["comparacion"]
    assert len(c["angulo_deg"]) == len(c["fem"]) == len(c["bem"]) == len(c["serie"])
    assert min(c["angulo_deg"]) >= -180.0 and max(c["angulo_deg"]) <= 180.0
    for clave in ("fem_s", "bem_s", "serie_s", "mapa_s", "total_s"):
        assert r["tiempos"][clave] >= 0.0
    assert r["malla"]["n_paneles_bem"] > 0


def test_el_campo_dispersado_satisface_la_condicion_de_radiacion():
    """El campo dispersado tiene que caer al alejarse; el **total** no, y por eso no se usa.

    El mapa que muestra la UI es ``|p_total|``, que lejos tiende a la incidente (|p| → 1) y por eso
    no decae: un test de decaimiento sobre el total falla por física correcta, no por un bug (lo
    comprobé: la esquina del mapa da −2.9 dB contra −3.8 dB «cerca», que es exactamente eso).
    La condición de radiación se chequea sobre ``p_scattered``, a dos radios.
    """
    from ultrasound_sim import materials
    from ultrasound_sim.bem import helmholtz2d as bem

    k = materials.WATER.wavenumber(scatterer.FREQUENCY)
    radio = 1.0 / abs(k)
    cerca = 1.6 * radio
    lejos = 4.0 * radio
    th = np.linspace(0.0, 2.0 * np.pi, 48, endpoint=False)
    p_cerca = bem.solve_exterior_neumann(
        np.column_stack([cerca * np.cos(th), cerca * np.sin(th)]),
        scatterer.FREQUENCY, materials.WATER, radio,
    ).p_scattered
    p_lejos = bem.solve_exterior_neumann(
        np.column_stack([lejos * np.cos(th), lejos * np.sin(th)]),
        scatterer.FREQUENCY, materials.WATER, radio,
    ).p_scattered
    ratio = float(np.linalg.norm(p_lejos) / np.linalg.norm(p_cerca))
    print(f"\n  |p_s(4a)| / |p_s(1.6a)| = {ratio:.3f}")
    assert ratio < 0.75, "el campo dispersado no decae al alejarse (¿falta la condición de radiación?)"
    assert ratio > 0.05, "decae demasiado: ¿se está apagando el campo por otra razón?"
