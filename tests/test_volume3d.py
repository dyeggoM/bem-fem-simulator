"""Tests de U4: escena 3D por rebanadas.

Dos modos y dos criterios:
- **puntual**: el campo tiene que reproducir el Green 3D (p = -iωρQ·e^{ikr}/(4πr)) dentro del 10 %
  en un anillo, con amplitud ≈ 1. D32 dio 8.4 % a ppw 8 en una caja de 1.1 λ.
- **plano**: R y T de la capa contra la matriz de transferencia, dentro del 2 % en amplitud (Fase 1).

El primer intento del modo puntual dio 48 % de error y **no bajaba al refinar**: la caja arrancaba en
x = 0 y la fuente (el nodo más cercano al origen) caía justo sobre la cara de la ABC. Se centraliza
la caja; acá queda fijado con una aserción sobre la posición de la fuente.
"""

from __future__ import annotations

import json

import numpy as np
import pytest

from ultrasound_sim.dashboard import volume3d

PUNTUAL = {
    "frecuencia_hz": 100.0e3,
    "ppw": 10,   # a ppw 8 el error queda en 9.9 %: pasa, pero roza el criterio y el test se vuelve frágil
    "alto_m": 22.0e-3,
    "fuente": "puntual",
    "capas": [{"material": "agua", "espesor_m": 22.0e-3}],
}

PLANO = {
    "frecuencia_hz": 100.0e3,
    "ppw": 10,
    "alto_m": 15.0e-3,
    "fuente": "plano",
    "capas": [
        {"material": "agua", "espesor_m": 22.5e-3},
        {"material": "hueso", "espesor_m": 7.5e-3},
        {"material": "agua", "espesor_m": 22.5e-3},
    ],
}


# --------------------------------------------------------------------------- validación


def test_material_desconocido():
    # dos capas: con una sola salta primero el mínimo de capas del modo plano, y lo que este test
    # fija es el mensaje del MATERIAL desconocido
    with pytest.raises(ValueError) as exc:
        volume3d.validar_volumen(dict(PLANO, capas=[
            {"material": "agua", "espesor_m": 20e-3},
            {"material": "kriptonita", "espesor_m": 1e-3},
        ]))
    assert "kriptonita" in str(exc.value)
    assert "agua" in str(exc.value)


def test_fuente_desconocida():
    with pytest.raises(ValueError):
        volume3d.validar_volumen(dict(PLANO, fuente="laser"))


def test_incidencia_oblicua_se_rechaza_con_la_razon():
    with pytest.raises(ValueError) as exc:
        volume3d.validar_volumen(dict(PLANO, angulo_deg=15.0))
    assert "incidencia normal" in str(exc.value)


def test_puntual_exige_medio_homogeneo():
    mezcla = dict(PUNTUAL, capas=[{"material": "agua", "espesor_m": 10e-3},
                                  {"material": "hueso", "espesor_m": 5e-3}])
    with pytest.raises(ValueError) as exc:
        volume3d.validar_volumen(mezcla)
    assert "homogéneo" in str(exc.value)


def test_malla_gigante_se_rechaza_con_que_bajar():
    grande = dict(PUNTUAL, ppw=20, alto_m=45e-3)
    with pytest.raises(ValueError) as exc:
        volume3d.estimar_volumen(grande)
    mensaje = str(exc.value)
    assert "tope interactivo" in mensaje
    assert "ppw" in mensaje


# --------------------------------------------------------------------------- costo


def test_estimar_no_resuelve_y_da_la_longitud_real():
    c = volume3d.estimar_volumen(PUNTUAL)
    assert c["malla"]["n_nodos"] > 1000
    assert "metricas" not in c
    # con fuente puntual la malla va centrada: la longitud es el rango, no xs[-1]
    assert c["malla"]["longitud_mm"] == pytest.approx(22.0, rel=0.05)
    assert c["tope_nodos"] > 0


# --------------------------------------------------------------------------- aceptación


def test_la_fuente_puntual_reproduce_el_green_3d():
    r = volume3d.resolver_volumen(PUNTUAL)
    m = r["metricas"]
    print(f"\n  error Green {m['error_green_pct']:.1f} % · amplitud {m['amplitud_green']:.3f} · "
          f"{r['malla']['n_nodos']} nodos · {r['tiempos']['solve_s']:.2f} s")
    assert m["error_green_pct"] < 10.0
    assert m["amplitud_green"] == pytest.approx(1.0, abs=0.06)


def test_la_fuente_no_queda_sobre_la_abc():
    """La regresión del 48 %: la fuente tiene que estar en el centro, no en la cara.

    Sin centrar, el nodo de la fuente coincidía con un nodo de la frontera y el error contra el Green
    no bajaba al refinar (medido: 48.7 % a 1.5 λ y 49.0 % a 2.5 λ, o sea que no era discretización).
    """
    from ultrasound_sim.meshing3d import boundary_faces

    p = volume3d.validar_volumen(PUNTUAL)
    xs, _ys, _zs, malla, _ = volume3d._malla(p)
    idx = int(np.argmin(np.linalg.norm(malla.nodes, axis=1)))
    caras = boundary_faces(malla.tets)
    assert idx not in set(caras.ravel().tolist()), "la fuente cayó sobre la frontera"
    centro = np.array([(xs[0] + xs[-1]) / 2.0, 0.0, 0.0])
    assert np.linalg.norm(malla.nodes[idx] - centro) < 0.1 * (xs[-1] - xs[0])


ESFERA = {
    "frecuencia_hz": 100.0e3,
    "ppw": 12,
    "alto_m": 30.0e-3,
    "fuente": "esfera",
    "ka": 1.0,
    "material_esfera": "hueso",
    "capas": [{"material": "agua", "espesor_m": 45.0e-3}],
}


def test_el_modo_esfera_compara_contra_anderson_y_sella():
    """Dispersión penetrable 3D contra la serie de ondas parciales de Anderson (1950)."""
    r = volume3d.resolver_volumen(ESFERA)
    m = r["metricas"]
    print(f"\n  ka={m['ka']}: error {m['error_serie_pct']:.2f} % contra Anderson · "
          f"{m['elementos_en_esfera']} elementos en la esfera · {m['puntos_comparados']} puntos")
    assert m["error_serie_pct"] < 5.0
    assert m["elementos_en_esfera"] > 20, "la malla no resuelve la esfera"
    assert m["puntos_comparados"] > 500
    assert r["sello"]["estado"] == "validado"
    ids = {c["id"] for c in r["sello"]["criterios"]}
    assert ids == {"volumen_esfera_serie"}


def test_el_modo_esfera_muestra_el_campo_dispersado_y_no_el_total():
    """Si se graficara el total, la incidente (|p| ≈ 1) taparía la perturbación.

    El síntoma medido fue un mapa casi uniforme con un rango de 0.3 dB; mostrando el dispersado el
    rango pasa a ser de decenas de dB.
    """
    r = volume3d.resolver_volumen(ESFERA)
    valores = [v for fila in r["rebanadas"]["xz"]["db"] for v in fila]
    assert max(valores) - min(valores) > 15.0, (
        "el mapa tiene rango dinámico de campo total, no de campo dispersado"
    )


def test_el_modo_esfera_requiere_que_la_malla_caiga_dentro():
    """Con ppw bajo y ka chico no entra ningún elemento en la esfera: se rechaza, no se resuelve."""
    chico = dict(ESFERA, ppw=6.0, ka=0.2)
    with pytest.raises(ValueError) as exc:
        volume3d.resolver_volumen(chico)
    assert "dentro de la esfera" in str(exc.value)


def test_el_modo_plano_da_ryt_dentro_del_2_por_ciento():
    r = volume3d.resolver_volumen(PLANO)
    m = r["metricas"]
    print(f"\n  R {m['R']:.4f} vs {m['R_analitico']:.4f} (abs {m['error_R_abs']:.4f}) · "
          f"T {m['T']:.4f} vs {m['T_analitico']:.4f} (abs {m['error_T_abs']:.4f}) · "
          f"energía {m['energia']:.4f}")
    assert m["error_R_abs"] < 0.02
    assert m["error_T_abs"] < 0.02
    assert m["energia"] == pytest.approx(1.0, abs=0.02)


# --------------------------------------------------------------------------- forma y serialización


def test_rebanadas_tienen_la_forma_de_sus_ejes():
    r = volume3d.resolver_volumen(PLANO)
    xy = r["rebanadas"]["xy"]
    assert len(xy["db"]) == len(xy["y_mm"]), "las filas son y"
    assert len(xy["db"][0]) == len(xy["x_mm"]), "las columnas son x"
    xz = r["rebanadas"]["xz"]
    assert len(xz["db"]) == len(xz["z_mm"])
    assert len(xz["db"][0]) == len(xz["x_mm"])
    assert len(r["rebanadas"]["interfaces_mm"]) == 2


def test_el_resultado_es_json_estricto():
    for datos in (PUNTUAL, PLANO):
        texto = json.dumps(volume3d.resolver_volumen(datos), allow_nan=False)
        assert len(texto) > 500
        json.loads(texto)


def test_las_rebanadas_no_son_uniformes():
    """Si el campo saliera constante el mapa no diría nada: tiene que haber estructura."""
    r = volume3d.resolver_volumen(PUNTUAL)
    valores = [v for fila in r["rebanadas"]["xy"]["db"] for v in fila]
    assert max(valores) - min(valores) > 6.0, "la rebanada no muestra estructura"
