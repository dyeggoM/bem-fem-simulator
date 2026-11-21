"""Tests de U6: trazabilidad del sello.

Lo que se fija acá:
1. Los tres estados existen y son excluyentes: ``validado``, ``fuera``, ``sin_referencia``.
2. ``sin_referencia`` **no** es un aprobado por omisión (el error clásico de este tipo de código es
   que un `None` termine contando como "cumple").
3. El resumen de varios criterios no promedia: si uno falla, el conjunto falla.
4. Cada criterio trae su fuente (ADR/fase), así el número no queda huérfano.
5. **Cada corrida de la UI trae su sello**, y ningún archivo JS del cliente tiene referencias ni
   umbrales: es la condición que hace que el sello sea auditable.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from ultrasound_sim.dashboard import validacion

ESTATICOS = Path(__file__).resolve().parents[1] / "src" / "ultrasound_sim" / "dashboard" / "static"


# --------------------------------------------------------------------------- estados


def test_dentro_de_tolerancia_es_validado():
    s = validacion.sello("escena_rt_r", 0.005)
    assert s["estado"] == "validado"
    assert s["margen"] == pytest.approx(0.005)


def test_fuera_de_tolerancia_no_es_validado():
    s = validacion.sello("escena_rt_r", 0.05)
    assert s["estado"] == "fuera"
    assert "FUERA" in s["texto"]


def test_el_limite_cuenta_como_validado():
    """`≤ tolerancia` y no `<`: en el borde exacto el criterio es que entra."""
    s = validacion.sello("escena_rt_r", 0.02)
    assert s["estado"] == "validado"


def test_sin_medido_no_es_un_aprobado_por_omision():
    """El error caro: un `None` que termina contando como que cumple."""
    s = validacion.sello("escena_rt_r", None)
    assert s["estado"] == "sin_referencia"
    assert s["estado"] != "validado"
    assert s["medido"] is None
    assert "SIN REFERENCIA" in s["texto"]


def test_referencia_relativa_usa_el_desvio_relativo():
    s = validacion.sello("cochlea_t_interior", 0.53187 * 1.005)  # +0.5 %
    assert s["estado"] == "validado"
    assert s["margen"] == pytest.approx(0.005, rel=1e-6)
    s2 = validacion.sello("cochlea_t_interior", 0.53187 * 1.02)  # +2 %
    assert s2["estado"] == "fuera"


# --------------------------------------------------------------------------- conjunto


def test_el_conjunto_no_promedia():
    """Un criterio que falla arrastra al conjunto, aunque el otro esté holgado."""
    r = validacion.sello_conjunto({"escena_rt_r": 0.001, "escena_rt_t": 0.5})
    assert r["estado"] == "fuera"
    assert len(r["criterios"]) == 2


def test_conjunto_sin_referencia_no_es_validado():
    r = validacion.sello_conjunto({"escena_rt_r": 0.001, "cochlea_t_interior": None})
    assert r["estado"] == "sin_referencia"


def test_conjunto_todo_dentro_es_validado():
    r = validacion.sello_conjunto({"escena_rt_r": 0.001, "escena_rt_t": 0.002})
    assert r["estado"] == "validado"
    assert "VALIDADO" in r["texto"]


def test_criterio_desconocido_falla_con_la_lista():
    with pytest.raises(KeyError) as exc:
        validacion.sello("inventado", 0.0)
    assert "escena_rt_r" in str(exc.value)


# --------------------------------------------------------------------------- catálogo


def test_cada_criterio_trae_su_fuente():
    for c in validacion.catalogo():
        assert c["fuente"], f"{c['id']} no dice de dónde sale"
        assert c["descripcion"]
        assert c["modo"] in ("abs", "rel")
        assert c["tolerancia"] > 0


def test_los_criterios_cubren_las_vistas():
    ids = {c["id"] for c in validacion.catalogo()}
    esperados = {"escena_rt_r", "escena_rt_t", "dispersor_fem_serie", "dispersor_bem_serie",
                 "cochlea_t_interior", "cochlea_balance", "volumen_green", "volumen_rt_r",
                 "volumen_rt_t", "barrido_d30"}
    assert esperados <= ids


# --------------------------------------------------------------------------- las corridas lo traen


def test_la_escena_trae_su_sello():
    from ultrasound_sim.dashboard import scene

    r = scene.resolver_escena({
        "frecuencia_hz": 100e3, "ppw": 12, "alto_m": 20e-3,
        "capas": [{"material": "agua", "espesor_m": 20e-3},
                  {"material": "hueso", "espesor_m": 10e-3},
                  {"material": "agua", "espesor_m": 20e-3}],
    })
    assert r["sello"]["estado"] == "validado"
    assert "Fuente:" in r["sello"]["texto"]


def test_el_dispersor_trae_su_sello():
    from ultrasound_sim.dashboard import scatterer

    r = scatterer.resolver_dispersor({"ka": 1.0})
    assert r["sello"]["estado"] == "validado"
    ids = {c["id"] for c in r["sello"]["criterios"]}
    assert ids == {"dispersor_fem_serie", "dispersor_bem_serie"}


def test_la_coclea_fuera_del_caso_base_da_sin_referencia_en_t_interior():
    """Con espesor 3 mm no hay referencia de t_interior: el sello tiene que decirlo."""
    from ultrasound_sim.dashboard import cochlea_field

    r = cochlea_field.resolver_cochlea({"thickness_mm": 3.0})
    por_id = {c["id"]: c for c in r["sello"]["criterios"]}
    assert por_id["cochlea_t_interior"]["estado"] == "sin_referencia"
    assert r["sello"]["estado"] in ("sin_referencia", "fuera")


# --------------------------------------------------------------------------- el cliente no decide

#: Números de referencia que NO pueden aparecer en el JS: si están, la validación quedó partida entre
#: servidor y cliente, que es exactamente lo que U6 vino a cerrar.
PROHIBIDOS = ("0.53187", "CRITERIO_AMPLITUD", "REF_T_INTERIOR", "CRITERIO_GREEN", "CRITERIO_DISP")


@pytest.mark.parametrize("archivo", ["escena.js", "cochlea.js", "dispersor.js", "volumen.js",
                                     "barridos.js"])
def test_el_js_del_cliente_no_tiene_la_referencia_ni_el_umbral(archivo):
    texto = (ESTATICOS / archivo).read_text(encoding="utf-8")
    for prohibido in PROHIBIDOS:
        assert prohibido not in texto, f"{archivo} todavía tiene {prohibido}: el sello no es del servidor"


@pytest.mark.parametrize("archivo", ["escena.js", "cochlea.js", "dispersor.js", "volumen.js",
                                     "barridos.js"])
def test_todas_las_vistas_usan_el_sello_del_servidor(archivo):
    texto = (ESTATICOS / archivo).read_text(encoding="utf-8")
    assert "mostrarSello(" in texto, f"{archivo} no pinta el sello del servidor"
    # y no decide por su cuenta comparando contra un número
    assert not re.search(r"<\s*0\.\d+\s*\)", texto.replace(" ", " ")), (
        f"{archivo} parece comparar contra un umbral propio"
    )
