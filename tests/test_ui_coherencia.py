"""Coherencia entre el cliente del dashboard y el servidor.

Son tests baratos y atajan una clase de error que ningún test de física ve: **que una vista exista en
la navegación pero no tenga su panel**, o que el cliente espere un campo que el servidor no manda. Los
dos casos dejan la pantalla en blanco sin una sola excepción, que es la peor forma de fallar.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

ESTATICOS = Path(__file__).resolve().parents[1] / "src" / "ultrasound_sim" / "dashboard" / "static"


@pytest.fixture(scope="module")
def html() -> str:
    return (ESTATICOS / "index.html").read_text(encoding="utf-8")


def test_cada_boton_de_vista_tiene_su_panel(html):
    """Un `data-vista` sin `vista<X>Main` hace que el click no muestre nada: el usuario ve la pestaña
    anterior y cree que la app se colgó (pasó con la vista de criterios)."""
    vistas = re.findall(r'data-vista="([a-z0-9_]+)"', html)
    ids = set(re.findall(r'id="([A-Za-z0-9_]+)"', html))
    assert vistas, "no se encontró ningún botón de vista"
    faltan = []
    for v in vistas:
        candidatos = {f"vista{v.capitalize()}Main", f"main{v.capitalize()}"}
        if not (candidatos & ids):
            faltan.append((v, sorted(candidatos)))
    print(f"\n  vistas: {vistas}")
    assert not faltan, f"vistas sin panel: {faltan}"


def test_cada_script_referenciado_existe(html):
    """Un `<script src>` que no está deja funciones sin definir y el error aparece lejos, al hacer
    click. Se verifica que el archivo exista en disco."""
    fuentes = re.findall(r'<script src="/static/([^"]+)"', html)
    assert fuentes, "la página no carga ningún script"
    faltan = [f for f in fuentes if not (ESTATICOS / f).exists()]
    print(f"\n  scripts: {fuentes}")
    assert not faltan, f"scripts referenciados que no existen: {faltan}"


def test_la_vista_de_criterios_lee_el_registro_del_servidor():
    """El contrato con `/api/validacion`: la ruta devuelve `{criterios: [...]}` **sin** campo `ok`, y
    el cliente tiene que aceptarlo así. Cuando el cliente exigía `ok`, la tabla salía vacía sin
    explicar nada."""
    js = (ESTATICOS / "criterios.js").read_text(encoding="utf-8")
    assert "/api/validacion" in js
    assert "r.criterios" in js
    # el guardia tiene que ser sobre `criterios`, no sobre `ok`
    assert "if (!r || !Array.isArray(r.criterios))" in js


def test_el_registro_no_declara_una_unidad_que_no_mide(html):
    """D49: el criterio de la esfera elástica decía 'L2 relativo' y la métrica que recibe es el
    máximo. La unidad tiene que decir lo que se mide."""
    from ultrasound_sim.dashboard import validacion

    criterio = validacion.CRITERIOS["volumen_esfera_elastica"]
    print(f"\n  unidad declarada: {criterio.unidad!r}")
    assert criterio.unidad != "L2 relativo", "declara L2 pero se le pasa un máximo"


def test_las_vistas_guardan_su_resultado_para_no_repetir_el_solve():
    """La persistencia es lo que permite volver a una pestaña sin pagar otro solve (la esfera precisa
    son ~30 s). Se comprueba que cada vista cara tenga su `vistaGuardar` y su `vistaCargar`."""
    caras = {"escena": "escena.js", "volumen": "volumen.js",
             "cochlea": "cochlea.js", "esfera": "esfera.js"}
    for vista, archivo in caras.items():
        js = (ESTATICOS / archivo).read_text(encoding="utf-8")
        assert f"vistaGuardar('{vista}'" in js, f"{archivo} no guarda"
        assert f"vistaCargar('{vista}')" in js, f"{archivo} no restaura"
