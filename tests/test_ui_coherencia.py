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


def test_el_tour_existe_tiene_pasos_y_se_puede_reabrir(html):
    """El tour es la única documentación que el usuario ve sin buscarla. Tiene que existir el archivo,
    estar cargado en la página, tener pasos y poder reabrirse (si sólo se viera la primera vez, el que
    lo cerró sin querer no lo encuentra más)."""
    js = (ESTATICOS / "tour.js").read_text(encoding="utf-8")
    assert 'src="/static/tour.js"' in html
    assert 'id="btnTour"' in html, "sin botón no se puede repetir el tour"
    assert "TOUR_PASOS" in js and js.count("titulo:") >= 5
    assert "tourAbrir" in js and "localStorage" in js
    print(f"\n  pasos del tour: {js.count('titulo:')}")


def test_los_tooltips_no_usan_el_title_nativo():
    """El `title` nativo se dibuja con los colores del sistema y en tema oscuro no se lee: por eso el
    texto se mueve a `data-tip` y se borra el `title`. Si alguien vuelve a dejar el `title`, salen los
    dos globos y el ilegible es el primero."""
    js = (ESTATICOS / "app.js").read_text(encoding="utf-8")
    css = (ESTATICOS / "style.css").read_text(encoding="utf-8")
    assert "activarTips" in js
    assert "removeAttribute('title')" in js
    assert ".globo-tip" in css
    # contraste: fondo claro y texto oscuro, no al revés
    assert "#f2f5f9" in css


def test_los_barridos_muestran_un_spinner_mientras_corren():
    """Un barrido de minutos con el panel vacío parece que no hizo nada. Cada barrido tiene que
    encender el spinner al lanzar y apagarlo al terminar, en éxito y en error."""
    js = (ESTATICOS / "app.js").read_text(encoding="utf-8")
    assert "window.cargando" in js and "terminarCargando" in js
    for archivo, contenedor in (("barridos.js", "graficoApertura"),
                                ("barridos.js", "graficoBarrido"),
                                ("esfera.js", "barridoEsf")):
        fuente = (ESTATICOS / archivo).read_text(encoding="utf-8")
        assert f"cargando('{contenedor}'" in fuente, f"{archivo} no enciende el spinner de {contenedor}"
        assert f"terminarCargando('{contenedor}')" in fuente, f"{archivo} no lo apaga"


def test_el_arrastre_de_fuentes_no_lo_come_el_zoom_de_plotly():
    """El síntoma que reportó el usuario: arrastrar hacía zoom. La causa era el `dragmode` de Plotly,
    que se queda con el gesto. Con el modo "mover" el dragmode tiene que estar apagado, y el arrastre
    tiene que existir también en modo arreglo (antes sólo había asas en modo "fuentes")."""
    js = (ESTATICOS / "app.js").read_text(encoding="utf-8")
    assert "dragmode: S.mover ? false : 'zoom'" in js
    assert "posicionesDelArreglo" in js, "en modo arreglo no hay nada que arrastrar"
    assert "btnMover" in js and "btnZoom" in js


def test_el_nav_no_necesita_scroll_horizontal():
    """Diez botones en una fila con `flex:1` desbordaban y había que scrollear para llegar a los
    últimos. La grilla de dos columnas los muestra todos."""
    css = (ESTATICOS / "style.css").read_text(encoding="utf-8")
    assert ".vistas { display:grid" in css and "grid-template-columns:1fr 1fr" in css
    assert ".acciones-header" in css, "tour y exportar van fuera del nav"


def test_la_accion_principal_queda_fija_y_visible():
    """El botón que corre la simulación estaba al final de una columna larga: había que scrollear
    para encontrarlo. Ahora la fila con la acción se mueve a una barra sticky al pie del panel."""
    js = (ESTATICOS / "app.js").read_text(encoding="utf-8")
    css = (ESTATICOS / "style.css").read_text(encoding="utf-8")
    html = (ESTATICOS / "index.html").read_text(encoding="utf-8")
    assert 'id="accionesVista"' in html
    assert "ubicarAcciones" in js
    assert "position: sticky" in css and ".acciones-vista" in css
    # se mueve el nodo, no se recrea: recrear perdería los listeners
    assert "appendChild(fila)" in js


def test_los_hoverlabels_de_plotly_son_claros():
    """El usuario reportó "tooltips muy oscuros": no eran sólo los `?`, también los hoverlabels de
    Plotly, que venían con fondo #1b2129."""
    js = (ESTATICOS / "app.js").read_text(encoding="utf-8")
    assert "hoverlabel" in js
    assert "bgcolor: '#f2f5f9'" in js
    assert "bgcolor: '#1b2129'" not in js


def test_los_estaticos_no_se_cachean():
    """`no-cache` permite guardar y sólo obliga a revalidar; sin ETag el navegador decide solo. Eso
    hacía que un arreglo de UI "no surtiera efecto" hasta recargar a mano."""
    fuente = (Path(__file__).resolve().parents[1] / "src" / "ultrasound_sim" / "dashboard"
              / "server.py").read_text(encoding="utf-8")
    assert "no-store" in fuente


def test_las_vistas_guardan_su_resultado_para_no_repetir_el_solve():
    """La persistencia es lo que permite volver a una pestaña sin pagar otro solve (la esfera precisa
    son ~30 s). Se comprueba que cada vista cara tenga su `vistaGuardar` y su `vistaCargar`."""
    caras = {"escena": "escena.js", "volumen": "volumen.js",
             "cochlea": "cochlea.js", "esfera": "esfera.js"}
    for vista, archivo in caras.items():
        js = (ESTATICOS / archivo).read_text(encoding="utf-8")
        assert f"vistaGuardar('{vista}'" in js, f"{archivo} no guarda"
        assert f"vistaCargar('{vista}')" in js, f"{archivo} no restaura"
