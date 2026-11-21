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
    # los botones deshabilitados (FEM en desarrollo, y las de «para el futuro») no tienen panel: no
    # se pueden abrir. Se excluyen del chequeo.
    vistas = re.findall(r'data-vista="([a-z0-9_]+)"(?![^>]*\bdisabled)', html)
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


def test_el_nav_muestra_el_problema_y_deja_lo_demas_para_el_futuro():
    """D54: el tablero pasa a un solo problema (BEM de fuentes, con FEM en desarrollo) y las otras
    cuatro vistas quedan marcadas «para el futuro». La que no está validada se muestra DESHABILITADA y
    con la etiqueta, no escondida ni habilitada."""
    html = (ESTATICOS / "index.html").read_text(encoding="utf-8")
    assert 'data-vista="fem"' in html and "en-desarrollo" in html
    assert 'class="proximamente"' in html
    for vista in ("cochlea", "barridos", "esfera", "criterios"):
        assert f'data-vista="{vista}" disabled' in html, f"{vista} debería estar deshabilitada"
    # el botón de FEM tiene que decir por qué
    assert "etiqueta-dev" in html


def test_la_vista_de_cadena_tiene_panel_main_y_script():
    """La cadena piel→cóclea es la vista donde se acoplan sólidos y fluidos: tiene que estar en el
    conmutador, con su panel y su script, o el botón abre una pestaña vacía."""
    html = (ESTATICOS / "index.html").read_text(encoding="utf-8")
    js_escena = (ESTATICOS / "escena.js").read_text(encoding="utf-8")
    assert 'data-vista="cadena"' in html
    assert 'id="vistaCadena"' in html and 'id="vistaCadenaMain"' in html
    assert "cadena: { panel: 'vistaCadena', main: 'vistaCadenaMain' }" in js_escena
    assert '<script src="/static/cadena.js">' in html
    # el aviso del explorador de campo ("el pico real está en…") no corresponde acá
    assert "vista-cadena" in (ESTATICOS / "style.css").read_text(encoding="utf-8")


def test_se_puede_cargar_un_stl_y_ver_su_seccion():
    """U9: la geometría propia entra como la sección del STL (el modelo es 2D)."""
    html = (ESTATICOS / "index.html").read_text(encoding="utf-8")
    js = (ESTATICOS / "cadena.js").read_text(encoding="utf-8")
    assert 'id="stlArchivo"' in html and 'accept=".stl' in html
    assert 'id="stlMaterial"' in html and 'id="stlPlano"' in html
    assert 'id="stlEjemplo"' in html, "hace falta un ejemplo para probar sin tener un archivo"
    assert '/api/jobs/stl' in js
    assert 'type: \'path\'' in js, "el contorno de la sección tiene que dibujarse"


def test_se_puede_cargar_el_arreglo_y_el_transductor_desde_el_panel():
    """U10: cargar un arreglo planeado y el catálogo de transductores dejan de ser 32 filas a mano."""
    html = (ESTATICOS / "index.html").read_text(encoding="utf-8")
    js = (ESTATICOS / "arreglos.js").read_text(encoding="utf-8")
    app = (ESTATICOS / "app.js").read_text(encoding="utf-8")
    assert 'id="selTransductor"' in html
    assert 'id="btnCargarArreglo"' in html and 'id="btnExportarArreglo"' in html
    assert 'id="archivoArreglo"' in html and 'accept=".csv,.json' in html
    assert '<script src="/static/arreglos.js">' in html
    # el puente con el estado privado de app.js
    assert "window.campoAplicarFuentes" in app
    assert "window.campoLeerFuentes" in app
    assert "window.campoAplicarTransductor" in app
    # y que un archivo rechazado NO toque lo que está en pantalla
    assert "no se toca lo que está en pantalla" in js


def test_la_estructura_se_arrastra_en_el_mapa():
    """El control de posicionamiento de la estructura es parte del pedido: tiene que haber un shape
    arrastrable, el listener del relayout registrado después del react (antes el div no era de Plotly y
    el handler no se registraba) y una señal visual de que se puede arrastrar."""
    js = (ESTATICOS / "cadena.js").read_text(encoding="utf-8")
    assert "edits: { shapePosition: e.finita }" in js
    assert "name: 'estructura'" in js
    assert "plotly_relayout" in js
    # registrado en pintar(), no en init()
    assert "gdMapa.__cadHandler" in js
    assert "style.cursor" in js and "grab" in js
    # las líneas de interfaz se recortan al alto de la pieza
    assert "const yTop = e.finita" in js


def test_la_cadena_guarda_su_corrida_con_la_clave_pelada():
    """`vistaGuardar` antepone 'usim.vista.' y envuelve en {fecha, resultado}: pasarle la clave con el
    prefijo ya puesto guardaba en 'usim.vista.usim.vista.cadena' y la restauración leía otro nivel, así
    que la vista volvía vacía sin ningún error."""
    js = (ESTATICOS / "cadena.js").read_text(encoding="utf-8")
    assert "const CLAVE = 'cadena';" in js
    assert "usim.vista.cadena'" not in js


def test_todo_lo_que_cambia_las_fuentes_del_arreglo_las_invalida():
    """Cambiar n o la apertura movía las posiciones; el FOCO, la FRECUENCIA y el MATERIAL cambian la
    FASE de cada elemento (`−k·|F − y|` con `k = 2π·f/c`). Mirar sólo `n` y `d` dejaba fases viejas
    guardadas y al pasar a «fuentes» se simulaba un arreglo distinto, sin aviso."""
    js = (ESTATICOS / "app.js").read_text(encoding="utf-8")
    assert "CLAVES_QUE_REGENERAN" in js
    for clave in ("'n'", "'d'", "'fx'", "'fy'", "'freq'", "'material'"):
        assert clave in js.split("CLAVES_QUE_REGENERAN = new Set(")[1].split(")")[0], (
            f"{clave} tiene que invalidar las fuentes guardadas"
        )
    # en «fuentes» se regenera al instante: si no, seguiría simulando fases viejas
    assert "generarDesdeArreglo();\n        S.fuentesMovidas = false;" in js


def test_los_valores_de_la_tabla_de_fuentes_se_pueden_leer_y_ubicar():
    """La tabla tenía 253 px para cinco columnas numéricas: el input quedaba en 52 px y `-29.06` se
    recortaba. Ahora el panel es más ancho, los spinners nativos de `type=number` no roban el campo, el
    número se alinea a la derecha, hay columna de fila y al editar se resalta el asa en el mapa."""
    css = (ESTATICOS / "style.css").read_text(encoding="utf-8")
    js = (ESTATICOS / "app.js").read_text(encoding="utf-8")
    html = (ESTATICOS / "index.html").read_text(encoding="utf-8")
    assert "width:380px" in css, "el panel tiene que dejar entrar los valores"
    assert "-webkit-inner-spin-button" in css, "los spinners nativos recortan el número"
    assert "text-align:right" in css
    assert 'class="c-idx"' in html and 'th-idx' in html
    assert "window.resaltarFuente" in js
    assert "tr.classList.add('editando')" in js


def test_la_tabla_de_fuentes_scrollea_sola_y_muestra_cuantas_hay():
    """Con 32 fuentes la tabla media 911 px, mas alta que la ventana: habia que arrastrar todo el panel
    para llegar a la fuente 30. Ahora tiene scroll propio y encabezado fijo, y un contador que dice
    cuantas hay y cuanto miden de apertura."""
    html = (ESTATICOS / "index.html").read_text(encoding="utf-8")
    css = (ESTATICOS / "style.css").read_text(encoding="utf-8")
    js = (ESTATICOS / "app.js").read_text(encoding="utf-8")
    assert 'class="tabla-scroll"' in html
    assert 'id="srcCuenta"' in html
    assert ".tabla-scroll" in css and "position: sticky" in css
    assert "$('srcCuenta')" in js and "apertura" in js


def test_la_frecuencia_se_sincroniza_entre_vistas():
    """Campo arrancaba en 250 kHz y escena/3D en 100 kHz: comparar obligaba a igualarlas a mano."""
    js = (ESTATICOS / "app.js").read_text(encoding="utf-8")
    assert "syncFrecuencia" in js
    for archivo in ("escena.js", "volumen.js"):
        fuente = (ESTATICOS / archivo).read_text(encoding="utf-8")
        assert "syncFrecuencia" in fuente, f"{archivo} no propaga la frecuencia"


def test_las_fuentes_movidasNo_se_descartan_al_volver_al_arreglo():
    """Bug reportado: "muevo las fuentes, vuelvo al arreglo y no quedan donde las moví".

    Eran tres fugas: (a) tocar n o d hacía `S.fuentes = null` en silencio, (b) el tab decidía qué se
    enviaba, así que en "arreglo" el campo se recalculaba ignorando lo movido, y (c) las asas y los
    marcadores se dibujaban con las posiciones del arreglo aunque el campo viniera de las fuentes.
    """
    js = (ESTATICOS / "app.js").read_text(encoding="utf-8")
    # (a) no se borra: se marca para regenerar al volver a «fuentes»
    assert "S.fuentes = null;   // el arreglo manda de nuevo" not in js
    assert "S.regenerarAlVolver" in js
    # (b) el envío depende de si hay fuentes movidas, no del tab
    assert "S.fuentes && (S.modo === 'manual' || S.fuentesMovidas)" in js
    # (c) asas y marcadores siguen la misma regla
    assert js.count("(S.modo === 'manual' || S.fuentesMovidas)") >= 3
    # y el aviso no vive en `#estado`, que la respuesta del campo reescribe
    assert 'id="avisoUX"' in (ESTATICOS / "index.html").read_text(encoding="utf-8")


def test_las_vistas_guardan_su_resultado_para_no_repetir_el_solve():
    """La persistencia es lo que permite volver a una pestaña sin pagar otro solve (la esfera precisa
    son ~30 s). Se comprueba que cada vista cara tenga su `vistaGuardar` y su `vistaCargar`."""
    caras = {"escena": "escena.js", "volumen": "volumen.js",
             "cochlea": "cochlea.js", "esfera": "esfera.js"}
    for vista, archivo in caras.items():
        js = (ESTATICOS / archivo).read_text(encoding="utf-8")
        assert f"vistaGuardar('{vista}'" in js, f"{archivo} no guarda"
        assert f"vistaCargar('{vista}')" in js, f"{archivo} no restaura"
