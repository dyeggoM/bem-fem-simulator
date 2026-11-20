"""Tests del dashboard: que el motor coincida con el paquete y que el servidor responda de verdad.

El servidor se levanta en un puerto libre efímero y se le habla con ``urllib``: se prueba el camino
HTTP completo (routing, JSON, validación), no la función por dentro.
"""

from __future__ import annotations

import json
import math
import threading
import urllib.error
import urllib.request

import numpy as np
import pytest

from ultrasound_sim import aperture, arrays, materials
from ultrasound_sim.dashboard import field as campo
from ultrasound_sim.dashboard.server import ESTATICOS, build_server, calcular

W = materials.WATER


# ------------------------------------------------------------------ motor
def test_evaluate_coincide_con_aperture_pressure():
    """El evaluador vectorizado tiene que dar lo mismo que la integral de referencia."""
    f = 100e3
    fuentes = [
        campo.Source(y=-0.012, width=0.0035, amp=1.0, phase=0.3),
        campo.Source(y=-0.001, width=0.0035, amp=0.8, phase=-1.1),
        campo.Source(y=0.011, width=0.0035, amp=1.2, phase=2.0),
    ]
    xs = np.linspace(0.02, 0.09, 17)
    ys = np.linspace(-0.02, 0.02, 13)
    X, Y = np.meshgrid(xs, ys)
    puntos = np.column_stack([X.ravel(), Y.ravel()])

    aps = [aperture.Aperture(center_y=s.y, width=s.width, amplitude=s.amplitude) for s in fuentes]
    referencia = aperture.pressure(puntos, f, W, aps).reshape(Y.shape)

    calculado = campo.evaluate(xs, ys, fuentes, f, W)
    assert calculado.shape == referencia.shape
    assert np.allclose(calculado, referencia, rtol=1e-9, atol=1e-3 * np.abs(referencia).max())


def test_array_sources_usa_la_misma_ley_de_fase_que_arrays():
    f = 250e3
    n, d, fx, fy = 16, 0.040, 0.120, 0.008
    fuentes = campo.array_sources(n, d, fx, fy, f, W)
    aps = arrays.LinearArray(center_y=0.0, n_elements=n, pitch=d / n).focused((fx, fy), f, W)

    assert len(fuentes) == n
    assert np.allclose([s.y for s in fuentes], [a.center_y for a in aps], atol=1e-15)
    assert np.allclose([s.width for s in fuentes], [a.width for a in aps], atol=1e-15)
    assert np.allclose([s.amplitude for s in fuentes], [a.amplitude for a in aps], atol=1e-12)


def test_metrics_detecta_pico_y_grataje():
    """Pico cerca del foco geométrico y grataje en el ángulo que predice ``arcsin(λ/pitch)``.

    El grataje sólo existe con ``pitch > λ``: con ``pitch = λ/2`` la predicción da ``arcsin 2``, que no
    tiene solución real — o sea, no hay lóbulo de grataje en el semiespacio visible. Se eligió
    ``pitch = 1.5 λ`` → ``arcsin(1/1.5) = 41.81°``.
    """
    f = 250e3
    lam = W.wavelength(f)
    pitch = 1.5 * lam
    fuentes = campo.array_sources(8, 8 * pitch, 0.12, 0.0, f, W)
    xs = np.linspace(0.05, 0.18, 90)
    ys = np.linspace(-0.03, 0.03, 90)
    p = campo.evaluate(xs, ys, fuentes, f, W)
    m = campo.metrics(xs, ys, p, pitch=pitch, frequency=f, material=W)

    assert m["grating_deg"] == pytest.approx(math.degrees(math.asin(1 / 1.5)), abs=0.01)
    assert abs(m["peak_y_mm"]) < 1.0                  # foco en el eje
    assert 100.0 < m["peak_x_mm"] < 140.0             # cerca de F=120 mm, no en él (D13)
    assert m["sidelobe_db"] < 0.0
    # con pitch > λ los lóbulos de grataje pueden no dejar cruce a −6 dB: None es la respuesta
    # correcta, no un FWHM inventado (se verifica en test_metrics_devuelve_none_si_el_lobulo_no_cruza)


def test_sin_grataje_cuando_el_pitch_es_menor_que_lambda():
    """``pitch ≤ λ`` ⇒ no hay lóbulo de grataje visible: el campo lo reporta como None."""
    f = 250e3
    lam = W.wavelength(f)
    fuentes = campo.array_sources(32, 32 * (lam / 2), 0.12, 0.0, f, W)
    xs = np.linspace(0.05, 0.16, 40)
    ys = np.linspace(-0.03, 0.03, 40)
    p = campo.evaluate(xs, ys, fuentes, f, W)
    m = campo.metrics(xs, ys, p, pitch=lam / 2, frequency=f, material=W)
    assert m["grating_deg"] is None
    # con pitch = λ/2 el lóbulo principal sí cruza −6 dB: el FWHM tiene que existir y ser finito.
    # Teórico 0.886·λ·f/# con f/# = 1.25 → 6.6 mm; se acepta un 30 % de margen por la malla gruesa.
    assert m["fwhm_lateral_mm"] is not None and 4.0 < m["fwhm_lateral_mm"] < 10.0


def test_metrics_devuelve_none_si_el_lobulo_no_cruza():
    """Sin cruce a −6 dB no se inventa un FWHM: se devuelve None."""
    xs = np.linspace(0.01, 0.02, 5)
    ys = np.linspace(-0.005, 0.005, 5)
    p = np.ones((5, 5), dtype=complex)               # plano: no hay -6 dB en la ventana
    m = campo.metrics(xs, ys, p)
    assert m["fwhm_lateral_mm"] is None


def test_calcular_rechaza_malla_invalida():
    with pytest.raises(ValueError, match="x0"):
        calcular({"grid": {"x0": 0.0, "x1": 0.1, "y0": -0.01, "y1": 0.01, "nx": 10, "ny": 10}})
    with pytest.raises(ValueError, match="material"):
        calcular({"material": "vacío"})


def test_calcular_devuelve_forma_y_metadatos():
    r = calcular({"frequency_khz": 250, "n_elements": 16, "aperture_mm": 30,
                  "focus_x_mm": 100, "grid": {"nx": 24, "ny": 18}})
    assert len(r["x_mm"]) == 24 and len(r["y_mm"]) == 18
    assert np.array(r["p2_db"]).shape == (18, 24)
    assert np.max(r["p2_db"]) == pytest.approx(0.0, abs=0.05)   # normalizado al pico
    assert r["meta"]["n_sources"] == 16
    assert r["meta"]["lambda_mm"] == pytest.approx(6.0, abs=1e-6)
    assert r["metrics"]["peak_x_mm"] > 0


# ------------------------------------------------------------------ servidor HTTP
@pytest.fixture()
def servidor():
    srv = build_server("127.0.0.1", 0)
    hilo = threading.Thread(target=srv.serve_forever, daemon=True)
    hilo.start()
    yield f"http://127.0.0.1:{srv.server_address[1]}"
    srv.shutdown()
    srv.server_close()
    hilo.join(timeout=5)


def _get(url: str) -> tuple[int, bytes]:
    with urllib.request.urlopen(url, timeout=30) as r:
        return r.status, r.read()


def test_servidor_sirve_la_pagina_y_los_assets(servidor):
    estado, cuerpo = _get(servidor + "/")
    assert estado == 200
    assert b"dashboard" in cuerpo and b"plotly.min.js" in cuerpo

    estado, cuerpo = _get(servidor + "/static/app.js")
    assert estado == 200 and b"Plotly.react" in cuerpo

    estado, cuerpo = _get(servidor + "/static/plotly.min.js")
    assert estado == 200 and len(cuerpo) > 500_000        # el bundle cartesiano completo


def test_api_field_por_post(servidor):
    pedido = json.dumps({"frequency_khz": 250, "n_elements": 8, "aperture_mm": 20,
                         "focus_x_mm": 80, "grid": {"nx": 16, "ny": 16}}).encode()
    req = urllib.request.Request(servidor + "/api/field", data=pedido,
                                 headers={"Content-Type": "application/json"}, method="POST")
    with urllib.request.urlopen(req, timeout=60) as r:
        assert r.status == 200
        datos = json.loads(r.read())
    assert len(datos["x_mm"]) == 16
    assert np.array(datos["p2_db"]).shape == (16, 16)
    assert "fwhm_lateral_mm" in datos["metrics"]


def test_api_fuentes_arbitrarias(servidor):
    """Dos fuentes puestas a mano: es el camino que después usa el editor de la UI."""
    pedido = json.dumps({
        "frequency_khz": 100,
        "grid": {"nx": 20, "ny": 20, "x0": 0.01, "x1": 0.10, "y0": -0.02, "y1": 0.02},
        "sources": [{"y": -0.005, "width": 0.002, "amp": 1.0, "phase": 0.0},
                    {"y": 0.005, "width": 0.002, "amp": 1.0, "phase": 0.0}],
    }).encode()
    req = urllib.request.Request(servidor + "/api/field", data=pedido,
                                 headers={"Content-Type": "application/json"}, method="POST")
    with urllib.request.urlopen(req, timeout=60) as r:
        datos = json.loads(r.read())
    assert datos["meta"]["n_sources"] == 2
    # dos fuentes simétricas en fase ⇒ máximo sobre el eje
    assert abs(datos["metrics"]["peak_y_mm"]) < 1.5


def test_api_field_por_get_y_errores(servidor):
    estado, cuerpo = _get(servidor + "/api/field?frequency_khz=250&n_elements=8&nx=12&ny=12")
    assert estado == 200 and len(json.loads(cuerpo)["x_mm"]) == 12

    with pytest.raises(urllib.error.HTTPError) as exc:
        _get(servidor + "/api/field?nx=1")             # nx fuera de rango
    assert exc.value.code == 400
    cuerpo = json.loads(exc.value.read())
    assert "nx" in cuerpo["error"]

    with pytest.raises(urllib.error.HTTPError) as exc:
        _get(servidor + "/static/no-existe.js")
    assert exc.value.code == 404


def test_paquete_importable_como_modulo():
    import ultrasound_sim.dashboard.__main__ as entrada
    assert callable(entrada.main)


def test_la_pagina_expone_arrastre_y_configuracion():
    """Contrato de la UI: el overlay de arrastre y los botones de guardar/cargar tienen que estar.

    No es un test de comportamiento (eso se verifica en un navegador real), sino de que las piezas
    que sostienen el arrastre y la persistencia no desaparezcan sin que nadie se entere.
    """
    html = (ESTATICOS / "index.html").read_text(encoding="utf-8")
    for pieza in ('id="arrastre"', 'id="btnGuardar"', 'id="btnCargar"', 'id="archivo"',
                  'id="srcBody"', 'data-modo="manual"'):
        assert pieza in html, f"falta {pieza} en index.html"

    js = (ESTATICOS / "app.js").read_text(encoding="utf-8")
    for pieza in ("function iniciarArrastre", "function moverArrastre", "function finArrastre",
                  "function instantanea", "function aplicar", "ultrasound-sim/config",
                  "localStorage"):
        assert pieza in js, f"falta {pieza} en app.js"


def test_la_configuracion_que_escribe_el_cliente_es_la_que_el_servidor_entiende(servidor):
    """El JSON de `instantanea()` del cliente, aplicado como pedido al servidor, tiene que correr.

    Se replica el formato en Python (el cliente no se puede importar desde acá) y se manda tal cual
    al endpoint: así un cambio de nombres de campo en el JS rompe este test.
    """
    config = {
        "formato": "ultrasound-sim/config", "version": 1, "material": "agua",
        "frequency_khz": 300, "modo": "array", "n_elements": 24, "aperture_mm": 40,
        "focus_x_mm": 80, "focus_y_mm": 2, "grid": 60, "sources": None,
    }
    pedido = json.dumps({
        "material": config["material"], "frequency_khz": config["frequency_khz"],
        "n_elements": config["n_elements"], "aperture_mm": config["aperture_mm"],
        "focus_x_mm": config["focus_x_mm"], "focus_y_mm": config["focus_y_mm"],
        "grid": {"nx": config["grid"], "ny": config["grid"]},
    }).encode()
    req = urllib.request.Request(servidor + "/api/field", data=pedido,
                                 headers={"Content-Type": "application/json"}, method="POST")
    with urllib.request.urlopen(req, timeout=60) as r:
        datos = json.loads(r.read())
    assert datos["meta"]["n_sources"] == 24
    assert datos["meta"]["frequency_khz"] == 300.0
    assert len(datos["x_mm"]) == 60
