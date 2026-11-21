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


# ------------------------------------------------------------------ los dos modelos de radiación
def test_el_campo_delantero_no_depende_del_baffle():
    """Identidad física que sostiene el dashboard: con baffle y sin baffle, el campo en x>0 es el
    mismo módulo con 180° de fase. La diferencia entre los dos modelos no está adelante: está en
    que sin baffle hay campo detrás.
    """
    from ultrasound_sim.dashboard import bem_field

    f = 250e3
    fuentes = campo.array_sources(16, 0.040, 0.090, 0.0, f, W)
    xs = np.linspace(0.055, 0.095, 9)
    ys = np.array([0.0, 0.006, -0.006])

    con_baffle = campo.evaluate(xs, ys, fuentes, f, W)
    sin_baffle = bem_field.free_space_field(xs, ys, bem_field.panels_from_sources(fuentes), f, W)

    desvio = np.abs(sin_baffle + con_baffle).max() / np.abs(con_baffle).max()
    assert desvio < 1e-6, f"p_libre debería ser -p_baffle exacto; desvío {desvio:.2e}"
    # y el módulo es idéntico, que es lo que se ve en pantalla
    assert np.allclose(np.abs(sin_baffle), np.abs(con_baffle), rtol=1e-6)


def test_modelo_libre_extiende_la_malla_al_semiespacio_trasero():
    r = calcular({"modelo": "libre", "frequency_khz": 250, "n_elements": 16, "aperture_mm": 40,
                  "focus_x_mm": 90, "grid": {"nx": 40, "ny": 30, "x0": -0.05, "x1": 0.08}})
    assert r["meta"]["modelo"] == "libre"
    assert r["x_mm"][0] < 0.0
    atras = [fila[i] for fila in r["p2_db"] for i, x in enumerate(r["x_mm"]) if x < 0.0]
    # el lóbulo trasero es tan fuerte como el delantero: no está a −40 dB, está cerca de 0
    assert max(atras) > -3.0, f"no hay campo trasero: máximo {max(atras):.1f} dB"


def test_modelo_baffle_rechaza_malla_negativa_y_modelo_desconocido():
    with pytest.raises(ValueError, match="baffle"):
        calcular({"modelo": "baffle", "grid": {"x0": -0.01, "x1": 0.10, "nx": 20, "ny": 20}})
    with pytest.raises(ValueError, match="modelo"):
        calcular({"modelo": "difraccion-exacta"})


def test_modelo_libre_por_http(servidor):
    pedido = json.dumps({"modelo": "libre", "frequency_khz": 250, "n_elements": 8,
                         "aperture_mm": 20, "focus_x_mm": 60,
                         "grid": {"nx": 30, "ny": 20, "x0": -0.04, "x1": 0.07}}).encode()
    req = urllib.request.Request(servidor + "/api/field", data=pedido,
                                 headers={"Content-Type": "application/json"}, method="POST")
    with urllib.request.urlopen(req, timeout=90) as r:
        datos = json.loads(r.read())
    assert datos["meta"]["modelo"] == "libre"
    assert datos["x_mm"][0] < 0.0
    assert len(datos["x_mm"]) == 30


def test_api_coclea_con_los_tres_parametros(servidor):
    """El corte se resuelve por HTTP y con los tres parámetros del panel.

    El valor de referencia es el de D29: t = 2 mm, c_s = 1900 m/s, sin atenuación → 0.5319 con la
    membrana puesta (sin ella sería 0.9722). Si el cableado de `alpha` o del hueso configurable se
    rompe, este número se mueve.
    """
    pedido = json.dumps({"thickness_mm": 2.0, "c_s": 1900.0, "alpha_np_m": 0.0}).encode()
    req = urllib.request.Request(servidor + "/api/cochlea", data=pedido,
                                 headers={"Content-Type": "application/json"}, method="POST")
    with urllib.request.urlopen(req, timeout=120) as r:
        d = json.loads(r.read())
    print(f"\n  t_interior = {d['t_interior']:.4f}  n0 = {d['t_interior_n0']:.4f}  "
          f"absorbe = {100 * d['absorbido']:.1f} %  ({d['meta']['seconds']} s)")
    assert d["t_interior"] == pytest.approx(0.5319, abs=2e-3)
    assert d["t_interior_n0"] == pytest.approx(0.1791, abs=2e-3)
    assert 0.5 < d["absorbido"] < 0.95      # la membrana absorbe la mayor parte de lo que entra
    assert d["meta"]["thickness_mm"] == 2.0

    # con atenuación la transmisión tiene que bajar
    pedido = json.dumps({"thickness_mm": 2.0, "c_s": 1900.0, "alpha_np_m": 15.31}).encode()
    req = urllib.request.Request(servidor + "/api/cochlea", data=pedido,
                                 headers={"Content-Type": "application/json"}, method="POST")
    with urllib.request.urlopen(req, timeout=120) as r:
        con_alpha = json.loads(r.read())
    assert con_alpha["t_interior"] < d["t_interior"]
    assert con_alpha["t_interior"] == pytest.approx(0.4808, abs=3e-3)


def test_api_coclea_rechaza_parametros_fuera_de_rango(servidor):
    for malo in ({"thickness_mm": 20}, {"c_s": 4000}, {"alpha_np_m": -1}):
        req = urllib.request.Request(servidor + "/api/cochlea", data=json.dumps(malo).encode(),
                                     headers={"Content-Type": "application/json"}, method="POST")
        with pytest.raises(urllib.error.HTTPError) as exc:
            urllib.request.urlopen(req, timeout=60)
        assert exc.value.code == 400


# ------------------------------------------------------------------ escena FEM (U1)


def _post(url: str, cuerpo: dict) -> dict:
    datos = json.dumps(cuerpo).encode("utf-8")
    req = urllib.request.Request(url, data=datos, headers={"content-type": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.loads(r.read())


ESCENA_HTTP = {
    "frecuencia_hz": 100.0e3,
    "ppw": 12,
    "alto_m": 20.0e-3,
    "capas": [
        {"material": "agua", "espesor_m": 20.0e-3},
        {"material": "hueso", "espesor_m": 10.0e-3},
        {"material": "agua", "espesor_m": 20.0e-3},
    ],
}


def test_api_scene_cost_no_resuelve(servidor):
    r = _post(servidor + "/api/scene/cost", ESCENA_HTTP)
    assert r["ok"] is True
    assert r["grilla"]["n_elementos"] > 0
    assert r["estimado_solve_s"] > 0
    assert "campo" not in r, "estimar no puede resolver"


def test_api_scene_solve_devuelve_el_campo_y_las_metricas(servidor):
    r = _post(servidor + "/api/scene/solve", ESCENA_HTTP)
    assert r["ok"] is True
    m = r["metricas"]
    assert m["error_R_abs"] < 0.02
    assert m["error_T_abs"] < 0.02
    assert m["energia"] == pytest.approx(1.0, abs=0.02)
    n_y = len(r["grilla"]["y_mm"])
    n_x = len(r["grilla"]["x_mm"])
    assert len(r["campo"]["db"]) == n_y
    assert len(r["campo"]["db"][0]) == n_x


def test_api_scene_solve_rechaza_incidencia_oblicua(servidor):
    req = urllib.request.Request(
        servidor + "/api/scene/solve",
        data=json.dumps(dict(ESCENA_HTTP, angulo_deg=30.0)).encode("utf-8"),
        headers={"content-type": "application/json"},
    )
    with pytest.raises(urllib.error.HTTPError) as exc:
        urllib.request.urlopen(req, timeout=60)
    assert exc.value.code == 400
    assert b"incidencia normal" in exc.value.read()


def test_api_materials_expone_la_escena_y_sus_topes(servidor):
    _estado, cuerpo = _get(servidor + "/api/materials")
    datos = json.loads(cuerpo)
    ids = {m["id"] for m in datos["scene_materials"]}
    assert {"agua", "tejido", "hueso", "perilinfa", "endolinfa"} <= ids
    assert datos["scene_limits"]["max_elementos"] > 0
    assert "abc1" in datos["scene_limits"]["abc"]


# ------------------------------------------------------------------ dispersor (U2)


def test_api_scatterer_cost_y_solve(servidor):
    costo = _post(servidor + "/api/scatterer/cost", {"ka": 1.0})
    assert costo["ok"] is True
    assert costo["malla"]["n_elementos"] > 0
    assert "metricas" not in costo

    r = _post(servidor + "/api/scatterer/solve", {"ka": 1.0})
    assert r["ok"] is True
    m = r["metricas"]
    assert m["fem_vs_serie_pct"] < 2.0
    assert m["bem_vs_serie_pct"] < 2.0
    assert len(r["mapa"]["db"]) == len(r["mapa"]["x_mm"])


def test_api_scatterer_solve_emite_json_sin_nan(servidor):
    """El cuerpo tiene que ser JSON estricto: con `NaN` el navegador no parsea nada."""
    datos = json.dumps({"ka": 1.0}).encode("utf-8")
    req = urllib.request.Request(
        servidor + "/api/scatterer/solve", data=datos, headers={"content-type": "application/json"}
    )
    with urllib.request.urlopen(req, timeout=120) as r:
        crudo = r.read()
    assert b"NaN" not in crudo, "el servidor emitió NaN: JSON.parse lo rechaza"
    json.loads(crudo)  # revienta si no es JSON estricto


def test_api_scatterer_rechaza_ecuacion_desconocida(servidor):
    req = urllib.request.Request(
        servidor + "/api/scatterer/solve",
        data=json.dumps({"ka": 1.0, "equation": "galerkin"}).encode("utf-8"),
        headers={"content-type": "application/json"},
    )
    with pytest.raises(urllib.error.HTTPError) as exc:
        urllib.request.urlopen(req, timeout=60)
    assert exc.value.code == 400


# ------------------------------------------------------------------ cóclea (U3)


def test_api_cochlea_field_cost(servidor):
    r = _post(servidor + "/api/cochlea/field/cost", {"thickness_mm": 2.0})
    assert r["ok"] is True
    assert r["malla"]["n_nodos"] > 10_000
    assert "metricas" not in r


def test_api_cochlea_field_devuelve_mapa_y_metricas(servidor):
    r = _post(servidor + "/api/cochlea/field", {"thickness_mm": 2.0, "c_s": 1900.0})
    assert r["ok"] is True
    m = r["metricas"]
    assert m["t_interior"] == pytest.approx(0.53187, rel=0.01)
    assert m["balance_energia"] < 0.01
    assert len(r["mapa"]["db"]) == len(r["mapa"]["x_mm"])


def test_api_cochlea_field_emite_json_sin_complejos(servidor):
    """`sol.as_dict()` trae `window_impedance` complejo: si no se sanea, el navegador no recibe nada."""
    datos = json.dumps({"thickness_mm": 2.0}).encode("utf-8")
    req = urllib.request.Request(
        servidor + "/api/cochlea/field", data=datos, headers={"content-type": "application/json"}
    )
    with urllib.request.urlopen(req, timeout=180) as r:
        crudo = r.read()
    assert b"NaN" not in crudo
    cuerpo = json.loads(crudo)
    assert cuerpo["detalle"]["window_impedance"]["re"] > 0
    assert cuerpo["detalle"]["window_impedance"]["im"] == 0.0


# ------------------------------------------------------------------ barridos y cola (U5)


def test_api_sweeps_lista_el_barrido_y_su_estimacion(servidor):
    _estado, cuerpo = _get(servidor + "/api/sweeps")
    datos = json.loads(cuerpo)
    barrido = datos["barridos"][0]
    assert barrido["id"] == "sensibilidad"
    assert barrido["puntos"] > 10
    assert barrido["estimado_s"] > 0
    assert "D30" in barrido["referencia"]


def test_api_jobs_lanza_y_se_sigue_por_id(servidor):
    """El ciclo completo por HTTP: encolar, consultar estado, y que el id desconocido dé 404."""
    datos = json.dumps({}).encode("utf-8")
    req = urllib.request.Request(
        servidor + "/api/jobs/sensibilidad", data=datos,
        headers={"content-type": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=30) as r:
        lanzado = json.loads(r.read())
    job_id = lanzado["job"]["id"]
    assert lanzado["ok"] is True
    assert lanzado["job"]["estado"] in ("en_cola", "corriendo")

    with urllib.request.urlopen(f"{servidor}/api/jobs/{job_id}", timeout=30) as r:
        estado = json.loads(r.read())
    assert estado["id"] == job_id
    assert 0.0 <= estado["progreso"] <= 1.0

    # se cancela para no dejar 90 s de barrido colgados en el suite
    req = urllib.request.Request(f"{servidor}/api/jobs/{job_id}/cancel", data=b"{}",
                                 headers={"content-type": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as r:
        cancelado = json.loads(r.read())
    assert cancelado["estado"] in ("cancelado", "corriendo")

    with pytest.raises(urllib.error.HTTPError) as exc:
        urllib.request.urlopen(f"{servidor}/api/jobs/noexiste", timeout=10)
    assert exc.value.code == 404


def test_api_jobs_lista_todos(servidor):
    _estado, cuerpo = _get(servidor + "/api/jobs")
    assert "jobs" in json.loads(cuerpo)


# ------------------------------------------------------------------ 3D por rebanadas (U4)


def test_api_volume3d_cost_y_solve(servidor):
    cuerpo = {
        "frecuencia_hz": 100.0e3, "ppw": 10, "alto_m": 22.0e-3, "fuente": "puntual",
        "capas": [{"material": "agua", "espesor_m": 22.0e-3}],
    }
    costo = _post(servidor + "/api/volume3d/cost", cuerpo)
    assert costo["ok"] is True
    assert costo["malla"]["n_nodos"] > 0
    assert costo["tope_nodos"] > 0
    assert "metricas" not in costo

    r = _post(servidor + "/api/volume3d/solve", cuerpo)
    assert r["ok"] is True
    assert r["metricas"]["error_green_pct"] < 10.0
    assert len(r["rebanadas"]["xy"]["db"]) == len(r["rebanadas"]["xy"]["y_mm"])


def test_api_volume3d_rechaza_oblicua(servidor):
    req = urllib.request.Request(
        servidor + "/api/volume3d/solve",
        data=json.dumps({"fuente": "plano", "angulo_deg": 20.0}).encode("utf-8"),
        headers={"content-type": "application/json"},
    )
    with pytest.raises(urllib.error.HTTPError) as exc:
        urllib.request.urlopen(req, timeout=60)
    assert exc.value.code == 400


# ------------------------------------------------------------------ trazabilidad (U6)


def test_api_validacion_publica_el_catalogo(servidor):
    _estado, cuerpo = _get(servidor + "/api/validacion")
    datos = json.loads(cuerpo)
    ids = {c["id"] for c in datos["criterios"]}
    assert {"escena_rt_r", "volumen_green", "barrido_d30"} <= ids
    for c in datos["criterios"]:
        assert c["fuente"], f"{c['id']} sin fuente"


def test_las_corridas_por_http_traen_sello(servidor):
    """El sello viaja en la respuesta: la UI no lo calcula."""
    r = _post(servidor + "/api/scene/solve", ESCENA_HTTP)
    assert "sello" in r
    assert r["sello"]["estado"] == "validado"
    assert "Fuente:" in r["sello"]["texto"]


def test_la_pagina_expone_la_vista_3d():
    html = (ESTATICOS / "index.html").read_text(encoding="utf-8")
    for pieza in ('data-vista="volumen"', 'id="vistaVolumen"', 'id="mapaVolumen"',
                  'id="corteVolumen"', 'id="volFuente"', 'volumen.js'):
        assert pieza in html, f"falta {pieza} en index.html"
    js = (ESTATICOS / "volumen.js").read_text(encoding="utf-8")
    for pieza in ("function volResolver", "function volDibujar", "/api/volume3d/solve",
                  "mostrarSello", "scaleanchor"):
        assert pieza in js, f"falta {pieza} en volumen.js"


def test_la_pagina_expone_la_vista_de_barridos():
    html = (ESTATICOS / "index.html").read_text(encoding="utf-8")
    for pieza in ('data-vista="barridos"', 'id="vistaBarridos"', 'id="barProgreso"',
                  'id="barRanking"', 'id="btnBarrido"', 'barridos.js'):
        assert pieza in html, f"falta {pieza} en index.html"
    js = (ESTATICOS / "barridos.js").read_text(encoding="utf-8")
    for pieza in ("/api/jobs/", "/api/sweeps", "function barPintar", "comparacion"):
        assert pieza in js, f"falta {pieza} en barridos.js"


def test_la_pagina_expone_la_vista_de_cochlea():
    html = (ESTATICOS / "index.html").read_text(encoding="utf-8")
    for pieza in ('data-vista="cochlea"', 'id="vistaCochlea"', 'id="mapaCochlea"',
                  'id="cocFieldT"', 'id="cocFieldBalance"', 'id="cocFieldMemb"', 'cochlea.js'):
        assert pieza in html, f"falta {pieza} en index.html"
    js = (ESTATICOS / "cochlea.js").read_text(encoding="utf-8")
    for pieza in ("function cocfResolver", "function cocfDibujar", "/api/cochlea/field",
                  "mostrarSello"):
        assert pieza in js, f"falta {pieza} en cochlea.js"


def test_la_pagina_expone_la_vista_de_dispersor():
    html = (ESTATICOS / "index.html").read_text(encoding="utf-8")
    for pieza in ('data-vista="dispersor"', 'id="vistaDispersor"', 'id="mapaDisp"',
                  'id="comparacionDisp"', 'id="disFemSerie"', 'dispersor.js'):
        assert pieza in html, f"falta {pieza} en index.html"
    js = (ESTATICOS / "dispersor.js").read_text(encoding="utf-8")
    for pieza in ("function disResolver", "function disDibujar", "/api/scatterer/solve",
                  "/api/scatterer/cost", "mostrarSello"):
        assert pieza in js, f"falta {pieza} en dispersor.js"


def test_la_pagina_expone_la_vista_de_escena():
    html = (ESTATICOS / "index.html").read_text(encoding="utf-8")
    for pieza in ('data-vista="escena"', 'id="vistaEscena"', 'id="capas"', 'id="btnEscena"',
                  'id="mapaEscena"', 'id="escR"', 'id="escErrR"', 'escena.js'):
        assert pieza in html, f"falta {pieza} en index.html"
    js = (ESTATICOS / "escena.js").read_text(encoding="utf-8")
    # el criterio ya no vive acá: lo evalúa el servidor y la vista lo pinta (U6)
    for pieza in ("function escResolver", "function escPintarCapas", "function escDibujar",
                  "function mostrarSello", "/api/scene/solve", "/api/scene/cost"):
        assert pieza in js, f"falta {pieza} en escena.js"


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
