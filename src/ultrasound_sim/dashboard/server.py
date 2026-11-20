"""Servidor del dashboard: HTTP mínimo sobre la biblioteca estándar.

Sin framework: ``ThreadingHTTPServer`` con cinco rutas. Se eligió así para que el dashboard no
agregue dependencias al paquete (``numpy`` y ``scipy`` ya están) y para no heredar el ciclo de vida de
un framework cuando lo que hace falta es un cálculo y un JSON.

Rutas
-----
``GET  /``              → ``static/index.html``
``GET  /static/<file>`` → assets (``app.js``, ``style.css``, ``plotly.min.js``)
``GET  /api/materials`` → materiales disponibles con sus propiedades
``POST /api/field``     → campo + métricas para un conjunto de fuentes
``GET  /api/field?...`` → atajo equivalente para probar desde el navegador
"""

from __future__ import annotations

import json
import time
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import numpy as np

from .. import cochlea, materials as mats
from . import bem_field
from . import field as campo

ESTATICOS = Path(__file__).resolve().parent / "static"
MATERIALES = {
    "agua": mats.WATER,
    "tejido": mats.SOFT_TISSUE,
    "hueso": mats.BONE_CORTICAL,
}
TIPOS = {
    ".html": "text/html; charset=utf-8",
    ".js": "application/javascript; charset=utf-8",
    ".css": "text/css; charset=utf-8",
    ".png": "image/png",
    ".svg": "image/svg+xml",
}

#: Modelos de radiación disponibles. ``baffle`` = integral de Rayleigh (semiespacio);
#: ``libre`` = lámina finita en espacio libre (capa simple), que además radia hacia atrás.
MODELOS = ("baffle", "libre")

MAX_NX = 400
MAX_NY = 400


def _float(d: dict, clave: str, defecto: float) -> float:
    v = d.get(clave, defecto)
    return float(defecto if v is None or v == "" else v)


def calcular(payload: dict) -> dict:
    """Campo y métricas a partir del pedido JSON. Separado del handler para poder testearlo."""
    t0 = time.perf_counter()

    nombre_mat = str(payload.get("material", "agua"))
    if nombre_mat not in MATERIALES:
        raise ValueError(f"material desconocido: {nombre_mat!r}; opciones: {sorted(MATERIALES)}")
    material = MATERIALES[nombre_mat]

    f_khz = _float(payload, "frequency_khz", 250.0)
    if not 1.0 <= f_khz <= 20_000.0:
        raise ValueError(f"frequency_khz fuera de rango [1, 20000]: {f_khz}")
    frecuencia = f_khz * 1e3

    modelo = str(payload.get("modelo", "baffle"))
    if modelo not in MODELOS:
        raise ValueError(f"modelo desconocido: {modelo!r}; opciones: {sorted(MODELOS)}")

    grilla = payload.get("grid") or {}

    def g(clave: str, defecto: float) -> float:
        """La grilla puede venir anidada (``{"grid": {...}}``) o suelta en el pedido (query string)."""
        if clave in grilla:
            return _float(grilla, clave, defecto)
        return _float(payload, clave, defecto)

    x0 = g("x0", 0.005)
    x1 = g("x1", 0.20)
    y0 = g("y0", -0.05)
    y1 = g("y1", 0.05)
    nx = int(g("nx", 120))
    ny = int(g("ny", 120))
    if not 2 <= nx <= MAX_NX or not 2 <= ny <= MAX_NY:
        raise ValueError(f"nx/ny deben estar en [2, {MAX_NX}]; recibido nx={nx}, ny={ny}")
    if not (x0 < x1):
        raise ValueError(f"la malla necesita x0 < x1; recibido x0={x0}, x1={x1}")
    if modelo == "baffle" and not (0.0 < x0):
        raise ValueError(
            "el modelo con baffle sólo está definido en x > 0; para ver la radiación trasera "
            f"elegí modelo='libre'. Recibido x0={x0}"
        )
    if not (y0 < y1):
        raise ValueError(f"la malla necesita y0 < y1; recibido y0={y0}, y1={y1}")

    fuentes_raw = payload.get("sources")
    pitch = None
    foco_x = None
    apertura = None
    if fuentes_raw:
        fuentes = [campo.Source.from_dict(s) for s in fuentes_raw]
        if len(fuentes) == 1:
            pitch = fuentes[0].width
    else:
        n = int(_float(payload, "n_elements", 32))
        apertura = _float(payload, "aperture_mm", 40.0) * 1e-3
        foco_x = _float(payload, "focus_x_mm", 120.0) * 1e-3
        foco_y = _float(payload, "focus_y_mm", 0.0) * 1e-3
        if not 1 <= n <= 512:
            raise ValueError(f"n_elements debe estar en [1, 512]; recibido {n}")
        fuentes = campo.array_sources(n, apertura, foco_x, foco_y, frecuencia, material)
        pitch = apertura / n

    xs = np.linspace(x0, x1, nx)
    ys = np.linspace(y0, y1, ny)

    if modelo == "libre":
        # Sin baffle: la lámina radia a los dos lados y el campo existe también en x < 0.
        # El punto sobre la lámina (x = 0) no es evaluable (H0 singular), así que se corre la
        # malla medio paso y ningún nodo cae exactamente sobre ella.
        paneles = bem_field.panels_from_sources(fuentes)
        if np.any(np.abs(xs) < 1e-9) and xs.size > 1:
            xs = xs + 0.5 * (xs[1] - xs[0])
        p = bem_field.free_space_field(xs, ys, paneles, frecuencia, material)
    else:
        p = campo.evaluate(xs, ys, fuentes, frecuencia, material)
    met = campo.metrics(xs, ys, p, pitch=pitch, frequency=frecuencia, material=material)

    intensidad = np.abs(p) ** 2
    pico = float(intensidad.max())
    db = 10.0 * np.log10(np.maximum(intensidad / pico, 1e-12))

    return {
        "x_mm": (xs * 1e3).round(4).tolist(),
        "y_mm": (ys * 1e3).round(4).tolist(),
        "p2_db": np.round(db, 2).tolist(),
        "metrics": met,
        "meta": {
            "material": nombre_mat,
            "frequency_khz": f_khz,
            "n_sources": len(fuentes),
            "pitch_mm": None if pitch is None else pitch * 1e3,
            "lambda_mm": material.wavelength(frecuencia) * 1e3,
            # ``D²/(4λ)``: más allá de esto la apertura no puede enfocar, por más que se le pida
            "near_field_mm": None if apertura is None
            else float(apertura ** 2 / (4.0 * material.wavelength(frecuencia))) * 1e3,
            "focus_x_mm": None if foco_x is None else foco_x * 1e3,
            "modelo": modelo,
            "nx": nx,
            "ny": ny,
            "quad": campo.DEFAULT_QUAD,
            "seconds": round(time.perf_counter() - t0, 3),
        },
    }


class Handler(BaseHTTPRequestHandler):
    server_version = "bem-fem-dashboard/1.0"

    def log_message(self, format: str, *args) -> None:  # noqa: A002 (firma de la clase base)
        return

    # ------------------------------------------------------------------ utilidades
    def _json(self, datos: dict, codigo: int = 200) -> None:
        cuerpo = json.dumps(datos, separators=(",", ":")).encode()
        self.send_response(codigo)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(cuerpo)))
        self.end_headers()
        self.wfile.write(cuerpo)

    def _error(self, mensaje: str, codigo: int = 400) -> None:
        self._json({"error": mensaje}, codigo)

    def _static(self, nombre: str) -> None:
        ruta = (ESTATICOS / nombre).resolve()
        if not str(ruta).startswith(str(ESTATICOS.resolve())) or not ruta.is_file():
            self._error("no encontrado", HTTPStatus.NOT_FOUND)
            return
        datos = ruta.read_bytes()
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", TIPOS.get(ruta.suffix, "application/octet-stream"))
        self.send_header("Content-Length", str(len(datos)))
        self.send_header("Cache-Control", "no-cache")
        self.end_headers()
        self.wfile.write(datos)

    # ------------------------------------------------------------------ rutas
    def do_GET(self) -> None:  # noqa: N802 (nombre impuesto por BaseHTTPRequestHandler)
        partes = urlparse(self.path)
        if partes.path in ("/", "/index.html"):
            self._static("index.html")
        elif partes.path.startswith("/static/"):
            self._static(partes.path[len("/static/"):])
        elif partes.path == "/api/materials":
            self._json({
                "materials": [
                    {"id": k, "c_m_s": m.c, "rho": m.rho, "alpha_1mhz_np_m": getattr(m, "alpha_1mhz", 0.0)}
                    for k, m in MATERIALES.items()
                ]
            })
        elif partes.path == "/api/cochlea":
            consulta: dict[str, object] = {k: v[0] for k, v in parse_qs(partes.query).items()}
            try:
                self._json(calcular_coclea(consulta))
            except (ValueError, KeyError) as exc:
                self._error(str(exc))
        elif partes.path == "/api/field":
            consulta = {k: v[0] for k, v in parse_qs(partes.query).items()}
            for clave in ("nx", "ny", "n_elements"):
                if clave in consulta:
                    consulta[clave] = int(float(str(consulta[clave])))
            try:
                self._json(calcular(consulta))
            except (ValueError, KeyError) as exc:
                self._error(str(exc))
        else:
            self._error("no encontrado", HTTPStatus.NOT_FOUND)

    def do_POST(self) -> None:  # noqa: N802
        ruta = urlparse(self.path).path
        if ruta not in ("/api/field", "/api/cochlea"):
            self._error("no encontrado", HTTPStatus.NOT_FOUND)
            return
        largo = int(self.headers.get("Content-Length") or 0)
        if largo > 8 * 1024 * 1024:
            self._error("pedido demasiado grande", HTTPStatus.REQUEST_ENTITY_TOO_LARGE)
            return
        try:
            payload = json.loads(self.rfile.read(largo) or b"{}")
        except json.JSONDecodeError as exc:
            self._error(f"JSON inválido: {exc}")
            return
        cuerpo = payload if isinstance(payload, dict) else {}
        try:
            self._json(calcular_coclea(cuerpo) if ruta == "/api/cochlea" else calcular(cuerpo))
        except (ValueError, KeyError) as exc:
            self._error(str(exc))
        except Exception as exc:  # el dashboard no debe morir por un pedido raro
            self._error(f"error interno: {type(exc).__name__}: {exc}", HTTPStatus.INTERNAL_SERVER_ERROR)


# --------------------------------------------------------------------------- cóclea (corte 2D)
#: Frecuencia y radio del corte con el que están calibrados D24/D29/D30; el modelo no es el mismo
#: que el del mapa (que es una apertura radiando), así que la UI lo presenta aparte.
COCLEA_FREQUENCY_HZ = 1.0e5
COCLEA_RADIUS_M = 20.0e-3


def calcular_coclea(payload: dict) -> dict:
    """Resuelve el corte de cápsula con los tres parámetros que el usuario puede editar.

    Devuelve ``t_interior`` y el armónico ``n = 0`` (D29 regla de D24) más la fracción de potencia
    que absorbe la membrana. Todo lo demás queda en el ballpark documentado: ``c_p = 3500 m/s``,
    ``ρ = 1950 kg/m³`` (D20, VERIFY) y la membrana de ``BASILAR_MEMBRANE_BALLPARK`` (VERIFY).
    """
    t0 = time.perf_counter()

    espesor_mm = _float(payload, "thickness_mm", 2.0)
    if not 0.5 <= espesor_mm <= 12.0:
        raise ValueError(f"thickness_mm debe estar en [0.5, 12]: recibido {espesor_mm}")
    c_s = _float(payload, "c_s", 1900.0)
    if not 300.0 <= c_s < 3500.0:
        raise ValueError(
            f"c_s debe estar en [300, 3500) m/s y ser menor que c_p = 3500 m/s; recibido {c_s}"
        )
    alpha = _float(payload, "alpha_np_m", 0.0)
    if not 0.0 <= alpha <= 200.0:
        raise ValueError(f"alpha_np_m debe estar en [0, 200] Np/m; recibido {alpha}")

    hueso = mats.ElasticSolid(
        name="hueso configurado en la UI", rho=1950.0, c_p=3500.0, c_s=c_s,
        provenance="editado por el usuario en el dashboard; c_p=3500 y rho=1950 son el ballpark "
                   "VERIFY de D20",
    )
    seccion = cochlea.build_cross_section(
        COCLEA_FREQUENCY_HZ, COCLEA_RADIUS_M, espesor_mm * 1e-3,
        ppw=20.0, pad_wavelengths=2.5, min_elements_across_wall=2, bone=hueso,
        window_half_angle=10.0, window_center_angle=0.0,
        window_model="impedance", window_impedance=complex(mats.WATER.rho * mats.WATER.c),
        membrane_angle=90.0,
    )
    sol = cochlea.solve_cross_section(
        seccion, alpha=alpha,
        membrane_impedance=cochlea.BASILAR_MEMBRANE_BALLPARK.impedance(COCLEA_FREQUENCY_HZ),
    )
    return {
        "t_interior": sol.t_interior,
        "t_interior_n0": sol.t_interior_complex,
        "absorbido": sol.power_membrane_fraction,
        "gdl": int(seccion.n_nodes),
        "meta": {
            "thickness_mm": espesor_mm,
            "c_s": c_s,
            "alpha_np_m": alpha,
            "frequency_khz": COCLEA_FREQUENCY_HZ / 1e3,
            "seconds": round(time.perf_counter() - t0, 2),
        },
    }


def build_server(host: str = "127.0.0.1", port: int = 8765) -> ThreadingHTTPServer:
    servidor = ThreadingHTTPServer((host, port), Handler)
    servidor.daemon_threads = True
    return servidor


def serve(host: str = "127.0.0.1", port: int = 8765) -> None:
    servidor = build_server(host, port)
    url = f"http://{host}:{servidor.server_address[1]}/"
    print(f"dashboard en {url}  (Ctrl+C para cortar)")
    try:
        servidor.serve_forever()
    except KeyboardInterrupt:
        print("\ncortado")
    finally:
        servidor.server_close()
