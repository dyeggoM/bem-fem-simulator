"""Dispersor cilíndrico rígido en agua: BEM, FEM y la serie analítica, en el mismo problema (U2).

Es la Fase 4 hecha interactiva: el mismo cilindro, resuelto por los tres caminos, con la
discrepancia entre ellos a la vista. La serie de Bessel no es un adorno —es lo que dice cuál de
los dos métodos se equivoca cuando no coinciden, que es justo cuando un simulador se vuelve
peligroso.

Los tres caminos están copiados del script validado `scripts/cross_validate.py` (mismos
parámetros por defecto: ka en {1,2}, ppw 20, R = 3a, ABC de orden 2, evaluación en r = 1.5a),
porque ese es el que tiene los números medidos en D12/D13 y el que corre en el suite.

Qué NO hace
-----------
- Sólo **cilindro rígido** (Neumann) en un medio homogéneo. Elipse, polígono y dispersor
  penetrable necesitan contornos BEM no circulares y la formulación de contraste — el mismo
  frente que D32 dejó abierto en 3D.
- Sólo agua: un dispersor con impedancia finita no es este problema.
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np
from scipy.interpolate import LinearNDInterpolator

from .. import materials, meshing
from ..analytic import scattering
from ..bem import helmholtz2d as bem
from ..fem import helmholtz2d as fem
from . import validacion

MATERIAL = materials.WATER
FREQUENCY = 100.0e3
ECUACIONES = ("standard", "burton_miller", "chief")
ABC_ORDERS = (1, 2)

_LIMITES = {
    "ka": (0.2, 6.0),
    "ppw": (8.0, 40.0),
    "r_out": (1.5, 6.0),
    "n_panels": (16, 512),
}

# costo: el FEM de la corona crece con (R/a)², y el mapa BEM con la cantidad de puntos
MAX_ELEMENTOS = 80_000
MAX_PUNTOS_MAPA = 8_100
PUNTOS_EVALUACION = 64
PUNTOS_MAPA = 70
#: El mapa se arma con cuadratura baja a propósito: los puntos del mapa están lejos de la
#: superficie (r > a) y el integrando es suave ahí. Medido contra `quad_order=12` en un
#: subconjunto: *0.00 %* de diferencia, y el costo baja de 11.3 s a 1.8 s (lo que lo vuelve
#: interactivo). Para las 64 curvas de comparación sí se usa la cuadratura por defecto.
QUAD_MAPA = 4


def _num(datos: dict[str, Any], clave: str, defecto: float) -> float:
    valor = datos.get(clave, defecto)
    try:
        numero = float(valor)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{clave} tiene que ser un número; llegó {valor!r}") from exc
    if not np.isfinite(numero):
        raise ValueError(f"{clave} no puede ser {numero}")
    return numero


def _en_rango(clave: str, valor: float) -> float:
    bajo, alto = _LIMITES[clave]
    if not (bajo <= valor <= alto):
        raise ValueError(f"{clave} = {valor:g} fuera de rango [{bajo:g}, {alto:g}]")
    return valor


def validar_dispersor(datos: dict[str, Any]) -> dict[str, Any]:
    """Normaliza y valida el pedido. Falla con mensaje, nunca con un `KeyError`."""
    if not isinstance(datos, dict):
        raise ValueError("el pedido tiene que ser un objeto JSON")

    frecuencia = _num(datos, "frecuencia_hz", FREQUENCY)
    if not (20.0e3 <= frecuencia <= 1.0e6):
        raise ValueError(f"frecuencia_hz = {frecuencia:g} fuera de [2e4, 1e6]")

    ka = _en_rango("ka", _num(datos, "ka", 1.0))
    ppw = _en_rango("ppw", _num(datos, "ppw", 20.0))
    r_out = _en_rango("r_out", _num(datos, "r_out", 3.0))

    n_panels = datos.get("n_panels")
    if n_panels is None or n_panels == "" or n_panels == "auto":
        n_panels = None
    else:
        n_panels = int(_en_rango("n_panels", float(n_panels)))

    abc_order = int(datos.get("abc_order", 2))
    if abc_order not in ABC_ORDERS:
        raise ValueError(f"abc_order desconocido: {abc_order}; opciones {list(ABC_ORDERS)}")

    ecuacion = str(datos.get("equation", "standard"))
    if ecuacion not in ECUACIONES:
        raise ValueError(f"ecuation desconocida: {ecuacion!r}; opciones {list(ECUACIONES)}")

    return {
        "frecuencia_hz": frecuencia,
        "ka": ka,
        "ppw": ppw,
        "r_out": r_out,
        "n_panels": n_panels,
        "abc_order": abc_order,
        "equation": ecuacion,
    }


def _radio(ka: float, frecuencia: float) -> float:
    k = abs(MATERIAL.wavenumber(frecuencia))
    return ka / k


def _malla_corona(radio: float, r_out: float, ppw: float, frecuencia: float):
    R = r_out * radio
    h = abs(MATERIAL.wavelength(frecuencia)) / ppw
    n_theta = max(16, int(np.ceil(2.0 * np.pi * R / h)))
    n_radial = max(2, int(np.ceil((R - radio) / h)))
    return meshing.polar_annulus_mesh(radio, R, n_radial, n_theta), R


def estimar_dispersor(datos: dict[str, Any]) -> dict[str, Any]:
    """Malla y costo estimado sin resolver nada."""
    p = validar_dispersor(datos)
    radio = _radio(p["ka"], p["frecuencia_hz"])
    malla, R = _malla_corona(radio, p["r_out"], p["ppw"], p["frecuencia_hz"])
    n_el = int(malla.n_elements)
    if n_el > MAX_ELEMENTOS:
        raise ValueError(
            f"la corona sale de {n_el} elementos; el tope interactivo es {MAX_ELEMENTOS}. "
            f"Bajá ppw (hoy {p['ppw']:g}) o R/a (hoy {p['r_out']:g})."
        )
    k = MATERIAL.wavenumber(p["frecuencia_hz"])
    paneles = p["n_panels"] or bem.panels_for_frequency(radio, k)
    return {
        "ok": True,
        "geo": {
            "radio_mm": radio * 1e3,
            "R_mm": R * 1e3,
            "longitud_onda_mm": abs(MATERIAL.wavelength(p["frecuencia_hz"])) * 1e3,
            "ka": p["ka"],
        },
        "malla": {"n_elementos": n_el, "n_nodos": int(malla.n_nodes)},
        "paneles_bem": int(paneles),
        "estimado_solve_s": round(0.6 * n_el / 1000.0 + 0.0004 * int(paneles) * PUNTOS_MAPA, 2),
    }


def _puntos_circulo(radio: float, n: int) -> np.ndarray:
    th = 2.0 * np.pi * np.arange(n) / n
    return np.column_stack([radio * np.cos(th), radio * np.sin(th)])


def _rel_l2(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


def _interpolar(nodos: np.ndarray, valores: np.ndarray, pts: np.ndarray) -> np.ndarray:
    """Interpolación P1 del campo complejo; exige que los puntos caigan dentro de la malla."""
    interp = LinearNDInterpolator(
        nodos, np.column_stack([np.real(valores), np.imag(valores)]), fill_value=np.nan
    )
    out = interp(pts)
    if np.any(np.isnan(out[:, 0])):
        raise RuntimeError(
            "hay puntos de evaluación fuera del dominio FEM (corona a <= r <= R); "
            "bajá r_out o evaluá más cerca del cilindro"
        )
    return out[:, 0] + 1j * out[:, 1]


def resolver_dispersor(datos: dict[str, Any]) -> dict[str, Any]:
    """Corre los tres caminos y devuelve campo, curvas de comparación y discrepancias."""
    t_inicio = time.perf_counter()
    p = validar_dispersor(datos)
    f = p["frecuencia_hz"]
    k = MATERIAL.wavenumber(f)
    radio = _radio(p["ka"], f)
    malla, R = _malla_corona(radio, p["r_out"], p["ppw"], f)
    if malla.n_elements > MAX_ELEMENTOS:
        raise ValueError(
            f"la corona sale de {malla.n_elements} elementos; el tope interactivo es "
            f"{MAX_ELEMENTOS}. Bajá ppw o R/a."
        )

    # ---- FEM en campo dispersado sobre la corona (receta de cross_validate.py) ----
    aristas = meshing.boundary_edges(malla.triangles)
    r_nodos = np.linalg.norm(malla.nodes[aristas[:, 0]], axis=1)
    exterior = np.isclose(r_nodos, R, atol=1e-12)
    interior = np.isclose(r_nodos, radio, atol=1e-12)

    p_inc = fem.plane_wave(f, MATERIAL, direction=(1.0, 0.0))
    aristas_int = aristas[interior]
    normales_int = meshing.outward_normals(malla.nodes, aristas_int)
    medios_int = meshing.edge_midpoints(malla.nodes, aristas_int)
    dp_dn_inc = p_inc.normal_derivative(medios_int, normales_int)

    k_el = np.full(malla.n_elements, k, dtype=complex)
    t0 = time.perf_counter()
    p_s = fem.solve_scattered(
        malla.nodes, malla.triangles, k_el,
        abc_edges=aristas[exterior], neumann_edges=aristas_int, dp_dn_incident=dp_dn_inc,
        abc_order=p["abc_order"], abc_center=(0.0, 0.0),
    )
    t_fem = time.perf_counter() - t0
    p_total_fem = fem.total_field(p_s, p_inc, malla.nodes)

    # ---- BEM ----
    pts_eval = _puntos_circulo(1.5 * radio, PUNTOS_EVALUACION)
    t0 = time.perf_counter()
    sol = bem.solve_exterior_neumann(
        pts_eval, f, MATERIAL, radio, n_panels=p["n_panels"], equation=p["equation"],
    )
    t_bem = time.perf_counter() - t0

    # ---- serie de referencia ----
    t0 = time.perf_counter()
    serie = scattering.cylinder_scattering(f, MATERIAL, radio, pts_eval, boundary="rigid")
    p_serie = serie["total"]
    t_serie = time.perf_counter() - t0

    p_fem_at = _interpolar(malla.nodes, p_total_fem, pts_eval)

    # el mapa se arma con BEM: evalúa cualquier punto del exterior sin interpolar
    x = np.linspace(-2.2 * R, 2.2 * R, PUNTOS_MAPA)
    gx, gy = np.meshgrid(x, x, indexing="xy")
    pts_mapa = np.column_stack([gx.ravel(), gy.ravel()])
    dentro = np.linalg.norm(pts_mapa, axis=1) < radio
    t0 = time.perf_counter()
    mapa_sol = bem.solve_exterior_neumann(
        pts_mapa[~dentro], f, MATERIAL, radio, n_panels=p["n_panels"],
        quad_order=QUAD_MAPA, equation=p["equation"],
    )
    t_mapa = time.perf_counter() - t0
    amplitud = np.full(pts_mapa.shape[0], np.nan)
    amplitud[~dentro] = np.abs(mapa_sol.p_total)
    amplitud = amplitud.reshape(PUNTOS_MAPA, PUNTOS_MAPA)
    pico = float(np.nanmax(amplitud)) or 1.0
    db = 20.0 * np.log10(np.maximum(amplitud / pico, 1e-6))
    dentro_mapa = dentro.reshape(PUNTOS_MAPA, PUNTOS_MAPA)
    # Dentro del cilindro no hay campo (el BEM sólo resuelve el exterior). Se manda `None`, que es
    # `null` en JSON: un NaN de Python se serializa como el literal `NaN`, que NO es JSON válido y
    # el navegador rechaza el cuerpo entero ("Unexpected token 'N'"). El mapa lo trata como hueco.
    db_mapa = [
        [None if dentro_mapa[i, j] else round(float(db[i, j]), 3) for j in range(PUNTOS_MAPA)]
        for i in range(PUNTOS_MAPA)
    ]

    th = np.degrees(np.arctan2(pts_eval[:, 1], pts_eval[:, 0]))
    orden = np.argsort(th)

    metricas = {
        "fem_vs_bem_pct": 100.0 * _rel_l2(p_fem_at, sol.p_total),
        "fem_vs_serie_pct": 100.0 * _rel_l2(p_fem_at, p_serie),
        "bem_vs_serie_pct": 100.0 * _rel_l2(sol.p_total, p_serie),
        "condicion": float(sol.condition_number),
        "residuo_bem": float(sol.residual),
        "max_bem": float(np.abs(sol.p_total).max()),
        "max_fem": float(np.abs(p_fem_at).max()),
        "max_serie": float(np.abs(p_serie).max()),
    }

    sello = validacion.sello_conjunto({
        # las métricas vienen en %, el criterio está en fracción
        "dispersor_fem_serie": metricas["fem_vs_serie_pct"] / 100.0,
        "dispersor_bem_serie": metricas["bem_vs_serie_pct"] / 100.0,
    })

    return {
        "ok": True,
        "sello": sello,
        "escena": p,
        "geo": {
            "radio_mm": radio * 1e3,
            "R_mm": R * 1e3,
            "ka": p["ka"],
            "frecuencia_hz": f,
            "longitud_onda_mm": abs(MATERIAL.wavelength(f)) * 1e3,
        },
        "malla": {
            "n_elementos": int(malla.n_elements),
            "n_nodos": int(malla.n_nodes),
            "n_paneles_bem": int(sol.n_panels),
            "n_aristas_interior": int(np.count_nonzero(interior)),
            "n_aristas_abc": int(np.count_nonzero(exterior)),
        },
        "mapa": {
            "x_mm": (x * 1e3).tolist(),
            "y_mm": (x * 1e3).tolist(),
            "db": db_mapa,
            "radio_mm": radio * 1e3,
        },
        "comparacion": {
            "angulo_deg": np.round(th[orden], 2).tolist(),
            "fem": np.round(np.abs(p_fem_at[orden]), 6).tolist(),
            "bem": np.round(np.abs(sol.p_total[orden]), 6).tolist(),
            "serie": np.round(np.abs(p_serie[orden]), 6).tolist(),
            "r_eval_mm": 1.5 * radio * 1e3,
        },
        "metricas": metricas,
        "tiempos": {
            "fem_s": round(t_fem, 4),
            "bem_s": round(t_bem, 4),
            "serie_s": round(t_serie, 6),
            "mapa_s": round(t_mapa, 4),
            "total_s": round(time.perf_counter() - t_inicio, 4),
        },
    }
