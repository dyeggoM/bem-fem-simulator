"""La cadena piel → hueso → cóclea resuelta por FEM y contrastada contra la referencia (U9).

Qué hace
--------
Arma un dominio 2D rectangular con agua a la entrada, las capas pedidas en el medio y agua a la salida,
incide una onda plana normal y resuelve con `fem.multicapa.solve_multicapa`. Devuelve el mapa de presión,
el corte axial y **el coeficiente de transmisión medido**, comparado contra la matriz de transferencia
cuando la comparación es válida.

Cuándo NO hay referencia
------------------------
Con alguna capa sólida elástica (hueso) hay cizalla y la fórmula de capas fluidas no vale. En ese caso
`t_teorico` sale `None` y la UI dice «sin referencia analítica», que es la verdad: no se compara contra un
número que no corresponde al mismo problema.
"""

from __future__ import annotations

import time
from typing import cast

import numpy as np

from .. import cadena as materiales
from ..materials import Material
from ..fem import helmholtz2d, multicapa
from ..meshing import boundary_edges, edge_midpoints, grid_tri_mesh

__all__ = ["resolver_cadena", "opciones_cadena"]

MAX_NODOS = 60_000


def opciones_cadena() -> dict:
    return {
        "medios": materiales.resumen_medios(),
        "aviso": (
            "La mayoría de estos valores son parámetros de MODELADO publicados, no mediciones de cada "
            "tejido: no existe un conjunto experimental completo validado de 0,1 a 1 MHz. Piel, grasa, "
            "músculo y líquidos entran como FLUIDOS porque no hay c_S citable para tejido blando; el "
            "hueso entra como sólido elástico."
        ),
    }


def _capas(datos: dict) -> list[tuple[str, float]]:
    crudas = datos.get("capas") or []
    if not crudas:
        raise ValueError("hacen falta 'capas': [{'medio': ..., 'espesor_mm': ...}, ...]")
    salida = []
    for i, capa in enumerate(crudas):
        medio = str(capa.get("medio", "")).strip()
        if not medio:
            raise ValueError(f"la capa {i + 1} no dice el medio")
        espesor = float(capa.get("espesor_mm", 0.0)) * 1e-3
        if not (espesor > 0):
            raise ValueError(f"la capa {i + 1} ({medio}) tiene espesor {capa.get('espesor_mm')!r} mm")
        salida.append((medio, espesor))
    return salida


def resolver_cadena(datos: dict) -> dict:
    t0 = time.perf_counter()
    capas = _capas(datos)
    f = float(datos.get("frecuencia_hz", 500e3))
    ppw = float(datos.get("ppw", 14.0))
    agua_mm = float(datos.get("agua_mm", 12.0))
    alto_mm = float(datos.get("alto_mm", 24.0))
    nx = int(datos.get("nx", 140))
    ny = int(datos.get("ny", 90))

    fluidos, solidos = materiales.medios_de([m for m, _ in capas] + ["agua"])
    # el agua tiene que estar como fluido; si el usuario la pone como capa, ya está
    agua = fluidos.get("agua")
    if agua is None:
        agua = cast(Material, materiales.MEDIOS_CADENA["agua"].material)
        fluidos["agua"] = agua

    lam = agua.wavelength(f)
    h = lam / ppw
    entrada = agua_mm * 1e-3
    espesor_total = sum(e for _, e in capas)
    salida = max(entrada, 8e-3)
    largo = entrada + espesor_total + salida

    n_largo = int(round(largo / h))
    n_alto = int(round(2 * (alto_mm * 1e-3) / h))
    if n_largo * n_alto > MAX_NODOS:
        raise ValueError(
            f"la malla saldría de ~{n_largo * n_alto} nodos y el tope interactivo es {MAX_NODOS}: "
            f"bajá los elementos por longitud de onda o el tamaño del dominio"
        )
    xs = np.linspace(0.0, largo, n_largo + 1)
    ys = np.linspace(-alto_mm * 1e-3, alto_mm * 1e-3, n_alto + 1)
    malla = grid_tri_mesh(xs, ys)
    nodos, tris = malla.nodes, malla.triangles
    t_malla = time.perf_counter()

    centro_x = nodos[tris][:, :, 0].mean(axis=1)
    medio_elem = np.full(centro_x.shape, "agua", dtype=object)
    corte = []
    x = entrada
    for nombre, espesor in capas:
        medio_elem[(centro_x >= x) & (centro_x < x + espesor)] = nombre
        corte.append((nombre, x, x + espesor))
        x += espesor

    ar = boundary_edges(tris)
    mid = edge_midpoints(nodos, ar)
    # ABC sólo en entrada y salida: en las paredes laterales la onda plana tiene dp/dn = 0 y una ABC ahí
    # refleja (D56)
    abc = ar[(np.abs(mid[:, 0]) < 1e-12) | (np.abs(mid[:, 0] - largo) < 1e-12)]

    onda = helmholtz2d.plane_wave(f, agua, (1.0, 0.0))
    r = multicapa.solve_multicapa(
        nodos, tris, medio_elem, fluidos, solidos, f, abc_edges=abc, onda=onda
    )
    t_solve = time.perf_counter()

    # ---- referencia analítica (sólo si no hay sólidos) ----------------------
    t_teorico = materiales.t_capas_analitico(capas, f) if not solidos else None
    hay_solido = bool(solidos)

    # ---- campo transmitido medido ------------------------------------------
    eje = np.abs(nodos[:, 1]) < 1e-12
    # con sólidos, los nodos del sólido no tienen gdl de presión: se toman sólo los que sí
    con_p = np.asarray(r.dof_p, dtype=int)[eje] >= 0
    x_m = nodos[eje, 0][con_p]
    p_m = np.abs(r.p_todos[eje])[con_p]
    # la zona de medición arranca 4 mm después de la última interfaz, para no leer el campo cercano
    zona = x_m > (entrada + sum(e for _, e in capas) + 0.004)
    t_medido = float(p_m[zona].mean()) if zona.any() else None

    # ---- mapa --------------------------------------------------------------
    X, Y = np.meshgrid(np.linspace(0, largo, nx), np.linspace(-alto_mm * 1e-3, alto_mm * 1e-3, ny))
    from matplotlib.tri import LinearTriInterpolator, Triangulation

    tri = Triangulation(nodos[:, 0], nodos[:, 1], tris)
    dentro = tri.get_trifinder()(X, Y) >= 0
    # El campo es un fasor: hay que interpolar PARTE REAL Y PARTE IMAGINARIA por separado.
    # `LinearTriInterpolator(tri, z_complejo)` tira la parte imaginaria con un ComplexWarning y el mapa
    # terminaba mostrando |Re(p)|, que no es la amplitud.
    re = np.asarray(LinearTriInterpolator(tri, np.real(r.p_todos))(X, Y), dtype=float)
    im = np.asarray(LinearTriInterpolator(tri, np.imag(r.p_todos))(X, Y), dtype=float)
    campo = np.hypot(re, im)
    campo = np.where(dentro, campo, np.nan)
    pico = np.nanmax(campo) if np.isfinite(campo).any() else 1.0
    db = 20.0 * np.log10(np.maximum(campo / (pico or 1.0), 1e-6))

    t_mapa = time.perf_counter()

    return {
        "ok": True,
        "capas": [
            {"medio": nombre, "etiqueta": materiales.MEDIOS_CADENA[nombre].etiqueta,
             "desde_mm": round(a * 1e3, 3), "hasta_mm": round(b * 1e3, 3),
             "espesor_mm": round((b - a) * 1e3, 3)}
            for nombre, a, b in corte
        ],
        "metricas": {
            "t_medido": t_medido,
            "t_teorico": t_teorico,
            "error_pct": (None if (t_teorico is None or t_medido is None or t_teorico == 0)
                          else round(100 * abs(t_medido - t_teorico) / t_teorico, 2)),
            "hay_solido": hay_solido,
            "referencia": ("matriz de transferencia de capas fluidas" if t_teorico is not None
                           else "sin referencia analítica: hay capa sólida y la fórmula de capas "
                                "fluidas no incluye cizalla"),
            "frecuencia_hz": f,
            "ppw": ppw,
            "nodos": int(nodos.shape[0]),
            "elementos": int(tris.shape[0]),
            "espesor_total_mm": round(espesor_total * 1e3, 3),
            "lambda_agua_mm": round(lam * 1e3, 3),
        },
        "grilla": {
            "x_mm": np.round(np.linspace(0, largo, nx) * 1e3, 4).tolist(),
            "y_mm": np.round(np.linspace(-alto_mm, alto_mm, ny), 4).tolist(),
            "db": np.round(db, 3).tolist(),
        },
        "corte_axial": {
            "x_mm": np.round(x_m * 1e3, 4).tolist(),
            "abs_p": np.round(p_m, 6).tolist(),
        },
        "medios": {n: {"etiqueta": m.etiqueta, "tipo": m.tipo, "fuente": m.fuente,
                       "verificado": m.verificado, "nota": m.nota}
                   for n, m in materiales.MEDIOS_CADENA.items()},
        "tiempos": {
            "malla_s": round(t_malla - t0, 3),
            "solve_s": round(t_solve - t_malla, 3),
            "pos_s": round(t_mapa - t_solve, 3),
        },
    }
