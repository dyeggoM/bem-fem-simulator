"""Cargar un STL y simular su **sección** en el campo (U9).

El simulador es 2D: un STL describe una superficie en 3D, así que lo que entra al modelo es el corte por
un plano —el mismo criterio con el que las aperturas se modelan como su corte lineal—. El contorno de la
sección pasa a ser la frontera de un sólido dentro del fluido, y en 2D esa frontera **sí** se puede
compartir entre las dos regiones (a diferencia de `mallar_stl`, donde dos superficies 3D que se tocan
rompen el mallador; ver D55).

Se mide contra dos casos que el test fija: un disco con las propiedades del fluido **no** perturba la onda
(|p| = 1) y un disco de acero hace sombra. El primero es el que distingue un acople bueno de uno roto.
"""

from __future__ import annotations

import base64
import binascii
import time

import numpy as np

from .. import meshing_volume, stl
from ..fem import helmholtz2d, multicapa
from ..materials import ElasticSolid, Material
from ..meshing import boundary_edges, edge_midpoints

__all__ = ["resolver_stl", "SOLIDOS_STL", "MAX_BYTES"]

MAX_BYTES = 24 * 1024 * 1024

#: Materiales que se pueden asignar a la sección, con su fuente.
SOLIDOS_STL: dict[str, dict] = {
    "acero": {
        "etiqueta": "acero",
        "material": ElasticSolid(name="acero", rho=7800.0, c_p=5900.0, c_s=3200.0, provenance="b"),
        "fuente": "valores de referencia del proyecto (acero estructural)",
        "verificado": True,
    },
    "hueso_cortical": {
        "etiqueta": "hueso cortical",
        "material": ElasticSolid(name="hueso", rho=1896.5, c_p=2471.0, c_s=1600.0,
                                 provenance="Pichardo 2011 / 2017"),
        "fuente": "ρ 1896,5 · c_L 2471 (Pichardo 2011) · c_S 1600 en el centro del rango 1545–1758 "
                  "(Pichardo 2017)",
        "verificado": False,
    },
    "aluminio": {
        "etiqueta": "aluminio",
        "material": ElasticSolid(name="aluminio", rho=2700.0, c_p=6320.0, c_s=3130.0, provenance="b"),
        "fuente": "valores de referencia del proyecto (aluminio 6061)",
        "verificado": True,
    },
    "transparente": {
        "etiqueta": "mismo fluido (control)",
        "material": ElasticSolid(name="transparente", rho=1000.0, c_p=1500.0, c_s=30.0,
                                 provenance="control"),
        "fuente": "control: el sólido tiene las propiedades del agua, así que no debe perturbar la onda",
        "verificado": True,
    },
}


def stl_de_ejemplo(radio: float = 5e-3, alto: float = 20e-3, n: int = 64) -> bytes:
    """Un cilindro cerrado en STL binario, para probar la vista sin tener un archivo a mano."""
    import struct

    ang = np.linspace(0, 2 * np.pi, n, endpoint=False)
    base_ = np.stack([radio * np.cos(ang), radio * np.sin(ang), np.full(n, -alto / 2)], axis=1)
    tapa = base_.copy()
    tapa[:, 2] = alto / 2
    tris = []
    for i in range(n):
        j = (i + 1) % n
        tris += [[base_[i], base_[j], tapa[j]], [base_[i], tapa[j], tapa[i]]]
    cb = np.array([0.0, 0.0, -alto / 2])
    ct = np.array([0.0, 0.0, alto / 2])
    for i in range(n):
        j = (i + 1) % n
        tris += [[cb, base_[j], base_[i]], [ct, tapa[i], tapa[j]]]

    salida = bytearray(84 + 50 * len(tris))
    struct.pack_into("<I", salida, 80, len(tris))
    off = 84
    for t in tris:
        a = np.asarray(t, dtype=float)
        nrm = np.cross(a[1] - a[0], a[2] - a[0])
        ln = np.linalg.norm(nrm)
        nrm = nrm / ln if ln > 0 else np.zeros(3)
        struct.pack_into("<3f", salida, off, *nrm)
        off += 12
        for v in a:
            struct.pack_into("<3f", salida, off, *v)
            off += 12
        struct.pack_into("<H", salida, off, 0)
        off += 2
    return bytes(salida)


def opciones_stl() -> dict:
    return {
        "solidos": [
            {"id": k, "etiqueta": v["etiqueta"], "fuente": v["fuente"], "verificado": v["verificado"]}
            for k, v in SOLIDOS_STL.items()
        ],
        "aviso": (
            "Se simula la SECCIÓN del STL por un plano (el modelo es 2D). El contorno de la sección es la "
            "frontera del sólido dentro del fluido; la onda plana entra desde la izquierda."
        ),
    }


def resolver_stl(datos: dict) -> dict:
    t0 = time.perf_counter()
    b64 = datos.get("datos_b64") or ""
    if not b64:
        raise ValueError("falta 'datos_b64' (el archivo STL en base64)")
    try:
        crudo = base64.b64decode(b64, validate=True)
    except (binascii.Error, ValueError) as exc:
        raise ValueError(f"los datos no son base64 válido: {exc}") from exc
    if len(crudo) > MAX_BYTES:
        raise ValueError(f"el STL pasa de {MAX_BYTES // 1024 // 1024} MB")
    nombre = str(datos.get("nombre") or "STL")

    malla3 = stl.leer_stl_bytes(crudo, nombre)
    stl.verificar_stl(malla3, nombre)

    plano = datos.get("plano") or {}
    eje = str(plano.get("eje", "z"))
    valor = float(plano.get("valor_mm", 0.0)) * 1e-3
    contornos = stl.seccion_por_plano(malla3, eje, valor)
    t_corte = time.perf_counter()

    f = float(datos.get("frecuencia_hz", 250e3))
    ppw = float(datos.get("ppw", 10.0))
    clave = str(datos.get("material", "hueso_cortical"))
    if clave not in SOLIDOS_STL:
        raise ValueError(f"material {clave!r} no disponible: {sorted(SOLIDOS_STL)}")
    solido = SOLIDOS_STL[clave]["material"]

    # dominio: se toma el tamaño de la sección y se le deja 2.5 secciones de agua a cada lado
    todos = np.concatenate(contornos)
    ancho_sec = float(todos[:, 0].max() - todos[:, 0].min())
    alto_sec = float(todos[:, 1].max() - todos[:, 1].min())
    if not (ancho_sec > 0 and alto_sec > 0):
        raise ValueError("la sección tiene extensión nula: revisá el plano de corte")
    W = ancho_sec + 4.0 * max(ancho_sec, 1e-3)
    H = alto_sec + 4.0 * max(alto_sec, 1e-3)

    from ..materials import WATER

    agua = WATER
    lam = agua.wavelength(f)
    tam = min(lam / ppw, max(ancho_sec, alto_sec) / 24.0)
    # `desde` es dónde va el centro de la sección: el dominio se centra en ella
    centro_x = 0.5 * float(todos[:, 0].min() + todos[:, 0].max())
    malla = meshing_volume.mallar_seccion(contornos, ancho=W, alto=H, tamano=tam, desde=centro_x)
    t_malla = time.perf_counter()

    nodos, tris, reg = malla.nodes, malla.tets, malla.regiones
    ar = boundary_edges(tris)
    mid = edge_midpoints(nodos, ar)
    # ABC sólo en entrada y salida (D56/D60): en las paredes laterales la onda plana tiene dp/dn = 0
    x0, x1 = float(nodos[:, 0].min()), float(nodos[:, 0].max())
    abc = ar[(np.abs(mid[:, 0] - x0) < 1e-9) | (np.abs(mid[:, 0] - x1) < 1e-9)]
    medio = np.where(reg == 0, "solido", "agua").astype(object)
    r = multicapa.solve_multicapa(nodos, tris, medio, {"agua": agua}, {"solido": solido}, f,
                                  abc_edges=abc, onda=helmholtz2d.plane_wave(f, agua, (1.0, 0.0)))
    t_solve = time.perf_counter()

    # mapa en la grilla del dibujo
    nx, ny = int(datos.get("nx", 150)), int(datos.get("ny", 110))
    xs = np.linspace(nodos[:, 0].min(), nodos[:, 0].max(), nx)
    ys = np.linspace(nodos[:, 1].min(), nodos[:, 1].max(), ny)
    X, Y = np.meshgrid(xs, ys)
    from matplotlib.tri import LinearTriInterpolator, Triangulation

    tri = Triangulation(nodos[:, 0], nodos[:, 1], tris)
    dentro = tri.get_trifinder()(X, Y) >= 0
    re = np.asarray(LinearTriInterpolator(tri, np.real(r.p_todos))(X, Y), dtype=float)
    im = np.asarray(LinearTriInterpolator(tri, np.imag(r.p_todos))(X, Y), dtype=float)
    campo = np.where(dentro, np.hypot(re, im), np.nan)
    pico = np.nanmax(campo) if np.isfinite(campo).any() else 1.0
    db = 20.0 * np.log10(np.maximum(campo / (pico or 1.0), 1e-6))

    # transmitido en una banda detrás del sólido
    x_max_sec = float(todos[:, 0].max())
    con_p = np.asarray(r.dof_p, dtype=int) >= 0
    banda = con_p & (np.abs(nodos[:, 1]) < 0.25 * alto_sec)
    x = nodos[banda, 0]
    p = np.abs(r.p_todos)[banda]
    detras = x > x_max_sec + 0.25 * ancho_sec
    transmitido = float(p[detras].mean()) if detras.any() else None

    return {
        "ok": True,
        "stl": {
            "nombre": nombre,
            "triangulos": malla3.n_triangulos,
            "resumen_malla": malla3.resumen(),
            "seccion": {
                "contornos": len(contornos),
                "puntos": int(sum(len(c) for c in contornos)),
                "ancho_mm": round(ancho_sec * 1e3, 3),
                "alto_mm": round(alto_sec * 1e3, 3),
                "area_mm2": round(sum(abs(stl.area_contorno(c)) for c in contornos) * 1e6, 3),
            },
        },
        "material": {
            "id": clave, "etiqueta": SOLIDOS_STL[clave]["etiqueta"],
            "fuente": SOLIDOS_STL[clave]["fuente"], "verificado": SOLIDOS_STL[clave]["verificado"],
        },
        "metricas": {
            "transmitido": transmitido,
            "frecuencia_hz": f,
            "ppw": ppw,
            "tamano_malla_mm": round(tam * 1e3, 4),
            "dominio_mm": [round(W * 1e3, 2), round(H * 1e3, 2)],
            "nodos": int(nodos.shape[0]),
            "elementos": int(tris.shape[0]),
        },
        "grilla": {
            "x_mm": np.round(xs * 1e3, 4).tolist(),
            "y_mm": np.round(ys * 1e3, 4).tolist(),
            "db": np.round(db, 3).tolist(),
            "contorno": [np.round(c[:, :2] * 1e3, 4).tolist() for c in contornos],
        },
        "tiempos": {
            "corte_s": round(t_corte - t0, 3),
            "malla_s": round(t_malla - t_corte, 3),
            "solve_s": round(t_solve - t_malla, 3),
        },
    }
