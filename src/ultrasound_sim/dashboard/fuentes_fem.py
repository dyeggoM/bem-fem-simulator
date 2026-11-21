"""El **mismo** problema de fuentes resuelto por FEM, para contrastarlo con el BEM del campo (U9).

Por qué existe
--------------
La vista `campo` calcula el campo de las fuentes por la integral de Rayleigh (BEM de capa simple,
baffle rígido infinito). Hasta ahora no había con qué contrastarlo: la escena de capas y el dispersor
resuelven **otros** problemas. Esto resuelve el mismo: las mismas fuentes, sobre el mismo baffle, en
una caja truncada con ABC.

El problema, escrito igual de los dos lados
------------------------------------------
Fluido en `x > 0`, baffle rígido en `x = 0` (Neumann homogéneo) salvo donde están las fuentes, donde
la velocidad normal está prescrita:

    dp/dx en x=0 = i*w*rho*v(y)   dentro de cada fuente
    dp/dx en x=0 = 0              en el resto del baffle

y la onda sale (ABC) por los tres bordes lejanos. El BEM de `dashboard.field` resuelve exactamente
este problema con la integral de Rayleigh; acá se resuelve con P1 + ABC de orden 1.

Por qué el ABC es de orden 1 y no de orden 2
--------------------------------------------
El de orden 2 exige un contorno **cerrado** (una circunferencia completa). Acá el dominio está abierto
por el baffle, así que no aplica. El error del orden 1 baja al agrandar el dominio, y se mide: es
parte de lo que compara esta vista.

El ensamblado: problema RADIANTE, no dispersor
----------------------------------------------
`helmholtz2d.solve_scattered` está escrito para un **dispersor**: impone la condicion sobre la onda
dispersada y resta el termino de la onda incidente. Aca no hay onda incidente, asi que usarla daba el
campo con la amplitud bien y la fase mal (D54). Lo que corresponde es `assemble(..., p_inc=None)` —que
da A y b = 0 y deja el baffle como Neumann natural— mas el dato de las fuentes en el lado derecho:

    forma debil:  integral sobre Gamma_N de (1/rho)*g*v ds  =>  b[nodo] += (1/rho)*g*L/2

con `g = dp/dn = +i*w*rho*v_s`. La normal saliente del dominio en `x = 0` es `-x` y la fuente empuja
hacia `+x`, de ahi el signo. Como `(1/rho)*g = i*w*v_s`, el dato **no depende de rho**.
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np
import scipy.sparse.linalg as spla

from ..fem import helmholtz2d
from ..materials import Material
from ..meshing import boundary_edges, edge_midpoints, grid_tri_mesh, outward_normals
from . import field as campo

__all__ = ["resolver_fuentes_fem", "estimar_fuentes_fem"]

#: Tope interactivo (mismo criterio que las otras vistas 2D).
MAX_NODOS = 60_000


def _material(datos: dict) -> Material:
    nombre = str(datos.get("material", "agua"))
    from .scene import MATERIALES

    if nombre not in MATERIALES:
        raise ValueError(f"material desconocido {nombre!r}")
    return MATERIALES[nombre]


def _fuentes(datos: dict) -> list[campo.Source]:
    crudas = datos.get("sources")
    if not crudas:
        raise ValueError("hacen falta 'sources' (las fuentes del campo)")
    return [campo.Source.from_dict(s) for s in crudas]


def _malla(fuentes, frecuencia, material, ppw: float, largo: float, alto: float):
    lam = material.wavelength(frecuencia)
    h = lam / ppw
    nx = int(np.clip(round(largo / h), 12, 300))
    ny = int(np.clip(round(2.0 * alto / h), 12, 300))
    # el paso en y incluye los bordes de las fuentes: si no, una fuente mas angosta que el paso
    # desaparece de la malla y el campo sale sin ella
    bordes = sorted({round(f.y - f.width / 2, 9) for f in fuentes}
                    | {round(f.y + f.width / 2, 9) for f in fuentes})
    ys = np.unique(np.concatenate([np.linspace(-alto, alto, ny + 1), np.clip(bordes, -alto, alto)]))
    xs = np.linspace(0.0, largo, nx + 1)
    malla = grid_tri_mesh(xs, ys)
    if malla.nodes.shape[0] > MAX_NODOS:
        raise ValueError(
            f"la malla sale de {malla.nodes.shape[0]} nodos y el tope interactivo es {MAX_NODOS}: "
            f"baja los elementos por longitud de onda o el largo del dominio"
        )
    return malla


def estimar_fuentes_fem(datos: dict) -> dict[str, Any]:
    """Costo antes de resolver, para que la UI no lance a ciegas."""
    fuentes = _fuentes(datos)
    f = float(datos.get("frecuencia_hz", 250e3))
    material = _material(datos)
    ppw = float(datos.get("ppw", 8.0))
    largo = float(datos.get("largo_mm", 160.0)) * 1e-3
    alto = float(datos.get("alto_mm", 60.0)) * 1e-3
    malla = _malla(fuentes, f, material, ppw, largo, alto)
    n = malla.nodes.shape[0]
    return {
        "ok": True,
        "n_nodos": int(n),
        "n_elementos": int(malla.triangles.shape[0]),
        "estimado_s": round(0.02 * (n / 3000.0) ** 1.2, 2),
        "fuentes": len(fuentes),
    }


def resolver_fuentes_fem(datos: dict) -> dict[str, Any]:
    """Resuelve el campo de las fuentes por FEM y lo compara con el BEM de la vista `campo`."""
    t0 = time.perf_counter()
    fuentes = _fuentes(datos)
    f = float(datos.get("frecuencia_hz", 250e3))
    material = _material(datos)
    ppw = float(datos.get("ppw", 8.0))
    largo = float(datos.get("largo_mm", 160.0)) * 1e-3
    alto = float(datos.get("alto_mm", 60.0)) * 1e-3
    quad = int(datos.get("quad", 8))
    nx = int(datos.get("nx", 120))
    ny = int(datos.get("ny", 120))

    malla = _malla(fuentes, f, material, ppw, largo, alto)
    t_malla = time.perf_counter()

    nodos = malla.nodes
    aristas = boundary_edges(malla.triangles)
    medios = edge_midpoints(nodos, aristas)

    # el baffle es el borde en x = 0; el resto es contorno lejano (ABC)
    en_baffle = np.abs(medios[:, 0]) < 1e-9
    abc_edges = aristas[~en_baffle]
    neumann = aristas[en_baffle]
    normales = outward_normals(nodos, abc_edges)

    omega = 2.0 * np.pi * f
    k = complex(material.wavenumber(f))
    rho = float(material.rho)

    # velocidad normal prescrita por arista del baffle
    v = np.zeros(neumann.shape[0], dtype=complex)
    for j, (cx, cy) in enumerate(edge_midpoints(nodos, neumann)):
        for s in fuentes:
            if abs(cy - s.y) <= 0.5 * s.width + 1e-12:
                v[j] = s.amp * np.exp(1j * s.phase)
                break

    k_el = np.full(malla.triangles.shape[0], k, dtype=complex)
    rho_el = np.full(malla.triangles.shape[0], rho)

    a, b = helmholtz2d.assemble(
        nodos, malla.triangles, k_el,
        abc_edges=abc_edges, p_inc=None,
        element_rho=rho_el, normals=normales,
    )
    largo_arista = np.linalg.norm(nodos[neumann[:, 1]] - nodos[neumann[:, 0]], axis=1)
    # b[nodo] += (1/rho)*g*L/2 con g = i*w*rho*v  =>  i*w*v*L/2
    contrib = 1j * omega * v * largo_arista / 2.0
    np.add.at(b, neumann[:, 0], contrib)
    np.add.at(b, neumann[:, 1], contrib)
    p = spla.spsolve(a.tocsc(), b)
    t_solve = time.perf_counter()

    # BEM de referencia en la MISMA grilla. `campo.evaluate` exige x > 0 estricto (el fluido esta
    # delante del baffle), asi que la grilla arranca un paso adentro.
    xs = np.linspace(largo / max(nx, 2), largo, nx)
    ys = np.linspace(-alto, alto, ny)
    p_bem = campo.evaluate(xs, ys, fuentes, f, material, quad=quad)
    t_bem = time.perf_counter()

    fem = _interpolar(nodos, malla.triangles, p, xs, ys)

    escala = np.abs(p_bem).max() or 1.0
    util = np.abs(p_bem) > 1e-3 * escala
    delta = np.abs(fem - p_bem)
    error = float(delta[util].max() / escala) if util.any() else None
    error_l2 = float(np.sqrt(np.sum(delta[util] ** 2) / max(np.sum(np.abs(p_bem[util]) ** 2), 1e-30)))
    razon = float(np.abs(fem).max() / escala)

    db = 20.0 * np.log10(np.maximum(np.abs(fem) / (np.abs(fem).max() or 1.0), 1e-6))
    db_bem = 20.0 * np.log10(np.maximum(np.abs(p_bem) / escala, 1e-6))

    return {
        "ok": True,
        "grilla": {
            "x_mm": np.round(xs * 1e3, 4).tolist(),
            "y_mm": np.round(ys * 1e3, 4).tolist(),
            "db": np.round(db, 3).tolist(),
            "db_bem": np.round(db_bem, 3).tolist(),
            "dif_db": np.round(db - db_bem, 3).tolist(),
        },
        "metricas": {
            "error_relativo": error,
            "error_l2_relativo": error_l2,
            "razon_amplitudes": razon,
            "fuentes": len(fuentes),
            "frecuencia_hz": f,
            "ppw": ppw,
            "nodos": int(nodos.shape[0]),
            "elementos": int(malla.triangles.shape[0]),
            "aristas_abc": int(abc_edges.shape[0]),
            "aristas_baffle": int(neumann.shape[0]),
            "aristas_fuente": int(np.count_nonzero(v)),
        },
        "tiempos": {
            "malla_s": round(t_malla - t0, 3),
            "solve_s": round(t_solve - t_malla, 3),
            "bem_s": round(t_bem - t_solve, 3),
        },
    }


def _interpolar(nodos, triangles, valores, xs, ys):
    """Interpola P1 sobre la malla en la grilla ``xs x ys``; fuera del dominio queda 0."""
    from matplotlib.tri import LinearTriInterpolator, Triangulation

    tri = Triangulation(nodos[:, 0], nodos[:, 1], triangles)
    xx, yy = np.meshgrid(xs, ys)
    dentro = tri.get_trifinder()(xx, yy) >= 0
    # parte real e imaginaria por separado: pasar un fasor complejo descarta la imaginaria
    re = np.asarray(LinearTriInterpolator(tri, np.real(valores))(xx, yy), dtype=float)
    im = np.asarray(LinearTriInterpolator(tri, np.imag(valores))(xx, yy), dtype=float)
    out = re + 1j * im
    return np.where(dentro, out, 0.0)
