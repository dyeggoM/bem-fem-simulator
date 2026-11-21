"""El **mismo** problema de fuentes resuelto por FEM, para contrastarlo con el BEM del campo (U9).

ESTADO: **NO VALIDADO — no está conectado a la UI**.

Lo que sí está verificado: el BEM de referencia (`dashboard.field.evaluate`) coincide con la solución
exacta de un monopolo 2D (`H₀(k·r)`) en módulo relativo 1.000 y fase 0.0° en cuatro distancias, así
que la referencia es buena.

Lo que no: este FEM, con las mismas fuentes y el mismo baffle, **no reproduce ese campo**. Medido con
una fuente de 2 mm en un dominio de 200×200 mm2: la razón de amplitudes va de 0.31 a 0.92 según el
punto y las fases difieren entre 45° y 165°, o sea que no es un factor global que se pueda corregir con
una constante. Sospechosos, en orden: (a) la forma débil de `solve_scattered` integra el dato de
Neumann con un signo/convención que no es la de un problema *radiante* (está escrita para un dispersor
con onda incidente), (b) el ABC-1 en un contorno rectangular abierto por el baffle, y (c) la malla no
resuelve el ancho de cada fuente (3 aristas para 2 mm).

Se deja el módulo en el repo porque la mitad del trabajo —la malla, el marcado de aristas del baffle,
la clasificación ABC/Neumann y la comparación contra el BEM— es correcta y reutilizable. La UI NO lo
expone hasta que cierre contra la referencia. Ver D54.
"""

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

    ∂p/∂x|₀ = i·ω·ρ·v(y)      dentro de cada fuente
    ∂p/∂x|₀ = 0              en el resto del baffle

y la onda sale (ABC) por los tres bordes lejanos. El BEM de `dashboard.field` resuelve exactamente
este problema con la integral de Rayleigh; acá se resuelve con P1 + ABC de orden 2.

El signo del dato
-----------------
La normal **saliente del dominio** en `x = 0` es `−x̂`, así que el dato de Neumann es
`∂p/∂n = −i·ω·ρ·v`. `solve_scattered` impone `∂p_s/∂n = −∂p_i/∂n` y suma `∫ dp_dn_incident·v`, de
modo que hay que pasarle `dp_dn_incident = +i·ω·ρ·v`. Se verifica contra el BEM: con el signo
invertido el campo sale en contrafase y la razón de amplitudes sigue siendo 1.
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np

from ..fem import helmholtz2d
from ..materials import WATER, Material
from ..meshing import boundary_edges, edge_midpoints, grid_tri_mesh, outward_normals
from . import field as campo

__all__ = ["resolver_fuentes_fem", "estimar_fuentes_fem"]

#: Tope interactivo (mismo criterio que las otras vistas 2D).
MAX_NODOS = 40_000


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
    nx = int(np.clip(round(largo / h), 12, 260))
    ny = int(np.clip(round(2.0 * alto / h), 12, 260))
    # el paso en y se ajusta para caer justo en los bordes de las fuentes: si no, una fuente de ancho
    # menor que el paso desaparece de la malla
    bordes = sorted({round(f.y - f.width / 2, 9) for f in fuentes}
                    | {round(f.y + f.width / 2, 9) for f in fuentes})
    ys = np.unique(np.concatenate([np.linspace(-alto, alto, ny + 1), np.clip(bordes, -alto, alto)]))
    xs = np.linspace(0.0, largo, nx + 1)
    malla = grid_tri_mesh(xs, ys)
    if malla.nodes.shape[0] > MAX_NODOS:
        raise ValueError(
            f"la malla sale de {malla.nodes.shape[0]} nodos y el tope interactivo es {MAX_NODOS}: "
            f"bajá los elementos por longitud de onda o el largo del dominio"
        )
    return malla


def estimar_fuentes_fem(datos: dict) -> dict[str, Any]:
    """Costo antes de resolver, para que la UI no lance a ciegas."""
    fuentes = _fuentes(datos)
    f = float(datos.get("frecuencia_hz", 250e3))
    material = _material(datos)
    ppw = float(datos.get("ppw", 10.0))
    largo = float(datos.get("largo_mm", 160.0)) * 1e-3
    alto = float(datos.get("alto_mm", 70.0)) * 1e-3
    malla = _malla(fuentes, f, material, ppw, largo, alto)
    n = malla.nodes.shape[0]
    return {
        "ok": True,
        "n_nodos": int(n),
        "n_elementos": int(malla.triangles.shape[0]),
        "estimado_s": round(0.15 * (n / 3000.0) ** 1.5, 2),
        "fuentes": len(fuentes),
    }


def resolver_fuentes_fem(datos: dict) -> dict[str, Any]:
    """Resuelve el campo de las fuentes por FEM y lo compara con el BEM de la vista `campo`."""
    t0 = time.perf_counter()
    fuentes = _fuentes(datos)
    f = float(datos.get("frecuencia_hz", 250e3))
    material = _material(datos)
    ppw = float(datos.get("ppw", 10.0))
    largo = float(datos.get("largo_mm", 160.0)) * 1e-3
    alto = float(datos.get("alto_mm", 70.0)) * 1e-3
    quad = int(datos.get("quad", 8))
    nx = int(datos.get("nx", 120))
    ny = int(datos.get("ny", 120))

    malla = _malla(fuentes, f, material, ppw, largo, alto)
    t_malla = time.perf_counter()

    nodos = malla.nodes
    aristas = boundary_edges(malla.triangles)
    medios = edge_midpoints(nodos, aristas)
    normales = outward_normals(nodos, aristas)

    # el baffle es el borde en x = 0; el resto es contorno lejano (ABC)
    en_baffle = np.abs(medios[:, 0]) < 1e-9
    abc_edges = aristas[~en_baffle]
    neumann = aristas[en_baffle]

    omega = 2.0 * np.pi * f
    k = complex(material.wavenumber(f))
    rho = float(material.rho)

    # velocidad normal prescrita en el punto medio de cada arista del baffle
    v = np.zeros(neumann.shape[0], dtype=complex)
    for j, (cx, cy) in enumerate(edge_midpoints(nodos, neumann)):
        for s in fuentes:
            if abs(cy - s.y) <= 0.5 * s.width + 1e-12:
                v[j] = s.amp * np.exp(1j * s.phase)
                break

    # ∂p/∂n = −i·ω·ρ·v con n saliente (−x̂) ⇒ dp_dn_incident = +i·ω·ρ·v.
    # El signo se verifica, no se deriva: con `+` el error era 117 % con razón de amplitudes 1.03, o
    # sea contrafase pura. Con `−` el campo queda alineado (ver la medición en el docstring y en D54).
    dato = -1j * omega * rho * v

    k_el = np.full(malla.triangles.shape[0], k, dtype=complex)
    rho_el = np.full(malla.triangles.shape[0], rho)

    p = helmholtz2d.solve_scattered(
        nodos, malla.triangles, k_el,
        abc_edges=abc_edges, neumann_edges=neumann, dp_dn_incident=dato,
        element_rho=rho_el, normals=normales[~en_baffle],
        # ABC de orden 1: el orden 2 exige un contorno cerrado (circunferencia completa) y acá el
        # dominio está abierto por el baffle. El error del ABC-1 baja al agrandar el dominio, y se
        # mide: es parte de lo que compara esta vista.
        abc_order=1,
        neumann_normals=normales[en_baffle],
    )
    t_solve = time.perf_counter()

    # BEM de referencia en la MISMA grilla. `campo.evaluate` exige x > 0 estricto (el fluido está
    # delante del baffle), así que la grilla arranca un paso adentro.
    xs = np.linspace(largo / max(nx, 2), largo, nx)
    ys = np.linspace(-alto, alto, ny)
    p_bem = campo.evaluate(xs, ys, fuentes, f, material, quad=quad)
    t_bem = time.perf_counter()

    # FEM interpolado a la grilla (búsqueda del triángulo contenedor, como en las otras vistas)
    fem = _interpolar(nodos, malla.triangles, p, xs, ys)

    escala = np.abs(p_bem).max() or 1.0
    util = np.abs(p_bem) > 1e-3 * escala
    error = float(np.abs(fem[util] - p_bem[util]).max() / escala) if util.any() else None
    delta = np.abs(fem[util] - p_bem[util]) if util.any() else np.array([0.0])
    error_l2 = float(np.sqrt(np.sum(delta ** 2) / max(np.sum(np.abs(p_bem[util]) ** 2), 1e-30)))
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
    """Interpola P1 sobre la malla en la grilla ``xs × ys`` (fuera del dominio ⇒ 1e-30)."""
    from matplotlib.tri import Triangulation

    tri = Triangulation(nodos[:, 0], nodos[:, 1], triangles)
    xx, yy = np.meshgrid(xs, ys)
    out = np.full(xx.shape, np.nan + 0j)
    mask = tri.get_trifinder()(xx, yy) >= 0
    from matplotlib.tri import LinearTriInterpolator

    interp = LinearTriInterpolator(tri, valores)
    real = interp(xx, yy)
    out = np.asarray(real, dtype=complex)
    return np.where(mask, out, 0.0)
