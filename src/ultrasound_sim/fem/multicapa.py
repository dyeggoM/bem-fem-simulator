"""Acople de **N regiones** en 2D: cadena de capas fluido/sólido (U9).

Por qué una función aparte
--------------------------
`solve_coupled` acopla **un** fluido con **un** sólido. La cadena piel -> hueso -> líquido necesita varios
materiales a la vez, y algunos son sólidos elásticos (piel, hueso) mientras otros son fluidos (perilinfa).

Cómo sale de la misma formulación, sin física nueva
---------------------------------------------------
Con **malla conforme** y elementos P1, las interfaces entre dos regiones del **mismo tipo** no necesitan
ningún bloque: los dos lados comparten los mismos gdl, así que

* la continuidad de ``u`` (sólido-sólido) es automática porque es la misma incógnita en el mismo nodo;
* el equilibrio de tracciones sale de la forma débil al ensamblar cada región con **su** rigidez;
* entre dos fluidos pasa lo mismo con ``p``, y la continuidad de ``(1/ρ)·∂p/∂n`` sale sola porque cada
  elemento lleva su ``1/ρ``.

La única interfaz que sí necesita bloque explícito es **fluido-sólido**, porque las incógnitas son
distintas (``p`` escalar y ``u`` vectorial) y hay que acoplarlas con ``C`` y ``ω²·Cᵀ``, exactamente el
mismo bloque de `coupled2d`.

Por eso esto puede ser una función nueva y no una reescritura: se reusan las primitivas de `coupled2d`
(`interface_edges`, la ABC, el término incidente) y ninguna fórmula cambia.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from ..materials import ElasticSolid, Material
from .assembly import p1_tri_local, weights_to_local
from .coupled2d import (
    _boundary_owner_map,
    abc_boundary_block,
    abc_incident_rhs,
    interface_edges,
)
from .elastic2d import p1_elastic_local

__all__ = ["Medio", "ResultadoMulticapa", "solve_multicapa", "medios_de_capas"]

_TWO_PI = 2.0 * np.pi


@dataclass(frozen=True)
class Medio:
    """Una región: su material y de qué tipo es."""

    nombre: str
    material: Material | ElasticSolid

    @property
    def es_solido(self) -> bool:
        return isinstance(self.material, ElasticSolid)


@dataclass(frozen=True)
class ResultadoMulticapa:
    p: np.ndarray | None            # presión en los nodos de fluido (None si no hay fluido)
    u: np.ndarray | None            # desplazamiento (N, 2) (None si no hay sólidos)
    dof_p: np.ndarray               # (N,) índice de gdl de presión o -1
    dof_u: np.ndarray               # (N, 2) índice de gdl de desplazamiento o -1
    p_todos: np.ndarray             # presión en TODOS los nodos (0 en los que no tienen gdl)
    u_todos: np.ndarray             # desplazamiento en todos los nodos
    medio_de_elemento: np.ndarray   # (M,) nombre del medio por elemento
    frecuencia: float
    #: |p| con la escala del solve, para poder comparar entre corridas
    presion_max: float


def medios_de_capas(capas: list[tuple[str, float]], fluidos: dict[str, Material],
                    solidos: dict[str, ElasticSolid], desde: float = 0.0) -> list[tuple[str, float]]:
    """Valida una lista de capas y devuelve ``[(medio, espesor)]`` comprobando que el medio exista."""
    salida = []
    for nombre, espesor in capas:
        if nombre not in fluidos and nombre not in solidos:
            raise ValueError(
                f"la capa {nombre!r} no está definida: medios disponibles "
                f"{sorted(set(fluidos) | set(solidos))}"
            )
        if not (espesor > 0):
            raise ValueError(f"la capa {nombre!r} tiene espesor {espesor!r} y tiene que ser positivo")
        salida.append((nombre, float(espesor)))
    return salida


def solve_multicapa(
    nodes,
    triangles,
    medio_de_elemento,
    fluidos: dict[str, Material],
    solidos: dict[str, ElasticSolid],
    frecuencia: float,
    *,
    abc_edges=None,
    onda=None,
    abc_order: int = 1,
    abc_center=None,
    neumann_edges=None,
    dp_dn=None,
) -> ResultadoMulticapa:
    """Resuelve el acople de N regiones con onda plana incidente (o sin fuente si ``onda is None``).

    Parameters
    ----------
    medio_de_elemento : ``(M,)`` nombre del medio de cada elemento.
    fluidos, solidos : diccionarios ``{nombre: material}``. Un nombre no puede estar en los dos.
    onda : objeto con ``.value(points)`` y ``.normal_derivative(points, normals)``
        (``fem.helmholtz2d.PlaneWave``). Entra como dato de borde en la ABC, sobre el campo total.
    neumann_edges, dp_dn : para excitar con **fuentes** en vez de una onda plana. ``dp_dn`` es
        ``∂p/∂n`` por arista de ``neumann_edges`` (la normal saliente del dominio), y entra en el lado
        derecho como ``∫ (1/ρ)·(∂p/∂n)·q ds`` con la regla del punto medio. Una fuente que empuja hacia
        adentro del dominio con velocidad ``v`` sobre un baffle de normal saliente ``−x̂`` tiene
        ``∂p/∂n = +i·ω·ρ·v``.
    """
    x = np.asarray(nodes, dtype=float)
    t = np.asarray(triangles, dtype=int)
    m = np.asarray(medio_de_elemento)
    if m.shape[0] != t.shape[0]:
        raise ValueError(f"medio_de_elemento tiene {m.shape[0]} valores y hay {t.shape[0]} elementos")

    repetidos = set(fluidos) & set(solidos)
    if repetidos:
        raise ValueError(f"estos nombres están como fluido y como sólido a la vez: {sorted(repetidos)}")
    desconocidos = sorted(set(m.tolist()) - set(fluidos) - set(solidos))
    if desconocidos:
        raise ValueError(f"medio_de_elemento usa medios no definidos: {desconocidos}")
    if not fluidos and not solidos:
        raise ValueError("hacen falta medios")

    n = x.shape[0]
    omega = _TWO_PI * float(frecuencia)
    k_f_ref = (next(iter(fluidos.values())).wavenumber(frecuencia) if fluidos else 1.0)

    # ---- gdl ---------------------------------------------------------------
    nodo_solido = np.zeros(n, dtype=bool)
    nodo_fluido = np.zeros(n, dtype=bool)
    for nombre, mat in solidos.items():
        sel = m == nombre
        if sel.any():
            nodo_solido[t[sel].ravel()] = True
    for nombre in fluidos:
        sel = m == nombre
        if sel.any():
            nodo_fluido[t[sel].ravel()] = True

    reps_u = np.flatnonzero(nodo_solido)
    reps_p = np.flatnonzero(nodo_fluido)
    n_u = 2 * reps_u.size
    n_p = reps_p.size
    size = n_u + n_p

    dof_u = -np.ones((n, 2), dtype=int)
    dof_u[reps_u, 0] = 2 * np.arange(reps_u.size)
    dof_u[reps_u, 1] = 2 * np.arange(reps_u.size) + 1
    dof_p = -np.ones(n, dtype=int)
    dof_p[reps_p] = n_u + np.arange(reps_p.size)

    filas: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    vals: list[np.ndarray] = []

    def _add(rdofs, cdofs, local):
        n_r = rdofs.shape[1]
        n_c = cdofs.shape[1]
        filas.append(np.repeat(rdofs, n_c, axis=1).ravel())
        cols.append(np.tile(cdofs, (1, n_r)).ravel())
        vals.append(np.asarray(local).reshape(-1))

    # ---- bloques sólidos, uno por material ---------------------------------
    for nombre, mat in solidos.items():
        sel = m == nombre
        if not sel.any():
            continue
        tris = t[sel]
        k_s, m_s, _ = p1_elastic_local(x, tris, mat)
        _add(dof_u[tris].reshape(-1, 6), dof_u[tris].reshape(-1, 6),
             k_s - mat.rho * omega ** 2 * m_s)

    # ---- bloques fluidos, uno por material ---------------------------------
    for nombre, mat in fluidos.items():
        sel = m == nombre
        if not sel.any():
            continue
        tris = t[sel]
        k_f = mat.wavenumber(frecuencia)
        m_f, k_l, _ = p1_tri_local(x, tris)
        _add(dof_p[tris], dof_p[tris], (k_l - k_f ** 2 * m_f) / mat.rho)

    # ---- acople fluido-sólido ---------------------------------------------
    es_solido = np.array([nombre in solidos for nombre in m.tolist()], dtype=bool)
    reg = np.where(es_solido, "solid", "fluid")
    if es_solido.any() and not es_solido.all():
        iface, iface_n = interface_edges(x, t, reg)
        if iface.shape[0]:
            a_i, b_i = iface[:, 0], iface[:, 1]
            largo = np.linalg.norm(x[b_i] - x[a_i], axis=1)
            ud = dof_u[iface].reshape(-1, 4)
            pd = dof_p[iface]
            if np.any(ud < 0) or np.any(pd < 0):
                raise RuntimeError("un nodo de interfaz no tiene los tres gdl")
            m_loc = weights_to_local(largo)
            c_local = np.empty((iface.shape[0], 4, 2), dtype=float)
            nx_e = iface_n[:, 0][:, None, None]
            ny_e = iface_n[:, 1][:, None, None]
            c_local[:, 0, :] = (nx_e * m_loc)[:, 0, :]
            c_local[:, 1, :] = (ny_e * m_loc)[:, 0, :]
            c_local[:, 2, :] = (nx_e * m_loc)[:, 1, :]
            c_local[:, 3, :] = (ny_e * m_loc)[:, 1, :]
            _add(ud, pd, c_local)
            _add(pd, ud, omega ** 2 * np.transpose(c_local, (0, 2, 1)))

    a = sp.coo_matrix(
        (np.concatenate(vals), (np.concatenate(filas), np.concatenate(cols))), shape=(size, size)
    ).tocsr()

    b = np.zeros(size, dtype=complex)

    # ---- ABC en la frontera de fluido -------------------------------------
    if n_p:
        bmap = _boundary_owner_map(t)
        if abc_edges is None:
            keep = [entrada[:2] for entrada in bmap.values() if not es_solido[entrada[2]]]
            abc = np.array(keep, dtype=int).reshape(-1, 2)
        else:
            abc = np.asarray(abc_edges, dtype=int).reshape(-1, 2)
        if abc.shape[0]:
            # la ABC usa el fluido "de referencia": en una cadena con líquidos distintos, el contorno
            # lejano está en el último medio y ése es el que tiene que absorber
            rho_abc = float(next(iter(fluidos.values())).rho)
            blk = abc_boundary_block(x, abc, k_f_ref, rho=rho_abc, abc_order=abc_order,
                                     abc_center=abc_center)
            coo = blk.tocoo()
            rd, cd = dof_p[coo.row], dof_p[coo.col]
            ok = (rd >= 0) & (cd >= 0)
            filas_b = rd[ok]
            cols_b = cd[ok]
            a = a + sp.coo_matrix((coo.data[ok], (filas_b, cols_b)), shape=(size, size)).tocsr()
            if onda is not None:
                rhs = abc_incident_rhs(x, abc, onda, k_f_ref, rho=rho_abc, abc_order=abc_order,
                                       abc_center=abc_center)
                has = dof_p >= 0
                np.add.at(b, dof_p[has], rhs[has])

    # ---- fuentes: dato de Neumann prescrito (baffle con un arreglo) ----------
    if neumann_edges is not None:
        ne = np.asarray(neumann_edges, dtype=int).reshape(-1, 2)
        if dp_dn is None:
            raise ValueError("neumann_edges requiere dp_dn (∂p/∂n por arista)")
        dato = np.asarray(dp_dn, dtype=complex).reshape(-1)
        if dato.size != ne.shape[0]:
            raise ValueError(
                f"dp_dn tiene {dato.size} valores y hay {ne.shape[0]} aristas de Neumann"
            )
        rd = dof_p[ne]
        if np.any(rd < 0):
            raise RuntimeError("una arista de Neumann no toca ningún nodo con gdl de presión")
        # ∫ (1/ρ)·(∂p/∂n)·q ds con regla del punto medio: cada nodo recibe (1/ρ)·g·L/2. El 1/ρ ya viene
        # incluido en dp_dn (quien llama pasa el dato físico dividido por ρ)
        largo_n = np.linalg.norm(x[ne[:, 1]] - x[ne[:, 0]], axis=1)
        contrib = dato * largo_n / 2.0
        np.add.at(b, rd[:, 0], contrib)
        np.add.at(b, rd[:, 1], contrib)

    sol = spla.spsolve(a.tocsc(), b)

    u = sol[:n_u].reshape(-1, 2) if n_u else None
    p = sol[n_u:] if n_p else None
    u_todos = np.zeros((n, 2), dtype=complex)
    if n_u:
        u_todos[reps_u] = u
    p_todos = np.zeros(n, dtype=complex)
    if n_p:
        p_todos[reps_p] = p
    return ResultadoMulticapa(
        p=p, u=u, dof_p=dof_p, dof_u=dof_u,
        p_todos=p_todos, u_todos=u_todos,
        medio_de_elemento=m, frecuencia=float(frecuencia),
        presion_max=float(np.abs(p).max()) if n_p else 0.0,
    )
