"""Acoplamiento fluido acústico / sólido elástico 2D en plano deformación (Fase 2).

Física (``e^{−iωt}``, D4), una sola malla, cada elemento marcado ``'fluid'`` o ``'solid'``, y ``n``
la normal **saliente del sólido** sobre la interfaz Γ::

    sólido:  ∇·σ + ρ_s·ω²·u = 0 ,     σ = λ·tr(ε)·I + 2·μ·ε
    fluido:  ∇²p + k²·p = 0 ,         k = ω/c_f + i·α
    Γ:       σ·n = −p·n        (tracción del fluido sobre el sólido)
             ∂p/∂n = ρ_f·ω²·(u·n)   (continuidad de velocidad normal, v_s = −iω·u)

Forma débil (identidades de Green, los dos signos salen **positivos**)::

    sólido:  ∫(σ:ε(v) − ρ_s·ω²·u·v) dA + ∫_Γ p·(n·v) ds = 0
    fluido:  ∫(∇p·∇q − k²·p·q)/ρ_f dA + ω²·∫_Γ (u·n)·q ds = 0

y en bloque, con ``C_ij = ∫_Γ (n·φ_i)·ψ_j ds`` (``φ_i`` base vectorial P1 del sólido, ``ψ_j`` base
escalar P1 del fluido)::

    ⎡ K_s − ρ_s·ω²·M_s        C      ⎤ ⎡u⎤   ⎡f_s⎤
    ⎢                                ⎥ ⎢ ⎥ = ⎢   ⎥
    ⎣ ω²·Cᵀ          K_f − k²·M_f/ρ_f⎦ ⎣p⎦   ⎣f_f⎦

Sobre el signo de ``C``
-----------------------
El plan escribe el bloque como ``[−C; −ω²·Cᵀ]``. Las dos versiones describen *el mismo problema
físico* salvo un signo global de ``C``, y ese signo global es exactamente la elección de hacia
dónde apunta la normal en ``C``: ``(n saliente, [+C; +ω²Cᵀ]) ≡ (n entrante, [−C; −ω²Cᵀ])``. Con
``n`` saliente y el bloque del plan el campo ``p`` sale idéntico y ``u`` sale con un desfase π
(equivale a ``σ·n = +p·n``): |T|, |R|, los ángulos y las resonancias no cambian, pero el
desplazamiento queda invertido. Acá se implementa la pareja consistente con las dos condiciones
físicas de arriba —``n`` saliente y ``+C``—, que es la misma matriz del plan con ``C`` definida
sobre la normal entrante; la verificación es contra `analytic.helmholtz1d.stack_response`
comparando el ``T`` **complejo** (fase incluida) y contra el balance de energía.

Grados de libertad (el punto delicado)
--------------------------------------
Una sola malla; cada nodo recibe gdl según los elementos que lo tocan:

===============================  ==========================================
nodo toca sólo elementos fluid  1 gdl  → ``p``
nodo toca sólo elementos solid  2 gdl  → ``(u_x, u_y)``
nodo de interfaz (los dos)      3 gdl  → ``(u_x, u_y, p)``
===============================  ==========================================

La numeración se hace **por nodo canónico** y es contigua: primero todos los gdl de desplazamiento
(2 por nodo sólido, en orden nodo mayor) y después todos los de presión. `CoupledDofs` guarda los
mapas ``u_dof`` ``(N, 2)`` y ``p_dof`` ``(N,)`` con ``-1`` donde el nodo no tiene ese gdl; las
matrices globales se ensamblan directamente sobre ese espacio (sin gdl "fantasma"), de modo que
la matriz de bloques tiene tamaño ``(n_u + n_p)²`` y no ``3N``.

Interfaz
--------
Se marcan como interfaz las **aristas interiores con un elemento de cada lado** (la malla tiene que
ser conforme: cada arista interior pertenece exactamente a dos triángulos). La normal es la
saliente del elemento **sólido**, o sea ``(t_y, −t_x)`` con ``t = p_j − p_i`` y ``(i, j)`` la arista
tal como la recorre el triángulo sólido CCW (la misma convención de `meshing.outward_normals`).

Fronteras
---------
- ``abc_edges``: condición absorbente sobre el campo **total**, con el término incidente para la
  onda que entra:

  - ``abc_order=1`` (default): ABC de primer orden ``∂p_s/∂n = i·k·p_s``, con el término incidente
    en las caras iluminadas (``n·d < 0``), igual que `fem.helmholtz2d.assemble`. Es válida en
    ductos con caras alineadas con la onda normal (D12); para contornos curvos no sirve — con la
    cápsula de `cochlea` mide 12 % L2 y 17–21 % en la cara de sombra, sin mejorar al agrandar el
    dominio (D18).
  - ``abc_order=2``: BGT-2 ``∂p_s/∂n = (i·k − 1/(2r) + i/(8k r²))·p_s + (i/(2k))·∂²_s p_s`` sobre un
    contorno circular de radio ``r`` respecto de ``abc_center`` (obligatorio), con el término
    tangencial ensamblado por ``∫∂²_s p_s·v = −∫∂_s p_s·∂_s v + [∂_s p_s·v]`` (el corchete se
    cancela en una cadena cerrada y se ensambla en los extremos de una cadena abierta). El término
    incidente entra en **toda** la ABC, no sólo en las caras iluminadas: la condición se impone de
    forma consistente sobre ``p_s``, así que la fuente espuria que la máscara deja en la sombra
    (``i·k·p_i·(d·n − 1)``) desaparece. Ver la documentación de ``solve_coupled`` para qué pasa si
    ``abc_center`` es ``None``, está desplazado o cae sobre el contorno.

  ``abc_edges=None`` ⇒ ABC en **toda** la frontera de los elementos de fluido.
- El resto de la frontera de fluido es Neumann natural (pared rígida, ``∂p/∂n = 0``).
- La frontera del sólido es libre de tracción (``σ·n = 0``) salvo donde se impongan gdl con
  ``dirichlet``.
- ``dirichlet``: gdl esenciales nulos, ``{'ux': nodos, 'uy': nodos, 'p': nodos}``. El caso de uso
  de la Fase 2 es el **rodillo** ``uy = 0`` en las paredes laterales de una placa: convierte el
  problema en deformación uniaxial (``ε_yy = 0``), que es la hipótesis con la que la capa elástica
  del stack analítico tiene velocidad ``c_p``. Sin ese rodillo una placa de caras libres transmite
  a la velocidad de placa ``c_p·√((1−2ν)/(1−ν)²)`` (≈8 % menor en acero) y la comparación contra
  el stack con ``c_p`` falla por física, no por código.
- ``periodic_pairs``: pares de nodos que se identifican (malla periódica en ``y``). Imprescindible
  para incidencia oblicua: con paredes rígidas la onda oblicua no satisface la condición de borde y
  las reflexiones contaminan el campo. Con ``H = λ_f/sin θ_i`` la onda plana incidente —y todas las
  transmitidas, porque ``k_y`` se conserva— son exactamente periódicas, así que la condición
  periódica es exacta y no una aproximación.

Convención de ``p`` devuelto: ``NaN`` donde el nodo no tiene gdl de fluido; ``u``: ``NaN`` donde no
hay sólido.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from ..materials import ElasticSolid, Material
from ..meshing import edge_midpoints, outward_normals
from .assembly import p1_tri_local, weights_to_local
from .elastic2d import p1_elastic_local

__all__ = [
    "CoupledDofs",
    "CoupledSolution",
    "abc_boundary_block",
    "abc_chain_ends",
    "abc_incident_rhs",
    "interface_edges",
    "segment_absorbed_power",
    "segment_impedance_block",
    "solve_coupled",
    "tangential_derivative",
]

_REGION_FLUID = "fluid"
_REGION_SOLID = "solid"
_TWO_PI = 2.0 * np.pi


# ---------------------------------------------------------------------------
# Estructuras de salida
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CoupledDofs:
    """Mapa de grados de libertad del sistema acoplado.

    Atributos
    ---------
    u : ``(N, 2)`` int; ``u[a, 0]``/``u[a, 1]`` son los gdl ``(u_x, u_y)`` del nodo ``a``, o ``-1``
        si el nodo no toca ningún elemento sólido.
    p : ``(N,)`` int; gdl de presión del nodo ``a``, o ``-1`` si no toca fluido.
    n_solid : cantidad de gdl de desplazamiento (2 por nodo sólido).
    n_pressure : cantidad de gdl de presión.
    size : ``n_solid + n_pressure``, el tamaño del sistema.
    """

    u: np.ndarray
    p: np.ndarray
    n_solid: int
    n_pressure: int
    size: int

    @property
    def interface_nodes(self) -> np.ndarray:
        """Nodos con los tres gdl (``u_x, u_y, p``)."""
        return np.flatnonzero((self.u[:, 0] >= 0) & (self.p >= 0))

    @property
    def solid_nodes(self) -> np.ndarray:
        return np.flatnonzero(self.u[:, 0] >= 0)

    @property
    def fluid_nodes(self) -> np.ndarray:
        return np.flatnonzero(self.p >= 0)

    def describe(self) -> str:
        return (
            f"gdl: {self.n_solid} de desplazamiento + {self.n_pressure} de presión "
            f"= {self.size}; nodos sólidos {self.solid_nodes.size}, "
            f"de fluido {self.fluid_nodes.size}, de interfaz {self.interface_nodes.size}"
        )


@dataclass(frozen=True)
class CoupledSolution:
    """Solución del sistema acoplado.

    Atributos
    ---------
    p : ``(N,)`` complejo, presión en los nodos; ``NaN+0j`` donde no hay gdl de presión.
    u : ``(N, 2)`` complejo, desplazamiento; ``NaN`` donde no hay gdl de desplazamiento.
    dof_map : :class:`CoupledDofs`.
    size : tamaño del sistema resuelto (``dof_map.size``).
    frequency : frecuencia [Hz] usada.
    wavenumber : ``k`` complejo del fluido [1/m] a esa frecuencia.
    residual : ``‖A·x − b‖ / ‖b‖`` del sistema **con signo y escala originales** (diagnóstico).
    interface_edges : ``(E, 2)`` aristas de interfaz, orientadas como las recorre el sólido CCW.
    interface_normals : ``(E, 2)`` normales salientes del sólido en esas aristas.
    """

    p: np.ndarray
    u: np.ndarray
    dof_map: CoupledDofs
    size: int
    frequency: float
    wavenumber: complex
    residual: float
    interface_edges: np.ndarray
    interface_normals: np.ndarray

    # alias: el encargo lo llama indistintamente dof_map / dofs
    @property
    def dofs(self) -> CoupledDofs:
        return self.dof_map

    def as_dict(self) -> dict[str, object]:
        return {
            "p": self.p,
            "u": self.u,
            "dof_map": self.dof_map,
            "size": self.size,
            "frequency": self.frequency,
            "residual": self.residual,
        }


# ---------------------------------------------------------------------------
# Validaciones y utilidades de malla
# ---------------------------------------------------------------------------


def _as_nodes(nodes) -> np.ndarray:
    x = np.asarray(nodes, dtype=float)
    if x.ndim != 2 or x.shape[1] != 2:
        raise ValueError(f"nodes debe tener forma (N, 2); recibido {x.shape}")
    if x.shape[0] < 3:
        raise ValueError("se requieren al menos 3 nodos")
    if not np.all(np.isfinite(x)):
        raise ValueError("nodes contiene valores no finitos")
    return x


def _as_triangles(triangles, n_nodes: int) -> np.ndarray:
    t = np.asarray(triangles, dtype=int)
    if t.ndim != 2 or t.shape[1] != 3:
        raise ValueError(f"triangles debe tener forma (M, 3); recibido {t.shape}")
    if t.shape[0] == 0:
        raise ValueError("se requiere al menos un triángulo")
    if np.any(t < 0) or np.any(t >= n_nodes):
        raise ValueError("triangles tiene índices de nodo fuera de rango")
    return t


def _as_regions(regions, n_elements: int) -> np.ndarray:
    reg = np.asarray(regions)
    if reg.ndim != 1:
        raise ValueError(f"regions debe ser un array 1D por elemento; recibido forma {reg.shape}")
    if reg.size != n_elements:
        raise ValueError(
            f"regions debe tener un valor por elemento ({n_elements}); recibido {reg.size}"
        )
    if reg.dtype.kind in "US":
        labels = np.array([str(v) for v in reg], dtype=object)
    else:
        labels = reg
    bad = sorted({str(v) for v in labels if str(v) not in (_REGION_FLUID, _REGION_SOLID)})
    if bad:
        raise ValueError(
            f"regions sólo admite {_REGION_FLUID!r} o {_REGION_SOLID!r}; valores inválidos: {bad}"
        )
    return np.array([str(v) for v in labels])


def _edge_owners(triangles: np.ndarray):
    """Aristas (orientadas por el triángulo CCW), sus claves y el elemento dueño de cada una."""
    t = np.asarray(triangles, dtype=int)
    m = t.shape[0]
    edges = np.vstack([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]])
    owner = np.tile(np.arange(m), 3)
    key = np.sort(edges, axis=1)
    uniq, inv, counts = np.unique(key, axis=0, return_inverse=True, return_counts=True)
    return edges, owner, uniq, np.ravel(inv), counts


def _first_two_occurrences(inv: np.ndarray, n_unique: int) -> tuple[np.ndarray, np.ndarray]:
    """Posiciones de la primera y la segunda ocurrencia de cada arista única."""
    order = np.argsort(inv, kind="stable")
    starts = np.searchsorted(inv[order], np.arange(n_unique))
    first = order[starts]
    second = order[np.minimum(starts + 1, order.size - 1)]
    return first, second


def abc_chain_ends(edges) -> tuple[np.ndarray, np.ndarray]:
    """Extremos libres de las cadenas que forman las aristas de la ABC.

    El término tangencial de ``abc_order=2`` sale de ``∫∂²_s p·v = −∫∂_s p·∂_s v + [∂_s p·v]``.
    En una **cadena cerrada** el corchete se cancela (cada nodo aparece en dos aristas y el término
    se recorre dos veces con signos opuestos); en una cadena **abierta**, los dos extremos dejan un
    término de borde que hay que ensamblar, o de lo contrario la condición de borde no es la que se
    cree estar imponiendo.

    Parameters
    ----------
    edges : ``(B, 2)`` aristas de la ABC.

    Returns
    -------
    (end_nodes, neighbours) : nodos que aparecen en **una sola** arista de la ABC y su vecino en esa
        arista. Con una cadena cerrada los dos arrays salen vacíos.

    Raises
    ------
    ValueError
        Si algún nodo aparece en más de dos aristas (la ABC se ramifica: no es una cadena y
        ``abc_order=2`` no está definido ahí).
    """
    e = np.asarray(edges, dtype=int).reshape(-1, 2)
    if e.shape[0] == 0:
        return np.zeros(0, dtype=int), np.zeros(0, dtype=int)
    nodes, counts = np.unique(e.reshape(-1), return_counts=True)
    if np.any(counts > 2):
        bad = nodes[counts > 2][:3]
        raise ValueError(
            "abc_order=2 necesita que las aristas de la ABC formen cadenas simples (cada nodo en "
            f"una o dos aristas); los nodos {bad.tolist()} aparecen en más de dos aristas"
        )
    ends = nodes[counts == 1]
    if ends.size == 0:
        return np.zeros(0, dtype=int), np.zeros(0, dtype=int)
    neighbours = np.empty(ends.size, dtype=int)
    for k, nd in enumerate(ends):
        row = int(np.flatnonzero(e == nd)[0] // 2)
        a, b = int(e[row, 0]), int(e[row, 1])
        neighbours[k] = b if a == int(nd) else a
    return ends.astype(int), neighbours


def _chain_local(lengths: np.ndarray, scale: complex) -> np.ndarray:
    """Bloque local (B, 2, 2) de ``scale·∫_e ∂_s φ_i·∂_s φ_j ds`` con ``φ`` P1 en la arista."""
    s = (np.asarray(scale, dtype=complex) / np.asarray(lengths, dtype=float))[:, None, None]
    local = np.empty((s.shape[0], 2, 2), dtype=complex)
    local[:, 0, 0] = s[:, 0, 0]
    local[:, 0, 1] = -s[:, 0, 0]
    local[:, 1, 0] = -s[:, 0, 0]
    local[:, 1, 1] = s[:, 0, 0]
    return local


def tangential_derivative(p_inc, points, tangent, step) -> np.ndarray:
    """``∂_s p_i = ∇p_i·t`` en ``points`` con ``t`` unitaria (``(M,)``).

    Usa ``p_inc.gradient`` si existe (exacto para :class:`fem.helmholtz2d.PlaneWave`) y, si no, una
    diferencia central de ``p_inc.value`` a lo largo de la tangente con paso ``step`` (el largo de la
    arista), que para un campo suave es de segundo orden.
    """
    pts = np.atleast_2d(np.asarray(points, dtype=float))
    tan = np.atleast_2d(np.asarray(tangent, dtype=float))
    grad = getattr(p_inc, "gradient", None)
    if grad is not None:
        return np.einsum("ij,ij->i", np.asarray(grad(pts), dtype=complex), tan)
    d = 0.5 * np.asarray(step, dtype=float).reshape(-1, 1)
    return (
        np.asarray(p_inc.value(pts + d * tan), dtype=complex).reshape(-1)
        - np.asarray(p_inc.value(pts - d * tan), dtype=complex).reshape(-1)
    ) / (2.0 * d.reshape(-1))


def _as_abc_center(abc_center) -> np.ndarray:
    """Valida y normaliza ``abc_center`` (obligatorio en BGT-2)."""
    if abc_center is None:
        raise ValueError(
            "abc_order=2 (BGT-2) necesita abc_center=(x0, y0), el centro de curvatura del "
            "contorno de la ABC: el término 1/(2r) sale de r = |x_mid − abc_center|. Sin él no hay "
            "forma de saber el radio y la condición quedaría mal sin avisar"
        )
    center = np.asarray(abc_center, dtype=float).reshape(-1)
    if center.size != 2 or not np.all(np.isfinite(center)):
        raise ValueError(f"abc_center debe ser un punto finito (2,); recibido {abc_center!r}")
    return center


def _bgt2_coefficients(nodes, edges, k, abc_center) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(coef, ct, lengths)`` de la ABC de orden 2 por arista.

    ``coef = i·k − 1/(2r) + i/(8k r²)`` es el factor de masa (``∂p_s/∂n = coef·p_s + (i/2k)∂²_s p_s``)
    y ``ct = i/(2k)`` el del término tangencial; ``r`` se evalúa en el punto medio de la arista
    respecto de ``abc_center``.
    """
    x = np.asarray(nodes, dtype=float)
    e = np.asarray(edges, dtype=int).reshape(-1, 2)
    lengths = np.linalg.norm(x[e[:, 1]] - x[e[:, 0]], axis=1)
    if np.any(lengths <= 0.0):
        raise ValueError("abc_order=2: hay aristas de la ABC de longitud nula")
    center = _as_abc_center(abc_center)
    r_nodes = np.linalg.norm(x[np.unique(e)] - center[None, :], axis=1)
    r_e = np.linalg.norm(edge_midpoints(x, e) - center[None, :], axis=1)
    if np.any(r_e <= 0.0) or np.any(r_nodes <= 0.0):
        raise ValueError(
            "abc_order=2: hay aristas de la ABC con r = |x − abc_center| = 0 (el centro de "
            "curvatura cae sobre el contorno); el término 1/(2r) diverge"
        )
    coef = 1j * k - 0.5 / r_e + 0.125j / (k * r_e**2)
    return coef, 0.5j / k, lengths


def abc_boundary_block(nodes, edges, k, rho=1.0, *, abc_order: int = 1, abc_center=None):
    """Bloque de la condición absorbente en la numeración de **nodos** (``N × N``).

    Es el término que la ABC aporta al lado izquierdo, ``−∮ (1/ρ)·B(p)·v ds``, con

    - ``abc_order=1``: ``B(p) = i·k·p`` (impedancia),
    - ``abc_order=2``: ``B(p) = coef·p + (i/2k)·∂²_s p`` con
      ``coef = i·k − 1/(2r) + i/(8k r²)`` (BGT-2, contorno circular de radio ``r``).

    El término tangencial se ensambla con ``∫∂²_s p·v ds = −∫∂_s p·∂_s v ds + [∂_s p·v]``: el
    corchete se cancela en una cadena cerrada de aristas y se agrega en los dos extremos de una
    cadena abierta (ver :func:`abc_chain_ends`).

    Se expone aparte de :func:`solve_coupled` para poder auditar el operador por separado (su
    respuesta a un modo ``exp(i·m·θ)`` tiene que ser el DtN exacto ``k·H_m′(kR)/H_m(kR)``).

    Parameters
    ----------
    nodes : ``(N, 2)``; ``edges`` : ``(B, 2)`` aristas de la ABC; ``k`` : número de onda complejo
    del fluido; ``rho`` : densidad del fluido (el bloque lleva el peso ``1/ρ``).
    """
    x = np.asarray(nodes, dtype=float)
    e = np.asarray(edges, dtype=int).reshape(-1, 2)
    n = x.shape[0]
    k = complex(k)
    lengths = np.linalg.norm(x[e[:, 1]] - x[e[:, 0]], axis=1)
    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    vals: list[np.ndarray] = []

    def push(r, c, v):
        rows.append(np.asarray(r, dtype=int).ravel())
        cols.append(np.asarray(c, dtype=int).ravel())
        vals.append(np.asarray(v, dtype=complex).ravel())

    if abc_order == 1:
        local = (-1j * k / float(rho)) * weights_to_local(lengths)
        push(np.repeat(e, 2, axis=1), np.tile(e, (1, 2)), local.reshape(-1))
    elif abc_order == 2:
        coef, ct, lengths = _bgt2_coefficients(x, e, k, abc_center)
        ct = ct / float(rho)
        local = (-coef / float(rho))[:, None, None] * weights_to_local(lengths) + _chain_local(
            lengths, ct
        )
        push(np.repeat(e, 2, axis=1), np.tile(e, (1, 2)), local.reshape(-1))
        end, nb = abc_chain_ends(e)
        if end.size:
            end_len = np.linalg.norm(x[nb] - x[end], axis=1)
            push(
                np.repeat(end[:, None], 2, axis=1),
                np.column_stack([end, nb]),
                np.column_stack([-ct / end_len, ct / end_len]),
            )
    else:
        raise ValueError(f"abc_order debe ser 1 o 2; recibido {abc_order!r}")

    return sp.coo_matrix(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
        shape=(n, n),
        dtype=complex,
    ).tocsr()


def abc_incident_rhs(nodes, edges, p_inc, k, rho=1.0, *, abc_order: int = 1, abc_center=None):
    """Término incidente de la ABC en la numeración de nodos (``(N,)`` complejo).

    Es ``b += ∮ (1/ρ)·(∂p_i/∂n − B(p_i))·v ds``, la parte conocida de la condición
    ``∂p/∂n = ∂p_i/∂n + B(p − p_i)``.

    - ``abc_order=1``: sólo en las caras iluminadas (``n·d < 0``), como en
      `fem.helmholtz2d.assemble`. En las de sombra el campo total ya es saliente y el término
      incidente metería una fuente espuria (D11/D12).
    - ``abc_order=2``: en **todas** las aristas. La condición se impone sobre ``p_s`` de forma
      consistente y ``∂p_i/∂n − B(p_i)`` —incluido el ``∂²_s p_i``, que se integra por partes y sólo
      necesita ``∂_s p_i``— es exacta para la onda plana en cualquier cara, así que no queda fuente
      espuria en la sombra.
    """
    x = np.asarray(nodes, dtype=float)
    e = np.asarray(edges, dtype=int).reshape(-1, 2)
    b = np.zeros(x.shape[0], dtype=complex)
    if e.shape[0] == 0:
        return b
    k = complex(k)
    lengths = np.linalg.norm(x[e[:, 1]] - x[e[:, 0]], axis=1)
    mid = edge_midpoints(x, e)
    nrm = outward_normals(x, e)

    if abc_order == 1:
        direction = getattr(p_inc, "direction", None)
        lit = (
            np.ones(e.shape[0], dtype=bool)
            if direction is None
            else (nrm @ np.asarray(direction, dtype=float)) < 0.0
        )
        f = p_inc.normal_derivative(mid, nrm) - 1j * k * p_inc.value(mid)
        contrib = np.where(lit, (f * lengths / 2.0) / float(rho), 0.0)
        np.add.at(b, e[:, 0], contrib)
        np.add.at(b, e[:, 1], contrib)
    elif abc_order == 2:
        coef, ct, lengths = _bgt2_coefficients(x, e, k, abc_center)
        ct = ct / float(rho)
        val = (p_inc.normal_derivative(mid, nrm) - coef * p_inc.value(mid)) * (
            lengths / 2.0
        ) / float(rho)
        np.add.at(b, e[:, 0], val)
        np.add.at(b, e[:, 1], val)
        # +ct·∫∂_s p_i·∂_s v ds = ct·∂_s p_i·(v_b − v_a)
        tangent = (x[e[:, 1]] - x[e[:, 0]]) / lengths[:, None]
        ds_inc = tangential_derivative(p_inc, mid, tangent, lengths)
        np.add.at(b, e[:, 0], -ct * ds_inc)
        np.add.at(b, e[:, 1], +ct * ds_inc)
        # extremos de cadenas abiertas: +ct·∂_u p_i (u = del extremo hacia su vecino)
        end, nb = abc_chain_ends(e)
        if end.size:
            end_len = np.linalg.norm(x[nb] - x[end], axis=1)
            t_end = (x[nb] - x[end]) / end_len[:, None]
            ds_end = tangential_derivative(p_inc, 0.5 * (x[end] + x[nb]), t_end, end_len)
            np.add.at(b, end, ct * ds_end)
    else:
        raise ValueError(f"abc_order debe ser 1 o 2; recibido {abc_order!r}")
    return b


def segment_impedance_block(nodes, edges, impedance, *, frequency, rho=1.0):
    """Bloque de presión de un **segmento con impedancia distribuida** ``p = Z·(v·n)``.

    Es la condición de una membrana delgada tratada como impedancia local (masa-resorte-resistencia
    por unidad de área), escrita de forma que reusa exactamente el operador de la ABC-1. Con
    ``e^{−iωt}`` y ``v_n = (1/(i·ω·ρ_f))·∂p/∂n`` (``n`` **saliente del fluido**, o sea hacia la
    membrana):

        p = Z·v_n   ⇒   ∂p/∂n = (i·ω·ρ_f/Z)·p ≡ i·k_w·p ,   k_w = ω·ρ_f/Z   [1/m]

    y el término débil es el mismo de siempre, ``−∫ (i·k_w/ρ_f)·p·q ds``, que es lo que devuelve
    :func:`abc_boundary_block` con ``k = k_w``. No hay una fórmula nueva que auditar: la identidad
    con la ABC-1 en el caso ``Z = ρ_f·c_f`` (``k_w = k_f``) es exacta término a término.

    **Límites (los dos son casos conocidos y por eso sirven de test)**: ``Z → ∞`` ⇒ ``k_w → 0`` ⇒ el
    término se anula y queda ``∂p/∂n = 0`` (pared rígida); ``Z → 0`` ⇒ ``k_w → ∞`` ⇒ ``p → 0``
    (pared blanda / presión nula). Los dos se alcanzan numéricamente porque ``Z`` sólo entra en el
    coeficiente, no en la incógnita.

    **Qué NO es**: no es la condición de salto de una placa elástica entre dos fluidos
    (``p⁺ − p⁻ = Z·v_n`` con ``v_n`` continuo). En una línea interior compartida por dos dominios,
    aplicar esto sobre las dos caras equivale a dos membranas locales independientes —el modelo
    "locally reacting" bilateral, que es el que da los dos límites pedidos: en el salto, ``Z → 0``
    daría transparencia, no ``p = 0``.

    Parameters
    ----------
    nodes : ``(N, 2)`` coordenadas [m]; ``edges : (B, 2)`` aristas del segmento (interiores o no).
    impedance : ``Z`` [Pa·s/m], complejo, ``Re Z > 0`` para que absorba.
    frequency : frecuencia [Hz].
    rho : densidad del fluido [kg/m³].

    Returns
    -------
    ``scipy.sparse`` CSR ``(N, N)`` compleja, lista para ensamblar en el bloque de presión. La
    potencia real que el segmento absorbe es ``½·∫|p|²·Re(1/Z) ds`` (positiva si ``Re Z > 0``),
    :func:`segment_absorbed_power`.
    """
    z = complex(impedance)
    if not (np.isfinite(z.real) and np.isfinite(z.imag)):
        raise ValueError(f"impedance debe ser finita; recibido {z!r}")
    if z.real <= 0.0:
        raise ValueError(
            "impedance necesita Re(Z) > 0 para que el segmento absorba potencia real "
            f"(½∫|p|²·Re(1/Z) ds); recibido {z!r}"
        )
    k_w = _TWO_PI * float(frequency) * float(rho) / z
    return abc_boundary_block(nodes, edges, k_w, rho=rho, abc_order=1)


def segment_absorbed_power(nodes, p, edges, impedance) -> float:
    """Potencia real **[W/m]** absorbida por un segmento con ``p = Z·(v·n)``.

    ``P = ½·∫ Re(p·v_n*) ds = ½·Re(1/Z)·∫|p|² ds``, con la masa P1 de arista
    (``L/6·[[2,1],[1,2]]``). Es cuadrática en ``p`` y no depende de diferencias finitas, así que es
    la medida **exacta** de lo que el ensamblado de :func:`segment_impedance_block` disipa. Positivo
    = potencia saliendo del fluido hacia la membrana (absorbida) cuando ``Re Z > 0``; con
    ``Re Z ≤ 0`` no tiene sentido físico y se rechaza.

    En una línea interior con dos caras hay que sumar **las dos** (cada una con su ``p`` local): la
    membrana recibe fluido por los dos lados.
    """
    z = complex(impedance)
    if not (np.isfinite(z.real) and np.isfinite(z.imag)):
        raise ValueError(f"impedance debe ser finita; recibido {z!r}")
    if z.real <= 0.0:
        raise ValueError(f"impedance necesita Re(Z) > 0; recibido {z!r}")
    x = np.asarray(nodes, dtype=float)
    e = np.asarray(edges, dtype=int).reshape(-1, 2)
    if e.shape[0] == 0:
        return 0.0
    lengths = np.linalg.norm(x[e[:, 1]] - x[e[:, 0]], axis=1)
    m_loc = weights_to_local(lengths)
    pe = np.asarray(p, dtype=complex)[e]
    dens = np.einsum("ei,ej,eij->e", pe, np.conj(pe), m_loc).real
    return float(0.5 * (1.0 / z).real * np.sum(dens))


def interface_edges(nodes, triangles, regions) -> tuple[np.ndarray, np.ndarray]:
    """Aristas interiores con un elemento de cada lado y normal saliente del sólido.

    Parameters
    ----------
    nodes : ``(N, 2)``.
    triangles : ``(M, 3)`` CCW.
    regions : ``(M,)`` con ``'fluid'``/``'solid'``.

    Returns
    -------
    (edges, normals) : ``(E, 2)`` aristas orientadas como las recorre el triángulo **sólido** CCW
        (``normals = (t_y, −t_x)/|t|`` es la saliente del sólido) y ``(E, 2)`` normales unitarias.
    """
    x = _as_nodes(nodes)
    t = _as_triangles(triangles, x.shape[0])
    reg = _as_regions(regions, t.shape[0])
    solid_e = reg == _REGION_SOLID
    if not np.any(solid_e):
        raise ValueError("regions no tiene ningún elemento 'solid'")
    if np.all(solid_e):
        raise ValueError("regions no tiene ningún elemento 'fluid'")

    edges, owner, _, inv, counts = _edge_owners(t)
    if np.any(counts > 2):
        bad = int(np.flatnonzero(counts > 2)[0])
        raise ValueError(
            f"malla no conforme: la arista {bad} pertenece a {counts[bad]} elementos; "
            "la interfaz se define sólo sobre aristas compartidas por exactamente dos triángulos"
        )
    first, second = _first_two_occurrences(inv, counts.size)
    interior = counts == 2
    o1, o2 = owner[first], owner[second]
    mixed = interior & (solid_e[np.minimum(o1, o2)] != solid_e[np.maximum(o1, o2)])
    solid_occ = np.where(solid_e[np.minimum(o1, o2)], first, second)[mixed]
    iface = edges[solid_occ]
    tangent = x[iface[:, 1]] - x[iface[:, 0]]
    normal = np.column_stack([tangent[:, 1], -tangent[:, 0]])
    lengths = np.linalg.norm(normal, axis=1)
    if np.any(lengths == 0.0):
        raise ValueError("arista de interfaz de longitud nula")
    return iface, normal / lengths[:, None]


def _boundary_owner_map(triangles: np.ndarray) -> dict[tuple[int, int], tuple[int, int, int]]:
    """``clave ordenada → (i, j, elemento dueño)`` con ``(i, j)`` orientada por el triángulo CCW."""
    edges, owner, _, inv, counts = _edge_owners(triangles)
    first, _ = _first_two_occurrences(inv, counts.size)
    out: dict[tuple[int, int], tuple[int, int, int]] = {}
    for u_i in np.flatnonzero(counts == 1):
        i, j = (int(v) for v in edges[first[u_i]])
        out[(min(i, j), max(i, j))] = (i, j, int(owner[first[u_i]]))
    return out


class _UnionFind:
    """Union-find mínimo para identificar nodos periódicos."""

    def __init__(self, n: int) -> None:
        self.parent = np.arange(n)

    def find(self, a: int) -> int:
        p = self.parent
        root = a
        while p[root] != root:
            root = p[root]
        while p[a] != root:  # compresión de camino
            p[a], a = root, p[a]
        return root

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[max(ra, rb)] = min(ra, rb)


def _resolve_frequency(frequency, p_inc, fluid: Material) -> float:
    """Frecuencia de trabajo: la explícita, o la deducida de ``p_inc.k``."""
    k_inc = getattr(p_inc, "k", None) if p_inc is not None else None
    if frequency is None:
        if k_inc is None:
            raise ValueError(
                "solve_coupled necesita frequency (o un p_inc con atributo .k del que deducirla, "
                "p. ej. fem.helmholtz2d.PlaneWave)"
            )
        frequency = complex(k_inc).real * fluid.c / _TWO_PI
    frequency = float(frequency)
    if not (frequency > 0.0):
        raise ValueError(f"frequency debe ser > 0; recibido {frequency!r}")
    if k_inc is not None:
        k_fluid = fluid.wavenumber(frequency)
        if not np.isclose(complex(k_inc), k_fluid, rtol=1e-6, atol=0.0):
            raise ValueError(
                f"el k de p_inc ({complex(k_inc)!r}) no corresponde a {fluid.name} a "
                f"{frequency:g} Hz ({k_fluid!r}); revisá la frecuencia o el material"
            )
    return frequency


# ---------------------------------------------------------------------------
# Ensamblado y resolución
# ---------------------------------------------------------------------------


def solve_coupled(
    nodes,
    triangles,
    regions,
    fluid,
    solid,
    abc_edges=None,
    p_inc=None,
    abc_order: int = 1,
    *,
    frequency=None,
    dirichlet=None,
    periodic_pairs=None,
    abc_center=None,
    robin_edges=None,
    robin_impedance=None,
) -> CoupledSolution:
    """Resuelve el sistema acoplado fluido/sólido 2D en plano deformación.

    Parameters
    ----------
    nodes : ``(N, 2)`` coordenadas.
    triangles : ``(M, 3)`` triángulos CCW.
    regions : ``(M,)`` ``'fluid'`` o ``'solid'`` por elemento (la placa es una **banda** de
        elementos sólidos; el resto, fluido).
    fluid : :class:`materials.Material` (acústico homogéneo).
    solid : :class:`materials.ElasticSolid` (isótropo, sin pérdidas).
    abc_edges : ``(B, 2)`` aristas de la ABC; ``None`` ⇒ toda la frontera de los elementos de
        fluido. Las aristas se reorientan contra la malla (el orden de los nodos da igual) y se
        valida que pertenezcan a elementos de fluido.
    p_inc : onda incidente (``.value(points)`` y ``.normal_derivative(points, normals)``), típicamente
        ``fem.helmholtz2d.plane_wave(...)``. Entra como dato de borde en las caras iluminadas
        (``n·d < 0``) de la ABC, sobre el campo total: ``b += ∮ (1/ρ)·(∂p_i/∂n − i·k·p_i)·q ds``.
        ``None`` ⇒ sin fuente (campo saliente).
    abc_order : orden de la condición absorbente en ``abc_edges``.

        - ``1`` (default, comportamiento de siempre): impedancia ``∂p_s/∂n = i·k·p_s``. El término
          incidente entra **sólo** en las caras iluminadas (``n·d < 0``). Es exacta en incidencia
          normal sobre una cara plana (ducto/placa) y tiene un piso de reflexión en contorno curvo
          (D12): con la cápsula de `cochlea` da 12 % L2 y 17–21 % en la cara de sombra, y **no**
          mejora al agrandar el dominio.
        - ``2`` (BGT-2 de Bayliss–Turkel–Gunzburger, el mismo de
          `fem.helmholtz2d.assemble_scattered`): sobre un contorno circular de radio ``r``

              ∂p_s/∂n = (i·k − 1/(2r) + i/(8k r²))·p_s + (i/(2k))·∂²_s p_s

          El término ``1/(2r)`` es la curvatura del frente cilíndrico y el ``∂²_s`` la parte
          angular; sin ellos la ABC-1 refleja la onda cilíndrica y ese error no baja con ``ppw``.
          La forma débil usa ``∫∂²_s p_s·v ds = −∫∂_s p_s·∂_s v ds + [∂_s p_s·v]``: el corchete se
          cancela si la ABC es una **cadena cerrada** (círculo completo) y se ensambla aparte en los
          dos extremos si es una **cadena abierta** (las dos caras de un ducto). Un nodo de la ABC
          que aparezca en más de dos aristas (ABC ramificada) se rechaza.
          **Requiere ``abc_center``** y sólo es correcto si la ABC es (aproximadamente) una
          circunferencia centrada ahí. Además, con ``abc_order=2`` el término incidente entra en
          **todas** las aristas de la ABC, no sólo en las iluminadas: la condición se impone sobre
          el campo dispersado de forma consistente, y ``∂p_i/∂n − B(p_i)`` (con el ``∂²_s p_i``
          calculado con ``∇p_i``) es exacto para la onda plana en cualquier cara. Eso es lo que
          elimina la fuente espuria ``i·k·p_i·(d·n − 1)`` que la máscara deja en la cara de sombra
          de un contorno curvo (D12) — y por eso el orden 2 **sí** baja el error de sombra en la
          cápsula mientras la ABC-1 no.
    abc_center : ``(2,)`` centro de curvatura del contorno de la ABC (**obligatorio con
        ``abc_order=2``**, se ignora con ``abc_order=1``). De él sale ``r_e = |x_mid − abc_center|``
        por arista. Qué pasa si se pasa mal:

        - ``None`` con ``abc_order=2`` ⇒ ``ValueError`` (no hay default: el radio es específico de
          la geometría y adivinarlo daría una condición silenciosamente equivocada).
        - centro **desplazado** (no es el centro del círculo) ⇒ no falla, pero ``r_e`` deja de ser
          el radio local y el término de curvatura ``−1/(2r_e) + i/(8k r_e²)`` queda mal: el
          resultado es peor que la ABC-1 en la zona donde el error de ``r_e`` es grande. Se degrada
          en silencio, así que conviene validar el centro (``|mid − c|`` constante).
        - centro **sobre el contorno** (alguna ``r_e = 0``) ⇒ ``ValueError`` (``1/r_e`` diverge).
        - contorno **no circular**: BGT-2 sigue siendo una condición de radiación de orden 2 pero
          con el radio equivocado; en una cara plana el límite correcto es ``r → ∞`` (allí la
          curvatura es cero y queda ``ik + (i/2k)∂²_s``, que coincide con la ABC-2 de una línea).
    frequency : frecuencia [Hz]. Obligatoria salvo que ``p_inc`` traiga ``.k``.
    dirichlet : ``{'ux': nodos, 'uy': nodos, 'p': nodos}`` con gdl esenciales nulos. El rodillo
        ``uy = 0`` en las paredes de la placa es lo que hace que la capa sólida del benchmark tenga
        velocidad ``c_p`` (ver el docstring del módulo).
    periodic_pairs : ``(K, 2)`` pares de nodos que se identifican (malla periódica). Comparten los
        tres gdl; sirve para incidencia oblicua sin paredes (ver el docstring del módulo).
    robin_edges : ``(B, 2)`` aristas de un **segmento interior** con impedancia distribuida
        ``p = Z·(v·n)`` (membrana delgada, p. ej. la basilar). Se ensambla el mismo bloque de borde
        de :func:`segment_impedance_block` —``−∫(i·k_w/ρ_f)·p·q ds`` con ``k_w = ω·ρ_f/Z``— en el
        bloque de presión, sin tocar el acoplamiento. A diferencia de ``abc_edges`` **no** tienen
        que ser frontera de la malla: la línea puede ser interior (en una malla partida, las dos
        caras coincidentes son nodos distintos y cada una lleva su término). Requiere
        ``robin_impedance``.
    robin_impedance : ``Z`` [Pa·s/m] del segmento; complejo con ``Re Z > 0`` para que absorba.
        ``Z → ∞`` reproduce la pared rígida (el término se anula) y ``Z → 0`` la de presión nula
        (``k_w → ∞`` fija ``p = 0`` en los nodos del segmento).

    Returns
    -------
    :class:`CoupledSolution` con ``p``, ``u``, ``dof_map``, ``size`` y diagnósticos.
    """
    x = _as_nodes(nodes)
    t = _as_triangles(triangles, x.shape[0])
    reg = _as_regions(regions, t.shape[0])
    n = x.shape[0]

    if not isinstance(fluid, Material):
        raise TypeError(f"fluid debe ser materials.Material; recibido {type(fluid).__name__}")
    if not isinstance(solid, ElasticSolid):
        raise TypeError(
            "solid debe ser materials.ElasticSolid (no se define otra clase de sólido); "
            f"recibido {type(solid).__name__}"
        )
    if abc_order not in (1, 2):
        raise ValueError(
            f"abc_order={abc_order!r} no está soportado: sólo 1 (impedancia) o 2 "
            "(Bayliss–Turkel–Gunzburger, que necesita abc_center)"
        )
    center_abc = _as_abc_center(abc_center) if abc_order == 2 else None

    frequency = _resolve_frequency(frequency, p_inc, fluid)
    omega = _TWO_PI * frequency
    k_f = fluid.wavenumber(frequency)

    solid_e = reg == _REGION_SOLID
    fluid_e = ~solid_e
    solid_tris = t[solid_e]
    fluid_tris = t[fluid_e]
    if solid_tris.shape[0] == 0 or fluid_tris.shape[0] == 0:
        raise ValueError(
            "solve_coupled requiere al menos un elemento de cada tipo "
            f"(solid: {solid_tris.shape[0]}, fluid: {fluid_tris.shape[0]})"
        )

    # ---- nodos canónicos (identificación periódica) ------------------------
    solid_node = np.zeros(n, dtype=bool)
    solid_node[solid_tris.ravel()] = True
    fluid_node = np.zeros(n, dtype=bool)
    fluid_node[fluid_tris.ravel()] = True

    uf = _UnionFind(n)
    if periodic_pairs is not None:
        pairs = np.asarray(periodic_pairs, dtype=int).reshape(-1, 2)
        if np.any(pairs < 0) or np.any(pairs >= n):
            raise ValueError("periodic_pairs tiene índices de nodo fuera de rango")
        if np.any(pairs[:, 0] == pairs[:, 1]):
            raise ValueError("periodic_pairs tiene un par con los dos nodos iguales")
        # los dos nodos de un par tienen que tener la misma disponibilidad de gdl
        bad = (solid_node[pairs[:, 0]] != solid_node[pairs[:, 1]]) | (
            fluid_node[pairs[:, 0]] != fluid_node[pairs[:, 1]]
        )
        if np.any(bad):
            a, b = pairs[bad][0]
            raise ValueError(
                f"periodic_pairs une los nodos {int(a)} y {int(b)}, que no tienen los mismos gdl "
                "(uno toca sólido y el otro no, o uno toca fluido y el otro no); los pares "
                "periódicos tienen que unir nodos del mismo tipo"
            )
        for a, b in pairs:
            uf.union(int(a), int(b))
    rep = np.array([uf.find(i) for i in range(n)], dtype=int)

    solid_rep = np.zeros(n, dtype=bool)
    solid_rep[rep[solid_node]] = True
    fluid_rep = np.zeros(n, dtype=bool)
    fluid_rep[rep[fluid_node]] = True

    reps_u = np.flatnonzero(solid_rep)
    reps_p = np.flatnonzero(fluid_rep)
    n_u = 2 * reps_u.size
    n_p = reps_p.size
    size = n_u + n_p

    u_dof_rep = -np.ones((n, 2), dtype=int)
    u_dof_rep[reps_u, 0] = 2 * np.arange(reps_u.size)
    u_dof_rep[reps_u, 1] = 2 * np.arange(reps_u.size) + 1
    p_dof_rep = -np.ones(n, dtype=int)
    p_dof_rep[reps_p] = n_u + np.arange(reps_p.size)
    u_dof = u_dof_rep[rep]
    p_dof = p_dof_rep[rep]
    dofs = CoupledDofs(u=u_dof, p=p_dof, n_solid=n_u, n_pressure=n_p, size=size)

    # ---- bloques -----------------------------------------------------------
    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    vals: list[np.ndarray] = []

    def _add(rdofs, cdofs, local):
        n_rows = rdofs.shape[1]
        n_cols = cdofs.shape[1]
        rows.append(np.repeat(rdofs, n_cols, axis=1).ravel())
        cols.append(np.tile(cdofs, (1, n_rows)).ravel())
        vals.append(np.asarray(local).reshape(-1))

    def _add_triplets(rdofs, cdofs, values):
        rows.append(np.asarray(rdofs, dtype=int).ravel())
        cols.append(np.asarray(cdofs, dtype=int).ravel())
        vals.append(np.asarray(values, dtype=complex).ravel())

    # (1,1) sólido: K_s − ρ_s·ω²·M_s
    k_s, m_s, _ = p1_elastic_local(x, solid_tris, solid)
    _add(
        u_dof[solid_tris].reshape(-1, 6),
        u_dof[solid_tris].reshape(-1, 6),
        k_s - solid.rho * omega**2 * m_s,
    )

    # (2,2) fluido: (K_f − k²·M_f)/ρ_f
    m_f, kf, _ = p1_tri_local(x, fluid_tris)
    _add(
        p_dof[fluid_tris],
        p_dof[fluid_tris],
        (kf - k_f**2 * m_f) / fluid.rho,
    )

    # (1,2) y (2,1) acoplamiento
    iface, iface_n = interface_edges(x, t, reg)
    if iface.shape[0] == 0:
        raise ValueError(
            "no hay aristas de interfaz: regions debe tener al menos un elemento de cada tipo "
            "compartiendo una arista"
        )
    a_i = iface[:, 0]
    b_i = iface[:, 1]
    lengths = np.linalg.norm(x[b_i] - x[a_i], axis=1)
    ud = u_dof[iface].reshape(-1, 4)  # [ux_a, uy_a, ux_b, uy_b]
    pd = p_dof[iface]  # (E, 2)
    if np.any(ud < 0) or np.any(pd < 0):
        raise RuntimeError(
            "un nodo de interfaz no tiene los tres gdl; revisá la marcación de regions"
        )
    m_loc = weights_to_local(lengths)  # ∫ψ_i ψ_j ds = L/6·[[2,1],[1,2]]
    c_local = np.empty((iface.shape[0], 4, 2), dtype=float)
    nx_e = iface_n[:, 0][:, None, None]
    ny_e = iface_n[:, 1][:, None, None]
    c_local[:, 0, :] = (nx_e * m_loc)[:, 0, :]
    c_local[:, 1, :] = (ny_e * m_loc)[:, 0, :]
    c_local[:, 2, :] = (nx_e * m_loc)[:, 1, :]
    c_local[:, 3, :] = (ny_e * m_loc)[:, 1, :]
    _add(ud, pd, c_local)                     # +C
    _add(pd, ud, omega**2 * np.transpose(c_local, (0, 2, 1)))  # +ω²·Cᵀ

    b_full = np.zeros(size, dtype=complex)

    # ---- ABC y término incidente ------------------------------------------
    bmap = _boundary_owner_map(t)

    if abc_edges is None:
        keep = [entry[:2] for entry in bmap.values() if fluid_e[entry[2]]]
        abc = np.array(keep, dtype=int).reshape(-1, 2)
    else:
        raw = np.asarray(abc_edges, dtype=int).reshape(-1, 2)
        keep = []
        for row, (i, j) in enumerate(raw):
            if int(i) == int(j):
                raise ValueError(f"abc_edges: la arista {row} tiene sus dos nodos iguales")
            entry = bmap.get((min(int(i), int(j)), max(int(i), int(j))))
            if entry is None:
                raise ValueError(
                    f"abc_edges: la arista {row} = ({int(i)}, {int(j)}) no es de frontera de la "
                    "malla; la ABC (y el signo del término incidente) sólo están definidos ahí"
                )
            if not fluid_e[entry[2]]:
                raise ValueError(
                    f"abc_edges: la arista {row} pertenece a un elemento 'solid'; la ABC de "
                    "impedancia es una condición del fluido acústico"
                )
            keep.append(entry[:2])
        abc = np.array(keep, dtype=int).reshape(-1, 2)

    if abc.shape[0]:
        blk = abc_boundary_block(
            x, abc, k_f, rho=fluid.rho, abc_order=abc_order, abc_center=center_abc
        )
        coo = blk.tocoo()
        row_dof = p_dof[coo.row]
        col_dof = p_dof[coo.col]
        keep = (row_dof >= 0) & (col_dof >= 0)
        if not np.all(keep):
            raise RuntimeError(
                "la ABC toca nodos sin gdl de presión; revisá que abc_edges no incluya nodos "
                "exclusivamente sólidos"
            )
        _add_triplets(row_dof, col_dof, coo.data)

        if p_inc is not None:
            rhs = abc_incident_rhs(
                x, abc, p_inc, k_f, rho=fluid.rho, abc_order=abc_order, abc_center=center_abc
            )
            has = p_dof >= 0
            np.add.at(b_full, p_dof[has], rhs[has])

    # ---- segmento interior con impedancia distribuida (membrana delgada) ----
    if robin_edges is not None:
        if robin_impedance is None:
            raise ValueError(
                "robin_edges requiere robin_impedance (Z [Pa·s/m] del segmento); sin Z no hay "
                "condición que imponer"
            )
        if periodic_pairs is not None:
            raise ValueError(
                "robin_edges con periodic_pairs no está soportado: los nodos identificados "
                "comparten gdl de presión y el término de borde se sumaría dos veces sobre la misma "
                "incógnita"
            )
        redges = np.asarray(robin_edges, dtype=int).reshape(-1, 2)
        if redges.size and (np.any(redges < 0) or np.any(redges >= n)):
            raise ValueError("robin_edges tiene índices de nodo fuera de rango")
        if np.any(redges[:, 0] == redges[:, 1]):
            raise ValueError("robin_edges tiene una arista con sus dos nodos iguales")
        if redges.shape[0]:
            keys = np.sort(redges, axis=1)
            if np.unique(keys, axis=0).shape[0] != keys.shape[0]:
                raise ValueError(
                    "robin_edges repite una arista: el término se ensamblaría dos veces y la "
                    "impedancia efectiva sería Z/2"
                )
            if abc.shape[0] and np.intersect1d(keys, np.sort(abc, axis=1)).size:
                raise ValueError(
                    "una arista está en abc_edges y en robin_edges a la vez: la ABC y la "
                    "impedancia del segmento son condiciones distintas sobre la misma superficie"
                )
            rlen = np.linalg.norm(x[redges[:, 1]] - x[redges[:, 0]], axis=1)
            if np.any(rlen <= 0.0):
                raise ValueError("robin_edges tiene una arista de longitud nula")
            blk = segment_impedance_block(
                x, redges, robin_impedance, frequency=frequency, rho=fluid.rho
            )
            coo = blk.tocoo()
            row_dof = p_dof[coo.row]
            col_dof = p_dof[coo.col]
            if np.any(row_dof < 0) or np.any(col_dof < 0):
                raise RuntimeError(
                    "el segmento con impedancia toca nodos sin gdl de presión; la membrana sólo "
                    "puede apoyarse en nodos de fluido (o de interfaz)"
                )
            _add_triplets(row_dof, col_dof, coo.data)

    a_full = sp.coo_matrix(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
        shape=(size, size),
        dtype=complex,
    ).tocsr()

    # ---- gdl esenciales nulos ---------------------------------------------
    fixed: list[np.ndarray] = []
    if dirichlet is not None:
        if not isinstance(dirichlet, dict):
            raise TypeError("dirichlet debe ser un dict {'ux': nodos, 'uy': nodos, 'p': nodos}")
        unknown = set(dirichlet) - {"ux", "uy", "p"}
        if unknown:
            raise ValueError(f"dirichlet tiene claves desconocidas: {sorted(unknown)}")
        for comp, key in ((0, "ux"), (1, "uy")):
            if key in dirichlet:
                idx = np.asarray(dirichlet[key], dtype=int).reshape(-1)
                if np.any(idx < 0) or np.any(idx >= n):
                    raise ValueError(f"dirichlet[{key!r}] tiene índices fuera de rango")
                dd = u_dof[idx, comp]
                if np.any(dd < 0):
                    badn = int(idx[dd < 0][0])
                    raise ValueError(f"dirichlet[{key!r}]: el nodo {badn} no tiene gdl de sólido")
                fixed.append(dd)
        if "p" in dirichlet:
            idx = np.asarray(dirichlet["p"], dtype=int).reshape(-1)
            if np.any(idx < 0) or np.any(idx >= n):
                raise ValueError("dirichlet['p'] tiene índices fuera de rango")
            dd = p_dof[idx]
            if np.any(dd < 0):
                badn = int(idx[dd < 0][0])
                raise ValueError(f"dirichlet['p']: el nodo {badn} no tiene gdl de presión")
            fixed.append(dd)
    fixed_dofs = np.unique(np.concatenate(fixed)) if fixed else np.zeros(0, dtype=int)
    free = np.setdiff1d(np.arange(size), fixed_dofs, assume_unique=True)
    if free.size == 0:
        raise ValueError("todos los gdl quedaron restringidos por dirichlet")

    a_red = a_full[free][:, free].tocsc()
    # Mismo ordenamiento que `cochlea._assemble_and_solve_arc` (MMD_AT_PLUS_A): el sistema es
    # estructuralmente simétrico y MMD es ~20× más rápido que el COLAMD por defecto, con la misma
    # solución. Tener el mismo en los dos caminos es lo que permite compararlos dígito a dígito.
    lu = spla.splu(a_red, permc_spec="MMD_AT_PLUS_A")
    sol = lu.solve(b_full[free])

    full = np.zeros(size, dtype=complex)
    full[free] = sol
    # residuo sobre las filas libres (en las filas de Dirichlet el residuo es la reacción, no un
    # error de la solución)
    resid_vec = a_full @ full - b_full
    scale = float(np.linalg.norm(b_full[free]))
    residual = float(np.linalg.norm(resid_vec[free]) / (scale if scale > 0.0 else 1.0))

    p = np.full(n, np.nan + 0.0j, dtype=complex)
    u = np.full((n, 2), np.nan + 0.0j, dtype=complex)
    has_p = p_dof >= 0
    p[has_p] = full[p_dof[has_p]]
    has_u = u_dof[:, 0] >= 0
    u[has_u, 0] = full[u_dof[has_u, 0]]
    u[has_u, 1] = full[u_dof[has_u, 1]]

    return CoupledSolution(
        p=p,
        u=u,
        dof_map=dofs,
        size=size,
        frequency=frequency,
        wavenumber=k_f,
        residual=residual,
        interface_edges=iface,
        interface_normals=iface_n,
    )
