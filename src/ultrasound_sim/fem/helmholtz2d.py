"""FEM P1 de Helmholtz 2D con ABC de primer orden y campo total (Fase 1).

Forma débil implementada (convención D4: ``e^{−iωt}``, ``∇²p + k²p = 0``, ``k = ω/c + i·α``)::

    ∫ (1/ρ)·(∇p·∇v − k² p v) dΩ − ∮_Γ (i·k/ρ)·p·v ds
        = ∮_Γ (1/ρ)·(∂p_i/∂n − i·k·p_i)·v ds

donde ``p`` es el campo **total** y ``p_i`` la onda plana incidente. Con densidad uniforme
(el valor por defecto, ``ρ = 1``) esto es literalmente la forma ``∇p·∇v − k²·p·v`` del plan; el
peso ``1/ρ`` sólo entra cuando se pasa ``element_rho``, y es necesario porque la ecuación fuerte
en un medio heterogéneo es ``∇·((1/ρ)·∇p) + (k²/ρ)·p = 0``: con ``ρ`` variable, ``∇²p + k²p = 0``
no vale aunque ``k = ω/c`` sea correcto por elemento. La condición de borde
``∂p_s/∂n = i·k·p_s`` sobre el campo dispersado ``p_s = p − p_i`` es la ABC de primer orden
(impedancia); en incidencia normal es exacta, en oblicua refleja una fracción
``(cosθ − 1)/(cosθ + 1)``.

El número de onda es **complejo y por elemento** (``element_k``), de modo que un medio atenuante
se modela simplemente con su ``Material.wavenumber(f)``. En las aristas de la ABC se usa el ``k``
del elemento dueño de la arista (cada arista de frontera pertenece a un único triángulo).

Caras iluminadas vs caras en sombra
-----------------------------------
El término ``∂p_i/∂n − i·k·p_i`` del lado derecho sólo se aplica en las aristas donde la onda
incidente **entra** al dominio (``n·d < 0``, ver ``illuminated``). En las caras de salida el campo
total es puramente saliente y no hay término incidente. Si se aplicara también ahí, en un medio
estratificado aparecería una fuente espuria —el ``k`` local (acero) no es el ``k`` de la onda
incidente (agua)— y el resultado mostraría una onda reflejada falsa con R ≈ 1. Es el error más
fácil de cometer con esta formulación.

Campo dispersado (``assemble_scattered``)
-----------------------------------------
La formulación de campo total de arriba sólo es consistente si la frontera exterior es alineada con
la onda (caja/ducto): en una cara de sombra **curva** el término ``∂p_i/∂n − i·k·p_i`` no se anula
porque ``n·d ≠ 1``, y el error no baja con ``ppw`` (es un error de modelo, no de discretización).

La alternativa, que es la que hay que usar en un contorno curvo (cilindro, esfera, malla polar), es
resolver el campo **dispersado** ``p_s = p − p_i``, que es puramente saliente en *todo* el contorno::

    ∫ (1/ρ)·(∇p_s·∇v − k² p_s v) dΩ − ∮_ABC (i·k/ρ)·p_s·v ds
        = − ∮_Γint (1/ρ)·(∂p_s/∂n)·v ds

El contorno ABC queda homogéneo (no hay máscara de iluminación: ``p_s`` sale por los 360°) y la onda
incidente entra **sólo** como dato de Neumann sobre la superficie interior rígida,
``∂p_total/∂n = 0 ⇒ ∂p_s/∂n = −∂p_i/∂n``, o sea ``b = −∮_Γint (1/ρ)·(∂p_i/∂n)·v ds``. El signo sale
de la integración por partes: el término de borde del volumen es ``−∮(1/ρ)∂p_s/∂n·v``; al pasar
``∂p_s/∂n = −∂p_i/∂n`` queda ``+∮(1/ρ)∂p_i/∂n·v`` del lado izquierdo, que es ``b`` con signo menos.

Ojo con la normal: ``∂p_i/∂n`` se evalúa con la normal **saliente del dominio FEM**, que en el
agujero del anillo apunta hacia el eje del cilindro (``n = −r̂``). Con la normal contraria el campo
dispersado sale con el signo invertido y el total da ~77 % de error.

``∂p_s/∂n = i·k·p_s`` es exacta sólo para una onda plana. Para una onda cilíndrica saliente sobre un
contorno circular de radio ``r`` vale ``∂p_s/∂n = (i·k − 1/(2r) + i/(8k r²))·p_s + (i/(2k))·∂²_s p_s``
(``abc_order=2``, Bayliss–Turkel–Gunzburger): el término ``1/(2r)`` es la curvatura del frente y el
``∂²_s`` la parte angular. Sin ellos la ABC de primer orden refleja ~6 % a ``kR = 3``, y ese error
tampoco baja con ``ppw``.

Notas de discretización
-----------------------
- La malla debe tener nodos sobre las interfaces de material; si un elemento cruza una interfaz,
  su ``k`` (constante por elemento) promedia dos medios y arruina el contraste. Con
  `meshing.rect_tri_mesh` alcanza con que la interfaz caiga sobre una línea de nodos.
- El término de borde se integra con la masa P1 de la arista (``k_e·S_e``), que para ``k``
  constante por arista es la discretización de Galerkin consistente.
- El lado derecho se integra con la regla del punto medio (exacta para ``v`` lineal sobre la
  arista cuando el integrando se evalúa en el punto medio).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from ..materials import Material
from ..meshing import boundary_edges as _boundary_edges
from ..meshing import edge_midpoints, outward_normals
from .assembly import p1_line_local, p1_tri_local, weights_to_local

__all__ = [
    "PlaneWave",
    "plane_wave",
    "assemble",
    "solve",
    "assemble_scattered",
    "solve_scattered",
    "total_field",
    "owner_elements",
]


# ---------------------------------------------------------------------------
# Onda plana incidente
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class PlaneWave:
    """Onda plana ``p_i(x) = A·exp(i·k·d·x)`` con ``d`` unitaria.

    Atributos
    ---------
    k : número de onda complejo del medio [1/m].
    direction : dirección de propagación unitaria.
    amplitude : amplitud compleja [Pa].
    """

    k: complex
    direction: np.ndarray
    amplitude: complex = 1.0

    @property
    def wavenumber_vector(self) -> np.ndarray:
        return self.k * self.direction

    def value(self, points) -> np.ndarray:
        """Presión incidente en ``points`` (``(M, 2)`` o ``(2,)``)."""
        pts = np.atleast_2d(np.asarray(points, dtype=float))
        phase = pts @ self.wavenumber_vector
        return self.amplitude * np.exp(1j * phase)

    # alias explícito para lectura del código de ensamblado
    evaluate = value

    def gradient(self, points) -> np.ndarray:
        """``∇p_i = i·k·d·p_i`` en ``points``, forma ``(M, 2)``."""
        pts = np.atleast_2d(np.asarray(points, dtype=float))
        return 1j * self.value(pts)[:, None] * self.wavenumber_vector[None, :]

    def normal_derivative(self, points, normals) -> np.ndarray:
        """``∂p_i/∂n = i·k·(d·n)·p_i`` con normales unitarias, forma ``(M,)``."""
        pts = np.atleast_2d(np.asarray(points, dtype=float))
        nrm = np.atleast_2d(np.asarray(normals, dtype=float))
        return 1j * self.value(pts) * (nrm @ self.wavenumber_vector)


def plane_wave(
    frequency: float,
    material: Material,
    direction,
    amplitude: complex = 1.0,
) -> PlaneWave:
    """Construye la onda plana incidente para un material y una dirección."""
    d = np.asarray(direction, dtype=float).reshape(-1)
    if d.size != 2:
        raise ValueError("direction debe ser un vector 2D")
    norm = np.linalg.norm(d)
    if norm == 0.0:
        raise ValueError("direction no puede ser nula")
    return PlaneWave(
        k=material.wavenumber(frequency),
        direction=d / norm,
        amplitude=complex(amplitude),
    )


# ---------------------------------------------------------------------------
# Utilidades de malla
# ---------------------------------------------------------------------------


def owner_elements(triangles) -> dict[tuple[int, int], int]:
    """Mapa ``arista ordenada → índice del triángulo`` (última ocurrencia gana)."""
    t = np.asarray(triangles, dtype=int)
    owners: dict[tuple[int, int], int] = {}
    for tri_i in range(t.shape[0]):
        a, b, c = t[tri_i]
        for i, j in ((a, b), (b, c), (c, a)):
            owners[(min(i, j), max(i, j))] = tri_i
    return owners


def _per_element(value, n_elements: int, name: str) -> np.ndarray:
    """Normaliza un parámetro escalar o por elemento a un array ``(M,)``."""
    if value is None:
        return np.ones(n_elements, dtype=float)
    arr = np.asarray(value)
    if arr.ndim == 0:
        return np.full(n_elements, arr.item(), dtype=arr.dtype)
    if arr.size != n_elements:
        raise ValueError(f"{name} debe ser escalar o tener un valor por elemento")
    return arr


def _edge_values(triangles, element_values: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """Valor (``k``, ``ρ``) del elemento dueño de cada arista dada."""
    owners = owner_elements(triangles)
    out = np.empty(edges.shape[0], dtype=element_values.dtype)
    for e, (i, j) in enumerate(edges):
        tri = owners.get((min(int(i), int(j)), max(int(i), int(j))))
        if tri is None:
            raise ValueError(f"arista {e} no pertenece a ningún triángulo")
        out[e] = element_values[tri]
    return out


def _edge_wavenumbers(triangles, element_k, edges) -> np.ndarray:
    """``k`` del elemento dueño de cada arista dada (compatibilidad)."""
    return _edge_values(triangles, _per_element(element_k, np.asarray(triangles).shape[0], "element_k"), edges)


def _illuminated_mask(normals: np.ndarray, p_inc) -> np.ndarray:
    """``True`` donde la onda incidente entra al dominio (``n·d < 0``).

    Si ``p_inc`` no expone ``direction`` (no es una onda plana) se considera iluminada toda la
    frontera de la ABC: es la elección correcta cuando el término incidente es sólo una fuente
    local, y hay que pasarla explícitamente en cualquier otro caso.
    """
    direction = getattr(p_inc, "direction", None)
    if direction is None:
        return np.ones(normals.shape[0], dtype=bool)
    return (normals @ np.asarray(direction, dtype=float)) < 0.0


def _as_boundary_edges(nodes, triangles, edges) -> np.ndarray:
    if edges is None:
        return _boundary_edges(np.asarray(triangles, dtype=int))
    return np.asarray(edges, dtype=int).reshape(-1, 2)


# ---------------------------------------------------------------------------
# Ensamblado y resolución
# ---------------------------------------------------------------------------


def assemble(
    nodes,
    triangles,
    element_k,
    abc_edges=None,
    p_inc: PlaneWave | None = None,
    edge_k=None,
    normals=None,
    element_rho=None,
    edge_rho=None,
    illuminated=None,
) -> tuple[sp.csr_matrix, np.ndarray]:
    """Ensambla el sistema ``A·p = b`` del campo total con ABC de primer orden.

    Parameters
    ----------
    nodes : ``(N, 2)`` coordenadas.
    triangles : ``(M, 3)`` triángulos CCW.
    element_k : ``(M,)`` (o escalar) número de onda complejo por elemento.
    abc_edges : ``(B, 2)`` aristas donde se aplica la ABC; ``None`` ⇒ toda la frontera.
    p_inc : :class:`PlaneWave` incidente (si es ``None`` no hay término fuente: campo saliente).
    edge_k : ``(B,)`` complejo por arista; si es ``None`` se toma el del elemento dueño.
    normals : ``(B, 2)`` normales salientes; si es ``None`` se calculan con `meshing.outward_normals`.
    element_rho : ``(M,)`` o escalar densidad por elemento [kg/m³]; **por defecto 1.0**, que deja
        la ecuación exactamente en la forma pedida ``∇p·∇v − k²·p·v``. Con densidad variable hay
        que pasar ``ρ`` porque la ecuación fuerte correcta es
        ``∇·((1/ρ)·∇p) + (k²/ρ)·p = 0``: el peso ``1/ρ`` es imprescindible cuando el contraste de
        impedancia viene de ``ρ`` y no sólo de ``c`` (agua/acero, p. ej.).
    edge_rho : ``(B,)`` densidad en cada arista de la ABC; ``None`` ⇒ la del elemento dueño.
    illuminated : ``(B,)`` bool; ``True`` en las aristas que la onda incidente **entra**
        (caras de incidencia) y ``False`` en las de sombra. ``None`` ⇒ se deduce de
        ``p_inc.direction`` con el criterio ``n·d < 0``. Es imprescindible que las caras de
        sombra queden en ``False``: en ellas el campo total es saliente y el término
        ``∂p_i/∂n − i·k·p_i`` no se anula si ``k`` local ≠ ``k`` del medio incidente (caso
        estratificado), y actúa como una fuente espuria que mete una onda reflejada falsa.

    Returns
    -------
    (A, b) : matriz dispersa compleja ``(N, N)`` y vector ``(N,)``.
    """
    x = np.asarray(nodes, dtype=float)
    t = np.asarray(triangles, dtype=int)
    n = x.shape[0]

    k_el = np.asarray(element_k, dtype=complex).reshape(-1)
    if k_el.size == 1:
        k_el = np.full(t.shape[0], k_el[0], dtype=complex)
    if k_el.size != t.shape[0]:
        raise ValueError("element_k debe ser escalar o tener un valor por elemento")

    rho_el = _per_element(element_rho, t.shape[0], "element_rho")
    if np.any(rho_el <= 0.0):
        raise ValueError("element_rho debe ser > 0")

    mass, stiff, _ = p1_tri_local(x, t)
    # peso 1/ρ (≡ 1 cuando no se pasa densidad, que es la forma pedida en el plan)
    local = (stiff - (k_el**2)[:, None, None] * mass) / rho_el[:, None, None]

    dofs = np.repeat(t, 3, axis=1).ravel()
    cols = np.tile(t, (1, 3)).ravel()
    a = sp.coo_matrix(
        (local.reshape(-1), (dofs, cols)), shape=(n, n), dtype=complex
    ).tocsr()
    b = np.zeros(n, dtype=complex)

    edges = _as_boundary_edges(x, t, abc_edges)
    if edges.shape[0] == 0:
        return a.tocsr(), b

    if edge_k is None:
        k_edge = _edge_wavenumbers(t, k_el, edges)
    else:
        k_edge = np.asarray(edge_k, dtype=complex).reshape(-1)
        if k_edge.size != edges.shape[0]:
            raise ValueError("edge_k debe tener un valor por arista de abc_edges")

    if normals is None:
        nrm = outward_normals(x, edges)
    else:
        nrm = np.asarray(normals, dtype=float).reshape(-1, 2)

    if edge_rho is None:
        rho_edge = _edge_values(t, rho_el, edges)
    else:
        rho_edge = np.asarray(edge_rho, dtype=float).reshape(-1)
        if rho_edge.size != edges.shape[0]:
            raise ValueError("edge_rho debe tener un valor por arista de abc_edges")
    w_edge = 1.0 / rho_edge

    lengths = np.linalg.norm(x[edges[:, 1]] - x[edges[:, 0]], axis=1)

    # Término de borde: −∫ (i·k/ρ)·p·v ds
    local_abc = (-1j * k_edge * w_edge)[:, None, None] * weights_to_local(lengths)
    rows = np.repeat(edges, 2, axis=1).ravel()
    cols_e = np.tile(edges, (1, 2)).ravel()
    a = a + sp.coo_matrix(
        (local_abc.reshape(-1), (rows, cols_e)), shape=(n, n), dtype=complex
    ).tocsr()

    # Lado derecho: ∫ (∂p_i/∂n − i·k·p_i)·v ds con regla del punto medio, sólo en las caras
    # iluminadas (en las de sombra el campo total ya es saliente y no hay término incidente).
    if p_inc is not None:
        if illuminated is None:
            illuminated = _illuminated_mask(nrm, p_inc)
        else:
            illuminated = np.asarray(illuminated, dtype=bool).reshape(-1)
            if illuminated.size != edges.shape[0]:
                raise ValueError("illuminated debe tener un valor por arista de abc_edges")
        mid = edge_midpoints(x, edges)
        f = p_inc.normal_derivative(mid, nrm) - 1j * k_edge * p_inc.value(mid)
        contrib = np.where(illuminated, w_edge * f * lengths / 2.0, 0.0)
        np.add.at(b, edges[:, 0], contrib)
        np.add.at(b, edges[:, 1], contrib)

    return a.tocsr(), b


def solve(
    nodes,
    triangles,
    element_k,
    abc_edges=None,
    p_inc: PlaneWave | None = None,
    edge_k=None,
    normals=None,
    element_rho=None,
    edge_rho=None,
    illuminated=None,
) -> np.ndarray:
    """Resuelve el Helmholtz 2D P1 y devuelve la presión total compleja en los nodos.

    Ver :func:`assemble` para el significado de los argumentos. ``element_rho`` es opcional:
    sin él la ecuación es exactamente ``∫(∇p·∇v − k²pv)`` (densidad uniforme), con él se agrega
    el peso ``1/ρ`` que hace falta cuando el contraste de impedancia es de densidad.
    """
    a, b = assemble(
        nodes, triangles, element_k, abc_edges=abc_edges, p_inc=p_inc,
        edge_k=edge_k, normals=normals, element_rho=element_rho, edge_rho=edge_rho,
        illuminated=illuminated,
    )
    return spla.spsolve(a.tocsc(), b)


# ---------------------------------------------------------------------------
# Campo dispersado p_s = p − p_i
# ---------------------------------------------------------------------------


def _boundary_edge_lookup(triangles) -> dict[tuple[int, int], tuple[int, int]]:
    """``arista ordenada → arista orientada`` tal como la recorre el triángulo CCW dueño."""
    lookup: dict[tuple[int, int], tuple[int, int]] = {}
    for i, j in _boundary_edges(np.asarray(triangles, dtype=int)):
        lookup[(min(int(i), int(j)), max(int(i), int(j)))] = (int(i), int(j))
    return lookup


def _oriented_boundary_edges(edges, lookup, name: str) -> np.ndarray:
    """Reorienta ``edges`` como las tiene la malla y valida que sean de frontera.

    La normal saliente de `meshing.outward_normals` depende del sentido en que se recorre la
    arista; si el usuario pasa la arista al revés, ``n`` cambia de signo y el dato de Neumann
    entra con el signo equivocado. Reorientar contra la malla elimina esa ambigüedad.
    """
    e = np.asarray(edges, dtype=int).reshape(-1, 2)
    out = np.empty_like(e)
    for row, (i, j) in enumerate(e):
        if int(i) == int(j):
            raise ValueError(f"{name}: la arista {row} tiene sus dos nodos iguales")
        key = (min(int(i), int(j)), max(int(i), int(j)))
        oriented = lookup.get(key)
        if oriented is None:
            raise ValueError(
                f"{name}: la arista {row} = {key} no es una arista de frontera de la malla; "
                "la normal saliente sólo está definida sobre la frontera"
            )
        out[row] = oriented
    return out


def _check_closed_chain(edges: np.ndarray, name: str) -> None:
    """Exige que las aristas formen una cadena cerrada (cada nodo, exactamente dos veces).

    El término tangencial de ``abc_order=2`` usa ``∫ ∂²_s p·v ds = −∫ ∂_s p·∂_s v ds``, igualdad
    que vale en una cadena cerrada (los términos de borde de la integración por partes se cancelan).
    """
    nodes, counts = np.unique(edges.reshape(-1), return_counts=True)
    if np.any(counts != 2):
        bad = nodes[counts != 2][:3]
        raise ValueError(
            f"{name}: abc_order=2 necesita que las aristas de la ABC formen una cadena cerrada "
            f"(cada nodo de la frontera aparece en dos aristas); los nodos {bad.tolist()} "
            "aparecen un número distinto de veces"
        )


def _chain_stiffness(n_dofs: int, edges: np.ndarray, lengths: np.ndarray,
                     scale=None) -> sp.csr_matrix:
    """Rigidez 1D ``∫ ∂_s u·∂_s v ds`` sobre la cadena de aristas de la ABC.

    ``scale`` es un factor complejo por arista (``i/(2·k·ρ)`` en la ABC de orden 2).
    """
    _, stiff, _ = p1_line_local(np.array([0.0, 1.0]), np.array([[0, 1]]))
    local = stiff[0][None, :, :] / lengths[:, None, None]
    if scale is not None:
        local = local * np.asarray(scale, dtype=complex).reshape(-1, 1, 1)
    rows = np.repeat(edges, 2, axis=1).ravel()
    cols = np.tile(edges, (1, 2)).ravel()
    return sp.coo_matrix(
        (local.reshape(-1), (rows, cols)), shape=(n_dofs, n_dofs), dtype=complex
    ).tocsr()


def assemble_scattered(
    nodes,
    triangles,
    element_k,
    abc_edges=None,
    neumann_edges=None,
    dp_dn_incident=None,
    element_rho=None,
    edge_rho=None,
    normals=None,
    abc_order: int = 1,
    abc_center=(0.0, 0.0),
    neumann_normals=None,
) -> tuple[sp.csr_matrix, np.ndarray]:
    r"""Ensambla ``A·p_s = b`` para el campo **dispersado** ``p_s = p − p_i``.

    Forma débil::

        ∫ (1/ρ)·(∇p_s·∇v − k²·p_s·v) dΩ − ∮_ABC (i·k/ρ)·p_s·v ds
            = − ∮_Γint (1/ρ)·(∂p_i/∂n)·v ds

    El contorno ABC es **homogéneo**: ``p_s`` es saliente en todo él, así que no hay máscara de
    iluminación/sombra (esa máscara es justamente lo que falla en un contorno curvo, ver el
    docstring del módulo). La onda incidente entra sólo por el dato de Neumann de la superficie
    interior, ``∂p_s/∂n = −∂p_i/∂n`` (pared rígida: ``∂p_total/∂n = 0``).

    Parameters
    ----------
    nodes, triangles, element_k : como en :func:`assemble`.
    abc_edges : ``(B, 2)`` aristas del contorno absorbente; ``None`` ⇒ no se agrega ABC.
    neumann_edges : ``(E, 2)`` aristas donde se impone el dato ``∂p_s/∂n = −∂p_i/∂n``
        (la superficie del dispersor). Se reorientan contra la malla, así que el orden de los dos
        nodos de cada fila es indistinto.
    dp_dn_incident : ``(E,)`` ``∂p_i/∂n`` en el punto medio de cada arista de ``neumann_edges``,
        con la normal **saliente del dominio FEM** (en el agujero del anillo, ``n = −r̂``). Se
        calcula con ``p_inc.normal_derivative(midpoints, outward_normals(nodes, neumann_edges))``.
    element_rho, edge_rho : como en :func:`assemble` (peso ``1/ρ``).
    normals : ``(B, 2)`` normales de la ABC; la ABC no las necesita para integrar (el signo lo
        fija el dominio), pero si se pasan se validan contra la malla para detectar un giro
        invertido. ``None`` ⇒ no se valida nada.
    abc_order : ``1`` (ABC de impedancia ``∂p_s/∂n = i·k·p_s``, la misma de :func:`assemble`) o
        ``2`` (Bayliss–Turkel–Gunzburger sobre un contorno circular:
        ``∂p_s/∂n = (i·k − 1/(2r) + i/(8k r²))·p_s + (i/(2k))·∂²_s p_s``). El orden 2 sólo es
        válido si la ABC es una circunferencia de radio ``r`` centrada en ``abc_center`` y sus
        aristas forman una cadena cerrada; en ese caso es notablemente más preciso que el orden 1
        (a ``kR = 3``: ~1 % contra ~6 %, y el orden 1 no mejora al refinar la malla).
    abc_center : centro de la circunferencia de la ABC (sólo ``abc_order=2``).
    neumann_normals : ``(E, 2)`` normales de ``neumann_edges``; no hacen falta para ensamblar (el
        dato ya viene en ``dp_dn_incident``) y sirven para auditar que no estén invertidas.

    Returns
    -------
    (A, b) : matriz dispersa compleja ``(N, N)`` y vector ``(N,)``, con ``p_s`` en los nodos.
    """
    x = np.asarray(nodes, dtype=float)
    t = np.asarray(triangles, dtype=int)
    n = x.shape[0]

    k_el = np.asarray(element_k, dtype=complex).reshape(-1)
    if k_el.size == 1:
        k_el = np.full(t.shape[0], k_el[0], dtype=complex)
    if k_el.size != t.shape[0]:
        raise ValueError(
            f"element_k debe ser escalar o tener un valor por elemento "
            f"({t.shape[0]} elementos, recibidos {k_el.size})"
        )

    rho_el = _per_element(element_rho, t.shape[0], "element_rho")
    if np.any(rho_el <= 0.0):
        raise ValueError("element_rho debe ser > 0")

    # Parte de volumen, idéntica a assemble() (ver allí la nota del peso 1/ρ).
    mass, stiff, _ = p1_tri_local(x, t)
    local = (stiff - (k_el**2)[:, None, None] * mass) / rho_el[:, None, None]
    dofs = np.repeat(t, 3, axis=1).ravel()
    cols = np.tile(t, (1, 3)).ravel()
    a = sp.coo_matrix(
        (local.reshape(-1), (dofs, cols)), shape=(n, n), dtype=complex
    ).tocsr()
    b = np.zeros(n, dtype=complex)

    if abc_order not in (1, 2):
        raise ValueError(f"abc_order debe ser 1 o 2; recibido {abc_order!r}")

    lookup = _boundary_edge_lookup(t)

    # ------------------------------------------------------------------ ABC
    if abc_edges is not None:
        abc = _oriented_boundary_edges(abc_edges, lookup, "abc_edges")
        if abc.shape[0]:
            if edge_rho is None:
                rho_abc = _edge_values(t, rho_el, abc)
            else:
                rho_abc = np.asarray(edge_rho, dtype=float).reshape(-1)
                if rho_abc.size != abc.shape[0]:
                    raise ValueError(
                        f"edge_rho debe tener un valor por arista de abc_edges "
                        f"({abc.shape[0]}, recibidos {rho_abc.size})"
                    )
            w_abc = 1.0 / rho_abc
            lengths = np.linalg.norm(x[abc[:, 1]] - x[abc[:, 0]], axis=1)
            if normals is not None:
                nrm = np.asarray(normals, dtype=float).reshape(-1, 2)
                if nrm.shape[0] != abc.shape[0]:
                    raise ValueError(
                        f"normals debe tener una normal por arista de abc_edges "
                        f"({abc.shape[0]}, recibidas {nrm.shape[0]})"
                    )
                # la ABC no usa la normal para nada (∂p_s/∂n = i·k·p_s vale sobre el contorno
                # saliente); se valida el signo sólo para que un error del usuario no pase
                # inadvertido, comparándola con la normal de la malla.
                check = np.einsum("ij,ij->i", nrm, outward_normals(x, abc))
                if np.any(check <= 0.0):
                    raise ValueError(
                        "normals tiene aristas con la normal invertida respecto de la malla "
                        "(n·n_malla ≤ 0); las normales de la ABC son salientes del dominio"
                    )

            if neumann_edges is not None:
                key_abc = {tuple(sorted(map(int, e))) for e in abc}
                key_neu = {tuple(sorted(map(int, e))) for e in
                           np.asarray(neumann_edges, dtype=int).reshape(-1, 2)}
                overlap = key_abc & key_neu
                if overlap:
                    raise ValueError(
                        "abc_edges y neumann_edges se solapan; una arista no puede ser ABC y "
                        f"dato de Neumann a la vez (p. ej. {sorted(overlap)[0]})"
                    )

            k_abc = np.asarray(_edge_wavenumbers(t, k_el, abc), dtype=complex)
            if abc_order == 1:
                coef = 1j * k_abc
            else:
                _check_closed_chain(abc, "abc_edges")
                mid = edge_midpoints(x, abc)
                center = np.asarray(abc_center, dtype=float).reshape(1, 2)
                r_e = np.linalg.norm(mid - center, axis=1)
                if np.any(r_e <= 0.0):
                    raise ValueError(
                        "abc_order=2 necesita un radio positivo para cada arista de la ABC; "
                        "¿abc_center cae sobre el contorno?"
                    )
                coef = 1j * k_abc - 0.5 / r_e + 0.125j / (k_abc * r_e**2)
                # término tangencial: −∫(1/ρ)(i/2k)∂²_s p_s·v = +∫(1/ρ)(i/2k)∂_s p_s·∂_s v
                # (cadena cerrada: el término de borde de la integración por partes se cancela)
                a = a + _chain_stiffness(n, abc, lengths, scale=0.5j * w_abc / k_abc)

            # término −∫ (i·k/ρ)·p_s·v ds
            local_abc = (-coef * w_abc)[:, None, None] * weights_to_local(lengths)
            rows = np.repeat(abc, 2, axis=1).ravel()
            cols_e = np.tile(abc, (1, 2)).ravel()
            a = a + sp.coo_matrix(
                (local_abc.reshape(-1), (rows, cols_e)), shape=(n, n), dtype=complex
            ).tocsr()

    # ------------------------------------------------- dato de Neumann (p_i)
    if neumann_edges is None:
        if dp_dn_incident is not None:
            raise ValueError("dp_dn_incident se pasó sin neumann_edges")
    else:
        neu = _oriented_boundary_edges(neumann_edges, lookup, "neumann_edges")
        if dp_dn_incident is None:
            raise ValueError(
                "neumann_edges requiere dp_dn_incident = ∂p_i/∂n sobre sus puntos medios "
                "(normal saliente del dominio)"
            )
        dp = np.asarray(dp_dn_incident, dtype=complex).reshape(-1)
        if dp.size != neu.shape[0]:
            raise ValueError(
                f"dp_dn_incident debe tener un valor por arista de neumann_edges "
                f"({neu.shape[0]}, recibidos {dp.size})"
            )
        if neumann_normals is not None:
            nrm_n = np.asarray(neumann_normals, dtype=float).reshape(-1, 2)
            if nrm_n.shape[0] != neu.shape[0]:
                raise ValueError(
                    f"neumann_normals debe tener una normal por arista de neumann_edges "
                    f"({neu.shape[0]}, recibidas {nrm_n.shape[0]})"
                )
            if np.any(np.einsum("ij,ij->i", nrm_n, outward_normals(x, neu)) <= 0.0):
                raise ValueError(
                    "neumann_normals tiene aristas con la normal invertida respecto de la malla; "
                    "∂p_i/∂n se define con la normal saliente del dominio FEM (en el agujero "
                    "del anillo, n = −r̂)"
                )
        rho_n = _edge_values(t, rho_el, neu)
        lengths_n = np.linalg.norm(x[neu[:, 1]] - x[neu[:, 0]], axis=1)
        # b = −∮ (1/ρ)·(∂p_i/∂n)·v ds, regla del punto medio (igual que assemble()).
        contrib = -(dp * lengths_n / 2.0) / rho_n
        np.add.at(b, neu[:, 0], contrib)
        np.add.at(b, neu[:, 1], contrib)

    return a.tocsr(), b


def solve_scattered(
    nodes,
    triangles,
    element_k,
    abc_edges=None,
    neumann_edges=None,
    dp_dn_incident=None,
    element_rho=None,
    edge_rho=None,
    normals=None,
    abc_order: int = 1,
    abc_center=(0.0, 0.0),
    neumann_normals=None,
) -> np.ndarray:
    """Resuelve el campo **dispersado** ``p_s`` en los nodos (ver :func:`assemble_scattered`).

    Devuelve sólo ``p_s``; para el campo físico hay que sumar la incidente con
    :func:`total_field`. Se eligió devolver ``p_s`` (y no un dict) para que la firma sea paralela a
    :func:`solve`.
    """
    a, b = assemble_scattered(
        nodes, triangles, element_k, abc_edges=abc_edges, neumann_edges=neumann_edges,
        dp_dn_incident=dp_dn_incident, element_rho=element_rho, edge_rho=edge_rho,
        normals=normals, abc_order=abc_order, abc_center=abc_center,
        neumann_normals=neumann_normals,
    )
    return spla.spsolve(a.tocsc(), b)


def total_field(p_scattered, p_inc, nodes=None) -> np.ndarray:
    """``p_total = p_s + p_i`` en los nodos.

    Parameters
    ----------
    p_scattered : ``(N,)`` salida de :func:`solve_scattered`.
    p_inc : :class:`PlaneWave` o callable ``f(points) -> (M,)``; la onda incidente.
    nodes : ``(N, 2)`` coordenadas donde evaluar ``p_i``; si se omite se usa ``p_inc`` sobre sus
        propios puntos (sólo válido si ``len(p_inc(nodes)) == len(p_scattered)``).
    """
    ps = np.asarray(p_scattered, dtype=complex).reshape(-1)
    if hasattr(p_inc, "value"):
        if nodes is None:
            raise ValueError("nodes es obligatorio cuando p_inc es una PlaneWave")
        pi = np.asarray(p_inc.value(nodes), dtype=complex).reshape(-1)
    elif callable(p_inc):
        if nodes is None:
            raise ValueError("nodes es obligatorio cuando p_inc es un callable")
        pi = np.asarray(p_inc(nodes), dtype=complex).reshape(-1)
    else:
        raise TypeError("p_inc debe ser un PlaneWave (con .value) o un callable")
    if pi.size != ps.size:
        raise ValueError(
            f"p_inc devolvió {pi.size} valores y p_scattered tiene {ps.size}; "
            "los nodos no coinciden"
        )
    return ps + pi
