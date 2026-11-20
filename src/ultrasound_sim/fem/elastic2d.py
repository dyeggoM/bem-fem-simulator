"""Elastodinámica 2D P1 en plano deformación (Fase 2).

Forma fuerte (isótropo, convención ``e^{−iωt}`` de `docs/decisions.md` D4)::

    ∇·σ + ρ_s·ω²·u = 0,      σ = λ·tr(ε)·I + 2·μ·ε,      ε = ½(∇u + ∇uᵀ)

Forma débil sin términos de borde (el acoplamiento con el fluido entra como tracción
``σ·n = −p·n`` en `fem.coupled2d`) ::

    ∫_Ω (σ(u):ε(v) − ρ_s·ω²·u·v) dA = 0

Elemento P1 **vectorial**: el desplazamiento es lineal por elemento, así que ``ε`` — y por lo
tanto ``σ`` — son constantes por elemento y las integrales volumétricas son exactas (no hay
cuadratura). Grados de libertad locales en orden **nodo mayor**::

    [u_x⁰, u_y⁰, u_x¹, u_y¹, u_x², u_y²]

Matriz de deformación ``B`` (3×6) con el orden de deformaciones ``(ε_xx, ε_yy, γ_xy)``, donde
``γ_xy = 2·ε_xy`` es la deformación de corte **de ingeniería** (la que hace ``σ_xy = μ·γ_xy``)::

        ⎡ ∂φ₀/∂x    0      ∂φ₁/∂x    0      ∂φ₂/∂x    0   ⎤
    B = ⎢    0    ∂φ₀/∂y      0    ∂φ₁/∂y      0    ∂φ₂/∂y ⎥
        ⎣ ∂φ₀/∂y  ∂φ₀/∂x  ∂φ₁/∂y  ∂φ₁/∂x  ∂φ₂/∂y  ∂φ₂/∂x ⎦

con lo cual ``K_e = A·Bᵀ·D·B`` y ``D`` en plano deformación (``ε_zz = 0``)::

        ⎡ λ+2μ   λ     0  ⎤
    D = ⎢  λ    λ+2μ   0  ⎥
        ⎣  0     0     μ  ⎦

``M_e`` es la masa **sin densidad** (``∫ φ_i·φ_j``), porque el término de la ecuación es
``ρ_s·ω²·M_s``: la densidad se aplica al ensamblar, igual que en `fem.assembly` (donde ``K`` y
``M`` también salen sin coeficientes). ``Σ_ij M_ij = 2·A`` para una traslación unitaria: sirve de
test.

Cara de tracción
----------------
`traction_vector` integra ``∫_Γ t·v ds`` con ``t`` constante por arista. Para una arista de
longitud ``L`` y una tracción uniforme, la carga consistente reparte ``t·L/2`` en cada nodo
(regla exacta para ``v`` lineal). ``normal_sign`` multiplica la tracción por ±1 y existe para que
el llamador pueda pasar la tracción *en el sentido de la normal* sin invertirla a mano.

Convención de normal (importante)
---------------------------------
``meshing.outward_normals`` devuelve la normal **saliente del dominio** de la arista tal como la
recorre el triángulo CCW dueño. En el sólido eso es la normal saliente del sólido, que es la que
usa el acoplamiento (``σ·n = −p·n``).
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from ..materials import ElasticSolid
from ..meshing import triangle_areas
from .assembly import triangle_element_areas, triangle_gradients

__all__ = [
    "plane_strain_matrix",
    "strain_matrix",
    "p1_elastic_local",
    "p1_elastic_matrices",
    "traction_vector",
]


# ---------------------------------------------------------------------------
# Validaciones
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
    if np.any((t[:, 0] == t[:, 1]) | (t[:, 1] == t[:, 2]) | (t[:, 0] == t[:, 2])):
        raise ValueError("hay triángulos con nodos repetidos (elemento degenerado)")
    return t


def _check_ccw(x: np.ndarray, t: np.ndarray) -> None:
    signed = triangle_areas(x, t)
    if np.any(signed <= 0.0):
        bad = int(np.flatnonzero(signed <= 0.0)[0])
        raise ValueError(
            "los triángulos deben ser CCW y no degenerados; el triángulo "
            f"{bad} tiene área con signo {signed[bad]!r} (meshing._ensure_ccw los orienta)"
        )


def _check_solid(solid) -> ElasticSolid:
    if not isinstance(solid, ElasticSolid):
        raise TypeError(
            "solid debe ser materials.ElasticSolid (la librería no define otra clase de sólido); "
            f"recibido {type(solid).__name__}"
        )
    if not (solid.mu > 0.0):
        raise ValueError(f"{solid.name}: μ = ρ·c_s² debe ser > 0")
    return solid


# ---------------------------------------------------------------------------
# Matrices locales
# ---------------------------------------------------------------------------


def plane_strain_matrix(solid: ElasticSolid) -> np.ndarray:
    """``D`` (3×3) de plano deformación en el orden ``(ε_xx, ε_yy, γ_xy)``."""
    lam, mu = solid.lam, solid.mu
    return np.array(
        [[lam + 2.0 * mu, lam, 0.0], [lam, lam + 2.0 * mu, 0.0], [0.0, 0.0, mu]],
        dtype=float,
    )


def strain_matrix(grads: np.ndarray) -> np.ndarray:
    """``B`` ``(m, 3, 6)`` a partir de los gradientes ``(m, 3, 2)`` de `assembly.triangle_gradients`.

    Fila 2 usa ``γ_xy = ∂u_x/∂y + ∂u_y/∂x`` (corte de ingeniería).
    """
    g = np.asarray(grads, dtype=float)
    if g.ndim != 3 or g.shape[1] != 3 or g.shape[2] != 2:
        raise ValueError(f"grads debe tener forma (m, 3, 2); recibido {g.shape}")
    m = g.shape[0]
    b = np.zeros((m, 3, 6), dtype=float)
    b[:, 0, 0::2] = g[:, :, 0]  # ε_xx = ∂u_x/∂x
    b[:, 1, 1::2] = g[:, :, 1]  # ε_yy = ∂u_y/∂y
    b[:, 2, 0::2] = g[:, :, 1]  # γ_xy = ∂u_x/∂y
    b[:, 2, 1::2] = g[:, :, 0]  #      + ∂u_y/∂x
    return b


def p1_elastic_local(nodes, triangles, solid) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Matrices locales ``(K_e, M_e, A)`` con formas ``(m,6,6)``, ``(m,6,6)``, ``(m,)``.

    ``K_e = A·Bᵀ·D·B`` (rigidez, exacta porque ``ε`` es constante por elemento) y
    ``M_e = A/12·(1 + I) ⊗ I₂`` (masa consistente **sin densidad**).
    """
    x = _as_nodes(nodes)
    t = _as_triangles(triangles, x.shape[0])
    _check_ccw(x, t)
    sol = _check_solid(solid)

    area = triangle_element_areas(x, t)
    grads = triangle_gradients(x, t)
    b = strain_matrix(grads)
    d = plane_strain_matrix(sol)
    stiff = area[:, None, None] * np.einsum("mki,kl,mlj->mij", b, d, b)

    m = area.size
    base = np.broadcast_to((area / 12.0)[:, None, None], (m, 3, 3)).copy()
    idx = np.arange(3)
    base[:, idx, idx] += (area / 12.0)[:, None]  # A/12·(1 + δ_ij)
    mass = np.zeros((m, 6, 6), dtype=float)
    mass[:, 0::2, 0::2] = base  # bloque u_x
    mass[:, 1::2, 1::2] = base  # bloque u_y
    return stiff, mass, area


def _accumulate(n_dofs: int, dofs: np.ndarray, local: np.ndarray) -> sp.csr_matrix:
    """Suma matrices locales ``(m, n, n)`` en una CSR ``(N, N)`` (dtype de ``local``)."""
    n = dofs.shape[1]
    rows = np.repeat(dofs, n, axis=1).ravel()
    cols = np.tile(dofs, (1, n)).ravel()
    return sp.coo_matrix(
        (np.asarray(local).reshape(-1), (rows, cols)), shape=(n_dofs, n_dofs)
    ).tocsr()


def p1_elastic_matrices(nodes, triangles, solid) -> tuple[sp.csr_matrix, sp.csr_matrix]:
    """Matrices globales ``(K, M)`` reales, ``(2N, 2N)``, para P1 vectorial en plano deformación.

    Orden global **nodo mayor**: el gdl ``2·a + α`` es la componente ``α`` (0: ``x``, 1: ``y``)
    del nodo ``a``. ``M`` no lleva densidad (el término de la ecuación es ``ρ_s·ω²·M``).

    Returns
    -------
    (K, M) : ``scipy.sparse.csr_matrix`` reales de forma ``(2N, 2N)``.
    """
    x = _as_nodes(nodes)
    t = _as_triangles(triangles, x.shape[0])
    stiff, mass, _ = p1_elastic_local(x, t, solid)
    loc = np.stack([2 * t, 2 * t + 1], axis=-1).reshape(t.shape[0], 6)
    n_dofs = 2 * x.shape[0]
    return _accumulate(n_dofs, loc, stiff), _accumulate(n_dofs, loc, mass)


# ---------------------------------------------------------------------------
# Tracción consistente
# ---------------------------------------------------------------------------


def traction_vector(
    nodes, edges, values, normal_sign: int = 1
) -> np.ndarray:
    """Vector global ``∫_Γ t·v ds`` con ``t`` constante por arista, forma ``(2N,)``.

    Para la arista ``(a, b)`` de longitud ``L`` y tracción ``t``::

        ∫_e t·φ_a ds = t·L/2 ,      ∫_e t·φ_b ds = t·L/2

    (exacto, porque ``φ`` es lineal sobre la arista y ``t`` constante).

    Parameters
    ----------
    nodes : ``(N, 2)`` coordenadas.
    edges : ``(B, 2)`` índices de nodo de cada arista cargada.
    values : ``(B, 2)`` (o ``(2,)``, que se repite) tracción ``[t_x, t_y]`` en Pa por arista.
    normal_sign : ``+1`` deja ``values``; ``−1`` lo invierte. Sirve para pasar una presión
        ``p`` como tracción ``−p·n`` sin cambiar el signo a mano. Sólo se acepta ``±1``.

    Returns
    -------
    ``(2N,)`` con el dtype de ``values`` (real o complejo según lo que se pase).
    """
    x = _as_nodes(nodes)
    e = np.asarray(edges, dtype=int)
    if e.ndim != 2 or e.shape[1] != 2:
        raise ValueError(f"edges debe tener forma (B, 2); recibido {e.shape}")
    if e.shape[0] == 0:
        raise ValueError("se requiere al menos una arista cargada")
    if np.any(e < 0) or np.any(e >= x.shape[0]):
        raise ValueError("edges tiene índices de nodo fuera de rango")
    if np.any(e[:, 0] == e[:, 1]):
        raise ValueError("hay aristas con sus dos nodos iguales")
    if int(normal_sign) not in (-1, 1):
        raise ValueError(f"normal_sign debe ser +1 o -1; recibido {normal_sign!r}")

    vals = np.asarray(values)
    if vals.ndim == 1:
        if vals.size != 2:
            raise ValueError(f"values (2,) debe ser [t_x, t_y]; recibido tamaño {vals.size}")
        vals = np.broadcast_to(vals, (e.shape[0], 2)).copy()
    elif vals.ndim == 2:
        if vals.shape != (e.shape[0], 2):
            raise ValueError(
                f"values debe tener forma (B, 2) = ({e.shape[0]}, 2); recibido {vals.shape}"
            )
    else:
        raise ValueError(f"values debe ser (2,) o (B, 2); recibido forma {vals.shape}")
    if not np.all(np.isfinite(vals)):
        raise ValueError("values contiene valores no finitos")

    lengths = np.linalg.norm(x[e[:, 1]] - x[e[:, 0]], axis=1)
    if np.any(lengths <= 0.0):
        raise ValueError("hay aristas de longitud nula")

    vals = int(normal_sign) * vals
    contrib = 0.5 * lengths[:, None] * vals  # (B, 2) = t·L/2
    out = np.zeros(2 * x.shape[0], dtype=vals.dtype)
    np.add.at(out, 2 * e[:, 0], contrib[:, 0])
    np.add.at(out, 2 * e[:, 0] + 1, contrib[:, 1])
    np.add.at(out, 2 * e[:, 1], contrib[:, 0])
    np.add.at(out, 2 * e[:, 1] + 1, contrib[:, 1])
    return out
