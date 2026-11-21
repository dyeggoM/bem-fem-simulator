"""Elasticidad lineal 3D con tetraedros P1 (Fase 6b).

Formulación
-----------
Desplazamiento ``u = (u_x, u_y, u_z)`` por nodo, tres grados de libertad por nodo en orden
``[u_x⁰, u_y⁰, u_z⁰, u_x¹, ...]``. Elemento P1: deformación **constante** por tetraedro, así que la
matriz ``B`` (6×12) y ``K_e = V·Bᵀ·D·B`` son exactas para campos lineales.

Voigt con **deformación de corte de ingeniería** ``γ = 2ε``:

    ε = (ε_xx, ε_yy, ε_zz, γ_yz, γ_xz, γ_xy)

que es lo que hace que ``D`` tenga ``μ`` (y no ``2μ``) en las tres últimas diagonales.

    D = [[λ+2μ, λ, λ, 0, 0, 0],
         [λ, λ+2μ, λ, 0, 0, 0],
         [λ, λ, λ+2μ, 0, 0, 0],
         [0, 0, 0, μ, 0, 0],
         [0, 0, 0, 0, μ, 0],
         [0, 0, 0, 0, 0, μ]]

Masa: consistente, ``M_e = ρV/20 · (2 en la diagonal de bloque, 1 fuera)``, que es la masa
consistente exacta del tetraedro P1 (no la masa concentrada: importa para las frecuencias).

Congruencia (lo que se verifica en los tests)
---------------------------------------------
- Los seis **modos de cuerpo rígido** (3 traslaciones + 3 rotaciones) tienen que estar en el núcleo de
  ``K``: si la matriz no los anula, hay un error en ``B`` o en el ensamblado.
- La **energía de una deformación uniforme** tiene que dar ``½·V·σ:ε`` exacto, porque P1 reproduce
  campos lineales sin error. Es un test de consistencia, no de convergencia: no debería tener error
  de discretización, sólo de redondeo.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from ..materials import ElasticSolid

__all__ = [
    "lame_parameters",
    "elasticity_matrix",
    "tet_shape_gradients",
    "strain_matrix",
    "p1_elastic3d_matrices",
    "rigid_body_modes",
]


def lame_parameters(solid: ElasticSolid) -> tuple[float, float]:
    """``(λ, μ)`` desde ``ρ, c_p, c_s`` del sólido."""
    if not isinstance(solid, ElasticSolid):
        raise TypeError("solid tiene que ser un ElasticSolid")
    mu = solid.rho * solid.c_s**2
    lam = solid.rho * (solid.c_p**2 - 2.0 * solid.c_s**2)
    if mu <= 0:
        raise ValueError(f"el sólido tiene c_s = {solid.c_s!r}: μ <= 0 no es elástico")
    return float(lam), float(mu)


def elasticity_matrix(solid: ElasticSolid) -> np.ndarray:
    """Matriz constitutiva ``D`` (6×6) en Voigt con corte de ingeniería."""
    lam, mu = lame_parameters(solid)
    d = np.zeros((6, 6))
    d[:3, :3] = lam
    d[0, 0] = d[1, 1] = d[2, 2] = lam + 2.0 * mu
    d[3, 3] = d[4, 4] = d[5, 5] = mu
    return d


def tet_shape_gradients(nodes: np.ndarray, tets: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Gradientes de las funciones de forma y volúmenes, por tetraedro.

    Returns
    -------
    (grads, volumes) : ``grads`` con forma ``(M, 4, 3)`` — ``grads[e, i]`` es ``∇N_i`` del elemento
        ``e`` — y ``volumes`` ``(M,)`` positivos.
    """
    p = nodes[tets]                                     # (M, 4, 3)
    jac = np.stack([p[:, 1] - p[:, 0], p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]], axis=1)  # (M,3,3)
    det = np.linalg.det(jac)
    if np.any(np.abs(det) < 1e-300):
        raise ValueError("hay tetraedros degenerados (determinante nulo)")
    inv = np.linalg.inv(jac)
    # ∇N_{i+1} es la COLUMNA i de J⁻¹, no la fila: J·∇N = e_i. Trasponer es lo que hace que los
    # modos de cuerpo rígido caigan en el núcleo de K y que la energía de una deformación uniforme
    # dé exacta; con la transpuesta mal, el error es de pocos por ciento y parece "convergencia".
    grads = np.empty((tets.shape[0], 4, 3))
    grads[:, 1:, :] = np.transpose(inv, (0, 2, 1))
    # ∇N₀ = −(∇N₁+∇N₂+∇N₃). Se suman los tres gradientes YA calculados en vez de sumar `inv` por un
    # eje: sobre un array (M, 3, 3) el `axis` correcto no es el mismo que sobre una matriz suelta, y
    # equivocarse ahí no da error — da gradientes que no suman cero, o sea traslaciones con
    # deformación y una energía 5× más grande.
    grads[:, 0, :] = -grads[:, 1:, :].sum(axis=1)
    volumes = np.abs(det) / 6.0
    return grads, volumes


def strain_matrix(grads: np.ndarray) -> np.ndarray:
    """``B`` (6×12) por elemento desde los gradientes ``(M, 4, 3)``.

    Columnas en orden ``[u_x⁰, u_y⁰, u_z⁰, u_x¹, ...]``.
    """
    m = grads.shape[0]
    b = np.zeros((m, 6, 12))
    for i in range(4):
        dx, dy, dz = grads[:, i, 0], grads[:, i, 1], grads[:, i, 2]
        b[:, 0, 3 * i + 0] = dx
        b[:, 1, 3 * i + 1] = dy
        b[:, 2, 3 * i + 2] = dz
        b[:, 3, 3 * i + 1] = dz      # γ_yz = ∂v/∂z + ∂w/∂y
        b[:, 3, 3 * i + 2] = dy
        b[:, 4, 3 * i + 0] = dz      # γ_xz = ∂u/∂z + ∂w/∂x
        b[:, 4, 3 * i + 2] = dx
        b[:, 5, 3 * i + 0] = dy      # γ_xy = ∂u/∂y + ∂v/∂x
        b[:, 5, 3 * i + 1] = dx
    return b


def p1_elastic3d_matrices_constantes(
    nodes: np.ndarray, tets: np.ndarray, lam: float, mu: float, rho: float,
) -> tuple[sp.csr_matrix, sp.csr_matrix]:
    """``(K, M)`` desde ``(λ, μ, ρ)`` sueltos, sin pasar por un :class:`ElasticSolid`.

    Hace falta para armar el operador del **fluido** con la misma malla y la misma cuadratura que el
    del sólido (``μ = 0``, ``λ = ρ₀c₀²``), que es lo que permite escribir el término de contraste como
    diferencia de matrices y que se anule exactamente cuando los dos materiales coinciden.
    """
    d = np.zeros((6, 6))
    d[:3, :3] = lam
    d[0, 0] = d[1, 1] = d[2, 2] = lam + 2.0 * mu
    d[3, 3] = d[4, 4] = d[5, 5] = mu
    return _matrices(nodes, tets, d, float(rho))


def p1_elastic3d_matrices(
    nodes: np.ndarray, tets: np.ndarray, solid: ElasticSolid
) -> tuple[sp.csr_matrix, sp.csr_matrix]:
    """Matrices globales de rigidez y masa ``(K, M)`` de forma ``(3N, 3N)``.

    ``K`` incorpora únicamente la rigidez elástica: la superficie libre de tracción es la condición
    natural del método, así que un contorno sin término de borde queda con ``σ·n = 0``.
    """
    return _matrices(nodes, tets, elasticity_matrix(solid), float(solid.rho))


def _matrices(nodes, tets, d: np.ndarray, rho: float) -> tuple[sp.csr_matrix, sp.csr_matrix]:
    """Núcleo del ensamblado, desde la matriz constitutiva ``D`` y la densidad."""
    nodes = np.asarray(nodes, dtype=float)
    tets = np.asarray(tets, dtype=np.int64)
    if tets.ndim != 2 or tets.shape[1] != 4:
        raise ValueError(f"tets debe tener forma (M, 4); recibida {tets.shape}")
    n = nodes.shape[0]
    grads, volumes = tet_shape_gradients(nodes, tets)
    b = strain_matrix(grads)

    # K_e = V·BᵀDB  (M, 12, 12) vectorizado
    db = d @ b                                  # (M, 6, 12)
    ke = volumes[:, None, None] * np.einsum("mki,mkj->mij", b, db)

    filas_local = (3 * tets[:, :, None] + np.arange(3)[None, None, :]).reshape(len(tets), 12)
    filas = np.repeat(filas_local[:, :, None], 12, axis=2).reshape(-1)
    cols = np.repeat(filas_local[:, None, :], 12, axis=1).reshape(-1)
    k = sp.coo_matrix((ke.reshape(-1), (filas, cols)), shape=(3 * n, 3 * n)).tocsr()

    masa_local = np.zeros((12, 12))
    for i in range(4):
        for j in range(4):
            masa_local[3 * i:3 * i + 3, 3 * j:3 * j + 3] = np.eye(3) * (2.0 if i == j else 1.0)
    me = rho * volumes[:, None, None] / 20.0 * masa_local[None, :, :]
    m = sp.coo_matrix((me.reshape(-1), (filas, cols)), shape=(3 * n, 3 * n)).tocsr()
    return k.tocsr(), m.tocsr()


def rigid_body_modes(nodes: np.ndarray) -> np.ndarray:
    """Los seis modos de cuerpo rígido como columnas de una matriz ``(3N, 6)``.

    Tres traslaciones y tres rotaciones infinitesimales alrededor de los ejes, con el centroide como
    origen (así las rotaciones no arrastran traslación).
    """
    nodes = np.asarray(nodes, dtype=float)
    n = nodes.shape[0]
    centro = nodes.mean(axis=0)
    x = nodes - centro
    modos = np.zeros((3 * n, 6))
    for k in range(3):
        modos[k::3, k] = 1.0
    # rotación alrededor de x: u = (0, -z, y); alrededor de y: (z, 0, -x); de z: (-y, x, 0)
    modos[1::3, 3] = -x[:, 2]
    modos[2::3, 3] = x[:, 1]
    modos[0::3, 4] = x[:, 2]
    modos[2::3, 4] = -x[:, 0]
    modos[0::3, 5] = -x[:, 1]
    modos[1::3, 5] = x[:, 0]
    return modos
