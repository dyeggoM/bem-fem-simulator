"""Elastodinámica 2D P1 en plano deformación: propiedades de las matrices y del vector de tracción.

Se verifica lo que puede verificarse sin comparar contra otra solución de la misma formulación:

- ``K`` y ``M`` simétricas; ``M`` semidefinida positiva y con la masa total correcta
  (``Σ M_ij = 2·A`` para la traslación unitaria, porque ``M`` **no** lleva densidad);
- modos de cuerpo rígido exactos en el núcleo de ``K`` (traslación y rotación infinitesimal: ambas
  tienen ``ε ≡ 0`` y por lo tanto energía nula — es la prueba de que el ``B`` y el ``D`` están bien
  armados, incluido el factor 2 del corte de ingeniería);
- patch test de deformación uniforme: la energía de un campo lineal vale exactamente
  ``A·εᵀ·D·ε`` (P1 reproduce el campo y ``ε`` es constante por elemento);
- velocidad de onda: en una banda con ``u_y = 0`` en las paredes, los autovalores de ``(K, M)``
  tienen que dar ``ω = c_p·k`` con ``k = n·π/L`` (deformación uniaxial ⇒ ``c_p``, no ``c_placa``);
- tracción consistente: fuerza y momento totales, y el signo de ``normal_sign``.

Los números medidos se imprimen: correr con ``-q -s`` para verlos.
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.linalg as sla

from ultrasound_sim import meshing
from ultrasound_sim.fem.elastic2d import (
    p1_elastic_local,
    p1_elastic_matrices,
    plane_strain_matrix,
    strain_matrix,
    traction_vector,
)
from ultrasound_sim.materials import ElasticSolid

STEEL = ElasticSolid(
    name="placa_benchmark_acero",
    rho=7850.0,
    c_p=5900.0,
    c_s=3200.0,
    provenance="parámetro de benchmark del test, no es material de librería",
)


def _square(n: int = 4, side: float = 1.0):
    return meshing.rect_tri_mesh(0.0, side, 0.0, side, n, n)


# ---------------------------------------------------------------------------
# Matrices
# ---------------------------------------------------------------------------


def test_matrices_shape_symmetry_and_mass():
    mesh = _square(3)
    k, m = p1_elastic_matrices(mesh.nodes, mesh.triangles, STEEL)
    n_dofs = 2 * mesh.n_nodes
    assert k.shape == (n_dofs, n_dofs) and m.shape == (n_dofs, n_dofs)
    assert np.allclose(k.toarray(), k.toarray().T, atol=1e-9 * np.abs(k.toarray()).max())
    assert np.allclose(m.toarray(), m.toarray().T, atol=1e-12 * np.abs(m.toarray()).max())

    # M no lleva densidad: para la traslación unitaria (todas las componentes a 1),
    # uᵀ·M·u = 2·∫ dA = 2·A_total
    ones = np.ones(n_dofs)
    total_area = mesh.areas().sum()
    assert np.isclose(ones @ (m @ ones), 2.0 * total_area, rtol=1e-12)
    print(f"\nM: uᵀMu(traslación) = {ones @ (m @ ones):.6f} vs 2·A = {2 * total_area:.6f}")

    # M semidefinida positiva
    eig = np.linalg.eigvalsh(m.toarray())
    assert eig.min() > -1e-10 * eig.max()


def test_rigid_body_modes_are_in_the_kernel():
    mesh = _square(3)
    k, _ = p1_elastic_matrices(mesh.nodes, mesh.triangles, STEEL)
    kd = k.toarray()
    x, y = mesh.nodes[:, 0], mesh.nodes[:, 1]
    scale = np.abs(kd).max()

    trans_x = np.zeros(2 * mesh.n_nodes)
    trans_y = np.zeros(2 * mesh.n_nodes)
    trans_x[0::2] = 1.0
    trans_y[1::2] = 1.0
    rot = np.zeros(2 * mesh.n_nodes)
    rot[0::2] = -y  # u = (−y, x): rotación infinitesimal, ε ≡ 0
    rot[1::2] = x

    for name, vec in (("traslación x", trans_x), ("traslación y", trans_y), ("rotación", rot)):
        residual = np.linalg.norm(kd @ vec) / (scale * np.linalg.norm(vec))
        print(f"K·u ({name}): residuo relativo = {residual:.3e}")
        assert residual < 1e-12


def test_uniform_strain_patch_test():
    """Energía de un campo lineal = A·εᵀ·D·ε (exacta para P1)."""
    mesh = _square(4, side=0.7)
    k, _, _ = p1_elastic_local(mesh.nodes, mesh.triangles, STEEL)
    x, y = mesh.nodes[:, 0], mesh.nodes[:, 1]
    a, b, c, d = 1.0e-3, -2.0e-3, 0.5e-3, 3.0e-3
    u = np.zeros((mesh.n_nodes, 2))
    u[:, 0] = a * x + b * y
    u[:, 1] = c * x + d * y
    local = u[mesh.triangles].reshape(mesh.n_elements, 6)
    energy = np.einsum("mi,mij,mj->", local, k, local)

    eps = np.array([a, d, b + c])  # (ε_xx, ε_yy, γ_xy)
    expected = mesh.areas().sum() * eps @ plane_strain_matrix(STEEL) @ eps
    print(f"\nenergía FEM = {energy:.10e}, A·εᵀDε = {expected:.10e}, "
          f"rel = {abs(energy - expected) / abs(expected):.2e}")
    assert np.isclose(energy, expected, rtol=1e-12)

    # B sobre el triángulo unitario reproduce (ε_xx, ε_yy, γ_xy) de un campo lineal
    grads = np.array([[[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]]])
    bmat = strain_matrix(grads)
    u_local = np.array([[0.0, 0.0, a, c, b, d]])  # u(0,0), u(1,0), u(0,1)
    got = bmat[0] @ u_local[0]
    print(f"B·u = {got} vs (ε_xx, ε_yy, γ_xy) = {eps}")
    assert np.allclose(got, eps, rtol=1e-13, atol=1e-18)


def test_longitudinal_wave_speed_is_cp():
    """Banda con u_y = 0 en las paredes: ω = c_p·n·π/L (deformación uniaxial)."""
    length, height = 1.0, 0.05
    nx, ny = 60, 2
    mesh = meshing.rect_tri_mesh(0.0, length, 0.0, height, nx, ny)
    k, m = p1_elastic_matrices(mesh.nodes, mesh.triangles, STEEL)
    kd, md = k.toarray(), m.toarray()

    x, y = mesh.nodes[:, 0], mesh.nodes[:, 1]
    tol = 1e-12 * height
    wall = (np.isclose(y, 0.0, atol=tol) | np.isclose(y, height, atol=tol))
    uy_dofs = 2 * np.flatnonzero(wall) + 1
    assert uy_dofs.size == 2 * (nx + 1)
    free = np.setdiff1d(np.arange(2 * mesh.n_nodes), uy_dofs)

    w2 = sla.eigh(kd[np.ix_(free, free)], md[np.ix_(free, free)], eigvals_only=True)
    # M no lleva densidad ⇒ el autovalor es ρ·ω²; la traslación rígida queda en ~0
    w = np.sqrt(np.sort(w2[w2 > 0.0]) / STEEL.rho)
    w = w[w > 1e2]
    expected = STEEL.c_p * np.pi / length
    rel = np.abs(w[:4] - expected * np.arange(1, 5)) / (expected * np.arange(1, 5))
    print(f"\nω₁ medida = {w[0]:.2f} rad/s, c_p·π/L = {expected:.2f} rad/s, "
          f"rel = {rel[0]:.2e}; tres modos: {np.round(rel[:3], 6)}")
    assert rel[0] < 5e-3
    # sin el rodillo el modo sería de placa (c_placa ≈ 0.92·c_p en acero): verifico que NO es ese
    nu = (STEEL.c_p**2 - 2 * STEEL.c_s**2) / (2 * (STEEL.c_p**2 - STEEL.c_s**2))
    c_plate = STEEL.c_p * np.sqrt((1 - 2 * nu) / (1 - nu) ** 2)
    assert abs(w[0] - c_plate * np.pi / length) > 0.05 * expected


# ---------------------------------------------------------------------------
# Tracción
# ---------------------------------------------------------------------------


def test_traction_vector_force_and_moment():
    mesh = _square(4, side=0.6)
    side = 0.6
    edges = meshing.boundary_edges(mesh.triangles)
    mid = 0.5 * (mesh.nodes[edges[:, 0]] + mesh.nodes[edges[:, 1]])
    right = edges[np.isclose(mid[:, 0], side, atol=1e-12)]
    assert right.shape[0] == 4

    t = np.array([3.0e5, -1.0e5])  # Pa
    f = traction_vector(mesh.nodes, right, t)
    total = np.array([f[0::2].sum(), f[1::2].sum()])
    print(f"\nfuerza total = {total} vs t·L = {t * side}")
    assert np.allclose(total, t * side, rtol=1e-12)

    x, y = mesh.nodes[:, 0], mesh.nodes[:, 1]
    moment = np.sum(x * f[1::2] - y * f[0::2])
    # ∫(x·t_y − y·t_x) ds sobre x = side:  x ≡ side, ∫y ds = side²/2
    expected_moment = side * t[1] * side - t[0] * (side * side / 2.0)
    print(f"momento = {moment:.6f} vs {expected_moment:.6f}")
    assert np.isclose(moment, expected_moment, rtol=1e-12)

    f_neg = traction_vector(mesh.nodes, right, t, normal_sign=-1)
    assert np.allclose(f_neg, -f)


def test_traction_vector_broadcast_and_validations():
    mesh = _square(1, side=1.0)
    edges = meshing.boundary_edges(mesh.triangles)[:1]
    f = traction_vector(mesh.nodes, edges, np.array([1.0, 0.0]))
    assert np.isclose(f[0::2].sum(), 1.0)  # |t|·L con L = 1

    with pytest.raises(ValueError):
        traction_vector(mesh.nodes, edges, np.array([1.0, 0.0, 0.0]))
    with pytest.raises(ValueError):
        traction_vector(mesh.nodes, edges, np.array([[1.0, 0.0], [1.0, 0.0]]))
    with pytest.raises(ValueError):
        traction_vector(mesh.nodes, edges, np.array([1.0, 0.0]), normal_sign=2)
    with pytest.raises(ValueError):  # arista con nodos repetidos
        traction_vector(mesh.nodes, np.array([[0, 0]]), np.array([1.0, 0.0]))
    with pytest.raises(ValueError):  # índice fuera de rango
        traction_vector(mesh.nodes, np.array([[0, mesh.n_nodes]]), np.array([1.0, 0.0]))


# ---------------------------------------------------------------------------
# Validaciones de entrada
# ---------------------------------------------------------------------------


def test_validations():
    mesh = _square(1)
    with pytest.raises(TypeError):
        p1_elastic_matrices(mesh.nodes, mesh.triangles, object())
    flipped = mesh.triangles[:, [0, 2, 1]]
    with pytest.raises(ValueError, match="CCW"):
        p1_elastic_matrices(mesh.nodes, flipped, STEEL)
    with pytest.raises(ValueError):
        p1_elastic_matrices(mesh.nodes, np.array([[0, 1, 99]]), STEEL)
    with pytest.raises(ValueError):
        p1_elastic_matrices(mesh.nodes[:, :1], mesh.triangles, STEEL)
    with pytest.raises(ValueError):
        strain_matrix(np.zeros((2, 3, 3)))
