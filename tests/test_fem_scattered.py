"""Fase 4 (pre): FEM de **campo dispersado** contra la serie de Bessel del cilindro rígido.

El FEM resuelve ``p_s = p − p_i`` en el anillo ``a <= r <= R`` con Neumann rígido en ``r = a`` y ABC
en ``r = R``; el campo físico se reconstruye como ``p_total = p_s + p_i`` y se compara con
``analytic.scattering.cylinder_scattering`` en un círculo de muestreo interior.

Por qué el campo dispersado y no el total: con la formulación de campo total
(``fem.helmholtz2d.solve``) el término ``∂p_i/∂n − i·k·p_i`` de las caras de sombra no se anula
sobre un contorno **curvo** (``n·d ≠ 1``), y el error —17–36 % en la cara de sombra según la malla,
18–19 % global— no baja con ``ppw`` porque es un error de modelo. Ese caso está medido acá abajo
(``test_total_field_formulation_fails_on_curved_boundary``) sólo como diagnóstico: no se exige que
pase.

Los números medidos se imprimen con ``print``: correr con ``-q -s`` para verlos.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.special import hankel1

from ultrasound_sim import meshing
from ultrasound_sim.analytic.scattering import (
    cylinder_radial_derivative,
    cylinder_scattering,
)
from ultrasound_sim.fem.assembly import p1_boundary_mass
from ultrasound_sim.fem.helmholtz2d import (
    assemble_scattered,
    plane_wave,
    solve,
    solve_scattered,
    total_field,
)
from ultrasound_sim.materials import WATER, Material

FREQUENCY = 1.0e5  # 100 kHz, el techo de las fases 0-5 (docs/decisions.md D8)
LAMBDA = WATER.wavelength(FREQUENCY)
PPW = 20
N_SAMPLES = 64
SAMPLE_RADIUS_RATIO = 1.5


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _annulus(radius: float, outer_ratio: float = 3.0, ppw: int = PPW):
    """Malla polar ``radius <= r <= outer_ratio·radius`` con sus dos fronteras separadas.

    ``ppw`` se aplica al paso radial y al arco del contorno exterior; ``n_radial`` se fuerza par
    para que ``r = 1.5a`` no caiga justo sobre una línea de nodos (el muestreo interpola).
    """
    outer_radius = outer_ratio * radius
    h = LAMBDA / ppw
    n_radial = max(2, int(round((outer_radius - radius) / h)))
    n_radial += n_radial % 2
    n_theta = max(12, int(round(2.0 * np.pi * outer_radius / h)))
    mesh = meshing.polar_annulus_mesh(radius, outer_radius, n_radial, n_theta)

    edges = meshing.boundary_edges(mesh.triangles)
    r_nodes = np.linalg.norm(mesh.nodes[edges], axis=2)  # (B, 2)
    tol = 1e-9 * outer_radius
    inner = edges[np.all(np.abs(r_nodes - radius) <= tol, axis=1)]
    outer = edges[np.all(np.abs(r_nodes - outer_radius) <= tol, axis=1)]
    if inner.size == 0 or outer.size == 0:
        raise AssertionError("la malla no separó las fronteras interior/exterior")
    return mesh, inner, outer


def _sample_circle(radius: float, n: int = N_SAMPLES) -> np.ndarray:
    theta = 2.0 * np.pi * np.arange(n) / n
    r = SAMPLE_RADIUS_RATIO * radius
    return np.column_stack([r * np.cos(theta), r * np.sin(theta)])


def _p1_interpolate(mesh, values: np.ndarray, points: np.ndarray) -> np.ndarray:
    """Interpola los valores nodales P1 en ``points`` (búsqueda del triángulo contenedor)."""
    tris = mesh.triangles
    p0 = mesh.nodes[tris[:, 0]]
    v0 = mesh.nodes[tris[:, 1]] - p0
    v1 = mesh.nodes[tris[:, 2]] - p0
    d00 = np.einsum("ij,ij->i", v0, v0)
    d01 = np.einsum("ij,ij->i", v0, v1)
    d11 = np.einsum("ij,ij->i", v1, v1)
    det = d00 * d11 - d01 * d01

    out = np.empty(points.shape[0], dtype=complex)
    for i, q in enumerate(points):
        v2 = q[None, :] - p0
        d20 = np.einsum("ij,ij->i", v2, v0)
        d21 = np.einsum("ij,ij->i", v2, v1)
        b1 = (d11 * d20 - d01 * d21) / det
        b2 = (d00 * d21 - d01 * d20) / det
        b0 = 1.0 - b1 - b2
        inside = np.flatnonzero((b0 >= -1e-9) & (b1 >= -1e-9) & (b2 >= -1e-9))
        if inside.size == 0:
            raise AssertionError(f"el punto {q} cayó fuera de la malla")
        t = inside[0]
        out[i] = (b0[t] * values[tris[t, 0]] + b1[t] * values[tris[t, 1]]
                  + b2[t] * values[tris[t, 2]])
    return out


def _solve_scattered_case(
    ka: float,
    outer_ratio: float = 3.0,
    ppw: int = PPW,
    material: Material = WATER,
    abc_order: int = 2,
    flip_incident_normal: bool = False,
):
    """Resuelve la dispersión por cilindro rígido y devuelve ``(mesh, radius, p_total)``."""
    k = material.wavenumber(FREQUENCY)
    radius = ka / k.real
    mesh, inner, outer = _annulus(radius, outer_ratio, ppw)
    incident = plane_wave(FREQUENCY, material, (1.0, 0.0))

    # Normal SALIENTE DEL DOMINIO FEM en el agujero: apunta al eje del cilindro (n = −r̂).
    n_inner = meshing.outward_normals(mesh.nodes, inner)
    mid = meshing.edge_midpoints(mesh.nodes, inner)
    dp_dn_incident = incident.normal_derivative(mid, n_inner)
    if flip_incident_normal:
        dp_dn_incident = -dp_dn_incident

    p_scattered = solve_scattered(
        mesh.nodes,
        mesh.triangles,
        np.full(mesh.n_elements, k, dtype=complex),
        abc_edges=outer,
        neumann_edges=inner,
        dp_dn_incident=dp_dn_incident,
        element_rho=material.rho,
        abc_order=abc_order,
    )
    return mesh, radius, total_field(p_scattered, incident, mesh.nodes)


def _errors_by_face(mesh, radius: float, p_total: np.ndarray, material: Material = WATER):
    """Error relativo L2 del campo total: global, cara iluminada (x<0) y cara de sombra (x>0)."""
    points = _sample_circle(radius)
    exact = cylinder_scattering(FREQUENCY, material, radius, points)["total"]
    fem = _p1_interpolate(mesh, p_total, points)
    error = fem - exact
    illuminated = points[:, 0] < 0.0

    def relative(mask) -> float:
        return float(np.linalg.norm(error[mask]) / np.linalg.norm(exact[mask]))

    return relative(np.ones_like(illuminated)), relative(illuminated), relative(~illuminated)


# ---------------------------------------------------------------------------
# Signo del dato de Neumann: ∂p_s/∂n = −∂p_i/∂n con la normal saliente del dominio
# ---------------------------------------------------------------------------


def test_incident_normal_derivative_sign_matches_series():
    """La normal del agujero es ``−r̂`` y el dato ``∂p_i/∂n`` cierra contra la serie de Bessel.

    Es la verificación que pide el invariante: si la normal se toma al revés, ``b`` entra con el
    signo contrario y el campo dispersado sale invertido (≈ 80 % de error en el total).
    """
    ka = 2.0
    k = WATER.wavenumber(FREQUENCY)
    radius = ka / k.real
    mesh, inner, _ = _annulus(radius, 3.0)

    inner_nodes = np.unique(inner)
    coords = mesh.nodes[inner_nodes]
    rhat_nodes = coords / np.linalg.norm(coords, axis=1, keepdims=True)

    edge_centres = meshing.edge_midpoints(mesh.nodes, inner)
    rhat_edges = edge_centres / np.linalg.norm(edge_centres, axis=1, keepdims=True)
    n_domain = meshing.outward_normals(mesh.nodes, inner)

    incident = plane_wave(FREQUENCY, WATER, (1.0, 0.0))
    dp_i_dr = incident.normal_derivative(coords, rhat_nodes)         # ∂p_i/∂r
    dp_total_dr = cylinder_radial_derivative(FREQUENCY, WATER, radius, coords)
    dp_s_dr = dp_total_dr - dp_i_dr                                  # ∂p_s/∂r
    dp_s_dn = -dp_s_dr                                               # ∂/∂n = −∂/∂r en el agujero
    dp_i_dn = incident.normal_derivative(coords, -rhat_nodes)        # ∂p_i/∂n

    deviation = float(np.abs(n_domain + rhat_edges).max())
    rigid_residual = float(np.abs(dp_total_dr).max())
    identity_error = float(np.abs(dp_s_dn + dp_i_dn).max() / np.abs(dp_i_dn).max())
    scale = float(np.abs(dp_i_dn).max())

    print("\n[dato de Neumann en r = a]")
    print(f"  max |n_FEM + r̂|           = {deviation:.2e}  (la normal del agujero es −r̂)")
    print(f"  max |∂p_total/∂r| en a    = {rigid_residual:.2e}  (pared rígida en la serie)")
    print(f"  max |∂p_i/∂n|             = {scale:.4e} Pa/m")
    print(f"  max |∂p_s/∂n + ∂p_i/∂n|   = {identity_error:.2e}  (identidad ∂p_s/∂n = −∂p_i/∂n)")

    assert deviation < 1e-12
    assert rigid_residual < 1e-9 * max(scale, 1.0)
    assert identity_error < 1e-9


def test_flipped_incident_normal_breaks_the_solution():
    """Control negativo: invertir el dato de Neumann arruina el contraste (guarda anti-regresión)."""
    mesh, radius, p_total = _solve_scattered_case(1.0, outer_ratio=2.0)
    good = _errors_by_face(mesh, radius, p_total)[0]

    mesh_f, radius_f, p_total_f = _solve_scattered_case(
        1.0, outer_ratio=2.0, flip_incident_normal=True
    )
    bad = _errors_by_face(mesh_f, radius_f, p_total_f)[0]

    print(f"\n[control de signo] error con ∂p_i/∂n correcto = {100 * good:.2f} %"
          f"   con el signo invertido = {100 * bad:.2f} %")

    assert bad > 10.0 * good
    assert bad > 0.2


# ---------------------------------------------------------------------------
# (1) Contraste con la serie: ka = 1 y ka = 2, R = 3a, ppw = 20
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("ka", [1.0, 2.0])
def test_scattered_field_matches_cylinder_series(ka):
    """Error L2 del campo total < 2 % global **y en cada cara** (iluminada y de sombra)."""
    mesh, radius, p_total = _solve_scattered_case(ka, outer_ratio=3.0, ppw=PPW, abc_order=2)
    total, lit, shadow = _errors_by_face(mesh, radius, p_total)

    print(f"\n[dispersión rígida] ka = {ka:g}, R = 3a, ppw = {PPW},"
          f" nodos = {mesh.n_nodes}, elementos = {mesh.n_elements}")
    print(f"  error L2 del campo total (r = 1.5a, {N_SAMPLES} puntos) = {100 * total:.3f} %")
    print(f"  cara iluminada (x < 0) = {100 * lit:.3f} %   cara de sombra (x > 0) = "
          f"{100 * shadow:.3f} %")

    assert total < 2.0e-2
    assert lit < 2.0e-2
    assert shadow < 2.0e-2
    # el error ya no depende de la cara (antes: 4.3 % vs 17-24 %)
    assert abs(shadow - lit) < 1.5e-2


# ---------------------------------------------------------------------------
# (2) La ABC de primer orden deja un error que NO depende de la cara
# ---------------------------------------------------------------------------


def test_first_order_abc_error_is_face_symmetric():
    """Con la ABC de impedancia (la de ``assemble``) el error residual es simétrico entre caras.

    Mide ~6 % a ``kR = 3`` y no baja con ``ppw``: es la reflexión de la ABC de primer orden frente a
    una onda **cilíndrica** (le falta el término de curvatura ``1/(2r)``). Lo que importa acá es que
    la asimetría iluminada/sombra desapareció: la formulación dispersada eliminó la fuente espuria
    de la cara de sombra.
    """
    mesh, radius, p_total = _solve_scattered_case(1.0, outer_ratio=3.0, ppw=PPW, abc_order=1)
    total, lit, shadow = _errors_by_face(mesh, radius, p_total)

    print(f"\n[ABC de primer orden, ka = 1, R = 3a] total = {100 * total:.2f} %,"
          f" iluminada = {100 * lit:.2f} %, sombra = {100 * shadow:.2f} %")

    assert shadow < 1.5 * lit + 2.0e-2       # antes la sombra era 4-6× la cara iluminada
    assert abs(shadow - lit) < 2.0e-2
    assert total < 1.0e-1


def test_first_order_abc_error_does_not_depend_on_ppw():
    """Refinar la malla no arregla el error de la ABC de primer orden: es reflexión, no ``h``.

    Deja separadas las dos causas: la formulación dispersada eliminó el error de *modelo* (la
    fuente espuria de la cara de sombra, que dependía de la cara); lo que queda con la ABC de
    impedancia es reflexión espuria del contorno, y por eso no baja al pasar de ppw = 20 a 40.
    """
    mesh_c, radius_c, p_c = _solve_scattered_case(1.0, outer_ratio=3.0, ppw=20, abc_order=1)
    coarse = _errors_by_face(mesh_c, radius_c, p_c)[0]
    mesh, radius, p_total = _solve_scattered_case(1.0, outer_ratio=3.0, ppw=40, abc_order=1)
    fine = _errors_by_face(mesh, radius, p_total)[0]

    print(f"\n[ABC de primer orden] ppw = 20: {100 * coarse:.2f} %   ppw = 40: {100 * fine:.2f} %"
          f"   (|Δ| = {100 * abs(fine - coarse):.2f} puntos)")

    assert abs(fine - coarse) < 0.15 * coarse


# ---------------------------------------------------------------------------
# (3) Convergencia al alejar el contorno absorbente
# ---------------------------------------------------------------------------


def test_error_decreases_when_abc_moves_away():
    """R = 2a, 3a, 4a con ka = 1 y ppw = 20: el error baja al alejar el contorno."""
    ratios = (2.0, 3.0, 4.0)
    errors = []
    print("\n[convergencia con la posición de la ABC] ka = 1, ppw = 20")
    for ratio in ratios:
        mesh, radius, p_total = _solve_scattered_case(1.0, outer_ratio=ratio, ppw=PPW)
        total, lit, shadow = _errors_by_face(mesh, radius, p_total)
        errors.append(total)
        print(f"  R = {ratio:g}a: {100 * total:.3f} %   (iluminada {100 * lit:.3f} %,"
              f" sombra {100 * shadow:.3f} %, nodos = {mesh.n_nodes})")

    assert errors[1] < errors[0]
    assert errors[2] < errors[1]


# ---------------------------------------------------------------------------
# (4) Medio con atenuación (k complejo)
# ---------------------------------------------------------------------------


def test_attenuating_medium_matches_series():
    """ka = 2 con pérdidas (α(100 kHz) = 20 Np/m): el total sigue a la serie con k complejo."""
    lossy = Material(
        name="lossy_test",
        c=1500.0,
        rho=1000.0,
        alpha0=200.0,        # Np/m a 1 MHz
        alpha_power=1.0,
        f_ref=1.0e6,
        provenance="test: α(1 MHz) = 200 Np/m lineal ⇒ 20 Np/m a 100 kHz; no es un material real",
    )
    alpha = lossy.attenuation(FREQUENCY)
    k = lossy.wavenumber(FREQUENCY)

    mesh, radius, p_total = _solve_scattered_case(
        2.0, outer_ratio=3.0, ppw=PPW, material=lossy, abc_order=2
    )
    total, lit, shadow = _errors_by_face(mesh, radius, p_total, material=lossy)

    print(f"\n[medio atenuante] ka = 2, α(100 kHz) = {alpha:.1f} Np/m, Im k = {k.imag:.1f} 1/m")
    print(f"  error total = {100 * total:.3f} %  (iluminada {100 * lit:.3f} %,"
          f" sombra {100 * shadow:.3f} %)")

    assert alpha > 1.0  # el test tiene sentido sólo si la atenuación pesa
    assert total < 2.0e-2


# ---------------------------------------------------------------------------
# (5) La ABC de orden 2 reproduce el DtN exacto de una circunferencia
# ---------------------------------------------------------------------------


def test_abc_order2_operator_matches_exact_dtn():
    """El bloque de borde de ``abc_order=2`` coincide con ``k·H_m⁽¹⁾′/H_m⁽¹⁾`` en ``r = R``.

    El operador de la ABC es ``A_abc = −λ·(masa de borde)``; su respuesta a un modo
    ``exp(i·m·θ)`` debe ser el factor ``λ = k·H_m′(kR)/H_m(kR)`` del contorno circular (DtN exacto).
    Es la justificación independiente de los coeficientes de curvatura: no están elegidos para
    que el test de dispersión pase, sino que son la expansión asintótica de ese cociente.
    """
    k = WATER.wavenumber(FREQUENCY)
    radius = 1.0 / k.real
    outer_radius = 3.0 * radius
    mesh, _, outer = _annulus(radius, 3.0)
    kk = np.full(mesh.n_elements, k, dtype=complex)
    no_edges = np.zeros((0, 2), dtype=int)

    a_volume, _ = assemble_scattered(mesh.nodes, mesh.triangles, kk, abc_edges=no_edges)
    a_abc, _ = assemble_scattered(
        mesh.nodes, mesh.triangles, kk, abc_edges=outer, abc_order=2
    )
    boundary = a_abc - a_volume
    mass = p1_boundary_mass(mesh.nodes, outer)

    nodes = np.unique(outer)
    theta = np.arctan2(mesh.nodes[nodes, 1], mesh.nodes[nodes, 0])
    z = k * outer_radius
    print(f"\n[ABC de orden 2 vs DtN exacto] kR = {abs(z):.2f}")
    # Sólo los modos que llevan campo a kR = 3 (|m| <= 1) se exigen: la condición es una
    # expansión asintótica en 1/(kR) y el error crece con |m| (a m = 2 ya es ~30 % en la parte real).
    tolerances = {0: 5.0e-2, 1: 8.0e-2, 2: None}
    for order in (0, 1, 2):
        v = np.zeros(mesh.n_nodes, dtype=complex)
        v[nodes] = np.exp(1j * order * theta)
        lam = complex(np.mean(-(boundary @ v)[nodes] / (mass @ v)[nodes]))
        exact = k * 0.5 * (hankel1(order - 1, z) - hankel1(order + 1, z)) / hankel1(order, z)
        error = float(abs(lam - exact) / abs(exact))
        print(f"  m = {order}: λ_ABC = {lam:.3f}   λ_exacto = {exact:.3f}"
              f"   error = {100 * error:.2f} %")
        if tolerances[order] is not None:
            assert error < tolerances[order]


# ---------------------------------------------------------------------------
# (6) Diagnóstico: la formulación de campo total falla en el contorno curvo
# ---------------------------------------------------------------------------


def test_total_field_formulation_fails_on_curved_boundary():
    """Deja registrado el bug que motiva la formulación dispersada (no es un criterio a cumplir)."""
    ka = 1.0
    k = WATER.wavenumber(FREQUENCY)
    radius = ka / k.real
    mesh, inner, outer = _annulus(radius, 3.0)
    incident = plane_wave(FREQUENCY, WATER, (1.0, 0.0))

    p_total = solve(
        mesh.nodes,
        mesh.triangles,
        np.full(mesh.n_elements, k, dtype=complex),
        abc_edges=outer,
        p_inc=incident,
        element_rho=WATER.rho,
    )
    total, lit, shadow = _errors_by_face(mesh, radius, p_total)
    print(f"\n[campo total, contorno curvo] ka = 1, R = 3a: total = {100 * total:.1f} %,"
          f" iluminada = {100 * lit:.1f} %, sombra = {100 * shadow:.1f} %"
          "   ← bug de la máscara de iluminación")

    # El síntoma: la cara de sombra es mucho peor que la iluminada y el global supera el 10 %.
    assert total > 0.10
    assert shadow > 2.0 * lit


# ---------------------------------------------------------------------------
# (7) Validación de la API
# ---------------------------------------------------------------------------


def test_assemble_scattered_validates_arguments():
    k = WATER.wavenumber(FREQUENCY)
    radius = 1.0 / k.real
    mesh, inner, outer = _annulus(radius, 2.0, ppw=8)
    kk = np.full(mesh.n_elements, k, dtype=complex)
    mid = meshing.edge_midpoints(mesh.nodes, inner)
    dp = np.ones(inner.shape[0], dtype=complex)

    with pytest.raises(ValueError, match="dp_dn_incident"):
        assemble_scattered(mesh.nodes, mesh.triangles, kk, abc_edges=outer,
                           neumann_edges=inner)

    with pytest.raises(ValueError, match="dp_dn_incident"):
        assemble_scattered(mesh.nodes, mesh.triangles, kk, abc_edges=outer,
                           neumann_edges=inner, dp_dn_incident=dp[:-1])

    with pytest.raises(ValueError, match="dp_dn_incident se pasó sin neumann_edges"):
        assemble_scattered(mesh.nodes, mesh.triangles, kk, abc_edges=outer,
                           dp_dn_incident=dp)

    with pytest.raises(ValueError, match="no es una arista de frontera"):
        assemble_scattered(mesh.nodes, mesh.triangles, kk, abc_edges=outer,
                           neumann_edges=np.array([[0, 3]]), dp_dn_incident=dp[:1])

    with pytest.raises(ValueError, match="abc_order"):
        assemble_scattered(mesh.nodes, mesh.triangles, kk, abc_edges=outer,
                           abc_order=3)

    # solapamiento ABC / Neumann
    with pytest.raises(ValueError, match="se solapan"):
        assemble_scattered(mesh.nodes, mesh.triangles, kk, abc_edges=outer,
                           neumann_edges=outer, dp_dn_incident=np.ones(len(outer)))

    with pytest.raises(ValueError, match="element_k"):
        assemble_scattered(mesh.nodes, mesh.triangles, np.ones(3, dtype=complex),
                           abc_edges=outer)

    # el orden de los dos nodos de cada arista de Neumann es indistinto (se reorienta contra la
    # malla): pasar las aristas al revés tiene que dar exactamente el mismo sistema.
    _, b_direct = assemble_scattered(mesh.nodes, mesh.triangles, kk, abc_edges=outer,
                                     neumann_edges=inner, dp_dn_incident=dp)
    _, b_reversed = assemble_scattered(mesh.nodes, mesh.triangles, kk, abc_edges=outer,
                                       neumann_edges=inner[:, ::-1], dp_dn_incident=dp)
    assert np.allclose(b_direct, b_reversed)
    assert mid.shape == (inner.shape[0], 2)
