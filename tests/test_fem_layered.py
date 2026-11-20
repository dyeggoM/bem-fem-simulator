"""Fase 1: FEM acústico 2D — onda plana, medios estratificados y ABC de primer orden.

Criterios de aceptación (docs/plan.md, D6):

(a) Dominio homogéneo: FEM contra onda plana analítica, error relativo L2 < 1 %.
(b) Dos semiespacios con incidencia normal: R medido contra el analítico, dentro del 2 % en
    amplitud.
(c) Reflexión espuria de la ABC en dominio homogéneo por debajo de −40 dB.

Los números medidos se imprimen con ``print``: correr con ``-q -s`` para verlos.
"""

from __future__ import annotations

import numpy as np

from ultrasound_sim import meshing
from ultrasound_sim.analytic.helmholtz1d import stack_response
from ultrasound_sim.fem.assembly import p1_tri_matrices
from ultrasound_sim.fem.helmholtz2d import plane_wave, solve
from ultrasound_sim.materials import WATER, Material, fluid_from_impedance

FREQUENCY = 1.0e5  # 100 kHz: λ = 15 mm en agua
LAMBDA = WATER.wavelength(FREQUENCY)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _steel_like() -> Material:
    """Material de prueba con Z = ρ·c explícito (no se agrega a `materials.py`)."""
    rho = 7850.0
    c = 5900.0
    return fluid_from_impedance("steel_test", impedance=rho * c, rho=rho)


def _strip(
    half_length_in_lambda: float,
    points_per_wavelength: int,
    height_in_lambda: float = 0.5,
):
    """Rectángulo centrado en x = 0 con la ABC sólo en las caras de entrada/salida.

    Las caras superior e inferior quedan con Neumann natural (paredes rígidas): con incidencia
    normal el campo es independiente de ``y`` y esas paredes son exactas, mientras que una ABC
    rasante sí reflejaría.
    """
    h = LAMBDA / points_per_wavelength
    nx = int(round(2.0 * half_length_in_lambda * LAMBDA / h))
    nx += nx % 2  # nodo exacto en x = 0 para que la interfaz caiga sobre una columna de nodos
    ny = max(2, int(round(height_in_lambda * LAMBDA / h)))
    xmin = -half_length_in_lambda * LAMBDA
    xmax = +half_length_in_lambda * LAMBDA
    mesh = meshing.rect_tri_mesh(xmin, xmax, 0.0, height_in_lambda * LAMBDA, nx, ny)

    edges = meshing.boundary_edges(mesh.triangles)
    coords = mesh.nodes[edges]
    left = edges[np.all(np.isclose(coords[:, :, 0], xmin), axis=1)]
    right = edges[np.all(np.isclose(coords[:, :, 0], xmax), axis=1)]
    return mesh, np.vstack([left, right])


def _l2_relative(mesh, p: np.ndarray, exact: np.ndarray) -> float:
    mass, _ = p1_tri_matrices(mesh.nodes, mesh.triangles)
    err = p - exact
    num = np.real(np.conj(err) @ (mass @ err))
    den = np.real(np.conj(exact) @ (mass @ exact))
    return float(np.sqrt(num / den))


def _decompose(x_coord: np.ndarray, values: np.ndarray, k: complex) -> tuple[complex, complex]:
    """Ajuste ``a·exp(i·k·x) + b·exp(−i·k·x)``; devuelve ``(a, b)``."""
    design = np.column_stack([np.exp(1j * k * x_coord), np.exp(-1j * k * x_coord)])
    a, b = np.linalg.lstsq(design, values, rcond=None)[0]
    return complex(a), complex(b)


# ---------------------------------------------------------------------------
# Ensamblado: propiedades básicas (mallas no uniformes incluidas)
# ---------------------------------------------------------------------------


def test_assembly_properties_on_non_uniform_meshes():
    from ultrasound_sim.fem.assembly import (
        p1_boundary_mass,
        p1_line_matrices,
        p1_tri_matrices,
    )

    # 1D: nodos NO equiespaciados; suma de la masa = longitud y K anula las constantes.
    nodes_1d = np.array([0.0, 0.1, 0.35, 0.4, 1.0])
    elements_1d = np.column_stack([np.arange(4), np.arange(1, 5)])
    mass_1d, stiff_1d = p1_line_matrices(nodes_1d, elements_1d)
    length = float(np.diff(nodes_1d).sum())

    mesh = meshing.rect_tri_mesh(0.0, 2.0, 0.0, 1.0, 7, 3)  # quads no cuadrados
    mass_2d, stiff_2d = p1_tri_matrices(mesh.nodes, mesh.triangles)
    edges = meshing.boundary_edges(mesh.triangles)
    boundary = p1_boundary_mass(mesh.nodes, edges)
    perimeter = float(
        np.linalg.norm(
            mesh.nodes[edges[:, 1]] - mesh.nodes[edges[:, 0]], axis=1
        ).sum()
    )
    area = float(mesh.areas().sum())

    print("\n[ensamblado] 1D no uniforme: suma M = %.12f (longitud %.12f), max|K·1| = %.2e"
          % (mass_1d.sum(), length, np.abs(stiff_1d @ np.ones(5)).max()))
    print("  2D: suma M = %.12f (área %.12f), suma S = %.12f (perímetro %.12f)"
          % (mass_2d.sum(), area, boundary.sum(), perimeter))

    assert abs(mass_1d.sum() - length) < 1e-12
    assert np.abs(stiff_1d @ np.ones(5)).max() < 1e-12
    assert abs(mass_2d.sum() - area) < 1e-10
    assert abs(stiff_2d @ np.ones(mesh.n_nodes)).max() < 1e-10
    assert abs(boundary.sum() - perimeter) < 1e-10
    assert np.allclose(mass_2d, mass_2d.T)


# ---------------------------------------------------------------------------
# (a) Homogéneo: FEM vs onda plana
# ---------------------------------------------------------------------------


def test_homogeneous_plane_wave():
    points_per_wavelength = 60
    mesh, abc = _strip(1.5, points_per_wavelength)
    k = WATER.wavenumber(FREQUENCY)
    incident = plane_wave(FREQUENCY, WATER, (1.0, 0.0))

    p = solve(
        mesh.nodes,
        mesh.triangles,
        np.full(mesh.n_elements, k, dtype=complex),
        abc_edges=abc,
        p_inc=incident,
    )
    exact = incident.value(mesh.nodes)
    error = _l2_relative(mesh, p, exact)

    print(f"\n[2D homogéneo] nodos = {mesh.n_nodes}, elementos = {mesh.n_elements},"
          f" ppw = {points_per_wavelength}")
    print(f"  error relativo L2 vs onda plana = {error:.4e}  ({100 * error:.3f} %)")
    print(f"  max |p_FEM − p_inc| = {np.abs(p - exact).max():.4e}"
          f"   (|p_inc| = {np.abs(exact).max():.3f})")

    assert error < 1.0e-2


# ---------------------------------------------------------------------------
# (b) Dos semiespacios: R medido vs analítico
# ---------------------------------------------------------------------------


def test_two_half_spaces_reflection_matches_impedance():
    steel = _steel_like()
    k1 = WATER.wavenumber(FREQUENCY)
    k2 = steel.wavenumber(FREQUENCY)

    z1 = WATER.impedance
    z2 = steel.impedance
    r_analytic = (z2 - z1) / (z2 + z1)
    t_analytic = 2.0 * z2 / (z1 + z2)

    points_per_wavelength = 20
    mesh, abc = _strip(1.5, points_per_wavelength)
    x = mesh.nodes

    # Material heterogéneo por elemento: la interfaz cae sobre la columna de nodos de x = 0.
    right_side = mesh.centroids()[:, 0] > 0.0
    element_k = np.where(right_side, k2, k1).astype(complex)
    # densidad por elemento: hace falta el peso 1/ρ de la ecuación heterogénea
    element_rho = np.where(right_side, steel.rho, WATER.rho)

    incident = plane_wave(FREQUENCY, WATER, (1.0, 0.0))
    p = solve(
        mesh.nodes,
        mesh.triangles,
        element_k,
        abc_edges=abc,
        p_inc=incident,
        element_rho=element_rho,
    )

    a_left, b_left = _decompose(x[x[:, 0] < 0.0, 0], p[x[:, 0] < 0.0], k1)
    a_right, b_right = _decompose(x[x[:, 0] > 0.0, 0], p[x[:, 0] > 0.0], k2)
    r_fem = b_left / a_left
    t_fem = a_right / a_left

    error = abs(abs(r_fem) - r_analytic) / r_analytic
    error_t = abs(abs(t_fem) - t_analytic) / t_analytic

    print(f"\n[2D dos semiespacios] ppw = {points_per_wavelength}, nodos = {mesh.n_nodes}")
    print(f"  Z₂/Z₁ = {z2 / z1:.3f}")
    print(f"  |R| FEM = {abs(r_fem):.6f}   analítico = {r_analytic:.6f}"
          f"   error = {100 * error:.3f} %   fase = {np.degrees(np.angle(r_fem)):+.2f}°")
    print(f"  |T| FEM = {abs(t_fem):.6f}   analítico = {t_analytic:.6f}"
          f"   error = {100 * error_t:.3f} %")
    print(f"  |a_izq| = {abs(a_left):.6f} (debería ser 1),  |b_der|/|a_der| = {abs(b_right / a_right):.2e}")

    assert error < 2.0e-2
    assert abs(a_left) > 0.99  # la onda incidente no se atenúa en el dominio sin pérdidas


def test_single_layer_slab_matches_analytic_1d():
    """R y T de **una capa** contra el analítico 1D, dentro del 2 % en amplitud (plan, Fase 1)."""
    slab = Material(name="slab_test", c=2000.0, rho=1100.0, provenance="test")
    k1 = WATER.wavenumber(FREQUENCY)
    k2 = slab.wavenumber(FREQUENCY)

    points_per_wavelength = 40
    h = LAMBDA / points_per_wavelength
    d = 6.0 * h  # espesor múltiplo del paso: las dos interfaces caen sobre columnas de nodos
    mesh, abc = _strip(1.5, points_per_wavelength)
    x = mesh.nodes

    cx = mesh.centroids()[:, 0]
    element_k = np.where(cx < 0.0, k1, np.where(cx < d, k2, k1)).astype(complex)
    element_rho = np.where(cx < 0.0, WATER.rho, np.where(cx < d, slab.rho, WATER.rho))

    incident = plane_wave(FREQUENCY, WATER, (1.0, 0.0))
    p = solve(
        mesh.nodes,
        mesh.triangles,
        element_k,
        abc_edges=abc,
        p_inc=incident,
        element_rho=element_rho,
    )

    a_left, b_left = _decompose(x[x[:, 0] < 0.0, 0], p[x[:, 0] < 0.0], k1)
    a_right, b_right = _decompose(x[x[:, 0] > d, 0], p[x[:, 0] > d], k1)
    r_fem = b_left / a_left
    t_fem = a_right / a_left

    analytic = stack_response(FREQUENCY, [WATER, slab, WATER], (d,))
    error_r = abs(abs(r_fem) - abs(analytic.R)) / abs(analytic.R)
    error_t = abs(abs(t_fem) - abs(analytic.T)) / abs(analytic.T)

    print(f"\n[2D una capa] espesor = {1e3 * d:.3f} mm (λ_capa = {1e3 * slab.wavelength(FREQUENCY):.1f} mm),"
          f" ppw = {points_per_wavelength}, nodos = {mesh.n_nodes}")
    print(f"  R: FEM {r_fem:.6f}  analítico {analytic.R:.6f}  error en amplitud = {100 * error_r:.3f} %")
    print(f"  T: FEM {t_fem:.6f}  analítico {analytic.T:.6f}  error en amplitud = {100 * error_t:.3f} %")
    print(f"  |b_derecha/a_derecha| = {abs(b_right / a_right):.2e} (onda reflejada espuria en la salida)")

    assert error_r < 2.0e-2
    assert error_t < 2.0e-2


# ---------------------------------------------------------------------------
# (c) Reflexión espuria de la ABC
# ---------------------------------------------------------------------------


def test_abc_spurious_reflection_below_minus_40_db():
    points_per_wavelength = 20
    mesh, abc = _strip(1.5, points_per_wavelength)
    k = WATER.wavenumber(FREQUENCY)
    incident = plane_wave(FREQUENCY, WATER, (1.0, 0.0))

    p = solve(
        mesh.nodes,
        mesh.triangles,
        np.full(mesh.n_elements, k, dtype=complex),
        abc_edges=abc,
        p_inc=incident,
    )

    # En un dominio homogéneo el campo debe ser sólo la onda incidente: la amplitud de la onda
    # que viaja hacia −x es la reflexión espuria de la ABC.
    x = mesh.nodes
    window = np.abs(x[:, 0]) <= 0.75 * LAMBDA
    a, b = _decompose(x[window, 0], p[window], k)

    reflection_db = 20.0 * np.log10(abs(b / a))
    exact = incident.value(x)
    print(f"\n[2D ABC] ppw = {points_per_wavelength}, nodos = {mesh.n_nodes}")
    print(f"  amplitud incidente ajustada |a| = {abs(a):.6f}")
    print(f"  amplitud reflejada ajustada |b| = {abs(b):.3e}")
    print(f"  reflexión espuria de la ABC = {reflection_db:.1f} dB  (criterio < −40 dB)")
    print(f"  error L2 del campo vs onda plana = {_l2_relative(mesh, p, exact):.3e}")

    assert reflection_db < -40.0
