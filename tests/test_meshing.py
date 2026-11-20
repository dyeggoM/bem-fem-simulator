"""Mallas: orientación CCW, normales salientes y nodos en las interfaces.

De estas propiedades dependen los signos de las integrales de borde de FEM y BEM: acá se fijan.
"""

import numpy as np
import pytest

from ultrasound_sim import meshing


def test_line_mesh_is_uniform_and_contiguous():
    m = meshing.line_mesh(0.0, 1.0, 4)
    assert m.n_elements == 4
    assert m.nodes.size == 5
    np.testing.assert_allclose(m.lengths, 0.25)
    np.testing.assert_array_equal(m.elements, [[0, 1], [1, 2], [2, 3], [3, 4]])


def test_layered_line_mesh_forces_nodes_on_interfaces():
    interfaces = [0.0, 0.002, 0.006, 0.01]
    m = meshing.layered_line_mesh(interfaces, max_step=0.0015)
    for x in interfaces:
        assert np.any(np.isclose(m.nodes, x, atol=1e-15)), f"falta nodo en la interfaz {x}"
    assert m.lengths.max() <= 0.0015 + 1e-15
    assert np.all(m.lengths > 0.0)


def test_rect_tri_mesh_is_ccw_and_covers_the_area():
    m = meshing.rect_tri_mesh(0.0, 2.0, 0.0, 1.0, nx=4, ny=3)
    assert m.n_nodes == 5 * 4
    assert m.n_elements == 2 * 4 * 3
    areas = m.areas()
    assert np.all(areas > 0.0), "todos los triángulos deben ser CCW"
    assert areas.sum() == pytest.approx(2.0 * 1.0)


def test_boundary_edges_of_a_rectangle_have_the_right_length():
    m = meshing.rect_tri_mesh(0.0, 2.0, 0.0, 1.0, nx=4, ny=3)
    edges = meshing.boundary_edges(m.triangles)
    # Perímetro = 2*(2+1) = 6
    length = np.sum(np.linalg.norm(m.nodes[edges[:, 0]] - m.nodes[edges[:, 1]], axis=1))
    assert length == pytest.approx(6.0)
    # Cada arista de frontera aparece una sola vez.
    keys = [tuple(sorted(e)) for e in edges]
    assert len(keys) == len(set(keys))


def test_outward_normals_point_out_of_the_domain():
    m = meshing.rect_tri_mesh(0.0, 2.0, 0.0, 1.0, nx=4, ny=3)
    edges = meshing.boundary_edges(m.triangles)
    normals = meshing.outward_normals(m.nodes, edges)
    mid = meshing.edge_midpoints(m.nodes, edges)
    center = np.array([1.0, 0.5])
    # En un rectángulo convexo, la normal saliente tiene producto positivo con (punto - centro).
    dots = np.sum(normals * (mid - center), axis=1)
    assert np.all(dots > 0.0)
    np.testing.assert_allclose(np.linalg.norm(normals, axis=1), 1.0)


def test_ccw_flip_repairs_a_cw_triangle():
    nodes = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    cw = np.array([[0, 2, 1]])
    fixed = meshing._ensure_ccw(nodes, cw)
    assert meshing.triangle_areas(nodes, fixed)[0] > 0.0


def test_polar_annulus_mesh_covers_the_annulus():
    m = meshing.polar_annulus_mesh(1.0, 3.0, n_radial=8, n_theta=96)
    assert np.all(m.areas() > 0.0), "todos los triángulos deben ser CCW"
    # Área analítica del anillo, con error de polígono despreciable a 96 tramos.
    assert m.areas().sum() == pytest.approx(np.pi * (3.0**2 - 1.0**2), rel=1e-3)
    radii = np.linalg.norm(m.nodes, axis=1)
    assert radii.min() == pytest.approx(1.0)
    assert radii.max() == pytest.approx(3.0)


def test_polar_annulus_has_two_boundaries_and_outward_normals():
    m = meshing.polar_annulus_mesh(2.0, 5.0, n_radial=4, n_theta=64)
    edges = meshing.boundary_edges(m.triangles)
    mid = meshing.edge_midpoints(m.nodes, edges)
    # Los nodos (no los puntos medios: la cuerda queda por dentro) marcan el radio del borde.
    r_nodes = np.linalg.norm(m.nodes[edges], axis=2)  # (E, 2)
    inner = np.all(np.isclose(r_nodes, 2.0, atol=1e-12), axis=1)
    outer = np.all(np.isclose(r_nodes, 5.0, atol=1e-12), axis=1)
    assert inner.sum() + outer.sum() == edges.shape[0], "sólo hay borde interior y exterior"
    assert inner.sum() == 64 and outer.sum() == 64
    # Perímetro de cada círculo (polígono inscripto: aproxima 2πr por abajo).
    for mask, radius in ((inner, 2.0), (outer, 5.0)):
        length = np.linalg.norm(
            m.nodes[edges[mask, 0]] - m.nodes[edges[mask, 1]], axis=1
        ).sum()
        assert length == pytest.approx(2.0 * np.pi * radius, rel=1e-3)
    # Normal saliente: hacia afuera en el borde exterior, hacia el eje en el interior.
    normals = meshing.outward_normals(m.nodes, edges)
    dot = np.sum(normals * mid, axis=1)  # producto con el vector posición
    assert np.all(dot[outer] > 0.0)
    assert np.all(dot[inner] < 0.0)
