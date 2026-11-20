"""Fase 4 del plan: contraste FEM vs BEM en el mismo problema acústico.

Este test es el que cierra la fase, y lo corre el orquestador sobre el árbol integrado: quien
implementó el FEM y quien implementó el BEM no firman el contraste entre ellos.

Problema común: onda plana sobre un cilindro rígido (Neumann) en agua.
- BEM: exterior del cilindro, 64 paneles (32 por longitud de onda a ka=2).
- FEM: corona a <= r <= 3a, campo dispersado, ABC de segundo orden en r=3a, dato Neumann en r=a.
- Serie de Bessel: verdad de referencia.

Criterio del plan: discrepancia < 2 % en campo, con mallas comparables, reportando error y tiempo.
"""

import time

import numpy as np
import pytest
from scipy.interpolate import LinearNDInterpolator

from ultrasound_sim import materials, meshing
from ultrasound_sim.analytic import scattering
from ultrasound_sim.bem import helmholtz2d as bem
from ultrasound_sim.fem import helmholtz2d as fem

FREQUENCY = 100.0e3
MATERIAL = materials.WATER
PPW = 20.0
R_OUT = 3.0
EVAL_POINTS = 64


def _eval_circle(radius, n=EVAL_POINTS):
    th = 2.0 * np.pi * np.arange(n) / n
    return np.column_stack([radius * np.cos(th), radius * np.sin(th)])


def _rel_l2(a, b):
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


def _annulus_fem(ka, abc_order=2):
    k = MATERIAL.wavenumber(FREQUENCY)
    a = ka / k.real
    R = R_OUT * a
    h = MATERIAL.wavelength(FREQUENCY) / PPW
    mesh = meshing.polar_annulus_mesh(
        a, R,
        n_radial=max(2, int(np.ceil((R - a) / h))),
        n_theta=max(16, int(np.ceil(2.0 * np.pi * R / h))),
    )
    edges = meshing.boundary_edges(mesh.triangles)
    r_nodes = np.linalg.norm(mesh.nodes[edges[:, 0]], axis=1)
    outer, inner = np.isclose(r_nodes, R, atol=1e-12), np.isclose(r_nodes, a, atol=1e-12)

    p_inc = fem.plane_wave(FREQUENCY, MATERIAL, direction=(1.0, 0.0))
    inner_edges = edges[inner]
    dp_dn_inc = p_inc.normal_derivative(
        meshing.edge_midpoints(mesh.nodes, inner_edges),
        meshing.outward_normals(mesh.nodes, inner_edges),
    )

    t0 = time.perf_counter()
    p_s = fem.solve_scattered(
        mesh.nodes, mesh.triangles, np.full(mesh.n_elements, k, dtype=complex),
        abc_edges=edges[outer], neumann_edges=inner_edges, dp_dn_incident=dp_dn_inc,
        abc_order=abc_order, abc_center=(0.0, 0.0),
    )
    t_fem = time.perf_counter() - t0
    return mesh, fem.total_field(p_s, p_inc, mesh.nodes), a, t_fem


def _interp_p1(mesh, p, pts):
    interp = LinearNDInterpolator(
        mesh.nodes, np.column_stack([p.real, p.imag]), fill_value=np.nan
    )
    v = interp(pts)
    assert not np.any(np.isnan(v)), "puntos de evaluación fuera del dominio FEM"
    return v[:, 0] + 1j * v[:, 1]


def _solve_both(ka, abc_order=2):
    mesh, p_fem, a, t_fem = _annulus_fem(ka, abc_order)
    pts = _eval_circle(1.5 * a)
    t0 = time.perf_counter()
    sol = bem.solve_exterior_neumann(pts, FREQUENCY, MATERIAL, a)
    t_bem = time.perf_counter() - t0
    true = scattering.cylinder_scattering(
        FREQUENCY, MATERIAL, a, pts, boundary="rigid"
    )["total"]
    return dict(
        fem=_interp_p1(mesh, p_fem, pts), bem=sol.p_total, true=true, pts=pts,
        a=a, mesh=mesh, t_fem=t_fem, t_bem=t_bem, bem_obj=sol,
    )


@pytest.mark.parametrize("ka", [1.0, 2.0])
def test_fem_and_bem_agree_to_within_two_percent(ka):
    r = _solve_both(ka)
    fem_bem = _rel_l2(r["fem"], r["bem"])
    fem_true = _rel_l2(r["fem"], r["true"])
    bem_true = _rel_l2(r["bem"], r["true"])
    print(
        f"\n  ka={ka:g}: |FEM-BEM|={fem_bem*100:.3f}%  |FEM-serie|={fem_true*100:.3f}%  "
        f"|BEM-serie|={bem_true*100:.3f}%  "
        f"(FEM {r['mesh'].n_elements} elem / {r['t_fem']*1e3:.1f} ms, "
        f"BEM {r['bem_obj'].n_panels} paneles / {r['t_bem']*1e3:.1f} ms)"
    )
    assert fem_bem < 0.02, f"FEM vs BEM = {fem_bem*100:.2f} % (criterio 2 %)"
    # Ninguno de los dos puede acertarle al otro y equivocarse: ambos contra la serie.
    assert fem_true < 0.02
    assert bem_true < 0.02


@pytest.mark.parametrize("ka", [1.0, 2.0])
def test_error_does_not_depend_on_the_side_of_the_scatterer(ka):
    """El bug de campo total daba 4 % en la cara iluminada y 24 % en la de sombra."""
    r = _solve_both(ka)
    th = np.degrees(np.arctan2(r["pts"][:, 1], r["pts"][:, 0]))
    lit = np.abs(th) > 90.0
    err = np.abs(r["fem"] - r["true"]) / np.abs(r["true"]).max()
    print(f"\n  ka={ka:g}: error medio iluminada={err[lit].mean()*100:.2f}%  "
          f"sombra={err[~lit].mean()*100:.2f}%")
    assert err[lit].mean() < 0.02
    assert err[~lit].mean() < 0.02


def test_first_order_abc_has_a_curvature_floor():
    """Registro del límite medido: la ABC de 1er orden en contorno curvo no baja de ~3 %.

    Si este número cambia, algo del ensamblado o de la geometría se movió: no es un criterio de
    aceptación, es un ancla contra deriva silenciosa (docs/decisions.md D12).
    """
    r = _solve_both(1.0, abc_order=1)
    err = _rel_l2(r["fem"], r["true"])
    print(f"\n  ka=1, ABC orden 1: |FEM-serie| = {err*100:.2f} % (piso de curvatura esperado ~6.5 %)")
    assert 0.03 < err < 0.10
    # Y con orden 2 sí entra en el criterio.
    r2 = _solve_both(1.0, abc_order=2)
    assert _rel_l2(r2["fem"], r2["true"]) < 0.02


def test_both_methods_scale_in_reasonable_time():
    """El contraste se corre en segundos: si un método se degrada, se ve acá."""
    r = _solve_both(2.0)
    assert r["t_fem"] < 5.0, f"FEM tardó {r['t_fem']:.2f} s"
    assert r["t_bem"] < 5.0, f"BEM tardó {r['t_bem']:.2f} s"
