"""Fase 4: contraste FEM vs BEM en el mismo problema acústico.

Problema: dispersión de una onda plana por un cilindro rígido (Neumann) en agua.

- **BEM** resuelve el exterior del cilindro de radio ``a`` (ecuación integral de frontera).
- **FEM** resuelve la corona ``a <= r <= R`` en CAMPO DISPERSADO: ABC en ``r = R`` y dato Neumann
  ``∂p_s/∂n = −∂p_i/∂n`` en ``r = a`` (superficie rígida). El campo total es ``p_s + p_i``.
- **Serie** (Jacobi-Anger + Hankel) es la verdad de referencia, para saber cuál de los dos se
  equivoca cuando no coinciden.

La formulación en campo total (``fem.solve``) NO sirve acá: en un contorno curvo la cara de sombra
recibe la ABC aplicada al campo total, sin el término incidente, y el error queda en 18 % sin
converger con la malla (medido: docs/decisions.md D12).

Se comparan los tres en el mismo conjunto de puntos (círculo de radio ``r_eval``), con el mismo
``k``. Se reporta error relativo L2 y tiempo de cada método.

Uso:
    .venv/bin/python scripts/cross_validate.py
    .venv/bin/python scripts/cross_validate.py --kas 1,2,5 --abc-order 2 --r-out 4
"""

from __future__ import annotations

import argparse
import time

import numpy as np
from scipy.interpolate import LinearNDInterpolator

from ultrasound_sim import materials, meshing
from ultrasound_sim.analytic import scattering
from ultrasound_sim.bem import helmholtz2d as bem
from ultrasound_sim.fem import helmholtz2d as fem

FREQUENCY = 100.0e3  # dentro del techo acordado (D8)
MATERIAL = materials.WATER


def eval_points(radius: float, n: int = 64) -> np.ndarray:
    th = 2.0 * np.pi * np.arange(n) / n
    return np.column_stack([radius * np.cos(th), radius * np.sin(th)])


def rel_l2(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


def fem_annulus_solution(ka: float, ppw: float, r_out: float, abc_order: int = 2):
    """FEM en campo dispersado sobre la corona. Devuelve (mesh, p_total, a, t, bordes)."""
    k = MATERIAL.wavenumber(FREQUENCY)
    a = ka / k.real
    R = r_out * a
    h = MATERIAL.wavelength(FREQUENCY) / ppw
    n_theta = max(16, int(np.ceil(2.0 * np.pi * R / h)))
    n_radial = max(2, int(np.ceil((R - a) / h)))
    mesh = meshing.polar_annulus_mesh(a, R, n_radial, n_theta)

    edges = meshing.boundary_edges(mesh.triangles)
    r_nodes = np.linalg.norm(mesh.nodes[edges[:, 0]], axis=1)
    outer = np.isclose(r_nodes, R, atol=1e-12)
    inner = np.isclose(r_nodes, a, atol=1e-12)

    k_el = np.full(mesh.n_elements, k, dtype=complex)
    p_inc = fem.plane_wave(FREQUENCY, MATERIAL, direction=(1.0, 0.0))

    # Dato Neumann en la superficie rígida: ∂p_s/∂n = −∂p_i/∂n con la normal saliente del dominio,
    # que en el agujero apunta hacia el eje.
    inner_edges = edges[inner]
    normals_inner = meshing.outward_normals(mesh.nodes, inner_edges)
    mid_inner = meshing.edge_midpoints(mesh.nodes, inner_edges)
    dp_dn_inc = p_inc.normal_derivative(mid_inner, normals_inner)

    t0 = time.perf_counter()
    p_s = fem.solve_scattered(
        mesh.nodes,
        mesh.triangles,
        k_el,
        abc_edges=edges[outer],
        neumann_edges=inner_edges,
        dp_dn_incident=dp_dn_inc,
        abc_order=abc_order,
        abc_center=(0.0, 0.0),
    )
    t_fem = time.perf_counter() - t0
    p_total = fem.total_field(p_s, p_inc, mesh.nodes)
    edges_count = (int(np.count_nonzero(inner)), int(np.count_nonzero(outer)))
    return mesh, p_total, a, t_fem, edges_count


def interpolate_fem(mesh, p: np.ndarray, pts: np.ndarray) -> np.ndarray:
    """Interpolación P1 (lineal por triángulo) del campo complejo en ``pts``."""
    interp = LinearNDInterpolator(
        mesh.nodes, np.column_stack([p.real, p.imag]), fill_value=np.nan
    )
    vals = interp(pts)
    if np.any(np.isnan(vals)):
        raise RuntimeError("puntos de evaluación fuera del dominio FEM")
    return vals[:, 0] + 1j * vals[:, 1]


def case(ka: float, ppw: float, r_out: float, n_panels: int | None,
         abc_order: int = 2) -> dict:
    mesh, p_fem, a, t_fem, edges_count = fem_annulus_solution(ka, ppw, r_out, abc_order)
    pts = eval_points(1.5 * a)

    t0 = time.perf_counter()
    sol = bem.solve_exterior_neumann(pts, FREQUENCY, MATERIAL, a, n_panels=n_panels)
    t_bem = time.perf_counter() - t0

    series = scattering.cylinder_scattering(FREQUENCY, MATERIAL, a, pts, boundary="rigid")
    p_true = series["total"]

    p_fem_at = interpolate_fem(mesh, p_fem, pts)
    p_bem_at = sol.p_total

    th = np.degrees(np.arctan2(pts[:, 1], pts[:, 0]))
    lit = np.abs(th) > 90.0  # cara iluminada: lado por donde entra la onda (x < 0)
    err_face = np.abs(p_fem_at - p_true) / np.abs(p_true).max()
    return dict(
        ka=ka, a=a, R=r_out * a, elems=mesh.n_elements, dofs=mesh.n_nodes,
        panels=sol.n_panels, cond=sol.condition_number,
        t_fem=t_fem, t_bem=t_bem,
        fem_bem=rel_l2(p_fem_at, p_bem_at),
        fem_series=rel_l2(p_fem_at, p_true),
        bem_series=rel_l2(p_bem_at, p_true),
        err_lit=float(err_face[lit].mean()),
        err_shadow=float(err_face[~lit].mean()),
        edges_count=edges_count,
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--kas", default="1,2")
    ap.add_argument("--ppw", type=float, default=20.0)
    ap.add_argument("--r-out", type=float, default=3.0)
    ap.add_argument("--panels", type=int, default=None)
    ap.add_argument("--abc-order", type=int, default=2, choices=(1, 2))
    args = ap.parse_args()

    kas = [float(x) for x in args.kas.split(",")]
    print(f"f = {FREQUENCY/1e3:g} kHz, {MATERIAL.name}, ppw={args.ppw:g}, R={args.r_out:g}a, "
          f"ABC orden {args.abc_order}, evaluación en r=1.5a (64 puntos)")
    print(" ka   a[mm]  R[mm]  elem   dofs  paneles  cond  t_FEM[s]  t_BEM[s]   "
          "|FEM-BEM|  |FEM-serie|  |BEM-serie|    ilum   sombra")
    for ka in kas:
        r = case(ka, args.ppw, args.r_out, args.panels, args.abc_order)
        print(
            f"{r['ka']:4.1f} {r['a']*1e3:6.3f} {r['R']*1e3:6.2f} {r['elems']:5d} "
            f"{r['dofs']:6d} {r['panels']:7d} {r['cond']:5.2f} {r['t_fem']:9.3f} "
            f"{r['t_bem']:8.3f}  {r['fem_bem']*100:7.3f}%  {r['fem_series']*100:8.3f}%  "
            f"{r['bem_series']*100:8.3f}%  {r['err_lit']*100:5.2f}%  {r['err_shadow']*100:5.2f}%"
        )


if __name__ == "__main__":
    main()
