"""Escalera de frecuencia: costo y precisión del FEM 2D, para decidir qué optimizar.

Tres barridos, cada uno responde una pregunta distinta:

A) **Auto-similar** (dominio medido en λ): a igual resolución en longitudes de onda, ¿cuánto cuesta?
   Con 8 λ de lado y 10 ppw la malla es idéntica en todas las frecuencias; lo que cambia es el
   tamaño FÍSICO del dominio que entra con ese presupuesto. Esa es la respuesta honesta a
   "¿cuánto cuesta pasar de 40 kHz a 1 MHz?".

B) **Precisión vs ppw** en un ducto a 100 kHz (paredes Neumann arriba/abajo, ABC sólo en entrada y
   salida, onda plana exacta). Se mide el error L2 contra la solución exacta y su orden de
   convergencia. El ducto es imprescindible: con ABC en los cuatro lados la onda incide rasante en
   las paredes superior/inferior y la ABC de primer orden se rompe (error > 50 %), lo que mide la
   frontera y no el elemento.

C) **Dominio físico fijo de 20 mm** (el tamaño de un modelo real de cóclea): cómo crece la malla con
   la frecuencia en el escenario que de verdad se va a simular, con tiempos de ensamblado y de
   resolución, y comparación `spsolve` vs factorizar una vez con `splu`.

Uso:
    .venv/bin/python scripts/bench_frequency.py
    .venv/bin/python scripts/bench_frequency.py --only B
"""

from __future__ import annotations

import argparse
import time

import numpy as np
import scipy.sparse.linalg as spla

from ultrasound_sim import materials, meshing
from ultrasound_sim.fem import helmholtz2d


def build_square(frequency, side: float, ppw: float, material):
    """Malla cuadrada de lado ``side`` [m] con ``ppw`` puntos por longitud de onda."""
    h = material.wavelength(frequency) / ppw
    n = max(2, int(np.ceil(side / h)))
    mesh = meshing.rect_tri_mesh(0.0, side, 0.0, side, n, n)
    k = material.wavenumber(frequency)
    return mesh, n


def duct_edges(nodes, triangles, side: float) -> np.ndarray:
    """Aristas de entrada y salida de un ducto (x=0 y x=side); arriba/abajo quedan Neumann."""
    edges = meshing.boundary_edges(triangles)
    pts = nodes[edges]  # (E, 2, 2)
    left = np.all(np.isclose(pts[:, :, 0], 0.0, atol=1e-12), axis=1)
    right = np.all(np.isclose(pts[:, :, 0], side, atol=1e-12), axis=1)
    return edges[left | right]


def relative_l2(p_fem: np.ndarray, p_exact: np.ndarray) -> float:
    return float(np.linalg.norm(p_fem - p_exact) / np.linalg.norm(p_exact))


def solve_duct(frequency, side, ppw, material):
    """Resuelve el ducto y devuelve (mesh, p, err, t_asm, t_spsolve, t_splu, nnz)."""
    mesh, _ = build_square(frequency, side, ppw, material)
    k = material.wavenumber(frequency)
    k_el = np.full(mesh.n_elements, k, dtype=complex)
    p_inc = helmholtz2d.plane_wave(frequency, material, direction=(1.0, 0.0))
    edges = duct_edges(mesh.nodes, mesh.triangles, side)

    t0 = time.perf_counter()
    a, b = helmholtz2d.assemble(
        mesh.nodes, mesh.triangles, k_el, abc_edges=edges, p_inc=p_inc
    )
    t_asm = time.perf_counter() - t0
    a = a.tocsc()

    t0 = time.perf_counter()
    p = spla.spsolve(a, b)
    t_spsolve = time.perf_counter() - t0

    t0 = time.perf_counter()
    lu = spla.splu(a)
    p_lu = lu.solve(b)
    t_splu = time.perf_counter() - t0

    agree = float(np.linalg.norm(p - p_lu) / max(np.linalg.norm(p_lu), 1e-300))
    err = relative_l2(p, p_inc.value(mesh.nodes))
    return dict(mesh=mesh, p=p, err=err, t_asm=t_asm, t_spsolve=t_spsolve,
                t_splu=t_splu, nnz=a.nnz, agree=agree, dofs=mesh.n_nodes,
                elems=mesh.n_elements)


def sweep_a(freqs, lambdas, ppw, material):
    print(f"\nA) auto-similar: lado = {lambdas:g} λ, {ppw:.0f} ppw — el costo es el mismo, "
          f"el dominio físico se encoge")
    print("f[kHz]  lambda[mm]  lado[mm]  elems   dofs   t_spsolve[s]")
    rows = []
    for f in freqs:
        lam = material.wavelength(f)
        r = solve_duct(f, lambdas * lam, ppw, material)
        print(f"{f/1e3:6.0f}  {lam*1e3:10.2f}  {lambdas*lam*1e3:8.1f}  {r['elems']:5d}  "
              f"{r['dofs']:5d}   {r['t_spsolve']:11.3f}   err={r['err']*100:6.3f}%")
        rows.append(r)
    return rows


def sweep_b(frequency, side, material, ppws=(5, 10, 20, 40)):
    print(f"\nB) precisión vs ppw en ducto de {side*1e3:g} mm a {frequency/1e3:g} kHz "
          f"(error L2 vs onda plana exacta)")
    print(" ppw   elems    dofs   err_L2      orden")
    prev = None
    for ppw in ppws:
        r = solve_duct(frequency, side, ppw, material)
        order = ""
        if prev is not None:
            k_el_ratio = np.log(prev["err"] / r["err"]) / np.log(2.0)
            order = f"{k_el_ratio:6.2f}"
        print(f"{ppw:4.0f}  {r['elems']:6d}  {r['dofs']:6d}  {r['err']*100:8.4f}%  {order}")
        prev = r


def sweep_c(freqs, side, ppw, material):
    print(f"\nC) dominio físico fijo de {side*1e3:g} mm, {ppw:.0f} ppw — el escenario real")
    print("f[kHz]  elems    dofs       nnz  t_asm[s]  t_spsolve[s]  t_splu[s]   "
          "err_L2   agree_spsolve_splu")
    rows = []
    for f in freqs:
        r = solve_duct(f, side, ppw, material)
        print(f"{f/1e3:6.0f}  {r['elems']:6d}  {r['dofs']:6d}  {r['nnz']:8d}   "
              f"{r['t_asm']:7.3f}  {r['t_spsolve']:11.3f}  {r['t_splu']:8.3f}  "
              f"{r['err']*100:7.3f}%   {r['agree']:.2e}")
        rows.append(r)
    return rows


def lu_reuse(freqs, side, ppw, material, n_sources=4):
    print(f"\nD) reutilización de LU (arreglo: misma matriz, {n_sources} rhs)")
    print("f[kHz]  elems   t_spsolve xN[s]   t_splu xN[s]   speedup")
    for f in freqs:
        mesh, _ = build_square(f, side, ppw, material)
        k = material.wavenumber(f)
        a, _ = helmholtz2d.assemble(
            mesh.nodes, mesh.triangles, np.full(mesh.n_elements, k), p_inc=None
        )
        a = a.tocsc()
        rng = np.random.default_rng(0)
        rhs = rng.standard_normal((mesh.n_nodes, n_sources)) + 1j * rng.standard_normal(
            (mesh.n_nodes, n_sources)
        )

        t0 = time.perf_counter()
        for j in range(n_sources):
            spla.spsolve(a, rhs[:, j])
        t_repeat = time.perf_counter() - t0

        t0 = time.perf_counter()
        lu = spla.splu(a)
        for j in range(n_sources):
            lu.solve(rhs[:, j])
        t_lu = time.perf_counter() - t0

        print(f"{f/1e3:6.0f}  {mesh.n_elements:5d}   {t_repeat:14.3f}   {t_lu:12.3f}   "
              f"{t_repeat/t_lu:6.2f}x")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", default="ABCD", help="subconjunto de sweeps: A,B,C,D")
    ap.add_argument("--ppw", type=float, default=10.0)
    ap.add_argument("--lambdas", type=float, default=8.0)
    ap.add_argument("--side", type=float, default=20e-3)
    args = ap.parse_args()
    which = args.only.upper()
    m = materials.WATER
    freqs = [20e3, 40e3, 100e3, 250e3, 500e3, 1e6]

    if "A" in which:
        sweep_a(freqs, args.lambdas, args.ppw, m)
    if "B" in which:
        sweep_b(100e3, 20e-3, m)
    if "C" in which:
        sweep_c([f for f in freqs if f >= 100e3], args.side, args.ppw, m)
    if "D" in which:
        lu_reuse([f for f in freqs if f >= 100e3], args.side, args.ppw, m)


if __name__ == "__main__":
    main()
