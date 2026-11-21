"""El acople de N regiones contra la teoría de tres medios (U9).

Qué se valida
-------------
`solve_multicapa` acopla varios materiales en una sola malla conforme. Lo que hay que demostrar no es que
"corra" sino que el campo transmitido sea el correcto:

1. **Caso transparente**: si la capa tiene las propiedades del medio, la transmisión es exactamente 1.
   Es el control que distingue un acople que no hace nada de uno que acopla bien.
2. **Una placa con contraste** contra el coeficiente de transmisión analítico de tres medios
   (Z1 = Z3 = agua):

       |t| = 1 / sqrt( cos²φ + ¼·(m + 1/m)²·sin²φ ) ,   φ = k₂·d ,  m = Z₂/Z₁

3. **Convergencia**: el error tiene que BAJAR al refinar la malla. Si no baja, el error no es de
   discretización y el número de arriba no significaría nada.

El error que se encontró armando esto (y que queda fijado abajo)
-----------------------------------------------------------------
La primera versión aplicaba la ABC a **toda** la frontera de fluido, incluidas las paredes laterales. Ahí
la onda plana tiene ``∂p/∂n = 0``, pero la ABC impone ``∂p/∂n = i·k·p``: metía una condición espuria que
reflejaba, y el caso base daba |p| entre 0.82 y 1.37 que **no mejoraba al refinar** (λ/12: 0.819; λ/20:
0.848). Con ABC sólo en entrada y salida el mismo caso da |p| = 1.003 (λ/10) y 1.002 (λ/16). Por eso los
tests de acá usan ABC sólo en ``x = 0`` y ``x = L``.
"""

from __future__ import annotations

import numpy as np
import pytest

from ultrasound_sim.fem import helmholtz2d, multicapa
from ultrasound_sim.materials import WATER, Material
from ultrasound_sim.meshing import boundary_edges, edge_midpoints, grid_tri_mesh

F = 250e3
LARGO = 0.06


def _placa_hueso() -> Material:
    """Hueso como **fluido** equivalente: sirve para comparar con la fórmula de tres medios (sin cizalla).

    Propiedades de `docs/materiales-cadena.md` (Pichardo 2011): ρ 1896,5 kg/m³, c_L 2471 m/s a ~0.8 MHz.
    """
    return Material(name="hueso_fluido", rho=1896.5, c=2471.0, alpha0=0.0)


def _correr(espesor: float, *, h_rel: float = 20.0, material_placa: Material | None = None):
    """Ducto agua | placa | agua con ABC sólo en entrada y salida."""
    lam = WATER.wavelength(F)
    h = lam / h_rel
    xs = np.linspace(0.0, LARGO, int(round(LARGO / h)) + 1)
    ys = np.linspace(-lam, lam, int(round(2 * lam / h)) + 1)
    malla = grid_tri_mesh(xs, ys)
    nodos, tris = malla.nodes, malla.triangles

    cx = nodos[tris][:, :, 0].mean(axis=1)
    placa = material_placa if material_placa is not None else _placa_hueso()
    medio = np.full(cx.shape, "agua", dtype=object)
    dentro = (cx >= 0.02) & (cx < 0.02 + espesor)
    medio[dentro] = "placa"

    ar = boundary_edges(tris)
    mid = edge_midpoints(nodos, ar)
    abc = ar[(np.abs(mid[:, 0]) < 1e-12) | (np.abs(mid[:, 0] - LARGO) < 1e-12)]

    onda = helmholtz2d.plane_wave(F, WATER, (1, 0))
    r = multicapa.solve_multicapa(
        nodos, tris, medio, {"agua": WATER, "placa": placa}, {}, F, abc_edges=abc, onda=onda
    )
    eje = (np.abs(nodos[:, 1]) < 1e-12) & (nodos[:, 0] > 0.02 + espesor + 0.008)
    return float(np.abs(r.p_todos[eje]).mean())


def _t_analitico(espesor: float) -> float:
    placa = _placa_hueso()
    z1 = WATER.rho * WATER.c
    z2 = placa.rho * placa.c
    phi = 2 * np.pi * F / placa.c * espesor
    cociente = z2 / z1
    return 1.0 / np.sqrt(np.cos(phi) ** 2 + 0.25 * (cociente + 1.0 / cociente) ** 2 * np.sin(phi) ** 2)


def test_el_caso_transparente_da_transmision_uno():
    """Capa con las propiedades del agua: nada que transmitir ni reflejar, |p| = 1."""
    sin_capa = _correr(0.0)
    con_capa = _correr(0.0025, material_placa=WATER)
    print(f"\n  sin capa |p| = {sin_capa:.5f} · capa del mismo material |p| = {con_capa:.5f}")
    assert abs(sin_capa - 1.0) < 0.01, "la onda plana sola no llega a 1: el dominio tiene un problema"
    assert abs(con_capa - 1.0) < 0.01, "una capa transparente tendría que transmitir sin atenuar"


def test_una_placa_con_contraste_reproduce_la_formula_de_tres_medios():
    espesor = 2.5e-3
    fem = _correr(espesor, h_rel=20.0)
    teo = _t_analitico(espesor)
    error = abs(fem - teo) / teo
    print(f"\n  placa de {espesor * 1e3:.1f} mm: |t| FEM = {fem:.4f} · analítico = {teo:.4f} "
          f"· error {100 * error:.2f} %")
    assert error < 0.02


def test_el_error_baja_al_refinar_la_malla():
    """Si el error no baja, no es discretización: sería un error de formulación disfrazado."""
    espesor = 2.5e-3
    teo = _t_analitico(espesor)
    errores = []
    for h_rel in (10.0, 20.0):
        fem = _correr(espesor, h_rel=h_rel)
        errores.append(abs(fem - teo) / teo)
    print(f"\n  error λ/10: {100 * errores[0]:.2f} % · λ/20: {100 * errores[1]:.2f} %")
    assert errores[1] < errores[0], "refinar la malla no mejoró el resultado"


def _t_matriz_de_transferencia(capas) -> float:
    """|t| exacto para N capas por matriz de transferencia (incidencia normal, sin pérdidas).

    Cada capa aporta ``M_j = [[cos φ, i·Z·sin φ], [i·sin φ/Z, cos φ]]`` que relaciona ``[p; u]`` entre sus
    dos caras. Con ``[p0; u0] = [1+R; (1-R)/Z1]`` en la entrada y ``p_N/u_N = Z_N`` en la salida:

        R = (β - α)/(α + β),  α = M11 - Z_N·M21,  β = (Z_N·M22 - M12)/Z1

    y la presión **transmitida** es ``p_N = M11·p0 + M12·u0``, NO ``1+R``.

    Esto se escribió mal la primera vez y el error vale la pena porque es sutil: ``1+R`` es la presión
    total en la **primera** interfaz, dentro del medio 1 (incidente + reflejada). Con una capa de cuarto
    de onda ``1+R`` puede pasar de 1 —es el transformador de impedancia, la presión en la cara de entrada
    sube— y tomarlo como transmitida daba |t| = 1.814 donde el FEM y la fórmula de tres medios dicen
    0.581. Un ``|t| > 1`` sin pérdidas es imposible: eso era la señal.
    """
    z1 = WATER.rho * WATER.c
    z_n = z1                                   # la cadena sale a agua en los dos lados
    m = np.eye(2, dtype=complex)
    for material, espesor in capas:
        z = material.rho * material.c
        phi = 2 * np.pi * F / material.c * espesor
        capa = np.array([[np.cos(phi), 1j * z * np.sin(phi)],
                         [1j * np.sin(phi) / z, np.cos(phi)]], dtype=complex)
        m = m @ capa
    alfa = m[0, 0] - z_n * m[1, 0]
    beta = (z_n * m[1, 1] - m[0, 1]) / z1
    reflejado = (beta - alfa) / (alfa + beta)
    p0 = 1.0 + reflejado
    u0 = (1.0 - reflejado) / z1
    return float(abs(m[0, 0] * p0 + m[0, 1] * u0))


def test_una_cadena_de_tres_materiales_coincide_con_la_matriz_de_transferencia():
    """piel 2 mm + hueso 5 mm contra el |t| exacto de dos capas.

    Ojo con la intuición: a 250 kHz el hueso tiene λ ≈ 9.9 mm, así que 5 mm es una **resonancia de media
    onda** y la placa transmite casi todo. Un test que exigiera "atenúa" pasaría por la razón equivocada
    el día que el código estuviera roto, por eso se compara contra el valor exacto.
    """
    piel = Material(name="piel_fluido", rho=1116.0, c=1537.0, alpha0=0.0)
    lam = WATER.wavelength(F)
    h = lam / 18.0
    xs = np.linspace(0.0, LARGO, int(round(LARGO / h)) + 1)
    ys = np.linspace(-lam, lam, int(round(2 * lam / h)) + 1)
    malla = grid_tri_mesh(xs, ys)
    nodos, tris = malla.nodes, malla.triangles
    cx = nodos[tris][:, :, 0].mean(axis=1)
    medio = np.full(cx.shape, "agua", dtype=object)
    medio[(cx >= 0.02) & (cx < 0.022)] = "piel"
    medio[(cx >= 0.022) & (cx < 0.027)] = "hueso"
    ar = boundary_edges(tris)
    mid = edge_midpoints(nodos, ar)
    abc = ar[(np.abs(mid[:, 0]) < 1e-12) | (np.abs(mid[:, 0] - LARGO) < 1e-12)]
    r = multicapa.solve_multicapa(
        nodos, tris, medio,
        {"agua": WATER, "piel": piel, "hueso": _placa_hueso()}, {}, F, abc_edges=abc,
        onda=helmholtz2d.plane_wave(F, WATER, (1, 0)),
    )
    eje = (np.abs(nodos[:, 1]) < 1e-12) & (nodos[:, 0] > 0.035)
    transmitido = float(np.abs(r.p_todos[eje]).mean())
    exacto = _t_matriz_de_transferencia([(piel, 2e-3), (_placa_hueso(), 5e-3)])
    error = abs(transmitido - exacto) / exacto
    print(f"\n  piel 2 mm + hueso 5 mm: |p| FEM = {transmitido:.4f} · matriz de transferencia "
          f"= {exacto:.4f} · error {100 * error:.2f} %")
    assert np.all(np.isfinite(r.p_todos))
    assert error < 0.03
