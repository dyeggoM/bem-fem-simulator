"""Frontera exterior del solver acoplado: ABC-1 (default) contra BGT-2 (``abc_order=2``).

Contexto (``docs/decisions.md`` D12 y D18, ``docs/plan.md`` Fase 5)
-----------------------------------------------------------------
El solver acoplado aplica sobre el campo **total** una ABC de primer orden en la frontera exterior
del fluido. Con un dispersor fuerte y un contorno **curvo** eso tiene un piso de error que no baja
al refinar ni al agrandar el dominio, porque es reflexión de la onda cilíndrica saliente y no
discretización. Medido en la Fase 5a (cápsula casi rígida contra la serie del cilindro rígido):
**12.0 % L2** y **17–21 % en la cara de sombra**.

Este archivo mide lo que cambia al subir la frontera a **BGT-2**
(``∂p_s/∂n = (i·k − 1/(2r) + i/(8kr²))·p_s + (i/2k)·∂²_s p_s``, el mismo operador que
`fem.helmholtz2d.assemble_scattered` con ``abc_order=2``), y deja fijado que el default sigue
siendo la ABC-1 con los números de la Fase 2 intactos.

Cuatro criterios, todos con el número medido impreso (correr con ``-q -s``):

1. Dispersión de la cápsula de `cochlea` (100 kHz, t = 2 mm, hueso 3500/1900/1950) contra la serie
   del cilindro rígido equivalente: el error L2 tiene que bajar de 12 % a < 5 % y el de la cara de
   sombra de 17–21 % a < 7 %. **Medido: 12.00 % → 8.59 % L2 y 16.93/21.87 % → 9.38/9.72 % de
   sombra** en la malla de la Fase 5a. El sesgo de sombra desaparece (queda al nivel de la cara
   iluminada) pero el criterio de 5 %/7 % **no se alcanza**: lo que queda es el piso de
   discretización del modelo acoplado, no la frontera (ver el test 1, el test 2 y el xfail).
2. Reflexión espuria en fluido homogéneo **sin dispersor** (misma malla y mismo contorno): error L2
   contra la onda plana exacta con ABC-1 vs BGT-2, y la mejora en dB. Medido: 20.55 % → 9.00 %
   (+7.2 dB) a ppw=20 y 18.75 % → 4.14 % (+13.1 dB) a ppw=30; la ABC-1 **no** converge y BGT-2 sí.
3. No-regresión: con ``abc_order=1`` los ``|T|`` de la placa de la Fase 2 son idénticos
   (acero 2 mm = 0.29345, hueso 2 mm = 0.79639); con BGT-2 en geometría plana (límite ``r → ∞``)
   coinciden dentro de 0.5 %.
4. Costo: mismo sistema (BGT-2 no agrega gdl), tiempo de ABC-1 vs BGT-2 en el caso de la cápsula.

El test (0) audita el bloque BGT-2 contra la implementación de referencia sin resolver nada: el
bloque que ensambla `coupled2d.abc_boundary_block` tiene que ser **el mismo** que el de
`fem.helmholtz2d.assemble_scattered` sobre la misma circunferencia, y su respuesta a un modo
``exp(i·m·θ)`` tiene que ser el DtN exacto ``k·H_m′(kR)/H_m(kR)``.

Dónde está el piso (medido, para no confundirlo con la frontera)
---------------------------------------------------------------
La ABC sola, en un anillo **acústico** con el mismo ``kR = 24.9`` (sin sólido, sin cavidad):

======================  =========  =========
ppw                     ABC-1      BGT-2
======================  =========  =========
10                       9.99 %     9.78 %
20                       2.84 %     2.52 %
30                       1.78 %     1.13 %
40                       1.56 %     0.64 %
======================  =========  =========

Eso es lo que aporta la frontera: con BGT-2 el error del anillo converge ``O(h²)`` hasta 0.64 % y
no tiene piso. El piso de 8–9 % de la cápsula es del **acoplamiento** (elementos P1, interfaz
curva): se mide con una cáscara de impedancia casada, donde no hay dispersor alguno, y da 9.00 % a
ppw=20 y 4.14 % a ppw=30.

Dos advertencias sobre la configuración
---------------------------------------
- La comparación contra la **serie del cilindro rígido** necesita una cáscara que se mueva muy poco:
  por eso el criterio 1 usa la cáscara casi rígida de `tests/test_cochlea.py` (ρ 200× el agua), que
  es el material con el que se midió el 12.0 %. Con el hueso del plan (3500/1900/1950) el blanco
  **no** es un cilindro rígido (transmite y convierte modo) y el error contra la serie sería el
  desvío del modelo respecto de otra física, no la frontera.
- La cáscara casi rígida sólo está bien con **3 elementos a través de la pared** (``h > t/3``):
  a ``h ≤ t/3`` el modelo salta a 32–36 % con **las dos** ABC y en **todas** las caras, o sea que
  cambia la mecánica discreta de la pared y no la condición de frontera (la cáscara casada, con el
  mismo refinamiento, converge 9.00 → 6.27 → 4.14 %). Por eso el criterio 1 se mide en la malla de
  la Fase 5a y no se refina.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
from scipy.special import hankel1

from ultrasound_sim import meshing
from ultrasound_sim.analytic.helmholtz1d import stack_response
from ultrasound_sim.analytic.scattering import cylinder_scattering
from ultrasound_sim.cochlea import build_cross_section
from ultrasound_sim.fem import helmholtz2d as fem
from ultrasound_sim.fem.assembly import p1_boundary_mass
from ultrasound_sim.fem.coupled2d import (
    abc_boundary_block,
    abc_chain_ends,
    solve_coupled,
)
from ultrasound_sim.materials import WATER, ElasticSolid, Material

FREQUENCY = 1.0e5  # 100 kHz: techo de las fases 0-5 (D8)
THICKNESS = 2.0e-3
LAMBDA = WATER.wavelength(FREQUENCY)

#: Cáscara ``200×`` el agua en densidad: la pared casi no se mueve, así que el blanco se comporta
#: como el cilindro **rígido** de radio ``a + t`` con el que se compara (es el material que usa
#: `tests/test_cochlea.py`; con ρ 20× el propio modelo aportaría ~20 % de error).
NEARLY_RIGID = ElasticSolid(
    name="cascara_casi_rigida",
    rho=200000.0,
    c_p=6000.0,
    c_s=3000.0,
    provenance="diagnóstico del test: impedancia 240× la del agua, no es material de librería",
)
#: Cáscara de impedancia casada (``Z = ρ·c_p = Z_agua``): el dominio es homogéneo a efectos
#: acústicos, así que la solución exacta del fluido es la onda plana incidente y todo lo que
#: aparezca es reflexión espuria. Es la validación de la cadena que ya usa `test_cochlea.py`.
MATCHED = ElasticSolid(
    name="cascara_casi_fluida",
    rho=1000.0,
    c_p=1500.0,
    c_s=100.0,
    provenance="diagnóstico del test: impedancia casada con el agua, no es material de librería",
)
STEEL = ElasticSolid(
    name="placa_benchmark_acero",
    rho=7850.0,
    c_p=5900.0,
    c_s=3200.0,
    provenance="parámetro de benchmark del test, no es material de librería",
)
BONE = ElasticSolid(
    name="placa_benchmark_hueso",
    rho=1950.0,
    c_p=3500.0,
    c_s=1900.0,
    provenance="parámetro de benchmark del test, no es material de librería",
)

#: ``|T|`` de la placa de 2 mm a 100 kHz con la convención de la Fase 2 (capa con ``c_p``),
#: medidos en `tests/test_coupled_fluid_solid.py`: no cambian con ninguna de las dos ABC.
PLATE_T_REFERENCE = {"acero": 0.29345, "hueso": 0.79639}


# ---------------------------------------------------------------------------
# Utilidades
# ---------------------------------------------------------------------------


def _capsule(bone, ppw=20.0, pad=2.5, a_over_t=10.0, min_wall=2):
    """Corte transversal de `cochlea.build_cross_section` con el material pedido."""
    return build_cross_section(
        FREQUENCY, a_over_t * THICKNESS, THICKNESS, ppw=ppw, pad_wavelengths=pad,
        bone=bone, min_elements_across_wall=min_wall,
    )


def _solve(sec, abc_order=1):
    """Resuelve el corte con el solver acoplado y devuelve ``(solución, segundos)``."""
    inc = fem.plane_wave(FREQUENCY, sec.fluid, (1.0, 0.0))
    kwargs = {} if abc_order == 1 else {"abc_center": (0.0, 0.0)}
    t0 = time.perf_counter()
    sol = solve_coupled(
        sec.nodes, sec.triangles, sec.regions, sec.fluid, sec.bone,
        abc_edges=sec.abc_edges, p_inc=inc, frequency=FREQUENCY,
        abc_order=abc_order, **kwargs,
    )
    return sol, time.perf_counter() - t0


def _rel(a, b):
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


def _ring_errors(sec, p, offsets=(0.6, 1.0)):
    """Error L2 del campo FEM contra la serie del cilindro rígido, por anillo y por cara.

    La cara iluminada es ``x < 0`` (la onda viaja en ``+x``), igual que en
    `tests/test_cochlea.py` y en la Fase 4.
    """
    radius = sec.a + sec.thickness
    out = []
    for offset in offsets:
        r = float(sec.radii[np.argmin(np.abs(sec.radii - (radius + offset * LAMBDA)))])
        idx = sec.ring_nodes(r)
        pts = sec.nodes[idx]
        series = cylinder_scattering(FREQUENCY, WATER, radius, pts, boundary="rigid")
        total = series["incident"] + series["scattered"]
        lit = pts[:, 0] < 0.0
        out.append(
            dict(
                r=r, offset=(r - radius) / LAMBDA,
                total=_rel(p[idx], total),
                lit=_rel(p[idx][lit], total[lit]),
                shadow=_rel(p[idx][~lit], total[~lit]),
            )
        )
    return out


def _duct(solid, thickness, ppw_target=40.0, n_lambda=3.0, ny=2):
    """Ducto de la Fase 2: banda sólida de espesor exacto ``thickness`` y ABC en las dos tapas."""
    lam = WATER.wavelength(FREQUENCY)
    n_plate = max(2, int(round(thickness * ppw_target / lam)))
    h = thickness / n_plate
    nx_f = int(round(n_lambda * lam / h))
    x1 = nx_f * h
    x2 = x1 + thickness
    length = x2 + nx_f * h
    height = ny * h
    mesh = meshing.rect_tri_mesh(0.0, length, 0.0, height, 2 * nx_f + n_plate, ny)
    centroids = mesh.centroids()
    regions = np.where((centroids[:, 0] > x1) & (centroids[:, 0] < x2), "solid", "fluid")
    edges = meshing.boundary_edges(mesh.triangles)
    mid = 0.5 * (mesh.nodes[edges[:, 0]] + mesh.nodes[edges[:, 1]])
    tol = 1e-9 * length
    abc = edges[np.isclose(mid[:, 0], 0.0, atol=tol) | np.isclose(mid[:, 0], length, atol=tol)]
    x, y = mesh.nodes[:, 0], mesh.nodes[:, 1]
    walls = np.flatnonzero(
        (x > x1 - tol) & (x < x2 + tol)
        & (np.isclose(y, 0.0, atol=tol) | np.isclose(y, height, atol=tol))
    )
    return dict(mesh=mesh, regions=regions, abc=abc, walls=walls, x1=x1, x2=x2)


def _plate_transmission(solid, thickness=THICKNESS, abc_order=1):
    """``|T|`` complejo del ducto de la Fase 2 (probe sobre la pared ``y = 0``)."""
    geo = _duct(solid, thickness)
    inc = fem.plane_wave(FREQUENCY, WATER, (1.0, 0.0))
    # En geometría plana el límite correcto de BGT-2 es r → ∞: la curvatura 1/(2r) → 0 y queda
    # ik + (i/2k)∂²_s, que es la ABC-2 de una recta. Con centro lejano se ensambla ese límite.
    kwargs = {} if abc_order == 1 else {"abc_center": (1.0e6, 0.0)}
    sol = solve_coupled(
        geo["mesh"].nodes, geo["mesh"].triangles, geo["regions"], WATER, solid,
        abc_edges=geo["abc"], p_inc=inc, dirichlet={"uy": geo["walls"]},
        frequency=FREQUENCY, abc_order=abc_order, **kwargs,
    )
    k = complex(sol.wavenumber)
    x, y = geo["mesh"].nodes[:, 0], geo["mesh"].nodes[:, 1]
    row = np.isclose(y, 0.0, atol=1e-12)
    sel = row & (x > geo["x2"] + 1e-9)
    T = complex((sol.p[sel] * np.exp(-1j * k * (x[sel] - geo["x1"]))).mean())
    ref = row & (x < geo["x1"] - 1e-9)
    R = complex(((sol.p[ref] * np.exp(-1j * k * x[ref]) - 1.0) * np.exp(2j * k * x[ref])).mean())
    return T, R, sol, geo


# ---------------------------------------------------------------------------
# (0) El bloque BGT-2 del solver acoplado es el de helmholtz2d (y el DtN exacto)
# ---------------------------------------------------------------------------


def test_bgt2_block_matches_reference_implementation():
    """``abc_boundary_block`` es, bit por bit, el bloque de `helmholtz2d.assemble_scattered`.

    Es la auditoría del "port": si el operador de la frontera difiere del que ya está validado en
    la Fase 4 (0.94 % contra la serie en un anillo), cualquier número de los tests de abajo mide
    otra cosa. Se compara sobre una circunferencia de ``kR = 3`` y ``ρ = 1`` (el peso ``1/ρ`` se
    chequea aparte por construcción).
    """
    k = WATER.wavenumber(FREQUENCY)
    radius = 1.0 / k.real
    outer_radius = 3.0 * radius
    h = LAMBDA / 20.0
    mesh = meshing.polar_annulus_mesh(
        radius, outer_radius,
        n_radial=max(2, int(np.ceil((outer_radius - radius) / h))),
        n_theta=max(16, int(np.ceil(2.0 * np.pi * outer_radius / h))),
    )
    edges = meshing.boundary_edges(mesh.triangles)
    r_first = np.linalg.norm(mesh.nodes[edges[:, 0]], axis=1)
    outer = edges[np.isclose(r_first, outer_radius, atol=1e-12)]
    kk = np.full(mesh.n_elements, k, dtype=complex)

    coupled = abc_boundary_block(
        mesh.nodes, outer, k, rho=1.0, abc_order=2, abc_center=(0.0, 0.0)
    )
    empty = np.zeros((0, 2), dtype=int)
    a_volume, _ = fem.assemble_scattered(mesh.nodes, mesh.triangles, kk, abc_edges=empty)
    a_with, _ = fem.assemble_scattered(
        mesh.nodes, mesh.triangles, kk, abc_edges=outer, abc_order=2, abc_center=(0.0, 0.0)
    )
    reference = (a_with - a_volume).tocsr()
    diff = abs(coupled - reference)
    rel = float(diff.max() / abs(reference).max())
    print(f"\n[bloque BGT-2] ‖B_acoplado − B_helmholtz2d‖ / ‖B‖ = {rel:.2e} "
          f"(max|B| = {abs(reference).max():.3f}, {outer.shape[0]} aristas)")

    # DtN exacto: la respuesta del bloque a exp(i·m·θ) es k·H_m′(kR)/H_m(kR)
    mass = p1_boundary_mass(mesh.nodes, outer)
    nodes_b = np.unique(outer)
    theta = np.arctan2(mesh.nodes[nodes_b, 1], mesh.nodes[nodes_b, 0])
    z = k * outer_radius
    errors = {}
    for order in (0, 1, 2):
        v = np.zeros(mesh.n_nodes, dtype=complex)
        v[nodes_b] = np.exp(1j * order * theta)
        lam = complex(np.mean(-(coupled @ v)[nodes_b] / (mass @ v)[nodes_b]))
        exact = k * 0.5 * (hankel1(order - 1, z) - hankel1(order + 1, z)) / hankel1(order, z)
        errors[order] = float(abs(lam - exact) / abs(exact))
        print(f"    m={order}: λ_ABC = {lam.real:8.2f}{lam.imag:+8.2f}i   "
              f"λ_exacto = {exact.real:8.2f}{exact.imag:+8.2f}i   error = {errors[order] * 100:5.2f} %")
    assert rel < 1e-12, "el bloque BGT-2 del solver acoplado no es el de helmholtz2d"
    # la condición es una expansión asintótica en 1/(kR): a kR=3 se exige hasta m=1
    assert errors[0] < 0.05
    assert errors[1] < 0.08


def test_abc_chain_ends_and_abc_center_validation():
    """Cadenas cerradas vs abiertas, y qué pasa si ``abc_center`` falta o cae en el contorno."""
    sec = _capsule(NEARLY_RIGID, ppw=10.0, pad=1.5)
    end, nb = abc_chain_ends(sec.abc_edges)
    print(f"\n[cadenas] cápsula: {sec.abc_edges.shape[0]} aristas, extremos libres = {end.size}")
    assert end.size == 0, "el contorno exterior de la cápsula es una cadena cerrada"

    geo = _duct(STEEL, THICKNESS)
    end, nb = abc_chain_ends(geo["abc"])
    print(f"[cadenas] ducto: {geo['abc'].shape[0]} aristas, extremos libres = {end.size} "
          f"(las esquinas donde la tapa se junta con la pared rígida)")
    assert end.size == 4
    pts = geo["mesh"].nodes[end]
    x, y = geo["mesh"].nodes[:, 0], geo["mesh"].nodes[:, 1]
    length = x.max()
    height = y.max()
    at_corners = (np.isclose(pts[:, 0], 0.0) | np.isclose(pts[:, 0], length)) & (
        np.isclose(pts[:, 1], 0.0) | np.isclose(pts[:, 1], height)
    )
    assert np.all(at_corners), "los extremos libres de la ABC tienen que ser las cuatro esquinas"

    inc = fem.plane_wave(FREQUENCY, WATER, (1.0, 0.0))
    nodes, tris, reg = sec.nodes, sec.triangles, sec.regions
    with pytest.raises(ValueError, match="abc_order"):
        solve_coupled(nodes, tris, reg, WATER, sec.bone, abc_edges=sec.abc_edges,
                      p_inc=inc, frequency=FREQUENCY, abc_order=2)
    with pytest.raises(ValueError, match="abc_center"):
        solve_coupled(nodes, tris, reg, WATER, sec.bone, abc_edges=sec.abc_edges,
                      p_inc=inc, frequency=FREQUENCY, abc_order=2,
                      abc_center=(sec.outer_radius, 0.0))  # el centro cae sobre el contorno
    with pytest.raises(ValueError, match="abc_order"):
        solve_coupled(nodes, tris, reg, WATER, sec.bone, abc_edges=sec.abc_edges,
                      p_inc=inc, frequency=FREQUENCY, abc_order=3)
    # una ABC ramificada no es una cadena: el término tangencial no está definido ahí
    with pytest.raises(ValueError, match="cadenas"):
        abc_chain_ends(np.array([[0, 1], [1, 2], [2, 3], [3, 0], [0, 2]]))


# ---------------------------------------------------------------------------
# (3) No-regresión: la placa de la Fase 2 con las dos fronteras
# ---------------------------------------------------------------------------


def test_plate_transmission_is_unchanged():
    """``|T|`` de la placa a 100 kHz: idéntico con ABC-1, y con BGT-2 dentro de 0.5 %.

    En geometría **plana** las dos condiciones son exactas en incidencia normal (una onda plana
    satisface ``∂p/∂n = i·k·p`` y también su versión de segundo orden con curvatura nula), así que
    los tres números no pueden moverse. Es el control de que la frontera nueva no toca el volumen,
    el acoplamiento ni el término incidente.
    """
    print(f"\n[placa] f = {FREQUENCY / 1e3:.0f} kHz, t = {THICKNESS * 1e3:.1f} mm")
    print(f"    {'material':>10} {'ABC':>4} {'|T|':>9} {'|T| ref':>9} {'|R|':>9} "
          f"{'residuo':>9} {'gdl':>6} {'ms':>6}")
    for name, solid in (("acero", STEEL), ("hueso", BONE)):
        layer = Material(name=f"{name}_capa", c=solid.c_p, rho=solid.rho, provenance="capa")
        analytic = stack_response(FREQUENCY, [WATER, layer, WATER], [THICKNESS])
        reference = PLATE_T_REFERENCE[name]
        values = {}
        for order in (1, 2):
            t0 = time.perf_counter()
            T, R, sol, _ = _plate_transmission(solid, abc_order=order)
            dt = (time.perf_counter() - t0) * 1e3
            values[order] = abs(T)
            print(f"    {name:>10} {order:4d} {abs(T):9.5f} {reference:9.5f} {abs(R):9.5f} "
                  f"{sol.residual:9.1e} {sol.size:6d} {dt:6.0f}")
            assert sol.residual < 1e-10
        # ABC-1 reproduce el número histórico de la Fase 2 (no-regresión exacta)
        assert abs(values[1] - reference) < 1e-5, (
            f"ABC-1 en {name}: |T| = {values[1]:.5f}, la Fase 2 medía {reference:.5f}")
        # BGT-2 en geometría plana: mismo resultado dentro de 0.5 %
        dev = abs(values[2] / values[1] - 1.0)
        print(f"    {name:>10} desvío BGT-2 vs ABC-1: {dev * 100:+.3f} %  "
              f"(analítico |T| = {abs(analytic.T):.5f})")
        assert dev < 5e-3
        assert abs(values[2] - abs(analytic.T)) / abs(analytic.T) < 5e-3


# ---------------------------------------------------------------------------
# (1) + (4) Cápsula: error de sombra con BGT-2 y costo
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_capsule_shadow_error_drops_with_bgt2():
    """Criterio 1 + 4: dispersión de la cápsula contra la serie, ABC-1 vs BGT-2, y su costo.

    Malla y material **exactamente** los de la Fase 5a (``ppw = 20``, ``pad = 2.5 λ``, cáscara casi
    rígida), que es de donde salen los números publicados. Medido acá:

    =========================  ==========  ==========  ==========
    métrica                    ABC-1       BGT-2       mejora
    =========================  ==========  ==========  ==========
    L2 (r = R+0.6λ)            12.00 %      8.59 %     +2.9 dB
    L2 (r = R+1.0λ)            13.66 %      8.69 %     +3.9 dB
    iluminada (R+0.6λ)          9.48 %      8.27 %     +1.2 dB
    sombra (R+0.6λ)            16.93 %      9.38 %     +5.1 dB
    sombra (R+1.0λ)            21.87 %      9.72 %     +7.0 dB
    =========================  ==========  ==========  ==========

    Lo que **sí** se verifica: la ABC-1 no se mueve ni un dígito (no-regresión), el sesgo
    iluminada/sombra desaparece (la sombra queda al nivel de la cara iluminada, 9.4 % vs 8.3 %), y
    BGT-2 no agrega grados de libertad.

    Lo que **no** se alcanza en esta malla es el criterio de 5 % / 7 %: el piso que queda (~8.5 %,
    igual en las dos caras) **no es la frontera**. Se mide en
    `test_spurious_reflection_in_homogeneous_fluid`: el mismo modelo acoplado con una cáscara de
    impedancia casada —sin dispersor alguno— da 9.00 % a esta malla y 4.14 % a ppw=30, o sea el
    mismo valor que el residuo de la cápsula y la misma tendencia. Refinar la cápsula casi rígida
    no sirve para bajar de ahí: con 4 elementos a través de la pared (``h ≤ t/3``, que es lo que
    pasa a ppw ≥ 23) el modelo salta a 32–36 % con **las dos** ABC, porque esa discretización de la
    pared grande cambia el modelo mecánico, no la frontera (la cáscara casada, con el mismo
    refinamiento, converge 9.00 → 6.27 → 4.14 %). Por eso el criterio de la fase queda registrado
    como xfail al final, con el número medido y lo que falta.
    """
    print("\n[cápsula casi rígida vs serie del cilindro rígido]")
    sec = _capsule(NEARLY_RIGID, ppw=20.0)
    print(f"  {sec.describe()}")
    results = {}
    for order in (1, 2):
        sol, seconds = _solve(sec, abc_order=order)
        errs = _ring_errors(sec, sol.p)
        results[order] = dict(sol=sol, seconds=seconds, errs=errs)
        print(f"    ABC-{order}: {sol.size} gdl, {seconds:6.1f} s, residuo "
              f"{sol.residual:.1e}")
        for e in errs:
            print(f"      r = R+{e['offset']:.2f}λ: L2 = {e['total'] * 100:5.2f} %  "
                  f"(iluminada {e['lit'] * 100:5.2f} %, sombra {e['shadow'] * 100:5.2f} %)")

    # --- no-regresión: los números de la Fase 5a con ABC-1 no se movieron --------------
    base1 = results[1]["errs"]
    assert abs(base1[0]["total"] - 0.1200) < 2e-3, (
        f"ABC-1 L2 = {base1[0]['total'] * 100:.2f} %, la Fase 5a medía 12.00 %")
    assert 0.160 < base1[0]["shadow"] < 0.180
    assert 0.210 < base1[1]["shadow"] < 0.225
    assert results[1]["sol"].residual < 1e-10

    # --- BGT-2: la sombra deja de estar dominada por la ABC ----------------------------
    base2 = results[2]["errs"]
    for e1, e2 in zip(base1, base2):
        print(f"  R+{e2['offset']:.2f}λ: sombra {e1['shadow'] * 100:5.2f} % → "
              f"{e2['shadow'] * 100:5.2f} %  "
              f"({20 * np.log10(e1['shadow'] / e2['shadow']):+5.1f} dB)")
        assert e2["shadow"] < 0.75 * e1["shadow"], "BGT-2 no bajó el error de sombra"
        assert e2["total"] < 0.90 * e1["total"], "BGT-2 no bajó el error global"
    # el sesgo de sombra desaparece: las dos caras quedan al mismo nivel (el piso)
    for e in base2:
        assert e["shadow"] < e["lit"] + 0.02, (
            f"queda sesgo de sombra: {e['shadow'] * 100:.2f} % vs {e['lit'] * 100:.2f} %")
    assert results[2]["sol"].residual < 1e-10

    # --- costo: mismo sistema, misma malla, mismo ensamblado de volumen ---------------
    cost1, cost2 = results[1], results[2]
    print(f"  costo (mismo sistema de {cost1['sol'].size} gdl): ABC-1 "
          f"{cost1['seconds']:.1f} s, BGT-2 {cost2['seconds']:.1f} s "
          f"(×{cost2['seconds'] / cost1['seconds']:.2f}, ruido de carga: el tiempo lo domina "
          f"`splu`, no el ensamblado)")
    assert cost2["sol"].size == cost1["sol"].size, (
        "BGT-2 no agrega grados de libertad: la condición es de borde")

    # --- el criterio de la fase (5 % L2 / 7 % sombra) en ESTA malla -------------------
    pytest.xfail(
        "criterio 5 %/7 % no alcanzado con BGT-2 en la malla de la Fase 5a: medido "
        f"L2 = {max(e['total'] for e in base2) * 100:.2f} % y sombra = "
        f"{max(e['shadow'] for e in base2) * 100:.2f} % (ABC-1: "
        f"{max(e['total'] for e in base1) * 100:.2f} % / "
        f"{max(e['shadow'] for e in base1) * 100:.2f} %). "
        "Lo que queda no es la frontera: con una cáscara de impedancia casada (sin dispersor) el "
        "mismo modelo acoplado da 9.00 % a ppw=20 y 4.14 % a ppw=30, y la ABC sola (anillo "
        "acústico, kR=24.9) baja a 0.64 % a ppw=40 con BGT-2. Falta bajar el piso de "
        "discretización del acoplamiento (P1, ~O(h²)) o refinar la pared, que a ≥4 elementos "
        "cambia el modelo mecánico de la cáscara casi rígida (32–36 %)."
    )


# ---------------------------------------------------------------------------
# (2) Reflexión espuria en fluido homogéneo sin dispersor
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_spurious_reflection_in_homogeneous_fluid():
    """Criterio 2: sin dispersor, todo lo que el campo se aparte de la onda plana es frontera.

    Cáscara de impedancia casada (``ρ = 1000, c_p = 1500`` → ``Z`` idéntica al agua): a efectos
    acústicos el dominio es homogéneo, la solución exacta del fluido es la onda plana incidente y
    el error L2 contra ella mide **reflexión espuria** de la ABC (más el piso de discretización del
    modelo acoplado, que es el mismo con las dos fronteras).

    Con ABC-1 la fuente espuria es el término incidente que la máscara no aplica en las caras de
    sombra: la condición queda ``∂p_s/∂n = i·k·p_s + i·k·p_i·(d·n − 1)``, y ese residuo no es
    pequeño en incidencia rasante. BGT-2 impone la condición sobre ``p_s`` en toda la frontera, con
    ``∂p_i/∂n − B(p_i)`` exacto para la onda plana, así que el término desaparece.

    Medido (error L2 en el fluido exterior):

    ==========  =========  =========  =========
    ppw         ABC-1      BGT-2      mejora
    ==========  =========  =========  =========
    20          20.55 %     9.00 %    +7.2 dB
    24          19.44 %     6.27 %    +9.8 dB
    30          18.75 %     4.14 %   +13.1 dB
    ==========  =========  =========  =========

    La firma de que la ABC-1 es un **defecto de modelo** y no de discretización: su error no baja al
    refinar (20.55 → 18.75 %), mientras BGT-2 converge (9.00 → 4.14 %, ``O(h²)``) hacia el piso del
    acoplamiento. Por eso la mejora crece con la malla.
    """
    print("\n[fluido homogéneo sin dispersor]")
    errors = {}
    for ppw in (20.0, 30.0):
        sec = _capsule(MATCHED, ppw=ppw)
        r_nodes = np.linalg.norm(sec.nodes, axis=1)
        inc = fem.plane_wave(FREQUENCY, sec.fluid, (1.0, 0.0))
        inc_nodes = inc.value(sec.nodes)
        outside = r_nodes > sec.a + sec.thickness  # fluido exterior: donde debe ser plana
        print(f"  {sec.describe()}")
        for order in (1, 2):
            sol, seconds = _solve(sec, abc_order=order)
            err = _rel(sol.p[outside], inc_nodes[outside])
            errors[(ppw, order)] = err
            print(f"    ABC-{order}: error L2 contra la onda plana = {err * 100:6.2f} %   "
                  f"({seconds:5.1f} s, residuo {sol.residual:.1e})")

    for ppw in (20.0, 30.0):
        gain = 20.0 * np.log10(errors[(ppw, 1)] / errors[(ppw, 2)])
        print(f"  ppw={ppw:2.0f}: ABC-1 {errors[(ppw, 1)] * 100:6.2f} % → BGT-2 "
              f"{errors[(ppw, 2)] * 100:6.2f} %   mejora {gain:+5.1f} dB")

    # con ABC-1 el error es grande y NO converge: es el defecto que motiva el orden 2
    assert errors[(20.0, 1)] > 0.15
    assert errors[(30.0, 1)] > 0.15
    assert errors[(30.0, 1)] > 0.85 * errors[(20.0, 1)], (
        "si la ABC-1 empieza a converger, el diagnóstico de D12 cambió")
    # BGT-2 baja al refinar y llega al piso de discretización del acoplamiento
    assert errors[(30.0, 2)] < 0.6 * errors[(20.0, 2)], "BGT-2 no converge con la malla"
    assert errors[(30.0, 2)] < 0.05
    gain = 20.0 * np.log10(errors[(30.0, 1)] / errors[(30.0, 2)])
    assert gain > 10.0, f"la mejora a ppw=30 es de sólo {gain:.1f} dB"
