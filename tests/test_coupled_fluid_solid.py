"""Fase 2 — placa elástica 2D en plano deformación acoplada a un fluido acústico.

Validación contra la solución analítica exacta (matriz de transferencia, incidencia normal) y
contra Snell con conversión de modo (incidencia oblicua, ángulos medidos por FFT 2D del
desplazamiento dentro del sólido).

Geometría
---------
Un ducto 2D ``0 <= x <= L``, ``0 <= y <= H`` con una **banda** de elementos sólidos de espesor ``d``
centrada en ``x``; el resto es fluido. ABC de primer orden en ``x = 0`` (entrada, con el término
incidente) y en ``x = L`` (salida); paredes ``y = 0, H`` rígidas para el fluido (Neumann natural) y
**rodillo** (``u_y = 0``) en la parte de la placa.

Por qué el rodillo y no cara libre: con ``u_y = 0`` en las paredes la placa está en deformación
uniaxial (``ε_yy = 0``) y la onda transmitida viaja a ``c_p = √((λ+2μ)/ρ)``, que es exactamente la
capa con la que se construye la referencia analítica
``helmholtz1d.stack_response(f, [agua, (c_p, ρ_s), agua], [d])``. Con caras libres la placa
transmite a la velocidad de placa ``c_p·√((1−2ν)/(1−ν)²)`` (≈8 % menor en acero) y la comparación no
sería contra el mismo problema físico. Para incidencia normal no hay conversión de modo: la
referencia de tres capas con ``c_p`` es exacta.

Incidencia oblicua
------------------
Con paredes rígidas una onda oblicua no satisface la condición de borde, así que la malla se hace
**periódica** en ``y`` con ``H = λ_f/sin θ_i``: entonces ``k_y·H = 2π`` y la onda incidente (y todas
las transmitidas, porque ``k_y`` se conserva) son exactamente periódicas — la condición periódica es
exacta, no una aproximación. Los ángulos se miden por FFT 2D del desplazamiento en la ventana de la
placa (picos en el espacio de k, con ventana de Hann en ``x`` y relleno de ceros).

Los números medidos se imprimen: correr con ``-q -s`` para verlos.
"""

from __future__ import annotations

import numpy as np
import pytest

from ultrasound_sim import meshing
from ultrasound_sim.analytic.helmholtz1d import stack_response
from ultrasound_sim.fem.coupled2d import interface_edges, solve_coupled
from ultrasound_sim.fem.helmholtz2d import plane_wave
from ultrasound_sim.materials import WATER, ElasticSolid, Material

FREQUENCY = 1.0e5  # 100 kHz: techo de las fases 0-5 (docs/decisions.md D8)

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


# ---------------------------------------------------------------------------
# Geometría y medición
# ---------------------------------------------------------------------------


def _layer(solid: ElasticSolid) -> Material:
    """Capa fluida equivalente para el stack analítico: sólo ``c_p`` y ``ρ`` (incidencia normal)."""
    return Material(
        name=f"{solid.name}_capa_equivalente",
        c=solid.c_p,
        rho=solid.rho,
        provenance="capa elástica equivalente en incidencia normal (Z = ρ·c_p)",
    )


def _duct(frequency, solid, thickness, ppw_target=40.0, n_lambda=3.0, ny=2):
    """Ducto con una banda sólida de espesor exactamente ``thickness``."""
    lam_f = WATER.wavelength(frequency)
    n_plate = max(2, int(round(thickness * ppw_target / lam_f)))
    h = thickness / n_plate  # el espesor de la placa es exacto: d = n_plate·h
    nx_f = int(round(n_lambda * lam_f / h))
    x1 = nx_f * h
    x2 = x1 + thickness
    length = x2 + nx_f * h
    nx = 2 * nx_f + n_plate
    height = ny * h
    mesh = meshing.rect_tri_mesh(0.0, length, 0.0, height, nx, ny)

    cen = mesh.centroids()
    regions = np.where((cen[:, 0] > x1) & (cen[:, 0] < x2), "solid", "fluid")

    edges = meshing.boundary_edges(mesh.triangles)
    mid = 0.5 * (mesh.nodes[edges[:, 0]] + mesh.nodes[edges[:, 1]])
    tol = 1e-9 * length
    abc = edges[np.isclose(mid[:, 0], 0.0, atol=tol) | np.isclose(mid[:, 0], length, atol=tol)]

    x, y = mesh.nodes[:, 0], mesh.nodes[:, 1]
    walls = np.flatnonzero(
        (x > x1 - tol) & (x < x2 + tol)
        & (np.isclose(y, 0.0, atol=tol) | np.isclose(y, height, atol=tol))
    )
    return dict(
        mesh=mesh, regions=regions, abc=abc, walls=walls, x1=x1, x2=x2, h=h,
        ppw=lam_f / h, length=length, height=height, n_plate=n_plate,
    )


def _solve_normal(frequency, solid, thickness, **kw):
    geo = _duct(frequency, solid, thickness, **kw)
    inc = plane_wave(frequency, WATER, (1.0, 0.0))
    sol = solve_coupled(
        geo["mesh"].nodes, geo["mesh"].triangles, geo["regions"], WATER, solid,
        abc_edges=geo["abc"], p_inc=inc, dirichlet={"uy": geo["walls"]}, frequency=frequency,
    )
    k = complex(sol.wavenumber)
    x, y = geo["mesh"].nodes[:, 0], geo["mesh"].nodes[:, 1]
    row = np.isclose(y, 0.0, atol=1e-12)

    def probe(mask, kind):
        xs = x[mask]
        if kind == "T":
            vals = sol.p[mask] * np.exp(-1j * k * (xs - geo["x1"]))
        else:
            vals = (sol.p[mask] * np.exp(-1j * k * xs) - 1.0) * np.exp(2j * k * xs)
        return complex(vals.mean()), float(np.max(np.abs(vals - vals.mean())))

    t, spread_t = probe(row & (x > geo["x2"] + 1e-9), "T")
    r, spread_r = probe(row & (x < geo["x1"] - 1e-9), "R")
    geo.update(T=t, R=r, Tspread=spread_t, Rspread=spread_r, sol=sol, k=k)
    return geo


# ---- FFT 2D ---------------------------------------------------------------


def _spectrum(components, h, nfft_factor=8):
    """``Σ_c |FFT2(c·Hann_x)|`` con relleno de ceros; devuelve ``(|F| normalizado, kx, ky)``."""
    ny, nx = components[0].shape
    wx = np.hanning(nx)[None, :]
    nfx, nfy = nfft_factor * nx, nfft_factor * ny
    total = np.zeros((nfy, nfx))
    for c in components:
        total += np.abs(np.fft.fft2(c * wx, s=(nfy, nfx)))
    mag = np.fft.fftshift(total / total.max())
    kx = np.fft.fftshift(2.0 * np.pi * np.fft.fftfreq(nfx, d=h))
    ky = np.fft.fftshift(2.0 * np.pi * np.fft.fftfreq(nfy, d=h))
    return mag, kx, ky


def _refine(mag, kx, ky, iy, ix):
    """Interpolación parabólica en ``log|F|`` de la posición del pico (en la malla rellenada)."""

    def par(values, grid):
        d = np.log(np.maximum(values, 1e-300))
        den = d[0] - 2.0 * d[1] + d[2]
        off = 0.0 if den == 0.0 else 0.5 * (d[0] - d[2]) / den
        off = float(np.clip(off, -0.5, 0.5))
        return grid[1] + off * (grid[2] - grid[1])

    kx_peak = par(mag[iy, ix - 1:ix + 2], kx[ix - 1:ix + 2])
    ky_peak = par(mag[iy - 1:iy + 2, ix], ky[iy - 1:iy + 2])
    return kx_peak, ky_peak


def _peaks(mag, kx, ky, n=8, threshold=0.02):
    """Máximos locales de ``|F|`` con su posición refinada, ordenados por amplitud."""
    out = []
    for iy in range(1, mag.shape[0] - 1):
        for ix in range(1, mag.shape[1] - 1):
            v = mag[iy, ix]
            if v > threshold and v >= mag[iy - 1:iy + 2, ix - 1:ix + 2].max():
                kxp, kyp = _refine(mag, kx, ky, iy, ix)
                out.append((float(v), float(kxp), float(kyp)))
    out.sort(reverse=True)
    return out[:n]


def _periodic_mesh(frequency, solid, thickness, angle, ppw_target=20.0, n_lambda=2.5):
    """Malla periódica en ``y`` con ``H = λ_f/sin(angle)`` (⇒ ``k_y·H = 2π``) y banda sólida."""
    lam_f = WATER.wavelength(frequency)
    height_want = lam_f / np.sin(angle)
    n_plate = max(2, int(round(thickness * ppw_target / lam_f)))
    h = thickness / n_plate
    ny = max(2, int(round(height_want / h)))
    height = ny * h
    nx_f = int(round(n_lambda * lam_f / h))
    x1 = nx_f * h
    x2 = x1 + thickness
    length = x2 + nx_f * h
    nx = 2 * nx_f + n_plate
    mesh = meshing.rect_tri_mesh(0.0, length, 0.0, height, nx, ny)

    cen = mesh.centroids()
    regions = np.where((cen[:, 0] > x1) & (cen[:, 0] < x2), "solid", "fluid")
    edges = meshing.boundary_edges(mesh.triangles)
    mid = 0.5 * (mesh.nodes[edges[:, 0]] + mesh.nodes[edges[:, 1]])
    tol = 1e-9 * length
    abc = edges[np.isclose(mid[:, 0], 0.0, atol=tol) | np.isclose(mid[:, 0], length, atol=tol)]

    nid = lambda i, j: j * (nx + 1) + i
    pairs = np.array([[nid(i, 0), nid(i, ny)] for i in range(nx + 1)], dtype=int)
    return dict(
        mesh=mesh, regions=regions, abc=abc, pairs=pairs, x1=x1, x2=x2, h=h, nx=nx, ny=ny,
        nx_f=nx_f, n_plate=n_plate, height=height, ppw=lam_f / h, length=length,
    )


def _measure_angles(frequency, solid, thickness, angle_deg, ppw_target=20.0):
    """Resuelve la incidencia oblicua periódica y devuelve ``(θ_p, θ_s, picos, geo, sol)``."""
    angle = np.deg2rad(angle_deg)
    geo = _periodic_mesh(frequency, solid, thickness, angle, ppw_target=ppw_target)
    inc = plane_wave(frequency, WATER, (np.cos(angle), np.sin(angle)))
    sol = solve_coupled(
        geo["mesh"].nodes, geo["mesh"].triangles, geo["regions"], WATER, solid,
        abc_edges=geo["abc"], p_inc=inc, periodic_pairs=geo["pairs"], frequency=frequency,
    )
    u = sol.u.reshape(geo["ny"] + 1, geo["nx"] + 1, 2)
    sl = (slice(0, geo["ny"]), slice(geo["nx_f"], geo["nx_f"] + geo["n_plate"] + 1))
    mag, kx, ky = _spectrum([u[sl][:, :, 0], u[sl][:, :, 1]], geo["h"])
    peaks = _peaks(mag, kx, ky)

    k_p = 2.0 * np.pi * frequency / solid.c_p
    k_s = 2.0 * np.pi * frequency / solid.c_s
    forward = [pk for pk in peaks if pk[1] > 0.0 and pk[0] > 0.05]

    def pick(target):
        cand = sorted(forward, key=lambda pk: abs(np.hypot(pk[1], pk[2]) - target))
        if not cand:
            raise AssertionError("no hay picos con k_x > 0 en el espectro del sólido")
        _, kx_p, ky_p = cand[0]
        return np.degrees(np.arctan2(abs(ky_p), abs(kx_p))), np.hypot(kx_p, ky_p)

    theta_p, k_p_meas = pick(k_p)
    theta_s, k_s_meas = pick(k_s)
    return theta_p, theta_s, peaks, geo, sol, (k_p_meas, k_s_meas)


# ---------------------------------------------------------------------------
# (0) Mapa de gdl
# ---------------------------------------------------------------------------


def test_dof_map_layout_and_masks():
    geo = _solve_normal(FREQUENCY, STEEL, 2.0e-3)
    sol, mesh = geo["sol"], geo["mesh"]
    dm = sol.dof_map
    print("\n" + dm.describe())

    # nodos de interfaz: dos líneas verticales de nodos (x = x1 y x = x2) sobre las dos filas
    iface_nodes = dm.interface_nodes
    n_rows = int(round(geo["height"] / geo["h"])) + 1
    assert iface_nodes.size == 2 * n_rows
    assert np.all(dm.u[np.setdiff1d(dm.solid_nodes, iface_nodes), 0] >= 0)
    assert np.all(dm.p[np.setdiff1d(dm.fluid_nodes, iface_nodes)] >= 0)
    # los nodos sólo-sólido no tienen presión y los sólo-fluido no tienen desplazamiento
    only_solid = np.setdiff1d(dm.solid_nodes, iface_nodes)
    only_fluid = np.setdiff1d(dm.fluid_nodes, iface_nodes)
    assert np.all(dm.p[only_solid] == -1)
    assert np.all(dm.u[only_fluid, 0] == -1)
    assert dm.size == dm.n_solid + dm.n_pressure == sol.size
    assert dm.n_solid == 2 * dm.solid_nodes.size
    assert dm.n_pressure == dm.fluid_nodes.size

    # NaN donde no hay gdl, finito donde sí
    assert np.all(np.isnan(sol.p[only_solid]))
    assert np.all(np.isfinite(sol.p[only_fluid]))
    assert np.all(np.isnan(sol.u[only_fluid, 0]))
    assert np.all(np.isfinite(sol.u[only_solid, 0]))
    assert sol.residual < 1e-10, sol.residual
    print(f"residuo relativo del sistema: {sol.residual:.2e}")
    assert isinstance(sol.as_dict()["p"], np.ndarray)


def test_interface_edges_are_tagged_with_the_solid_normal():
    geo = _duct(FREQUENCY, STEEL, 2.0e-3)
    iface, normals = interface_edges(geo["mesh"].nodes, geo["mesh"].triangles, geo["regions"])
    x = geo["mesh"].nodes[iface].mean(axis=1)
    left = iface[np.isclose(x[:, 0], geo["x1"], atol=1e-9)]
    right = iface[np.isclose(x[:, 0], geo["x2"], atol=1e-9)]
    assert left.shape[0] == right.shape[0] == 2  # 2 elementos en y
    nl = normals[np.isclose(x[:, 0], geo["x1"], atol=1e-9)]
    nr = normals[np.isclose(x[:, 0], geo["x2"], atol=1e-9)]
    assert np.allclose(nl, [-1.0, 0.0])  # saliente del sólido en su cara izquierda
    assert np.allclose(nr, [1.0, 0.0])


# ---------------------------------------------------------------------------
# (a) Incidencia normal: |T| contra la fórmula analítica de tres capas
# ---------------------------------------------------------------------------


def test_normal_incidence_transmission_table():
    print("\n(a) incidencia normal, f = 100 kHz")
    header = (f"{'material':12s} {'d[mm]':>6s} {'d/λp':>6s} {'d/λs':>6s} {'ppw':>5s} "
              f"{'|T| FEM':>9s} {'|T| anal':>9s} {'err %':>7s} {'fase(T)°':>9s}")
    print(header)
    print("-" * len(header))
    worst = 0.0
    for solid in (STEEL, BONE):
        for thickness in (2.0e-3, 5.0e-3):
            geo = _solve_normal(FREQUENCY, solid, thickness)
            ref = stack_response(FREQUENCY, [WATER, _layer(solid), WATER], [thickness])
            err = abs(abs(geo["T"]) - abs(ref.T)) / abs(ref.T)
            phase = np.degrees(np.angle(geo["T"] / ref.T))
            worst = max(worst, err)
            print(f"{solid.name:12s} {thickness * 1e3:6.1f} "
                  f"{thickness / (solid.c_p / FREQUENCY):6.3f} "
                  f"{thickness / (solid.c_s / FREQUENCY):6.3f} {geo['ppw']:5.1f} "
                  f"{abs(geo['T']):9.5f} {abs(ref.T):9.5f} {err * 100:7.3f} {phase:9.3f}")
            assert err < 0.02, (
                f"{solid.name}, d={thickness * 1e3:.1f} mm: |T| FEM={abs(geo['T']):.5f} vs "
                f"analítico {abs(ref.T):.5f} → {err * 100:.2f} % > 2 %"
            )
            assert abs(phase) < 10.0
            assert geo["sol"].residual < 1e-10
    print(f"peor error relativo de |T| en la tabla: {worst * 100:.3f} %")


# ---------------------------------------------------------------------------
# (b) Incidencia oblicua: Snell con conversión de modo
# ---------------------------------------------------------------------------


def test_snell_mode_conversion_bone_20deg():
    """Hueso a 20°: se transmiten los dos modos (θ_p = 52.94°, θ_s = 25.67°)."""
    angle_i = 20.0
    thickness = 4.0 * BONE.c_p / FREQUENCY  # 4 λ_p: ventana suficiente para resolver k_x
    theta_p, theta_s, peaks, geo, sol, k_meas = _measure_angles(
        FREQUENCY, BONE, thickness, angle_i, ppw_target=20.0
    )
    sin_i = np.sin(np.deg2rad(angle_i))
    pred_p = np.degrees(np.arcsin(BONE.c_p / WATER.c * sin_i))
    pred_s = np.degrees(np.arcsin(BONE.c_s / WATER.c * sin_i))
    k_y = 2.0 * np.pi / geo["height"]

    print(f"\n(b) hueso, θ_i = {angle_i:.1f}°, d = {thickness * 1e3:.1f} mm "
          f"({thickness / (BONE.c_p / FREQUENCY):.1f} λ_p), ppw = {geo['ppw']:.1f}, "
          f"gdl = {sol.size}, residuo = {sol.residual:.1e}")
    print(f"    k_y = 2π/H = {k_y:.4f} rad/m;  k_p = {k_meas[0]:.3f} (exacto "
          f"{2 * np.pi * FREQUENCY / BONE.c_p:.3f}), k_s = {k_meas[1]:.3f} (exacto "
          f"{2 * np.pi * FREQUENCY / BONE.c_s:.3f})")
    print(f"    picos (|F|/max, kx, ky): "
          + "; ".join(f"({v:.2f}, {kx:.1f}, {ky:.1f})" for v, kx, ky in peaks[:6]))
    print(f"    θ_p: medido {theta_p:7.2f}°  predicho {pred_p:7.2f}°  "
          f"error {theta_p - pred_p:+.2f}°")
    print(f"    θ_s: medido {theta_s:7.2f}°  predicho {pred_s:7.2f}°  "
          f"error {theta_s - pred_s:+.2f}°")

    assert abs(theta_p - pred_p) < 1.0, f"θ_p: {theta_p:.2f}° vs {pred_p:.2f}°"
    assert abs(theta_s - pred_s) < 1.0, f"θ_s: {theta_s:.2f}° vs {pred_s:.2f}°"
    # los dos modos tienen que estar separados en |k| (identificación no ambigua)
    assert k_meas[1] > 1.3 * k_meas[0]


def test_steel_has_no_longitudinal_wave_at_20deg():
    """A 20° el acero está por encima del ángulo crítico del modo P (sin θ_p > 1).

    No es una limitación del código: ``sin θ_p = (c_p/c_f)·sin θ_i = 2.69 > 1`` ⇒ la onda P es
    evanescente en el sólido y no hay ángulo que medir. Se verifica la aritmética de Snell.
    """
    sin_i = np.sin(np.deg2rad(20.0))
    ratio_p = STEEL.c_p / WATER.c * sin_i
    ratio_s = STEEL.c_s / WATER.c * sin_i
    critical = np.degrees(np.arcsin(WATER.c / STEEL.c_p))
    print(f"\n(b') acero a 20°: sin θ_p = {ratio_p:.3f} > 1 (ángulo crítico {critical:.2f}°); "
          f"sin θ_s = {ratio_s:.3f} → θ_s = {np.degrees(np.arcsin(ratio_s)):.2f}°")
    assert ratio_p > 1.0
    assert 0.0 < ratio_s < 1.0


def test_snell_mode_conversion_steel_10deg():
    """Acero a 10°, por debajo del crítico (14.7°): se transmiten P y S."""
    angle_i = 10.0
    thickness = 4.0 * STEEL.c_p / FREQUENCY
    theta_p, theta_s, peaks, geo, sol, k_meas = _measure_angles(
        FREQUENCY, STEEL, thickness, angle_i, ppw_target=16.0
    )
    sin_i = np.sin(np.deg2rad(angle_i))
    pred_p = np.degrees(np.arcsin(STEEL.c_p / WATER.c * sin_i))
    pred_s = np.degrees(np.arcsin(STEEL.c_s / WATER.c * sin_i))
    print(f"\n(b'') acero, θ_i = {angle_i:.1f}°, d = {thickness * 1e3:.1f} mm "
          f"({thickness / (STEEL.c_p / FREQUENCY):.1f} λ_p), ppw = {geo['ppw']:.1f}, "
          f"gdl = {sol.size}, residuo = {sol.residual:.1e}")
    print(f"    |k| medidos: {k_meas[0]:.2f} (P, exacto {2 * np.pi * FREQUENCY / STEEL.c_p:.2f}), "
          f"{k_meas[1]:.2f} (S, exacto {2 * np.pi * FREQUENCY / STEEL.c_s:.2f})")
    print(f"    picos: " + "; ".join(f"({v:.2f}, {kx:.1f}, {ky:.1f})" for v, kx, ky in peaks[:6]))
    print(f"    θ_p: medido {theta_p:7.2f}°  predicho {pred_p:7.2f}°  "
          f"error {theta_p - pred_p:+.2f}°")
    print(f"    θ_s: medido {theta_s:7.2f}°  predicho {pred_s:7.2f}°  "
          f"error {theta_s - pred_s:+.2f}°")
    assert abs(theta_p - pred_p) < 1.0, f"θ_p: {theta_p:.2f}° vs {pred_p:.2f}°"
    assert abs(theta_s - pred_s) < 1.0, f"θ_s: {theta_s:.2f}° vs {pred_s:.2f}°"


# ---------------------------------------------------------------------------
# (c) Resonancia de espesor
# ---------------------------------------------------------------------------


def test_thickness_resonance_sweep():
    print("\n(c) resonancia de espesor (|T| máximo en f₁ = c_p/2d)")
    for solid, thickness in ((STEEL, 2.0e-3), (STEEL, 5.0e-3),
                             (BONE, 2.0e-3), (BONE, 5.0e-3)):
        f1 = solid.c_p / (2.0 * thickness)
        freqs = np.linspace(0.85, 1.15, 41) * f1
        trans = np.array([abs(_solve_normal(f, solid, thickness, ppw_target=30.0)["T"])
                          for f in freqs])
        i = int(np.argmax(trans))
        assert 0 < i < freqs.size - 1, "el pico cayó en el borde del barrido"
        y = np.log(trans[i - 1:i + 2])
        off = 0.5 * (y[0] - y[2]) / (y[0] - 2.0 * y[1] + y[2])
        f_peak = freqs[i] + off * (freqs[1] - freqs[0])
        err = (f_peak - f1) / f1
        print(f"    {solid.name:20s} d = {thickness * 1e3:.1f} mm | f₁ = c_p/(2d) = "
              f"{f1 / 1e6:.4f} MHz | pico medido = {f_peak / 1e6:.4f} MHz | "
              f"error = {err * 100:+.3f} % | |T|(pico) = {trans[i]:.4f}")
        assert abs(err) < 0.02, f"{solid.name}: pico {f_peak:.1f} Hz vs {f1:.1f} Hz"
        assert trans[i] > 0.98  # a la resonancia de media onda la placa es transparente


# ---------------------------------------------------------------------------
# (d) Energía
# ---------------------------------------------------------------------------


def test_energy_balance():
    """|R|² + |T|² = 1 (mismo fluido a los dos lados ⇒ sin factor de impedancia)."""
    print("\n(d) balance de energía (incidencia normal)")
    for solid in (STEEL, BONE):
        geo = _solve_normal(FREQUENCY, solid, 2.0e-3)
        total = abs(geo["R"]) ** 2 + abs(geo["T"]) ** 2
        print(f"    {solid.name:12s} |R| = {abs(geo['R']):.6f}  |T| = {abs(geo['T']):.6f}  "
              f"|R|²+|T|² = {total:.6f}  (desvío {total - 1:+.2e})")
        assert abs(total - 1.0) < 1e-2, f"{solid.name}: |R|²+|T|² = {total:.6f}"


# ---------------------------------------------------------------------------
# (e) Condición de interfaz (fija la convención de signos)
# ---------------------------------------------------------------------------


def test_interface_state_matches_analytic():
    """En la cara izquierda, ``(p, u_x)`` del FEM contra el estado analítico del stack.

    En la capa equivalente (``c_p``, ``ρ_s``) el segundo componente del estado es
    ``w = (1/ρ)·dp/dx``; la continuidad de velocidad normal da ``w = ω²·u_x``, así que
    ``w/ω²`` **es** el desplazamiento normal del sólido. Comparar las dos componentes contra el
    estado analítico en ``x = 0`` (la cara izquierda de la placa) fija la convención de signos del
    acoplamiento: un signo global invertido en ``C`` dejaría ``p`` idéntico y ``u`` con un desfase
    π, que este test detecta. Además ``σ_xx = −p`` en la cara, que es la otra condición de
    interfaz.
    """
    from ultrasound_sim.analytic.helmholtz1d import state_at

    geo = _solve_normal(FREQUENCY, STEEL, 2.0e-3, ppw_target=60.0)
    sol, mesh = geo["sol"], geo["mesh"]
    omega = 2.0 * np.pi * FREQUENCY
    x = mesh.nodes[:, 0]
    row = np.isclose(mesh.nodes[:, 1], 0.0, atol=1e-12)
    face = np.flatnonzero(row & np.isclose(x, geo["x1"], atol=1e-9))[0]

    state = state_at(np.array([0.0]), FREQUENCY, [WATER, _layer(STEEL), WATER], [2.0e-3])[0]
    p_ref = state[0]
    u_ref = state[1] / omega**2

    p_fem = sol.p[face]
    u_fem = sol.u[face, 0]
    rel_p = abs(p_fem - p_ref) / abs(p_ref)
    rel_u = abs(u_fem - u_ref) / abs(u_ref)
    phase_u = np.degrees(np.angle(u_fem / u_ref))
    print(f"\n(e) cara izquierda: p FEM = {p_fem:.6f} vs analítico {p_ref:.6f} (rel {rel_p:.2e})")
    print(f"    u_x FEM = {u_fem:.4e} vs w/ω² = {u_ref:.4e} "
          f"(rel {rel_u:.2e}, fase {phase_u:+.3f}°)")
    assert rel_p < 5e-2
    assert rel_u < 5e-2, f"u_x difiere del analítico en {rel_u * 100:.2f} %"
    assert abs(phase_u) < 2.0  # un signo invertido daría 180°
    # σ_xx = (λ+2μ)·du_x/dx = −p  ⇒  el gradiente del desplazamiento es negativo con p > 0
    assert (sol.p[face].real > 0) and (np.sign((sol.u[face + 1, 0] - u_fem).real) < 0
                                       or np.sign((u_fem - sol.u[face - 1, 0]).real) < 0)


# ---------------------------------------------------------------------------
# Validaciones
# ---------------------------------------------------------------------------


def test_validations():
    geo = _duct(FREQUENCY, STEEL, 2.0e-3)
    nodes, tris, reg = geo["mesh"].nodes, geo["mesh"].triangles, geo["regions"]
    inc = plane_wave(FREQUENCY, WATER, (1.0, 0.0))

    with pytest.raises(ValueError, match="regions"):
        solve_coupled(nodes, tris, np.where(reg == "solid", "hueso", "agua"), WATER, STEEL,
                      frequency=FREQUENCY)
    with pytest.raises(ValueError, match="un valor por elemento"):
        solve_coupled(nodes, tris, reg[:-1], WATER, STEEL, frequency=FREQUENCY)
    with pytest.raises(TypeError, match="ElasticSolid"):
        solve_coupled(nodes, tris, reg, WATER, WATER, frequency=FREQUENCY)
    with pytest.raises(TypeError, match="Material"):
        solve_coupled(nodes, tris, reg, STEEL, STEEL, frequency=FREQUENCY)
    with pytest.raises(ValueError, match="frequency"):
        solve_coupled(nodes, tris, reg, WATER, STEEL)
    with pytest.raises(ValueError, match="abc_order"):
        solve_coupled(nodes, tris, reg, WATER, STEEL, abc_order=2, frequency=FREQUENCY)
    with pytest.raises(ValueError, match="no corresponde"):
        solve_coupled(nodes, tris, reg, WATER, STEEL, frequency=2.0 * FREQUENCY, p_inc=inc)
    with pytest.raises(ValueError, match="solid"):
        # una arista de la frontera de la placa (pared superior del sólido) no puede ser ABC
        edges = meshing.boundary_edges(tris)
        mid = 0.5 * (nodes[edges[:, 0]] + nodes[edges[:, 1]])
        on_plate_top = edges[
            np.isclose(mid[:, 1], geo["height"], atol=1e-12)
            & (mid[:, 0] > geo["x1"] + 1e-9) & (mid[:, 0] < geo["x2"] - 1e-9)
        ]
        solve_coupled(nodes, tris, reg, WATER, STEEL, abc_edges=on_plate_top, frequency=FREQUENCY)
    with pytest.raises(ValueError, match="no es de frontera"):
        solve_coupled(nodes, tris, reg, WATER, STEEL, abc_edges=np.array([[0, geo["mesh"].triangles.shape[0]]]),
                      frequency=FREQUENCY)
    with pytest.raises(ValueError, match="dirichlet"):
        solve_coupled(nodes, tris, reg, WATER, STEEL, frequency=FREQUENCY,
                      dirichlet={"uy": [0]})
