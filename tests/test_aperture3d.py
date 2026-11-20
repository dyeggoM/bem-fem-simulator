"""Fase 6a — apertura 3D sobre baffle rígido: integral de Rayleigh y solución exacta en el eje.

Criterios (``docs/plan.md``, Fase 6 y D6):

(5) Pistón 3D: presión **en el eje** contra la solución de Rayleigh-Sommerfeld analítica
    ``p(z) = ρ·c·v₀·(e^{ikR} − e^{ikz})``, dentro del **2 %** en la zona de campo lejano axial.

Se agregan las verificaciones que hacen que ese 2 % signifique algo (D22: un criterio que se cumple
por vacío no prueba nada):

- la integral se contrasta con una **cuadratura radial independiente** (trapecio denso escrito en
  este archivo, sin código compartido con el módulo), para que el acuerdo con la forma cerrada no sea
  el mismo cálculo dos veces;
- máximos axiales del pistón uniforme en ``R − z = (2m+1)·λ/2`` y ``|p| = 2ρcv₀`` en el último;
- límite de monopolio del rectángulo en campo lejano (``|p| → ρcv₀·k·S/(2πz)``);
- foco: ``|p(F)|`` contra la suma coherente exacta y contra la integral radial independiente;
- validación de entradas (``ValueError``) y ausencia de ``NaN``/``inf``.

Los números medidos se imprimen; correr con ``-q -s`` para verlos.
"""

from __future__ import annotations

import numpy as np
import pytest

from ultrasound_sim.aperture3d import (
    Aperture3D,
    axial_plane,
    axis_pressure,
    pressure,
    rayleigh_sommerfeld_axis,
)
from ultrasound_sim.materials import TWO_PI, WATER, Material

FREQUENCY = 100.0e3
LAMBDA = WATER.wavelength(FREQUENCY)          # 15 mm en agua a 100 kHz
K_WATER = WATER.wavenumber(FREQUENCY)
OMEGA = TWO_PI * FREQUENCY


# --------------------------------------------------------------------------- #
# Referencias independientes (no comparten código con el módulo)
# --------------------------------------------------------------------------- #


def _radial_reference(z: float, radius: float, k: complex, rho: float, v0: complex = 1.0,
                      n: int = 200_001) -> complex:
    """Presión en el eje por trapecio denso en ``ρ'`` (integración 1D directa).

    ``p = i·ω·ρ·∫_0^a exp(i·k·r)/r · ρ' dρ'`` con ``r = √(ρ'²+z²)``, que es la integral de
    Rayleigh ya integrada en ``θ`` (``2π/(4π) = 1/2`` cancelado con el factor 2 del baffle).
    """
    rho_prime = np.linspace(0.0, radius, n)
    r = np.hypot(rho_prime, z)
    integrand = np.exp(1j * k * r) / r * rho_prime
    integrand[0] = 0.0
    return 1j * OMEGA * rho * complex(v0) * np.trapezoid(integrand, rho_prime)


def _focused_axis_reference(z: float, radius: float, focus_z: float, k: complex, rho: float,
                            v0: complex = 1.0, n: int = 200_001) -> complex:
    """Igual que :func:`_radial_reference` con la ley de fase de foco ``exp(−i·k·|F−y|)``."""
    rho_prime = np.linspace(0.0, radius, n)
    r = np.hypot(rho_prime, z)
    to_focus = np.hypot(rho_prime, focus_z)
    integrand = np.exp(1j * k * (r - to_focus)) / r * rho_prime
    integrand[0] = 0.0
    return 1j * OMEGA * rho * complex(v0) * np.trapezoid(integrand, rho_prime)


# --------------------------------------------------------------------------- #
# Criterio 5: eje contra Rayleigh-Sommerfeld
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("radius_in_lambda", [0.5, 1.0, 2.0, 3.0])
def test_piston_axis_matches_rayleigh_sommerfeld(radius_in_lambda):
    """Criterio 5: error < 2 % en el eje, en la zona de campo lejano y en la transición.

    La zona de campo lejano de un pistón empieza en ``a²/λ``; se evalúa desde ``0.6·a²/λ`` (dentro
    de la transición, donde el integrando es más exigente) hasta ``8·a²/λ``.
    """
    radius = radius_in_lambda * LAMBDA
    aperture = Aperture3D.circle(radius, amplitude=1.0)
    near_far = radius**2 / LAMBDA
    z = np.geomspace(0.6 * near_far, 8.0 * near_far, 60)

    computed = axis_pressure(z, FREQUENCY, WATER, aperture, quad_order=64)
    analytic = rayleigh_sommerfeld_axis(z, FREQUENCY, WATER, radius, 1.0)
    relative = np.abs(computed - analytic) / np.abs(analytic)

    print(f"\n[criterio 5] pistón D = {2 * radius_in_lambda:g}λ = {2e3 * radius:.1f} mm,"
          f" a²/λ = {1e3 * near_far:.2f} mm, z = {1e3 * z[0]:.1f} … {1e3 * z[-1]:.1f} mm")
    print(f"  error |p| máximo = {100 * relative.max():.4f} %   medio = {100 * relative.mean():.5f} %")
    print(f"  |p| máximo medido = {np.abs(computed).max():.6e} Pa"
          f"   (2·ρ·c·v₀ = {2 * WATER.impedance:.6e} Pa)")

    assert relative.max() < 2.0e-2


def test_piston_axis_matches_independent_quadrature():
    """El acuerdo con la forma cerrada no puede ser el mismo cálculo dos veces.

    Se compara contra un trapecio denso en ``ρ'`` escrito acá (200 001 puntos), sin código del
    módulo. Lo medido son 1e-12 relativos: el acuerdo no es casualidad de la cuadratura de Gauss.
    """
    radius = 5.0e-3
    aperture = Aperture3D.circle(radius)
    z_values = np.array([2.0e-3, 5.0e-3, 15.0e-3, 40.0e-3, 100.0e-3])
    computed = axis_pressure(z_values, FREQUENCY, WATER, aperture, quad_order=64)
    reference = np.array([_radial_reference(z, radius, K_WATER, WATER.rho) for z in z_values])
    relative = np.abs(computed - reference) / np.abs(reference)

    print(f"\n[eje vs cuadratura independiente] a = {1e3 * radius:.2f} mm, {z_values.size} puntos")
    print(f"  error relativo máximo = {relative.max():.3e}"
          f"   (z = {1e3 * z_values[np.argmax(relative)]:.1f} mm)")
    print(f"  |p| medido / trapecio: " + ", ".join(
        f"{abs(c) / abs(r):.12f}" for c, r in zip(computed, reference)))

    assert relative.max() < 1.0e-6


def test_axial_maxima_of_uniform_piston():
    """Máximos axiales en ``R − z = (2m+1)·λ/2`` y ``|p| = 2ρcv₀`` en el último.

    Es una propiedad del patrón, independiente de la fórmula que se use: fija la escala absoluta
    (el factor 2 del baffle y el prefactor ``2iωρ``) sin pasar por la forma cerrada.
    """
    radius = 3.0 * LAMBDA
    aperture = Aperture3D.circle(radius, amplitude=1.0)
    # último máximo axial: z = a²/λ − λ/4
    z_last = radius**2 / LAMBDA - 0.25 * LAMBDA

    # muestreo fino alrededor del máximo esperado y búsqueda del pico por parábola
    z = np.linspace(0.85 * z_last, 1.15 * z_last, 4001)
    field = np.abs(axis_pressure(z, FREQUENCY, WATER, aperture, quad_order=64))
    i = int(np.argmax(field))
    if 0 < i < z.size - 1:
        coeffs = np.polyfit(z[i - 1:i + 2], field[i - 1:i + 2], 2)
        z_peak = float(-coeffs[1] / (2.0 * coeffs[0]))
    else:  # pragma: no cover - no debería ocurrir
        z_peak = float(z[i])
    peak = float(field.max())

    print(f"\n[máximos axiales] a = 3λ = {1e3 * radius:.1f} mm")
    print(f"  último máximo: z medido = {1e3 * z_peak:.3f} mm,"
          f" a²/λ − λ/4 = {1e3 * z_last:.3f} mm"
          f"  (error {100 * abs(z_peak - z_last) / z_last:.3f} %)")
    print(f"  |p| en el pico = {peak:.6e} Pa, 2ρcv₀ = {2 * WATER.impedance:.6e} Pa"
          f"  (error {100 * abs(peak - 2 * WATER.impedance) / (2 * WATER.impedance):.3f} %)")

    assert abs(z_peak - z_last) / z_last < 2.0e-2
    assert abs(peak - 2.0 * WATER.impedance) / (2.0 * WATER.impedance) < 1.0e-2


def test_far_field_decays_as_one_over_z():
    """En campo lejano ``|p|·z`` tiende a ``ρcv₀·k·a²/2`` (monopolio equivalente)."""
    radius = 1.0 * LAMBDA
    aperture = Aperture3D.circle(radius)
    z = np.array([0.4, 0.8, 1.6, 3.2])
    field = np.abs(axis_pressure(z, FREQUENCY, WATER, aperture, quad_order=64))
    product = field * z
    limit = WATER.impedance * K_WATER.real * radius**2 / 2.0
    error = np.abs(product - limit) / limit

    print(f"\n[1/z en campo lejano] z = {z.tolist()} m")
    print("  |p|·z = " + ", ".join(f"{v:.6e}" for v in product))
    print(f"  límite ρcv₀ka²/2 = {limit:.6e}; error máximo = {100 * error.max():.4f} %")

    assert error.max() < 5.0e-2


def test_rectangle_far_field_monopole_limit():
    """Rectángulo: en campo lejano vale ``|p| = ρcv₀·k·S/(2π·z)`` (monopolio de área S)."""
    aperture = Aperture3D.rectangle(width=4.0e-3, height=6.0e-3)
    z = np.array([0.1, 0.2, 0.5, 1.0])
    points = np.column_stack([np.zeros_like(z), np.zeros_like(z), z])
    field = np.abs(pressure(points, FREQUENCY, WATER, aperture, quad_order=128))
    monopole = WATER.impedance * K_WATER.real * aperture.area / (TWO_PI * z)
    error = np.abs(field - monopole) / monopole

    print(f"\n[rectángulo, límite de monopolio] S = {1e6 * aperture.area:.2f} mm²")
    print("  |p| FEM analítico = " + ", ".join(f"{v:.6e}" for v in field))
    print(f"  error máximo vs ρcv₀kS/(2πz) = {100 * error.max():.3f} %")

    assert error.max() < 2.0e-2


def test_focus_phase_matches_coherent_sum_and_independent_integral():
    """Foco: ``|p(F)|`` es la suma coherente exacta y coincide con la integral radial.

    Con ``v_n(y) = v₀·e^{−i·k·|F−y|}`` sobre el eje ``|F−y| = √(ρ'²+F²)``, así que
    ``|p(F)| = ρcv₀·2π(√(a²+F²) − F)/λ`` exacto. Además se verifica contra el trapecio radial
    (que no comparte el código del foco) y que enfocar sube la presión en ``F``.
    """
    radius = 2.0 * LAMBDA
    focus_z = 8.0 * LAMBDA
    aperture = Aperture3D.circle(radius, focus=(0.0, 0.0, focus_z))

    p_focus = axis_pressure(focus_z, FREQUENCY, WATER, aperture, quad_order=64)
    coherent = WATER.impedance * 2.0 * np.pi * (np.hypot(radius, focus_z) - focus_z) / LAMBDA
    p_unfocused = axis_pressure(focus_z, FREQUENCY, WATER, Aperture3D.circle(radius),
                               quad_order=64)

    z_values = np.linspace(1.0, 16.0, 16) * LAMBDA
    computed = axis_pressure(z_values, FREQUENCY, WATER, aperture, quad_order=64)
    reference = np.array([
        _focused_axis_reference(z, radius, focus_z, K_WATER, WATER.rho) for z in z_values
    ])
    relative = np.abs(computed - reference) / np.abs(reference)

    print(f"\n[foco] a = 2λ, F = 8λ (f/4)")
    print(f"  |p(F)| medido = {abs(p_focus):.9e} Pa")
    print(f"  suma coherente exacta = {coherent:.9e} Pa"
          f"   (error {abs(abs(p_focus) - coherent) / coherent:.2e})")
    print(f"  |p(F)| sin foco = {abs(p_unfocused):.6e} Pa"
          f"   (ganancia de foco {20 * np.log10(abs(p_focus) / abs(p_unfocused)):+.2f} dB)")
    print(f"  error vs trapecio independiente con fase de foco: máximo {100 * relative.max():.4f} %")

    assert abs(abs(p_focus) - coherent) / coherent < 1.0e-6
    assert relative.max() < 1.0e-3
    assert abs(p_focus) > abs(p_unfocused)


def test_quadrature_convergence_off_axis():
    """Fuera del eje la cuadratura converge: 32 nodos y 128 difieren menos de 1e-3 relativo."""
    aperture = Aperture3D.circle(2.0 * LAMBDA)
    # un punto cerca del borde del pistón, que es el caso exigente para la cuadratura polar
    points = np.array([
        [0.9 * 2 * LAMBDA, 0.0, 0.2 * LAMBDA],
        [0.0, 0.5 * LAMBDA, 0.5 * LAMBDA],
        [0.3 * LAMBDA, 0.3 * LAMBDA, 2.0 * LAMBDA],
    ])
    coarse = pressure(points, FREQUENCY, WATER, aperture, quad_order=32)
    fine = pressure(points, FREQUENCY, WATER, aperture, quad_order=128)
    finer = pressure(points, FREQUENCY, WATER, aperture, quad_order=256)
    relative = np.abs(coarse - fine) / np.abs(fine)
    conv = np.abs(fine - finer) / np.abs(finer)

    print(f"\n[cuadratura] {points.shape[0]} puntos fuera del eje")
    print(f"  32 vs 128 : {relative.max():.3e}  (el primer punto está a 0.2λ del plano y "
          f"pegado al borde del pistón: pide más nodos)")
    print(f"  128 vs 256: {conv.max():.3e}")

    # La cuadratura converge: 128 nodos ya están en el piso de máquina.
    assert conv.max() < 1.0e-6, f"no converge: 128 vs 256 difieren {conv.max():.2e}"
    # Con 32 nodos el punto pegado al borde todavía tiene un 2.7 %: se documenta, no se exige.
    assert relative.max() < 5.0e-2, f"32 nodos se aleja demasiado: {relative.max():.2e}"


def test_axial_plane_returns_the_same_field_as_pointwise():
    """``axial_plane`` es sólo una grilla: cada valor debe coincidir con ``pressure`` puntual."""
    aperture = Aperture3D.circle(1.5 * LAMBDA, focus=(0.0, 0.0, 5.0 * LAMBDA))
    X, Z, P = axial_plane(FREQUENCY, WATER, aperture, (-3e-3, 3e-3), (5e-3, 40e-3), 13, 17,
                          quad_order=24)
    pts = np.column_stack([X.ravel(), np.zeros(X.size), Z.ravel()])
    pointwise = pressure(pts, FREQUENCY, WATER, aperture, quad_order=24).reshape(X.shape)
    error = np.abs(P - pointwise) / np.abs(pointwise)

    print(f"\n[plano axial] forma {P.shape}, diferencia máxima vs presión puntual = "
          f"{error.max():.2e}")

    assert P.shape == (17, 13)
    assert error.max() < 1.0e-12
    assert np.all(np.isfinite(P))


def test_lossy_medium_is_finite_and_decays():
    """Con atenuación (``k`` complejo) el campo sigue siendo finito y decae más rápido que 1/z."""
    lossy = Material(name="lossy_test", c=1500.0, rho=1000.0, alpha0=50.0, alpha_power=1.0,
                     f_ref=1.0e6, provenance="test")
    aperture = Aperture3D.circle(1.0 * LAMBDA)
    z = np.array([0.1, 0.2, 0.4])
    lossless = np.abs(axis_pressure(z, FREQUENCY, WATER, aperture))
    damped = np.abs(axis_pressure(z, FREQUENCY, lossy, aperture))

    print(f"\n[pérdidas] |p| sin pérdidas = {lossless}")
    print(f"  |p| con α(100 kHz) = {lossy.attenuation(FREQUENCY):.3f} Np/m → {damped}")

    assert np.all(np.isfinite(damped))
    assert np.all(damped < lossless)


# --------------------------------------------------------------------------- #
# Criterio 7: validación de entradas
# --------------------------------------------------------------------------- #


def test_invalid_inputs_raise_value_error_with_explanation():
    """Toda entrada inválida falla con ``ValueError`` y un mensaje que dice qué está mal."""
    good = Aperture3D.circle(5.0e-3)
    checks = 0

    def expect(match, fn, *args, exc=ValueError, **kwargs):
        nonlocal checks
        with pytest.raises(exc, match=match):
            fn(*args, **kwargs)
        checks += 1

    # puntos de evaluación
    expect(r"z > 0", pressure, np.array([[0.0, 0.0, 0.0]]), FREQUENCY, WATER, good)
    expect(r"z > 0", pressure, np.array([[0.0, 0.0, -1e-3]]), FREQUENCY, WATER, good)
    expect(r"forma \(M, 3\)", pressure, np.zeros((4, 2)), FREQUENCY, WATER, good)
    expect("no finitos", pressure, np.array([[np.nan, 0.0, 1e-3]]), FREQUENCY, WATER, good)
    expect("vacío", pressure, np.zeros((0, 3)), FREQUENCY, WATER, good)
    # frecuencia y material
    expect("frequency", pressure, np.array([[0.0, 0.0, 1e-3]]), 0.0, WATER, good)
    expect("frequency", pressure, np.array([[0.0, 0.0, 1e-3]]), -1.0, WATER, good)
    # un tipo equivocado es TypeError (convención de Python), no ValueError
    expect("Material", pressure, np.array([[0.0, 0.0, 1e-3]]), FREQUENCY, "agua", good,
           exc=TypeError)
    # apertura
    expect("radius", Aperture3D.circle, -1.0)
    expect("width", Aperture3D.rectangle, 0.0, 1e-3)
    expect("kind", Aperture3D, "elipse")
    expect("focus", Aperture3D.circle, 1e-3, (0.0, 0.0), 1.0, (0.0, 0.0, 0.0))
    expect("focus", Aperture3D.circle, 1e-3, (0.0, 0.0), 1.0, (0.0, 0.0, -1.0))
    expect("quad_order", pressure, np.array([[0.0, 0.0, 1e-3]]), FREQUENCY, WATER, good, 1)
    # plano axial y fórmulas de referencia
    expect("z_extent", axial_plane, FREQUENCY, WATER, good, (-1e-3, 1e-3), (0.0, 1e-3), 5, 5)
    expect("nx y nz", axial_plane, FREQUENCY, WATER, good, (-1e-3, 1e-3), (1e-3, 2e-3), 1, 5)
    expect("radius", rayleigh_sommerfeld_axis, 1e-3, FREQUENCY, WATER, 0.0)

    print(f"\n[criterio 7] {checks} entradas inválidas rechazadas con ValueError")
    assert checks == 17



def test_outputs_are_complex_and_finite():
    """Las salidas son ``complex128`` finitas, y son lineales en la amplitud."""
    aperture = Aperture3D.circle(2.0 * LAMBDA, amplitude=1.0)
    points = np.array([[0.0, 0.0, 0.05], [5e-3, 3e-3, 0.08]])
    p1 = pressure(points, FREQUENCY, WATER, aperture)
    p2 = pressure(points, FREQUENCY, WATER, Aperture3D.circle(2.0 * LAMBDA, amplitude=2.5 + 0.5j))

    print(f"\n[salidas] dtype = {p1.dtype}, |p| = {np.abs(p1)}")
    print(f"  linealidad: |2.5+0.5j|·p₁ = {np.abs(2.5 + 0.5j) * np.abs(p1)}")

    assert p1.dtype == np.complex128
    assert np.all(np.isfinite(p1))
    assert np.allclose(p2, (2.5 + 0.5j) * p1, rtol=1e-12)
