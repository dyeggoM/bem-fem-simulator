"""Fase 0: Helmholtz 1D — analítico de referencia y FEM P1 contra analítico.

Criterios de aceptación (docs/plan.md, D6):

1. Interfaz única agua/acero: ``R`` contra ``(Z₂−Z₁)/(Z₂+Z₁)`` con error < 1e-10.
2. Tres capas: error relativo L2 del campo FEM vs analítico < 1 % con ≥ 200 elementos.
3. Caso con atenuación (``alpha0 ≠ None``): el campo transmitido decae como ``exp(−α·d)``.

Los números medidos se imprimen con ``print``: correr con ``-q -s`` para verlos.
"""

from __future__ import annotations

import numpy as np

from ultrasound_sim.analytic.helmholtz1d import pressure, stack_response
from ultrasound_sim.fem.assembly import p1_line_matrices
from ultrasound_sim.fem.helmholtz1d import extract_coefficients, solve
from ultrasound_sim.materials import WATER, Material, fluid_from_impedance

FREQUENCY = 1.0e5  # 100 kHz: λ = 15 mm en agua


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _steel_like() -> Material:
    """Material de prueba construido acá a partir de Z = ρ·c explícito.

    No se agrega a `materials.py`: es un caso de test, no un dato de la librería
    (ver D5 — nada de números plausibles sin procedencia).
    """
    rho = 7850.0  # kg/m³, acero
    c = 5900.0  # m/s, acero
    return fluid_from_impedance("steel_test", impedance=rho * c, rho=rho)


def _l2_relative(x: np.ndarray, p: np.ndarray, exact: np.ndarray) -> float:
    """Error relativo en norma L2 usando la masa P1 (norma exacta del interpolante)."""
    elements = np.column_stack([np.arange(x.size - 1), np.arange(1, x.size)])
    mass, _ = p1_line_matrices(x, elements)
    err = p - exact
    num = np.real(np.conj(err) @ (mass @ err))
    den = np.real(np.conj(exact) @ (mass @ exact))
    return float(np.sqrt(num / den))


# ---------------------------------------------------------------------------
# 1. Interfaz única contra la fórmula de impedancias
# ---------------------------------------------------------------------------


def test_single_interface_matches_impedance_formula():
    steel = _steel_like()
    response = stack_response(FREQUENCY, [WATER, steel], ())

    z1 = WATER.impedance
    z2 = steel.impedance
    r_expected = (z2 - z1) / (z2 + z1)
    t_expected = 2.0 * z2 / (z1 + z2)

    print(f"\n[interfaz única] Z1 = {z1:.6e}  Z2 = {z2:.6e}  Z2/Z1 = {z2 / z1:.3f}")
    print(f"  R analítico (fórmula) = {r_expected:.15f}")
    print(f"  R stack_response      = {response.R.real:.15f}  (|err| = {abs(response.R - r_expected):.2e})")
    print(f"  T analítico (fórmula) = {t_expected:.15f}")
    print(f"  T stack_response      = {response.T.real:.15f}  (|err| = {abs(response.T - t_expected):.2e})")

    assert abs(response.R - r_expected) < 1e-10
    assert abs(response.T - t_expected) < 1e-10


def test_single_interface_field_is_continuous_and_energy_conserving():
    """Continuidad de p y de (1/ρ)·dp/dx en la interfaz, y balance de energía."""
    steel = _steel_like()
    media = [WATER, steel]
    eps = 1e-9

    p_left = pressure(-eps, FREQUENCY, media, ())
    p_right = pressure(+eps, FREQUENCY, media, ())
    # estado (p, (1/ρ)·dp/dx) a cada lado de la interfaz
    from ultrasound_sim.analytic.helmholtz1d import state_at

    w_left = state_at(-eps, FREQUENCY, media, ())[0, 1]
    w_right = state_at(+eps, FREQUENCY, media, ())[0, 1]
    response = stack_response(FREQUENCY, media, ())

    # Balance de intensidad entre medios distintos: |R|² + (Z₁/Z₂)·|T|² = 1
    z1, z2 = WATER.impedance, steel.impedance
    budget = abs(response.R) ** 2 + (z1 / z2) * abs(response.T) ** 2
    print(f"\n[continuidad] salto de p = {abs(p_left - p_right):.2e}")
    print(f"  salto de (1/ρ)dp/dx = {abs(w_left - w_right):.2e}  (relativo {abs(w_left - w_right) / abs(w_right):.1e})")
    print(f"  |R|² + (Z₁/Z₂)·|T|² = {budget:.15f}")

    assert abs(p_left - p_right) < 1e-5
    assert abs(w_left - w_right) < 1e-3 * abs(w_right)
    assert abs(budget - 1.0) < 1e-12


# ---------------------------------------------------------------------------
# 2. Tres capas: FEM P1 vs analítico
# ---------------------------------------------------------------------------


def test_fem_three_layers_matches_analytic():
    # Materiales explícitos de test: contrastes moderados y distintos ρ/c.
    m1 = Material(name="m1_test", c=2500.0, rho=1200.0, provenance="test")
    m2 = Material(name="m2_test", c=1750.0, rho=1300.0, provenance="test")
    m3 = Material(name="m3_test", c=3200.0, rho=1600.0, provenance="test")
    media = [WATER, m1, m2, m3, WATER]
    thicknesses = (4.0e-3, 5.0e-3, 3.0e-3)

    points_per_wavelength = 80
    x, p = solve(
        FREQUENCY, media, thicknesses, points_per_wavelength=points_per_wavelength
    )
    exact = pressure(x, FREQUENCY, media, thicknesses)
    n_elements = x.size - 1
    error = _l2_relative(x, p, exact)

    analytic = stack_response(FREQUENCY, media, thicknesses)
    fem = extract_coefficients(x, p, FREQUENCY, media, thicknesses)

    print(f"\n[3 capas FEM vs analítico] n_elements = {n_elements} (ppw = {points_per_wavelength})")
    print(f"  error relativo L2 = {error:.4e}  ({100 * error:.3f} %)")
    print(f"  R: FEM {fem.R:.10f}  analítico {analytic.R:.10f}  |Δ| = {abs(fem.R - analytic.R):.2e}")
    print(f"  T: FEM {fem.T:.10f}  analítico {analytic.T:.10f}  |Δ| = {abs(fem.T - analytic.T):.2e}")

    assert n_elements >= 200
    assert error < 1.0e-2
    assert abs(fem.R - analytic.R) < 1.0e-2
    assert abs(fem.T - analytic.T) < 1.0e-2


def test_fem_three_layers_converges_with_mesh():
    """El error baja al refinar: confirma que no hay una cancelación accidental."""
    m1 = Material(name="m1_test", c=2500.0, rho=1200.0, provenance="test")
    m2 = Material(name="m2_test", c=1750.0, rho=1300.0, provenance="test")
    m3 = Material(name="m3_test", c=3200.0, rho=1600.0, provenance="test")
    media = [WATER, m1, m2, m3, WATER]
    thicknesses = (4.0e-3, 5.0e-3, 3.0e-3)

    errors = []
    for ppw in (20, 40, 80):
        x, p = solve(FREQUENCY, media, thicknesses, points_per_wavelength=ppw)
        errors.append(_l2_relative(x, p, pressure(x, FREQUENCY, media, thicknesses)))

    print("\n[convergencia 1D] ppw → error L2: " + ", ".join(
        f"{ppw}:{e:.3e}" for ppw, e in zip((20, 40, 80), errors)
    ))
    assert errors[0] > errors[1] > errors[2]


# ---------------------------------------------------------------------------
# 3. Atenuación
# ---------------------------------------------------------------------------


def test_attenuation_transmitted_field_decays_exponentially():
    """Con ``alpha0 ≠ None`` el campo transmitido decae como ``exp(−α·d)``."""
    lossy = Material(
        name="lossy_test",
        c=2000.0,
        rho=1400.0,
        alpha0=800.0,  # Np/m a f_ref
        alpha_power=1.0,
        f_ref=1.0e6,
        provenance="test",
    )
    d = 5.0e-3
    media = [WATER, lossy, lossy]
    thicknesses = (d,)

    alpha_expected = lossy.attenuation(FREQUENCY)
    print(f"\n[atenuación] α(100 kHz) = {alpha_expected:.6f} Np/m "
          f"({lossy.attenuation_db_per_m(FREQUENCY):.3f} dB/m), k = {lossy.wavenumber(FREQUENCY)}")

    # (a) la referencia analítica decae exactamente como exp(−α·Δx)
    x1, x2 = d + 2.0e-3, d + 9.0e-3
    ratio = abs(pressure(x2, FREQUENCY, media, thicknesses)) / abs(
        pressure(x1, FREQUENCY, media, thicknesses)
    )
    expected_ratio = float(np.exp(-alpha_expected * (x2 - x1)))
    print(f"  (a) |p(x2)|/|p(x1)| analítico = {ratio:.15f}   exp(−α·Δx) = {expected_ratio:.15f}")
    assert abs(ratio - expected_ratio) < 1e-12

    # (b) FEM: tasa de decaimiento medida en el semiespacio atenuante
    ppw = 80
    x, p = solve(FREQUENCY, media, thicknesses, points_per_wavelength=ppw)
    in_medium = x > d + 1.0e-3
    alpha_fit = float(
        -np.polyfit(x[in_medium], np.log(np.abs(p[in_medium])), deg=1)[0]
    )
    relative_alpha = abs(alpha_fit - alpha_expected) / alpha_expected
    print(f"  (b) ppw = {ppw} (n = {x.size - 1} elementos)")
    print(f"      α medida por ajuste del campo FEM = {alpha_fit:.6f} Np/m"
          f"  (error relativo {100 * relative_alpha:.4f} %)")
    assert relative_alpha < 5.0e-3

    # (c) T del FEM contra el analítico, y decaimiento medido en dos profundidades
    analytic = stack_response(FREQUENCY, media, thicknesses)
    fem = extract_coefficients(x, p, FREQUENCY, media, thicknesses)
    print(f"  (c) T: FEM {fem.T:.10f}  analítico {analytic.T:.10f}"
          f"  |Δ|/|T| = {abs(fem.T - analytic.T) / abs(analytic.T):.3e}")
    for depth in (1.0e-3, 4.0e-3, 8.0e-3):
        node = int(np.argmin(np.abs(x - (d + depth))))
        measured = abs(p[node])
        predicted = abs(analytic.T) * float(np.exp(-alpha_expected * x[node]))
        print(f"      x−L = {1e3 * depth:4.1f} mm: |p_FEM| = {measured:.6f}"
              f"  |T|·exp(−α·x) = {predicted:.6f}")
        assert abs(measured - predicted) < 2.0e-2 * predicted

    assert abs(fem.T - analytic.T) / abs(analytic.T) < 1.0e-2


def test_field_scales_with_incident_amplitude():
    """El campo y la fuente son lineales en ``p_inc``; ``R`` y ``T`` no dependen de él."""
    media = [WATER, Material(name="m1_test", c=2500.0, rho=1200.0, provenance="test"), WATER]
    thicknesses = (3.0e-3,)

    x1, p1 = solve(FREQUENCY, media, thicknesses, points_per_wavelength=40, p_inc=1.0)
    x3, p3 = solve(FREQUENCY, media, thicknesses, points_per_wavelength=40, p_inc=3.0)
    r1 = extract_coefficients(x1, p1, FREQUENCY, media, thicknesses, p_inc=1.0)
    r3 = extract_coefficients(x3, p3, FREQUENCY, media, thicknesses, p_inc=3.0)

    analytic = stack_response(FREQUENCY, media, thicknesses)
    scale_error = float(np.abs(p3 - 3.0 * p1).max() / np.abs(p3).max())
    print(f"\n[linealidad] max|p(3)−3·p(1)|/max|p(3)| = {scale_error:.3e}")
    print(f"  R(p_inc=1) = {r1.R:.12f}   R(p_inc=3) = {r3.R:.12f}")
    print(f"  |R(1) − R_analítico| = {abs(r1.R - analytic.R):.2e},"
          f" |R(3) − R_analítico| = {abs(r3.R - analytic.R):.2e}")

    assert scale_error < 1e-12
    assert abs(r1.R - r3.R) < 1e-11
    assert abs(pressure(5.0e-3, FREQUENCY, media, thicknesses, p_inc=2.0)
               - 2.0 * pressure(5.0e-3, FREQUENCY, media, thicknesses)) < 1e-12


def test_lossy_wavenumber_follows_repo_convention():
    """k = ω/c + i·α con α ≥ 0: exp(+i·k·x) decae con x (convención D4)."""
    lossy = Material(
        name="lossy_test", c=1500.0, rho=1000.0, alpha0=100.0, f_ref=1.0e6, provenance="test"
    )
    k = lossy.wavenumber(1.0e5)
    print(f"\n[convención] k = {k}  (parte imaginaria {k.imag:+.3f})")
    assert k.imag > 0.0
    assert abs(k.real - 2.0 * np.pi * 1.0e5 / 1500.0) < 1e-9

    alpha = lossy.attenuation(1.0e5)
    media = [lossy, lossy]
    values = pressure(np.array([0.0, 0.01]), 1.0e5, media, ())
    decay = abs(values[1] / values[0])
    expected = float(np.exp(-alpha * 0.01))
    print(f"  |p(1 cm)|/|p(0)| = {decay:.10f}   exp(−α·1 cm) = {expected:.10f}  (α = {alpha} Np/m)")
    assert abs(decay - expected) < 1e-12
