"""Contraste BEM propio vs serie exacta para un cilindro rígido (Fase 3).

Criterio de aceptación (docs/plan.md, Fase 3): error relativo L2 del campo **dispersado**,
evaluado sobre un círculo ``r = 2a``, menor al 2 % contra la serie de Bessel, para ``ka = 1, 2, 5``
y también con atenuación (``k`` complejo).

Además se verifican las dos piezas que el caso rígido no ejercita por sí solo:

* el auto-término ``S_ii`` (única integral singular de la formulación) contra una cuadratura
  numérica fina independiente;
* un campo fabricado con fuente dentro del cuerpo, que valida la capa simple ``S`` y la ecuación
  integral con dato de Neumann no nulo.
"""

import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss
from scipy.special import hankel1

from ultrasound_sim.analytic.scattering import cylinder_scattering
from ultrasound_sim.bem.helmholtz2d import (
    _layer_matrices,
    _self_panel_single_layer,
    circle_boundary,
    panels_for_frequency,
    solve_exterior_neumann,
)
from ultrasound_sim.materials import WATER, Material

RADIUS = 0.01
KA_VALUES = (1.0, 2.0, 5.0)

#: Tolerancia de aceptación de la fase 3.
TOLERANCE = 0.02


def frequency_for(ka: float, material: Material = WATER) -> float:
    return ka * material.c / (2.0 * np.pi * RADIUS)


def circle_points(radius: float, n: int = 64) -> np.ndarray:
    theta = 2.0 * np.pi * np.arange(n) / n
    return radius * np.column_stack([np.cos(theta), np.sin(theta)])


def relative_error(a, b) -> float:
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


@pytest.mark.parametrize("ka", KA_VALUES)
def test_rigid_cylinder_matches_series(ka):
    """Campo dispersado BEM vs serie exacta sobre ``r = 2a`` (error L2 < 2 %)."""
    frequency = frequency_for(ka)
    k = WATER.wavenumber(frequency)
    points = circle_points(2.0 * RADIUS)

    solution = solve_exterior_neumann(points, frequency, WATER, RADIUS)
    reference = cylinder_scattering(frequency, WATER, RADIUS, points)["scattered"]

    error = relative_error(solution.p_scattered, reference)
    error_total = relative_error(solution.p_total,
                                 cylinder_scattering(frequency, WATER, RADIUS,
                                                     points)["total"])
    wavelength = 2.0 * np.pi / k.real
    per_wavelength = wavelength / solution.boundary.panel_length

    print(f"\n  ka={ka:<4g} paneles={solution.n_panels:<4d} incógnitas={solution.n_panels:<4d} "
          f"paneles/λ={per_wavelength:.1f} cond={solution.condition_number:.2f} "
          f"error_dispersado={error:.3e} error_total={error_total:.3e}")
    assert error < TOLERANCE


def test_attenuated_medium_matches_series():
    """Con ``k`` complejo (α0 ≠ None) el BEM y la serie siguen coincidiendo dentro del 2 %."""
    ka = 2.0
    # α0 = 400 Np/m a 1 MHz con n = 1 ⇒ α ≈ 19 Np/m en ka=2 (Im k / Re k ≈ 10 %).
    # Es un parámetro de test deliberadamente fuerte, no un valor físico de tejido.
    lossy = Material(
        name="water_lossy_test",
        c=WATER.c,
        rho=WATER.rho,
        alpha0=400.0,
        alpha_power=1.0,
        f_ref=1.0e6,
        provenance="parámetro de test (ejercita k complejo), no es un dato físico",
    )
    frequency = frequency_for(ka, lossy)
    k = lossy.wavenumber(frequency)
    points = circle_points(2.0 * RADIUS)

    solution = solve_exterior_neumann(points, frequency, lossy, RADIUS)
    reference = cylinder_scattering(frequency, lossy, RADIUS, points)["scattered"]
    error = relative_error(solution.p_scattered, reference)

    print(f"\n  ka={ka:<4g} k={k:.3f} 1/m (Im/Re={k.imag / k.real:.1%}) "
          f"paneles={solution.n_panels} cond={solution.condition_number:.2f} "
          f"error_dispersado={error:.3e}")
    assert k.imag > 1.0, "el test no estaría ejercitando k complejo"
    assert error < TOLERANCE


def test_self_panel_single_layer_matches_numerical_integral():
    """``S_ii`` (extracción de la singularidad logarítmica) vs cuadratura fina independiente.

    La referencia numérica usa la sustitución ``u = v·|v|``, que quita la singularidad logarítmica
    de ``H₀`` en ``s = 0`` y converge con muchos nodos de Gauss-Legendre.
    """
    for ka in (0.5, 2.0, 5.0, 8.0):
        k = WATER.wavenumber(frequency_for(ka))
        boundary = circle_boundary(RADIUS, 64)
        length = boundary.panel_length

        order = 3000
        v, w = leggauss(order)
        u = v * np.abs(v)
        jacobian = 2.0 * np.abs(v)
        reference = (length / 2.0) * np.sum(
            w * jacobian * 0.25j * hankel1(0, k * (length / 2.0) * np.abs(u))
        )
        value = _self_panel_single_layer(k, length)
        error = abs(value - reference) / abs(reference)

        print(f"\n  ka={ka:<5g} L={length:.4e} m  S_ii={value:.6e}  "
              f"referencia={reference:.6e}  error_rel={error:.2e}")
        assert error < 1e-5


def test_manufactured_point_source_validates_single_layer():
    """Campo fabricado (fuente dentro del cuerpo): valida ``S`` y la ecuación con ``q ≠ 0``.

    ``p = (i/4)H₀⁽¹⁾(k|x − x₀|)`` con ``x₀`` dentro del cilindro satisface Helmholtz homogéneo y
    radiación en el fluido, así que ``(½I − D)p = −S·q`` debe reproducirlo sin usar la onda plana.
    El caso rígido (``q = 0``) nunca toca ``S``; este test sí.
    """
    ka = 2.0
    frequency = frequency_for(ka)
    k = WATER.wavenumber(frequency)
    source = np.array([0.3 * RADIUS, 0.2 * RADIUS])

    def field(x):
        r = np.hypot(*(x - source).T)
        return 0.25j * hankel1(0, k * r)

    def normal_derivative(y, normals):
        d = y - source
        r = np.hypot(d[:, 0], d[:, 1])
        return 0.25j * (-k) * hankel1(1, k * r) * np.sum(d * normals, axis=1) / r

    points = circle_points(2.0 * RADIUS)
    errors = []
    for n_panels in (128, 256):
        boundary = circle_boundary(RADIUS, n_panels)
        p_boundary = field(boundary.midpoints)
        q_boundary = normal_derivative(boundary.midpoints, boundary.normals)
        S, D = _layer_matrices(boundary.midpoints, boundary, k, handle_self_panel=True)

        solved = np.linalg.solve(0.5 * np.eye(n_panels) - D, -(S @ q_boundary))
        S_eval, D_eval = _layer_matrices(points, boundary, k)
        evaluated = D_eval @ solved - S_eval @ q_boundary
        errors.append(relative_error(evaluated, field(points)))
        print(f"\n  fabricado paneles={n_panels}: error en la frontera="
              f"{relative_error(solved, p_boundary):.3e} error en r=2a={errors[-1]:.3e}")

    assert errors[0] < 5e-3
    # convergencia O(h²): duplicar paneles reduce el error ~4x
    assert errors[1] < errors[0] / 2.5


def test_scattered_formulation_matches_total_formulation():
    """Las dos incógnitas posibles (``p`` total y ``p_s``) dan el mismo campo dispersado.

    Verificación cruzada de la ecuación integral: en la formulación ``"scattered"`` el dato de
    Neumann ``q_s = −q_inc`` entra multiplicando a ``S``, así que un error de signo o de escala en
    la capa simple se ve acá y no en el caso rígido (donde ``q = 0``).
    """
    ka = 2.0
    frequency = frequency_for(ka)
    points = circle_points(2.0 * RADIUS)

    total = solve_exterior_neumann(points, frequency, WATER, RADIUS)
    scattered = solve_exterior_neumann(points, frequency, WATER, RADIUS,
                                       formulation="scattered")
    difference = relative_error(scattered.p_scattered, total.p_scattered)

    print(f"\n  ka={ka:<4g} paneles={total.n_panels} diferencia entre formulaciones="
          f"{difference:.3e} (error de discretización del BEM: ~2.7e-3)")
    assert difference < TOLERANCE


def test_panels_for_frequency_is_power_of_two_and_grows_with_ka():
    """La malla recomendada es ``2^m``, respeta el mínimo y crece con ``ka``."""
    counts = [panels_for_frequency(RADIUS, WATER.wavenumber(frequency_for(ka)))
              for ka in (0.5, 1.0, 2.0, 5.0, 20.0)]
    print(f"\n  paneles recomendados por ka: {counts}")
    assert all(c >= 64 for c in counts)
    assert all(c & (c - 1) == 0 for c in counts)
    assert counts == sorted(counts)
    assert counts[-1] > counts[0]
