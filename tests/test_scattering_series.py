"""Verificación de la serie exacta de dispersión por un cilindro (Fase 3).

La serie se valida contra la física, no contra sí misma:

* sobre un cilindro **rígido** la velocidad normal del campo total debe anularse
  (``∂p/∂n = 0``) en todos los puntos de la superficie. Un error de signo en los coeficientes
  ``A_n`` o en la fase ``iⁿ`` rompe esta condición de inmediato;
* la serie incidente reconstruida con Jacobi-Anger debe coincidir con la onda plana exacta;
* el truncado debe estar convergido: duplicar ``n_terms`` no puede mover el campo.
"""

import numpy as np
import pytest

from ultrasound_sim.analytic.scattering import (
    cylinder_radial_derivative,
    cylinder_scattered_field,
    cylinder_scattering,
)
from ultrasound_sim.materials import WATER

#: Radio de trabajo [m]. Es arbitrario: la física depende solo de ``ka``.
RADIUS = 0.01

#: Valores de ``ka = Re(k)·a`` barridos (régimen adimensional, docs/plan.md).
KA_VALUES = (0.5, 1.0, 2.0, 5.0, 8.0)


def frequency_for(ka: float, material=WATER) -> float:
    """Frecuencia [Hz] que da ``ka`` con ese material y ese radio."""
    return ka * material.c / (2.0 * np.pi * RADIUS)


def surface_points(n: int = 512) -> np.ndarray:
    """``n`` puntos equiespaciados sobre ``r = RADIUS``."""
    theta = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    return RADIUS * np.column_stack([np.cos(theta), np.sin(theta)])


@pytest.mark.parametrize("ka", KA_VALUES)
def test_rigid_surface_has_zero_normal_velocity(ka):
    """``|∂p/∂n| / (|k|·|p_inc|) < 1e-6`` sobre la superficie del cilindro rígido."""
    frequency = frequency_for(ka)
    k = WATER.wavenumber(frequency)
    points = surface_points()

    dp_dr = cylinder_radial_derivative(frequency, WATER, RADIUS, points)
    residual = float(np.max(np.abs(dp_dr)) / abs(k))  # p_inc = 1 Pa

    print(f"\n  ka={ka:<4g} k={k:.3f} 1/m  residual |dp/dn|/(|k|p_inc) = {residual:.3e}")
    assert residual < 1e-6


@pytest.mark.parametrize("angle", [0.0, 0.7, -1.9])
def test_incident_series_matches_plane_wave(angle):
    """La parte incidente de Jacobi-Anger es la onda plana exacta (fija la convención de signo)."""
    frequency = frequency_for(2.0)
    points = surface_points(64)
    result = cylinder_scattering(frequency, WATER, RADIUS, points, angle=angle)

    k = result["k"]
    exact = np.exp(1j * k * (points[:, 0] * np.cos(angle) + points[:, 1] * np.sin(angle)))
    error = np.linalg.norm(result["incident"] - exact) / np.linalg.norm(exact)

    print(f"\n  ángulo={angle:+.1f} rad  error serie incidente vs exp(ik·r·cos(θ−α)) = {error:.3e}")
    assert error < 1e-12


@pytest.mark.parametrize("ka", (1.0, 5.0))
def test_convergence_doubling_n_terms(ka):
    """Duplicar ``n_terms`` mueve el campo total menos de ``1e-9`` relativo."""
    frequency = frequency_for(ka)
    points = 2.0 * RADIUS * np.column_stack([
        np.cos(np.linspace(0.0, 2.0 * np.pi, 128, endpoint=False)),
        np.sin(np.linspace(0.0, 2.0 * np.pi, 128, endpoint=False)),
    ])
    base = cylinder_scattering(frequency, WATER, RADIUS, points)
    refined = cylinder_scattering(frequency, WATER, RADIUS, points,
                                  n_terms=2 * base["n_terms"])

    change = np.linalg.norm(base["total"] - refined["total"]) / np.linalg.norm(base["total"])
    print(f"\n  ka={ka:<4g} n_terms {base['n_terms']} → {2 * base['n_terms']}: "
          f"cambio relativo = {change:.3e}")
    assert change < 1e-9


@pytest.mark.parametrize("ka", (1.0, 3.0))
def test_pressure_release_has_zero_pressure(ka):
    """En la variante Dirichlet el campo total se anula sobre la superficie."""
    frequency = frequency_for(ka)
    points = surface_points()
    result = cylinder_scattering(frequency, WATER, RADIUS, points, boundary="pressure_release")

    residual = float(np.max(np.abs(result["total"])))
    print(f"\n  ka={ka:<4g} |p_total| máx sobre la superficie (release) = {residual:.3e} Pa")
    assert residual < 1e-10


def test_scattered_field_matches_total_minus_incident():
    """``total = incident + scattered`` y ``cylinder_scattered_field`` devuelve el mismo campo."""
    frequency = frequency_for(2.0)
    points = 1.7 * surface_points(40)
    result = cylinder_scattering(frequency, WATER, RADIUS, points, angle=0.3)

    identity = np.max(np.abs(result["total"] - result["incident"] - result["scattered"]))
    shortcut = cylinder_scattered_field(frequency, WATER, RADIUS, points, angle=0.3)
    diff = np.max(np.abs(shortcut - result["scattered"]))

    print(f"\n  |total − incidente − dispersado| máx = {identity:.3e}; atajo vs dict = {diff:.3e}")
    assert identity < 1e-12
    assert diff == 0.0


def test_points_inside_cylinder_are_rejected():
    """La serie solo vale afuera del cilindro: ``r < a`` debe fallar explícitamente."""
    frequency = frequency_for(1.0)
    with pytest.raises(ValueError):
        cylinder_scattering(frequency, WATER, RADIUS, np.array([[0.0, 0.5 * RADIUS]]))
