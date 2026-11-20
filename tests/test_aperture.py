"""Fase 4B — apertura plana sobre baffle rígido: ceros del pistón, cuadratura, 1/√r, linealidad.

Criterios (``docs/plan.md``, Fase 4B):

(a) ceros del patrón de un pistón uniforme en ``sinθ = m·λ/D``, medidos por interpolación del mínimo
    de ``|p|``, dentro de 0.5°;
(b) convergencia de cuadratura: contra una cuadratura 10× más fina, error < 1e-6 en campo lejano y
    < 1e-3 en campo cercano;
(c) el campo lejano decae como ``1/√r`` (2D) dentro del 3 %, para ``r`` entre 20λ y 40λ;
(d) linealidad: duplicar la amplitud duplica ``p`` (error < 1e-12).

Todos los tests imprimen el número medido; las aserciones son los criterios, no "cerca de".
"""

import numpy as np
import pytest

from ultrasound_sim.aperture import Aperture, far_field_pattern, pressure
from ultrasound_sim.materials import WATER, Material

FREQUENCY = 100.0e3
LAMBDA = WATER.wavelength(FREQUENCY)          # 15 mm en agua a 100 kHz

#: Radio de evaluación del campo lejano: con D = 2λ..3λ, ``r >> D²/λ`` (6 cm..9 cm) en 5 m.
FAR_RADIUS = 5.0


def parabola_minimum(x: np.ndarray, y: np.ndarray, neighbours: int = 2) -> float:
    """Vértice del ajuste cuadrático local alrededor del mínimo muestreado de ``y(x)``.

    Un ajuste sobre toda la ventana sesga el vértice (``|sinc|²`` sólo es parabólico muy cerca del
    cero), así que se ajustan los ``2·neighbours+1`` puntos que rodean al mínimo muestreado.
    """
    index = int(np.argmin(y))
    lo = max(index - neighbours, 0)
    hi = min(index + neighbours + 1, x.size)
    if hi - lo < 3:
        raise AssertionError("no hay suficientes muestras alrededor del mínimo")
    coeffs = np.polyfit(x[lo:hi], y[lo:hi], 2)
    if coeffs[0] <= 0.0:
        raise AssertionError("el ajuste no es un mínimo (curvatura <= 0)")
    return float(-coeffs[1] / (2.0 * coeffs[0]))


@pytest.mark.parametrize("widths_in_lambda", [2.0, 3.0])
def test_piston_zeros_at_sin_theta_m_lambda_over_D(widths_in_lambda):
    """Ceros del pistón uniforme: ``sinθ_m = m·λ/D``, medidos por el mínimo de ``|p|``.

    ``D = 2λ`` (30 mm a 100 kHz) es el caso del plan: un único cero visible, en 30°. Con ``D = 3λ``
    entran dos (``sinθ = 1/3`` y ``2/3``), así que la regla ``m·λ/D`` se contrasta con dos datos.
    """
    width = widths_in_lambda * LAMBDA
    aperture = Aperture(0.0, width)
    orders = [m for m in range(1, 10) if m * LAMBDA / width < 0.999]
    print(f"\n  pistón D = {widths_in_lambda:g}λ = {width * 1e3:.1f} mm, f = {FREQUENCY / 1e3:.0f} kHz, "
          f"r = {FAR_RADIUS:g} m")

    errors = []
    for m in orders:
        expected = float(np.arcsin(m * LAMBDA / width))
        span = np.deg2rad(3.0)
        theta = np.linspace(expected - span, expected + span, 1201)
        field = far_field_pattern(theta, FREQUENCY, WATER, [aperture], radius=FAR_RADIUS)
        magnitude_sq = np.abs(field) ** 2
        measured = parabola_minimum(theta, magnitude_sq)
        error_deg = abs(np.rad2deg(measured - expected))
        on_axis = abs(pressure([FAR_RADIUS, 0.0], FREQUENCY, WATER, [aperture])[0]) ** 2
        depth_db = 10.0 * np.log10(magnitude_sq.min() / on_axis)
        errors.append(error_deg)
        print(f"    m={m}: teórico sinθ={m * LAMBDA / width:.4f} -> θ={np.rad2deg(expected):7.4f}°   "
              f"medido θ={np.rad2deg(measured):7.4f}°   error={error_deg:.4f}°   "
              f"profundidad={depth_db:6.1f} dB")

    # simetría: el patrón es par en θ (apertura centrada), el cero de -θ debe coincidir
    theta = np.deg2rad(np.linspace(-(np.rad2deg(np.arcsin(orders[0] * LAMBDA / width)) + 3.0),
                                   -(np.rad2deg(np.arcsin(orders[0] * LAMBDA / width)) - 3.0), 1201))
    field = far_field_pattern(theta, FREQUENCY, WATER, [aperture], radius=FAR_RADIUS)
    mirrored = parabola_minimum(theta, np.abs(field) ** 2)
    print(f"    cero simétrico m=1: θ={np.rad2deg(mirrored):7.4f}° (esperado "
          f"{-np.rad2deg(np.arcsin(orders[0] * LAMBDA / width)):.4f}°)")

    assert max(errors) < 0.5, f"cero fuera de 0.5°: error máximo {max(errors):.4f}°"
    assert abs(np.rad2deg(mirrored) + np.rad2deg(np.arcsin(orders[0] * LAMBDA / width))) < 0.5


def test_quadrature_converges_against_10x_finer():
    """Cuadratura: ``quad_order=16, refinement=2`` vs 10× más fina (orden 160, refinement 20).

    Criterio: error relativo < 1e-6 en campo lejano y < 1e-3 en campo cercano (puntos a 0.2–2 mm de
    la apertura, donde el núcleo es casi singular y la cuadratura fija sola fallaría).
    """
    width = 2.0 * LAMBDA
    aperture = Aperture(0.0, width)
    h_base = width / 16.0
    far_points = np.array([[0.5, 0.0], [0.5, 0.02], [0.2, 0.01], [1.0, -0.03], [0.05, 0.0]])
    near_points = np.array([[0.002, 0.0], [0.001, 0.005], [0.0005, -0.002], [0.0002, 0.014]])
    print(f"\n  D = {width * 1e3:.1f} mm, subsegmento base = {h_base * 1e3:.3f} mm "
          f"(umbral de cercanía {2 * h_base * 1e3:.3f} mm)")

    worst = {}
    for label, points, tolerance in (("lejano", far_points, 1e-6), ("cercano", near_points, 1e-3)):
        coarse = pressure(points, FREQUENCY, WATER, [aperture])
        fine = pressure(points, FREQUENCY, WATER, [aperture], quad_order=160, refinement=20.0)
        error = np.abs(coarse - fine) / np.abs(fine)
        worst[label] = float(error.max())
        print(f"    {label}: error relativo máximo = {error.max():.3e}  (criterio {tolerance:g})")
        for p, e, value in zip(points, error, fine):
            print(f"       x={p[0] * 1e3:7.4f} mm y={p[1] * 1e3:8.4f} mm  "
                  f"|p|={abs(value):.6e} Pa  err={e:.3e}")
        assert worst[label] < tolerance

    # la cuadratura reportada no debe depender del tamaño del bloque vectorizado
    blocked = pressure(np.vstack([far_points, near_points]), FREQUENCY, WATER, [aperture])
    print(f"    (evaluación en bloque de los {len(far_points) + len(near_points)} puntos: "
          f"error del bloque vs fina = "
          f"{np.max(np.abs(blocked - np.concatenate([pressure(far_points, FREQUENCY, WATER, [aperture]), pressure(near_points, FREQUENCY, WATER, [aperture])])) / np.abs(blocked)):.3e})")

    # contraste: sin graduación (refinement → 0 desactiva la subdivisión) la Gauss fija falla
    naive = pressure(near_points, FREQUENCY, WATER, [aperture], refinement=1e-9)
    fine_near = pressure(near_points, FREQUENCY, WATER, [aperture], quad_order=160, refinement=20.0)
    naive_error = np.abs(naive - fine_near) / np.abs(fine_near)
    print(f"    contraste sin graduación (Gauss fijo, 16 nodos por apertura): error máximo = "
          f"{naive_error.max():.2e}  ⇒ la subdivisión es lo que lleva el campo cercano a 1e-16")
    assert naive_error.max() > 1e-3


def test_far_field_decays_as_inverse_sqrt_r():
    """En 2D el campo lejano decae como ``1/√r``: ``|p|·√r`` constante dentro del 3 %."""
    aperture = Aperture(0.0, 2.0 * LAMBDA)
    radii = LAMBDA * np.linspace(20.0, 40.0, 41)
    amplitude = np.array([abs(pressure([r, 0.0], FREQUENCY, WATER, [aperture])[0]) for r in radii])
    scaled = amplitude * np.sqrt(radii)
    scale = float(np.sum(amplitude / np.sqrt(radii)) / np.sum(1.0 / radii))   # LS de |p| = C/√r
    deviation = np.abs(amplitude / (scale / np.sqrt(radii)) - 1.0)
    print(f"\n  pistón D=2λ, eje x, r = {radii[0] * 1e3:.1f}..{radii[-1] * 1e3:.1f} mm "
          f"({20:g}λ..{40:g}λ)")
    print(f"    C (LS) = {scale:.6e} Pa·√m;  |p|·√r relativo: min = {scaled.min() / scaled.mean():.6f}, "
          f"max = {scaled.max() / scaled.mean():.6f}")
    print(f"    desvío máximo respecto de C/√r = {deviation.max() * 100:.4f} %  "
          f"(criterio 3 %);  |p·√r|(20λ)/|p·√r|(40λ) - 1 = "
          f"{scaled[0] / scaled[-1] - 1:+.4%}")
    assert deviation.max() < 0.03


def test_linearity_in_amplitude():
    """``p ∝ v_n``: duplicar la amplitud duplica la presión (error < 1e-12) y ``i`` rota 90°."""
    points = np.array([[0.5, 0.0], [0.3, 0.01], [0.002, 0.0], [0.02, 0.005]])
    width = 2.0 * LAMBDA
    p1 = pressure(points, FREQUENCY, WATER, [Aperture(0.0, width, 1.0)])
    p2 = pressure(points, FREQUENCY, WATER, [Aperture(0.0, width, 2.0)])
    p_neg = pressure(points, FREQUENCY, WATER, [Aperture(0.0, width, -1.0)])
    p_i = pressure(points, FREQUENCY, WATER, [Aperture(0.0, width, 1.0j)])
    error = float(np.max(np.abs(2.0 * p1 - p2)) / np.max(np.abs(p2)))
    error_neg = float(np.max(np.abs(p1 + p_neg)) / np.max(np.abs(p1)))
    error_i = float(np.max(np.abs(1.0j * p1 - p_i)) / np.max(np.abs(p1)))
    print(f"\n  duplicar amplitud: error relativo = {error:.3e}   (-amplitud: {error_neg:.3e})   "
          f"(i·amplitud: {error_i:.3e})")
    assert error < 1e-12
    assert error_neg < 1e-12
    assert error_i < 1e-12

    # superposición: dos aperturas contiguas suman exactamente igual que sumar sus campos
    left = Aperture(-LAMBDA, LAMBDA)
    right = Aperture(+LAMBDA, LAMBDA)
    both = pressure(points, FREQUENCY, WATER, [left, right])
    separately = pressure(points, FREQUENCY, WATER, [left]) + pressure(points, FREQUENCY, WATER, [right])
    print(f"  superposición (2 aperturas): error relativo = "
          f"{np.max(np.abs(both - separately)) / np.max(np.abs(both)):.3e}")
    assert np.max(np.abs(both - separately)) / np.max(np.abs(both)) < 1e-12


def test_attenuation_enters_through_complex_wavenumber():
    """Con ``k`` complejo el campo decae como ``exp(−α·Δr)`` además del ``1/√r`` (D4)."""
    lossy = Material(name="lossy_test", c=WATER.c, rho=WATER.rho, alpha0=40.0, alpha_power=1.0,
                     f_ref=1.0e6, provenance="parámetro de test (α=4 Np/m a 100 kHz), no es un dato físico")
    alpha = lossy.attenuation(FREQUENCY)
    k = lossy.wavenumber(FREQUENCY)
    aperture = Aperture(0.0, 2.0 * LAMBDA)
    r1, r2 = 0.30, 0.60
    p1 = abs(pressure([r1, 0.0], FREQUENCY, lossy, [aperture])[0])
    p2 = abs(pressure([r2, 0.0], FREQUENCY, lossy, [aperture])[0])
    measured = p2 / p1
    predicted = np.exp(-alpha * (r2 - r1)) / np.sqrt(r2 / r1)
    print(f"\n  α = {alpha:.3f} Np/m, k = {k.real:.2f} + {k.imag:.2f}i 1/m")
    print(f"    |p(0.60)|/|p(0.30)| medido = {measured:.6f}   "
          f"predicho (1/√(r2/r1)·exp(-αΔr)) = {predicted:.6f}   "
          f"error = {abs(measured / predicted - 1):.3e}")
    assert k.imag > 0.0
    assert abs(measured / predicted - 1.0) < 0.02


def test_points_on_or_behind_the_baffle_are_rejected():
    """El modelo sólo vale en ``x > 0``: sobre el baffle o detrás tiene que fallar con mensaje."""
    aperture = Aperture(0.0, 2.0 * LAMBDA)
    with pytest.raises(ValueError, match="x > 0"):
        pressure([0.0, 0.0], FREQUENCY, WATER, [aperture])
    with pytest.raises(ValueError, match="x > 0"):
        pressure(np.array([[0.01, 0.0], [-0.01, 0.0]]), FREQUENCY, WATER, [aperture])
    with pytest.raises(ValueError, match=r"\|θ\| < π/2"):
        far_field_pattern(np.deg2rad(90.0), FREQUENCY, WATER, [aperture])

    # un punto suelto (2,) y una lista de puntos dan lo mismo
    single = pressure([0.1, 0.005], FREQUENCY, WATER, [aperture])
    listed = pressure(np.array([[0.1, 0.005]]), FREQUENCY, WATER, [aperture])
    print(f"\n  punto suelto vs lista: |Δp| = {abs(single[0] - listed[0]):.3e} Pa")
    assert single.shape == (1,)
    assert single[0] == pytest.approx(listed[0], rel=1e-15)
