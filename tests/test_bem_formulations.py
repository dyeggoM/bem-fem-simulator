"""Elección de ECUACIÓN en el BEM: ``standard`` / ``burton_miller`` / ``chief`` (Fase 4B).

Contexto (``docs/decisions.md`` D11.3): la ecuación integral directa de segunda especie

    (½I − D)·p = p_inc − S·q

tiene resonancias espurias en las frecuencias propias del problema interior de Dirichlet
(``J_n(ka) = 0``). La Fase 4B barre ``ka``, así que cruza esas frecuencias y necesita una ecuación
que no sea singular ahí. Acá se verifica que Burton-Miller y CHIEF lo logran, que no cambian nada
fuera de las resonancias y que la ecuación estándar sigue dando exactamente lo de antes.

**Qué se mide y qué no.** El criterio de aceptación es el error relativo L2 del campo *dispersado*
en ``r = 2a`` contra la serie exacta (el mismo de la Fase 3), con la malla que recomienda
``panels_for_frequency``.

Sobre el fallo de la ecuación estándar hay un detalle que conviene tener escrito, porque los
números no son los que uno espera de la teoría: en un **círculo** con onda plana la resonancia
**no** degrada el campo. La condición de compatibilidad de Fredholm es
``⟨p_inc, ψ⟩ = 0`` con ``ψ = e^{inθ}`` el vector nulo, y ``⟨p_inc, e^{inθ}⟩ ∝ J_n(ka)``, que se anula
exactamente en la resonancia: el modo espurio no se excita. Lo que sí explota es el
**condicionamiento** (``cond(½I − D)`` = 3.4e2 con 64 paneles y 2.1e4 con 512, creciendo como
``1/h``), y cerca del cero de ``J₂`` el modo espurio sí se excita lo suficiente para que el error
sea ``8.9e-2`` y **no converja** (``8.6e-2``, ``8.9e-2``, ``9.1e-2``, ``9.2e-2`` con 64…512
paneles). Para una geometría sin esa simetría (o con la fuente descentrada) el fallo es el clásico
``O(1)``; con este BEM circular no es observable, así que los tests documentan los números medidos
en vez de exigir un valor que este problema no produce.
"""

import numpy as np
import pytest

from ultrasound_sim.analytic.scattering import cylinder_scattering
from ultrasound_sim.bem.helmholtz2d import (
    DEFAULT_CHIEF_POINTS,
    EQUATIONS,
    solve_exterior_neumann,
)
from ultrasound_sim.materials import WATER

RADIUS = 0.01

#: Tolerancia de aceptación (la misma de la Fase 3): 2 % en el campo dispersado.
TOLERANCE = 0.02

#: Ceros de ``J_n`` dentro del barrido de la Fase 4B (``J₀`` y ``J₂``).
RESONANT_KA = (
    (2.4048255577, "J0"),
    (5.1356223, "J2"),
)


def frequency_for(ka: float) -> float:
    """Frecuencia [Hz] que da ``ka`` para el cilindro de prueba en agua."""
    return ka * WATER.c / (2.0 * np.pi * RADIUS)


def circle_points(radius: float, n: int = 64) -> np.ndarray:
    theta = 2.0 * np.pi * np.arange(n) / n
    return radius * np.column_stack([np.cos(theta), np.sin(theta)])


def relative_error(a, b) -> float:
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


def solve_at(ka: float, equation: str = "standard", points=None, **kwargs):
    """Resuelve a ``ka`` con la malla recomendada y devuelve ``(solución, error vs serie)``."""
    frequency = frequency_for(ka)
    if points is None:
        points = circle_points(2.0 * RADIUS)
    solution = solve_exterior_neumann(points, frequency, WATER, RADIUS,
                                      equation=equation, **kwargs)
    reference = cylinder_scattering(frequency, WATER, RADIUS, points)["scattered"]
    return solution, relative_error(solution.p_scattered, reference)


# --------------------------------------------------------------------------- #
# 1) resonancias: la estándar se rompe, Burton-Miller y CHIEF no
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("ka,modo", RESONANT_KA)
def test_resonant_wavenumbers(ka, modo):
    """En los ceros de ``J_n`` BM y CHIEF siguen dentro del 2 % y la estándar no.

    Se imprimen los tres números por ``ka``. La estándar se mide contra el campo (donde el círculo
    la perdona, ver el docstring del módulo) y contra el condicionamiento (donde no).
    """
    print(f"\n  ka = {ka:.10g}  (cero de {modo})")
    results = {}
    for equation in EQUATIONS:
        solution, error = solve_at(ka, equation)
        results[equation] = (solution, error)
        print(f"    {equation:14s} error_campo={error:.3e}  cond={solution.condition_number:9.3e}  "
              f"filas={solution.n_equations:<4d} residuo={solution.residual:.1e}")

    standard, error_standard = results["standard"]
    bm, error_bm = results["burton_miller"]
    chief, error_chief = results["chief"]

    # BM y CHIEF: el criterio de aceptación, en la resonancia.
    assert error_bm < TOLERANCE, f"burton_miller: error {error_bm:.3e} >= {TOLERANCE}"
    assert error_chief < TOLERANCE, f"chief: error {error_chief:.3e} >= {TOLERANCE}"

    # La resonancia está ahí: el operador estándar se vuelve (casi) singular y BM no.
    # Medido: 3.4e2 vs 6.8 (J0) y 6.7e2 vs 7.3 (J2) con 64 y 128 paneles.
    assert standard.condition_number > 20.0 * bm.condition_number, (
        "se esperaba que la ecuación estándar quedara mucho peor condicionada que Burton-Miller; "
        f"cond(standard)={standard.condition_number:.3e}, cond(BM)={bm.condition_number:.3e}"
    )

    # Chasquido de la resonancia en el campo: en el cero de J2 el modo espurio se excita y el
    # error supera la tolerancia (y no converge con la malla). En el de J0, no (ver docstring).
    if modo == "J2":
        assert error_standard > TOLERANCE
        assert error_standard > 10.0 * error_bm


# --------------------------------------------------------------------------- #
# 2) barrido de ka: ninguna formulación se cae
# --------------------------------------------------------------------------- #


@pytest.mark.slow
def test_frequency_sweep():
    """Barrido ``ka ∈ [0.5, 8]`` en 60 puntos: peor error de cada ecuación contra la serie.

    Criterio: Burton-Miller y CHIEF < 2 % en **todo** el barrido. La estándar se reporta como
    referencia (con la malla recomendada su peor error es ~1.7e-2, en parte porque los 60 puntos no
    caen justo sobre un cero de ``J_n``; su problema real es el condicionamiento, ver test 1).
    """
    ka_values = np.linspace(0.5, 8.0, 60)
    points = circle_points(2.0 * RADIUS, 48)
    worst = {equation: (np.nan, np.nan) for equation in EQUATIONS}

    for ka in ka_values:
        frequency = frequency_for(float(ka))
        reference = cylinder_scattering(frequency, WATER, RADIUS, points)["scattered"]
        for equation in EQUATIONS:
            solution = solve_exterior_neumann(points, frequency, WATER, RADIUS,
                                              equation=equation)
            error = relative_error(solution.p_scattered, reference)
            if not np.isfinite(worst[equation][1]) or error > worst[equation][1]:
                worst[equation] = (float(ka), error)

    print()
    for equation in EQUATIONS:
        ka_worst, error_worst = worst[equation]
        print(f"  {equation:14s} peor ka={ka_worst:.4f}  error={error_worst:.3e}")

    assert worst["burton_miller"][1] < TOLERANCE
    assert worst["chief"][1] < TOLERANCE


# --------------------------------------------------------------------------- #
# 3) no-regresión: fuera de las resonancias las tres coinciden
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("ka", (1.0, 2.0, 5.0))
def test_formulations_agree_off_resonance(ka):
    """``ka = 1, 2, 5``: las tres ecuaciones coinciden dentro del error de discretización.

    Y la estándar reproduce exactamente los números de la Fase 3 (``1.3e-3``, ``2.7e-3`` y
    ``2.4e-3`` de error contra la serie); además el camino por defecto es bit a bit el mismo que
    ``equation="standard"``.
    """
    solutions = {}
    for equation in EQUATIONS:
        solutions[equation] = solve_at(ka, equation)

    standard, error_standard = solutions["standard"]
    bm, error_bm = solutions["burton_miller"]
    chief, error_chief = solutions["chief"]

    diffs = {
        "standard-BM": relative_error(bm.p_scattered, standard.p_scattered),
        "standard-CHIEF": relative_error(chief.p_scattered, standard.p_scattered),
        "BM-CHIEF": relative_error(chief.p_scattered, bm.p_scattered),
    }
    print(f"\n  ka={ka:<4g} paneles={standard.n_panels}  "
          f"error standard={error_standard:.3e} BM={error_bm:.3e} CHIEF={error_chief:.3e}")
    print("    diferencias entre formulaciones: "
          + "  ".join(f"{k}={v:.3e}" for k, v in diffs.items()))

    # Las tres ecuaciones resuelven el mismo problema: tienen que coincidir al nivel del error de
    # discretización (medido: 3e-4 … 3e-3, dominado por la parte hipersingular de BM).
    assert all(d < TOLERANCE for d in diffs.values())
    assert error_standard < TOLERANCE
    assert error_bm < TOLERANCE
    assert error_chief < TOLERANCE

    # El default no cambió: sin argumentos y con equation="standard" es el mismo cálculo.
    points = circle_points(2.0 * RADIUS)
    default = solve_exterior_neumann(points, frequency_for(ka), WATER, RADIUS)
    assert np.array_equal(default.p_scattered, standard.p_scattered)
    assert np.array_equal(default.p_boundary, standard.p_boundary)
    assert default.equation == "standard"
    assert default.n_equations == default.n_panels


# --------------------------------------------------------------------------- #
# 4) CHIEF: con 4 puntos interiores alcanza
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("ka", (1.0, 2.0, 5.0, 2.4048255577))
def test_chief_four_interior_points_are_enough(ka):
    """``M = 4`` (default) y ``M = 8`` puntos interiores dan el mismo campo (< 1e-3 relativo)."""
    four, error_four = solve_at(ka, "chief")
    eight, error_eight = solve_at(ka, "chief", chief_points=8)

    change = relative_error(four.p_scattered, eight.p_scattered)
    change_boundary = relative_error(four.p_boundary, eight.p_boundary)
    print(f"\n  ka={ka:<12g} M=4 vs M=8: campo={change:.3e} frontera={change_boundary:.3e}  "
          f"error campo: M=4 {error_four:.3e} / M=8 {error_eight:.3e}")

    assert four.n_equations == four.n_panels + DEFAULT_CHIEF_POINTS
    assert eight.n_equations == eight.n_panels + 8
    assert change < 1e-3
    assert change_boundary < 1e-3
    # Los dos M resuelven el problema físico con la misma precisión.
    assert error_four < TOLERANCE
    assert error_eight < TOLERANCE


# --------------------------------------------------------------------------- #
# 5) costo: condicionamiento y tiempo
# --------------------------------------------------------------------------- #


def test_cost_and_conditioning_at_ka5():
    """Costo y condicionamiento de las tres ecuaciones a ``ka = 5`` (se reporta, no se exige tiempo).

    El condicionamiento de la estándar crece como ``1/h`` en las resonancias; BM y CHIEF lo
    mantienen acotado (ver test 1). Fuera de las resonancias los tres condicionamientos son
    parecidos, y BM/CHIEF cuestan ~2x por la matriz hipersingular.
    """
    import time

    ka = 5.0
    frequency = frequency_for(ka)
    points = circle_points(2.0 * RADIUS)

    print()
    for equation in EQUATIONS:
        start = time.perf_counter()
        solution = solve_exterior_neumann(points, frequency, WATER, RADIUS,
                                          equation=equation)
        elapsed = time.perf_counter() - start
        print(f"  ka={ka:g} {equation:14s} cond={solution.condition_number:8.3f}  "
              f"filas={solution.n_equations:<4d} residuo={solution.residual:.1e}  "
              f"t={elapsed * 1e3:7.1f} ms")

        assert solution.residual < 1e-3
        if equation == "chief":
            assert solution.n_equations == solution.n_panels + DEFAULT_CHIEF_POINTS
        else:
            assert solution.n_equations == solution.n_panels


# --------------------------------------------------------------------------- #
# validación de argumentos
# --------------------------------------------------------------------------- #


def test_argument_validation():
    """Errores claros (en español) con ``equation``/``chief_points``/``bm_eta`` inválidos."""
    frequency = frequency_for(1.0)
    points = circle_points(2.0 * RADIUS)
    kwargs = dict(points_eval=points, frequency=frequency, material=WATER, radius=RADIUS)

    with pytest.raises(ValueError, match="equation debe ser uno de"):
        solve_exterior_neumann(**kwargs, equation="burton-miller")

    with pytest.raises(ValueError, match="solo se usa con equation='chief'"):
        solve_exterior_neumann(**kwargs, equation="standard", chief_points=4)

    with pytest.raises(ValueError, match="solo se usa con equation='burton_miller'"):
        solve_exterior_neumann(**kwargs, equation="chief", bm_eta=1j)

    with pytest.raises(ValueError, match="chief_points debe ser"):
        solve_exterior_neumann(**kwargs, equation="chief", chief_points=0)

    with pytest.raises(ValueError, match="bm_eta debe ser un escalar complejo"):
        solve_exterior_neumann(**kwargs, equation="burton_miller", bm_eta=0.0)

    # Y los valores válidos no levantan nada.
    for equation in EQUATIONS:
        solution = solve_exterior_neumann(**kwargs, equation=equation)
        assert solution.equation == equation
