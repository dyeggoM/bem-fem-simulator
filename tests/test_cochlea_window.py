"""Fase 5b-1 — la ventana redonda como salida de energía del corte transversal.

Qué se prueba y qué se mide
---------------------------
La Fase 5a midió la transmisión **dentro de una cavidad cerrada** y por eso no convergía al ``T`` de
una placa (D18). Acá se abre la **ventana redonda** sobre la cara interna de la cápsula
(``r = a``) y se mide si eso alcanza para que la transmisión deje de ser estacionaria.

Hay dos modelos, y el negativo del primero es parte del resultado:

``pressure_release`` (lo que pidió el encargo)
    Parche de **presión nula** (``dirichlet={'p': nodos del arco}``) sobre el arco. Es la
    idealización declarada de la membrana flexible (masa e impedancia despreciables). Es una **pared
    blanda**: ``R = −1``, no absorbe. Medido a 10°: ``P_ventana/P_inc = +3.1e−3`` real contra
    ``−2.1e−2`` reactiva (×6.7), y la dispersión de ``t_interior`` en ±5 % de frecuencia **no** baja
    del 17.8 % de la cavidad cerrada a menos de 10 % (12.9 % a 10°, 16.5 % a 5°, 17.2 % a 20 %).

``open`` (diagnóstico, mismo arco)
    Se **quitan** los elementos sólidos del arco: el fluido pasa y la ABC exterior absorbe potencia
    real. Es el límite de membrana transparente. Con 20° la dispersión baja a 8.6 % (std) / 11.6 %
    (semirrango), pero 20° es cuatro veces un arco fisiológico y la cavidad sigue resonando: el
    amortiguador de banda ancha de la cóclea real es la membrana basilar (5b-2).

Criterios del encargo, uno por test
-----------------------------------
1. **Dispersión**: ``test_dispersion_pressure_release_window`` (xfail: medido 12.9 % a 10°) y
   ``test_dispersion_open_window`` (8.6 % a 20°).
2. **Energía**: ``test_energy_balance`` — la potencia neta que entra por el contorno exterior y la
   que sale por la ventana se miden por separado; con ``p = 0`` el balance cierra con **todos los
   términos en cero** y eso es el diagnóstico, no un aprobado.
3. **Tamaño**: ``test_window_size_table`` — con ``p = 0`` la presión al interior *baja* al agrandar
   el arco (no sube); en el modelo ``open`` la potencia por la ventana sí crece monótona con el
   área (0.21 % → 1.12 % → 8.84 % de ``P_inc``).
4. **No-regresión**: ``test_no_regression_closed_cavity`` — ``window_half_angle=None`` reproduce
   ``t_interior = 1.02554`` (a/t = 10, ppw = 20, t = 2 mm), residuo ``1.7e−13``.

Notas de método (para que los números sean comparables)
-------------------------------------------------------
- El barrido de frecuencia usa **la misma malla** (construida a 100 kHz y barrida con
  ``dataclasses.replace(section, frequency=f)``). Si se remalla en cada frecuencia, ``n_theta``
  cambia (492 → 508 entre 95 y 105 kHz) y con él el número de nodos del arco a 5° (**13 → 15**,
  +15 % de área de ventana): un artefacto mayor que la dispersión que se quiere medir.
- El anillo ``r = a`` es la interfaz y no sirve como anillo central de una diferencia radial (sus
  vecinos exteriores son nodos sólo sólidos). El flujo de la ventana se evalúa en el anillo de malla
  justo por dentro, ``r = 18.667 mm`` (``a − 1.33 mm``); en ``r = a`` la parte real es exactamente 0
  porque ``p = 0`` en el parche.
- El elemento de área del parche que se imprime es ``2·half·a·t`` (arco por espesor), no el área
  lateral exacta.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from ultrasound_sim.cochlea import (
    arc_power_flux,
    build_cross_section,
    intercepted_power,
    radial_power_flux,
    solve_cross_section,
    wall_injection_power,
)
from ultrasound_sim.materials import WATER

FREQUENCY = 1.0e5
THICKNESS = 2.0e-3
RADIUS = 20.0e-3
A_OVER_T = RADIUS / THICKNESS  # 10: el caso de 5a que da t_interior = 1.02554

#: Malla de 5a: ppw=20 y h = t/2 ⇒ 42 169 nodos, 45 181 gdl, residuo ~1.7e-13.
PPW = 20.0
PAD = 2.5
MEAW = 2

#: Barrido de ±5 % alrededor de 100 kHz (5 puntos).
SWEEP = (0.95, 0.975, 1.0, 1.025, 1.05)

_CACHE = {}


def _section(half=None, model="pressure_release"):
    """Malla a 100 kHz (cacheada: construir es barato, pero el barrido la reusa 5 veces)."""
    key = (half, model)
    if key not in _CACHE:
        if half is None:
            _CACHE[key] = build_cross_section(
                FREQUENCY, RADIUS, THICKNESS,
                ppw=PPW, pad_wavelengths=PAD, min_elements_across_wall=MEAW,
            )
        else:
            _CACHE[key] = build_cross_section(
                FREQUENCY, RADIUS, THICKNESS,
                ppw=PPW, pad_wavelengths=PAD, min_elements_across_wall=MEAW,
                window_half_angle=half, window_model=model,
            )
    return _CACHE[key]


def _solve(half=None, model="pressure_release", ratio=1.0, abc_order=1):
    sec = _section(half, model)
    if ratio != 1.0:
        sec = replace(sec, frequency=FREQUENCY * ratio)
    key = ("sol", half, model, ratio, abc_order)
    if key not in _CACHE:
        _CACHE[key] = solve_cross_section(
            sec, window_half_angle=half, window_model=model, abc_order=abc_order
        )
    return _CACHE[key]


def _dispersion(values):
    v = np.asarray(values, dtype=float)
    return float(v.std() / v.mean()), float((v.max() - v.min()) / (2.0 * v.mean()))


def _evaluation_radius(section):
    """Anillo de malla donde se mide el flujo de la ventana (el central del último trío)."""
    inside = section.radii[section.radii < section.a - 1e-12]
    return float(inside[-2])


# ---------------------------------------------------------------------------
# (4) No-regresión: sin ventana, los números de 5a
# ---------------------------------------------------------------------------


def test_no_regression_closed_cavity():
    """``window_half_angle=None`` tiene que reproducir la Fase 5a exacta (criterio 4)."""
    sec = _section()
    sol = _solve()
    print("\n(4) no-regresión, cavidad cerrada")
    print("    " + sec.describe())
    print("    " + sol.describe())

    assert sol.window_half_angle is None
    assert sol.window_node_count == 0
    assert sol.power_window_real_fraction == 0.0
    assert sol.power_window_reactive_fraction == 0.0
    # el número de 5a (a/t=10, ppw=20, t=2 mm, pad 2.5, 2 elementos en la pared)
    assert abs(sol.t_interior - 1.0255) < 5e-5, f"t_interior = {sol.t_interior:.6f}"
    assert sol.solution.residual < 1e-10, sol.solution.residual
    assert sol.solution.size == 45181
    print(f"    t_interior = {sol.t_interior:.6f}  (5a: 1.02554)   "
          f"residuo = {sol.solution.residual:.1e}   gdl = {sol.solution.size}")

    # pasar la ventana explícitamente como None no cambia nada
    same = solve_cross_section(sec, window_half_angle=None)
    assert same.t_interior == sol.t_interior
    assert np.allclose(same.p, sol.p, rtol=0.0, atol=0.0, equal_nan=True)


# ---------------------------------------------------------------------------
# (2) Energía: qué mide de verdad el parche p = 0
# ---------------------------------------------------------------------------


def test_energy_balance():
    """Criterio 2: el balance cierra, pero con ``p = 0`` **todos** los términos son cero.

    Se miden por separado:

    - potencia neta que entra por el contorno exterior (``r = R_out + 0.6·λ``): ``2e−14`` de
      ``P_inc`` sin ventana, ``1.9e−6`` con el parche a 10° (piso numérico del solver);
    - potencia neta almacenada/transportada en la cavidad (``r = 0.75·a``): ``2e−15`` / ``−2.0e−5``;
    - potencia que cruza el **arco de la ventana**: real ``+3.1e−3`` contra reactiva ``−2.1e−2``
      (×6.7). En ``r = a`` exacto la parte real es 0 por construcción (``p = 0`` en el parche).

    Es decir: ``|P_in − P_out| < 5 %`` se cumple, pero porque **nada** transporta potencia real.
    Eso es lo que hay que reportar del criterio, y por eso el test lo verifica como tal.

    Se imprime además la ruta del lado sólido (``wall_injection_power``), que **no** sirve para el
    balance: su parte real da ``−0.575`` cuando el valor exacto es 0 (error ``O(h)`` de la tracción
    P1) y no se cancela entre las dos caras del anillo.
    """
    sec = _section()
    closed = _solve()
    p_inc = intercepted_power(sec)
    lam = WATER.wavelength(FREQUENCY)
    r_out = sec.a + sec.thickness + 0.6 * lam

    print("\n(2) energía (criterio 2), P en unidades de P_inc")
    print(f"    {'caso':>12} {'P_net(R_out)':>14} {'P_cav(0.75a)':>14} "
          f"{'P_vent.real':>13} {'P_vent.react':>13}")
    rows = {}
    for label, sol in (("cerrada", closed), ("p=0 10°", _solve(10.0))):
        s = sol.section
        p_tot = radial_power_flux(s, sol.p, r_out) / intercepted_power(s)
        rows[label] = (p_tot, sol.power_cavity_fraction,
                       sol.power_window_real_fraction, sol.power_window_reactive_fraction)
        print(f"    {label:>12} {p_tot:+14.2e} {sol.power_cavity_fraction:+14.2e} "
              f"{sol.power_window_real_fraction:+13.2e} "
              f"{sol.power_window_reactive_fraction:+13.2e}")

    for label, (p_tot, p_cav, p_win, p_react) in rows.items():
        assert abs(p_tot) < 1e-3, f"{label}: entra {p_tot:+.2e} de P_inc por el contorno"
        assert abs(p_cav) < 1e-3, f"{label}: la cavidad no almacena ({p_cav:+.2e})"

    # la ventana p=0 sí transporta menos potencia real que reactiva (×6.7)
    p_win = rows["p=0 10°"][2]
    p_react = rows["p=0 10°"][3]
    assert abs(p_win) < 5e-2, p_win
    assert abs(p_react) > 5.0 * abs(p_win), (p_win, p_react)
    print(f"    |P_real|/|P_reactiva| en la ventana = {abs(p_win / p_react):.3f}  "
          "⇒ pared blanda, no salida")

    # y en r = a la parte real es exactamente 0 (p = 0 nodal)
    sol = _solve(10.0)
    ring = sol.section.window_nodes()
    assert np.allclose(sol.p[ring], 0.0, atol=1e-14)
    print(f"    max|p| en los {ring.size} nodos del parche = {np.abs(sol.p[ring]).max():.1e} Pa")

    # auditoría del lado sólido: la ruta NO sirve (error O(h) de la tracción P1)
    p_solid, q_solid = wall_injection_power(sec, closed.solution)
    print(f"    lado sólido (diagnóstico): p·(u·n) real = {p_solid / p_inc:+.3f} P_inc cuando el "
          f"valor exacto es 0; reactiva = {q_solid / p_inc:+.3f}  ⇒ no sirve para el balance")


def test_arc_quadrature_matches_parseval():
    """La cuadratura directa por sector tiene que coincidir con la de Parseval del círculo.

    Con ``half = 180°`` la suma por nodos es el círculo completo: se compara contra
    :func:`radial_power_flux` en el mismo anillo. Se hace sobre el caso con ventana, donde el flujo
    neto **no** es cero (``−3.5e−5`` de ``P_inc``): rel. ``4.6e−12``. En la cavidad cerrada los dos
    valores son ruido de máquina (``3.2e−15`` vs ``2.9e−15``), así que no informan.
    """
    sol = _solve(10.0)
    sec = sol.section
    r = _evaluation_radius(sec)
    p_inc = intercepted_power(sec)
    arc, _ = arc_power_flux(sec, sol.p, 180.0, radius=r)
    par = radial_power_flux(sec, sol.p, r)
    print("\ncuadratura del sector vs Parseval")
    print(f"    r = {r * 1e3:.3f} mm:  sector = {arc / p_inc:+.6e}   "
          f"Parseval = {par / p_inc:+.6e}   rel = {abs(arc / par - 1):.2e}")
    assert abs(arc / par - 1.0) < 1e-6

    closed = _solve()
    a_c, _ = arc_power_flux(closed.section, closed.p, 180.0, radius=r)
    p_c = radial_power_flux(closed.section, closed.p, r)
    print(f"    cerrada: sector = {a_c / p_inc:+.2e}  Parseval = {p_c / p_inc:+.2e}  (los dos ~0)")
    assert abs(a_c / p_inc) < 1e-12 and abs(p_c / p_inc) < 1e-12


# ---------------------------------------------------------------------------
# (1) Dispersión en frecuencia
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_dispersion_pressure_release_window():
    """Criterio 1 con el parche ``p = 0``: la dispersión **no** baja de 10 %.

    Misma malla, 5 frecuencias (±5 %). Medido (std/media / semirrango): cerrada 17.8 % / 22.8 %;
    ventana 5° 16.5 % / 21.9 %; 10° 12.9 % / 18.9 %; 20° 17.2 % / 23.5 %. La ventana *desintoniza*
    (baja el nivel medio: 1.0255 → 0.9374 con 10°) pero no amortigua: sigue siendo una cavidad
    cerrada, sólo que con una pared blanda en un arco (D18). El test verifica el hecho medido y
    falla explícitamente si alguien lo arregla sin querer.
    """
    print("\n(1) dispersión de t_interior, barrido ±5 % (misma malla)")
    print(f"    {'caso':>10} {'t(95k)':>8} {'t(97.5k)':>9} {'t(100k)':>8} "
          f"{'t(102.5k)':>10} {'t(105k)':>8} {'std/media':>10} {'semirr':>7}")
    results = {}
    for half in (None, 5.0, 10.0, 20.0):
        ts = [_solve(half, ratio=r).t_interior for r in SWEEP]
        sd, rng = _dispersion(ts)
        results[half] = (ts, sd, rng)
        name = "cerrada" if half is None else f"p=0 {half:g}°"
        print(f"    {name:>10} " + " ".join(f"{t:8.4f}" for t in ts) +
              f" {sd * 100:9.1f}% {rng * 100:6.1f}%")

    sd_closed = results[None][1]
    assert sd_closed > 0.15, f"la referencia cerrada ya no dispersa: {sd_closed:.3f}"
    for half in (5.0, 10.0, 20.0):
        assert results[half][1] >= 0.10, (
            f"a {half:g}° la dispersión bajó a {results[half][1] * 100:.1f} %: "
            "si esto pasa, el criterio 1 se cumple y hay que sacar el xfail"
        )
    pytest.xfail(
        "criterio 1 (<10 %) NO alcanzado con el parche p=0: dispersión "
        f"{results[10.0][1] * 100:.1f} % a 10° (std/media) contra {sd_closed * 100:.1f} % de la "
        "cavidad cerrada. Una pared blanda no absorbe: p=0 es R=−1, la potencia real que la cruza "
        f"es {results[10.0][2]:.0e} de P_inc. Desintoniza el nivel (-8.6 % a 10°) pero no el Q."
    )


@pytest.mark.slow
def test_dispersion_open_window():
    """El mismo arco **abierto** (sin pared) sí reduce la dispersión, pero no gratis.

    Medido: 5° 14.6 % / 19.6 %; 10° 12.9 % / 17.0 %; **20° 8.6 % / 11.6 %**. Con 20° de semiángulo
    (40° de circunferencia, ~7 mm de abertura en a = 20 mm: cuatro veces la ventana redonda real) el
    criterio < 10 % se cumple con std/media y **no** con semirrango. La conclusión del 5b-1 queda
    entonces: la ventana sola no alcanza para medir una transmisión no estacionaria; el amortiguador
    de banda ancha que falta es la membrana basilar (5b-2).
    """
    print("\n(1b) dispersión con la ventana abierta (misma malla por tamaño)")
    print(f"    {'caso':>10} {'t(95k)':>8} {'t(97.5k)':>9} {'t(100k)':>8} "
          f"{'t(102.5k)':>10} {'t(105k)':>8} {'std/media':>10} {'semirr':>7}")
    results = {}
    for half in (5.0, 10.0, 20.0):
        ts = [_solve(half, "open", ratio=r).t_interior for r in SWEEP]
        sd, rng = _dispersion(ts)
        results[half] = (ts, sd, rng)
        print(f"    {'open ' + format(half, 'g') + '°':>10} " + " ".join(f"{t:8.4f}" for t in ts) +
              f" {sd * 100:9.1f}% {rng * 100:6.1f}%")

    assert results[20.0][1] < 0.10, f"20°: {results[20.0][1] * 100:.1f} %"
    assert results[5.0][1] > results[10.0][1] > results[20.0][1], "la dispersión no baja con el área"
    print(f"    con 20° (4× el arco fisiológico): {results[20.0][1] * 100:.1f} % (std/media), "
          f"{results[20.0][2] * 100:.1f} % (semirrango)")


# ---------------------------------------------------------------------------
# (3) Tamaño de la ventana
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_window_size_table():
    """Criterio 3: la tabla de tamaños, con los dos modelos.

    Lo que se mide a 100 kHz (``r = 18.667 mm``):

    ===========  =========  ==========  ==========  ===========  ===========
    modelo       semiáng.   área[mm²]   t_int       P_vent real  P_vent reac
    ===========  =========  ==========  ==========  ===========  ===========
    cerrada      —          0           1.02554     0            0
    p = 0        5°         6.98        1.00559     +7.0e−4      −6.0e−3
    p = 0        10°        13.96       0.93739     +3.1e−3      −2.1e−2
    p = 0        20°        27.93       0.83006     +3.1e−3      −1.3e−1
    open         5°         6.98        0.97909     +2.1e−3      +1.9e−2
    open         10°        13.96       0.96394     +1.1e−2      +1.6e−2
    open         20°        27.93       0.90596     +8.8e−2      +3.6e−3
    ===========  =========  ==========  ==========  ===========  ===========

    Con ``p = 0`` la transmisión al interior **baja** monótonamente al agrandar el arco (el parche
    fija un nodo de presión y cortocircuita el campo): el criterio de monotonía creciente **no se
    cumple** en esa dirección y el test lo deja registrado. Lo que sí es monótono es la **potencia
    que cruza la ventana** en el modelo ``open`` (0.21 % → 1.12 % → 8.84 % de ``P_inc``): es la
    prueba de que la ventana funciona de salida sólo cuando es una abertura, no una pared blanda.
    """
    print("\n(3) tamaño de la ventana (100 kHz), r_eval = "
          f"{_evaluation_radius(_section()) * 1e3:.3f} mm")
    print(f"    {'modelo':>18} {'semiáng':>8} {'área[mm²]':>10} {'t_int':>9} "
          f"{'P_vent.real':>13} {'P_vent.reac':>13} {'residuo':>9}")
    table = {}
    for model in ("pressure_release", "open"):
        for half in (None, 5.0, 10.0, 20.0):
            if half is None and model == "open":
                continue
            sol = _solve(half, model) if half else _solve()
            area = 0.0 if half is None else np.radians(2 * half) * RADIUS * THICKNESS
            table[(model, half)] = sol
            name = "cerrada" if half is None else f"{model} {half:g}°"
            print(f"    {name:>18} {('—' if half is None else format(half, 'g') + '°'):>8} "
                  f"{area * 1e6:10.3f} {sol.t_interior:9.5f} "
                  f"{sol.power_window_real_fraction:+13.2e} "
                  f"{sol.power_window_reactive_fraction:+13.2e} "
                  f"{sol.solution.residual:9.1e}")
            assert sol.solution.residual < 1e-10

    print(f"    a/t = {A_OVER_T:.0f}, a = {RADIUS * 1e3:.0f} mm, t = {THICKNESS * 1e3:.1f} mm; "
          "área = 2·semiángulo·a·t")

    # p = 0: la presión al interior BAJA con el área (opuesto al criterio)
    t_p0 = [table[("pressure_release", h)].t_interior for h in (5.0, 10.0, 20.0)]
    assert t_p0[0] > t_p0[1] > t_p0[2], f"p=0 dejó de ser monótono decreciente: {t_p0}"
    print(f"    p=0  : t_int = {['%.5f' % t for t in t_p0]} → decreciente, "
          "el criterio de crecimiento monótono NO se cumple")

    # open: la potencia por la ventana SÍ crece monótona con el área
    p_open = [table[("open", h)].power_window_real_fraction for h in (5.0, 10.0, 20.0)]
    assert p_open[0] < p_open[1] < p_open[2], f"la potencia por la ventana no crece: {p_open}"
    print(f"    open : P_vent/P_inc = {['%.4f' % p for p in p_open]} → creciente; "
          "es la ventana como salida (no como adorno)")
