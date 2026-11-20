"""Fase 5b-1b — la ventana redonda como **impedancia resistiva** (``p = Z_r·v_n``).

Por qué existe este archivo
---------------------------
La Fase 5b-1 modeló la ventana como parche ``p = 0`` y midió el negativo (D22): una pared blanda
(``R = −1``) **no transporta potencia real** (``P = p·(v·n) = 0`` en la superficie) y la dispersión
de ``t_interior`` en ±5 % de frecuencia apenas bajaba de 17.8 % a 12.9 %. La condición correcta es

    p = Z_r·(v·n)      con      v_n = (1/(i·ω·ρ_f))·∂p/∂n
    ⇒  ∂p/∂n = (i·ω·ρ_f/Z_r)·p  ≡  i·k_w·p ,   k_w = ω·ρ_f/Z_r

o sea **el mismo término de borde que la ABC-1 con ``k → k_w``**, ``−∫ (i·k_w/ρ_f)·p·q ds`` sobre las
aristas del arco (``r = a``, interfaz fluido/sólido, no frontera de la malla — por eso el caso
``'impedance'`` ensambla el sistema acoplado en ``cochlea._assemble_and_solve_arc`` en vez de llamar
a ``solve_coupled``).

Tres identidades que se usan como test (ninguna es aproximación)
----------------------------------------------------------------
1. **El ensamblado no inventa nada**: con el arco apagado tiene que dar lo mismo que
   ``solve_cross_section`` (Fase 5a) — medido: diferencia 0.0 exacta en ``p`` y en ``t_interior``.
2. **``Z_r = ρ_f·c_f`` ⇒ ``k_w = k_f``**: la condición coincide término a término con la ABC-1 sobre
   ese arco. Medido: ``t_interior`` y ``P_ventana`` idénticos (diferencia 0.0) contra
   ``solve_arc_condition(condition='abc1')``.
3. **La potencia sale y el balance cierra**: ``−P_net(R_out) = P_ventana(borde)``. Medido a 10° y
   ``Z_r = ρc``: ``−P_net = 2.068e−2`` contra ``P_ventana = 2.101e−2`` de ``P_inc`` (**1.6 %**).

Números medidos (a = 20 mm, t = 2 mm, ppw = 20, a/t = 10, 100 kHz, malla de 45 181 gdl)
---------------------------------------------------------------------------------------
======================  ==========  ===========  ===========  =================
caso                    t_interior  P_vent.real  P_vent.borde  dispersión ±5 %
======================  ==========  ===========  ===========  =================
cerrada (5a)            1.025542    —            —            17.8 %
``p = 0`` 10°           0.937392    +3.08e−3     —            12.9 %
impedancia 10° Z=ρc     0.972180    +1.70e−2     +2.10e−2     15.5 %
impedancia 20° Z=ρc     0.884500    +1.07e−1     +1.18e−1     10.7 %
impedancia 20° Z=ρc/2   0.837710    +9.43e−2     +1.25e−1      9.6 %
======================  ==========  ===========  ===========  =================

**El criterio de <10 % no se alcanza con un arco fisiológico en la cara de sombra** y el test lo
deja medido (``xfail``): con 10° la absorción es sólo 2–3 % de ``P_inc`` porque la ventana está donde
``|p|`` es bajo (``|p|_rms ≈ 0.38`` de la incidente) y la cavidad sigue siendo un resonador de Q alto
(D18). La ventana **sí amortigua** — el pico de la resonancia a 98 kHz baja de ``1.443`` a ``1.207``
con 10° y ``Z_r = ρc`` (``test_fine_sweep_shows_damping``) — pero no lo suficiente para un barrido
de ±5 % que la cruza. Sí se cumple con 20° y ``Z_r ≤ ρc`` (9.6 %) y — diagnóstico, no física de la
ventana redonda — con 10° sobre la cara **iluminada** (7.2 %, ahí es donde va la ventana **oval**).
Falta lo mismo que en 5b-1: **la membrana basilar** (5b-2), que es el amortiguador de banda ancha de
la cóclea real.

``Z_r`` es un supuesto declarado: Sim-Datos no tiene impedancia citable de la ventana redonda (D20),
así que en vez de fijar un valor se barre 0.5 / 1 / 2 × ρc y se reporta la sensibilidad.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from ultrasound_sim.cochlea import (
    _assemble_and_solve_arc,
    _measure_solution,
    build_cross_section,
    intercepted_power,
    radial_power_flux,
    solve_arc_condition,
    solve_cross_section,
    window_arc_edges,
    window_boundary_power,
)
from ultrasound_sim.fem.coupled2d import abc_boundary_block
from ultrasound_sim.materials import WATER

FREQUENCY = 1.0e5
THICKNESS = 2.0e-3
RADIUS = 20.0e-3
A_OVER_T = RADIUS / THICKNESS  # 10: el caso de 5a (t_interior = 1.025542)

PPW = 20.0
PAD = 2.5
MEAW = 2

#: Impedancia característica del fluido: el valor que hace ``k_w = k_f`` (identidad con la ABC).
Z0 = WATER.impedance

#: Barrido de ±5 % alrededor de 100 kHz (5 puntos, la misma malla: ver D22 y test_cochlea_window).
SWEEP = (0.95, 0.975, 1.0, 1.025, 1.05)

_CACHE: dict = {}


def _section(half=None, model="pressure_release", impedance=None):
    key = ("sec", half, model, impedance)
    if key not in _CACHE:
        _CACHE[key] = build_cross_section(
            FREQUENCY,
            RADIUS,
            THICKNESS,
            ppw=PPW,
            pad_wavelengths=PAD,
            min_elements_across_wall=MEAW,
            window_half_angle=half,
            window_model=model,
            window_impedance=impedance,
        )
    return _CACHE[key]


def _solve(half=None, model="pressure_release", impedance=None, ratio=1.0):
    key = ("sol", half, model, impedance, ratio)
    if key not in _CACHE:
        sec = _section(half, model, impedance)
        if ratio != 1.0:
            sec = replace(sec, frequency=FREQUENCY * ratio)
        _CACHE[key] = solve_cross_section(sec)
    return _CACHE[key]


def _dispersion(values):
    """``(std/media, semirrango/media)`` — las dos que reporta ``test_cochlea_window``."""
    v = np.asarray(values, dtype=float)
    return float(v.std() / v.mean()), float((v.max() - v.min()) / (2.0 * v.mean()))


def _arc_pressure_stats(section, p, half=10.0, center=0.0):
    """``(∫|p|² ds, L, |p|_rms)`` sobre el arco, con la masa P1 de arista."""
    edges = window_arc_edges(section, half, center)
    x = section.nodes
    pe = np.asarray(p, dtype=complex)[edges]
    lengths = np.linalg.norm(x[edges[:, 1]] - x[edges[:, 0]], axis=1)
    m_loc = np.array([[2.0, 1.0], [1.0, 2.0]]) / 6.0
    m_loc = lengths[:, None, None] * m_loc[None, :, :]
    dens = float(np.einsum("ei,ej,eij->", pe, np.conj(pe), m_loc).real)
    return dens, float(lengths.sum()), float(np.sqrt(dens / lengths.sum()))


# ---------------------------------------------------------------------------
# (0) el ensamblado del arco no introduce ni un factor ni un signo propios
# ---------------------------------------------------------------------------


def test_assembled_system_matches_reference_solver():
    """Con el arco apagado, el ensamblado nuevo tiene que ser el de ``solve_coupled`` exacto.

    ``solve_coupled`` no acepta términos de borde en aristas interiores, así que el caso
    ``'impedance'`` reensambla el sistema acoplado en ``cochlea._assemble_and_solve_arc`` con los
    bloques públicos de ``fem``. Si esa copia tuviera un factor, un signo o un bloque distinto, este
    test lo muestra: la cavidad cerrada tiene que coincidir **dígito a dígito** con la Fase 5a
    (``t_interior = 1.025542``, residuo ``1.7e−13``, 45 181 gdl).
    """
    sec = _section()
    ref = _solve()
    sol = _assemble_and_solve_arc(sec, arc_edges=np.zeros((0, 2), dtype=int), k_arc=0.0)
    mine = _measure_solution(
        sec, sol, half=None, center_angle=0.0, model="pressure_release",
        amplitude=1.0, direction=(1.0, 0.0),
    )

    dp = np.nanmax(np.abs(mine.p - ref.p)) / np.nanmax(np.abs(ref.p))
    print("\n(0) ensamblado propio vs solve_coupled (cavidad cerrada)")
    print(f"    t_interior: propio {mine.t_interior!r}  referencia {ref.t_interior!r}  "
          f"Δ = {mine.t_interior - ref.t_interior:.3e}")
    print(f"    max|Δp|/max|p| = {dp:.3e}   residuo propio = {sol.residual:.1e}   "
          f"gdl = {sol.size}")

    assert mine.t_interior == ref.t_interior
    assert np.array_equal(np.nan_to_num(mine.p), np.nan_to_num(ref.p))
    assert sol.size == ref.solution.size == 45181
    assert sol.residual < 1e-10

    # y el default de Z_r es ρ·c (el que da la identidad con la ABC)
    assert _section(10.0, "impedance").window_impedance == pytest.approx(Z0)


# ---------------------------------------------------------------------------
# (1) potencia real no nula y SALIENDO
# ---------------------------------------------------------------------------


def test_window_power_is_real_and_outgoing():
    """Criterio 1: ``P_vent.real/P_inc > 0`` con ``Z_r = ρc``, y del orden esperado.

    Se mide de tres maneras que tienen que contar la misma historia:

    - **borde (exacta)**: ``½∫|p|²·Re(1/Z_r) ds`` sobre las aristas del arco
      (:func:`window_boundary_power`), ``+2.10e−2`` de ``P_inc``;
    - **flujo radial un anillo por dentro** (:func:`arc_power_flux`, la misma vara con la que se
      midió ``p = 0``: ``+3.08e−3``): ``+1.70e−2``, o sea **5.5×** el parche de presión nula;
    - **balance**: la potencia neta que entra por ``R_out`` (``−2.07e−2``) es exactamente la que
      absorbe la ventana (``+2.10e−2``), 1.6 %.

    El orden de magnitud no es casual: ``P = ½·⟨|p|²⟩_arco·L/Z_r`` con ``|p|_rms ≈ 0.36`` de la
    incidente y ``L = 6.98 mm`` de arco. Eso se verifica contra la integral, así que el factor ½ y el
    ``1/Z`` quedan clavados.
    """
    sec = _section(10.0, "impedance")
    imp = solve_cross_section(sec)
    p0 = _solve(10.0)

    p_inc = intercepted_power(sec)
    lam = WATER.wavelength(FREQUENCY)
    r_out = sec.a + sec.thickness + 0.6 * lam
    p_net = radial_power_flux(imp.section, imp.p, r_out) / p_inc

    dens, arc_len, p_rms = _arc_pressure_stats(sec, imp.p)
    p_boundary_abs = 0.5 * (1.0 / Z0) * dens  # [W/m] = ½∫|p|²/Z ds
    estimate_abs = 0.5 * p_rms**2 * arc_len / Z0

    print("\n(1) potencia de la ventana, Z_r = ρ·c = %.3f MRayl, ±10°" % (Z0 / 1e6))
    print(f"    arco: {window_arc_edges(sec, 10.0, 0.0).shape[0]} aristas, "
          f"L = {arc_len * 1e3:.2f} mm, |p|_rms = {p_rms:.3f}·|p_inc|")
    print(f"    P_ventana/P_inc: borde = {imp.power_window_boundary_fraction:+.4e}   "
          f"flujo radial = {imp.power_window_real_fraction:+.4e}   "
          f"(p=0 10°: {p0.power_window_real_fraction:+.4e})")
    print(f"    reactiva = {imp.power_window_reactive_fraction:+.3e}   "
          f"⇒ |real|/|reactiva| = "
          f"{abs(imp.power_window_real_fraction / imp.power_window_reactive_fraction):.1f}")
    print(f"    balance: −P_net(R_out)/P_inc = {-p_net:+.4e}  vs  P_ventana(borde) = "
          f"{imp.power_window_boundary_fraction:+.4e}   "
          f"rel = {abs(-p_net / imp.power_window_boundary_fraction - 1.0) * 100:.2f} %")
    print(f"    ∫|p|²ds = {dens:.4e}  ⇒  ½·∫|p|²/Z = {p_boundary_abs / p_inc:+.4e} de P_inc   "
          f"(estimación ½·⟨|p|²⟩·L/Z = {estimate_abs / p_inc:+.4e})")

    # (a) positiva y saliendo: el balance dice que la potencia ENTRA al dominio y la absorbe la
    #     ventana; con el signo invertido el balance no cerraría.
    assert imp.power_window_boundary_fraction > 0.0
    assert imp.power_window_real_fraction > 0.0
    assert -p_net > 0.0, "la potencia neta tiene que ENTRAR al dominio (y salir por la ventana)"
    assert abs(-p_net / imp.power_window_boundary_fraction - 1.0) < 0.05

    # (b) orden esperado: 1e−2 de P_inc, no el 3e−3 (ni el ruido) del parche p = 0
    assert 1e-2 < imp.power_window_boundary_fraction < 1e-1
    assert imp.power_window_real_fraction > 3.0 * p0.power_window_real_fraction
    assert abs(p0.power_window_reactive_fraction) > 3.0 * abs(p0.power_window_real_fraction)
    # el modelo de impedancia no impone ningún nodo con p = 0: actúa sobre 26 aristas
    assert imp.window_node_count == 0
    assert imp.window_edge_count == 26
    assert imp.window_impedance == pytest.approx(Z0)

    # (c) el factor ½ y el 1/Z quedan clavados: la integral de borde y la estimación con |p|_rms
    #     coinciden (la masa P1 es exacta para esta cuadratura).
    assert abs(p_boundary_abs / estimate_abs - 1.0) < 1e-12


# ---------------------------------------------------------------------------
# (2) dispersión en frecuencia
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_dispersion_impedance_window():
    """Criterio 2: ¿baja la dispersión de ``t_interior`` en ±5 % por debajo de 10 %?

    **Medido: no, con un arco fisiológico en la cara de sombra.** Misma malla, 5 frecuencias:

    ==============  =============================================  ========  =======
    caso            95 / 97.5 / 100 / 102.5 / 105 kHz             std/media semirr.
    ==============  =============================================  ========  =======
    cerrada         0.741 / … (17.8 % global)                      17.8 %    22.8 %
    ``p = 0`` 10°   0.895 / 0.809 / 0.937 / 1.004 / 1.064          12.9 %    18.9 %
    impedancia 10°  0.782 / 0.724 / 0.972 / 1.108 / 0.994          15.5 %    21.0 %
    impedancia 20°  0.838 / 0.705 / 0.885 / 0.985 / 0.908          10.7 %    16.2 %
    ==============  =============================================  ========  =======

    La impedancia **sí amortigua** (el pico de resonancia a 98 kHz baja de ``t = 1.44`` en la
    cavidad cerrada a ``1.21``), pero con 10° absorbe sólo 2–3 % de ``P_inc`` y no alcanza para
    bajar de 10 % la dispersión de un barrido que cruza esa resonancia. El ``p = 0`` da un número
    menor **sin absorber**: desintona el modo (mueve el pico), que es otra cosa (D22).

    Lo que falta no es el modelo de la ventana —está validado en los otros tests— sino **área
    absorbente en una banda ancha**: la membrana basilar (5b-2). El criterio se cumple con 20° y
    ``Z_r ≤ ρc`` (9.6 %) y, como diagnóstico, con 10° sobre la cara **iluminada** (7.2 %; allí el
    ``|p|`` es alto, pero ésa es la ventana **oval**, no la redonda).
    """
    print("\n(2) dispersión de t_interior, barrido ±5 % (misma malla)")
    print(f"    {'caso':>18} " + " ".join(f"{r:>8.3f}" for r in SWEEP) +
          f" {'std/med':>8} {'semirr':>7} {'P_vent(borde)':>14}")

    table = {}
    cases = [
        ("cerrada", dict(half=None, model="pressure_release", imp=None)),
        ("p=0 10°", dict(half=10.0, model="pressure_release", imp=None)),
        ("impedancia 10°", dict(half=10.0, model="impedance", imp=None)),
        ("impedancia 20°", dict(half=20.0, model="impedance", imp=None)),
        ("impedancia 20° Z/2", dict(half=20.0, model="impedance", imp=0.5 * Z0)),
    ]
    for name, kw in cases:
        ts, pw = [], []
        for r in SWEEP:
            sol = _solve(kw["half"], kw["model"], kw["imp"], ratio=r)
            ts.append(sol.t_interior)
            pw.append(sol.power_window_boundary_fraction)
        sd, rng = _dispersion(ts)
        table[name] = (ts, sd, rng, float(np.mean(pw)))
        print(f"    {name:>18} " + " ".join(f"{t:8.4f}" for t in ts) +
              f" {sd * 100:7.1f}% {rng * 100:6.1f}% {np.mean(pw):14.3e}")

    # la cavidad cerrada ya dispersa; la referencia no cambió
    assert table["cerrada"][1] > 0.15
    # con 20° y Z ≤ ρc el criterio SÍ se cumple
    assert table["impedancia 20° Z/2"][1] < 0.10
    assert table["impedancia 20° Z/2"][1] < table["impedancia 20°"][1] < table["cerrada"][1]
    print(f"    ⇒ 20° con Z_r = ρc/2: {table['impedancia 20° Z/2'][1] * 100:.1f} % (std/media) "
          f"< 10 %: el criterio se alcanza con 4× el arco fisiológico y Z_r bajo")

    # y con 10° NO se cumple: se deja medido, no se relaja
    sd10 = table["impedancia 10°"][1]
    assert sd10 > 0.10
    pytest.xfail(
        f"criterio 2 (<10 %) NO alcanzado con un arco fisiológico: dispersión {sd10 * 100:.1f} % "
        f"a 10° (std/media) contra {table['cerrada'][1] * 100:.1f} % de la cavidad cerrada. La "
        f"ventana absorbe {np.mean(table['impedancia 10°'][3]) * 100:.1f} % de P_inc (potencia real, "
        "balance cerrado al 1.6 %), pero el barrido cruza la resonancia de la cavidad (pico a "
        "98 kHz: t = 1.44 cerrada → 1.21 con la ventana) y 3 % de absorción no baja el Q. Falta el "
        "amortiguador de banda ancha: la membrana basilar (5b-2)."
    )


@pytest.mark.slow
def test_fine_sweep_shows_damping():
    """La ventana impedancia **amortigua** la resonancia; el parche ``p = 0`` sólo la desintona.

    Barrido fino de 11 puntos (±5 %, paso 1 %) sobre la misma malla. El pico de la cavidad está en
    ~98 kHz y el barrido de 5 puntos del criterio 2 lo salta entre dos muestras, así que la
    dispersión sola no distingue "amortiguar" de "mover el modo". Medido:

    ==============  ===========================================  =========  ============
    caso            t en 95…105 kHz (11 puntos)                pico       P_vent/P_inc
    ==============  ===========================================  =========  ============
    cerrada         0.741 … 1.443 (98 k) … 1.031               1.443      —
    ``p = 0`` 10°   0.895 … 1.133 (98 k) … 1.017               1.133      3.1e−3 (no absorbe)
    impedancia 10°  0.782 … 1.207 (98 k) … 0.994               1.207      2.1e−2 (absorbe)
    ==============  ===========================================  =========  ============

    El pico baja 16 % con la ventana resistiva y **además** la ventana transporta 2.1 % de ``P_inc``
    con balance cerrado; el ``p = 0`` baja el pico a 1.133 sin transportar nada, o sea por
    desintonización (pared blanda, ``R = −1``, D22). Las dos cosas bajan el número del barrido
    grueso: por eso el criterio 2 se mide con la potencia y el balance, no con ``t_interior`` solo.
    """
    fine = np.linspace(0.95, 1.05, 11)
    print("\n(2b) barrido fino: ¿amortigua o desintona?")
    print(f"    {'caso':>14} " + " ".join(f"{r:6.3f}" for r in fine) + f"  {'pico':>6}")
    peaks = {}
    for name, half, model in (("cerrada", None, "pressure_release"),
                              ("p = 0 10°", 10.0, "pressure_release"),
                              ("impedancia 10°", 10.0, "impedance")):
        ts = []
        for r in fine:
            ts.append(_solve(half, model, None, ratio=r).t_interior)
        ts = np.asarray(ts)
        peaks[name] = float(ts.max())
        print(f"    {name:>14} " + " ".join(f"{t:6.3f}" for t in ts) +
              f"  {ts.max():6.3f}")

    p_c, p_i = peaks["cerrada"], peaks["impedancia 10°"]
    print(f"    pico: cerrada {p_c:.3f} → impedancia {p_i:.3f} "
          f"({(p_i / p_c - 1.0) * 100:+.1f} %) con P_ventana = "
          f"{_solve(10.0, 'impedance').power_window_boundary_fraction * 100:.1f} % de P_inc")
    assert p_i < 0.9 * p_c, "la ventana tiene que bajar el pico de la resonancia"


# ---------------------------------------------------------------------------
# (3) monotonía con el tamaño
# ---------------------------------------------------------------------------

def test_power_grows_with_window_size():
    """Criterio 3: 5°/10°/20° ⇒ la potencia real por la ventana crece (y ``t_interior`` baja).

    Es la diferencia cualitativa con ``p = 0``, que **decrecía** con el área (1.006 / 0.937 / 0.830
    para 5°/10°/20°: una pared blanda grande desintoniza más, no absorbe más — D22).
    """
    p_inc = None
    rows = {}
    print("\n(3) tamaño de la ventana, Z_r = ρc")
    print(f"    {'half':>6} {'aristas':>8} {'L [mm]':>8} {'P_vent(borde)':>14} "
          f"{'P_vent(radial)':>15} {'t_interior':>11} {'P_net(R_out)':>13}")
    for half in (5.0, 10.0, 20.0):
        sol = _solve(half, "impedance")
        sec = sol.section
        p_inc = intercepted_power(sec)
        lam = WATER.wavelength(FREQUENCY)
        p_net = radial_power_flux(sec, sol.p, sec.a + sec.thickness + 0.6 * lam) / p_inc
        edges = window_arc_edges(sec, half, 0.0)
        arc_len = float(
            np.sum(np.linalg.norm(sec.nodes[edges[:, 1]] - sec.nodes[edges[:, 0]], axis=1))
        )
        rows[half] = sol
        print(f"    {half:6.1f} {edges.shape[0]:8d} {arc_len * 1e3:8.2f} "
              f"{sol.power_window_boundary_fraction:14.4e} "
              f"{sol.power_window_real_fraction:15.4e} {sol.t_interior:11.4f} "
              f"{p_net:13.4e}")
        assert sol.power_window_boundary_fraction > 0.0
        # el balance de energía cierra al 1.6 % en los tres tamaños
        assert abs(-p_net / sol.power_window_boundary_fraction - 1.0) < 0.05

    pb = [rows[h].power_window_boundary_fraction for h in (5.0, 10.0, 20.0)]
    pr = [rows[h].power_window_real_fraction for h in (5.0, 10.0, 20.0)]
    tt = [rows[h].t_interior for h in (5.0, 10.0, 20.0)]
    print(f"    borde : {pb[0]:.3e} → {pb[1]:.3e} → {pb[2]:.3e}   (creciente: "
          f"{pb[0] < pb[1] < pb[2]})")
    print(f"    radial: {pr[0]:.3e} → {pr[1]:.3e} → {pr[2]:.3e}   (creciente: "
          f"{pr[0] < pr[1] < pr[2]})")
    assert pb[0] < pb[1] < pb[2]
    assert pr[0] < pr[1] < pr[2]
    # más área ⇒ más energía sale ⇒ menos presión media adentro
    assert tt[0] > tt[1] > tt[2]


# ---------------------------------------------------------------------------
# (4) control de consistencia con la ABC
# ---------------------------------------------------------------------------


def test_impedance_matches_abc_on_the_arc():
    """Criterio 4: con ``Z_r = ρ_f·c_f``, ``k_w = ω·ρ/Z_r = k_f`` ⇒ es la ABC sobre ese arco.

    Con ``Z_r = ρc``:

    - el coeficiente del término débil es **idéntico** al de la ABC-1 (``abc_boundary_block`` con
      ``k = k_f``): matriz igual a máquina;
    - la solución es la misma que ``solve_arc_condition(condition='abc1')``: ``t_interior`` y
      ``P_ventana`` idénticos (diferencia **0.0**). Ésa es la identidad que anuncia D22 y el test
      que atrapa un factor o un signo.

    Además se compara contra **BGT-2 sobre el arco** (``condition='abc2'``). **No es la misma
    condición**: agrega los términos de curvatura (``−1/(2a) + i/(8k a²)``) y el tangencial
    ``∂²_s``, que en un arco de ``ka ≈ 8`` aportan ``1/(2ka) ≈ 6 %``. Medido: ``t_interior``
    ``0.9722`` contra ``0.9513`` (**2.2 %**). Ese número se reporta como el peso de los términos de
    curvatura y no como un error de factor, que no existe: la ABC-1 coincide **exacta** (matriz y
    solución), y ésa es la identidad que anuncia D22 (``k_w = k_f`` ⇔ ``Z_r = ρc``).
    """
    z = Z0
    k = WATER.wavenumber(FREQUENCY)
    k_w = 2.0 * np.pi * FREQUENCY * WATER.rho / z
    sec = _section(10.0, "impedance")
    edges = window_arc_edges(sec, 10.0, 0.0)

    b_imp = abc_boundary_block(sec.nodes, edges, k_w, rho=WATER.rho)
    b_abc = abc_boundary_block(sec.nodes, edges, k, rho=WATER.rho)
    d_mat = abs(b_imp - b_abc).max() / abs(b_abc).max()

    imp = solve_cross_section(sec)
    ab1 = solve_arc_condition(sec, half_angle=10.0, condition="abc1")
    ab2 = solve_arc_condition(sec, half_angle=10.0, condition="abc2")

    # la misma integral de borde sobre los tres campos (es la potencia que la condición saca)
    def _bwin(sol):
        return window_boundary_power(sec, sol.p, impedance=Z0, half_angle=10.0,
                                     center_angle=0.0) / intercepted_power(sec)

    b_imp, b_ab1, b_ab2 = _bwin(imp), _bwin(ab1), _bwin(ab2)
    d1 = abs(imp.t_interior - ab1.t_interior) / ab1.t_interior
    d2 = abs(imp.t_interior - ab2.t_interior) / ab2.t_interior
    print("\n(4) consistencia con la ABC sobre el arco (Z_r = ρc ⇒ k_w = k_f = %.4f 1/m)"
          % k.real)
    print(f"    k_w = {k_w:.6f}  vs  k_f = {k:.6f}   max|B_imp − B_abc1|/max|B| = {d_mat:.2e}")
    print(f"    t_interior: impedance {imp.t_interior:.6f}   abc1 {ab1.t_interior:.6f}   "
          f"abc2 {ab2.t_interior:.6f}")
    print(f"    Δ vs ABC-1 = {d1:.2e}   Δ vs ABC-2 = {d2 * 100:.2f} % "
          f"(BGT-2 no es la misma condición: agrega curvatura y tangencial)")
    print(f"    P_ventana(borde)/P_inc: impedance {b_imp:+.4e}   abc1 {b_ab1:+.4e}   "
          f"abc2 {b_ab2:+.4e}")

    # la matriz es la misma: no hay factor ni signo propios
    assert d_mat < 1e-14
    # y la solución también
    assert d1 < 1e-10, f"impedance y ABC-1 sobre el arco difieren {d1:.2e}"
    assert np.allclose(np.nan_to_num(imp.p), np.nan_to_num(ab1.p), rtol=0.0,
                       atol=1e-10 * np.nanmax(np.abs(ab1.p)))
    assert imp.power_window_real_fraction == pytest.approx(
        ab1.power_window_real_fraction, rel=1e-12
    )
    assert b_imp == pytest.approx(b_ab1, rel=1e-12)
    # BGT-2 es otra condición: se reporta el peso de sus términos de curvatura (≈2 %) sin
    # presentarlo como identidad.
    assert 0.0 < d2 < 0.05


# ---------------------------------------------------------------------------
# (5) sensibilidad a Z_r
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_impedance_sensitivity():
    """Criterio 5: ``Z_r = 0.5 / 1 / 2 × ρc`` ⇒ ``t_interior`` y dispersión de cada uno.

    Tendencia esperada: **más resistencia = menos salida**. Medido a 100 kHz (arco de 10°):
    ``P_ventana/P_inc = 3.27e−2 / 3.04e−2 / 2.33e−2`` y ``t_interior = 0.958 / 0.972 / 0.989``
    (menos absorción ⇒ la cavidad vuelve hacia su 1.0255 cerrado). La dispersión sigue la misma
    dirección: ``14.7 / 15.5 / 16.2 %`` — amortiguar ayuda, pero con 10° no alcanza (criterio 2).

    ``Z_r`` no tiene valor citable (D20): esto es la sensibilidad, no una calibración. Que el orden
    sea monótono y de la magnitud de ``½∫|p|²/Z ds`` es lo que se puede afirmar.
    """
    print("\n(5) sensibilidad a Z_r (arco de 10°, 5 frecuencias)")
    print(f"    {'Z_r':>10} {'t(100k)':>9} {'P_vent(borde)':>14} {'P_vent(radial)':>15} "
          f"{'std/med':>8} {'semirr':>7}")
    rows = {}
    for fac in (0.5, 1.0, 2.0):
        z = fac * Z0
        ts, pb, pr = [], [], []
        for r in SWEEP:
            sol = _solve(10.0, "impedance", z, ratio=r)
            ts.append(sol.t_interior)
            pb.append(sol.power_window_boundary_fraction)
            pr.append(sol.power_window_real_fraction)
        sd, rng = _dispersion(ts)
        rows[fac] = dict(t100=ts[2], pb=float(np.mean(pb)), pr=float(np.mean(pr)), sd=sd)
        print(f"    {fac:4.2f}·ρc   {ts[2]:9.4f} {np.mean(pb):14.4e} {np.mean(pr):15.4e} "
              f"{sd * 100:7.1f}% {rng * 100:6.1f}%   (Z = {z / 1e6:.3f} MRayl)")

    pb = [rows[f]["pb"] for f in (0.5, 1.0, 2.0)]
    pr = [rows[f]["pr"] for f in (0.5, 1.0, 2.0)]
    tt = [rows[f]["t100"] for f in (0.5, 1.0, 2.0)]
    sd = [rows[f]["sd"] for f in (0.5, 1.0, 2.0)]
    print(f"    P_vent(borde) : {pb[0]:.3e} → {pb[1]:.3e} → {pb[2]:.3e}  ⇒ más Z, menos salida")
    print(f"    t_interior    : {tt[0]:.4f} → {tt[1]:.4f} → {tt[2]:.4f}  ⇒ más Z, más cerca "
          f"de la cavidad cerrada (1.0255)")
    print(f"    dispersión    : {sd[0] * 100:.1f}% → {sd[1] * 100:.1f}% → {sd[2] * 100:.1f}%")

    assert pb[0] > pb[1] > pb[2]
    assert pr[0] > pr[1] > pr[2]
    assert tt[0] < tt[1] < tt[2]
    assert sd[0] < sd[1] < sd[2]
    # y todas siguen siendo "salida real", no el ruido del parche p = 0 (3e−3)
    assert all(v > 1e-2 for v in pb)
    assert all(v > 1e-3 for v in pr)


# ---------------------------------------------------------------------------
# (6) no-regresión
# ---------------------------------------------------------------------------


def test_no_regression_closed_and_pressure_release():
    """Criterio 6: los números de 5a y de 5b-1 no se mueven.

    ``window_model=None``/``'pressure_release'`` siguen yendo por ``solve_coupled`` y por
    ``dirichlet``, sin tocar el ensamblado nuevo: ``t_interior = 1.025542`` con la ventana apagada
    (a/t = 10, ppw = 20, 45 181 gdl, residuo ``1.7e−13``) y ``0.937392`` con ``p = 0`` a 10°.
    """
    sec = _section()
    closed = solve_cross_section(sec)
    same = solve_cross_section(sec, window_model=None, window_half_angle=None)
    p0 = solve_cross_section(sec, window_half_angle=10.0, window_model="pressure_release")
    p0_again = _solve(10.0)

    print("\n(6) no-regresión")
    print(f"    cerrada       t_interior = {closed.t_interior:.6f}  (5a: 1.025542)   "
          f"residuo = {closed.solution.residual:.1e}   gdl = {closed.solution.size}")
    print(f"    p = 0 10°     t_interior = {p0.t_interior:.6f}  (5b-1: 0.937392)   "
          f"P_vent.real = {p0.power_window_real_fraction:+.3e}")

    assert abs(closed.t_interior - 1.025542) < 5e-6
    assert closed.solution.size == 45181
    assert closed.solution.residual < 1e-10
    assert closed.window_impedance is None
    assert closed.power_window_real_fraction == 0.0
    assert p0.t_interior == closed.t_interior or abs(p0.t_interior - 0.937392) < 5e-6
    assert p0.window_model == "pressure_release"
    assert p0.window_node_count > 0
    # el parche p=0 sigue sin transportar potencia real (D22): la real es ruido comparada con la
    # reactiva, y la de borde ni se calcula
    assert abs(p0.power_window_real_fraction) < 5e-2
    assert abs(p0.power_window_reactive_fraction) > 3.0 * abs(p0.power_window_real_fraction)
    assert p0.power_window_boundary_fraction == 0.0

    # pasar la ventana como None explícito no cambia nada (mismo camino)
    assert np.allclose(np.nan_to_num(same.p), np.nan_to_num(closed.p), rtol=0.0, atol=0.0)
    assert same.t_interior == closed.t_interior
    # y la sección con ventana pero modelo pressure_release sigue dando lo mismo
    assert p0_again.t_interior == p0.t_interior


# ---------------------------------------------------------------------------
# validaciones del modelo
# ---------------------------------------------------------------------------


def test_impedance_validation():
    """``Re Z_r > 0`` obligatorio y arco obligatorio: una parte reactiva pura no absorbe nada."""
    from ultrasound_sim.cochlea import build_cross_section as build

    with pytest.raises(ValueError, match="Re\\(Z_r\\) > 0"):
        build(FREQUENCY, RADIUS, THICKNESS, ppw=PPW, pad_wavelengths=PAD,
              min_elements_across_wall=MEAW, window_half_angle=10.0,
              window_model="impedance", window_impedance=1j * Z0)

    sec = _section(10.0, "impedance")
    with pytest.raises(ValueError, match="sin ventana"):
        solve_cross_section(replace(sec, window_half_angle=None))
    with pytest.raises(ValueError, match="window_model"):
        solve_cross_section(replace(sec, window_model="absorber"))

    # el default es ρc y se puede sobreescribir sin remallar
    other = solve_cross_section(sec, window_impedance=2.0 * Z0)
    print("\nvalidación: Z_r = 2ρc sobre la misma malla ⇒ t_interior = "
          f"{other.t_interior:.6f} (ρc: {solve_cross_section(sec).t_interior:.6f})")
    assert other.window_impedance == pytest.approx(2.0 * Z0)
    assert other.t_interior != solve_cross_section(sec).t_interior
