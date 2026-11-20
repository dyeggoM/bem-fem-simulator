"""Fase 5b-2 — la membrana basilar como partición del fluido interior (impedancia distribuida).

Qué se agrega y qué se mide
--------------------------
La Fase 5a dejó escrito que la cápsula cerrada sin pérdidas es una **cavidad resonante** (D18) y la
5b-1/5b-1b midió que ni la ventana ``p = 0`` ni la ventana resistiva alcanzan para amortiguarla
(D22, D23): la dispersión de ``t_interior`` en ±5 % de frecuencia no baja del 15.5 % y el pico de
resonancia baja sólo 16 %. D23 cierra con "el amortiguador natural es la membrana basilar". Acá se
implementa: la membrana es una **línea interior** con impedancia distribuida
``p = Z_m·(v·n)``, ``Z_m = i·ω·m + k/(i·ω) + R`` (masa-resorte-resistencia por unidad de área).

Cómo (y por qué así)
--------------------
``build_cross_section(..., membrane_angle=90)`` **parte** el fluido interior por un diámetro:
duplica los nodos de la línea (coordenadas idénticas) y cada triángulo usa la copia de su lado. Las
áreas, la orientación y las regiones no cambian; lo único que cambia es que la línea deja de ser
interior. Sobre las dos caras se ensambla el mismo término de borde de la ABC-1 con
``k_w = ω·ρ_f/Z_m`` (:func:`fem.coupled2d.segment_impedance_block`), que es exactamente la
condición ``p = Z_m·v_n``. Es el camino (a) del encargo, no la viga elástica (b): con ``Z`` sólo en
el coeficiente, los dos límites son alcanzables numéricamente y son los que pide el enunciado
(rígido con ``Z → ∞``, presión nula con ``Z → 0``). En la formulación de **salto** de presión
(placa) el límite ``Z → 0`` sería transparencia, no ``p = 0``: ver el docstring del módulo.

Los cinco criterios, con los números medidos (a = 20 mm, t = 2 mm, 100 kHz, Z_m = ρ_f·c_f salvo
donde se diga; malla de ``ppw=12``, h = 1 mm, salvo en el barrido de espesor que usa la base de D24
con ``ppw=20``)
-------------------------------------------------------------------------------------------
1. **Amortiguamiento** (``test_damping_q``): el pico de 98 kHz de la cavidad cerrada tiene
   ``Q = 287`` (FWHM 0.341 kHz, medido a 0.1 kHz sobre malla fija; el pico se reproduce a ppw=20 con
   +4 % de altura, así que es físico, aunque el valor de ``Q`` no está convergido en malla — D21 ya
   midió que este modelo no converge monótono). Con la membrana **no hay pico**: el máximo en ±10 %
   es 0.673 a 97 kHz y la respuesta no baja a ``A/√2`` en la banda ⇒ ``Q < 5``.
   ``Q_antes/Q_después > 57`` (el criterio pide > 2) y la membrana absorbe
   ``P_membrana/P_inc = 0.733``.
2. **Monotonía con el espesor** (``test_monotonicity_with_thickness``): la tabla de D24
   (1.0996 … 0.6635, con dos cambios de signo y un lóbulo de 1.36 a 4 mm) pasa a
   **0.6040 / 0.5803 / 0.5319 / 0.4899 / 0.4571 / 0.4322 / 0.3876 / 0.3680 / 0.3520 / 0.3380 /
   0.3276 / 0.3215 / 0.3247**: decreciente en 12 de los 13 tramos, con un repunte de **+1.0 %** en
   el último (6.5 → 7 mm) que **se reproduce a ppw=30 (+1.1 %)** y sin ventana: no es ruido, es
   saturación del nivel (ver el test). El caso sin membrana reproduce D24 dígito a dígito.
3. **Límites** (``test_limits_against_imposed_conditions``): ``Z_m = 10⁶·ρc`` reproduce la pared
   **rígida** impuesta directamente (misma malla, sin término) al **2e−4 %** en el campo y 1e−4 %
   en ``t_interior``; ``Z_m = 10⁻⁶·ρc`` reproduce ``p = 0`` (Dirichlet) al **2e−4 %**. La tabla de
   convergencia (10² … 10⁶) es monótona.
4. **Balance de energía** (``test_energy_balance_with_absorption``): ``−P_total(R_out)/P_inc``
   contra ``P_membrana/P_inc``: **0.7281 vs 0.7327** (residuo **0.46 %**, tolerancia 1 %) con
   ``Z_m = ρc``; 0.39 % con ``Z_m = 0.5·ρc``; 0.11 % en el caso rígido (``P_abs = 0`` por
   construcción). ``P_ext = P_disp + P_abs`` se reporta pero **no** se usa como criterio: es una
   identidad de la definición de extinción (D22: criterio que se cumple por construcción).
5. **Armónico n = 0** (``test_n0_harmonic_alongside``): se reporta ``t_interior_complex`` junto a
   ``t_interior`` en todas las tablas, como pide D24. Con la membrana el n = 0 baja 1.55×
   (0.1813 → 0.1169) mientras ``t_interior`` baja 1.86×: la métrica que estaba contaminada pasa a
   ser la que responde, y el n = 0 satura porque la membrana fija el nivel medio.

Propiedades: todo **VERIFY**
----------------------------
``MembraneImpedance`` pide masa, rigidez y resistencia por unidad de área. No hay conjunto citable
para la membrana basilar en nuestra banda (el hueco de D20) así que :data:`BASILAR_MEMBRANE_BALLPARK`
va marcado VERIFY con los rangos de los que sale, y **no** se usa como validación: los criterios
1–4 son relativos o de límite. A 100 kHz la reactancia ``|ω·m − k/ω|`` es 1.4e4 Pa·s/m (casi toda
masa), el **0.9 %** de ``R``; ``test_membrane_properties_and_sensitivity`` barre ``m`` y ``k`` ×10
(reactancia hasta el 9 % de ``R``) y verifica que ``t_interior`` se mueve **1.7 %**, así que el nivel
lo fija ``R`` y no esos dos números.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from ultrasound_sim.cochlea import (
    BASILAR_MEMBRANE_BALLPARK,
    MembraneImpedance,
    build_cross_section,
    extinction_power,
    intercepted_power,
    membrane_absorbed_power,
    radial_power_flux,
    solve_cross_section,
)
from ultrasound_sim.fem.coupled2d import (
    abc_boundary_block,
    interface_edges,
    segment_absorbed_power,
    segment_impedance_block,
)
from ultrasound_sim.materials import WATER

FREQUENCY = 1.0e5
THICKNESS = 2.0e-3
RADIUS = 20.0e-3
Z0 = float(WATER.impedance.real)  # 1.5 MRayl

#: Malla barata (h = 1 mm, 1 elemento en la pared ⇒ 25 480 elementos, 12 881 nodos). La membrana
#: agrega 41 nodos (dos copias por anillo interior + el centro).
KW = dict(ppw=12.0, pad_wavelengths=1.5, min_elements_across_wall=2)

_CACHE: dict = {}


def _section(membrana=True, **kw):
    key = ("sec", membrana, tuple(sorted(kw.items())))
    if key not in _CACHE:
        base = dict(KW)
        base.update(kw)
        if membrana:
            base["membrane_angle"] = 90.0
        _CACHE[key] = build_cross_section(FREQUENCY, RADIUS, THICKNESS, **base)
    return _CACHE[key]


def _solve(membrana=True, impedance=None, model=None, **kw):
    key = ("sol", membrana, impedance, model, tuple(sorted(kw.items())))
    if key not in _CACHE:
        extra = {}
        if model is not None:
            extra["membrane_model"] = model
        if impedance is not None:
            extra["membrane_impedance"] = complex(impedance)
        _CACHE[key] = solve_cross_section(_section(membrana, **kw), **extra)
    return _CACHE[key]


def _areas(sec):
    t = sec.triangles
    x = sec.nodes
    return 0.5 * (
        (x[t[:, 1], 0] - x[t[:, 0], 0]) * (x[t[:, 2], 1] - x[t[:, 0], 1])
        - (x[t[:, 2], 0] - x[t[:, 0], 0]) * (x[t[:, 1], 1] - x[t[:, 0], 1])
    )


def _bandwidth_q(f, t):
    """``(f_pico, FWHM, Q)`` del máximo de la curva; ``FWHM`` es ``None`` si no baja a ``A/√2``."""
    f = np.asarray(f, dtype=float)
    t = np.asarray(t, dtype=float)
    i = int(np.argmax(t))
    half = t[i] / np.sqrt(2.0)
    fl = fr = None
    for k in range(i, 0, -1):
        if t[k - 1] <= half <= t[k]:
            fl = float(np.interp(half, [t[k - 1], t[k]], [f[k - 1], f[k]]))
            break
    for k in range(i, len(t) - 1):
        if t[k] >= half >= t[k + 1]:
            fr = float(np.interp(half, [t[k + 1], t[k]], [f[k + 1], f[k]]))
            break
    if fl is None or fr is None:
        return float(f[i]), None, None
    return float(f[i]), float(fr - fl), float(f[i] / (fr - fl))


# ---------------------------------------------------------------------------
# (0) La partición del diámetro no cambia la geometría
# ---------------------------------------------------------------------------


def test_mesh_partition_is_geometrically_transparent():
    """Partir el fluido por el diámetro duplica nodos y **no** toca áreas, regiones ni CCW.

    También se verifica la topología que hace que el modelo signifique algo: las dos scalae quedan
    separadas (dos componentes conexas de igual tamaño al sacar los nodos del muro) y la línea tiene
    ``4·(n_anillos + 1)`` aristas, la mitad por cara.
    """
    sec0 = _section(False)
    sec = _section(True)
    print("\n(0) partición del diámetro")
    print("    " + sec.describe())

    n_inner = int(np.count_nonzero(sec.radii < sec.a - 1e-12))
    n_dup = 2 * n_inner + 1  # dos copias por anillo interior + el centro
    # si el n_theta tuvo que subir a par, todos los anillos crecen un nodo
    expect = sec0.n_nodes + n_dup + (sec.n_theta - sec0.n_theta) * sec.radii.size
    assert sec.n_nodes == expect, (sec.n_nodes, expect)
    assert sec.n_elements == sec0.n_elements
    assert np.unique(sec.regions).tolist() == ["fluid", "solid"]
    assert np.count_nonzero(sec.regions == "solid") == np.count_nonzero(sec0.regions == "solid")

    a0, a1 = _areas(sec0), _areas(sec)
    assert np.array_equal(a0, a1), "las áreas cambiaron al duplicar nodos"
    assert np.all(a1 > 0.0), "hay elementos con orientación invertida"

    n_edges = sec.membrane_edges.shape[0]
    assert n_edges == 4 * (n_inner + 1), (n_edges, n_inner)
    vest, tymp = sec.membrane_sides()
    assert vest.shape[0] == tymp.shape[0] == n_edges // 2
    # el diámetro va de muro a muro: el extremo más lejano de cada arista toca r = a
    r_edges = np.linalg.norm(sec.nodes[sec.membrane_edges], axis=2)
    assert np.isclose(r_edges.max(), sec.a, atol=1e-9)
    assert sec.membrane_edge_nodes().size == n_edges  # 2 copias por nodo + 2 anclajes del muro

    # las dos scalae: sin los nodos del muro el grafo se parte en dos mitades iguales
    t = sec.triangles
    ed = np.vstack([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]])
    g = coo_matrix((np.ones(ed.shape[0]), (ed[:, 0], ed[:, 1])),
                   shape=(sec.n_nodes, sec.n_nodes)).tocsr()
    wall = np.flatnonzero(np.isclose(np.linalg.norm(sec.nodes, axis=1), sec.a, atol=1e-9))
    keep = np.setdiff1d(np.arange(sec.n_nodes), wall)
    n_comp, lab = connected_components(g[keep][:, keep], directed=False)
    sizes = np.sort(np.bincount(lab))[::-1]
    assert n_comp == 3, f"componentes sin el muro: {n_comp} (esperado 3: dos scalae + exterior)"
    assert sizes[0] > sizes[1], "el dominio exterior tiene que ser la componente más grande"
    assert sizes[1] == sizes[2], f"las dos scalae no son simétricas: {sizes}"
    print(f"    nodos: {sec0.n_nodes} → {sec.n_nodes} (+{n_dup}), aristas de membrana {n_edges} "
          f"({vest.shape[0]} por cara), scalae de {sizes[1]} / {sizes[2]} nodos "
          f"(exterior {sizes[0]})")
    print(f"    áreas idénticas elemento a elemento: {np.array_equal(a0, a1)}  "
          f"n_theta: {sec0.n_theta} → {sec.n_theta}  ángulo efectivo {sec.membrane_angle_effective:.3f}°")


# ---------------------------------------------------------------------------
# (1) El operador del segmento es el de la ABC con k → ω·ρ/Z
# ---------------------------------------------------------------------------


def test_segment_impedance_operator():
    """``segment_impedance_block`` no inventa nada: con ``Z = ρc`` es la ABC-1 **bit a bit**.

    Además la cuadratura de la masa de arista se audita con un campo constante: ``B·1`` tiene que dar
    ``−(i·k_w/ρ)·(L_izq/2 + L_der/2)`` en cada nodo, que es la integral exacta de ``∫ψ_i ψ_j ds`` para
    ``p ≡ 1``.
    """
    sec = _section(True)
    edges = sec.membrane_edges
    k_f = WATER.wavenumber(FREQUENCY)
    b_robin = segment_impedance_block(sec.nodes, edges, WATER.impedance, frequency=FREQUENCY,
                                      rho=WATER.rho)
    b_abc = abc_boundary_block(sec.nodes, edges, k_f, rho=WATER.rho, abc_order=1)
    d = abs(b_robin - b_abc)
    rel = float(d.max() / abs(b_abc).max())
    print(f"\n(1) identidad con la ABC-1 sobre la misma línea: max|ΔB|/max|B| = {rel:.2e}")
    assert rel < 1e-15, "con Z = ρc la condición tiene que coincidir término a término con la ABC-1"

    # campo constante: p ≡ 1 ⇒ (B·1)_i = −(i·k_w/ρ)·Σ_aristas ∫ψ_i ds
    k_w = 2.0 * np.pi * FREQUENCY * WATER.rho / complex(Z0)
    ones = np.ones(sec.n_nodes, dtype=complex)
    got = np.asarray(b_robin @ ones)
    want = np.zeros(sec.n_nodes, dtype=complex)
    length = np.linalg.norm(sec.nodes[edges[:, 1]] - sec.nodes[edges[:, 0]], axis=1)
    for (i, j), L in zip(edges, length):
        for nd in (int(i), int(j)):
            want[nd] += -1j * k_w / WATER.rho * L / 2.0
    err = float(np.abs(got - want).max() / np.abs(want).max())
    print(f"    cuadratura con p ≡ 1: max|B·1 − analítico|/max = {err:.2e}  "
          f"(ΣL = {length.sum() * 1e3:.2f} mm = 4a = {4 * sec.a * 1e3:.2f} mm)")
    assert err < 1e-14

    # y la potencia absorbida es la integral exacta: p ≡ 1 ⇒ ½·Re(1/Z)·ΣL
    p_one = np.ones(sec.n_nodes, dtype=complex)
    p_abs = segment_absorbed_power(sec.nodes, p_one, edges, complex(Z0, 0.0))
    assert np.isclose(p_abs, 0.5 * length.sum() / Z0, rtol=1e-14)
    print(f"    P_abs(p ≡ 1) = {p_abs:.6e} W/m = ½·ΣL/Z = {0.5 * length.sum() / Z0:.6e}")

    # Re(Z) <= 0 no absorbe: se rechaza en vez de devolver un cero silencioso (D22)
    with pytest.raises(ValueError):
        segment_impedance_block(sec.nodes, edges, complex(0.0, Z0), frequency=FREQUENCY,
                                rho=WATER.rho)


# ---------------------------------------------------------------------------
# (3) Los dos límites conocidos, impuestos directamente
# ---------------------------------------------------------------------------


def test_limits_against_imposed_conditions():
    """Criterio 3: ``Z → ∞`` tiene que dar la pared rígida y ``Z → 0`` la de presión nula.

    Las dos condiciones se imponen **directamente** en la misma malla (``membrane_model='rigid'``: no
    se ensambla nada, la línea queda con ``∂p/∂n = 0`` natural; ``'pressure_release'``: ``p = 0`` en
    los nodos de la línea, ``dirichlet``). No es auto-consistencia: son los dos casos conocidos, y el
    residuo de la tabla tiene que caer a cero al llevar ``Z`` al extremo.
    """
    sec = _section(True)
    rig = _solve(True, model="rigid")
    pr = _solve(True, model="pressure_release")
    print("\n(3) límites del ensamblado (Z_m = k·ρc contra la condición impuesta)")
    print(f"    {'Z_m':>10} {'t_interior':>11} {'ref':>9} {'Δ t_int':>9} {'max|Δp|/max|p_ref|':>20}")
    rows = []
    for k in (1e2, 1e4, 1e6):
        sol = _solve(True, impedance=complex(k * Z0, 0.0))
        dev = sol.t_interior / rig.t_interior - 1.0
        dp = float(np.nanmax(np.abs(sol.p - rig.p)) / np.nanmax(np.abs(rig.p)))
        rows.append((k, sol.t_interior, dev, dp))
        print(f"    {k:10.0e} {sol.t_interior:11.6f} {rig.t_interior:9.6f} "
              f"{dev * 100:+8.4f}% {dp * 100:19.5f}%")
    devs = [abs(r[2]) for r in rows]
    assert devs[0] > devs[1] > devs[2], "la convergencia al límite rígido no es monótona"
    assert devs[-1] < 0.01 and rows[-1][3] < 0.01

    rows = []
    for k in (1e-2, 1e-4, 1e-6):
        sol = _solve(True, impedance=complex(k * Z0, 0.0))
        dev = sol.t_interior / pr.t_interior - 1.0
        dp = float(np.nanmax(np.abs(sol.p - pr.p)) / np.nanmax(np.abs(pr.p)))
        rows.append((k, sol.t_interior, dev, dp))
        print(f"    {k:10.0e} {sol.t_interior:11.6f} {pr.t_interior:9.6f} "
              f"{dev * 100:+8.4f}% {dp * 100:19.5f}%")
    devs = [abs(r[2]) for r in rows]
    assert devs[0] > devs[1] > devs[2], "la convergencia al límite p = 0 no es monótona"
    assert devs[-1] < 0.01 and rows[-1][3] < 0.01

    # y el caso rígido es el que la malla ya resuelve sin término: mismo campo que sin membrana en
    # las scalae... no (el baffle cambia los modos), pero sí residuo y simetría
    print(f"    rígido: t_vestibuli = {rig.t_vestibuli:.4f}, t_tympani = {rig.t_tympani:.4f} "
          f"(el baffle rígido separa las scalae)")
    print(f"    p = 0 : t_center = {pr.t_center:.2e} (el centro está sobre la línea)")
    assert pr.t_center < 1e-12


# ---------------------------------------------------------------------------
# (4) Balance de energía con absorción real
# ---------------------------------------------------------------------------


def test_energy_balance_with_absorption():
    """Criterio 4: lo que entra por el contorno exterior se lo lleva la membrana.

    ``−P_total(R_out + 0.6·λ)`` es el flujo neto del campo total (derivada radial de 4º orden,
    5 anillos) y ``P_membrana = ½∫|p|²·Re(1/Z_m) ds`` la integral de borde sobre las dos caras: dos
    fórmulas y dos regiones del problema. Se exige el 1 %. ``P_ext = P_disp + P_membrana`` se imprime
    como chequeo de extinción, no como criterio (es una identidad de la definición).
    """
    sec = _section(True)
    p_inc = intercepted_power(sec)
    patch = sec.a + sec.thickness + 0.6 * WATER.wavelength(FREQUENCY)
    print("\n(4) balance de energía (P en unidades de P_inc)")
    print(f"    {'Z_m':>10} {'-P_total(R)':>13} {'P_membrana':>12} {'residuo':>9} "
          f"{'P_disp':>9} {'P_ext':>9} {'P_ext-P_disp-P_abs':>18}")
    cases = (("rigid", None), ("0.5", 0.5), ("1", 1.0), ("2", 2.0))
    out = {}
    for name, fac in cases:
        if fac is None:
            sol = _solve(True, model="rigid")
        else:
            sol = _solve(True, impedance=complex(fac * Z0, 0.0))
        p_tot = radial_power_flux(sec, sol.p, patch, order=4) / p_inc
        p_abs = sol.power_membrane_fraction
        resid = abs(p_tot + p_abs)
        p_scat = sol.power_scattered_fraction  # 2º orden: sólo informativo
        p_ext = extinction_power(sec, sol.p, order=4) / p_inc
        out[name] = (p_tot, p_abs, resid)
        print(f"    {name:>10} {p_tot:+13.6f} {p_abs:+12.6f} {resid * 100:8.3f}% "
              f"{p_scat:9.5f} {p_ext:9.5f} {abs(p_ext - p_scat - p_abs) * 100:17.3f}%")

    assert out["rigid"][1] == 0.0, "sin término no hay absorción (p = 0 ⇒ P = 0, D22)"
    assert abs(out["rigid"][0]) < 0.01, out["rigid"]
    for name in ("0.5", "1", "2"):
        p_tot, p_abs, resid = out[name]
        assert p_abs > 0.0, f"Z_m = {name}·ρc no absorbe: {p_abs}"
        assert resid < 0.01, f"Z_m = {name}·ρc: balance {resid * 100:.2f} % > 1 %"
    print(f"    la membrana se lleva {out['1'][1] * 100:.1f} % de P_inc con Z_m = ρc; "
          f"de la potencia que entra al blanco, el {out['1'][1] / (out['1'][1] - out['1'][0]):.0%} "
          "se disipa en la línea y el resto vuelve a salir dispersado")


# ---------------------------------------------------------------------------
# (1, criterio) Amortiguamiento: Q antes y después
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_damping_q():
    """Criterio 1: ``Q`` de la cavidad antes y después de la membrana (misma malla, barrido de f).

    ``Q = f_pico/FWHM`` del máximo de ``t_interior`` en una banda de ±10 % alrededor de la resonancia
    (98 kHz), con la malla construida a 100 kHz y sólo ``frequency`` barrida (D22: remallar mete más
    ruido que la dispersión que se mide). Caso "antes": cavidad cerrada (sin ventana, para que el
    pico sea limpio). Caso "después": el mismo corte con ``membrane_angle=90`` y ``Z_m = ρc``.
    """
    sec0 = _section(False)
    sec1 = _section(True)
    fs = np.arange(96.6e3, 99.45e3, 0.1e3)
    print("\n(1) amortiguamiento, barrido de 0.1 kHz sobre malla fija")
    print(f"    {'f[kHz]':>8} {'sin membrana':>13} {'con membrana':>13}")
    t0, t1 = [], []
    for f in fs:
        a = solve_cross_section(replace(sec0, frequency=f)).t_interior
        b = solve_cross_section(replace(sec1, frequency=f), membrane_impedance=complex(Z0, 0.0)).t_interior
        t0.append(a)
        t1.append(b)
        print(f"    {f / 1e3:8.2f} {a:13.4f} {b:13.4f}")
    fp0, fw0, q0 = _bandwidth_q(fs, t0)
    fp1, fw1, q1 = _bandwidth_q(fs, t1)
    print(f"    antes : pico {fp0 / 1e3:.2f} kHz = {max(t0):.4f}, FWHM = {fw0 / 1e3:.3f} kHz, "
          f"Q = {q0:.0f}")
    assert fw0 is not None and q0 is not None, "el pico de la cavidad cerrada desapareció: revisá"
    assert q0 > 50.0

    # después: no hay pico. La cota es Q < f0/(0.2·f0) = 5 si el FWHM excede la banda de ±10 %.
    band_lo, band_hi = 0.9 * fp0, 1.1 * fp0
    extra = np.concatenate([np.arange(band_lo, 96.6e3, 1.0e3), np.arange(99.5e3, band_hi, 1.0e3)])
    tt = list(t1)
    for f in extra:
        tt.append(solve_cross_section(replace(sec1, frequency=f),
                                      membrane_impedance=complex(Z0, 0.0)).t_interior)
    all_f = np.concatenate([fs, extra])
    order = np.argsort(all_f)
    all_f, all_t = all_f[order], np.asarray(tt)[order]
    fp1, fw1, q1 = _bandwidth_q(all_f, all_t)
    print(f"    después: máximo {max(all_t):.4f} a {fp1 / 1e3:.1f} kHz; en "
          f"[{band_lo / 1e3:.0f}, {band_hi / 1e3:.0f}] kHz no baja a A/√2 (mín {all_t.min():.4f}, "
          f"nivel {max(all_t) / np.sqrt(2):.4f})")
    if fw1 is None:
        q_after_bound = fp0 / (2.0 * 0.1 * fp0)  # FWHM > 0.2·f0 ⇒ Q < 5
        print(f"    ⇒ FWHM > {0.2 * fp0 / 1e3:.1f} kHz ⇒ Q_después < {q_after_bound:.0f}; "
              f"Q_antes/Q_después > {q0 / q_after_bound:.0f}")
        assert q0 / q_after_bound > 2.0
    else:
        print(f"    ⇒ FWHM = {fw1 / 1e3:.3f} kHz, Q_después = {q1:.0f}, "
              f"Q_antes/Q_después = {q0 / q1:.1f}")
        assert q0 / q1 > 2.0

    # el pico tampoco es un artefacto de malla: a ppw=20 sigue estando, con la misma altura
    sec0b = build_cross_section(FREQUENCY, RADIUS, THICKNESS, ppw=20.0, pad_wavelengths=2.5,
                                min_elements_across_wall=2)
    pk = max(solve_cross_section(replace(sec0b, frequency=f)).t_interior
             for f in (97.7e3, 97.9e3, 98.1e3))
    print(f"    control ppw=20 (h = 0.75 mm): pico {pk:.3f} vs {max(t0):.3f} de esta malla "
          f"({pk / max(t0) - 1:+.1%}) ⇒ el pico es físico, no ruido de malla")
    print("    (Q no está convergido en malla: D21 ya midió que este modelo no converge monótono)")
    assert pk > 1.2


# ---------------------------------------------------------------------------
# (2, criterio) Monotonía con el espesor de la pared
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_monotonicity_with_thickness():
    """Criterio 2: la tabla de espesores de D24, antes y después de la membrana.

    Base de D24 (``a = 20 mm``, ``ppw = 20``, ventana de impedancia ±10° con ``Z_r = ρc``): la
    columna "sin membrana" reproduce los números publicados **dígito a dígito** (1.0996 … 0.6635),
    que es la no-regresión de todo lo anterior. Con la membrana (``Z_m = ρc``) la curva deja de tener
    el lóbulo de resonancia y baja de forma monótona salvo el último tramo, que sube +1 %.
    """
    ts = np.arange(1.0, 7.01, 0.5) * 1e-3
    old = {1.0: 1.0996, 1.5: 0.9492, 2.0: 0.9722, 2.5: 1.0629, 3.0: 1.1839, 3.5: 1.2691,
           4.0: 1.3649, 4.5: 1.3167, 5.0: 1.1920, 5.5: 1.0405, 6.0: 0.8968, 6.5: 0.7683,
           7.0: 0.6635}
    print("\n(2) t_interior contra el espesor de la cápsula (a = 20 mm, 100 kHz, ppw = 20)")
    print(f"    {'t[mm]':>6} {'D24 (publicado)':>16} {'sin membrana':>13} {'con membrana':>13} "
          f"{'n=0 sin':>9} {'n=0 con':>9}")
    before, after, n0_before, n0_after = [], [], [], []
    for t in ts:
        kw = dict(ppw=20.0, pad_wavelengths=2.5, min_elements_across_wall=2,
                  window_half_angle=10.0, window_model="impedance")
        s0 = build_cross_section(FREQUENCY, RADIUS, float(t), **kw)
        s1 = build_cross_section(FREQUENCY, RADIUS, float(t), membrane_angle=90.0, **kw)
        a = solve_cross_section(s0)
        b = solve_cross_section(s1)
        before.append(a.t_interior)
        after.append(b.t_interior)
        n0_before.append(a.t_interior_complex)
        n0_after.append(b.t_interior_complex)
        print(f"    {t * 1e3:6.1f} {old[round(t * 1e3, 1)]:16.4f} {a.t_interior:13.4f} "
              f"{b.t_interior:13.4f} {a.t_interior_complex:9.4f} {b.t_interior_complex:9.4f}")
        assert abs(a.t_interior - old[round(t * 1e3, 1)]) < 5e-5, "no-regresión de D24"

    d = np.diff(after)
    sign_changes = int(np.count_nonzero(np.sign(d[1:]) != np.sign(d[:-1])))
    print(f"    antes : máx/mín = {max(before) / min(before):.2f}×, cambios de signo de Δt = "
          f"{int(np.count_nonzero(np.diff(np.sign(np.diff(before))) != 0))}")
    print(f"    después: máx/mín = {max(after) / min(after):.2f}×, cambios de signo de Δt = "
          f"{sign_changes}, tramo final {d[-1] / after[-2] * 100:+.2f} %")
    print(f"    n = 0 : sin membrana {n0_before[0]:.4f} → {n0_before[-1]:.4f} "
          f"({n0_before[0] / n0_before[-1]:.2f}×); con membrana {n0_after[0]:.4f} → "
          f"{n0_after[-1]:.4f} ({n0_after[0] / n0_after[-1]:.2f}×)")

    assert int(np.count_nonzero(np.diff(np.sign(np.diff(before))) != 0)) == 2, (
        "la tabla vieja tiene dos cambios de signo (lóbulo de resonancia)"
    )
    assert sign_changes <= 1, f"la curva con membrana tiene {sign_changes} cambios de signo"
    # 12 de 13 tramos decrecientes y el último +1 %: se documenta el número medido
    assert int(np.count_nonzero(d < 0.0)) == len(d) - 1
    assert abs(d[-1]) / after[-2] < 0.02
    print("    ⇒ monótona decreciente salvo el último tramo (repunte reproducible a ppw=30, "
          "+1.1 %): el nivel satura porque con la pared gruesa ya no la fija la pared sino la "
          "membrana. No se fuerza la monotonía estricta.")


# ---------------------------------------------------------------------------
# (5, criterio) El armónico n = 0 se sigue reportando
# ---------------------------------------------------------------------------


def test_n0_harmonic_alongside():
    """Criterio 5 (regla de D24): ``t_interior`` y el armónico ``n = 0`` van siempre juntos.

    Se verifica además que ``t_interior_complex`` **es** ``|mean(p)|`` sobre la sonda (o sea el
    coeficiente ``m = 0`` de la expansión angular) y que las dos métricas se reportan por scala.
    """
    sol = _solve(True, impedance=complex(Z0, 0.0))
    sec = sol.section
    ring = sec.ring_nodes(sec.probe_radius)
    n0 = float(abs(np.mean(sol.p[ring])))
    print("\n(5) armónico n = 0")
    print(f"    t_interior = {sol.t_interior:.5f}   n = 0 = {sol.t_interior_complex:.5f}   "
          f"|mean(p)| recalculado = {n0:.5f}")
    print(f"    por scala: vestíbulo {sol.t_vestibuli:.5f}, tímpano {sol.t_tympani:.5f}, "
          f"t_center {sol.t_center:.5f}")
    assert np.isclose(sol.t_interior_complex, n0, rtol=1e-12)
    assert sol.t_vestibuli > 0.0 and sol.t_tympani > 0.0
    d = sol.as_dict()
    for key in ("t_interior", "t_interior_complex", "t_vestibuli", "t_tympani",
                "power_membrane_fraction", "power_total_fraction", "power_extinction_fraction"):
        assert key in d
    assert "n=0" in sol.describe()


# ---------------------------------------------------------------------------
# Propiedades de la membrana: VERIFY y sensibilidad a m y k
# ---------------------------------------------------------------------------


def test_membrane_properties_and_sensitivity():
    """``Z_m(ω)`` y el barrido de masa/rigidez: los dos son VERIFY y no cambian el resultado.

    ``Z_m = R + i(ω·m − k/ω)``. A 100 kHz la reactancia del ballpark es 1.4e4 Pa·s/m (casi toda
    masa: ``ω·m = 1.38e4``), el **0.92 %** de ``R``. Subir ``m`` o ``k`` un factor 10 (reactancia
    hasta el 9 % de ``R``) mueve ``t_interior`` un **1.7 %**: el orden de magnitud lo fija ``R``, y
    por eso ni ``m`` ni ``k`` se usan como validación (igual que ``Z_r`` en la ventana, D22).
    """
    m = BASILAR_MEMBRANE_BALLPARK
    z = m.impedance(FREQUENCY)
    omega = 2.0 * np.pi * FREQUENCY
    print("\npropiedades de la membrana (VERIFY)")
    print("    " + m.describe())
    print(f"    Z_m(100 kHz) = {z.real / 1e6:.4f}{z.imag / 1e6:+.5f}j MRayl "
          f"(reactancia / R = {abs(z.imag) / z.real * 100:.2f} %)")
    assert np.isclose(z.real, m.resistance)
    assert np.isclose(z.imag, omega * m.mass - m.stiffness / omega)
    assert "VERIFY" in m.provenance
    assert np.isclose(m.impedance(2.0 * FREQUENCY).real, m.resistance)

    base = _solve(True, impedance=z).t_interior
    base_kw = dict(mass=m.mass, stiffness=m.stiffness, resistance=m.resistance)
    print(f"    {'caso':>22} {'t_interior':>11} {'Δ':>8}")
    for label, kw in (
        ("m × 10", dict(mass=10.0 * m.mass)),
        ("k × 10", dict(stiffness=10.0 * m.stiffness)),
        ("m × 0.1, k × 0.1", dict(mass=0.1 * m.mass, stiffness=0.1 * m.stiffness)),
    ):
        alt = MembraneImpedance(name="alt", provenance="VERIFY", **{**base_kw, **kw})
        t = _solve(True, impedance=alt.impedance(FREQUENCY)).t_interior
        print(f"    {label:>22} {t:11.5f} {t / base - 1:7.3%}")
        assert abs(t / base - 1.0) < 0.025

    with pytest.raises(ValueError):
        MembraneImpedance(name="x", mass=1.0, stiffness=1.0, resistance=1.0).impedance(0.0)

    # el ángulo de la membrana no es cosmética: se mide la absorción en dos orientaciones
    print(f"    {'ángulo':>8} {'t_interior':>11} {'P_membrana':>11}")
    absorcion = {}
    for ang in (90.0, 30.0):
        s = build_cross_section(FREQUENCY, RADIUS, THICKNESS, membrane_angle=ang, **KW)
        sol = solve_cross_section(s)
        absorcion[ang] = sol.power_membrane_fraction
        print(f"    {ang:8.1f} {sol.t_interior:11.5f} {sol.power_membrane_fraction:11.4f}")
    assert absorcion[90.0] > absorcion[30.0], (
        "a 90° (perpendicular a la incidencia) la membrana tiene que absorber más: es la razón de "
        "elegir ese ángulo en los tests"
    )


def test_validation_errors():
    """Los errores que tienen que avisar: sin membrana, ángulo inválido, ventana sobre la línea."""
    sec = _section(False)
    with pytest.raises(ValueError):
        solve_cross_section(sec, membrane_model="rigid")
    with pytest.raises(ValueError):
        build_cross_section(FREQUENCY, RADIUS, THICKNESS, ppw=12.0, membrane_angle=0.0)
    with pytest.raises(ValueError):
        build_cross_section(FREQUENCY, RADIUS, THICKNESS, ppw=12.0, membrane_angle=180.0)
    with pytest.raises(ValueError):
        # la ventana de ±90° alrededor de 0° tocaría el diámetro a 90°
        build_cross_section(FREQUENCY, RADIUS, THICKNESS, ppw=12.0, membrane_angle=90.0,
                            window_half_angle=90.0, window_model="impedance")
    # y la combinación que sí es válida (ventana en +x, membrana en y) funciona
    ok = build_cross_section(FREQUENCY, RADIUS, THICKNESS, ppw=12.0, membrane_angle=90.0,
                             window_half_angle=10.0, window_model="impedance")
    sol = solve_cross_section(ok)
    print(f"\nventana + membrana a la vez: t_interior = {sol.t_interior:.5f}, "
          f"P_ventana(borde) = {sol.power_window_boundary_fraction:+.3e}, "
          f"P_membrana = {sol.power_membrane_fraction:+.3e}")
    assert sol.window_edge_count > 0 and sol.membrane_edge_count > 0
