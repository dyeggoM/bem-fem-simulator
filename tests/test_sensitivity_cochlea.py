"""Sensibilidad global del corte de cápsula coclear (``scripts/sensitivity_cochlea.py``).

Qué se verifica (los cuatro criterios del encargo)
--------------------------------------------------
1. **El caso base del script reproduce el modelo.** ``t_interior`` del ensamblado propio (con
   ``D`` complejo, que a ``α = 0`` es el mismo sistema) contra ``solve_cross_section`` del repo, y
   contra el valor ancla ``0.972180`` (D22) con tolerancia del 1 %.
2. **Monotonía con el espesor.** El criterio que pide la física («más hueso, menos pasa») se
   verifica donde es cierto —``|T|`` de la placa equivalente y el armónico n=0 de la cavidad— y se
   deja **medido y explicado** que ``mean(|p|)`` **no** es monótona en un rango 1–7 mm: es una
   amplitud de onda estacionaria de una cavidad cerrada (D18), no un coeficiente de transmisión.
   El test marca ``xfail`` con los números, como hace ``test_cochlea_impedance`` con el criterio
   de dispersión que no se alcanza (D22).
3. **La atenuación siempre reduce** ``t_interior``, para las tres leyes de ``attenuation.py``, y
   el orden sigue a ``α``.
4. **El ranking se imprime y se guarda**, y en la grilla reducida el parámetro más influyente es
   el **espesor** (si no lo fuera, el test falla: sería un hallazgo, no un detalle).

Costos (medidos en este equipo)
-------------------------------
Un solve son ~2 s con ``MMD_AT_PLUS_A`` (43.5 s con el COLAMD que trae scipy por defecto). El
archivo completo corre en ~1 min porque los escenarios se cachean entre tests; el barrido grande
vive en el script (``--step 0.5``), no acá.
"""

from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

from ultrasound_sim.cochlea import build_cross_section, solve_cross_section

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "sensitivity_cochlea.py"

#: Valores medidos con el script completo (ppw=20, paso 0.5 mm); se citan en los mensajes.
T_BASE = 0.972180
T_THICKNESS_1_TO_7_MM = (1.099587, 0.972180, 1.183929, 1.364903, 1.191993, 0.896818, 0.663523)


def test_el_barrido_con_membrana_reproduce_la_tabla_de_d29():
    """``--membrane`` tiene que dar los números de D29/D30, no parecerse a los de D24.

    A t = 2 mm y ppw = 20: ``0.972180`` sin membrana (ancla de D22) y ``0.531870`` con la membrana
    basilar. Si el bloque de borde se cableara mal —o no se aplicara— el valor volvería a ~0.97 y
    este test lo detecta; el ranking de D30 entero depende de que este número sea el correcto.
    """
    from ultrasound_sim.cochlea import BASILAR_MEMBRANE_BALLPARK

    # setattr: el módulo se carga con importlib, así que el chequeador de tipos no ve sus globales
    setattr(SC, "MEMBRANE_IMPEDANCE", None)
    sin = SC.evaluate(SC.Scenario(), ppw=20.0)["t_interior"]
    setattr(SC, "MEMBRANE_IMPEDANCE", BASILAR_MEMBRANE_BALLPARK.impedance(SC.FREQUENCY_HZ))
    try:
        con = SC.evaluate(SC.Scenario(), ppw=20.0)
    finally:
        setattr(SC, "MEMBRANE_IMPEDANCE", None)   # no contaminar los demás tests

    print(f"\n  t = 2 mm: sin membrana {sin:.5f} | con membrana {con['t_interior']:.5f} "
          f"| n=0 {con['t_interior_complex']:.5f}")
    assert abs(sin - T_BASE) < 0.01 * T_BASE, f"el caso base sin membrana se movió: {sin:.5f}"
    assert con["t_interior"] == pytest.approx(0.53187, rel=2e-3), (
        f"con membrana esperaba 0.5319 (D29), medí {con['t_interior']:.5f}"
    )
    assert con["t_interior_complex"] == pytest.approx(0.17913, rel=2e-3)


def _load_script():
    """Carga ``scripts/sensitivity_cochlea.py`` como módulo (no es un paquete)."""
    name = "sensitivity_cochlea"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, SCRIPT)
    if spec is None or spec.loader is None:
        raise ImportError(f"no se pudo cargar {SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module  # dataclass necesita el módulo registrado
    spec.loader.exec_module(module)
    return module


SC = _load_script()


def _t(**kw) -> float:
    """``t_interior`` de un escenario (usa el caché del script)."""
    return float(SC.evaluate(SC.Scenario(**kw))["t_interior"])


def _monotone_decreasing(values) -> bool:
    v = list(values)
    return all(b < a for a, b in zip(v, v[1:]))


# --------------------------------------------------------------------------- #
# (1) el caso base reproduce el modelo del repo
# --------------------------------------------------------------------------- #


def test_base_case_matches_solve_cross_section():
    """Criterio 1: el ensamblado del script da el mismo ``t_interior`` que el solver del repo.

    Dos comprobaciones independientes:

    - contra ``solve_cross_section(section)`` con **la misma malla y los mismos parámetros**
      (tolerancia 1e-6 relativa: con ``α = 0`` el sistema es el mismo y sólo cambia el pivoteo de
      SuperLU, así que la diferencia tiene que ser de redondeo);
    - contra el ancla documentada ``0.972180`` (D22, ``tests/test_cochlea_impedance.py``) con el
      1 % que pide el encargo, para que el número base no se mueva en silencio.
    """
    ours = SC.evaluate(SC.Scenario())
    section = build_cross_section(
        SC.FREQUENCY_HZ,
        SC.RADIUS_M,
        SC.BASE_THICKNESS_M,
        ppw=SC.PPW,
        pad_wavelengths=SC.PAD_WAVELENGTHS,
        min_elements_across_wall=SC.MIN_ELEMENTS_ACROSS_WALL,
        bone=SC.Scenario().bone(),
        window_half_angle=SC.WINDOW_HALF_ANGLE_DEG,
        window_center_angle=SC.WINDOW_CENTER_ANGLE_DEG,
        window_model="impedance",
    )
    ref = solve_cross_section(section)

    d_rel = abs(ours["t_interior"] / ref.t_interior - 1.0)
    d_anchor = abs(ours["t_interior"] / T_BASE - 1.0)
    print("\n(1) caso base")
    print(f"    script  t_interior = {ours['t_interior']:.9f}")
    print(f"    repo    t_interior = {ref.t_interior:.9f}   Δ rel = {d_rel:.2e}")
    print(f"    ancla   t_interior = {T_BASE:.6f}            Δ rel = {d_anchor:.2e}")
    print(
        f"    gdl = {ours['gdl']}   nodos = {ours['nodes']}   "
        f"residuo repo = {ref.solution.residual:.1e}"
    )

    assert d_rel < 1e-6, f"el ensamblado del script difiere del solver del repo en {d_rel:.2e}"
    assert d_anchor < 0.01, f"el caso base se movió {100 * d_anchor:.2f} % del ancla 0.972180"
    assert section.window_impedance == pytest.approx(SC.Z_FLUID)


# --------------------------------------------------------------------------- #
# (2) monotonía con el espesor
# --------------------------------------------------------------------------- #


def test_thickness_sweep_monotonicity_is_a_standing_wave_not_a_transmission():
    """Criterio 2: «más hueso, menos pasa». Se cumple en la placa y en el armónico n=0; **no** en
    ``mean(|p|)``, y eso es el hallazgo del análisis.

    Medido (ppw=20, ventana 10°, ``Z_r = ρc``, sin atenuación):

    ==========  =========  ==========  ================
    t [mm]      mean|p|    |mean p|    |T| placa
    ==========  =========  ==========  ================
    1           1.099587   0.386602    0.932700
    2           0.972180   0.254662    0.795900
    3           1.183929   0.193594    0.669200
    4           1.364903   0.165707    0.574600
    5           1.191993   0.122031    0.508600
    6           0.896818   0.083105    0.464500
    7           0.663523   0.076085    0.436800
    ==========  =========  ==========  ================

    ``|mean(p)|`` y ``|T|`` bajan siempre (transmisión); ``mean(|p|)`` sube entre 1 y 4 mm porque
    la cavidad cerrada es un resonador y su nivel medio es una **amplitud de onda estacionaria**
    (D18). No es malla: a ppw=30 los mismos tres puntos (2/4/6 mm) dan 0.9805 / 1.2938 / 0.8533,
    o sea el mismo patrón con el mismo signo de las diferencias.

    Por eso el test **afirma** la monotonía donde la física la exige y deja el criterio tal como
    se pidió en ``xfail`` con los números, en vez de relajarlo a escondidas.
    """
    thicknesses = (1.0e-3, 2.0e-3, 3.0e-3, 4.0e-3, 5.0e-3, 6.0e-3, 7.0e-3)
    res = [SC.evaluate(SC.Scenario(thickness_m=t)) for t in thicknesses]
    mean_abs = [r["t_interior"] for r in res]
    breathing = [r["t_interior_complex"] for r in res]
    plate = [r["plate_T"] for r in res]

    print("\n(2) espesor 1…7 mm")
    print("      t [mm]   mean|p|   |mean p|   |T| placa")
    for t, a, b, c in zip(thicknesses, mean_abs, breathing, plate):
        print(f"      {t * 1e3:5.0f}   {a:8.6f}  {b:8.6f}   {c:8.6f}")
    print(f"    mean|p|    monótona decreciente: {_monotone_decreasing(mean_abs)}")
    print(f"    |mean p|   monótona decreciente: {_monotone_decreasing(breathing)}")
    print(f"    |T| placa  monótona decreciente: {_monotone_decreasing(plate)}")

    # (a) lo que la física exige, y que el solver cumple: una placa más gruesa transmite menos
    assert _monotone_decreasing(plate), "|T| de la placa tiene que bajar con el espesor"
    # (b) el armónico que se comporta como transmisión: la «respiración» n=0 de la cavidad
    assert _monotone_decreasing(breathing), "|mean(p)| tiene que bajar con el espesor"
    # (c) y el número medido, para que la explicación no se separe de los datos
    assert not _monotone_decreasing(mean_abs)
    assert mean_abs[3] > mean_abs[1]  # 4 mm por encima de 2 mm: el pico de la resonancia
    assert np.allclose(mean_abs, T_THICKNESS_1_TO_7_MM, rtol=0.0, atol=5e-6)

    pytest.xfail(
        "criterio 2 en su forma literal NO se cumple: mean(|p|) en r=a/2 no es monótona en el "
        f"espesor (1…7 mm: {', '.join(f'{v:.4f}' for v in mean_abs)}; sube 1→4 mm y baja 4→7 mm). "
        "No es un artefacto de malla (ppw=30 da el mismo patrón) ni un error de ensamblado "
        "(criterio 1 reproduce el solver al dígito): la cavidad es cerrada y sin pérdidas, así que "
        "mean(|p|) es una amplitud de onda estacionaria (D18). La monotonía que pide la física se "
        "verifica en |T| de la placa (0.933→0.437) y en el armónico n=0 (0.387→0.076), que es la "
        "magnitud que se comporta como transmisión."
    )


# --------------------------------------------------------------------------- #
# (3) la atenuación reduce siempre
# --------------------------------------------------------------------------- #


def test_attenuation_always_reduces_interior_pressure():
    """Criterio 3: las tres leyes de ``attenuation.py`` bajan ``t_interior`` respecto de sin
    pérdidas, y más ``α`` ⇒ menos presión interior.

    ``α(100 kHz)`` por ley (extrapolación declarada: D20 no tiene dato bajo 270 kHz): White 2006
    ``3.78 Np/m``, Pichardo 2011 ``6.24 Np/m`` (interpolación log-log de sus puntos, por debajo del
    primero), Pinton 2012 ``15.31 Np/m`` (lineal). Medido: ``0.9524 / 0.9418 / 0.9191`` contra
    ``0.9722`` sin pérdidas — o sea **−2.0 / −3.1 / −5.5 %**.

    La atenuación se aplica con módulos de Lamé complejos y el mismo ``α`` para P y S (el material
    no lleva ``α``); las leyes son coeficientes *efectivos* de cráneo, así que esto mide la
    magnitud de la pérdida, no su mecanismo (D19).
    """
    base = _t()
    alphas, values = {}, {}
    for law in ("white_2006", "pichardo_2011", "pinton_2012"):
        scn = SC.Scenario(law=law)
        alphas[law] = scn.alpha()
        values[law] = _t(law=law)

    print("\n(3) atenuación (t = 2 mm, 100 kHz)")
    print(f"    sin atenuación   α = 0.00 Np/m   t_interior = {base:.6f}")
    for law in values:
        print(
            f"    {law:>15}  α = {alphas[law]:6.2f} Np/m   t_interior = {values[law]:.6f}   "
            f"Δ = {100 * (values[law] / base - 1):+.2f} %"
        )

    # (a) las tres reducen
    for law, v in values.items():
        assert v < base, f"{law} no reduce t_interior ({v:.6f} vs {base:.6f})"
    # (b) el orden sigue a α: más atenuación, menos presión media adentro
    by_alpha = sorted(values, key=lambda name: alphas[name])
    ordered = [values[name] for name in by_alpha]
    assert _monotone_decreasing(ordered), f"el orden no sigue a α: {by_alpha} → {ordered}"
    # (c) los α son los de las leyes del repo (no números inventados acá)
    assert alphas["pinton_2012"] == pytest.approx(15.31, abs=0.02)
    assert alphas["white_2006"] == pytest.approx(3.78, abs=0.02)
    assert 0.0 < alphas["white_2006"] < alphas["pichardo_2011"] < alphas["pinton_2012"]
    # (d) y el efecto es chico: no compite con el espesor (criterio 4)
    assert max(abs(v / base - 1.0) for v in values.values()) < 0.10


# --------------------------------------------------------------------------- #
# (4) ranking: el espesor manda
# --------------------------------------------------------------------------- #


@pytest.mark.slow
def test_thickness_dominates_ranking(tmp_path):
    """Criterio 4: en la grilla reducida, el parámetro más influyente es el **espesor**.

    ``analyse(include_grid=False)`` corre el OAT completo (con espesor cada 2 mm para abaratar) y
    devuelve el ranking por rango relativo ``(máx−mín)/base``. Medido en la grilla completa:
    espesor 72.1 %, ``c_p`` 16.1 %, atenuación 5.5 %, ``Z_r`` 3.2 %, ``c_s`` 2.2 %, ``ρ`` 1.2 %.

    El test **falla** si otro parámetro pasa al espesor: eso sería un hallazgo que hay que
    reportar, no un test a relajar (por eso no lleva tolerancia amplia: el margen medido es 4.5×).
    """
    report = SC.analyse(
        thickness_values=[1.0e-3, 3.0e-3, 5.0e-3, 7.0e-3],
        include_grid=False,
        verbose=False,
    )
    names = [r[0] for r in report.ranking]
    ranges = {r[0]: r[3] for r in report.ranking}

    print("\n(4) ranking OAT (grilla reducida, espesor 1/3/5/7 mm)")
    for name, lo, hi, rng in report.ranking:
        print(f"    {name:>12}  mín {lo:.6f}  máx {hi:.6f}  rango {100 * rng:6.2f} %")

    assert names[0] == "espesor", f"el más influyente es {names[0]!r}, no el espesor: {names}"
    assert ranges["espesor"] > 0.20, "el espesor tiene que mover t_interior más de 20 %"
    # margen contra el segundo: medido 4.5× en la grilla completa; acá se exige 2×
    assert ranges["espesor"] > 2.0 * ranges[names[1]]

    # el ranking se guarda, con las columnas que pide el encargo
    path = SC.write_csv(report, tmp_path / "sensitivity_cochlea.csv")
    with path.open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    print(f"    CSV: {path.name}, {len(rows)} filas")
    assert list(rows[0]) == ["parametro", "valor", "t_interior", "notas"]
    assert rows[-1]["parametro"] == "ranking"
    assert any(r["parametro"] == "espesor" for r in rows)
    assert all(r["notas"] is not None for r in rows)


def test_csv_schema(tmp_path):
    """El CSV tiene exactamente las cuatro columnas del encargo (sin ``solve`` de por medio)."""
    report = SC.Report(
        base={"t_interior": 0.5},
        rows=(SC.Row("espesor", "2.0 mm", 0.5, "nota; con punto y coma"),),
        grid=(SC.Row("espesor_x_atenuacion", "2.0 mm + pinton_2012", 0.48, "grilla"),),
        ranking=(("espesor", 0.4, 0.6, 0.4),),
    )
    path = SC.write_csv(report, tmp_path / "s.csv")
    with path.open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    assert [r["parametro"] for r in rows] == [
        "espesor",
        "espesor_x_atenuacion",
        "ranking",
    ]
    assert all(set(r) == {"parametro", "valor", "t_interior", "notas"} for r in rows)
    assert rows[0]["t_interior"] == "0.500000"
    assert "espesor=" in rows[2]["valor"]
