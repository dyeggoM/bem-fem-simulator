"""Atenuación en el sólido por módulos complejos: que la librería y el script den lo mismo.

El script de sensibilidad (``scripts/sensitivity_cochlea.py``) tiene su **propio** ensamblado acoplado
con ``D`` compleja desde antes de que la librería supiera de ``alpha``. Ahora ``solve_cross_section``
acepta ``alpha`` [Np/m], así que hay dos implementaciones independientes del mismo damping y este
archivo las obliga a coincidir: si una cambia un signo o una convención, el test lo ve.

Los dos criterios
-----------------
1. ``alpha = 0`` no cambia nada: mismo camino de código, valor idéntico al de no pasarlo.
2. Con la ley de Pinton (``α = 15.31 Np/m`` a 100 kHz, D19) la librería y el script coinciden dentro
   de ``1e-6`` relativo, con el mismo mallado y las mismas condiciones.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

from ultrasound_sim import materials
from ultrasound_sim.cochlea import build_cross_section, solve_cross_section

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "sensitivity_cochlea.py"


def _load_script():
    name = "sensitivity_cochlea"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, SCRIPT)
    if spec is None or spec.loader is None:
        raise ImportError(f"no se pudo cargar {SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


SC = _load_script()


def _seccion(scn, ppw: float = 20.0):
    """El mismo ``build_cross_section`` que usa ``SC.evaluate``, para comparar mallados idénticos."""
    return build_cross_section(
        SC.FREQUENCY_HZ,
        SC.RADIUS_M,
        scn.thickness_m,
        ppw=ppw,
        pad_wavelengths=SC.PAD_WAVELENGTHS,
        min_elements_across_wall=SC.MIN_ELEMENTS_ACROSS_WALL,
        bone=scn.bone(),
        window_half_angle=SC.WINDOW_HALF_ANGLE_DEG,
        window_center_angle=SC.WINDOW_CENTER_ANGLE_DEG,
        window_model="impedance",
        window_impedance=scn.impedance(),
    )


def test_alpha_cero_no_cambia_nada():
    """``alpha = 0`` tiene que dar exactamente lo mismo que no pasarlo (misma rama, no redondeo)."""
    scn = SC.Scenario()
    sec = _seccion(scn)
    sin = solve_cross_section(sec)
    con = solve_cross_section(sec, alpha=0.0)
    print(f"\n  t_interior sin alpha = {sin.t_interior:.12f} | con alpha=0 = {con.t_interior:.12f}")
    assert con.t_interior == sin.t_interior
    assert con.t_interior_complex == sin.t_interior_complex


def test_la_libreria_reproduce_el_ensamblado_del_script_con_atenuacion():
    """Dos ensamblados independientes, la misma ley de atenuación: tienen que coincidir."""
    scn = SC.Scenario(law="pinton_2012")
    alpha = float(scn.alpha())
    assert alpha > 0.0, "la ley de Pinton tiene que dar α > 0 a 100 kHz"

    del_script = SC.evaluate(scn, ppw=20.0)["t_interior"]
    sec = _seccion(scn)
    de_la_libreria = solve_cross_section(sec, alpha=alpha).t_interior

    rel = abs(de_la_libreria - del_script) / abs(del_script)
    print(f"\n  α = {alpha:.2f} Np/m   script = {del_script:.6f}   librería = {de_la_libreria:.6f}   "
          f"Δ rel = {rel:.2e}")
    assert rel < 1e-6, f"los dos ensamblados discrepan: {rel:.2e}"


def test_la_atenuacion_reduce_la_transmision_y_es_monotona():
    """Más α, menos pasa: la atenuación no puede aumentar la transmisión."""
    scn = SC.Scenario()
    sec = _seccion(scn)
    valores = [solve_cross_section(sec, alpha=a).t_interior for a in (0.0, 5.0, 15.31, 40.0)]
    print("\n  α [Np/m] = 0 / 5 / 15.31 / 40 -> t_interior = "
          + " / ".join(f"{v:.5f}" for v in valores))
    assert all(b < a for a, b in zip(valores, valores[1:])), f"no decrece: {valores}"


def test_alpha_invalido_se_rechaza():
    sec = _seccion(SC.Scenario())
    with pytest.raises(ValueError, match="alpha"):
        solve_cross_section(sec, alpha=-1.0)
