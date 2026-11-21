"""El arreglo y las fuentes generadas tienen que dar EL MISMO campo.

Es la pregunta que aparece en cuanto se usan las dos pestañas: «¿estamos usando los mismos elementos?».
Antes de este test la respuesta no estaba fijada en ningún lado, y las dos rutas podían divergir en
silencio —de hecho divergían, no en la fórmula sino en el estado del cliente (D64)—.

La ruta del servidor (`n_elements` + `aperture_mm`) y la de fuentes explícitas tienen que coincidir al
bit: es el mismo problema escrito de dos formas.
"""

from __future__ import annotations

import numpy as np

from ultrasound_sim.dashboard.server import calcular

GRILLA = {"nx": 100, "ny": 80, "x_max_mm": 150, "y_max_mm": 60}
BASE = {
    "material": "agua", "frequency_khz": 250.0, "modelo": "baffle", "grid": GRILLA,
    "n_elements": 16, "aperture_mm": 24.0, "focus_x_mm": 90.0, "focus_y_mm": 0.0,
}


def _fuentes_como_el_cliente(n, apertura_mm, fx_mm, fy_mm, f_khz, c=1500.0):
    """Replica `generarDesdeArreglo` del cliente: paso = D/n, y el foco en radianes."""
    k = 2 * np.pi * (f_khz * 1e3) / c
    paso = apertura_mm / n
    salida = []
    for j in range(n):
        y = -apertura_mm / 2 + paso * (j + 0.5)
        fase = -k * float(np.hypot(fx_mm, fy_mm - y)) / 1e3
        salida.append({"y": y * 1e-3, "width": paso * 1e-3, "amp": 1.0, "phase": fase})
    return salida


def test_el_arreglo_y_sus_fuentes_dan_el_mismo_campo():
    a = calcular(dict(BASE))
    b = calcular({**BASE, "n_elements": 0,
                  "sources": _fuentes_como_el_cliente(16, 24.0, 90.0, 0.0, 250.0)})
    da = np.array(a["p2_db"])
    db = np.array(b["p2_db"])
    dif = float(np.max(np.abs(da - db)))
    print(f"\n  diferencia máxima entre arreglo y fuentes: {dif:.3e} dB · "
          f"pico {a['metrics']['p_max']:.6f} vs {b['metrics']['p_max']:.6f}")
    assert dif < 1e-6, "el arreglo y sus fuentes no describen el mismo problema"
    assert abs(a["metrics"]["p_max"] - b["metrics"]["p_max"]) / a["metrics"]["p_max"] < 1e-9


def test_cambiar_el_foco_cambia_las_fases_y_el_campo():
    """El foco entra en la fase, no en la posición: mover el foco tiene que mover el pico."""
    cerca = calcular({**BASE, "focus_x_mm": 55.0})
    lejos = calcular({**BASE, "focus_x_mm": 90.0})
    f_cerca = cerca["metrics"]["peak_x_mm"]
    f_lejos = lejos["metrics"]["peak_x_mm"]
    print(f"\n  foco 55 mm → pico en {f_cerca:.1f} mm · foco 90 mm → pico en {f_lejos:.1f} mm")
    assert f_cerca < f_lejos, "acercar el foco tiene que acercar el máximo"


def test_la_frecuencia_cambia_las_fases_del_mismo_arreglo():
    a = calcular({**BASE, "frequency_khz": 250.0})
    b = calcular({**BASE, "frequency_khz": 400.0})
    fa = _fuentes_como_el_cliente(16, 24.0, 90.0, 0.0, 250.0)
    fb = _fuentes_como_el_cliente(16, 24.0, 90.0, 0.0, 400.0)
    diferencia_fases = max(abs(x["phase"] - y["phase"]) for x, y in zip(fa, fb))
    print(f"\n  250 kHz vs 400 kHz: la misma fuente cambia {np.degrees(diferencia_fases):.1f}°")
    assert diferencia_fases > 0.1, "la frecuencia tiene que entrar en la fase de foco"
    assert a["metrics"]["peak_x_mm"] != b["metrics"]["peak_x_mm"]
