"""El BEM de referencia está validado contra la solución exacta (D54).

Este test fija la pieza sobre la que se apoya todo el contraste de métodos: el campo que muestra la
vista `campo` tiene que ser la solución exacta de un monopolo 2D cuando la fuente es mucho más chica
que la longitud de onda. Si esto se rompe, cualquier "el FEM no coincide" deja de tener sentido —
primero hay que saber cuál de los dos es la referencia.

Lo que se valida es la **forma**: que `p_BEM / H₀(k·r)` sea el mismo número complejo en todas las
distancias. Un factor complejo global (la escala física y el signo de convención) no se juzga, porque
no cambia ni el decaimiento ni la fase relativa, que es lo que hace que dos métodos "resuelvan
igual".
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.special import hankel1

from ultrasound_sim.dashboard import field
from ultrasound_sim.materials import WATER

F = 250e3
RADIOS_MM = (30.0, 60.0, 100.0, 140.0)


def test_el_campo_de_una_fuente_chica_es_la_hankel_en_todo_el_rango():
    k = complex(WATER.wavenumber(F))
    # fuente de 0.5 mm = λ/12
    fuentes = [field.Source(y=0.0, width=0.0005, amp=1.0, phase=0.0)]

    cocientes = []
    for r_mm in RADIOS_MM:
        r = r_mm * 1e-3
        y = -0.6 * r                      # fuera del eje: no se apoya en el borde del baffle
        x = np.sqrt(r**2 - y**2)
        p_bem = complex(field.evaluate(np.array([x]), np.array([y]), fuentes, F, WATER)[0, 0])
        cocientes.append(p_bem / hankel1(0, k * r))

    cocientes = np.array(cocientes)
    normalizados = cocientes / cocientes[0]
    modulos = np.abs(normalizados)
    fases = np.degrees(np.angle(normalizados))
    print("\n  |cociente|/|cociente(30mm)|:", np.round(modulos, 6))
    print("  fase relativa (°):", np.round(fases, 4))
    print(f"  factor global: |{abs(cocientes[0]):.1f}| · {np.degrees(np.angle(cocientes[0])):+.1f}°")

    assert np.allclose(modulos, 1.0, atol=2e-4), "el decaimiento no es el de H₀(k·r)"
    assert np.allclose(fases, 0.0, atol=0.05), "la fase no es la de H₀(k·r)"


def test_el_ancho_de_la_fuente_entra_en_el_campo():
    """Contraprueba: si `evaluate` ignorara el ancho, el test de arriba pasaría igual. Se compara la
    MISMA posición con dos anchos muy distintos (λ/12 contra 3λ): el campo tiene que cambiar."""
    lam = WATER.wavelength(F)
    x, y = 0.03, 0.0
    fina = [field.Source(y=0.0, width=lam / 12, amp=1.0, phase=0.0)]
    ancha = [field.Source(y=0.0, width=3 * lam, amp=1.0, phase=0.0)]
    p_fina = complex(field.evaluate(np.array([x]), np.array([y]), fina, F, WATER)[0, 0])
    p_ancha = complex(field.evaluate(np.array([x]), np.array([y]), ancha, F, WATER)[0, 0])
    cambio = abs(p_ancha / p_fina - 1.0)
    print(f"\n  mismo punto, ancho λ/12 vs 3λ: |cambio| = {100 * cambio:.1f} %")
    assert cambio > 0.05, "el ancho de la fuente no está entrando en el cálculo"
