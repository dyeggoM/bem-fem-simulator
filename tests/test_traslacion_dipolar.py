"""La traslación dipolar contra la ley analítica `3/(2g+1)` en el límite de Rayleigh (D50).

Es el chequeo más fuerte que tiene el acoplamiento fluido-sólido: en `k₀a ≪ 1` la amplitud de
oscilación del cuerpo tiene forma cerrada (Bruus 2012, ecs. 22 y 26) y es **independiente de la
malla**, así que no hay nada que ajustar ni tolerancia que negociar.

De paso fija el defecto de la formulación total+restricción: no aproxima la traslación, la impone en
cero. Donde la física pide 0.18 (acero) devuelve 0.0000.
"""

from __future__ import annotations

import numpy as np
import pytest

from ultrasound_sim import materials as mats
from ultrasound_sim.fem import coupled3d
from ultrasound_sim.meshing3d_sphere import sphere_shell_tet_mesh

F = 100.0e3
AGUA = mats.WATER
KA = 0.1


@pytest.fixture(scope="module")
def escena():
    k0 = 2 * np.pi * F / float(AGUA.c)
    a = KA / k0
    mall = sphere_shell_tet_mesh(a, 4.0 * a, a / 4.5)
    return mall


def _traslacion(r, datos) -> np.ndarray:
    """Desplazamiento del centro de masa: proyección sobre el modo de traslación con la matriz de
    masa. Promediar nodos daría otro número porque la malla no es uniforme."""
    u = r["u"][datos["nodos_solidos"]].reshape(-1)
    mm = datos["M"]
    salida = []
    for d in range(3):
        t = np.zeros_like(u)
        t[d::3] = 1.0
        salida.append((t @ (mm @ u)) / (t @ (mm @ t)))
    return np.array(salida)


@pytest.mark.parametrize(
    "nombre, rho, cp, cs",
    [
        ("acero", 7800.0, 5900.0, 3200.0),
        ("hueso", 1950.0, 3500.0, 1900.0),
    ],
)
def test_la_traslacion_sigue_la_ley_3_sobre_2g_mas_1(escena, nombre, rho, cp, cs):
    """`|u_cm|/|u_i| = 3/(2g+1)` dentro del 5 %: acero g=7.8 → 0.1807, hueso g=1.95 → 0.6122.

    Se mide la formulación dispersada (la correcta) y la total sin restricción, que tienen que
    coincidir: el operador de contraste no puede sesgar la respuesta dipolar.
    """
    sol = mats.ElasticSolid(name=nombre, rho=rho, c_p=cp, c_s=cs, provenance="bench")
    g = rho / float(AGUA.rho)
    predicho = 3.0 / (2.0 * g + 1.0)

    r_dis = coupled3d.resolver(escena, AGUA, sol, F, abc_order=2, formulacion="dispersada")
    r_tot = coupled3d.resolver(escena, AGUA, sol, F, abc_order=2, formulacion="total",
                               fijar_rigidos=False)
    datos = r_dis["datos"]
    # referencia: el mismo proyector aplicado a la incidente, que viene como (n_solid, 3)
    u_i = r_dis["u_incidente"].reshape(-1)
    tt = np.zeros_like(u_i)
    tt[2::3] = 1.0
    ref = abs((tt @ (datos["M"] @ u_i)) / (tt @ (datos["M"] @ tt)))

    disp = abs(_traslacion(r_dis, datos)[2]) / ref
    tot = abs(_traslacion(r_tot, r_tot["datos"])[2]) / ref
    print(f"\n  {nombre}: predicho {predicho:.4f} · dispersada {disp:.4f} · total {tot:.4f}")
    assert abs(disp - predicho) / predicho < 0.05
    assert abs(tot - predicho) / predicho < 0.05
    assert abs(disp - tot) < 1e-3, "el operador de contraste sesgó la traslación"


def test_la_restriccion_anula_la_traslacion(escena):
    """El defecto queda fijado como número, no como comentario: con `fijar_rigidos=True` la
    traslación es cero donde la ley pide 0.18 (acero)."""
    sol = mats.ElasticSolid(name="acero", rho=7800.0, c_p=5900.0, c_s=3200.0, provenance="bench")
    r = coupled3d.resolver(escena, AGUA, sol, F, abc_order=2, formulacion="total",
                           fijar_rigidos=True)
    datos = r["datos"]
    u_i = r["u_incidente"].reshape(-1)
    tt = np.zeros_like(u_i)
    tt[2::3] = 1.0
    ref = abs((tt @ (datos["M"] @ u_i)) / (tt @ (datos["M"] @ tt)))
    con = abs(_traslacion(r, datos)[2]) / ref
    print(f"\n  acero con restricción: {con:.5f} (la ley pide 0.1807)")
    assert con < 1e-3
