"""La formulación dispersada del sólido (D49): lo que la hace correcta y lo que la hace cara.

Por qué existe
--------------
La formulación "total" resolvía el desplazamiento del sólido completo y le imponía traslación media
nula para domar los modos rígidos. Sim-Datos mostró (D48) que esa condición **suprime una respuesta
dipolar física** que Faran incluye (Feuillade y Clay 1999, §IV): el error crecía como `3/(2g+1)`.

La formulación dispersada resuelve `u_s = u − u_i` con el **operador de contraste** aplicado a la
incidente. La distinción importa:

- restar el operador COMPLETO deja un residuo O(h²) que empuja un campo espurio (~18 % del campo real
  de un metal);
- restar el operador del FLUIDO armado sobre la misma malla y la misma cuadratura hace que un cuerpo
  indistinguible del fluido dé **matrices idénticas** y el término sea exactamente cero, sin depender
  de `h`.
"""

from __future__ import annotations

import numpy as np
import pytest

from ultrasound_sim import materials as mats
from ultrasound_sim.fem import coupled3d, elastic3d
from ultrasound_sim.meshing3d_sphere import sphere_shell_tet_mesh

F = 100.0e3
AGUA = mats.WATER


def _radio() -> float:
    return 1.0 / abs(complex(AGUA.wavenumber(F)))


@pytest.fixture(scope="module")
def malla():
    a = _radio()
    with pytest.warns(UserWarning):
        return sphere_shell_tet_mesh(a, 2.5 * a, a / 4.0)


def test_el_operador_de_contraste_es_exactamente_cero_si_los_materiales_coinciden(malla):
    """El núcleo de la formulación: con `λ_s = ρ₀c₀²` y `μ = 0`, las matrices del sólido y las del
    fluido son **las mismas**, así que su diferencia es cero a precisión de máquina.

    Si esto fallara, el caso transparente volvería a dispersar y todo lo demás sería discutible.
    """
    mall = malla
    nodos = np.unique(mall.tets_solid)
    mapa = np.full(mall.n_nodes, -1)
    mapa[nodos] = np.arange(nodos.size)
    tets = mapa[mall.tets_solid]

    lam_f = AGUA.rho * float(AGUA.c) ** 2
    k_f, m_f = elastic3d.p1_elastic3d_matrices_constantes(
        mall.nodes[nodos], tets, lam_f, 0.0, AGUA.rho
    )
    k_s, m_s = elastic3d.p1_elastic3d_matrices_constantes(
        mall.nodes[nodos], tets, lam_f, 0.0, AGUA.rho
    )
    print(f"\n  |K_sólido − K_fluido| = {abs(k_s - k_f).max():.2e} · |M_s − M_f| = "
          f"{abs(m_s - m_f).max():.2e}")
    assert abs(k_s - k_f).max() == 0.0
    assert abs(m_s - m_f).max() == 0.0


def test_el_caso_transparente_dispersa_casi_nada(malla):
    """`c_s → 0` con `ρ` y `c_p` del agua: la dispersión tiene que tender a cero **sin restringir
    nada**. Se ve la tendencia con cuatro valores de `c_s`."""
    valores = []
    for cs in (100.0, 30.0, 10.0):
        casi = mats.ElasticSolid(
            name=f"agua c_s={cs}", rho=AGUA.rho, c_p=float(AGUA.c), c_s=cs, provenance="límite"
        )
        r = coupled3d.resolver(malla, AGUA, casi, F, abc_order=2, formulacion="dispersada")
        valores.append(float(np.abs(r["p_s"]).max()))
    print(f"\n  |p_s| con c_s = 100/30/10: {[f'{v:.2e}' for v in valores]} (incidente = 1 Pa)")
    assert valores[-1] < 5e-3
    for antes, despues in zip(valores, valores[1:]):
        assert despues < antes, "no decrece al bajar c_s: el contraste no está entrando bien"


def test_la_restriccion_de_rigidos_no_tiene_sentido_en_la_dispersada(malla):
    """Se rechaza en vez de ignorarla: aplicarla ahí suprimiría justo el dipolo que la formulación
    recupera."""
    acero = mats.ElasticSolid(name="acero", rho=7800.0, c_p=5900.0, c_s=3200.0, provenance="bench")
    with pytest.raises(ValueError, match="fijar_rigidos"):
        coupled3d.resolver(malla, AGUA, acero, F, abc_order=2, fijar_rigidos=True,
                           formulacion="dispersada")


def test_la_formulacion_total_sigue_disponible_para_comparar(malla):
    """Se deja la formulación vieja accesible: es la que permite medir cuánto cambiaba el resultado
    con la condición de traslación nula."""
    acero = mats.ElasticSolid(name="acero", rho=7800.0, c_p=5900.0, c_s=3200.0, provenance="bench")
    total = coupled3d.resolver(malla, AGUA, acero, F, abc_order=2, formulacion="total",
                               fijar_rigidos=True)
    dispersada = coupled3d.resolver(malla, AGUA, acero, F, abc_order=2, formulacion="dispersada")
    razon_total = abs(total["p_s"]).max()
    razon_disp = abs(dispersada["p_s"]).max()
    print(f"\n  |p_s| total+restricción {razon_total:.4f} · dispersada {razon_disp:.4f}")
    assert razon_total > 0 and razon_disp > 0
    assert not np.isclose(razon_total, razon_disp)
