"""Validación de la serie de Faran (1951) para la esfera elástica sumergida.

El problema de esta serie es que **no hay tabla de referencia a mano**, así que la validación tiene
que ser estructural. Las cinco que se usan:

1. **Condición unitaria por coeficiente**: para un dispersor sin pérdidas cada `Cₙ` cumple
   `Re Cₙ = −|Cₙ|²`. Se mide a **precisión de máquina**, y no se puede cumplir por casualidad: es una
   relación cuadrática entre la parte real y el módulo de cada coeficiente.
2. **Teorema óptico**: `σ_ext = σ_sca` con `σ_ext = (4π/k₀)Im A(0)`.
3. **Límite de fluido**: si la esfera es indistinguible del medio (`ρ` y `c_p` iguales, `c_s → 0`),
   la dispersión tiene que irse a cero, y hacerlo como `c_s²`.
4. **Aproximación al rígido**: cuanto más impedante es el sólido, más cerca tiene que quedar de la
   esfera rígida de `sphere.py` — y el orden entre materiales tiene que ser monótono.
5. **Residuos de las tres condiciones de borde** en `r = a`.

Motivo de cada umbral: los números salen de lo medido, no de una meta. La unidad de `Cₙ` no existe
(es un cociente de amplitudes), las secciones están en m² y los residuos son relativos.
"""

from __future__ import annotations

import numpy as np
import pytest

from ultrasound_sim import materials as mats
from ultrasound_sim.analytic import faran, sphere

F = 100.0e3
AGUA = mats.WATER
C0 = float(np.real(AGUA.c))
Z0 = AGUA.rho * C0

ACERO = mats.ElasticSolid(
    name="acero (benchmark)", rho=7800.0, c_p=5900.0, c_s=3200.0,
    provenance="benchmark paramétrico, no material de librería (D5)",
)
ALUMINIO = mats.ElasticSolid(
    name="aluminio (benchmark)", rho=2700.0, c_p=6320.0, c_s=3130.0,
    provenance="benchmark paramétrico, no material de librería (D5)",
)
TUNGSTENO = mats.ElasticSolid(
    name="tungsteno (benchmark)", rho=19300.0, c_p=5200.0, c_s=2800.0,
    provenance="benchmark paramétrico, no material de librería (D5)",
)
HUESO = mats.ElasticSolid(
    name="hueso (benchmark)", rho=1950.0, c_p=3500.0, c_s=1900.0,
    provenance="ballpark VERIFY de D20; la esfera elástica no tiene dato citable todavía",
)


def _radio(ka: float) -> float:
    return ka / abs(complex(AGUA.wavenumber(F)))


def _rigido(ka: float, n: int) -> np.ndarray:
    a = _radio(ka)
    return sphere.fluid_sphere_scattering(
        F, AGUA, AGUA, a, np.array([[0.0, 0.0, 1.5 * a]]), boundary="rigid"
    )["coefficients"][:n]


# --------------------------------------------------------------------------- 1. unitariedad


@pytest.mark.parametrize("ka", [0.5, 1.0, 2.0, 3.0])
def test_condicion_unitaria_a_precision_de_maquina(ka):
    """`Re Cₙ = −|Cₙ|²` para un dispersor sin pérdidas. Es el chequeo más fuerte: no depende de
    ninguna tabla, y una relación cuadrática entre la parte real y el módulo de cada coeficiente no
    se cumple por casualidad."""
    C = faran.coeficientes(F, AGUA, ACERO, _radio(ka))["C"]
    error = np.abs(np.real(C) + np.abs(C) ** 2).max()
    print(f"\n  ka={ka}: max |Re C + |C|²| = {error:.2e}")
    assert error < 1e-12


# --------------------------------------------------------------------------- 2. teorema óptico


@pytest.mark.parametrize("material", [ACERO, HUESO, ALUMINIO])
def test_el_teorema_optico_cierra(material):
    """`σ_ext = σ_sca`: conservación de energía. El desvío residual (~1e-5) es la cuadratura de
    `σ_sca`, no la serie."""
    a = _radio(1.0)
    sec = faran.secciones(faran.coeficientes(F, AGUA, material, a), a)
    print(f"\n  {material.name}: σ_sca {sec['sigma_sca']:.6e} m² · σ_ext {sec['sigma_ext']:.6e} m² "
          f"· desvío {100 * sec['desvio_relativo']:.3e} %")
    assert sec["desvio_relativo"] < 1e-4


# --------------------------------------------------------------------------- 3. límite de fluido


def test_solido_indistinguible_del_medio_no_dispersa():
    """`ρ` y `c_p` iguales al agua y `c_s → 0`: la dispersión tiene que caer como `c_s²`.

    Es el límite de contraste nulo, y sí es numéricamente usable: la singularidad está en `c_s → 0`
    *exacto*, pero la convergencia se ve antes de llegar. Ojo: la esfera tiene que ser de verdad
    indistinguible — con las mismas `ρ`/`c_p` pero `c_s` grande, la esfera es otro material y
    dispersa (medido: `|C₀| = 4.9e-2` con `c_s = 500`).
    """
    valores = []
    for c_s in (100.0, 50.0, 25.0):
        sol = mats.ElasticSolid(
            name=f"agua con c_s={c_s}", rho=AGUA.rho, c_p=C0, c_s=c_s, provenance="límite"
        )
        valores.append(np.abs(faran.coeficientes(F, AGUA, sol, _radio(1.0))["C"][0]))
    print(f"\n  |C₀| para c_s = 100/50/25: {[f'{v:.2e}' for v in valores]}")
    assert valores[0] < 5e-3
    for antes, despues in zip(valores, valores[1:]):
        assert despues < antes / 3.0, "no cae como c_s²: el potencial de corte no está entrando bien"


# --------------------------------------------------------------------------- 4. límite rígido


def test_cuanto_mas_impedante_mas_cerca_del_rigido():
    """Aproximación monótona a la esfera rígida con la impedancia `ρ·c_p`.

    Cruza contra otro camino de código (Anderson con contorno rígido). Los valores medidos a ka=1:
    aluminio (Z/Z₀ = 11.4) → 3.9e-3, acero (30.7) → 2.3e-3, tungsteno (66.9) → 1.2e-3. **No** son
    cero: un sólido elástico no es rígido, y el acero a ka=3 se aparta 2.9 %.
    """
    distancias = []
    for material in (ALUMINIO, ACERO, TUNGSTENO):
        C = faran.coeficientes(F, AGUA, material, _radio(1.0))["C"]
        d = np.abs(C - _rigido(1.0, len(C))).max()
        distancias.append(d)
        print(f"\n  {material.name:22s} Z/Z₀ = {material.rho * material.c_p / Z0:5.1f} → "
              f"|C − C_rígido| = {d:.3e}")
    assert all(a > b for a, b in zip(distancias, distancias[1:])), "el orden por impedancia no es monótono"
    assert distancias[2] < 2e-3


def test_el_acero_coincide_con_el_rigido_en_el_limite_de_ka_chico():
    """A `ka → 0` la longitud de onda es mucho mayor que la esfera: todo sólido se ve rígido."""
    C = faran.coeficientes(F, AGUA, ACERO, _radio(0.2))["C"]
    d = np.abs(C - _rigido(0.2, len(C))).max()
    print(f"\n  ka=0.2: |C_acero − C_rígido| = {d:.3e}")
    assert d < 1e-3


# --------------------------------------------------------------------------- 5. condiciones de borde


@pytest.mark.parametrize("material", [ACERO, HUESO])
def test_las_tres_condiciones_de_borde_se_satisfacen(material):
    """Desplazamiento normal, tensión normal y tensión tangencial nula, después de resolver."""
    d = faran.residuos(F, AGUA, material, _radio(1.0))
    print(f"\n  {material.name}: u_r {d['u_r']:.1e} · σ_rr {d['sigma_rr']:.1e} · σ_rθ {d['sigma_rt']:.1e}")
    assert d["u_r"] < 1e-9
    assert d["sigma_rr"] < 1e-9
    assert d["sigma_rt"] < 1e-9


# --------------------------------------------------------------------------- 6. física


def test_la_serie_elastica_no_es_la_fluida():
    """Con las mismas `c_p`/`ρ`, el sólido tiene que dispersar distinto que el fluido: si diera
    exactamente lo mismo, las ondas de corte no estarían entrando en el sistema."""
    a = _radio(1.0)
    fluido = mats.Material(name="hueso como fluido", c=HUESO.c_p, rho=HUESO.rho, provenance="test")
    amp_elast = faran.amplitud_lejana(faran.coeficientes(F, AGUA, HUESO, a), [np.pi], a)[0]
    p_fluido = sphere.fluid_sphere_scattering(
        F, AGUA, fluido, a, np.array([[0.0, 0.0, -1.5 * a]])
    )["scattered"][0]
    print(f"\n  |A(π)| elástico {abs(amp_elast):.4e} · esfera fluida |p_s| {abs(p_fluido):.4e}")
    assert abs(abs(amp_elast) - abs(p_fluido)) > 1e-4


def test_validacion_de_argumentos():
    with pytest.raises(ValueError):
        faran.coeficientes(0.0, AGUA, ACERO, 1e-3)
    with pytest.raises(ValueError):
        faran.coeficientes(F, AGUA, ACERO, -1e-3)
    with pytest.raises(TypeError):
        faran.coeficientes(F, AGUA, AGUA, 1e-3)   # la esfera tiene que ser ElasticSolid

# --------------------------------------------------------------------------- campo interior


def test_el_campo_interior_cumple_la_condicion_de_continuidad():
    """Validación del desplazamiento **dentro** del sólido contra la condición de borde.

    `campo_interior` se valida contra lo que ya está validado: en `r = a` la componente radial tiene
    que igualar a la del fluido, `(1/(ρ₀ω²))·∂p/∂r`. Si los potenciales interiores estuvieran mal
    escalados —o con un orden de Bessel corrido— esto no cierra.

    Hace falta porque en el eje `Pₙ(1) = 1` y `dPₙ/dθ = 0`: mirar sólo ahí **esconde todos los
    órdenes `n ≥ 1`**, que es justo donde aparece la deformación frente a la traslación.
    """
    a = _radio(1.0)
    coef = faran.coeficientes(F, AGUA, ACERO, a)
    k0 = coef["k0"]
    omega = 2 * np.pi * F
    th = np.radians(np.array([20.0, 50.0, 90.0, 130.0, 160.0]))
    dirs = np.stack([np.sin(th), np.zeros_like(th), np.cos(th)], axis=1)

    def punto(rr):
        return rr * a * dirs

    u = faran.campo_interior(coef, punto(1.0))
    u_r_solido = np.einsum("ij,ij->i", u, dirs)

    dr = 1e-6 * a
    p_abajo = faran.campo_dispersado(coef, punto(1 - 1e-6), a) + np.exp(1j * k0 * punto(1 - 1e-6)[:, 2])
    p_arriba = faran.campo_dispersado(coef, punto(1 + 1e-6), a) + np.exp(1j * k0 * punto(1 + 1e-6)[:, 2])
    u_r_fluido = (p_arriba - p_abajo) / (2 * dr) / (AGUA.rho * omega**2)

    razon = np.abs(u_r_solido / u_r_fluido)
    print(f"\n  razón u_r(sólido)/u_r(fluido) en r=a: {np.round(razon, 6)}")
    assert np.allclose(razon, 1.0, atol=1e-4)


def test_el_interior_tiene_deformacion_y_no_solo_traslacion():
    """El campo interior de un sólido elástico no puede ser un movimiento de cuerpo rígido: si
    quitando la parte rígida queda ~0, los potenciales P y S no están entrando."""
    a = _radio(1.0)
    coef = faran.coeficientes(F, AGUA, ACERO, a)
    th = np.radians(np.array([30.0, 60.0, 90.0, 120.0, 150.0]))
    pts = np.stack([(0.7 * a) * np.sin(th), np.zeros_like(th), (0.7 * a) * np.cos(th)], axis=1)
    u = faran.campo_interior(coef, pts)

    x = pts - pts.mean(axis=0)
    R = np.zeros((3 * len(pts), 6))
    for k in range(3):
        R[k::3, k] = 1.0
    R[1::3, 3], R[2::3, 3] = -x[:, 2], x[:, 1]
    R[0::3, 4], R[2::3, 4] = x[:, 2], -x[:, 0]
    R[0::3, 5], R[1::3, 5] = -x[:, 1], x[:, 0]
    coef_r, *_ = np.linalg.lstsq(R, u.reshape(-1), rcond=None)
    deformacion = np.abs(u - (R @ coef_r).reshape(-1, 3)).mean()
    print(f"\n  |u| medio {np.abs(u).mean():.3e} · deformación {deformacion:.3e} "
          f"({100*deformacion/np.abs(u).mean():.0f} % del total)")
    assert deformacion > 0.2 * np.abs(u).mean()
