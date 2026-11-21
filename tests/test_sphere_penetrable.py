"""Validación de la esfera fluida de Anderson (1950) y del dispersado penetrable 3D.

Anderson: *Sound Scattering from a Fluid Sphere*, JASA 22:426 — solución por ondas parciales (la
"Mie acústica"), válida para todo `ka` (la restricción `ka ≪ 1` es de Rayleigh, no de esta serie).

Estos tests fijan la referencia ANTES de comparar el FEM contra ella, porque una referencia sin
validar es peor que ninguna: si la serie está mal, el FEM "coincide" con un error.
"""

from __future__ import annotations

import numpy as np
import pytest

from ultrasound_sim import materials as mats
from ultrasound_sim.analytic import sphere
from ultrasound_sim.fem import helmholtz3d as h3
from ultrasound_sim.meshing3d import boundary_faces, box_tet_mesh_from_spacing
from scipy.interpolate import LinearNDInterpolator

F = 100.0e3
WATER = mats.WATER
HUESO_FLUIDO = mats.Material(
    name="hueso como fluido (test)", c=3500.0, rho=1950.0,
    provenance="benchmark de dispersión; no es un material de librería",
)


def _esfera_de_puntos(radio: float, n_theta: int = 48, n_phi: int = 20) -> np.ndarray:
    th = np.linspace(0.0, np.pi, n_theta)
    ph = np.linspace(0.0, 2.0 * np.pi, n_phi, endpoint=False)
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    return np.column_stack([
        (radio * np.sin(TH) * np.cos(PH)).ravel(),
        (radio * np.sin(TH) * np.sin(PH)).ravel(),
        (radio * np.cos(TH)).ravel(),
    ])


def _radio_para(ka: float) -> float:
    return ka / abs(WATER.wavenumber(F))


# --------------------------------------------------------------------------- coherencia de la serie


def test_esfera_casada_no_dispersa():
    """ρ y c iguales al medio: Aₙ = 0. Si no, hay un error de signo o de factor en la fórmula."""
    a = _radio_para(1.0)
    r = sphere.fluid_sphere_scattering(F, WATER, WATER, a, _esfera_de_puntos(1.5 * a))
    print(f"\n  |Aₙ| máximo = {np.abs(r['coefficients']).max():.2e}")
    assert np.abs(r["coefficients"]).max() < 1e-12
    assert np.abs(r["scattered"]).max() < 1e-12


def test_el_limite_de_densidad_infinita_es_el_rigido():
    """ρ₁ → ∞ tiene que dar la condición rígida (∂ᵣp = 0), no la blanda (p = 0).

    Fue el chequeo que encontró el cociente de densidades invertido: con ρ₁/ρ₀ el límite daba
    blando, y como la presión sigue siendo continua por construcción, el error sólo se veía en el
    salto de flujo (73.7 % constante en ka).
    """
    a = _radio_para(1.0)
    pts = _esfera_de_puntos(1.5 * a)
    rigido = sphere.fluid_sphere_scattering(F, WATER, WATER, a, pts, boundary="rigid")
    denso = mats.Material(name="casi rígido", c=1500.0, rho=1e9, provenance="límite ρ→∞")
    limite = sphere.fluid_sphere_scattering(F, WATER, denso, a, pts)
    err = np.linalg.norm(limite["scattered"] - rigido["scattered"]) / np.linalg.norm(rigido["scattered"])
    print(f"\n  ρ→∞ vs rígido explícito: {100 * err:.4f} %")
    # 1e-4 y no 1e-6: el límite es numérico (ρ = 1e9), no algebraico
    assert err < 1e-4


def test_la_derivada_normal_se_anula_en_la_superficie_rigida():
    a = _radio_para(1.0)
    r = sphere.fluid_sphere_scattering(F, WATER, WATER, a, _esfera_de_puntos(a), boundary="rigid")
    ratio = np.abs(r["d_total"]).max() / np.abs(r["total"]).max()
    print(f"\n  |∂ᵣp| / |p| en r = a: {ratio:.2e} 1/m")
    assert ratio < 1e-9


@pytest.mark.parametrize("ka", [0.5, 1.0, 3.0, 5.0])
def test_continuidad_de_presion_y_flujo_en_la_superficie(ka):
    """Las dos condiciones de la interfaz, evaluadas sobre r = a con los tres campos."""
    a = _radio_para(ka)
    r = sphere.fluid_sphere_scattering(F, WATER, HUESO_FLUIDO, a, _esfera_de_puntos(a))
    salto_p = np.abs(r["total"] - r["interior"]).max() / np.abs(r["total"]).max()
    salto_flujo = (np.abs(r["d_total"] / r["rho_out"] - r["d_interior"] / r["rho_in"]).max()
                   / np.abs(r["d_total"] / r["rho_out"]).max())
    print(f"\n  ka={r['ka']:.2f}: salto de p {100 * salto_p:.1e} % · salto de flujo "
          f"{100 * salto_flujo:.1e} %")
    assert salto_p < 1e-9
    assert salto_flujo < 1e-9


def test_puntos_dentro_de_la_esfera_se_rechazan():
    a = _radio_para(1.0)
    adentro = np.array([[0.0, 0.0, 0.5 * a]])
    with pytest.raises(ValueError) as exc:
        sphere.fluid_sphere_scattering(F, WATER, HUESO_FLUIDO, a, adentro)
    assert "dentro de la esfera" in str(exc.value)


def test_boundary_desconocida():
    a = _radio_para(1.0)
    with pytest.raises(ValueError):
        sphere.fluid_sphere_scattering(F, WATER, HUESO_FLUIDO, a, _esfera_de_puntos(2 * a),
                                       boundary="elastica")


def test_radio_y_frecuencia_invalidos():
    with pytest.raises(ValueError):
        sphere.fluid_sphere_scattering(0.0, WATER, HUESO_FLUIDO, 1e-3, [[0, 0, 1]])
    with pytest.raises(ValueError):
        sphere.fluid_sphere_scattering(F, WATER, HUESO_FLUIDO, -1e-3, [[0, 0, 1]])


# --------------------------------------------------------------------------- el FEM contra Anderson


def _fem_dispersado(ka: float, ppw: float):
    """FEM 3D del dispersado con el término de contraste; devuelve (malla, total, a)."""
    k0 = WATER.wavenumber(F)
    a = ka / abs(k0)
    lam = abs(WATER.wavelength(F))
    R = 3.0 * a
    malla = box_tet_mesh_from_spacing((-R,) * 3, (R,) * 3, lam / ppw)
    centro = malla.nodes[malla.tets].mean(axis=1)
    dentro = np.linalg.norm(centro, axis=1) < a
    k_el = np.where(dentro, HUESO_FLUIDO.wavenumber(F), k0)
    rho_el = np.where(dentro, HUESO_FLUIDO.rho, WATER.rho)
    p_inc = h3.plane_wave(F, WATER, (0.0, 0.0, 1.0))
    p_s = h3.solve_scattered(
        malla.nodes, malla.tets, k_el, abc_faces=boundary_faces(malla.tets),
        element_rho=rho_el, background_k=k0, background_rho=WATER.rho, p_inc=p_inc,
    )
    return malla, h3.total_field(p_s, p_inc, malla.nodes), a, int(dentro.sum())


def _error_contra_anderson(malla, total, a, ka: float) -> float:
    pts = _esfera_de_puntos(1.5 * a)
    serie = sphere.fluid_sphere_scattering(F, WATER, HUESO_FLUIDO, a, pts)
    interp = LinearNDInterpolator(
        malla.nodes, np.column_stack([np.real(total), np.imag(total)]), fill_value=np.nan
    )
    fem = interp(pts)
    ok = ~np.isnan(fem[:, 0])
    valores = fem[ok, 0] + 1j * fem[ok, 1]
    return float(np.linalg.norm(valores - serie["total"][ok]) / np.linalg.norm(serie["total"][ok]))


@pytest.mark.parametrize("ka", [1.0, 2.0])
def test_fem_penetrable_3d_contra_anderson(ka):
    """El criterio: el campo total del FEM dentro del 5 % de la serie sobre r = 1.5a."""
    malla, total, a, n_dentro = _fem_dispersado(ka, ppw=12)
    err = _error_contra_anderson(malla, total, a, ka)
    print(f"\n  ka={ka}: {malla.n_nodes} nodos, {n_dentro} elementos en la esfera, error {100 * err:.2f} %")
    assert n_dentro > 0, "ningún elemento cayó dentro de la esfera: la malla no la representa"
    assert err < 0.05, f"FEM fuera del 5 % de Anderson: {100 * err:.1f} %"


def test_convergencia_del_fem_penetrable():
    """Refinar tiene que acercar a la serie: es lo que distingue una solución de una coincidencia."""
    errores = []
    for ppw in (8, 16):
        malla, total, a, _ = _fem_dispersado(1.0, ppw=ppw)
        errores.append(_error_contra_anderson(malla, total, a, 1.0))
    print(f"\n  ppw 8: {100 * errores[0]:.2f} % · ppw 16: {100 * errores[1]:.2f} %")
    assert errores[1] < errores[0], f"no converge al refinar: {errores}"
    assert errores[1] < 0.03


def test_sin_contraste_el_dispersado_es_nulo():
    """Esfera casada con el medio: no hay fuente de contraste, así que p_s = 0 (y no "casi cero")."""
    k0 = WATER.wavenumber(F)
    lam = abs(WATER.wavelength(F))
    malla = box_tet_mesh_from_spacing((-lam / 2,) * 3, (lam / 2,) * 3, lam / 10)
    p_inc = h3.plane_wave(F, WATER, (0.0, 0.0, 1.0))
    p_s = h3.solve_scattered(
        malla.nodes, malla.tets, np.full(malla.tets.shape[0], k0),
        abc_faces=boundary_faces(malla.tets), element_rho=WATER.rho,
        background_k=k0, background_rho=WATER.rho, p_inc=p_inc,
    )
    assert np.abs(p_s).max() < 1e-9
