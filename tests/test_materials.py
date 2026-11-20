"""Materiales: convención de signos, unidades y procedencia.

Fija D4/D5 de docs/decisions.md. Si estos tests cambian, cambia la física del repo entero.
"""

import math

import numpy as np
import pytest

from ultrasound_sim import materials as mat


def test_water_values_come_from_the_plan():
    assert mat.WATER.c == 1500.0
    assert mat.WATER.rho == 1000.0
    assert mat.WATER.impedance == pytest.approx(1.5e6)
    assert mat.WATER.provenance.startswith("plan:")


def test_fluids_are_flagged_when_the_plan_has_no_numbers():
    # El plan no fija c, rho ni alpha de tejido blando / hueso: deben estar marcados, no silenciosos.
    assert mat.SOFT_TISSUE.name in mat.NEEDS_VERIFICATION
    assert mat.BONE_CORTICAL.name in mat.NEEDS_VERIFICATION
    assert "VERIFY" in mat.SOFT_TISSUE.provenance
    assert "ballpark" in mat.BONE_CORTICAL.provenance
    # Y los que sí vienen del plan no están marcados.
    assert mat.WATER.name not in mat.NEEDS_VERIFICATION


def test_attenuation_free_medium_has_real_wavenumber():
    f = 1.0e6
    k = mat.WATER.wavenumber(f)
    assert k.imag == pytest.approx(0.0)
    assert k.real == pytest.approx(2.0 * math.pi * f / 1500.0)


def test_power_law_attenuation_and_sign_convention():
    m = mat.Material(name="t", c=1500.0, rho=1000.0, alpha0=10.0, alpha_power=1.0, f_ref=1.0e6)
    assert m.attenuation(1.0e6) == pytest.approx(10.0)
    assert m.attenuation(2.0e6) == pytest.approx(20.0)

    n = mat.Material(name="n", c=1500.0, rho=1000.0, alpha0=10.0, alpha_power=2.0, f_ref=1.0e6)
    assert n.attenuation(2.0e6) == pytest.approx(40.0)

    k = m.wavenumber(1.0e6)
    assert k.imag == pytest.approx(10.0), "parte imaginaria positiva: exp(+ikx) debe decaer"


def test_plane_wave_decays_with_distance_under_attenuation():
    """exp(i k x) con k = kr + i*alpha decae como exp(-alpha x): la convención de D4."""
    m = mat.Material(name="t", c=1500.0, rho=1000.0, alpha0=5.0, f_ref=1.0e6)
    k = m.wavenumber(1.0e6)
    x = 0.7
    ratio = abs(np.exp(1j * k * x))
    assert ratio == pytest.approx(math.exp(-5.0 * x), rel=1e-12)
    assert ratio < 1.0


def test_wavelength_matches_the_plan_numbers():
    # Del plan: 75 mm @ 20 kHz, 15 mm @ 100 kHz, 1.5 mm @ 1 MHz en agua.
    assert mat.WATER.wavelength(20e3) == pytest.approx(0.075, rel=1e-12)
    assert mat.WATER.wavelength(100e3) == pytest.approx(0.015, rel=1e-12)
    assert mat.WATER.wavelength(1e6) == pytest.approx(0.0015, rel=1e-12)


def test_db_conversion_is_the_standard_8_686():
    m = mat.Material(name="t", c=1500.0, rho=1000.0, alpha0=1.0, f_ref=1.0e6)
    assert m.attenuation_db_per_m(1.0e6) == pytest.approx(8.685889638, rel=1e-9)


def test_reflection_at_a_single_interface_is_the_impedance_formula():
    """R = (Z2 - Z1)/(Z2 + Z1): el patrón contra el que se medirá el FEM."""
    z1, z2 = mat.WATER.impedance, 4.0e6
    r = (z2 - z1) / (z2 + z1)
    assert r == pytest.approx(0.45454545454545453, rel=1e-12)


def test_invalid_materials_are_rejected():
    with pytest.raises(ValueError):
        mat.Material(name="bad", c=0.0, rho=1000.0)
    with pytest.raises(ValueError):
        mat.Material(name="bad", c=1500.0, rho=-1.0)
    with pytest.raises(ValueError):
        mat.Material(name="bad", c=1500.0, rho=1000.0, alpha0=-1.0)
    with pytest.raises(ValueError):
        mat.ElasticSolid(name="bad", rho=2000.0, c_p=2000.0, c_s=3000.0)


def test_elastic_solid_lamy_parameters():
    bone = mat.ElasticSolid(name="bone", rho=1950.0, c_p=3500.0, c_s=1900.0)
    assert bone.mu == pytest.approx(1950.0 * 1900.0**2)
    assert bone.lam == pytest.approx(1950.0 * (3500.0**2 - 2 * 1900.0**2))
    # K = lambda + 2mu/3
    assert bone.bulk_modulus == pytest.approx(bone.lam + 2 * bone.mu / 3)
