"""Importador de transductores: catálogo TOML → geometría → aperturas.

Lo que se prueba acá no es la física (esa está en test_aperture/test_array_focus) sino que el
importador no invente ni esconda nada: campo faltante falla con nombre, geometría no soportada
falla en vez de aproximarse, y lo no verificado queda marcado.
"""

import numpy as np
import pytest

from ultrasound_sim import aperture, materials, transducers
from ultrasound_sim.transducers import CatalogError, TransducerSpec

WATER = materials.WATER

CATALOGO_VALIDO = """
[transducers.mono]
name = "Mono de prueba"
frequency_hz = 100000.0
n_elements = 1
pitch_m = 0.03
source = "test"

[transducers.arr]
name = "Arreglo de prueba"
frequency_hz = 100000.0
n_elements = 8
pitch_m = 0.0075
focus_distance_m = 0.04
source = "datasheet de mentira para el test"
verified = true
notes = "sólo test"
"""


def test_parse_catalog_reads_every_entry():
    cat = transducers.parse_catalog(CATALOGO_VALIDO)
    assert set(cat) == {"mono", "arr"}
    arr = cat["arr"]
    assert arr.n_elements == 8
    assert arr.pitch_m == pytest.approx(0.0075)
    assert arr.focus_distance_m == pytest.approx(0.04)
    assert arr.verified is True
    assert cat["mono"].verified is False, "verified debe ser False por defecto"


def test_missing_required_field_names_the_field_and_the_entry():
    texto = """
[transducers.roto]
name = "sin frecuencia"
n_elements = 4
pitch_m = 0.01
source = "x"
"""
    with pytest.raises(CatalogError) as exc:
        transducers.parse_catalog(texto)
    msg = str(exc.value)
    assert "roto" in msg and "frequency_hz" in msg


def test_empty_source_is_rejected():
    """Procedencia obligatoria (D5): sin fuente, el número no entra."""
    texto = """
[transducers.sin_fuente]
name = "x"
frequency_hz = 100000.0
n_elements = 1
pitch_m = 0.01
source = "   "
"""
    with pytest.raises(CatalogError, match="source"):
        transducers.parse_catalog(texto)


def test_kerf_is_refused_instead_of_approximated():
    texto = """
[transducers.con_kerf]
name = "x"
frequency_hz = 100000.0
n_elements = 4
pitch_m = 0.01
element_width_m = 0.008
source = "datasheet"
"""
    with pytest.raises(CatalogError, match="kerf"):
        transducers.parse_catalog(texto)


@pytest.mark.parametrize(
    "override, needle",
    [
        ("frequency_hz = -1.0", "frequency_hz"),
        ("frequency_hz = 0.0", "frequency_hz"),
        ("pitch_m = 0.0", "pitch_m"),
        ("n_elements = 0", "n_elements"),
        ("focus_distance_m = -0.02", "focus_distance_m"),
        ("amplitude = 0.0", "amplitude"),
    ],
)
def test_invalid_values_are_rejected(override, needle):
    base = {
        "name": '"x"',
        "frequency_hz": "100000.0",
        "n_elements": "4",
        "pitch_m": "0.01",
        "source": '"d"',
    }
    key = override.split("=")[0].strip()
    base[key] = override.split("=", 1)[1].strip()
    texto = "[transducers.t]\n" + "\n".join(f"{k} = {v}" for k, v in base.items()) + "\n"
    with pytest.raises(CatalogError, match=needle):
        transducers.parse_catalog(texto)


def test_unknown_field_is_reported():
    texto = """
[transducers.t]
name = "x"
frequency_hz = 100000.0
n_elements = 4
pitch_m = 0.01
source = "d"
frequncy_hz = 5
"""
    with pytest.raises(CatalogError, match="desconocidos"):
        transducers.parse_catalog(texto)


def test_shipped_catalog_loads_and_marks_what_is_not_verified():
    cat = transducers.load_catalog()
    assert cat, "el catálogo del repo no debería estar vacío"
    for key, spec in cat.items():
        assert spec.source.strip(), f"{key} sin procedencia"
        # Regla: nada del repo se declara verificado sin datasheet. Si algún día se carga uno real,
        # este test se actualiza a propósito, no por accidente.
        assert spec.verified is False, f"{key} está marcado verificado sin datasheet en el repo"


def test_derived_geometry_matches_the_wavelength():
    cat = transducers.parse_catalog(CATALOGO_VALIDO)
    arr = cat["arr"]
    lam = WATER.wavelength(1e5)
    assert arr.wavelength(WATER) == pytest.approx(lam)
    assert arr.aperture_width_m == pytest.approx(8 * 0.0075)
    assert arr.aperture_over_wavelength(WATER) == pytest.approx(0.06 / lam)
    assert arr.pitch_over_wavelength(WATER) == pytest.approx(0.0075 / lam)
    assert arr.f_number == pytest.approx(0.04 / 0.06)


def test_grating_lobe_only_when_pitch_is_coarse():
    cat = transducers.parse_catalog(CATALOGO_VALIDO)
    lam = WATER.wavelength(1e5)  # 15 mm
    fino = TransducerSpec(key="f", name="f", frequency_hz=1e5, n_elements=8,
                          pitch_m=lam / 4, source="test")
    grueso = TransducerSpec(key="g", name="g", frequency_hz=1e5, n_elements=8,
                            pitch_m=lam, source="test")
    assert transducers.grating_lobe_angle_deg(fino, WATER) is None
    assert transducers.grating_lobe_angle_deg(grueso, WATER) == pytest.approx(90.0)
    # a λ/2 el seno da 2.0 > 1: no hay lóbulo de grataje real
    medio = TransducerSpec(key="m", name="m", frequency_hz=1e5, n_elements=8,
                          pitch_m=lam / 2, source="test")
    assert transducers.grating_lobe_angle_deg(medio, WATER) is None
    assert cat["arr"].n_elements == 8  # sanity


def test_apertures_are_focused_at_the_declared_distance():
    cat = transducers.parse_catalog(CATALOGO_VALIDO)
    aps = cat["arr"].apertures(WATER)
    assert len(aps) == 8
    # La fase de cada elemento es exactamente exp(-i k |F - y_j|) con F = (focus_distance, 0).
    # Comparar fases envueltas con np.angle sería ambiguo: se compara el fasor.
    k = WATER.wavenumber(1e5).real
    esperado = np.array([np.exp(-1j * k * np.hypot(0.04, a.center_y)) for a in aps])
    np.testing.assert_allclose([a.amplitude for a in aps], esperado, atol=1e-12)
    assert np.allclose(np.abs([a.amplitude for a in aps]), 1.0)


def test_catalog_transducer_reproduces_the_4b_reference_field():
    """El baseline del repo (λ/8 a 100 kHz) tiene que dar el FWHM del estudio de Fase 4B."""
    cat = transducers.load_catalog()
    spec = cat["referencia_lambda8_100k"]
    aps = spec.apertures(WATER, focus_point=(0.04, 0.0))
    ys = np.linspace(-20e-3, 20e-3, 401)
    pts = np.column_stack([np.full(ys.size, 0.0295), ys])  # plano del máximo axial real (D13)
    p = aperture.pressure(pts, spec.frequency_hz, WATER, aps)
    intens = np.abs(p) ** 2
    ipk = int(intens.argmax())
    half = intens[ipk] / 2

    def cruce(side):
        if side < 0:
            sy, si = ys[: ipk + 1][::-1], intens[: ipk + 1][::-1]
        else:
            sy, si = ys[ipk:], intens[ipk:]
        idx = int(np.where(si <= half)[0][0])
        return sy[idx - 1] + (half - si[idx - 1]) * (sy[idx] - sy[idx - 1]) / (si[idx] - si[idx - 1])

    fwhm = cruce(1) - cruce(-1)
    print(f"\n  referencia λ/8: pico en y={ys[ipk]*1e3:+.3f} mm, FWHM={fwhm*1e3:.3f} mm "
          f"(estudio 4B: 8.766 mm)")
    assert ys[ipk] == pytest.approx(0.0, abs=0.5e-3)
    assert fwhm * 1e3 == pytest.approx(8.766, abs=0.15)
