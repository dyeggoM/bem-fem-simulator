"""El lector STL tiene que distinguir un sólido de una superficie rota (U9).

Del lector depende que el mallado no invente: un STL con un agujero malla igual y el solve devuelve un
número, que es peor que un error. Estos tests fijan las tres cosas que el lector promete: leer binario y
ASCII, detectar que la malla no es cerrada, y detectar normales invertidas.
"""

from __future__ import annotations

import struct

import numpy as np
import pytest

from ultrasound_sim import stl


def _escribir_binario(ruta, malla: stl.MallaSTL):
    """STL binario desde una malla (lo que exportaría cualquier CAD)."""
    with open(ruta, "wb") as f:
        f.write(b"\0" * 80)
        f.write(struct.pack("<I", malla.n_triangulos))
        for tri, n in zip(malla.triangulos, malla.normales):
            f.write(struct.pack("<3f", *n))
            for v in tri:
                f.write(struct.pack("<3f", *v))
            f.write(struct.pack("<H", 0))


def _escribir_ascii(ruta, malla: stl.MallaSTL):
    lineas = ["solid prueba"]
    for tri, n in zip(malla.triangulos, malla.normales):
        lineas.append(f"  facet normal {n[0]:.6e} {n[1]:.6e} {n[2]:.6e}")
        lineas.append("    outer loop")
        for v in tri:
            lineas.append(f"      vertex {v[0]:.6e} {v[1]:.6e} {v[2]:.6e}")
        lineas.append("    endloop")
        lineas.append("  endfacet")
    lineas.append("endsolid prueba")
    ruta.write_text("\n".join(lineas), encoding="utf-8")


def test_el_cubo_generado_es_cerrado_y_con_volumen_positivo(tmp_path):
    malla = stl.cubo_stl(0.01)          # 1 cm
    detalle = stl.verificar_stl(malla, "cubo")
    print("\n  cubo:", detalle)
    assert abs(detalle["volumen_m3"] - 1e-6) < 1e-12      # (0.01 m)^3
    assert detalle["aristas_una_cara"] == 0


@pytest.mark.parametrize("formato", ["binario", "ascii"])
def test_lee_los_dos_formatos_y_da_lo_mismo(tmp_path, formato):
    malla = stl.cubo_stl(0.01)
    ruta = tmp_path / ("cubo_bin.stl" if formato == "binario" else "cubo_asc.stl")
    (_escribir_binario if formato == "binario" else _escribir_ascii)(ruta, malla)

    leida = stl.leer_stl(ruta)
    print(f"\n  {formato}: {leida.n_triangulos} triángulos · origen {leida.origen}")
    assert leida.n_triangulos == malla.n_triangulos
    np.testing.assert_allclose(leida.triangulos, malla.triangulos, atol=1e-9)
    assert abs(stl.volumen_con_signo(leida) - 1e-6) < 1e-12


def test_un_stl_abierto_se_rechaza_con_el_motivo(tmp_path):
    """Al cubo se le saca una cara: el lector tiene que decir «agujero», no mallar."""
    malla = stl.cubo_stl(0.01)
    abierta = stl.MallaSTL(triangulos=malla.triangulos[:-1], normales=malla.normales[:-1])
    cerrada, detalle = stl.escerrada(abierta.triangulos)
    print("\n  sin una cara:", detalle)
    assert not cerrada
    assert detalle["aristas_una_cara"] > 0
    with pytest.raises(ValueError, match="agujero"):
        stl.verificar_stl(abierta, "abierto")


def test_un_stl_con_las_normales_invertidas_se_rechaza(tmp_path):
    """Mismo sólido con las caras dadas vuelta: cerrado, pero de volumen negativo."""
    malla = stl.cubo_stl(0.01)
    invertida = stl.MallaSTL(triangulos=malla.triangulos[:, ::-1], normales=-malla.normales)
    cerrada, _ = stl.escerrada(invertida.triangulos)
    print(f"\n  invertida: cerrada={cerrada} · volumen {stl.volumen_con_signo(invertida):.3e} m³")
    assert cerrada, "invertir las caras no abre la superficie"
    assert stl.volumen_con_signo(invertida) < 0
    with pytest.raises(ValueError, match="hacia adentro"):
        stl.verificar_stl(invertida, "invertida")


def test_un_archivo_truncado_no_se_confunde_con_uno_valido(tmp_path):
    malla = stl.cubo_stl(0.01)
    ruta = tmp_path / "cubo.stl"
    _escribir_binario(ruta, malla)
    datos = ruta.read_bytes()
    roto = tmp_path / "roto.stl"
    roto.write_bytes(datos[:-7])          # se cortan 7 bytes: el tamaño ya no cierra
    with pytest.raises(ValueError):
        stl.leer_stl(roto)
