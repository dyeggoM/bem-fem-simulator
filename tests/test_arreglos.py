"""Cargar un arreglo de fuentes desde archivo tiene que fallar cuando el archivo no describe un arreglo.

El objetivo de estos tests es que un archivo mal formado no entre nunca al solver: 32 filas en metros,
dos elementos pisados o una amplitud negativa producen un campo que parece válido y no lo es.
"""

from __future__ import annotations

import numpy as np
import pytest

from ultrasound_sim import arreglos


CSV_BUENO = """y_mm,ancho_mm,amp,fase_deg
-3.75,3.75,1.0,0.0
0.00,3.75,1.0,0.0
3.75,3.75,1.0,-90.0
"""


def test_csv_lee_y_convierte_a_metros():
    fuentes = arreglos.cargar_arreglo(CSV_BUENO)
    print(f"\n  {len(fuentes)} fuentes · y = {[round(f['y'] * 1e3, 3) for f in fuentes]} mm")
    assert len(fuentes) == 3
    np.testing.assert_allclose([f["y"] for f in fuentes], [-3.75e-3, 0.0, 3.75e-3])
    np.testing.assert_allclose([f["width"] for f in fuentes], 3.75e-3)
    assert all(f["amp"] == 1.0 for f in fuentes)
    # la columna dice grados: −90° tiene que entrar como −π/2 rad
    np.testing.assert_allclose(fuentes[2]["phase"], -np.pi / 2)


def test_acepta_radianes_en_archivos_viejos():
    """Los archivos con `fase_rad` (o el alias pelado `phase`) siguen significando radianes."""
    viejo = "y_mm,ancho_mm,amp,fase_rad\n0,2,1,1.5708\n"
    nuevo = "y_mm,ancho_mm,amp,phase\n0,2,1,1.5708\n"
    np.testing.assert_allclose(arreglos.cargar_arreglo(viejo)[0]["phase"], 1.5708)
    np.testing.assert_allclose(arreglos.cargar_arreglo(nuevo)[0]["phase"], 1.5708)


def test_la_fase_se_envuelve_a_180_grados():
    """Una fase de foco vale −k·r y en grados puede ser −5674.6°: es la misma fase módulo 360, pero
    ilegible. Se guarda envuelta a (−180, 180]."""
    texto = "y_mm,ancho_mm,amp,fase_deg\n0,2,1,-5674.6\n"
    fase = arreglos.cargar_arreglo(texto)[0]["phase"]
    print(f"\n  −5674.6° → {fase:.4f} rad → {fase * 180 / np.pi:.2f}°")
    assert -np.pi <= fase <= np.pi
    np.testing.assert_allclose(fase * 180 / np.pi, 85.4, atol=0.05)


def test_las_dos_columnas_de_fase_a_la_vez_usan_grados():
    """Si vienen las dos, manda la que declara la unidad: `fase_deg`."""
    texto = "y_mm,ancho_mm,amp,fase_deg,fase_rad\n0,2,1,90,0\n"
    np.testing.assert_allclose(arreglos.cargar_arreglo(texto)[0]["phase"], np.pi / 2)


def test_json_acepta_las_dos_formas():
    lista = '[{"y_mm": -2, "ancho_mm": 4, "amp": 1, "fase_rad": 0},' \
            ' {"y_mm": 2, "ancho_mm": 4, "amp": 1, "fase_rad": 0.5}]'
    objeto = '{"fuentes": [{"y_mm": -2, "ancho_mm": 4, "amp": 1, "fase_rad": 0}]}'
    assert len(arreglos.cargar_arreglo(lista)) == 2
    assert len(arreglos.cargar_arreglo(objeto)) == 1


def test_acepta_los_nombres_en_ingles_que_exporta_el_sistema():
    """El exportador escribe la cabecera en inglés; volver a cargar lo exportado tiene que funcionar."""
    fuentes = arreglos.cargar_arreglo(CSV_BUENO)
    texto = arreglos.exportar_arreglo(fuentes)
    print(f"\n  exportado:\n{texto}")
    assert texto.splitlines()[0] == "y_mm,ancho_mm,amp,fase_deg"
    ida = arreglos.cargar_arreglo(texto)
    np.testing.assert_allclose([f["y"] for f in ida], [f["y"] for f in fuentes], atol=1e-12)


def test_dos_elementos_solapados_se_rechazan_diciendo_cuales():
    texto = """y_mm,ancho_mm,amp,fase_rad
0,4,1,0
2,4,1,0
"""
    with pytest.raises(arreglos.ErrorArreglo, match="se solapan"):
        arreglos.cargar_arreglo(texto)


def test_una_amplitud_negativa_se_rechaza():
    texto = """y_mm,ancho_mm,amp,fase_rad
0,2,1,0
4,2,-1,0
"""
    with pytest.raises(arreglos.ErrorArreglo, match="amplitud negativa"):
        arreglos.cargar_arreglo(texto)


def test_un_ancho_cero_se_rechaza():
    texto = """y_mm,ancho_mm,amp,fase_rad
0,0,1,0
"""
    with pytest.raises(arreglos.ErrorArreglo, match="ancho_mm"):
        arreglos.cargar_arreglo(texto)


def test_falta_una_columna_y_se_dice_cual():
    texto = """y_mm,amp,fase_rad
0,1,0
"""
    with pytest.raises(arreglos.ErrorArreglo, match="ancho_mm"):
        arreglos.cargar_arreglo(texto)


def test_un_archivo_vacio_o_sin_fuentes_falla():
    with pytest.raises(arreglos.ErrorArreglo, match="vacío"):
        arreglos.cargar_arreglo("")
    with pytest.raises(arreglos.ErrorArreglo, match="ninguna fila"):
        arreglos.cargar_arreglo("y_mm,ancho_mm,amp,fase_rad\n")
    with pytest.raises(arreglos.ErrorArreglo, match="ninguna fuente"):
        arreglos.cargar_arreglo("[]")
