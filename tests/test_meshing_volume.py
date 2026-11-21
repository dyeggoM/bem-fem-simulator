"""STL -> malla de volumen, y la cadena de capas con interfaz compartida (U9).

Lo que estos tests fijan es lo que hace que el número que sale del solve signifique algo: que el volumen
mallado sea el que se pidió, que no haya nodos duplicados (una malla con duplicados tiene la interfaz
rota y el acople entre regiones no se arma), y que cada tetraedro sepa a qué capa pertenece.
"""

from __future__ import annotations

import struct

import numpy as np
import pytest

from ultrasound_sim import meshing_volume, stl


def _escribir_stl_binario(ruta, malla: stl.MallaSTL):
    with open(ruta, "wb") as f:
        f.write(b"\0" * 80)
        f.write(struct.pack("<I", malla.n_triangulos))
        for tri, n in zip(malla.triangulos, malla.normales):
            f.write(struct.pack("<3f", *n))
            for v in tri:
                f.write(struct.pack("<3f", *v))
            f.write(struct.pack("<H", 0))


def test_una_pieza_da_el_volumen_exacto(tmp_path):
    ruta = tmp_path / "placa.stl"
    _escribir_stl_binario(ruta, stl.cubo_stl(0.010, centro=(0.005, 0, 0)))
    malla = meshing_volume.mallar_stl({"placa": ruta}, 0.003, fino=0.002)
    r = malla.resumen()
    print(f"\n  {r}")
    assert r["n_tets"] > 100
    assert abs(r["volumen_mm3"] - 1000.0) < 1.0            # cubo de 1 cm
    assert r["regiones"] == {"placa": r["n_tets"]}


def test_la_cadena_de_capas_es_conforme_y_volumen_exacto():
    """piel 2 mm + hueso 5 mm + liquido 4 mm en una seccion de 20x20 mm: 4400 mm3."""
    capas = [("piel", 2e-3), ("hueso", 5e-3), ("liquido", 4e-3)]
    malla = meshing_volume.mallar_capas(capas, ancho=0.020, alto=0.020, tamano=0.0012, fino=0.0008)
    r = malla.resumen()
    print(f"\n  {r}")

    assert abs(r["volumen_mm3"] - 4400.0) < 20.0
    # sin duplicados: es lo que hace que la interfaz sea una sola superficie y el acople se pueda armar
    unicos = np.unique(np.round(malla.nodes, 10), axis=0).shape[0]
    assert unicos == malla.n_nodes, "la malla tiene nodos duplicados: la interfaz está rota"

    # cada capa tiene que caer en su rango de x, en el orden pedido
    limites = np.cumsum([0.0] + [e for _, e in capas]) * 1e3
    for i, (nombre, espesor) in enumerate(capas):
        sel = malla.regiones == i
        assert sel.any(), f"la capa {nombre} quedó sin tetraedros"
        xs = malla.nodes[malla.tets[sel]].mean(axis=1)[:, 0] * 1e3
        print(f"    {nombre}: {int(sel.sum())} tets · x [{xs.min():.2f}, {xs.max():.2f}] mm "
              f"· pedido [{limites[i]:.1f}, {limites[i + 1]:.1f}]")
        assert xs.min() >= limites[i] - 1e-6
        assert xs.max() <= limites[i + 1] + 1e-6
        assert malla.nombres[i] == nombre


def test_dos_stl_sueltos_fallan_diciendo_por_que(tmp_path):
    """No se mallan dos piezas que se tocan: la interfaz no se puede garantizar y hay que decirlo."""
    a = tmp_path / "a.stl"
    b = tmp_path / "b.stl"
    _escribir_stl_binario(a, stl.cubo_stl(0.010, centro=(0.005, 0, 0)))
    _escribir_stl_binario(b, stl.cubo_stl(0.010, centro=(0.015, 0, 0)))
    with pytest.raises(NotImplementedError, match="mallar_capas"):
        meshing_volume.mallar_stl({"a": a, "b": b}, 0.0025)


def test_el_mallador_rechaza_un_stl_abierto_antes_de_mallar(tmp_path):
    malla = stl.cubo_stl(0.010)
    abierta = stl.MallaSTL(triangulos=malla.triangulos[:-1], normales=malla.normales[:-1])
    ruta = tmp_path / "abierto.stl"
    _escribir_stl_binario(ruta, abierta)
    with pytest.raises(ValueError, match="no es cerrada"):
        meshing_volume.mallar_stl({"roto": ruta}, 0.003)
