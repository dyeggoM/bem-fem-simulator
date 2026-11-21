"""STL → sección → malla → acople: lo que hay que medir para poder decir que sirve.

Tres cosas, en orden de confianza:

1. **El corte** de un cilindro por su plano medio tiene que dar un círculo con el área exacta.
2. **La malla** tiene que ser conforme (sin nodos duplicados) y con las dos regiones.
3. **El acople**, que es lo que importa: un disco con las propiedades del fluido **no perturba la onda**
   (|p| = 1), y un disco de acero **hace sombra**. Sin el caso transparente, "hay sombra" también lo daría
   un acople roto.

El error que costó encontrar esto (y que está anotado en D60): la ABC iba en las cuatro caras del dominio.
En las laterales la onda plana tiene `∂p/∂n = 0` y la ABC impone `∂p/∂n = i·k·p`: refleja. Es el mismo
error de D56, repetido. Con ABC sólo en entrada y salida el caso transparente pasa de 1.15 a 1.016.
"""

from __future__ import annotations

import numpy as np
import pytest

from ultrasound_sim import meshing_volume, stl
from ultrasound_sim.fem import helmholtz2d, multicapa
from ultrasound_sim.materials import WATER, ElasticSolid
from ultrasound_sim.meshing import boundary_edges, edge_midpoints

F = 250e3
ANCHO, ALTO = 0.040, 0.030


def cilindro(radio=0.005, alto=0.020, n=64):
    ang = np.linspace(0, 2 * np.pi, n, endpoint=False)
    base = np.stack([radio * np.cos(ang), radio * np.sin(ang), np.full(n, -alto / 2)], axis=1)
    tapa = base.copy()
    tapa[:, 2] = alto / 2
    tris = []
    for i in range(n):
        j = (i + 1) % n
        tris += [[base[i], base[j], tapa[j]], [base[i], tapa[j], tapa[i]]]
    cb = np.array([0, 0, -alto / 2])
    ct = np.array([0, 0, alto / 2])
    for i in range(n):
        j = (i + 1) % n
        tris += [[cb, base[j], base[i]], [ct, tapa[i], tapa[j]]]
    tri = np.array(tris, dtype=float)
    nrm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    nrm /= np.linalg.norm(nrm, axis=1, keepdims=True)
    return stl.MallaSTL(triangulos=tri, normales=nrm, origen="generado")


def test_el_corte_de_un_cilindro_da_un_circulo():
    contornos = stl.seccion_por_plano(cilindro(), "z", 0.0)
    assert len(contornos) == 1
    c = contornos[0]
    area = abs(stl.area_contorno(c))
    exacta = np.pi * 0.005 ** 2
    print(f"\n  sección: {len(c)} puntos · área {area * 1e6:.4f} mm² · exacta {exacta * 1e6:.4f} "
          f"· error {100 * abs(area - exacta) / exacta:.3f} %")
    assert abs(area - exacta) / exacta < 0.01
    # centrada en el origen y con el radio correcto
    assert abs(c.mean(axis=0)).max() < 1e-9
    assert abs(np.linalg.norm(c, axis=1).max() - 0.005) < 1e-9


def test_un_plano_que_no_corta_falla_diciendo_que():
    with pytest.raises(ValueError, match="no corta"):
        stl.seccion_por_plano(cilindro(alto=0.010), "z", 0.5)


def _malla_y_resolver(material, tam=0.0006):
    contornos = stl.seccion_por_plano(cilindro(), "z", 0.0)
    malla = meshing_volume.mallar_seccion(contornos, ancho=ANCHO, alto=ALTO, tamano=tam, desde=0.012)
    nodos, tris, reg = malla.nodes, malla.tets, malla.regiones
    ar = boundary_edges(tris)
    mid = edge_midpoints(nodos, ar)
    # ABC SÓLO en entrada y salida (ver el docstring y D60). Los bordes se leen de la malla: el dominio
    # se centra en la sección, así que no empieza en x = 0
    x0, x1 = float(nodos[:, 0].min()), float(nodos[:, 0].max())
    abc = ar[(np.abs(mid[:, 0] - x0) < 1e-9) | (np.abs(mid[:, 0] - x1) < 1e-9)]
    medio = np.where(reg == 0, "disco", "agua").astype(object)
    r = multicapa.solve_multicapa(nodos, tris, medio, {"agua": WATER}, {"disco": material}, F,
                                  abc_edges=abc, onda=helmholtz2d.plane_wave(F, WATER, (1.0, 0.0)))
    con_p = np.asarray(r.dof_p, dtype=int) >= 0
    banda = con_p & (np.abs(nodos[:, 1]) < 0.002)
    x = nodos[banda, 0]
    p = np.abs(r.p_todos)[banda]
    # delante y detrás del disco (centro en 0.012 m, radio 5 mm)
    return malla, {"delante": float(p[x < 0.004].mean()), "detras": float(p[x > 0.021].mean()),
                   "max": float(np.abs(r.p_todos).max())}


def test_la_malla_de_la_seccion_es_conforme():
    contornos = stl.seccion_por_plano(cilindro(), "z", 0.0)
    malla = meshing_volume.mallar_seccion(contornos, ancho=ANCHO, alto=ALTO, tamano=0.0008, desde=0.012)
    unicos = np.unique(np.round(malla.nodes, 12), axis=0).shape[0]
    print(f"\n  {malla.resumen()}")
    assert unicos == malla.n_nodes, "nodos duplicados: la interfaz está rota"
    assert set(np.unique(malla.regiones).tolist()) == {0, 1}, "faltan regiones (sólido o fluido)"
    assert abs(malla.resumen()["area_mm2"] - ANCHO * ALTO * 1e6) < 1.0


def test_un_disco_transparente_no_perturba_la_onda():
    """El control que distingue un acople bueno de uno roto."""
    transparente = ElasticSolid(name="t", rho=WATER.rho, c_p=1500.0, c_s=30.0, provenance="test")
    _, m = _malla_y_resolver(transparente, tam=0.0005)
    print(f"\n  transparente: delante {m['delante']:.4f} · detrás {m['detras']:.4f} · máx {m['max']:.4f}")
    assert abs(m["detras"] - 1.0) < 0.05, "un disco del mismo material que el fluido tendría que transmitir 1"


def test_un_disco_de_acero_hace_sombra():
    acero = ElasticSolid(name="a", rho=7800.0, c_p=5900.0, c_s=3200.0, provenance="literatura")
    _, m = _malla_y_resolver(acero, tam=0.0005)
    print(f"\n  acero: delante {m['delante']:.4f} · detrás {m['detras']:.4f} · máx {m['max']:.4f}")
    assert m["detras"] < 0.8, "un cilindro de acero tiene que bloquear parte de la onda"
    assert m["delante"] > 1.05, "delante tiene que haber incidente más reflejada"
