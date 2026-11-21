"""STL -> malla de tetraedros con una región por pieza (U9).

Qué hace
--------
Toma uno o más STL **cerrados** y devuelve una malla de volumen con cada tetraedro etiquetado con la
región a la que pertenece. Es la pieza que hoy no existe: `meshing3d_sphere` sabe de una esfera y
`grid_tri_mesh` de cajas, pero no de una geometría propia.

Cómo decide a qué región va cada tetraedro
------------------------------------------
Por **ángulo sólido** (winding number): se suma el ángulo sólido subtendido por cada triángulo del STL
visto desde el centroide del tetraedro. Vale 4π si el punto está adentro y 0 si está afuera, sin depender
de la orientación de las caras ni de que el sólido sea convexo. Se descarta antes por bounding box, así
que el costo es proporcional a los tetraedros que están cerca de cada pieza y no al producto total.

Por qué no se usan los grupos físicos de gmsh
---------------------------------------------
Con `occ.fragment` (necesario para que dos piezas que se tocan compartan la interfaz) los tags de los
volúmenes cambian y hay que volver a mapearlos. El ángulo sólido es independiente de eso: etiqueta
geometría, no tags.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from . import stl

__all__ = ["MallaVolumen", "mallar_stl", "etiquetar_regiones"]


@dataclass
class MallaVolumen:
    nodes: np.ndarray                       # (N, 3) en metros
    tets: np.ndarray                        # (M, 4)
    regiones: np.ndarray                    # (M,) int -> índice en `nombres`
    nombres: list[str] = field(default_factory=list)
    origen: str = ""

    @property
    def n_nodes(self) -> int:
        return int(self.nodes.shape[0])

    @property
    def es_2d(self) -> bool:
        return self.tets.shape[1] == 3

    def resumen(self) -> dict:
        conteo = {n: int(np.count_nonzero(self.regiones == i)) for i, n in enumerate(self.nombres)}
        base = {
            "n_nodos": self.n_nodes,
            "n_elementos": int(self.tets.shape[0]),
            "regiones": conteo,
            "dim": 2 if self.es_2d else 3,
            "origen": self.origen,
        }
        if self.es_2d:
            base["area_mm2"] = round(_area(self.nodes, self.tets) * 1e6, 3)
        else:
            # `n_tets` se mantiene para no romper a quien ya lo leía
            base["n_tets"] = int(self.tets.shape[0])
            base["volumen_mm3"] = round(_volumen(self.nodes, self.tets) * 1e9, 2)
        return base


def _area(nodes: np.ndarray, tris: np.ndarray) -> float:
    p = nodes[tris]
    # triángulos en 2D: el producto cruz en el plano es el determinante 2x2
    a = p[:, 1] - p[:, 0]
    b = p[:, 2] - p[:, 0]
    det = a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]
    return float(np.abs(det).sum() / 2.0)


def _volumen(nodes: np.ndarray, tets: np.ndarray) -> float:
    p = nodes[tets]
    return float(np.abs(np.einsum("ij,ij->i", p[:, 1] - p[:, 0],
                                  np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]))).sum() / 6.0)


def _angulo_solido(punto: np.ndarray, tri: np.ndarray) -> float:
    """Ángulo sólido de un triángulo visto desde `punto` (formula de Van Oosterom-Strackee)."""
    a = tri[:, 0] - punto
    b = tri[:, 1] - punto
    c = tri[:, 2] - punto
    la = np.linalg.norm(a, axis=1)
    lb = np.linalg.norm(b, axis=1)
    lc = np.linalg.norm(c, axis=1)
    num = np.einsum("ij,ij->i", a, np.cross(b, c))
    den = (la * lb * lc + np.einsum("ij,ij->i", a, b) * lc
           + np.einsum("ij,ij->i", b, c) * la + np.einsum("ij,ij->i", c, a) * lb)
    return float(2.0 * np.arctan2(num, den).sum())


def etiquetar_regiones(nodes, tets, mallas: dict[str, stl.MallaSTL]) -> tuple[np.ndarray, list[str]]:
    """Devuelve ``(region_por_tet, nombres)`` por ángulo sólido, con ``-1`` para lo que no cae adentro."""
    nombres = list(mallas)
    regiones = np.full(tets.shape[0], -1, dtype=int)
    centroides = nodes[tets].mean(axis=1)

    cajas = {}
    for nombre, m in mallas.items():
        lo, hi = m.bbox()
        cajas[nombre] = (lo - 1e-12, hi + 1e-12)

    for i, c in enumerate(centroides):
        for j, nombre in enumerate(nombres):
            lo, hi = cajas[nombre]
            if np.any(c < lo) or np.any(c > hi):
                continue
            if abs(_angulo_solido(c, mallas[nombre].triangulos)) > 2.0 * np.pi:
                regiones[i] = j
                break
    return regiones, nombres


def mallar_capas(capas: list[tuple[str, float]], *, ancho: float, alto: float,
                 desde: float = 0.0, tamano: float, fino: float | None = None) -> MallaVolumen:
    """Malla de una **cadena de placas** contiguas, con la interfaz compartida por los dos volúmenes.

    Es el caso piel -> hueso -> cóclea: cada capa es una placa de espesor propio, apiladas en x. Se
    construyen como geometría analítica y se fragmentan con occ, que es lo que garantiza que la
    interfaz sea **una sola superficie** y la malla salga conforme por construcción. Con STL la
    interfaz no se puede garantizar (ver `mallar_stl`).

    Parameters
    ----------
    capas : list of ``(nombre, espesor_m)``, en el orden en que la onda las atraviesa.
    ancho, alto : tamaño transversal de las placas [m].
    desde : posición del comienzo de la primera capa en x [m].
    """
    import gmsh

    if not capas:
        raise ValueError("hacen falta capas")
    if any(e <= 0 for _, e in capas):
        raise ValueError("todas las capas necesitan espesor positivo")

    t0 = time.perf_counter()
    gmsh.initialize(interruptible=False)
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("capas")
        occ = gmsh.model.occ

        x = float(desde)
        cajas = []
        for _nombre, espesor in capas:
            cajas.append(occ.addBox(x, -alto / 2.0, -ancho / 2.0, float(espesor), alto, ancho))
            x += float(espesor)
        occ.synchronize()

        if len(cajas) > 1:
            occ.fragment([(3, t) for t in cajas], [], removeObject=True, removeTool=True)
            occ.synchronize()

        gmsh.option.setNumber("Mesh.MeshSizeMax", float(tamano))
        gmsh.option.setNumber("Mesh.MeshSizeMin", float(tamano))
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)
        gmsh.option.setNumber("Mesh.Optimize", 1)
        if fino is not None:
            superficies = [tag for dim, tag in gmsh.model.getEntities(2)]
            campo = gmsh.model.mesh.field.add("Distance")
            gmsh.model.mesh.field.setNumbers(campo, "SurfacesList", superficies)
            gmsh.model.mesh.field.setNumber(campo, "Sampling", 200)
            umbral = gmsh.model.mesh.field.add("Threshold")
            gmsh.model.mesh.field.setNumber(umbral, "InField", campo)
            gmsh.model.mesh.field.setNumber(umbral, "SizeMin", float(fino))
            gmsh.model.mesh.field.setNumber(umbral, "SizeMax", float(tamano))
            gmsh.model.mesh.field.setNumber(umbral, "DistMin", float(fino))
            gmsh.model.mesh.field.setNumber(umbral, "DistMax", float(tamano))
            gmsh.model.mesh.field.setAsBackgroundMesh(umbral)
            gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
            gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
            gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

        gmsh.model.mesh.generate(3)

        tag_nodos, coords, _ = gmsh.model.mesh.getNodes()
        orden = np.argsort(tag_nodos)
        nodes = np.asarray(coords, dtype=float).reshape(-1, 3)[orden]
        mapa = {int(t): i for i, t in enumerate(np.asarray(tag_nodos)[orden])}

        tets_lista = []
        regiones_lista = []
        # se pide volumen por volumen: `getElements(3)` sin tag no dice a qué volumen pertenece cada
        # bloque, y la región es justamente eso
        for _dim, tag_vol in gmsh.model.getEntities(3):
            tipos, _tags, conectividad = gmsh.model.mesh.getElements(3, tag=tag_vol)
            for tipo, conn in zip(tipos, conectividad):
                if tipo != 4:
                    continue
                bloques = np.asarray(conn, dtype=int).reshape(-1, 4)
                tets_lista.append(np.vectorize(mapa.get)(bloques))
                regiones_lista.append(np.full(bloques.shape[0], int(tag_vol)))
        if not tets_lista:
            raise ValueError("gmsh no generó tetraedros: revisá el tamaño de malla")
        tets = np.concatenate(tets_lista, axis=0)
        crudas = np.concatenate(regiones_lista, axis=0)
    finally:
        gmsh.finalize()

    nombres = [n for n, _ in capas]
    orden_tags = sorted(set(crudas.tolist()))
    if len(orden_tags) != len(nombres):
        raise ValueError(
            f"gmsh devolvió {len(orden_tags)} volúmenes y hay {len(nombres)} capas: la geometría no "
            f"salió como se pidió"
        )
    # el fragment conserva los tags de los volúmenes de entrada, así que se ordenan y se mapean
    mapa_tags = {tag: i for i, tag in enumerate(orden_tags)}
    regiones = np.array([mapa_tags[t] for t in crudas], dtype=int)
    # y el orden de los tags no tiene por qué coincidir con el de las capas: se verifica por el
    # centroide en x, que es lo que dice cuál es cuál
    centros = np.array([nodes[tets[regiones == i]].mean(axis=1)[:, 0].mean() for i in range(len(nombres))])
    orden_geo = np.argsort(centros)
    remapeo = np.empty_like(orden_geo)
    remapeo[orden_geo] = np.arange(len(nombres))
    regiones = remapeo[regiones]

    tele = time.perf_counter() - t0
    return MallaVolumen(nodes=nodes, tets=tets, regiones=regiones, nombres=nombres,
                        origen=f"gmsh capas · {tele:.1f} s")


def mallar_seccion(contornos: list[np.ndarray], *, ancho: float, alto: float,
                   tamano: float, desde: float = 0.0) -> MallaVolumen:
    """Malla 2D de una **sección de STL** dentro de un dominio rectangular de fluido.

    En 2D sí se puede tener la interfaz compartida: el contorno de la sección es **la misma curva** que
    el borde del agujero del fluido (`occ.cut`), así que la malla sale conforme por construcción. Es la
    diferencia con `mallar_stl`, donde dos superficies en 3D no comparten cara y el mallador se cae.

    El dominio se centra en el origen del contorno; ``desde`` corre el rectángulo en x.
    """
    import gmsh

    if not contornos:
        raise ValueError("hacen falta contornos (ver stl.seccion_por_plano)")
    if tamano <= 0:
        raise ValueError("tamano tiene que ser > 0")

    t0 = time.perf_counter()
    gmsh.initialize(interruptible=False)
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("seccion")
        occ = gmsh.model.occ

        superficies_solido = []
        for contorno in contornos:
            c = np.asarray(contorno, dtype=float)
            if len(c) < 3:
                raise ValueError("un contorno de la sección tiene menos de 3 puntos")
            # `desde` es dónde va el CENTRO de la sección dentro del dominio
            c = c - c.mean(axis=0) + np.array([float(desde), 0.0])
            puntos = [occ.addPoint(float(x), float(y), 0.0) for x, y in c]
            lineas = [occ.addLine(puntos[i], puntos[(i + 1) % len(puntos)]) for i in range(len(puntos))]
            loop = occ.addCurveLoop(lineas)
            superficies_solido.append(occ.addPlaneSurface([loop]))
        occ.synchronize()

        # addRectangle(x, y, z, dx, dy): el z va en 0 (malla 2D) y la caja se centra en `desde`, que es
        # donde está el centro de la sección; si la caja arrancara en 0 la sección podría quedar a
        # caballo del borde y el recorte sale degenerado
        caja = occ.addRectangle(float(desde) - ancho / 2, -alto / 2, 0.0, ancho, alto)
        occ.synchronize()
        # el fluido es la caja menos las secciones: la frontera del agujero es la MISMA curva que el
        # borde del sólido, así que la malla queda conforme
        fluido, _ = occ.cut([(2, caja)], [(2, s) for s in superficies_solido],
                            removeObject=True, removeTool=False)
        occ.synchronize()
        if not fluido:
            raise ValueError(
                "la sección quedó fuera del dominio: el rectángulo tiene que contenerla entera"
            )

        gmsh.option.setNumber("Mesh.MeshSizeMax", float(tamano))
        gmsh.option.setNumber("Mesh.MeshSizeMin", float(tamano))
        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.model.mesh.generate(2)

        tag_nodos, coords, _ = gmsh.model.mesh.getNodes()
        orden = np.argsort(tag_nodos)
        nodos3 = np.asarray(coords, dtype=float).reshape(-1, 3)[orden]
        mapa = {int(t): i for i, t in enumerate(np.asarray(tag_nodos)[orden])}

        triangulos = []
        regiones = []
        # una superficie por región: la del sólido y la del fluido
        for etiqueta, tags in (
            (0, [(2, s) for s in superficies_solido]),
            (1, [dim_tag for dim_tag in fluido]),
        ):
            for dim, tag in tags:
                tipos, _t, conn = gmsh.model.mesh.getElements(dim, tag=tag)
                for tipo, c in zip(tipos, conn):
                    if tipo != 2:           # 2 = triángulo
                        continue
                    bloques = np.asarray(c, dtype=int).reshape(-1, 3)
                    triangulos.append(np.vectorize(mapa.get)(bloques))
                    regiones.append(np.full(bloques.shape[0], etiqueta))
        if not triangulos:
            raise ValueError("gmsh no generó triángulos")
    finally:
        gmsh.finalize()

    nodos = nodos3[:, :2]
    tris = np.concatenate(triangulos, axis=0)
    reg = np.concatenate(regiones, axis=0)
    # gmsh no garantiza el sentido de giro: el sólido puede salir horario y `elastic2d` lo rechaza
    from .meshing import _ensure_ccw

    tris = _ensure_ccw(nodos, tris)
    tele = time.perf_counter() - t0
    return MallaVolumen(nodes=nodos, tets=tris, regiones=reg, nombres=["seccion", "fluido"],
                        origen=f"gmsh sección · {tele:.1f} s")


def mallar_stl(rutas: dict[str, str | Path], tamano: float, *,
               fino: float | None = None, algoritmo: int = 1) -> MallaVolumen:
    """Malla de volumen desde STL cerrados, una región por archivo.

    Parameters
    ----------
    rutas : ``{nombre_region: ruta_stl}``. Cada archivo tiene que ser un sólido cerrado.
    tamano : tamaño de malla objetivo [m] en el interior.
    fino : tamaño sobre las superficies (si es ``None``, ``tamano / 2``).
    """
    import gmsh

    mallas = {}
    for nombre, ruta in rutas.items():
        m = stl.leer_stl(ruta)
        stl.verificar_stl(m, nombre)
        mallas[nombre] = m

    t0 = time.perf_counter()
    gmsh.initialize(interruptible=False)
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("stl")
        occ = gmsh.model.occ

        # gmsh no importa STL como geometría BREP: hay que mergear la malla discreta, clasificar sus
        # superficies y recién ahí crear la geometría. Con eso occt puede hacer el volumen.
        def _geometria_desde_stl(ruta, nombre):
            gmsh.merge(str(ruta))
            entidades = gmsh.model.getEntities(2)
            if not entidades:
                raise ValueError(f"{nombre}: el STL no aportó ninguna superficie")
            gmsh.model.mesh.classifySurfaces(
                40.0 * np.pi / 180.0,          # ángulo entre caras que separa parches
                boundary=True,                 # crear las curvas de frontera
                forReparametrization=True,
                curveAngle=np.pi,
            )
            gmsh.model.mesh.createGeometry()
            gmsh.model.geo.synchronize()
            superficies = [tag for dim, tag in gmsh.model.getEntities(2)]
            if not superficies:
                raise ValueError(f"{nombre}: no se pudo construir geometría desde el STL")
            return superficies

        volumenes = {}
        for nombre, ruta in rutas.items():
            superficies = _geometria_desde_stl(ruta, nombre)
            loop = gmsh.model.geo.addSurfaceLoop(superficies)
            volumenes[nombre] = gmsh.model.geo.addVolume([loop])
        gmsh.model.geo.synchronize()

        gmsh.option.setNumber("Mesh.MeshSizeMax", float(tamano))
        gmsh.option.setNumber("Mesh.MeshSizeMin", float(tamano))
        gmsh.option.setNumber("Mesh.Algorithm3D", int(algoritmo))
        gmsh.option.setNumber("Mesh.Optimize", 1)
        if fino is not None:
            gmsh.option.setNumber("Mesh.MeshSizeMin", float(fino))
            superficies = [tag for dim, tag in gmsh.model.getEntities(2)]
            campo = gmsh.model.mesh.field.add("Distance")
            gmsh.model.mesh.field.setNumbers(campo, "SurfacesList", superficies)
            gmsh.model.mesh.field.setNumber(campo, "Sampling", 200)
            umbral = gmsh.model.mesh.field.add("Threshold")
            gmsh.model.mesh.field.setNumber(umbral, "InField", campo)
            gmsh.model.mesh.field.setNumber(umbral, "SizeMin", float(fino))
            gmsh.model.mesh.field.setNumber(umbral, "SizeMax", float(tamano))
            gmsh.model.mesh.field.setNumber(umbral, "DistMin", float(fino))
            gmsh.model.mesh.field.setNumber(umbral, "DistMax", float(tamano))
            gmsh.model.mesh.field.setAsBackgroundMesh(umbral)
            gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
            gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
            gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

        if len(rutas) > 1:
            raise NotImplementedError(
                "mallar_stl acepta UNA pieza por malla. Dos STL que se tocan tienen caras coincidentes "
                "y el mallador 3D se cae («A segment and a facet intersect»); fusionarlas después de "
                "mallar rompe la correspondencia con la geometría («Could not recover boundary mesh»). "
                "Para una cadena de capas usá `mallar_capas`, que construye las placas como geometría "
                "analítica y las fragmenta, y así la interfaz queda compartida por los dos volúmenes."
            )
        gmsh.model.mesh.generate(3)

        tag_nodos, coords, _ = gmsh.model.mesh.getNodes()
        orden = np.argsort(tag_nodos)
        nodes = np.asarray(coords, dtype=float).reshape(-1, 3)[orden]
        mapa = {int(t): i for i, t in enumerate(np.asarray(tag_nodos)[orden])}

        tipos, tags, conectividad = gmsh.model.mesh.getElements(3)
        tets_lista = []
        for tipo, conn in zip(tipos, conectividad):
            if tipo != 4:                       # 4 = tetraedro
                continue
            bloques = np.asarray(conn, dtype=int).reshape(-1, 4)
            tets_lista.append(np.vectorize(mapa.get)(bloques))
        if not tets_lista:
            raise ValueError("gmsh no generó tetraedros: revisá el tamaño de malla contra la geometría")
        tets = np.concatenate(tets_lista, axis=0)
    finally:
        gmsh.finalize()

    tele = time.perf_counter() - t0
    regiones, nombres = etiquetar_regiones(nodes, tets, mallas)
    huerfanos = int(np.count_nonzero(regiones < 0))
    if huerfanos:
        # no se descartan en silencio: se quedan con la región 0 y el resumen lo dice
        regiones[regiones < 0] = 0
    malla = MallaVolumen(nodes=nodes, tets=tets, regiones=regiones, nombres=nombres,
                         origen=f"gmsh · {tele:.1f} s")
    malla.resumen()["huerfanos"] = huerfanos
    return malla
