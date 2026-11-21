"""Malla 3D de una esfera elástica sumergida en una cáscara de fluido (Fase 6b).

Construcción con Gmsh
---------------------
El mallador propio de este paquete (`meshing3d`) sólo sabe de cajas estructuradas, y para la esfera
eléctrica hace falta una malla **conforme sobre una superficie curva con dos regiones**. Se probó
armarla a mano con una icósfera y capas de prismas y el split de prismas en tetraedros **no es
conforme** (medido: 2 240 caras internas en vez de 320, cada una apareciendo una sola vez). En vez de
inventar el split, se usa Gmsh, que es lo que el plan preveía para esta fase.

Dos volúmenes:

- ``r < radius``  : la esfera (sólido elástico), elementos ``tets_solid``.
- ``radius < r < outer_radius`` : el fluido, elementos ``tets_fluid``.

La interfaz se comparte nodo a nodo (Gmsh fractura el sólido con la cáscara), así que el acoplamiento
fluido-sólido se puede ensamblar sin multiplicadores de Lagrange ni proyecciones. Las caras de
interfaz se **calculan de la malla** (caras que aparecen una vez en cada región) en vez de pedírselas
a Gmsh: así la malla se valida contra sí misma y no contra el bookkeeping del mallador.

Unidades: SI (m). ``spacing`` es el tamaño de elemento objetivo; ``refine_interface`` permite achicar
los elementos cerca de la interfaz, que es donde está la curvatura.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["SphereShellMesh", "sphere_shell_tet_mesh"]


@dataclass(frozen=True)
class SphereShellMesh:
    """Malla de esfera + cáscara de fluido con la interfaz conforme.

    Atributos
    ---------
    nodes : ``(N, 3)`` coordenadas [m].
    tets_solid : ``(Ms, 4)`` tetraedros de ``r < radius`` (índices globales de nodo).
    tets_fluid : ``(Mf, 4)`` tetraedros de ``radius < r < outer_radius``.
    interface_faces : ``(Fi, 3)`` triángulos de la superficie ``r = radius``, orientados con la normal
        saliente del sólido (o sea, hacia el fluido).
    radius, outer_radius : radios [m].
    """

    nodes: np.ndarray
    tets_solid: np.ndarray
    tets_fluid: np.ndarray
    interface_faces: np.ndarray
    radius: float
    outer_radius: float

    @property
    def n_nodes(self) -> int:
        return int(self.nodes.shape[0])

    def volume(self, tets: np.ndarray) -> float:
        p = self.nodes[tets]
        return float(np.sum(np.einsum("ij,ij->i", np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]),
                                      p[:, 3] - p[:, 0]) / 6.0))

    def solid_volume(self) -> float:
        return self.volume(self.tets_solid)

    def fluid_volume(self) -> float:
        return self.volume(self.tets_fluid)

    def interface_normals(self) -> np.ndarray:
        """Normal (unitaria, saliente del sólido) de cada cara de interfaz."""
        p = self.nodes[self.interface_faces]
        n = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
        norm = np.linalg.norm(n, axis=1)
        return n / np.where(norm > 0, norm, 1.0)[:, None]

    def interface_areas(self) -> np.ndarray:
        p = self.nodes[self.interface_faces]
        n = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
        return 0.5 * np.linalg.norm(n, axis=1)


def _tet_faces(tets: np.ndarray) -> np.ndarray:
    """Las cuatro caras de cada tetraedro, ``(4M, 3)``."""
    return np.vstack([
        tets[:, [0, 1, 2]], tets[:, [0, 1, 3]], tets[:, [0, 2, 3]], tets[:, [1, 2, 3]],
    ])


def _faces_once(faces: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Caras que aparecen una sola vez, y su conteo (para detectar no conformidad)."""
    key = np.sort(faces, axis=1)
    uniq, inverse, counts = np.unique(key, axis=0, return_inverse=True, return_counts=True)
    return uniq[counts == 1], counts


def sphere_shell_tet_mesh(
    radius: float,
    outer_radius: float,
    spacing: float,
    *,
    refine_interface: float | None = None,
) -> SphereShellMesh:
    """Malla esfera + cáscara de fluido con Gmsh.

    Parámetros
    ----------
    radius : radio de la esfera sólida [m].
    outer_radius : radio de la frontera exterior del fluido [m].
    spacing : tamaño de elemento objetivo [m].
    refine_interface : si se da, tamaño de elemento en la interfaz [m]; útil porque la curvatura
        está en ``r = radius`` y el ABC en ``r = outer_radius``.
    """
    if not (radius > 0 and outer_radius > radius):
        raise ValueError(f"hacen falta 0 < radius < outer_radius; llegó {radius!r}, {outer_radius!r}")
    if spacing <= 0:
        raise ValueError(f"spacing tiene que ser > 0; llegó {spacing!r}")

    import gmsh

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("esfera_sumergida")
        occ = gmsh.model.occ

        solido_vol = occ.addSphere(0.0, 0.0, 0.0, radius)
        fluido_vol = occ.addSphere(0.0, 0.0, 0.0, outer_radius)
        # fragment deja la superficie de la esfera compartida por los dos volúmenes
        occ.fragment([(3, fluido_vol)], [(3, solido_vol)])
        occ.synchronize()

        gmsh.option.setNumber("Mesh.MeshSizeMax", spacing)
        gmsh.option.setNumber("Mesh.MeshSizeMin", spacing)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)      # Delaunay, robusto para esferas
        gmsh.option.setNumber("Mesh.Optimize", 1)

        if refine_interface is not None:
            # tamaño chico sobre la superficie de la esfera, crecimiento hacia afuera
            superficies = [tag for dim, tag in gmsh.model.getEntities(2)
                           if _surface_is_sphere(gmsh, dim, tag, radius)]
            campo = gmsh.model.mesh.field.add("Distance")
            gmsh.model.mesh.field.setNumbers(campo, "SurfacesList", superficies)
            gmsh.model.mesh.field.setNumber(campo, "Sampling", 200)
            umbral = gmsh.model.mesh.field.add("Threshold")
            gmsh.model.mesh.field.setNumber(umbral, "InField", campo)
            gmsh.model.mesh.field.setNumber(umbral, "SizeMin", refine_interface)
            gmsh.model.mesh.field.setNumber(umbral, "SizeMax", spacing)
            gmsh.model.mesh.field.setNumber(umbral, "DistMin", radius * 0.25)
            gmsh.model.mesh.field.setNumber(umbral, "DistMax", outer_radius - radius)
            gmsh.model.mesh.field.setAsBackgroundMesh(umbral)
            gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
            gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
            gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

        gmsh.model.mesh.generate(3)

        tag_nodos, coords, _ = gmsh.model.mesh.getNodes()
        tag_nodos = np.asarray(tag_nodos, dtype=np.int64)
        order = np.argsort(tag_nodos)
        tag_nodos = tag_nodos[order]
        nodes = np.asarray(coords, dtype=float).reshape(-1, 3)[order]
        indice = {int(t): i for i, t in enumerate(tag_nodos)}

        tipos, etiquetas, nodos_por_el = gmsh.model.mesh.getElements(3)
        tets_crudos = []
        for i, tipo in enumerate(tipos):     # `tipos` es un ndarray, no una lista
            if int(tipo) != 4:               # 4 = tetraedro de 4 nodos
                continue
            tets_crudos.append(np.asarray(nodos_por_el[i], dtype=np.int64).reshape(-1, 4))
        tets = np.vstack(tets_crudos) if tets_crudos else np.zeros((0, 4), dtype=np.int64)
        tets = np.vectorize(indice.__getitem__)(tets).astype(np.int64)

        # separar por región usando el centroide: la esfera es r < radius
        centros = nodes[tets].mean(axis=1)
        dentro = np.linalg.norm(centros, axis=1) < radius
        tets_solid = tets[dentro]
        tets_fluid = tets[~dentro]

        caras_solido, _ = _faces_once(_tet_faces(tets_solid))
        caras_fluido, _ = _faces_once(_tet_faces(tets_fluid))
        # la interfaz son las caras de frontera del sólido que también son frontera del fluido
        clave_s = {tuple(sorted(f)) for f in caras_solido}
        interfaz = np.array([f for f in caras_fluido if tuple(sorted(f)) in clave_s], dtype=np.int64)

        if interfaz.size == 0:
            raise RuntimeError("la malla no tiene interfaz compartida: revisar el fragment de Gmsh")

        # orientar con la normal saliente del sólido, para que el signo del acoplamiento sea fijo
        p = nodes[interfaz]
        n = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
        centros_cara = p.mean(axis=1)
        afuera = np.einsum("ij,ij->i", n, centros_cara) < 0
        interfaz[afuera] = interfaz[afuera][:, [0, 2, 1]]

        return SphereShellMesh(
            nodes=nodes, tets_solid=tets_solid, tets_fluid=tets_fluid,
            interface_faces=interfaz, radius=float(radius), outer_radius=float(outer_radius),
        )
    finally:
        gmsh.finalize()


def _surface_is_sphere(gmsh, dim: int, tag: int, radius: float, tol: float = 1e-6) -> bool:
    """Si una superficie entera está sobre ``r = radius`` (por su caja envolvente)."""
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
    extension = max(xmax - xmin, ymax - ymin, zmax - zmin) / 2.0
    return abs(extension - radius) <= tol * max(radius, 1e-30)
