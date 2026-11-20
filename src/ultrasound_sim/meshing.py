"""Mallas simples y deterministas (1D y 2D trianguladas), sin mallador externo.

Todo lo que devuelve este módulo es geometría pura: no sabe de física. Las mallas 2D son
estructuradas (quads partidos en dos triángulos), con orientación **CCW** garantizada, porque de
eso dependen los signos de las integrales de borde en `fem/` y `bem/`.

Para geometrías reales (CT, cóclea) la malla vendrá de Gmsh en fases posteriores; estas funciones
existen para las fases 0–4 y para los tests.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class LineMesh:
    """Malla 1D: ``nodes`` (n+1,) y ``elements`` (n, 2) con ``elements[i] = (i, i+1)``."""

    nodes: np.ndarray
    elements: np.ndarray

    @property
    def n_elements(self) -> int:
        return int(self.elements.shape[0])

    @property
    def lengths(self) -> np.ndarray:
        return np.abs(np.diff(self.nodes))


@dataclass(frozen=True)
class TriMesh:
    """Malla 2D: ``nodes`` (N, 2) y ``triangles`` (M, 3) con orientación CCW."""

    nodes: np.ndarray
    triangles: np.ndarray

    @property
    def n_nodes(self) -> int:
        return int(self.nodes.shape[0])

    @property
    def n_elements(self) -> int:
        return int(self.triangles.shape[0])

    def areas(self) -> np.ndarray:
        return triangle_areas(self.nodes, self.triangles)

    def centroids(self) -> np.ndarray:
        return self.nodes[self.triangles].mean(axis=1)


def line_mesh(x0: float, x1: float, n_elements: int) -> LineMesh:
    """Malla 1D uniforme de ``n_elements`` tramos entre ``x0`` y ``x1``."""
    if n_elements < 1:
        raise ValueError("n_elements debe ser >= 1")
    if not x1 > x0:
        raise ValueError("se requiere x1 > x0")
    nodes = np.linspace(float(x0), float(x1), n_elements + 1)
    elements = np.column_stack([np.arange(n_elements), np.arange(1, n_elements + 1)])
    return LineMesh(nodes=nodes, elements=elements.astype(int))


def layered_line_mesh(interfaces: list[float] | np.ndarray, max_step: float) -> LineMesh:
    """Malla 1D con nodos forzados en cada interfaz de material.

    Cada tramo entre interfaces consecutivas se subdivide en partes iguales de paso ``<= max_step``.
    Imprescindible en FEM con medios por capas: un elemento que cruza una interfaz promedia dos
    materiales y arruina el contraste.
    """
    xs = np.asarray(interfaces, dtype=float)
    if xs.ndim != 1 or xs.size < 2:
        raise ValueError("se requieren al menos dos interfaces")
    if np.any(np.diff(xs) <= 0.0):
        raise ValueError("las interfaces deben ser estrictamente crecientes")
    if not max_step > 0.0:
        raise ValueError("max_step debe ser > 0")

    nodes = [xs[0]]
    for xa, xb in zip(xs[:-1], xs[1:]):
        n = max(1, int(np.ceil((xb - xa) / max_step)))
        nodes.extend(np.linspace(xa, xb, n + 1)[1:].tolist())
    nodes_arr = np.asarray(nodes, dtype=float)
    elements = np.column_stack(
        [np.arange(nodes_arr.size - 1), np.arange(1, nodes_arr.size)]
    )
    return LineMesh(nodes=nodes_arr, elements=elements.astype(int))


def rect_tri_mesh(
    xmin: float, xmax: float, ymin: float, ymax: float, nx: int, ny: int
) -> TriMesh:
    """Malla de un rectángulo con ``nx × ny`` quads, cada uno partido en 2 triángulos CCW."""
    if nx < 1 or ny < 1:
        raise ValueError("nx y ny deben ser >= 1")
    if not (xmax > xmin and ymax > ymin):
        raise ValueError("se requiere xmax > xmin y ymax > ymin")

    xs = np.linspace(float(xmin), float(xmax), nx + 1)
    ys = np.linspace(float(ymin), float(ymax), ny + 1)
    grid_x, grid_y = np.meshgrid(xs, ys, indexing="xy")
    nodes = np.column_stack([grid_x.ravel(), grid_y.ravel()])

    def nid(i: int, j: int) -> int:  # i: columna, j: fila
        return j * (nx + 1) + i

    triangles: list[tuple[int, int, int]] = []
    for j in range(ny):
        for i in range(nx):
            a, b = nid(i, j), nid(i + 1, j)
            c, d = nid(i + 1, j + 1), nid(i, j + 1)
            triangles.append((a, b, c))
            triangles.append((a, c, d))
    tris = np.asarray(triangles, dtype=int)
    return TriMesh(nodes=nodes, triangles=_ensure_ccw(nodes, tris))


def polar_annulus_mesh(
    inner_radius: float, outer_radius: float, n_radial: int, n_theta: int
) -> TriMesh:
    """Malla de una corona circular ``inner_radius <= r <= outer_radius``.

    Es la malla que pide el contraste FEM/BEM de la Fase 4: el FEM resuelve el anillo
    ``a <= r <= R`` (Neumann en ``r = a``, ABC en ``r = R``) y el BEM resuelve el exterior del
    cilindro de radio ``a``; ambos comparan el mismo campo en puntos del anillo.

    ``n_radial`` tramos radiales y ``n_theta`` tramos angulares (periódicos: el nodo ``j = n_theta``
    es el mismo que ``j = 0``, así que no se duplica y la malla no tiene costura).
    """
    if not (0.0 < inner_radius < outer_radius):
        raise ValueError("se requiere 0 < inner_radius < outer_radius")
    if n_radial < 1 or n_theta < 3:
        raise ValueError("se requiere n_radial >= 1 y n_theta >= 3")

    radii = np.linspace(float(inner_radius), float(outer_radius), n_radial + 1)
    theta = 2.0 * np.pi * np.arange(n_theta) / n_theta
    rr, tt = np.meshgrid(radii, theta, indexing="ij")  # (n_radial+1, n_theta)
    nodes = np.column_stack([(rr * np.cos(tt)).ravel(), (rr * np.sin(tt)).ravel()])

    def nid(i: int, j: int) -> int:
        return i * n_theta + (j % n_theta)

    triangles: list[tuple[int, int, int]] = []
    for i in range(n_radial):
        for j in range(n_theta):
            a, b = nid(i, j), nid(i, j + 1)
            c, d = nid(i + 1, j + 1), nid(i + 1, j)
            triangles.append((a, b, c))
            triangles.append((a, c, d))
    tris = np.asarray(triangles, dtype=int)
    return TriMesh(nodes=nodes, triangles=_ensure_ccw(nodes, tris))


def triangle_areas(nodes: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    """Área con signo de cada triángulo (positiva si es CCW)."""
    p0 = nodes[triangles[:, 0]]
    p1 = nodes[triangles[:, 1]]
    p2 = nodes[triangles[:, 2]]
    return 0.5 * (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p2[:, 0] - p0[:, 0]) * (p1[:, 1] - p0[:, 1])
    )


def _ensure_ccw(nodes: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    """Voltea los triángulos con área negativa para que todos sean CCW."""
    tris = triangles.copy()
    negative = triangle_areas(nodes, tris) < 0.0
    if np.any(negative):
        tris[negative] = tris[negative][:, [0, 2, 1]]
    return tris


def boundary_edges(triangles: np.ndarray) -> np.ndarray:
    """Aristas de frontera de una malla triangular, orientadas según el triángulo CCW.

    Cada fila ``(i, j)`` recorre la arista en el sentido del triángulo dueño; para un triángulo
    CCW el normal **saliente** es ``(t_y, −t_x)`` con ``t = p_j − p_i``. Las aristas interiores
    (compartidas por dos triángulos) no aparecen.
    """
    edges = np.vstack(
        [
            triangles[:, [0, 1]],
            triangles[:, [1, 2]],
            triangles[:, [2, 0]],
        ]
    )
    key = np.sort(edges, axis=1)
    _, inverse, counts = np.unique(
        key, axis=0, return_inverse=True, return_counts=True
    )
    on_boundary = counts[inverse] == 1
    return edges[on_boundary]


def outward_normals(nodes: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """Normal unitaria saliente de cada arista de frontera (malla CCW)."""
    p0 = nodes[edges[:, 0]]
    p1 = nodes[edges[:, 1]]
    tangent = p1 - p0
    normal = np.column_stack([tangent[:, 1], -tangent[:, 0]])
    norm = np.linalg.norm(normal, axis=1, keepdims=True)
    if np.any(norm == 0.0):
        raise ValueError("arista de longitud nula en la frontera")
    return normal / norm


def edge_midpoints(nodes: np.ndarray, edges: np.ndarray) -> np.ndarray:
    return 0.5 * (nodes[edges[:, 0]] + nodes[edges[:, 1]])
