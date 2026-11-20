"""Mallador de cajas estructuradas en tetraedros y utilidades de frontera 3D (Fase 6a).

Todo lo que devuelve este módulo es geometría pura: no sabe de física. Es el equivalente 3D de
`meshing.py` y existe para que el núcleo acústico 3D tenga una malla verificable sin mallador
externo (Gmsh llega en la Fase 6b, con la geometría de CT).

Construcción
------------
Una caja ``[xmin, xmax] × [ymin, ymax] × [zmin, zmax]`` se parte en ``nx × ny × nz`` celdas
hexaédricas y cada celda en **6 tetraedros** por la subdivisión de Kuhn (Freudenthal): los seis
caminos que van del vértice ``000`` al ``111`` sumando un eje por vez. Las diagonales de cada cara
del cubo coinciden con las de la celda vecina (todas usan la diagonal ``000–111`` proyectada sobre
la cara), así que la malla es **conforme**: no quedan caras colgantes ni nodos huérfanos, que es
justo lo que rompe el ensamblado P1.

Orientación
-----------
Cada tetraedro se reordena para que su volumen con signo sea **positivo**
(``det[v₁−v₀, v₂−v₀, v₃−v₀] > 0``). De eso depende el signo de la normal saliente de
:func:`boundary_faces` y, a través suyo, el signo de todos los términos de borde del FEM. La
orientación se verifica en la construcción y en :func:`check_tet_mesh`, no se supone.

Frontera
--------
:func:`boundary_faces` devuelve las caras (triángulos) que pertenecen a un único tetraedro, con la
orientación del tetraedro dueño, que para un tetraedro positivo es la que deja la normal
``(p₁−p₀)×(p₂−p₀)`` **saliente** del dominio. :func:`boundary_edges` devuelve las aristas de esa
superficie (cada una una sola vez).

Unidades: SI (m). Los espaciados se piden en metros o como número de celdas por eje.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass

import numpy as np

__all__ = [
    "TetMesh",
    "tet_signed_volumes",
    "check_tet_mesh",
    "box_tet_mesh",
    "box_tet_mesh_from_spacing",
    "boundary_faces",
    "face_areas",
    "face_normals",
    "face_centroids",
    "boundary_edges",
    "faces_on_plane",
    "unique_rows_appearing_once",
]

#: Los seis caminos ``000 → 111`` de la subdivisión de Kuhn, en índices locales ``dx + 2·dy + 4·dz``.
_KUHN_PATHS = np.array(list(itertools.permutations([1, 2, 4])), dtype=np.int64)

#: Caras locales de un tetraedro positivo, orientadas con la normal saliente.
_LOCAL_FACES = np.array([[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]], dtype=np.int64)


@dataclass(frozen=True)
class TetMesh:
    """Malla 3D: ``nodes`` ``(N, 3)`` y ``tets`` ``(M, 4)`` con orientación positiva.

    Atributos
    ---------
    nodes : coordenadas de los nodos [m], forma ``(N, 3)``.
    tets : conectividad, forma ``(M, 4)``; ``tets[i]`` son los cuatro nodos del tetraedro ``i``
        en orden tal que su volumen con signo es positivo.
    """

    nodes: np.ndarray
    tets: np.ndarray

    @property
    def n_nodes(self) -> int:
        return int(self.nodes.shape[0])

    @property
    def n_elements(self) -> int:
        return int(self.tets.shape[0])

    def signed_volumes(self) -> np.ndarray:
        """Volumen con signo de cada tetraedro [m³] (positivo en una malla bien orientada)."""
        return tet_signed_volumes(self.nodes, self.tets)

    def volumes(self) -> np.ndarray:
        """Volumen (positivo) de cada tetraedro [m³]."""
        return np.abs(self.signed_volumes())

    def centroids(self) -> np.ndarray:
        """Centroide de cada tetraedro [m], forma ``(M, 3)``."""
        return self.nodes[self.tets].mean(axis=1)

    def boundary_faces(self) -> np.ndarray:
        """Caras de frontera orientadas (ver :func:`boundary_faces`)."""
        return boundary_faces(self.tets)


# --------------------------------------------------------------------------- #
# Geometría de elemento
# --------------------------------------------------------------------------- #


def tet_signed_volumes(nodes, tets) -> np.ndarray:
    """Volumen con signo de cada tetraedro [m³].

    ``V = det[v₁−v₀, v₂−v₀, v₃−v₀] / 6``; positivo si ``(v₁−v₀, v₂−v₀, v₃−v₀)`` es dextrógira.
    """
    x = np.asarray(nodes, dtype=float)
    t = np.asarray(tets, dtype=np.int64)
    p = x[t]
    matrix = np.stack([p[:, 1] - p[:, 0], p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]], axis=-1)
    return np.linalg.det(matrix) / 6.0


def check_tet_mesh(nodes, tets) -> None:
    """Valida forma, finitud e **orientación positiva** de una malla tetraédrica.

    Lanza ``ValueError`` explicando qué está mal: forma incorrecta, índices fuera de rango, nodos
    repetidos dentro de un tetraedro, valores no finitos o volumen con signo no positivo (que es
    a la vez un tetraedro degenerado y una orientación invertida).
    """
    x = np.asarray(nodes, dtype=float)
    t = np.asarray(tets)

    if x.ndim != 2 or x.shape[1] != 3:
        raise ValueError(f"nodes debe tener forma (N, 3); recibida {x.shape}")
    if x.shape[0] < 4:
        raise ValueError(f"se requieren al menos 4 nodos; recibidos {x.shape[0]}")
    if not np.all(np.isfinite(x)):
        raise ValueError("nodes contiene valores no finitos (NaN o inf)")

    if t.ndim != 2 or t.shape[1] != 4:
        raise ValueError(f"tets debe tener forma (M, 4); recibida {t.shape}")
    if t.shape[0] < 1:
        raise ValueError("tets está vacío")
    if not np.issubdtype(t.dtype, np.integer):
        raise ValueError(f"tets debe ser entero; recibido dtype {t.dtype}")
    if t.min() < 0 or t.max() >= x.shape[0]:
        raise ValueError(
            f"tets tiene índices fuera de rango: mínimo {int(t.min())}, máximo {int(t.max())}, "
            f"nodos disponibles {x.shape[0]}"
        )
    repeated = (t[:, 0] == t[:, 1]) | (t[:, 1] == t[:, 2]) | (t[:, 2] == t[:, 3]) | (
        t[:, 0] == t[:, 2]
    ) | (t[:, 0] == t[:, 3]) | (t[:, 1] == t[:, 3])
    if np.any(repeated):
        i = int(np.flatnonzero(repeated)[0])
        raise ValueError(f"el tetraedro {i} repite nodos: {t[i].tolist()}")

    signed = tet_signed_volumes(x, t)
    bad = signed <= 0.0
    if np.any(bad):
        i = int(np.flatnonzero(bad)[0])
        raise ValueError(
            f"el tetraedro {i} = {t[i].tolist()} tiene volumen con signo {signed[i]:.3e} m³ "
            "(≤ 0): está degenerado o mal orientado. Los tets deben ser positivos "
            "(v₁−v₀, v₂−v₀, v₃−v₀ dextrógiras)"
        )


# --------------------------------------------------------------------------- #
# Mallado de cajas
# --------------------------------------------------------------------------- #


def _box_tet_mesh(xs: np.ndarray, ys: np.ndarray, zs: np.ndarray) -> TetMesh:
    """Malla de una caja a partir de las coordenadas de sus planos de nodos (por eje)."""
    nx, ny, nz = xs.size - 1, ys.size - 1, zs.size - 1
    if min(nx, ny, nz) < 1:
        raise ValueError("se requiere al menos una celda por eje")

    gx, gy, gz = np.meshgrid(xs, ys, zs, indexing="ij")
    nodes = np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()]).astype(float)

    sx = (ny + 1) * (nz + 1)
    sy = nz + 1
    i, j, k = np.meshgrid(
        np.arange(nx, dtype=np.int64),
        np.arange(ny, dtype=np.int64),
        np.arange(nz, dtype=np.int64),
        indexing="ij",
    )
    first = i * sx + j * sy + k  # nodo (0,0,0) de cada celda

    local_bit = np.arange(8, dtype=np.int64)  # dx + 2·dy + 4·dz
    off = (local_bit % 2) * sx + ((local_bit // 2) % 2) * sy + (local_bit // 4)
    corners = first.reshape(-1, 1) + off.reshape(1, -1)  # (n_cells, 8)

    # (n_cells, 6, 4): cada camino 000 → p1 → p1+p2 → 111
    tets = np.empty((corners.shape[0], 6, 4), dtype=np.int64)
    tets[:, :, 0] = corners[:, :1]
    tets[:, :, 1] = corners[:, _KUHN_PATHS[:, 0]]
    tets[:, :, 2] = corners[:, _KUHN_PATHS[:, 0] + _KUHN_PATHS[:, 1]]
    tets[:, :, 3] = corners[:, 7:8]
    tets = tets.reshape(-1, 4)

    # Orientación positiva (los seis caminos no la comparten: el determinante cambia de signo).
    signed = tet_signed_volumes(nodes, tets)
    negative = signed < 0.0
    if np.any(negative):
        tets[negative] = tets[negative][:, [0, 1, 3, 2]]
    if np.any(tet_signed_volumes(nodes, tets) <= 0.0):
        raise RuntimeError("malla degenerada: algún tetraedro quedó con volumen nulo")

    return TetMesh(nodes=nodes, tets=tets)


def box_tet_mesh(
    xmin: float, xmax: float, ymin: float, ymax: float, zmin: float, zmax: float,
    nx: int, ny: int, nz: int,
) -> TetMesh:
    """Malla tetraédrica de una caja con ``nx × ny × nz`` celdas (espaciado propio por eje).

    Parameters
    ----------
    xmin, xmax, ymin, ymax, zmin, zmax : límites de la caja [m].
    nx, ny, nz : número de celdas por eje (``>= 1``). Cada celda da 6 tetraedros, así que la malla
        tiene ``6·nx·ny·nz`` elementos y ``(nx+1)(ny+1)(nz+1)`` nodos.

    Returns
    -------
    :class:`TetMesh` con orientación positiva verificada.
    """
    bounds = np.array([xmin, xmax, ymin, ymax, zmin, zmax], dtype=float)
    if not np.all(np.isfinite(bounds)):
        raise ValueError("los límites de la caja deben ser finitos")
    if not (xmax > xmin and ymax > ymin and zmax > zmin):
        raise ValueError(
            f"se requiere xmax > xmin, ymax > ymin, zmax > zmin; recibido "
            f"x=[{xmin}, {xmax}], y=[{ymin}, {ymax}], z=[{zmin}, {zmax}]"
        )
    counts = (nx, ny, nz)
    if any(int(c) < 1 for c in counts) or any(int(c) != c for c in counts):
        raise ValueError(f"nx, ny, nz deben ser enteros >= 1; recibido {counts}")

    xs = np.linspace(float(xmin), float(xmax), int(nx) + 1)
    ys = np.linspace(float(ymin), float(ymax), int(ny) + 1)
    zs = np.linspace(float(zmin), float(zmax), int(nz) + 1)
    return _box_tet_mesh(xs, ys, zs)


def box_tet_mesh_from_spacing(lower, upper, spacing) -> TetMesh:
    """Malla de una caja con espaciado objetivo ``h`` [m] por eje.

    ``nx = ceil(Lx/hx)`` y análogos, así que el espaciado real es ``<= spacing`` en cada eje.
    Es la forma natural de pedir "ppw puntos por longitud de onda":
    ``spacing = λ / ppw``.

    Parameters
    ----------
    lower, upper : esquinas opuestas de la caja [m], cada una de longitud 3.
    spacing : paso objetivo [m]; escalar (isótropo) o secuencia de tres ``(hx, hy, hz)``.
    """
    lo = np.asarray(lower, dtype=float).reshape(-1)
    hi = np.asarray(upper, dtype=float).reshape(-1)
    if lo.size != 3 or hi.size != 3:
        raise ValueError(f"lower y upper deben ser de longitud 3; recibido {lo.size} y {hi.size}")
    if not (np.all(np.isfinite(lo)) and np.all(np.isfinite(hi))):
        raise ValueError("lower y upper deben ser finitos")
    lengths = hi - lo
    if np.any(lengths <= 0.0):
        raise ValueError(f"se requiere upper > lower en cada eje; longitudes {lengths.tolist()}")

    h = np.asarray(spacing, dtype=float).reshape(-1)
    if h.size == 1:
        h = np.full(3, h[0])
    if h.size != 3:
        raise ValueError(f"spacing debe ser escalar o de longitud 3; recibido {h.size}")
    if np.any(~np.isfinite(h)) or np.any(h <= 0.0):
        raise ValueError(f"spacing debe ser finito y > 0; recibido {h.tolist()}")

    counts = np.maximum(1, np.ceil(lengths / h)).astype(int)
    mesh = box_tet_mesh(
        lo[0], hi[0], lo[1], hi[1], lo[2], hi[2], int(counts[0]), int(counts[1]), int(counts[2])
    )
    return mesh


# --------------------------------------------------------------------------- #
# Frontera
# --------------------------------------------------------------------------- #


def unique_rows_appearing_once(rows: np.ndarray, n_nodes: int) -> tuple[np.ndarray, np.ndarray]:
    """Filas (de 3 columnas) que aparecen una sola vez, y su conteo, en orden de aparición.

    Se codifica cada fila ordenada como un único entero (``n0 + N·n1 + N²·n2`` con
    ``N = n_nodes + 1``) y se busca en un ``argsort`` 1D en vez de un ``np.unique(axis=0)``, que
    sobre millones de caras es mucho más caro en tiempo y memoria.

    Returns
    -------
    (mask, counts) : ``mask`` booleano ``(F,)`` con las filas que aparecen una vez y ``counts``
        ``(F,)`` con el número de apariciones de cada fila.
    """
    r = np.asarray(rows, dtype=np.int64)
    if r.ndim != 2 or r.shape[1] != 3:
        raise ValueError(f"rows debe tener forma (F, 3); recibida {r.shape}")
    n = int(n_nodes) + 1
    if n**3 > 2**62:
        # Camino lento pero sin desbordar: sólo se alcanza con > 2·10⁶ nodos.
        sorted_rows = np.sort(r, axis=1)
        _, inverse, counts = np.unique(sorted_rows, axis=0, return_inverse=True,
                                       return_counts=True)
        return counts[inverse] == 1, counts[inverse]

    sorted_rows = np.sort(r, axis=1)
    key = sorted_rows[:, 0] + n * (sorted_rows[:, 1] + n * sorted_rows[:, 2])
    order = np.argsort(key, kind="stable")
    sk = key[order]
    counts_sorted = np.ones(sk.size, dtype=np.int64)
    same = np.empty(sk.size, dtype=bool)
    same[0] = False
    same[1:] = sk[1:] == sk[:-1]
    # arranque de cada grupo: cuenta el tamaño del grupo
    starts = np.flatnonzero(~same)
    sizes = np.diff(np.append(starts, sk.size))
    counts_sorted[:] = np.repeat(sizes, sizes)
    counts = np.empty_like(counts_sorted)
    counts[order] = counts_sorted
    return counts == 1, counts


def boundary_faces(tets) -> np.ndarray:
    """Caras de frontera de una malla tetraédrica, orientadas hacia afuera.

    Cada tetraedro aporta sus cuatro caras con la orientación que, para un tetraedro **positivo**,
    deja la normal ``(p₁−p₀)×(p₂−p₀)`` saliente del elemento. Las caras interiores aparecen dos
    veces (una por cada tetraedro dueño, con orientaciones opuestas) y se descartan; las de
    frontera aparecen una sola vez y conservan la orientación saliente.

    Returns
    -------
    ``(F, 3)`` índices de nodo de las caras de frontera [—].
    """
    t = np.asarray(tets, dtype=np.int64)
    if t.ndim != 2 or t.shape[1] != 4:
        raise ValueError(f"tets debe tener forma (M, 4); recibida {t.shape}")
    all_faces = t[:, _LOCAL_FACES].reshape(-1, 3)
    n_nodes = int(t.max()) + 1
    once, _ = unique_rows_appearing_once(all_faces, n_nodes)
    return all_faces[once]


def face_areas(nodes, faces) -> np.ndarray:
    """Área [m²] de cada cara triangular."""
    x = np.asarray(nodes, dtype=float)
    f = np.asarray(faces, dtype=np.int64)
    a, b, c = x[f[:, 0]], x[f[:, 1]], x[f[:, 2]]
    return 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)


def face_normals(nodes, faces) -> np.ndarray:
    """Normal unitaria de cada cara [—], según el orden de sus vértices.

    Para las caras que devuelve :func:`boundary_faces` la normal es **saliente** del dominio.
    """
    x = np.asarray(nodes, dtype=float)
    f = np.asarray(faces, dtype=np.int64)
    a, b, c = x[f[:, 0]], x[f[:, 1]], x[f[:, 2]]
    n = np.cross(b - a, c - a)
    norm = np.linalg.norm(n, axis=1, keepdims=True)
    if np.any(norm == 0.0):
        i = int(np.flatnonzero(norm[:, 0] == 0.0)[0])
        raise ValueError(f"la cara {i} es degenerada (área nula): no tiene normal")
    return n / norm


def face_centroids(nodes, faces) -> np.ndarray:
    """Centroide de cada cara [m], forma ``(F, 3)``."""
    x = np.asarray(nodes, dtype=float)
    f = np.asarray(faces, dtype=np.int64)
    return x[f].mean(axis=1)


def boundary_edges(faces) -> np.ndarray:
    """Aristas de la superficie de frontera, sin repetir y ordenadas por fila.

    Devuelve ``(E, 2)`` con cada arista una sola vez (``e[0] < e[1]``). Es la conectividad de la
    superficie: sirve para cadenas, para graficar la frontera y para auditar que las caras forman
    una superficie cerrada (cada arista de una superficie cerrada pertenece a exactamente dos
    caras).
    """
    f = np.asarray(faces, dtype=np.int64)
    if f.ndim != 2 or f.shape[1] != 3:
        raise ValueError(f"faces debe tener forma (F, 3); recibida {f.shape}")
    edges = np.vstack([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]])
    edges = np.sort(edges, axis=1)
    n_nodes = int(edges.max()) + 1
    # clave 1D para deduplicar aristas (2 columnas: n0 + N·n1)
    key = edges[:, 0] + n_nodes * edges[:, 1]
    _, index = np.unique(key, return_index=True)
    return edges[np.sort(index)]


def faces_on_plane(nodes, faces, axis: int, coordinate: float, tol: float = 1e-12) -> np.ndarray:
    """Máscara de las caras cuyos **tres** vértices están sobre el plano ``x_axis = coordinate``.

    Es el selector que hace falta para armar la lista de caras de la ABC (o de un dato de
    Neumann) en una caja.

    Parameters
    ----------
    axis : ``0``, ``1`` o ``2`` para ``x``, ``y`` o ``z``.
    coordinate : posición del plano [m].
    tol : tolerancia absoluta [m] con la que se compara la coordenada.
    """
    x = np.asarray(nodes, dtype=float)
    f = np.asarray(faces, dtype=np.int64)
    if axis not in (0, 1, 2):
        raise ValueError(f"axis debe ser 0, 1 o 2; recibido {axis!r}")
    on = np.abs(x[f][:, :, int(axis)] - float(coordinate)) <= float(tol)
    return np.all(on, axis=1)
