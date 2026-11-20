"""Ensamblado de elementos finitos P1 (lineales y triangulares), numpy puro.

Convenciones (comunes a todo `fem/`):

- La incógnita es la presión ``p`` en los nodos; las funciones de forma P1 son lineales por
  elemento y continuas en los nodos (nodal, no conforme para el gradiente).
- Mallas 1D: ``nodes`` ``(N,)`` y ``elements`` ``(m, 2)``.
- Mallas 2D: ``nodes`` ``(N, 2)`` y ``triangles`` ``(m, 3)`` con orientación **CCW**
  (la garantiza `meshing.py`); las áreas se toman positivas.
- Todo el ensamblado usa ``np.add.at`` sobre índices de grado de libertad, así que funciona
  igual con mallas no uniformes ni estructuradas.
- Las matrices locales se exponen por separado (``p1_line_local``, ``p1_tri_local``) para que
  los ensambladores con coeficientes variables (p. ej. ``k`` por elemento) las reutilicen sin
  duplicar fórmulas.

Términos que se construyen aquí:

- ``M``: masa (``∫ u·v``), útil para normas L2 y para el término ``−k²·M`` de Helmholtz.
- ``K``: rigidez (``∫∇u·∇v``).
- ``S``: masa de frontera (``∫_Γ u·v·ds``) para el término de borde/ABC.
"""

from __future__ import annotations

import numpy as np

from ..meshing import triangle_areas

__all__ = [
    "line_element_lengths",
    "triangle_element_areas",
    "triangle_gradients",
    "p1_line_local",
    "p1_tri_local",
    "p1_line_matrices",
    "p1_tri_matrices",
    "p1_boundary_mass",
    "weights_to_local",
]


# ---------------------------------------------------------------------------
# Geometría de elemento
# ---------------------------------------------------------------------------


def line_element_lengths(nodes, elements) -> np.ndarray:
    """Longitud (positiva) de cada elemento de línea."""
    x = np.asarray(nodes, dtype=float)
    e = np.asarray(elements, dtype=int)
    return np.abs(x[e[:, 1]] - x[e[:, 0]])


def triangle_element_areas(nodes, triangles) -> np.ndarray:
    """Área (positiva) de cada triángulo."""
    return np.abs(triangle_areas(np.asarray(nodes, dtype=float), np.asarray(triangles, dtype=int)))


def triangle_gradients(nodes, triangles) -> np.ndarray:
    """Gradientes (constantes por elemento) de las tres funciones de forma P1.

    Devuelve ``(m, 3, 2)``: ``grad[:, i, :]`` es ``∇φ_i`` en el elemento ``i``.
    Fórmula estándar ``∇φ_i = (b_i, c_i) / (2A)`` con
    ``b_i = y_{i+1} − y_{i+2}`` y ``c_i = x_{i+2} − x_{i+1}`` (índices cíclicos).
    """
    xyz = np.asarray(nodes, dtype=float)[np.asarray(triangles, dtype=int)]  # (m, 3, 2)
    y = xyz[:, :, 1]
    x = xyz[:, :, 0]
    b = np.column_stack([y[:, 1] - y[:, 2], y[:, 2] - y[:, 0], y[:, 0] - y[:, 1]])
    c = np.column_stack([x[:, 2] - x[:, 1], x[:, 0] - x[:, 2], x[:, 1] - x[:, 0]])
    area = np.abs(triangle_areas(np.asarray(nodes, dtype=float), np.asarray(triangles, dtype=int)))
    two_a = 2.0 * area
    if np.any(two_a <= 0.0):
        raise ValueError("triángulo degenerado (área nula)")
    return np.stack([b, c], axis=-1) / two_a[:, None, None]


# ---------------------------------------------------------------------------
# Matrices locales
# ---------------------------------------------------------------------------


def p1_line_local(nodes, elements) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Matrices locales 1D: ``(M_e, K_e, h)`` con formas ``(m,2,2)``, ``(m,2,2)``, ``(m,)``.

    ``M_e = h/6·[[2,1],[1,2]]``, ``K_e = 1/h·[[1,−1],[−1,1]]``.
    """
    h = line_element_lengths(nodes, elements)
    m = h.size
    mass = np.zeros((m, 2, 2))
    mass[:, 0, 0] = 2.0 * h / 6.0
    mass[:, 0, 1] = h / 6.0
    mass[:, 1, 0] = h / 6.0
    mass[:, 1, 1] = 2.0 * h / 6.0
    stiff = np.zeros((m, 2, 2))
    stiff[:, 0, 0] = 1.0 / h
    stiff[:, 0, 1] = -1.0 / h
    stiff[:, 1, 0] = -1.0 / h
    stiff[:, 1, 1] = 1.0 / h
    return mass, stiff, h


def p1_tri_local(nodes, triangles) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Matrices locales 2D: ``(M_e, K_e, A)`` con formas ``(m,3,3)``, ``(m,3,3)``, ``(m,)``.

    ``M_e = A/12·(1 + I)`` (fila suma = A, consistente), ``K_e[i,j] = A·∇φ_i·∇φ_j``.
    """
    area = triangle_element_areas(nodes, triangles)
    grads = triangle_gradients(nodes, triangles)
    m = area.size
    mass = np.broadcast_to(area[:, None, None] / 12.0, (m, 3, 3)).copy()
    idx = np.arange(3)
    mass[:, idx, idx] += (area / 12.0)[:, None]  # diagonal doble: A/12·(1 + δ_ij)
    stiff = area[:, None, None] * np.einsum("mik,mjk->mij", grads, grads)
    return mass, stiff, area


# ---------------------------------------------------------------------------
# Ensamblado
# ---------------------------------------------------------------------------


def _accumulate(n_dofs: int, dofs: np.ndarray, local: np.ndarray) -> np.ndarray:
    """Suma matrices locales (``(m,n,n)``) en una matriz global densa ``(N,N)``."""
    n = dofs.shape[1]
    rows = np.repeat(dofs, n, axis=1).ravel()
    cols = np.tile(dofs, (1, n)).ravel()
    out = np.zeros((n_dofs, n_dofs), dtype=local.dtype)
    np.add.at(out, (rows, cols), np.asarray(local).reshape(-1))
    return out


def p1_line_matrices(nodes, elements) -> tuple[np.ndarray, np.ndarray]:
    """Matrices globales 1D ``(M, K)`` reales para elementos P1 de línea."""
    x = np.asarray(nodes, dtype=float)
    e = np.asarray(elements, dtype=int)
    mass, stiff, _ = p1_line_local(x, e)
    n = x.size
    return _accumulate(n, e, mass), _accumulate(n, e, stiff)


def p1_tri_matrices(nodes, triangles) -> tuple[np.ndarray, np.ndarray]:
    """Matrices globales 2D ``(M, K)`` reales para elementos P1 triangulares."""
    x = np.asarray(nodes, dtype=float)
    t = np.asarray(triangles, dtype=int)
    mass, stiff, _ = p1_tri_local(x, t)
    n = x.shape[0]
    return _accumulate(n, t, mass), _accumulate(n, t, stiff)


def p1_boundary_mass(nodes, edges) -> np.ndarray:
    """Matriz de masa 1D sobre aristas de frontera 2D: ``∫_e u·v·ds``.

    Sólo depende de la longitud de cada arista, de modo que sirve igual para aristas rectas
    de una malla triangular. Para el término ``−i·k·∫ p v ds`` se reemplaza por
    ``Σ_e k_e·S_e`` (ver `fem.helmholtz2d`).
    """
    x = np.asarray(nodes, dtype=float)
    e = np.asarray(edges, dtype=int)
    lengths = np.linalg.norm(x[e[:, 1]] - x[e[:, 0]], axis=1)
    mass, _, _ = p1_line_local(np.array([0.0, 1.0]), np.array([[0, 1]]))
    local = lengths[:, None, None] * mass[0][None, :, :]
    return _accumulate(x.shape[0], e, local)


def weights_to_local(lengths: np.ndarray) -> np.ndarray:
    """Matrices locales ``(m,2,2)`` de masa de borde a partir de longitudes."""
    m = lengths.size
    local = np.zeros((m, 2, 2))
    local[:, 0, 0] = 2.0 * lengths / 6.0
    local[:, 0, 1] = lengths / 6.0
    local[:, 1, 0] = lengths / 6.0
    local[:, 1, 1] = 2.0 * lengths / 6.0
    return local
