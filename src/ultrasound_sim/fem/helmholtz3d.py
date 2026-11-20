"""FEM P1 de Helmholtz 3D en tetraedros, con ABC-1 y BGT-2 sobre caras planas (Fase 6a).

Forma débil implementada (convención D4: ``e^{−iωt}``, ``∇²p + k²p = 0``, ``k = ω/c + i·α``)::

    ∫ (1/ρ)·(∇p·∇v − k² p v) dΩ + ∫_Γ (1/ρ)·[−c·p·v + d·∇_s p·∇_s v] dS
        = ∫_Γ (1/ρ)·(∂p_i/∂n − c_i·p_i)·v dS   (sólo caras iluminadas)

donde ``∂p/∂n = c·p + d·∇²_s p`` es la condición de frontera (``c = i·k``, ``d = 0`` para ABC-1;
ver :func:`abc_coefficients`) y ``∇_s`` es el gradiente tangencial sobre la cara. El peso ``1/ρ``
sólo entra cuando se pasa ``element_rho``; sin él (``ρ = 1``) la ecuación queda exactamente en la
forma ``∇p·∇v − k²·p·v`` del plan.

Por qué el término tangencial va con **+**
------------------------------------------
La condición se impone sustituyendo ``∂p/∂n`` en el término de borde de la integración por partes,
``+∮ (1/ρ)·∂p/∂n·v``, y pasándolo al lado izquierdo (es el mismo signo que en 2D, `fem.helmholtz2d`)::

    ∫_Γ (1/ρ)·(c·p + d·∇²_s p)·v dS = ∫_Γ (1/ρ)·c·p·v dS − ∫_Γ (1/ρ)·d·∇_s p·∇_s v dS

(la segunda igualdad es integración por partes **sobre la superficie**, cuyo término de contorno se
anula si la superficie es cerrada). De ahí que el bloque de masa de borde entre con signo ``−`` y el
de rigidez tangencial con signo ``+``, igual que en el 2D.

ABC de segundo orden (BGT) en 3D
--------------------------------
En 2D el repo usa la forma truncada de Bayliss–Turkel–Gunzburger para un contorno circular
(``i·k − 1/(2r) + i/(8k r²)`` más el término ``(i/2k)·∂²_s``). En 3D el operador aniquilador de dos
términos es ``B₂ = (∂_r − i·k + 3/r)(∂_r − i·k + 1/r)``, que se anula exactamente sobre
``p = e^{ikr}·(1/r + A/r²)`` para cualquier ``A`` (y por conmutación, sobre cualquier armónico
angular). Combinando ``B₂ p = 0`` con la ecuación de Helmholtz
``∂_r²p + (2/r)·∂_r p + k²p + ∇²_s p = 0`` (``∇²_s = Δ_S/r²``) y truncando en ``O(1/r²)`` queda::

    ∂p/∂n = (i·k − 1/r − i/(2·k·r²))·p + (i/(2·k))·∇²_s p

que es la condición implementada. El término ``−1/r`` es la curvatura del frente esférico (en 2D es
``−1/(2r)``, la del frente cilíndrico) y es el que arregla la incidencia oblicua: sin él, una onda
esférica que llega a 50° se refleja ~25 %. El error residual de la forma truncada es
``O(1/(k·r)³)``, medido en ``tests/test_helmholtz3d.py`` contra la onda esférica exacta.

``abc_order=2`` exige ``abc_center`` (el centro del frente esférico) y que las caras de la ABC
formen una **superficie cerrada**: si no, el término de contorno de la integración por partes sobre
la superficie no se anula y la condición deja de ser consistente (se rechaza con ``ValueError``).

Campo total vs campo dispersado
-------------------------------
:func:`assemble` resuelve el campo **total** con la ABC aplicada sobre él y la onda incidente
entrando por el término ``∂p_i/∂n − i·k·p_i`` en las caras **iluminadas** (``n·d < 0``). Vale en
cajas/ductos alineados con la onda (incidencia normal), que es el uso de la Fase 6a: ahí la ABC es
exacta y el error es sólo de discretización.

:func:`assemble_scattered` resuelve ``p_s = p − p_i``, saliente en **todo** el contorno (ABC
homogénea, sin máscara), con la onda incidente como dato de Neumann en la superficie del dispersor.
Es la formulación que hay que usar con contornos curvos o con obstáculos (D12, D21).

Notas de discretización
-----------------------
- La malla debe tener nodos sobre las interfaces de material: ``k`` y ``ρ`` son constantes por
  elemento, y un tetraedro que cruza una interfaz promedia dos medios. Con
  `meshing3d.box_tet_mesh` alcanza con que la interfaz caiga sobre un plano de nodos.
- El término de borde usa la masa P1 de la cara (``A/12·(1+I)``), que es la discretización de
  Galerkin consistente para ``k`` constante por cara.
- El lado derecho se integra con la regla del punto medio (exacta para ``v`` lineal sobre un
  triángulo plano).
- La fuente puntual es un delta **lumpado en el nodo más cercano**: el campo lejos de la fuente
  converge con ``h``, pero cerca de ella el valor nodal no representa el singular ``1/r``. Para
  validar contra ``e^{ikr}/(4πr)`` hay que medir a distancia (ver el test del criterio 3).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from ..materials import TWO_PI, Material
from ..meshing3d import (
    _LOCAL_FACES,
    boundary_edges,
    check_tet_mesh,
    face_areas,
    face_centroids,
    face_normals,
    tet_signed_volumes,
)

__all__ = [
    "PlaneWave3D",
    "PointSource",
    "plane_wave",
    "point_source",
    "tet_gradients",
    "p1_tet_local",
    "p1_tet_matrices",
    "face_mass",
    "surface_stiffness",
    "abc_coefficients",
    "assemble",
    "solve",
    "assemble_scattered",
    "solve_scattered",
    "total_field",
]


# --------------------------------------------------------------------------- #
# Ondas y fuentes
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class PlaneWave3D:
    """Onda plana ``p_i(x) = A·exp(i·k·d·x)`` con ``d`` unitaria.

    Atributos
    ---------
    k : número de onda complejo del medio [1/m].
    direction : dirección de propagación unitaria [—], forma ``(3,)``.
    amplitude : amplitud compleja [Pa].
    """

    k: complex
    direction: np.ndarray
    amplitude: complex = 1.0

    @property
    def wavenumber_vector(self) -> np.ndarray:
        """``k·d`` [1/m]."""
        return complex(self.k) * self.direction

    def value(self, points) -> np.ndarray:
        """Presión incidente [Pa] en ``points`` (``(M, 3)`` o ``(3,)``)."""
        pts = np.atleast_2d(np.asarray(points, dtype=float))
        return complex(self.amplitude) * np.exp(1j * (pts @ self.wavenumber_vector))

    # alias para leer el código de ensamblado igual que en 2D
    evaluate = value

    def gradient(self, points) -> np.ndarray:
        """``∇p_i = i·k·d·p_i`` [Pa/m] en ``points``, forma ``(M, 3)``."""
        pts = np.atleast_2d(np.asarray(points, dtype=float))
        return 1j * self.value(pts)[:, None] * self.wavenumber_vector[None, :]

    def normal_derivative(self, points, normals) -> np.ndarray:
        """``∂p_i/∂n = i·k·(d·n)·p_i`` [Pa/m] con normales unitarias."""
        pts = np.atleast_2d(np.asarray(points, dtype=float))
        nrm = np.atleast_2d(np.asarray(normals, dtype=float))
        return 1j * self.value(pts) * (nrm @ self.wavenumber_vector)


def plane_wave(frequency: float, material: Material, direction, amplitude: complex = 1.0
               ) -> PlaneWave3D:
    """Construye la onda plana incidente ``p_i = A·exp(i·k·d·x)``.

    Parameters
    ----------
    frequency : [Hz].
    material : medio por el que viaja la onda incidente.
    direction : vector 3D (no nulo); se normaliza.
    amplitude : amplitud compleja de presión [Pa].
    """
    if not isinstance(material, Material):
        raise TypeError("material debe ser un materials.Material")
    _positive_frequency(frequency)
    d = np.asarray(direction, dtype=float).reshape(-1)
    if d.size != 3:
        raise ValueError(f"direction debe ser un vector 3D; recibido tamaño {d.size}")
    if not np.all(np.isfinite(d)):
        raise ValueError("direction contiene valores no finitos")
    norm = float(np.linalg.norm(d))
    if norm == 0.0:
        raise ValueError("direction no puede ser nula")
    return PlaneWave3D(
        k=material.wavenumber(float(frequency)), direction=d / norm, amplitude=complex(amplitude)
    )


@dataclass(frozen=True)
class PointSource:
    """Fuente monopolar en espacio libre: ``p(x) = −i·ω·ρ·Q·e^{i·k·r}/(4π·r)``.

    El signo sale de la convención D4 (``e^{−iωt}``): con ``Q`` el caudal volumétrico inyectado
    [m³/s], la ecuación es ``∇·((1/ρ)∇p) + (k²/ρ)p = i·ω·Q·δ(x−x0)`` y su solución es la de arriba.
    Con ``Q`` real la presión va desfasada ``−90°`` respecto del caudal, que es lo correcto para un
    monopolo.

    Atributos
    ---------
    position : posición de la fuente [m], forma ``(3,)``.
    volume_velocity : caudal volumétrico complejo ``Q`` [m³/s].
    k : número de onda del medio [1/m].
    omega : pulsación ``2π·f`` [rad/s].
    rho : densidad del medio [kg/m³].
    """

    position: np.ndarray
    volume_velocity: complex
    k: complex
    omega: float
    rho: float

    def value(self, points) -> np.ndarray:
        """Presión [Pa] en ``points`` (``(M, 3)`` o ``(3,)``); diverge en la posición de la fuente."""
        pts = np.atleast_2d(np.asarray(points, dtype=float))
        r = np.linalg.norm(pts - self.position[None, :], axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            green = np.exp(1j * self.k * r) / (4.0 * np.pi * r)
        return -1j * self.omega * self.rho * self.volume_velocity * green

    def normal_derivative(self, points, normals) -> np.ndarray:
        """``∂p/∂n`` [Pa/m] con normales unitarias (``∂G/∂r = (i·k − 1/r)·G``)."""
        pts = np.atleast_2d(np.asarray(points, dtype=float))
        nrm = np.atleast_2d(np.asarray(normals, dtype=float))
        delta = pts - self.position[None, :]
        r = np.linalg.norm(delta, axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            green = np.exp(1j * self.k * r) / (4.0 * np.pi * r)
            dgreen_dr = (1j * self.k - 1.0 / r) * green
        return -1j * self.omega * self.rho * self.volume_velocity * dgreen_dr * (
            np.einsum("ij,ij->i", delta, nrm) / r
        )


def point_source(frequency: float, material: Material, position, volume_velocity: complex = 1.0
                 ) -> PointSource:
    """Construye la fuente puntual (monopolo) de referencia analítica.

    Parameters
    ----------
    frequency : [Hz].
    material : medio (de ahí salen ``k``, ``ρ`` y la atenuación).
    position : posición [m], forma ``(3,)``.
    volume_velocity : caudal volumétrico complejo [m³/s].
    """
    if not isinstance(material, Material):
        raise TypeError("material debe ser un materials.Material")
    f = _positive_frequency(frequency)
    pos = np.asarray(position, dtype=float).reshape(-1)
    if pos.size != 3:
        raise ValueError(f"position debe tener 3 componentes; recibido {pos.size}")
    if not np.all(np.isfinite(pos)):
        raise ValueError("position contiene valores no finitos")
    if not np.isfinite(np.abs(complex(volume_velocity))):
        raise ValueError("volume_velocity debe ser finito")
    return PointSource(
        position=pos,
        volume_velocity=complex(volume_velocity),
        k=material.wavenumber(f),
        omega=TWO_PI * f,
        rho=float(material.rho),
    )


# --------------------------------------------------------------------------- #
# Validación
# --------------------------------------------------------------------------- #


def _positive_frequency(frequency) -> float:
    try:
        f = float(frequency)
    except (TypeError, ValueError):
        raise ValueError(f"frequency debe ser un número [Hz]; recibido {frequency!r}") from None
    if not np.isfinite(f) or f <= 0.0:
        raise ValueError(f"frequency debe ser finita y > 0 [Hz]; recibido {f!r}")
    return f


def _per_element(value, n_elements: int, name: str) -> np.ndarray:
    """Normaliza un parámetro escalar o por elemento a un array ``(M,)``."""
    if value is None:
        return np.ones(n_elements, dtype=float)
    arr = np.asarray(value).reshape(-1)
    if arr.size == 1:
        return np.full(n_elements, arr[0], dtype=arr.dtype)
    if arr.size != n_elements:
        raise ValueError(f"{name} debe ser escalar o tener un valor por elemento ({n_elements}); "
                         f"recibido {arr.size}")
    return arr


def _element_k_array(element_k, n_elements: int) -> np.ndarray:
    k = np.asarray(element_k, dtype=complex).reshape(-1)
    if k.size == 1:
        k = np.full(n_elements, k[0], dtype=complex)
    if k.size != n_elements:
        raise ValueError(f"element_k debe ser escalar o tener un valor por elemento ({n_elements}); "
                         f"recibido {k.size}")
    if not (np.all(np.isfinite(k.real)) and np.all(np.isfinite(k.imag))):
        raise ValueError("element_k contiene valores no finitos (NaN o inf)")
    if np.any(k.real <= 0.0):
        raise ValueError("element_k debe tener parte real > 0 (número de onda de propagación)")
    if np.any(k.imag < 0.0):
        raise ValueError(
            "element_k debe tener parte imaginaria >= 0: con la convención e^{−iωt} la atenuación "
            "es k = ω/c + i·α (α > 0)"
        )
    return k


def _element_rho_array(element_rho, n_elements: int) -> np.ndarray:
    rho = _per_element(element_rho, n_elements, "element_rho").astype(float)
    if not np.all(np.isfinite(rho)):
        raise ValueError("element_rho contiene valores no finitos")
    if np.any(rho <= 0.0):
        raise ValueError("element_rho debe ser > 0 [kg/m³]")
    return rho


# --------------------------------------------------------------------------- #
# Elemento P1 de tetraedro
# --------------------------------------------------------------------------- #


def tet_gradients(nodes, tets) -> tuple[np.ndarray, np.ndarray]:
    """Gradientes (constantes) de las cuatro funciones de forma P1 y volumen [m³].

    Con ``J = [v₁−v₀, v₂−v₀, v₃−v₀]`` (columnas), las coordenadas baricéntricas son
    ``λ_{1..3} = J⁻¹(x − v₀)``, de modo que ``∇λ_{1..3}`` son las **filas** de ``J⁻¹`` y
    ``∇λ₀ = −(∇λ₁+∇λ₂+∇λ₃)``.

    Returns
    -------
    (grads, volumes) : ``grads`` ``(M, 4, 3)`` con ``grads[m, i] = ∇φ_i`` [1/m] y ``volumes``
        ``(M,)`` en [m³] (positivos).
    """
    x = np.asarray(nodes, dtype=float)
    t = np.asarray(tets, dtype=np.int64)
    p = x[t]
    jmat = np.stack([p[:, 1] - p[:, 0], p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]], axis=-1)
    det = np.linalg.det(jmat)
    volumes = det / 6.0
    if np.any(volumes <= 0.0):
        i = int(np.flatnonzero(volumes <= 0.0)[0])
        raise ValueError(
            f"el tetraedro {i} tiene volumen {volumes[i]:.3e} m³ (≤ 0): malla degenerada o mal "
            "orientada; usar meshing3d.check_tet_mesh para localizarlo"
        )
    inv = np.linalg.inv(jmat)
    g1, g2, g3 = inv[:, 0, :], inv[:, 1, :], inv[:, 2, :]
    g0 = -(g1 + g2 + g3)
    return np.stack([g0, g1, g2, g3], axis=1), volumes


def p1_tet_local(nodes, tets) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Matrices locales P1 de tetraedro: ``(M_e, K_e, V)`` con formas ``(m,4,4)`` y ``(m,)``.

    ``M_e = V/20·(1 + I)`` (diagonal ``V/10``, fuera de la diagonal ``V/20``) y
    ``K_e[i,j] = V·∇φ_i·∇φ_j``. ``M_e`` es consistente: su suma es el volumen y
    ``Σ_j M_e[i,j] = V/4``.
    """
    grads, volumes = tet_gradients(nodes, tets)
    m = volumes.size
    mass = np.broadcast_to(volumes[:, None, None] / 20.0, (m, 4, 4)).copy()
    idx = np.arange(4)
    mass[:, idx, idx] += (volumes / 20.0)[:, None]
    stiff = volumes[:, None, None] * np.einsum("mik,mjk->mij", grads, grads)
    return mass, stiff, volumes


def _accumulate(n_dofs: int, dofs: np.ndarray, local: np.ndarray) -> sp.csr_matrix:
    """Suma matrices locales ``(m,n,n)`` en una matriz dispersa global ``(N,N)``."""
    n = dofs.shape[1]
    rows = np.repeat(dofs, n, axis=1).ravel()
    cols = np.tile(dofs, (1, n)).ravel()
    return sp.coo_matrix(
        (np.asarray(local).reshape(-1), (rows, cols)), shape=(n_dofs, n_dofs), dtype=complex
    ).tocsr()


def p1_tet_matrices(nodes, tets) -> tuple[sp.csr_matrix, sp.csr_matrix]:
    """Matrices globales ``(M, K)`` reales de masa y rigidez P1 en tetraedros."""
    x = np.asarray(nodes, dtype=float)
    t = np.asarray(tets, dtype=np.int64)
    mass, stiff, _ = p1_tet_local(x, t)
    n = x.shape[0]
    return _accumulate(n, t, mass), _accumulate(n, t, stiff)


# --------------------------------------------------------------------------- #
# Bloques de frontera
# --------------------------------------------------------------------------- #


def face_mass(nodes, faces) -> np.ndarray:
    """Matriz de masa P1 sobre caras triangulares: ``∫_f u·v dS``, forma ``(F, 3, 3)``.

    ``A/12·(1 + I)``: la misma matriz que la masa de un triángulo 2D. Su suma es el área de la cara.
    """
    x = np.asarray(nodes, dtype=float)
    f = np.asarray(faces, dtype=np.int64)
    areas = face_areas(x, f)
    m = areas.size
    local = np.broadcast_to(areas[:, None, None] / 12.0, (m, 3, 3)).copy()
    idx = np.arange(3)
    local[:, idx, idx] += (areas / 12.0)[:, None]
    return local


def _surface_gradients(nodes, faces) -> np.ndarray:
    """Gradiente tangencial de las tres funciones de forma P1 sobre una cara plana.

    ``∇_s φ_i = n × (P_{i+2} − P_{i+1}) / (2A)`` con ``n`` la normal unitaria de la cara y los
    índices cíclicos (es el gradiente 2D rotado al plano de la cara). Cumple
    ``Σ_i ∇_s φ_i = 0`` y ``∇_s φ_i · (P_j − P_i) = δ_ij``.

    Returns
    -------
    ``(F, 3, 3)`` con ``[:, i, :] = ∇_s φ_i`` [1/m].
    """
    x = np.asarray(nodes, dtype=float)
    f = np.asarray(faces, dtype=np.int64)
    a, b, c = x[f[:, 0]], x[f[:, 1]], x[f[:, 2]]
    cross = np.cross(b - a, c - a)
    two_area = np.linalg.norm(cross, axis=1)
    if np.any(two_area <= 0.0):
        i = int(np.flatnonzero(two_area <= 0.0)[0])
        raise ValueError(f"la cara {i} es degenerada (área nula)")
    n_hat = cross / two_area[:, None]
    g0 = np.cross(n_hat, c - b) / two_area[:, None]
    g1 = np.cross(n_hat, a - c) / two_area[:, None]
    g2 = np.cross(n_hat, b - a) / two_area[:, None]
    return np.stack([g0, g1, g2], axis=1)


def surface_stiffness(nodes, faces) -> np.ndarray:
    """Rigidez tangencial ``∫_f ∇_s u·∇_s v dS`` por cara, forma ``(F, 3, 3)`` [m⁰].

    Es el bloque que acompaña al término ``d·∇²_s p`` de la ABC de orden 2 tras integrar por partes
    sobre la superficie.
    """
    grads = _surface_gradients(nodes, faces)
    areas = face_areas(nodes, faces)
    return areas[:, None, None] * np.einsum("fik,fjk->fij", grads, grads)


def abc_coefficients(k, radius=None, order: int = 1) -> tuple[complex, complex]:
    """Coeficientes ``(c, d)`` de la condición ``∂p/∂n = c·p + d·∇²_s p``.

    ``order=1``: ``c = i·k``, ``d = 0`` (impedancia plana; exacta en incidencia normal).
    ``order=2``: ``c = i·k − 1/r − i/(2·k·r²)``, ``d = i/(2·k)`` (BGT-2 sobre un frente esférico
    de radio ``r``; ver el docstring del módulo).

    Parameters
    ----------
    k : número de onda complejo [1/m].
    radius : radio ``r`` del frente esférico [m]; obligatorio para ``order=2`` y debe ser > 0.
    order : ``1`` o ``2``.
    """
    if order not in (1, 2):
        raise ValueError(f"abc_order debe ser 1 o 2; recibido {order!r}")
    k_arr = np.asarray(k, dtype=complex)
    if not (np.all(np.isfinite(k_arr.real)) and np.all(np.isfinite(k_arr.imag))) \
            or np.any(k_arr.real <= 0.0):
        raise ValueError(f"k debe ser finito con parte real > 0; recibido {k!r}")
    if order == 1:
        c = 1j * k_arr
        d = np.zeros_like(k_arr)
        shape = k_arr.shape
    else:
        if radius is None:
            raise ValueError("abc_order=2 requiere radius (el radio del frente esférico) [m]")
        r = np.asarray(radius, dtype=float)
        if np.any(~np.isfinite(r)) or np.any(r <= 0.0):
            raise ValueError(f"radius debe ser finito y > 0 [m]; recibido {radius!r}")
        shape = np.broadcast_shapes(k_arr.shape, r.shape)
        k_arr = np.broadcast_to(k_arr, shape)
        r = np.broadcast_to(r, shape)
        c = 1j * k_arr - 1.0 / r - 0.5j / (k_arr * r**2)
        d = 0.5j / k_arr
    c = np.broadcast_to(c, shape)
    d = np.broadcast_to(d, shape)
    if shape == ():
        return complex(c), complex(d)
    return np.array(c), np.array(d)



# --------------------------------------------------------------------------- #
# Utilidades de malla para el ensamblado
# --------------------------------------------------------------------------- #


def _encode_faces(faces: np.ndarray, n_nodes: int) -> np.ndarray:
    """Clave entera única de cada cara (tripleta ordenada), para búsquedas vectorizadas."""
    f = np.sort(np.asarray(faces, dtype=np.int64), axis=1)
    n = int(n_nodes) + 1
    return f[:, 0] + n * (f[:, 1] + n * f[:, 2])


def _face_owner(tets: np.ndarray, faces: np.ndarray) -> np.ndarray:
    """Índice del tetraedro dueño de cada cara (cualquiera de los dos si es interior)."""
    t = np.asarray(tets, dtype=np.int64)
    f = np.asarray(faces, dtype=np.int64)
    n_nodes = int(t.max()) + 1
    all_faces = t[:, _LOCAL_FACES].reshape(-1, 3)
    keys = _encode_faces(all_faces, n_nodes)
    order = np.argsort(keys, kind="stable")
    sorted_keys = keys[order]
    query = _encode_faces(f, n_nodes)
    pos = np.searchsorted(sorted_keys, query)
    pos = np.clip(pos, 0, sorted_keys.size - 1)
    found = sorted_keys[pos] == query
    if not np.all(found):
        i = int(np.flatnonzero(~found)[0])
        raise ValueError(
            f"la cara {i} = {np.sort(f[i]).tolist()} no pertenece a ningún tetraedro de la malla"
        )
    return (order[pos] // 4).astype(np.int64)


def _face_values(tets, element_values: np.ndarray, faces: np.ndarray, name: str) -> np.ndarray:
    """Valor (``k`` o ``ρ``) del elemento dueño de cada cara."""
    values = np.asarray(element_values)
    if values.size == 1:
        return np.full(faces.shape[0], values[0], dtype=values.dtype)
    if np.all(values == values[0]):
        return np.full(faces.shape[0], values[0], dtype=values.dtype)
    return values[_face_owner(tets, faces)]


def _as_faces(faces, tets: np.ndarray, name: str) -> np.ndarray | None:
    if faces is None:
        return None
    f = np.asarray(faces, dtype=np.int64)
    if f.ndim == 1 and f.size == 0:
        return f.reshape(0, 3)
    if f.ndim != 2 or f.shape[1] != 3:
        raise ValueError(f"{name} debe tener forma (F, 3); recibida {f.shape}")
    if f.size and (f.min() < 0 or f.max() > int(np.max(tets))):
        raise ValueError(f"{name} tiene índices de nodo fuera de rango")
    return f


def _validate_abc_surface_closed(faces: np.ndarray, name: str) -> None:
    """Exige que las caras formen una superficie cerrada (cada arista en dos caras)."""
    edges = boundary_edges(faces)  # únicas
    f = np.sort(np.asarray(faces, dtype=np.int64), axis=1)
    e = np.vstack([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]])
    e = np.sort(e, axis=1)
    n_nodes = int(e.max()) + 1
    key = e[:, 0] + n_nodes * e[:, 1]
    _, counts = np.unique(key, return_counts=True)
    if edges.shape[0] and np.any(counts != 2):
        raise ValueError(
            f"{name}: abc_order=2 necesita que las caras de la ABC formen una superficie cerrada "
            "(cada arista pertenece a exactamente dos caras); hay "
            f"{int(np.count_nonzero(counts != 2))} aristas que no cumplen"
        )


def _add_point_sources(b, nodes, point_sources, volume_velocities, frequency) -> None:
    """Suma ``−i·ω·Q`` en el nodo más cercano a cada fuente (delta lumpado)."""
    x = np.asarray(nodes, dtype=float)
    pos = np.asarray(point_sources, dtype=float)
    if pos.ndim == 1:
        pos = pos.reshape(1, -1)
    if pos.ndim != 2 or pos.shape[1] != 3:
        raise ValueError(f"point_sources debe tener forma (P, 3); recibida {pos.shape}")
    if not np.all(np.isfinite(pos)):
        raise ValueError("point_sources contiene valores no finitos")
    q = np.asarray(volume_velocities, dtype=complex).reshape(-1)
    if q.size == 1 and pos.shape[0] > 1:
        q = np.full(pos.shape[0], q[0], dtype=complex)
    if q.size != pos.shape[0]:
        raise ValueError(f"point_volume_velocities debe tener un valor por fuente ({pos.shape[0]}); "
                         f"recibidos {q.size}")
    if not np.all(np.isfinite(q.real)) or not np.all(np.isfinite(q.imag)):
        raise ValueError("point_volume_velocities contiene valores no finitos")
    if frequency is None:
        raise ValueError("point_sources requiere frequency [Hz] para armar el término −i·ω·Q")
    omega = TWO_PI * _positive_frequency(frequency)

    for p, qi in zip(pos, q):
        d = np.linalg.norm(x - p[None, :], axis=1)
        k = int(np.argmin(d))
        if d[k] > 0.0:
            # el delta lumpado sólo es defendible si la fuente cae sobre un nodo
            raise ValueError(
                f"la fuente puntual {p.tolist()} no coincide con ningún nodo (el más cercano está "
                f"a {d[k]:.3e} m). El delta se lumpa en un nodo de la malla: ubicar la fuente "
                "sobre un nodo o refinar"
            )
        b[k] += -1j * omega * qi


# --------------------------------------------------------------------------- #
# Ensamblado de campo total
# --------------------------------------------------------------------------- #


def assemble(
    nodes,
    tets,
    element_k,
    abc_faces=None,
    p_inc: PlaneWave3D | None = None,
    element_rho=None,
    face_k=None,
    face_rho=None,
    normals=None,
    illuminated=None,
    abc_order: int = 1,
    abc_center=None,
    point_sources=None,
    point_volume_velocities=None,
    frequency=None,
) -> tuple[sp.csr_matrix, np.ndarray]:
    """Ensambla ``A·p = b`` para el campo **total** con ABC sobre caras planas.

    Parameters
    ----------
    nodes : ``(N, 3)`` coordenadas [m].
    tets : ``(M, 4)`` tetraedros con orientación positiva (`meshing3d`).
    element_k : ``(M,)`` (o escalar) número de onda complejo [1/m] por elemento.
    abc_faces : ``(F, 3)`` caras donde se aplica la ABC; ``None`` ⇒ sin ABC (frontera rígida).
    p_inc : onda plana incidente; ``None`` ⇒ sin término incidente (campo saliente o fuente).
    element_rho : ``(M,)`` o escalar densidad [kg/m³]; por defecto ``1.0`` (deja la ecuación en la
        forma ``∇p·∇v − k²pv``). Con densidad variable hay que pasarla: la ecuación fuerte correcta
        es ``∇·((1/ρ)∇p) + (k²/ρ)p = 0``.
    face_k : ``(F,)`` número de onda por cara de la ABC; ``None`` ⇒ el del elemento dueño.
    face_rho : ``(F,)`` densidad por cara de la ABC; ``None`` ⇒ la del elemento dueño.
    normals : ``(F, 3)`` normales de las caras de la ABC; si se pasan se validan contra la malla
        (una normal invertida no cambia la ABC —ésta no la usa— pero sí el término incidente).
    illuminated : ``(F,)`` bool; ``True`` donde la onda incidente **entra** al dominio
        (``n·d < 0``). ``None`` ⇒ se deduce de ``p_inc.direction``. Es imprescindible que las caras
        de sombra queden en ``False``: ahí el campo total es saliente y el término incidente, con
        ``k`` local distinto del de la onda, actúa como fuente espuria (D11.2).
    abc_order : ``1`` o ``2`` (ver :func:`abc_coefficients`).
    abc_center : centro del frente esférico [m], ``(3,)``; obligatorio con ``abc_order=2``.
    point_sources : ``(P, 3)`` posiciones de fuentes monopolares [m]; deben coincidir con nodos.
    point_volume_velocities : ``(P,)`` caudales volumétricos complejos [m³/s].
    frequency : [Hz]; obligatorio si hay fuentes puntuales (para ``ω``).

    Returns
    -------
    (A, b) : matriz dispersa compleja ``(N, N)`` y vector ``(N,)``.
    """
    x = np.asarray(nodes, dtype=float)
    t = np.asarray(tets, dtype=np.int64)
    check_tet_mesh(x, t)
    n = x.shape[0]

    k_el = _element_k_array(element_k, t.shape[0])
    rho_el = _element_rho_array(element_rho, t.shape[0])

    mass, stiff, _ = p1_tet_local(x, t)
    local = (stiff - (k_el**2)[:, None, None] * mass) / rho_el[:, None, None]
    a = _accumulate(n, t, local)
    b = np.zeros(n, dtype=complex)

    faces = _as_faces(abc_faces, t, "abc_faces")
    if faces is not None and faces.shape[0] > 0:
        a = a + _abc_block(x, t, faces, k_el, rho_el, face_k, face_rho, normals, abc_order,
                           abc_center)
        if p_inc is not None:
            b += _incident_rhs(x, t, faces, k_el, rho_el, face_k, face_rho, normals,
                               illuminated, p_inc)
    elif p_inc is not None:
        raise ValueError(
            "p_inc requiere abc_faces: el término incidente se integra sobre las caras de la ABC "
            "(sin ellas no hay dónde imponer la condición de radiación)"
        )

    if point_sources is not None:
        _add_point_sources(b, x, point_sources, point_volume_velocities, frequency)
    elif point_volume_velocities is not None:
        raise ValueError("point_volume_velocities se pasó sin point_sources")

    return a.tocsr(), b


def _abc_block(x, t, faces, k_el, rho_el, face_k, face_rho, normals, abc_order, abc_center):
    """Bloque de la ABC sobre las caras dadas (masa de borde y, en orden 2, rigidez tangencial)."""
    n = x.shape[0]
    k_f = (np.asarray(face_k, dtype=complex).reshape(-1) if face_k is not None
           else _face_values(t, k_el, faces, "k"))
    if k_f.size != faces.shape[0]:
        raise ValueError(f"face_k debe tener un valor por cara ({faces.shape[0]}); "
                         f"recibido {k_f.size}")
    if not (np.all(np.isfinite(k_f.real)) and np.all(np.isfinite(k_f.imag))):
        raise ValueError("face_k contiene valores no finitos")
    if np.any(k_f.real <= 0.0):
        raise ValueError("face_k debe tener parte real > 0")

    if face_rho is not None:
        rho_f = np.asarray(face_rho, dtype=float).reshape(-1)
        if rho_f.size != faces.shape[0]:
            raise ValueError(f"face_rho debe tener un valor por cara ({faces.shape[0]}); "
                             f"recibido {rho_f.size}")
        if np.any(~np.isfinite(rho_f)) or np.any(rho_f <= 0.0):
            raise ValueError("face_rho debe ser finito y > 0 [kg/m³]")
    else:
        rho_f = _face_values(t, rho_el, faces, "rho").astype(float)

    if normals is not None:
        nrm = np.asarray(normals, dtype=float).reshape(-1, 3)
        if nrm.shape[0] != faces.shape[0]:
            raise ValueError(f"normals debe tener una normal por cara ({faces.shape[0]}); "
                             f"recibidas {nrm.shape[0]}")
        if np.any(np.einsum("ij,ij->i", nrm, face_normals(x, faces)) <= 0.0):
            raise ValueError(
                "normals tiene caras con la normal invertida respecto de la malla "
                "(n·n_malla ≤ 0); las normales de la ABC son salientes del dominio"
            )

    if abc_order == 2:
        if abc_center is None:
            raise ValueError("abc_order=2 requiere abc_center (centro del frente esférico) [m]")
        center = np.asarray(abc_center, dtype=float).reshape(-1)
        if center.size != 3 or not np.all(np.isfinite(center)):
            raise ValueError(f"abc_center debe ser un punto 3D finito; recibido {abc_center!r}")
        _validate_abc_surface_closed(faces, "abc_faces")
        mid = face_centroids(x, faces)
        r_f = np.linalg.norm(mid - center[None, :], axis=1)
        if np.any(r_f <= 0.0):
            raise ValueError("abc_center cae sobre alguna cara de la ABC (radio nulo)")
        c_f, d_f = abc_coefficients(k_f, r_f, order=2)
    else:
        if abc_center is not None:
            raise ValueError("abc_center sólo tiene sentido con abc_order=2")
        c_f, d_f = abc_coefficients(k_f, None, order=1)

    w = 1.0 / rho_f
    local = (-c_f * w)[:, None, None] * face_mass(x, faces)
    a = _accumulate(n, faces, local)
    if abc_order == 2:
        a = a + _accumulate(n, faces, (d_f * w)[:, None, None] * surface_stiffness(x, faces))
    return a


def _incident_rhs(x, t, faces, k_el, rho_el, face_k, face_rho, normals, illuminated, p_inc) -> np.ndarray:
    """Lado derecho ``∫ (1/ρ)(∂p_i/∂n − i·k·p_i)·v dS`` en las caras iluminadas."""
    n = x.shape[0]
    b = np.zeros(n, dtype=complex)
    k_f = (np.asarray(face_k, dtype=complex).reshape(-1) if face_k is not None
           else _face_values(t, k_el, faces, "k"))
    rho_f = (np.asarray(face_rho, dtype=float).reshape(-1) if face_rho is not None
             else _face_values(t, rho_el, faces, "rho").astype(float))
    nrm = (np.asarray(normals, dtype=float).reshape(-1, 3) if normals is not None
           else face_normals(x, faces))

    if illuminated is None:
        direction = getattr(p_inc, "direction", None)
        illum = (np.ones(faces.shape[0], dtype=bool) if direction is None
                 else (np.asarray(nrm) @ np.asarray(direction, dtype=float)) < 0.0)
    else:
        illum = np.asarray(illuminated, dtype=bool).reshape(-1)
        if illum.size != faces.shape[0]:
            raise ValueError(f"illuminated debe tener un valor por cara ({faces.shape[0]}); "
                             f"recibidos {illum.size}")

    mid = face_centroids(x, faces)
    areas = face_areas(x, faces)
    f = p_inc.normal_derivative(mid, nrm) - 1j * k_f * p_inc.value(mid)
    contrib = np.where(illum, f * areas / 3.0 / rho_f, 0.0).astype(complex)
    for j in range(3):
        np.add.at(b, faces[:, j], contrib)
    return b


def solve(
    nodes,
    tets,
    element_k,
    abc_faces=None,
    p_inc: PlaneWave3D | None = None,
    element_rho=None,
    face_k=None,
    face_rho=None,
    normals=None,
    illuminated=None,
    abc_order: int = 1,
    abc_center=None,
    point_sources=None,
    point_volume_velocities=None,
    frequency=None,
) -> np.ndarray:
    """Resuelve el Helmholtz 3D P1 y devuelve la presión total compleja ``(N,)`` [Pa].

    Factoriza con ``scipy.sparse.linalg.splu`` (D10: ~17× más rápido que ``spsolve`` en 2D; en 3D
    la diferencia es mayor). Garantiza que la salida es finita: si el sistema sale con ``NaN`` o
    ``inf`` levanta ``RuntimeError`` en vez de devolverlos.

    Ver :func:`assemble` para el significado de los argumentos.
    """
    a, b = assemble(
        nodes, tets, element_k, abc_faces=abc_faces, p_inc=p_inc, element_rho=element_rho,
        face_k=face_k, face_rho=face_rho, normals=normals, illuminated=illuminated,
        abc_order=abc_order, abc_center=abc_center, point_sources=point_sources,
        point_volume_velocities=point_volume_velocities, frequency=frequency,
    )
    return _solve_sparse(a, b)


def _solve_sparse(a: sp.csr_matrix, b: np.ndarray) -> np.ndarray:
    """Factoriza y resuelve con ``splu``, con control de finitud y de residuo.

    Orden de intentos (medido en esta máquina, caja de 6 253 gdl: **32 s → 1.8 s**):

    1. ``permc_spec="MMD_AT_PLUS_A"`` + ``SymmetricMode`` + ``diag_pivot_thresh=0``. La matriz del
       FEM P1 es **estructuralmente simétrica** (rigidez, masa y bloque de borde), y con el
       reordenamiento simétrico SuperLU no paga el pivoteo dinámico: 17× menos tiempo para el mismo
       residuo (1e-14). Es el camino rápido.
    2. Si el residuo relativo de (1) supera ``1e-8`` —posible con un pivoteo que no conviene
       omitir— se rehace con el pivoteo por defecto de SuperLU, que es más lento pero robusto.

    La salida se valida siempre: primero finitud y después residuo relativo ``‖A·p − b‖/‖b‖``.
    """
    a = a.tocsr()
    b = np.asarray(b, dtype=complex).reshape(-1)
    if b.size != a.shape[0]:
        raise ValueError(f"b tiene {b.size} componentes y A es {a.shape[0]}×{a.shape[1]}")

    fast = dict(
        permc_spec="MMD_AT_PLUS_A",
        diag_pivot_thresh=0.0,
        options=dict(SymmetricMode=True),
    )
    robust = dict(permc_spec="COLAMD", diag_pivot_thresh=1.0)

    factorizations: list[tuple[str, dict]] = [("rápida (simétrica)", fast),
                                              ("robusta (pivoteo completo)", robust)]
    norm_b = float(np.linalg.norm(b))
    if norm_b == 0.0:
        norm_b = 1.0

    last = ""
    for name, kwargs in factorizations:
        try:
            lu = spla.splu(a.tocsc(), **kwargs)
            x = np.asarray(lu.solve(b), dtype=complex).reshape(-1)
        except RuntimeError as exc:
            last = f"{name}: {exc}"
            continue
        if not np.all(np.isfinite(x.real)) or not np.all(np.isfinite(x.imag)):
            last = f"{name}: la solución tiene valores no finitos"
            continue
        residual = float(np.linalg.norm(a @ x - b)) / norm_b
        if residual > 1e-8:
            last = f"{name}: residuo relativo {residual:.2e} > 1e-8"
            continue
        return x

    raise RuntimeError(
        f"no se pudo resolver el sistema de {a.shape[0]} gdl ({a.nnz} no nulos). "
        f"Último intento — {last}. El sistema está mal escalado o es singular: revisar que el "
        "dominio tenga condiciones de frontera en todas las caras (ABC o Neumann) y que k no sea nulo"
    )



# --------------------------------------------------------------------------- #
# Campo dispersado p_s = p − p_i
# --------------------------------------------------------------------------- #


def assemble_scattered(
    nodes,
    tets,
    element_k,
    abc_faces=None,
    neumann_faces=None,
    dp_dn_incident=None,
    element_rho=None,
    face_rho=None,
    normals=None,
    abc_order: int = 1,
    abc_center=None,
    neumann_normals=None,
    point_sources=None,
    point_volume_velocities=None,
    frequency=None,
) -> tuple[sp.csr_matrix, np.ndarray]:
    r"""Ensambla ``A·p_s = b`` para el campo **dispersado** ``p_s = p − p_i``.

    Forma débil::

        ∫ (1/ρ)·(∇p_s·∇v − k²·p_s·v) dΩ + ∫_ABC (1/ρ)·[−c·p_s·v + d·∇_s p_s·∇_s v] dS
            = − ∫_Γint (1/ρ)·(∂p_i/∂n)·v dS

    El contorno ABC es **homogéneo** (``p_s`` es saliente en todo él: no hay máscara de
    iluminación). La onda incidente entra sólo como dato de Neumann sobre la superficie del
    dispersor, ``∂p_s/∂n = −∂p_i/∂n`` (pared rígida: ``∂p_total/∂n = 0``).

    Parameters
    ----------
    abc_faces, abc_order, abc_center, face_rho, normals : como en :func:`assemble`.
    neumann_faces : ``(E, 3)`` caras del dispersor donde se impone ``∂p_s/∂n = −∂p_i/∂n``.
    dp_dn_incident : ``(E,)`` ``∂p_i/∂n`` [Pa/m] en el centroide de cada cara de
        ``neumann_faces``, con la normal **saliente del dominio FEM** (en el hueco de un anillo,
        hacia el eje). Se reorienta y valida contra la malla.
    neumann_normals : ``(E, 3)`` normales de esas caras; se auditan contra la malla.

    Returns
    -------
    (A, b) : ``p_s`` en los nodos.
    """
    x = np.asarray(nodes, dtype=float)
    t = np.asarray(tets, dtype=np.int64)
    check_tet_mesh(x, t)
    n = x.shape[0]

    k_el = _element_k_array(element_k, t.shape[0])
    rho_el = _element_rho_array(element_rho, t.shape[0])

    mass, stiff, _ = p1_tet_local(x, t)
    local = (stiff - (k_el**2)[:, None, None] * mass) / rho_el[:, None, None]
    a = _accumulate(n, t, local)
    b = np.zeros(n, dtype=complex)

    faces = _as_faces(abc_faces, t, "abc_faces")
    if faces is not None and faces.shape[0] > 0:
        a = a + _abc_block(x, t, faces, k_el, rho_el, None, face_rho, normals, abc_order,
                           abc_center)

    neu = _as_faces(neumann_faces, t, "neumann_faces")
    if neu is not None:
        if dp_dn_incident is None:
            raise ValueError(
                "neumann_faces requiere dp_dn_incident = ∂p_i/∂n sobre los centroides "
                "(normal saliente del dominio)"
            )
        dp = np.asarray(dp_dn_incident, dtype=complex).reshape(-1)
        if dp.size != neu.shape[0]:
            raise ValueError(f"dp_dn_incident debe tener un valor por cara de neumann_faces "
                             f"({neu.shape[0]}); recibidos {dp.size}")
        if not (np.all(np.isfinite(dp.real)) and np.all(np.isfinite(dp.imag))):
            raise ValueError("dp_dn_incident contiene valores no finitos")
        if neumann_normals is not None:
            nrm_n = np.asarray(neumann_normals, dtype=float).reshape(-1, 3)
            if nrm_n.shape[0] != neu.shape[0]:
                raise ValueError(f"neumann_normals debe tener una normal por cara "
                                 f"({neu.shape[0]}); recibidas {nrm_n.shape[0]}")
            if np.any(np.einsum("ij,ij->i", nrm_n, face_normals(x, neu)) <= 0.0):
                raise ValueError(
                    "neumann_normals tiene caras con la normal invertida respecto de la malla; "
                    "∂p_i/∂n se define con la normal saliente del dominio FEM"
                )
        if faces is not None and faces.shape[0]:
            key_abc = set(map(tuple, np.sort(faces, axis=1).tolist()))
            key_neu = set(map(tuple, np.sort(neu, axis=1).tolist()))
            overlap = key_abc & key_neu
            if overlap:
                raise ValueError(
                    "abc_faces y neumann_faces se solapan; una cara no puede ser ABC y dato de "
                    f"Neumann a la vez (p. ej. {sorted(overlap)[0]})"
                )
        rho_n = (np.asarray(face_rho, dtype=float).reshape(-1) if face_rho is not None
                 else _face_values(t, rho_el, neu, "rho").astype(float))
        areas = face_areas(x, neu)
        # b = −∮ (1/ρ)(∂p_i/∂n)·v dS, regla del punto medio
        contrib = -(dp * areas / 3.0) / rho_n
        for j in range(3):
            np.add.at(b, neu[:, j], contrib)
    elif dp_dn_incident is not None:
        raise ValueError("dp_dn_incident se pasó sin neumann_faces")

    if point_sources is not None:
        _add_point_sources(b, x, point_sources, point_volume_velocities, frequency)
    elif point_volume_velocities is not None:
        raise ValueError("point_volume_velocities se pasó sin point_sources")

    return a.tocsr(), b


def solve_scattered(
    nodes,
    tets,
    element_k,
    abc_faces=None,
    neumann_faces=None,
    dp_dn_incident=None,
    element_rho=None,
    face_rho=None,
    normals=None,
    abc_order: int = 1,
    abc_center=None,
    neumann_normals=None,
    point_sources=None,
    point_volume_velocities=None,
    frequency=None,
) -> np.ndarray:
    """Resuelve el campo **dispersado** ``p_s`` en los nodos (ver :func:`assemble_scattered`)."""
    a, b = assemble_scattered(
        nodes, tets, element_k, abc_faces=abc_faces, neumann_faces=neumann_faces,
        dp_dn_incident=dp_dn_incident, element_rho=element_rho, face_rho=face_rho,
        normals=normals, abc_order=abc_order, abc_center=abc_center,
        neumann_normals=neumann_normals, point_sources=point_sources,
        point_volume_velocities=point_volume_velocities, frequency=frequency,
    )
    return _solve_sparse(a, b)


def total_field(p_scattered, p_inc, nodes=None) -> np.ndarray:
    """``p_total = p_s + p_i`` en los nodos.

    Parameters
    ----------
    p_scattered : ``(N,)`` salida de :func:`solve_scattered`.
    p_inc : :class:`PlaneWave3D`, :class:`PointSource` o callable ``f(points) -> (M,)``.
    nodes : ``(N, 3)`` coordenadas donde evaluar ``p_i``.
    """
    ps = np.asarray(p_scattered, dtype=complex).reshape(-1)
    if hasattr(p_inc, "value"):
        if nodes is None:
            raise ValueError("nodes es obligatorio cuando p_inc es una onda o fuente")
        pi = np.asarray(p_inc.value(nodes), dtype=complex).reshape(-1)
    elif callable(p_inc):
        if nodes is None:
            raise ValueError("nodes es obligatorio cuando p_inc es un callable")
        pi = np.asarray(p_inc(nodes), dtype=complex).reshape(-1)
    else:
        raise TypeError("p_inc debe tener .value (onda/fuente) o ser un callable")
    if pi.size != ps.size:
        raise ValueError(f"p_inc devolvió {pi.size} valores y p_scattered tiene {ps.size}; "
                         "los nodos no coinciden")
    return ps + pi
