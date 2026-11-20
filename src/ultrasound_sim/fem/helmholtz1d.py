"""FEM P1 de Helmholtz 1D (Fase 0) con radiación saliente en ambos extremos.

Forma débil (convención D4: ``e^{−iωt}``)::

    ∫ (1/ρ)·(p' v' − k² p v) dx = [ (1/ρ)·p' n v ]_∂Ω

El peso ``1/ρ`` es el de la ecuación heterogénea ``∇·((1/ρ)·∇p) + (k²/ρ)·p = 0``, que por tramos
homogéneos se lee ``(1/ρ)·p'' + (k²/ρ)·p = 0``. Su condición natural en la interfaz es la
continuidad de ``(1/ρ)·p'``, es decir de la velocidad normal (``u = (1/(i·ω·ρ))·p'``); sin el peso,
un contraste de densidad da un R/T incorrecto. Con ``ρ`` uniforme se reduce a
``∫(p' v' − k²p v)`` y el resultado no cambia.

En 1D la condición de radiación de primer orden es **exacta** (no hay ángulos oblicuos), de modo
que el único error que queda es el de discretización P1.

Extremos
--------
Con ``p = p_i + p_r`` en ``x = x₀`` y onda incidente ``p_i = p_inc·exp(+i·k₁·x)``:

- extremo izquierdo, normal saliente ``n = −1``:  ``∂p/∂n = i·k₁·(p − 2·p_i)``
- extremo derecho,     normal saliente ``n = +1``:  ``∂p/∂n = i·k₂·p``

Con ``q_j = i·k_j/ρ_j`` (el mismo ``q`` de `analytic.helmholtz1d`), el sistema queda con la
diagonal compleja ``−q`` en los dos nodos extremos y un término conocido a la derecha::

    (K̃ − Σ_e (k_e²/ρ_e)·M_e − q₁·e₀e₀ᵀ − q₂·e_N e_Nᵀ) p = −2·q₁·p_i(x₀)·e₀

La malla sale de `meshing.layered_line_mesh` con nodos forzados en las interfaces, así que cada
elemento pertenece a un único material (no hay elementos que promedien dos medios).

Extracción de coeficientes
--------------------------
En la región izquierda (``x < 0``) el campo es una superposición de onda incidente y reflejada;
un ajuste por mínimos cuadrados de ``a·exp(i·k₁·x) + b·exp(−i·k₁·x)`` sobre los nodos de esa
región da ``R = b/a`` (referido a ``x = 0``). En la región derecha (``x > L``) el mismo
procedimiento con ``a·exp(i·k₂·x) + b·exp(−i·k₂·x)`` da ``T = a/p_inc``, con la misma definición
de ``T`` que `analytic.helmholtz1d` (referida a ``x = 0``).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..analytic.helmholtz1d import _validate
from ..materials import Material
from ..meshing import LineMesh, layered_line_mesh
from .assembly import p1_line_local

__all__ = ["WaveCoefficients", "solve", "extract_coefficients", "build_mesh"]


@dataclass(frozen=True)
class WaveCoefficients:
    """Coeficientes extraídos de un campo 1D por descomposición de onda.

    ``R`` y ``T`` están referidos a ``x = 0`` (igual que `analytic.helmholtz1d`); ``a_left``
    y ``b_right`` son las amplitudes crudas del ajuste, útiles para diagnosticar.
    """

    R: complex
    T: complex
    a_left: complex
    b_left: complex
    a_right: complex
    b_right: complex


def build_mesh(
    frequency: float,
    media,
    thicknesses,
    points_per_wavelength: float = 10.0,
    pad_wavelengths: float = 2.0,
) -> tuple[LineMesh, np.ndarray]:
    """Malla 1D del stack con regiones de amortiguación (padding) a ambos lados.

    Se agrega ``pad_wavelengths`` longitudes de onda a cada lado; esas regiones hacen falta para
    poder descomponer la onda (ajustar incidente/reflejada) dentro del dominio.

    El paso se toma como el mínimo entre los ``λ_i/points_per_wavelength`` de todos los tramos,
    porque `meshing.layered_line_mesh` acepta un único ``max_step``. Es conservador: los tramos
    con λ grande quedan sobre-refinados.
    """
    media_t, th = _validate(media, thicknesses)
    if th.size == 0:
        raise ValueError("se requiere al menos una capa (thicknesses no vacío)")
    if not points_per_wavelength > 0.0:
        raise ValueError("points_per_wavelength debe ser > 0")

    k = np.array([m.wavenumber(frequency) for m in media_t], dtype=complex)
    lambda_ = 2.0 * np.pi / np.real(k)
    length = float(th.sum())

    pad_left = pad_wavelengths * lambda_[0]
    pad_right = pad_wavelengths * lambda_[-1]

    # Longitudes de onda por tramo: [pad izq, capas..., pad der]
    lam_per_segment = np.concatenate([[lambda_[0]], lambda_[1:-1], [lambda_[-1]]])
    max_step = float(np.min(lam_per_segment) / points_per_wavelength)

    interfaces = np.concatenate([[-pad_left], [0.0], np.cumsum(th), [length + pad_right]])
    return layered_line_mesh(interfaces.tolist(), max_step), interfaces


def _element_k(mesh: LineMesh, interfaces: np.ndarray, k: np.ndarray) -> np.ndarray:
    """Número de onda complejo de cada elemento, por pertenencia de su punto medio a un tramo."""
    mid = 0.5 * (mesh.nodes[mesh.elements[:, 0]] + mesh.nodes[mesh.elements[:, 1]])
    segment = np.searchsorted(interfaces, mid, side="right") - 1
    segment = np.clip(segment, 0, k.size - 1)
    return k[segment]


def solve(
    frequency: float,
    media,
    thicknesses,
    points_per_wavelength: float = 10.0,
    p_inc: float = 1.0,
    pad_wavelengths: float = 2.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Resuelve el Helmholtz 1D por FEM P1 con ABC de primer orden en ambos extremos.

    Parameters
    ----------
    frequency : frecuencia [Hz].
    media : secuencia de :class:`Material` de longitud ``len(thicknesses) + 2``.
    thicknesses : espesores de las capas [m]; las capas ocupan ``x = 0 .. sum(thicknesses)``.
    points_per_wavelength : nodos por longitud de onda del tramo más fino.
    p_inc : amplitud de la onda incidente ``p_inc·exp(+i·k₁·x)`` [Pa].
    pad_wavelengths : extensión de las regiones izquierda/derecha, en longitudes de onda.

    Returns
    -------
    (x, p) : coordenadas nodales ``(N,)`` y presión compleja ``(N,)``.
    """
    mesh, interfaces = build_mesh(
        frequency, media, thicknesses, points_per_wavelength, pad_wavelengths
    )
    media_t, th = _validate(media, thicknesses)
    k = np.array([m.wavenumber(frequency) for m in media_t], dtype=complex)
    rho = np.array([m.rho for m in media_t], dtype=float)
    k_el = _element_k(mesh, interfaces, k)
    rho_el = _element_k(mesh, interfaces, rho)

    mass, stiff, _ = p1_line_local(mesh.nodes, mesh.elements)
    n = mesh.nodes.size
    dofs = mesh.elements

    # (1/ρ)·(K − k²·M) por elemento: la forma heterogénea correcta
    local = (stiff - (k_el**2)[:, None, None] * mass) / rho_el[:, None, None]
    a = np.zeros((n, n), dtype=complex)
    rows = np.repeat(dofs, 2, axis=1).ravel()
    cols = np.tile(dofs, (1, 2)).ravel()
    np.add.at(a, (rows, cols), local.reshape(-1))

    k_left = k[0]
    k_right = k[-1]
    q_left = 1j * k_left / rho[0]
    q_right = 1j * k_right / rho[-1]
    x = np.asarray(mesh.nodes, dtype=float)

    # Términos de frontera: −q en los dos nodos extremos.
    a[0, 0] -= q_left
    a[-1, -1] -= q_right

    # Lado derecho: −2·q₁·p_i(x₀) en el nodo izquierdo.
    p_i_left = p_inc * np.exp(1j * k_left * x[0])
    b = np.zeros(n, dtype=complex)
    b[0] = -2.0 * q_left * p_i_left

    p = np.linalg.solve(a, b)
    return x, p


def _fit_two_waves(
    x: np.ndarray, p: np.ndarray, k: complex
) -> tuple[complex, complex]:
    """Mínimos cuadrados de ``a·exp(i·k·x) + b·exp(−i·k·x)`` sobre los nodos dados."""
    if x.size < 3:
        raise ValueError("se necesitan al menos 3 nodos para separar las dos ondas")
    design = np.column_stack([np.exp(1j * k * x), np.exp(-1j * k * x)])
    a, b = np.linalg.lstsq(design, p, rcond=None)[0]
    return complex(a), complex(b)


def extract_coefficients(
    x,
    p,
    frequency: float,
    media,
    thicknesses,
    p_inc: float = 1.0,
) -> WaveCoefficients:
    """Extrae ``R`` y ``T`` del campo FEM por descomposición de onda a cada lado del stack.

    Usa los nodos de ``x < 0`` para ``R`` y los de ``x > L`` para ``T``, con ``L = sum(thicknesses)``.
    """
    x = np.asarray(x, dtype=float)
    p = np.asarray(p, dtype=complex)
    media_t, th = _validate(media, thicknesses)
    k = np.array([m.wavenumber(frequency) for m in media_t], dtype=complex)
    length = float(th.sum())

    left = x < 0.0
    right = x > length
    if np.count_nonzero(left) < 3 or np.count_nonzero(right) < 3:
        raise ValueError(
            "hacen falta al menos 3 nodos en cada región exterior; "
            "aumentá pad_wavelengths o points_per_wavelength"
        )

    a_l, b_l = _fit_two_waves(x[left], p[left], k[0])
    a_r, b_r = _fit_two_waves(x[right], p[right], k[-1])
    return WaveCoefficients(
        R=b_l / a_l,
        T=a_r / p_inc,
        a_left=a_l,
        b_left=b_l,
        a_right=a_r,
        b_right=b_r,
    )
