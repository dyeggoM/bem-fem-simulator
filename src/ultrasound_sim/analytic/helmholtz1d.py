"""Helmholtz 1D exacto en medios estratificados, incidencia normal.

Referencia analítica para los tests de FEM 1D (Fase 0) y para los casos 2D que son
efectivamente 1D (incidencia normal, campo independiente de ``y``).

Convención del repo (`docs/decisions.md` D4):

- tiempo ``exp(-i·ω·t)``, ondas salientes ``∝ exp(+i·k·r)``, ``∇²p + k²p = 0``;
- ``k = ω/c + i·α`` con parte imaginaria positiva: ``exp(+i·k·x)`` decae con ``x``;
- unidades SI; la API recibe Hz.

Geometría
---------
Las capas ocupan ``x = 0 .. L`` con ``L = sum(thicknesses)``; el medio ``media[0]`` es el
semiespacio ``x < 0`` y ``media[-1]`` el semiespacio ``x > L``. Es decir
``len(media) == len(thicknesses) + 2``.

Método
------
Matriz de transferencia 2×2 por capa sobre el vector de estado ``(p, w)`` con
``w = (1/ρ)·dp/dx``. La continuidad de presión y de velocidad normal
(``u_n = (1/(i·ω·ρ))·dp/dx``, con ``i·ω`` común a ambos lados) equivale a la continuidad de
``(p, w)``. Para una capa de espesor ``d`` y número de onda ``k``::

    [ cos(kd)          ρ·sin(kd)/k ]
    [ −k·sin(kd)/ρ     cos(kd)     ]

cuyo determinante es 1 (unimodular, como corresponde a un problema de transmisión).

Coeficientes
------------
``R`` y ``T`` son coeficientes de amplitud de **presión** referidos a la interfaz izquierda
``x = 0``:

    x < 0           p = p_inc·[ exp(i·k₁·x) + R·exp(−i·k₁·x) ]
    0 ≤ x ≤ L       campo dentro de las capas (continuo)
    x > L           p = p_inc·T·exp(i·k₂·x)

``T`` incluye por lo tanto la fase (y la atenuación) acumulada en el espesor del stack. Para
una interfaz única (``L = 0``) esta definición coincide con la fórmula clásica de impedancias
``R = (Z₂−Z₁)/(Z₂+Z₁)``, ``T = 2·Z₂/(Z₁+Z₂)`` cuando no hay atenuación.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..materials import Material

__all__ = ["StackResponse", "stack_response", "pressure", "state_at"]


@dataclass(frozen=True)
class StackResponse:
    """Coeficientes de un stack estratificado a una frecuencia dada.

    Atributos
    ---------
    R : coeficiente de reflexión de presión en ``x = 0`` (amplitud, sin ``p_inc``).
    T : coeficiente de transmisión de presión referido a ``x = 0``
        (``p(x>L) = T·p_inc·exp(i·k_derecha·x)``).
    frequency : frecuencia [Hz].
    total_thickness : espesor total del stack [m].
    k : números de onda complejos de cada medio, en orden (``len(media)``).
    impedance : impedancias ``ρ·c`` reales de cada medio [Pa·s/m].
    """

    R: complex
    T: complex
    frequency: float
    total_thickness: float
    k: tuple[complex, ...]
    impedance: tuple[float, ...]

    # -- comodidades ---------------------------------------------------------

    @property
    def k_left(self) -> complex:
        return self.k[0]

    @property
    def k_right(self) -> complex:
        return self.k[-1]

    @property
    def reflected_power_fraction(self) -> float:
        """Fracción de intensidad reflejada en el semiespacio izquierdo (sin pérdidas)."""
        return float(abs(self.R) ** 2)

    def as_dict(self) -> dict[str, object]:
        return {
            "R": self.R,
            "T": self.T,
            "frequency": self.frequency,
            "total_thickness": self.total_thickness,
        }


def _validate(media, thicknesses) -> tuple[tuple[Material, ...], np.ndarray]:
    media_t = tuple(media)
    if not all(isinstance(m, Material) for m in media_t):
        raise TypeError("media debe ser una secuencia de Material")
    th = np.asarray(thicknesses, dtype=float).reshape(-1)
    if media_t.__len__() != th.size + 2:
        raise ValueError(
            "se requiere len(media) == len(thicknesses) + 2 "
            f"(media={len(media_t)}, thicknesses={th.size})"
        )
    if th.size and np.any(th <= 0.0):
        raise ValueError("los espesores deben ser estrictamente positivos")
    return media_t, th


def layer_matrix(k: complex, rho: float, thickness: float) -> np.ndarray:
    """Matriz de transferencia ``(p, w)`` de una capa homogénea de espesor ``thickness``."""
    kd = k * thickness
    c = np.cos(kd)
    s = np.sin(kd)
    return np.array(
        [[c, rho * s / k], [-k * s / rho, c]],
        dtype=complex,
    )


@dataclass(frozen=True)
class _Internal:
    """Estado interno: números de onda, interfaces y estados ``(p, w)`` en cada interfaz."""

    k: np.ndarray
    x_iface: np.ndarray  # (n_capas + 1,) posiciones de las interfaces, x_iface[0] = 0
    states: np.ndarray  # (n_capas + 1, 2) estados (p, w) en cada interfaz
    media: tuple[Material, ...]
    thicknesses: np.ndarray
    R: complex
    T: complex


def _propagate(frequency: float, media, thicknesses, p_inc: float = 1.0) -> _Internal:
    media_t, th = _validate(media, thicknesses)
    if not frequency > 0.0:
        raise ValueError("la frecuencia debe ser > 0")

    k = np.array([m.wavenumber(frequency) for m in media_t], dtype=complex)
    rho = np.array([m.rho for m in media_t], dtype=float)
    k_left, k_right = k[0], k[-1]
    q_left = 1j * k_left / rho[0]
    q_right = 1j * k_right / rho[-1]

    # Producto de matrices de transferencia desde x = 0 hasta la interfaz j.
    n_layers = th.size
    cumulative = np.empty((n_layers + 1, 2, 2), dtype=complex)
    cumulative[0] = np.eye(2, dtype=complex)
    for j in range(n_layers):
        cumulative[j + 1] = layer_matrix(k[j + 1], rho[j + 1], th[j]) @ cumulative[j]
    m_total = cumulative[-1]

    length = float(th.sum())

    # Incógnitas (R, P) con P = amplitud del estado transmitido en x = L.
    #   fila 1:  m11·p0 + m12·w0 = P
    #   fila 2:  m21·p0 + m22·w0 = q_right·P
    # con p0 = p_inc·(1 + R) y w0 = q_left·p_inc·(1 − R).
    a1 = m_total[0, 0] - m_total[0, 1] * q_left  # coeficiente de R en fila 1
    b1 = m_total[0, 0] + m_total[0, 1] * q_left
    a2 = m_total[1, 0] - m_total[1, 1] * q_left
    b2 = m_total[1, 0] + m_total[1, 1] * q_left

    det = a2 - q_right * a1
    R = (b1 * q_right - b2) / det
    P = (a2 * b1 - a1 * b2) / det

    # Estado en x = 0 y estados en cada interfaz.
    p0 = p_inc * (1.0 + R)
    w0 = q_left * p_inc * (1.0 - R)
    states = cumulative @ np.array([p0, w0], dtype=complex)

    # T referido a x = 0: p(x>L) = T·p_inc·exp(i·k_right·x) ⇒ en x = L vale T·p_inc·exp(i·k·L).
    T = P / (p_inc * np.exp(1j * k_right * length))

    x_iface = np.concatenate([[0.0], np.cumsum(th)])
    return _Internal(
        k=k,
        x_iface=x_iface,
        states=states,
        media=media_t,
        thicknesses=th,
        R=complex(R),
        T=complex(T),
    )


def stack_response(frequency: float, media, thicknesses) -> StackResponse:
    """Coeficientes ``R`` y ``T`` de un stack estratificado con incidencia normal.

    Parameters
    ----------
    frequency : frecuencia [Hz].
    media : secuencia de :class:`Material` de longitud ``len(thicknesses) + 2``
        (semiespacio izquierdo, capas, semiespacio derecho).
    thicknesses : espesores de las capas [m].

    Returns
    -------
    StackResponse con ``R`` y ``T`` complejos referidos a la interfaz izquierda ``x = 0``.
    """
    internal = _propagate(frequency, media, thicknesses, p_inc=1.0)
    return StackResponse(
        R=internal.R,
        T=internal.T,
        frequency=float(frequency),
        total_thickness=float(internal.thicknesses.sum()),
        k=tuple(complex(v) for v in internal.k),
        impedance=tuple(float(m.impedance) for m in internal.media),
    )


def state_at(x: float | np.ndarray, frequency: float, media, thicknesses) -> np.ndarray:
    """Estado ``(p, w)`` complejo en ``x`` (una fila por punto, forma ``(..., 2)``)."""
    internal = _propagate(frequency, media, thicknesses, p_inc=1.0)
    xs = np.atleast_1d(np.asarray(x, dtype=float))
    out = np.empty((xs.size, 2), dtype=complex)

    length = float(internal.thicknesses.sum())
    k = internal.k
    rho = np.array([m.rho for m in internal.media], dtype=float)
    q_left = 1j * k[0] / rho[0]
    q_right = 1j * k[-1] / rho[-1]

    left = xs < 0.0
    right = xs > length if length > 0.0 else xs >= 0.0
    inside = ~(left | right)

    # Semiespacio izquierdo: incidente + reflejada.
    if np.any(left):
        xl = xs[left]
        p_inc_l = np.exp(1j * k[0] * xl)
        p_out = p_inc_l + internal.R * np.exp(-1j * k[0] * xl)
        w_out = q_left * (p_inc_l - internal.R * np.exp(-1j * k[0] * xl))
        out[left] = np.column_stack([p_out, w_out])

    # Semiespacio derecho: transmitida referida a x = 0.
    if np.any(right):
        xr = xs[right]
        p_out = internal.T * np.exp(1j * k[-1] * xr)
        out[right] = np.column_stack([p_out, q_right * p_out])

    # Capas: propago desde la interfaz j-ésima con la matriz de la capa.
    if np.any(inside):
        xi = xs[inside]
        idx = np.flatnonzero(inside)
        layer = np.searchsorted(internal.x_iface, xi, side="right") - 1
        layer = np.clip(layer, 0, internal.thicknesses.size - 1)
        p_out = np.empty(xi.size, dtype=complex)
        w_out = np.empty(xi.size, dtype=complex)
        for j in np.unique(layer):
            sel = layer == j
            d = xi[sel] - internal.x_iface[j]
            kk = k[j + 1]
            p_j, w_j = internal.states[j]
            p_out[sel] = p_j * np.cos(kk * d) + w_j * rho[j + 1] * np.sin(kk * d) / kk
            w_out[sel] = -(kk / rho[j + 1]) * p_j * np.sin(kk * d) + w_j * np.cos(kk * d)
        out[idx, 0] = p_out
        out[idx, 1] = w_out

    return out


def pressure(
    x_query,
    frequency: float,
    media,
    thicknesses,
    p_inc: float = 1.0,
) -> np.ndarray:
    """Campo de presión complejo en ``x_query`` (array o escalar) para un stack estratificado.

    La onda incidente es ``p_inc·exp(+i·k₁·x)`` (viaja hacia +x). Devuelve un array de la misma
    forma que ``x_query``.
    """
    xs = np.asarray(x_query, dtype=float)
    flat = np.atleast_1d(xs).reshape(-1)
    # Propagación a amplitud incidente unitaria y escalado al final: así no se cuenta dos veces
    # p_inc (los estados internos ya vienen escalados por la amplitud incidente).
    internal = _propagate(frequency, media, thicknesses, p_inc=1.0)
    values = _evaluate(internal, flat) * p_inc
    if np.ndim(x_query) == 0:
        return np.asarray(values[0])
    return values.reshape(xs.shape)


def _evaluate(internal: _Internal, xs: np.ndarray) -> np.ndarray:
    """Presión compleja para amplitud incidente unitaria."""
    k = internal.k
    rho = np.array([m.rho for m in internal.media], dtype=float)
    length = float(internal.thicknesses.sum())
    out = np.empty(xs.size, dtype=complex)

    left = xs < 0.0
    if np.any(left):
        xl = xs[left]
        out[left] = np.exp(1j * k[0] * xl) + internal.R * np.exp(-1j * k[0] * xl)

    if length > 0.0:
        right = xs > length
    else:
        right = xs >= 0.0
    if np.any(right):
        out[right] = internal.T * np.exp(1j * k[-1] * xs[right])

    inside = ~(left | right)
    if np.any(inside):
        xi = xs[inside]
        layer = np.searchsorted(internal.x_iface, xi, side="right") - 1
        layer = np.clip(layer, 0, internal.thicknesses.size - 1)
        local = np.empty(xi.size, dtype=complex)
        for j in np.unique(layer):
            sel = layer == j
            d = xi[sel] - internal.x_iface[j]
            kk = k[j + 1]
            p_j, w_j = internal.states[j]
            local[sel] = p_j * np.cos(kk * d) + w_j * rho[j + 1] * np.sin(kk * d) / kk
        out[inside] = local

    return out
