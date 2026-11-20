"""Dispersión acústica por un cilindro infinito en 2D: serie exacta de referencia.

Convención del repo (``docs/decisions.md`` D4), obligatoria aquí porque todo lo demás se mide
contra esta solución:

- tiempo ``exp(-i·ω·t)``;
- Helmholtz ``∇²p + k²p = 0``;
- ``k = ω/c + i·α`` con parte imaginaria **positiva**;
- ondas salientes ``∝ exp(+i·k·r)`` ⇒ en 2D la onda saliente es ``H₀⁽¹⁾(k·r)``;
- unidades SI; la API pública recibe Hz y devuelve Pa.

Geometría
---------
Cilindro circular infinito de radio ``a`` centrado en el origen, eje paralelo a ``z``, onda plana
incidente que viaja en el plano ``xy`` (incidencia perpendicular al eje, el caso 2D puro)::

    p_inc(x, y) = P0 · exp(i·k·(x·cos α + y·sin α))

Serie
-----
Con ``x = r·cos θ``, ``y = r·sin θ`` la identidad de Jacobi-Anger con esta convención es::

    exp(i·k·r·cos(θ − α)) = Σ_n iⁿ · J_n(k·r) · exp(i·n·(θ − α))

y el campo total fuera del cilindro (``r ≥ a``) se escribe término a término::

    p_total = Σ_n [ iⁿ·J_n(k·r) + A_n·H_n⁽¹⁾(k·r) ] · exp(i·n·(θ − α))

Los coeficientes ``A_n`` salen de la condición de contorno en ``r = a``; con ``n̂ = r̂`` (normal
saliente del cuerpo, que apunta **hacia el fluido**):

- ``boundary="rigid"`` (Neumann, pared rígida) ⇒ ``∂p/∂r|_{r=a} = 0``::

      A_n = −iⁿ · J_n′(k·a) / H_n⁽¹⁾′(k·a)

- ``boundary="pressure_release"`` (Dirichlet, presión nula) ⇒ ``p|_{r=a} = 0``::

      A_n = −iⁿ · J_n(k·a) / H_n⁽¹⁾(k·a)

La derivada de las funciones de Bessel se calcula con ``f_n′ = (f_{n−1} − f_{n+1})/2``, válida para
``J`` y para ``H⁽¹⁾``.

Con ``Im(k) > 0`` los argumentos ``k·r`` son complejos: ``J_n`` y ``H_n⁽¹⁾`` de SciPy los aceptan
directamente (``k`` entra analíticamente en toda la construcción, la atenuación es solo un ``k``
complejo). El decaimiento ``exp(−α·r)`` aparece solo, no hay que agregarlo a mano.

Convergencia
------------
La serie se trunca en ``|n| ≤ N``. Si no se pasa ``n_terms``, ``N`` se duplica hasta que la
contribución del último par de órdenes (``±N``) caiga por debajo de ``1e-12`` relativo a la amplitud
incidente; el ``N`` usado se devuelve en ``result["n_terms"]``.
"""

from __future__ import annotations

import numpy as np
from scipy.special import hankel1, jv

from ..materials import Material

__all__ = [
    "BOUNDARY_CONDITIONS",
    "cylinder_scattering",
    "cylinder_scattered_field",
    "cylinder_radial_derivative",
]

#: Condiciones de contorno soportadas en la superficie del cilindro.
BOUNDARY_CONDITIONS = ("rigid", "pressure_release")

#: Tolerancia relativa por defecto para elegir ``n_terms`` automáticamente.
DEFAULT_REL_TOL = 1e-12

#: Techo de seguridad del truncado (``ka`` realista nunca lo necesita: N ≈ ka + 15).
_MAX_TERMS = 512


# --------------------------------------------------------------------------- #
# utilidades internas
# --------------------------------------------------------------------------- #


def _as_points(points) -> np.ndarray:
    """Convierte ``points`` a un array ``(M, 2)`` de floats."""
    pts = np.asarray(points, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 2:
        raise ValueError(f"points debe tener forma (M, 2); recibido {pts.shape}")
    if not np.all(np.isfinite(pts)):
        raise ValueError("points contiene valores no finitos")
    return pts


def _check_boundary(boundary: str) -> str:
    if boundary not in BOUNDARY_CONDITIONS:
        raise ValueError(
            f"boundary debe ser uno de {BOUNDARY_CONDITIONS}; recibido {boundary!r}"
        )
    return boundary


def _check_arguments(frequency, material, radius) -> complex:
    """Valida y devuelve ``k`` [1/m] del medio a ``frequency``."""
    if not isinstance(material, Material):
        raise TypeError("material debe ser un materials.Material")
    if not (radius > 0.0):
        raise ValueError(f"radius debe ser > 0; recibido {radius!r}")
    if not (frequency >= 0.0):
        raise ValueError(f"frequency debe ser >= 0; recibido {frequency!r}")
    return material.wavenumber(frequency)


def _orders(n_terms: int) -> np.ndarray:
    return np.arange(-n_terms, n_terms + 1)


def _coefficients(k: complex, radius: float, n_terms: int, boundary: str):
    """Coeficientes ``A_n`` de la serie dispersada para ``|n| ≤ n_terms``.

    Devuelve ``(orders, A)``. Los órdenes donde el cociente es numéricamente indefinido
    (``H_n⁽¹⁾′`` desborda para ``|n| ≫ |k·a|``) se anulan: su contribución verdadera es
    despreciable frente a ``1e-200`` y solo aparecen si el usuario pide un truncado absurdo.
    """
    n = _orders(n_terms)
    z = k * radius

    if boundary == "rigid":
        num = 0.5 * (jv(n - 1, z) - jv(n + 1, z))          # J_n′(k·a)
        den = 0.5 * (hankel1(n - 1, z) - hankel1(n + 1, z))  # H_n⁽¹⁾′(k·a)
    else:  # pressure_release
        num = jv(n, z)
        den = hankel1(n, z)

    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        ratio = num / den
    ratio = np.where(np.isfinite(ratio), ratio, 0.0)
    return n, -np.power(1j, n) * ratio


def _tail_magnitude(k: complex, radius: float, r_values: np.ndarray,
                    n_terms: int, boundary: str) -> float:
    """Magnitud del aporte del par de órdenes ``±n_terms`` en los puntos de mayor ``r``.

    ``|A_{-n}| = |A_n|`` (``J_{-n} = (−1)ⁿJ_n`` y ``H_{-n}⁽¹⁾ = (−1)ⁿH_n⁽¹⁾``), así que alcanza con
    mirar el orden positivo.
    """
    _, coeff = _coefficients(k, radius, n_terms, boundary)
    amplitude = abs(coeff[-1])  # orden +n_terms
    z = k * r_values
    with np.errstate(invalid="ignore", over="ignore"):
        sca = amplitude * np.abs(hankel1(n_terms, z))
        inc = np.abs(jv(n_terms, z))
    sca = np.where(np.isfinite(sca), sca, np.inf)
    return float(max(np.max(sca), np.max(inc)))


def _auto_n_terms(k: complex, radius: float, r_max: float, boundary: str,
                  rel_tol: float, amplitude: float) -> int:
    """Elige ``n_terms`` duplicando hasta que el último par de órdenes aporte < ``rel_tol``."""
    r_values = np.array([max(r_max, radius)])
    n_terms = max(8, int(np.ceil(abs(k) * radius)) + 6)
    while n_terms <= _MAX_TERMS:
        if _tail_magnitude(k, radius, r_values, n_terms, boundary) <= rel_tol * amplitude:
            return n_terms
        n_terms *= 2
    raise RuntimeError(
        f"la serie no converge con {_MAX_TERMS} términos (ka={abs(k) * radius:.3g}); "
        "revisar los parámetros"
    )


def _evaluate_series(k: complex, n: np.ndarray, coefficients: np.ndarray,
                     r: np.ndarray, theta: np.ndarray, angle: float):
    """Evalúa sobre ``(r, θ)`` la serie incidente y la dispersada (Jacobi-Anger + Hankel)."""
    z = k * r
    phase = np.exp(1j * np.outer(n, theta - angle))
    p_inc = (np.power(1j, n)[:, None] * jv(n[:, None], z[None, :]) * phase).sum(axis=0)
    with np.errstate(invalid="ignore", over="ignore"):
        terms = coefficients[:, None] * hankel1(n[:, None], z[None, :]) * phase
    # un 0·∞ aislado no debe contaminar la suma (ver _coefficients)
    p_sca = np.where(np.isfinite(terms), terms, 0.0).sum(axis=0)
    return p_inc, p_sca


# --------------------------------------------------------------------------- #
# API pública
# --------------------------------------------------------------------------- #


def cylinder_scattering(frequency: float, material: Material, radius: float, points,
                        angle: float = 0.0, boundary: str = "rigid",
                        p_inc: float = 1.0, n_terms: int | None = None) -> dict:
    """Campo total (incidente + dispersado) de un cilindro infinito en ``points``.

    Parámetros
    ----------
    frequency : frecuencia [Hz].
    material : medio exterior (``materials.Material``); de ahí salen ``c`` y ``α(f)``.
    radius : radio del cilindro [m].
    points : ``(M, 2)`` posiciones de evaluación [m], con ``r ≥ radius``.
    angle : ángulo de incidencia ``α`` [rad]; la onda viaja según ``(cos α, sin α)``.
    boundary : ``"rigid"`` (Neumann, ``∂p/∂n = 0``) o ``"pressure_release"`` (``p = 0``).
    p_inc : amplitud ``P0`` de la onda plana incidente [Pa].
    n_terms : truncado ``|n| ≤ n_terms``; ``None`` lo elige automáticamente.

    Devuelve
    --------
    dict con
        ``total``      : ``(M,)`` complejo, ``p_inc + p_sca`` [Pa]
        ``incident``   : ``(M,)`` complejo, onda plana incidente [Pa]
        ``scattered``  : ``(M,)`` complejo, campo dispersado [Pa]
        ``k``          : número de onda complejo [1/m] usado
        ``n_terms``    : truncado efectivo
        ``orders``     : ``(2N+1,)`` órdenes ``n``
        ``coefficients`` : ``(2N+1,)`` coeficientes ``A_n``
        ``boundary``   : condición de contorno usada
    """
    boundary = _check_boundary(boundary)
    k = _check_arguments(frequency, material, radius)
    pts = _as_points(points)

    r = np.hypot(pts[:, 0], pts[:, 1])
    if np.any(r < radius * (1.0 - 1e-12)):
        raise ValueError(
            "hay puntos dentro del cilindro (r < radius): la serie solo vale en el exterior "
            f"(r mínimo = {r.min():.6g} m, radius = {radius:.6g} m)"
        )
    theta = np.arctan2(pts[:, 1], pts[:, 0])

    scale = abs(float(p_inc)) if p_inc != 0.0 else 1.0
    if n_terms is None:
        n_terms = _auto_n_terms(k, radius, float(r.max()), boundary, DEFAULT_REL_TOL, scale)
    elif n_terms < 0:
        raise ValueError(f"n_terms debe ser >= 0; recibido {n_terms!r}")

    n, coeff = _coefficients(k, radius, int(n_terms), boundary)
    p_i, p_s = _evaluate_series(k, n, coeff, r, theta, angle)

    return {
        "total": p_inc * (p_i + p_s),
        "incident": p_inc * p_i,
        "scattered": p_inc * p_s,
        "k": k,
        "n_terms": int(n_terms),
        "orders": n,
        "coefficients": coeff,
        "boundary": boundary,
    }


def cylinder_scattered_field(frequency: float, material: Material, radius: float, points,
                             angle: float = 0.0, boundary: str = "rigid",
                             p_inc: float = 1.0, n_terms: int | None = None) -> np.ndarray:
    """Solo el campo dispersado ``p_sca`` en ``points``; atajo de :func:`cylinder_scattering`."""
    result = cylinder_scattering(frequency, material, radius, points, angle=angle,
                                 boundary=boundary, p_inc=p_inc, n_terms=n_terms)
    return result["scattered"]


def cylinder_radial_derivative(frequency: float, material: Material, radius: float, points,
                               angle: float = 0.0, boundary: str = "rigid",
                               p_inc: float = 1.0, n_terms: int | None = None) -> np.ndarray:
    r"""``∂p_total/∂r`` de la serie, evaluado en ``points`` (``r ≥ radius``).

    Sirve para verificar la condición de contorno: sobre la superficie de un cilindro rígido el
    campo total debe cumplir ``∂p/∂r|_{r=a} = 0``. Con ``n̂ = r̂`` esta es exactamente la velocidad
    normal (``∂p/∂n = ρ·ω²·u_n`` en el fluido), así que ``|∂p/∂n| / (|k|·|p_inc|)`` es el residual
    adimensional de la condición de Neumann.

    Derivando la serie término a término (la derivada conmuta con la suma porque ambas series
    convergen uniformemente en ``r ≥ radius``)::

        ∂p/∂r = Σ_n [ iⁿ·k·J_n′(k·r) + A_n·k·H_n⁽¹⁾′(k·r) ] · exp(i·n·(θ − α))
    """
    boundary = _check_boundary(boundary)
    k = _check_arguments(frequency, material, radius)
    pts = _as_points(points)

    r = np.hypot(pts[:, 0], pts[:, 1])
    if np.any(r < radius * (1.0 - 1e-12)):
        raise ValueError("hay puntos dentro del cilindro (r < radius)")
    theta = np.arctan2(pts[:, 1], pts[:, 0])

    scale = abs(float(p_inc)) if p_inc != 0.0 else 1.0
    if n_terms is None:
        n_terms = _auto_n_terms(k, radius, float(r.max()), boundary, DEFAULT_REL_TOL, scale)
    elif n_terms < 0:
        raise ValueError(f"n_terms debe ser >= 0; recibido {n_terms!r}")

    n, coeff = _coefficients(k, radius, int(n_terms), boundary)
    z = k * r
    j_prime = 0.5 * (jv(n[:, None] - 1, z[None, :]) - jv(n[:, None] + 1, z[None, :]))
    h_prime = 0.5 * (hankel1(n[:, None] - 1, z[None, :])
                     - hankel1(n[:, None] + 1, z[None, :]))
    phase = np.exp(1j * np.outer(n, theta - angle))
    terms = k * (np.power(1j, n)[:, None] * j_prime + coeff[:, None] * h_prime) * phase
    return p_inc * np.where(np.isfinite(terms), terms, 0.0).sum(axis=0)
