"""Dispersión de una onda plana por una esfera **elástica** sumergida: serie de Faran (1951).

Referencia: Faran, *Sound Scattering by Solid Cylinders and Spheres*, J. Acoust. Soc. Am. 23:405–418
(1951). Es la versión con ondas P y S dentro del sólido de la serie de ondas parciales; la de esfera
fluida es Anderson (1950) y está en `analytic/sphere.py`.

Derivación
----------
Sólido (``r < a``): ``λ, μ, ρ₁``, ``k_α = ω/c_p``, ``k_β = ω/c_s``, con ``c_p² = (λ+2μ)/ρ₁`` y
``c_s² = μ/ρ₁``. Fluido (``r > a``): ``ρ₀, c₀, k₀``. Convención ``e^{−iωt}``.

Potenciales, con ``Pₙ = Pₙ(cosθ)``:

    Φ = Σ Aₙ jₙ(k_α r) Pₙ              (longitudinal)
    Ψ = Σ Bₙ jₙ(k_β r) Pₙ              (corte)
    u = ∇Φ + ∇×(∇×(Ψ r̂))

En el fluido, ``p = p_i + p_s`` con ``p_i = Σ jₙ(k₀r)Pₙ`` y ``p_s = Σ Cₙ hₙ⁽¹⁾(k₀r)Pₙ`` (se omite el
factor ``(2n+1)iⁿ``, que es común y se repone al final).

Los desplazamientos y las tensiones que usa el sistema salieron de **sympy**, calculando el
rotacional-rotacional en esféricas con simetría axial y las deformaciones
``ε_rr, ε_θθ, ε_φφ, ε_rθ`` con ``σ = λ·tr(ε)I + 2μ·ε``:

    u_r = Φ_r + n(n+1)·Ψ/r²
    u_θ = (Φ_θ + Ψ_rθ)/r

    σ_rr/Pₙ = −λn(n+1)f/r² + 2λ f′/r + (λ+2μ)f″ − 4μn(n+1)g/r³ + 2μn(n+1)g′/r²
    σ_rθ/Pₙ′ = −2μ f/r² + 2μ f′/r + μn(n+1)g/r³ − 2μ g′/r² + μ g″/r

con ``f(r) = jₙ(k_α r)`` y ``g(r) = jₙ(k_β r)`` (las primas son derivadas en ``r``). En ``σ_rr`` se usó
la ecuación de Legendre ``P″ + cotθ·P′ = −n(n+1)P`` para que quede proporcional a ``Pₙ``, y en
``σ_rθ`` para que quede proporcional a ``Pₙ′``.

Condiciones en ``r = a`` (tres ecuaciones por orden ``n``, una más que en Anderson):

    1. u_r(sólido) = u_r(fluido),  con u_r(fluido) = (1/(ρ₀ω²))·∂p/∂r
    2. σ_rr = −p
    3. σ_rθ = 0

Validación (sin tabla de referencia)
------------------------------------
- **Teorema óptico**: para un dispersor sin pérdidas ``σ_ext = σ_sca``, con
  ``σ_ext = (4π/k₀)·Im A(0)`` y ``σ_sca = ∫|A|²dΩ``. Es conservación de energía: si los coeficientes
  están mal, no cierra. Ver `amplitud_lejana`.
- **Residuos** de las tres condiciones en ``r = a`` después de resolver.
- **Contraste bajo**: con la esfera casi igual al fluido, ``Cₙ → 0``.

No se valida contra el límite ``μ → 0`` de `sphere.py`: es físicamente correcto pero
**numéricamente singular** (``k_β → ∞`` y ``jₙ(k_β a)`` oscila sin límite).
"""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.polynomial import legendre as _leg
from scipy.special import spherical_jn, spherical_yn

from ..materials import ElasticSolid, Material

_I = 1j


# --------------------------------------------------------------------------- funciones esféricas


def _h1(n: int, x: complex) -> complex:
    return complex(spherical_jn(n, x) + _I * spherical_yn(n, x))


def _h1_deriv(n: int, x: complex) -> complex:
    """``hₙ⁽¹⁾′(x) = hₙ₋₁⁽¹⁾(x) − (n+1)/x·hₙ⁽¹⁾(x)``, con ``h₀′ = −h₁`` para ``n = 0``.

    Recurrencia y no diferencia finita: a ``k₀a`` grande la diferencia finita pierde precisión justo
    donde el término de corte domina.
    """
    if n == 0:
        return -_h1(1, x)
    return _h1(n - 1, x) - (n + 1) / x * _h1(n, x)


def _j_deriv(n: int, x: complex) -> complex:
    """``jₙ′(x) = jₙ₋₁(x) − (n+1)/x·jₙ(x)``, con ``j₀′ = −j₁`` para ``n = 0``.

    El caso ``n = 0`` va aparte porque la recurrencia pide ``j₋₁``, que no existe: `scipy` devuelve
    inf ahí y el NaN se propaga a todo el sistema.
    """
    if n == 0:
        return complex(-spherical_jn(1, x))
    return complex(spherical_jn(n - 1, x) - (n + 1) / x * spherical_jn(n, x))


def _j_second(n: int, x: complex) -> complex:
    """``jₙ″`` desde la ecuación de Bessel esférica: ``j″ = −(2/x)j′ − (1 − n(n+1)/x²)j``."""
    return -(2.0 / x) * _j_deriv(n, x) - (1.0 - n * (n + 1) / (x * x)) * complex(spherical_jn(n, x))


def _n_terms(k0: complex, radio: float) -> int:
    ka = abs(k0) * radio
    return int(min(max(int(ka + 4.0 * ka ** (1.0 / 3.0) + 2.0) + 12, 12), 4000))


def _coeficiente(
    n: int, radio: float, k0: complex, k_alpha: complex, k_beta: complex,
    rho0: float, omega: float, lam: float, mu: float,
) -> tuple[complex, complex, complex]:
    """Resuelve el 3×3 de ``(Aₙ, Bₙ, Cₙ)`` para un orden ``n``."""
    r = radio
    xa, xb, x0 = k_alpha * r, k_beta * r, k0 * r

    f = complex(spherical_jn(n, xa))
    f1 = k_alpha * _j_deriv(n, xa)
    f2 = k_alpha**2 * _j_second(n, xa)
    g = complex(spherical_jn(n, xb))
    g1 = k_beta * _j_deriv(n, xb)
    g2 = k_beta**2 * _j_second(n, xb)

    jn = complex(spherical_jn(n, x0))
    jn1 = k0 * _j_deriv(n, x0)
    hn = _h1(n, x0)
    hn1 = k0 * _h1_deriv(n, x0)

    nn = n * (n + 1)
    # fila 1: u_r del sólido menos u_r del fluido (dividido por Pₙ).
    # OJO: `hn1` ya es ∂h/∂r = k₀·h′(k₀a); el factor 1/(ρ₀ω²) va una sola vez.
    # Poner `k0/(ρ₀ω²)·hn1` duplica k₀ y el sistema entero sale como si la esfera fuera rígida.
    fila1 = np.array([f1, nn * g / r**2, -hn1 / (rho0 * omega**2)], dtype=complex)
    # fila 2: σ_rr del sólido + p (dividido por Pₙ)
    fila2 = np.array([
        -lam * nn * f / r**2 + 2.0 * lam * f1 / r + (lam + 2.0 * mu) * f2,
        -4.0 * mu * nn * g / r**3 + 2.0 * mu * nn * g1 / r**2,
        hn,
    ], dtype=complex)
    # fila 3: σ_rθ del sólido (dividido por Pₙ′)
    fila3 = np.array([
        -2.0 * mu * f / r**2 + 2.0 * mu * f1 / r,
        mu * nn * g / r**3 - 2.0 * mu * g1 / r**2 + mu * g2 / r,
        0.0,
    ], dtype=complex)

    matriz = np.vstack([fila1, fila2, fila3])
    # término independiente: lo que aporta la incidente sola
    rhs = np.array([
        jn1 / (rho0 * omega**2),
        -jn,
        0.0,
    ], dtype=complex)
    sol = np.linalg.solve(matriz, rhs)
    return complex(sol[0]), complex(sol[1]), complex(sol[2])


# --------------------------------------------------------------------------- API pública


def _check(frequency: float, radio: float, medio, solido) -> None:
    if not np.isfinite(frequency) or frequency <= 0.0:
        raise ValueError(f"frequency tiene que ser > 0; llegó {frequency!r}")
    if not np.isfinite(radio) or radio <= 0.0:
        raise ValueError(f"radius tiene que ser > 0; llegó {radio!r}")
    if not isinstance(medio, Material):
        raise TypeError("medium tiene que ser un Material")
    if not isinstance(solido, ElasticSolid):
        raise TypeError("sphere tiene que ser un ElasticSolid (la esfera es sólida)")


def _propiedades(solido: ElasticSolid) -> tuple[float, float]:
    lam = solido.rho * (solido.c_p**2 - 2.0 * solido.c_s**2)
    mu = solido.rho * solido.c_s**2
    if mu <= 0:
        raise ValueError("el sólido tiene c_s <= 0: no es un sólido elástico")
    return float(lam), float(mu)


def coeficientes(
    frequency: float, medio: Material, solido: ElasticSolid, radio: float,
) -> dict[str, Any]:
    """``Cₙ`` (dispersado), ``Aₙ``, ``Bₙ`` y ``ka`` para la esfera elástica."""
    _check(frequency, radio, medio, solido)
    lam, mu = _propiedades(solido)
    k0 = complex(medio.wavenumber(frequency))
    omega = 2.0 * np.pi * frequency
    k_alpha = omega / solido.c_p + 1j * abs(complex(medio.wavenumber(frequency)).imag)
    k_beta = omega / solido.c_s

    n_max = _n_terms(k0, radio)
    A = np.zeros(n_max + 1, dtype=complex)
    B = np.zeros(n_max + 1, dtype=complex)
    C = np.zeros(n_max + 1, dtype=complex)
    for n in range(n_max + 1):
        A[n], B[n], C[n] = _coeficiente(
            n, radio, k0, k_alpha, k_beta, float(medio.rho), omega, lam, mu
        )
    return {
        "C": C, "A": A, "B": B, "ka": float(abs(k0) * radio), "n_terms": n_max,
        "k0": k0, "k_alpha": k_alpha, "k_beta": k_beta, "lambda": lam, "mu": mu,
    }


def amplitud_lejana(coef: dict[str, Any], theta, radio: float) -> np.ndarray:
    """Amplitud de dispersión ``A(θ)``: ``p_s ≈ A(θ)·e^{ik₀r}/r``."""
    C = np.asarray(coef["C"])
    k0 = coef["k0"]
    th = np.atleast_1d(np.asarray(theta, dtype=float))
    salida = np.zeros(th.size, dtype=complex)
    for n, cn in enumerate(C):
        if cn == 0.0:
            continue
        coef_leg = np.zeros(n + 1)
        coef_leg[n] = 1.0
        salida += (2 * n + 1) * cn * _leg.legval(np.cos(th), coef_leg)
    return salida / (_I * k0)


def secciones(coef: dict[str, Any], radio: float, n_theta: int = 360) -> dict[str, float]:
    """``σ_sca`` (integrada) y ``σ_ext`` (teorema óptico). Tienen que coincidir si no hay pérdidas."""
    k0 = coef["k0"]
    th = np.linspace(0.0, np.pi, n_theta)
    amp = amplitud_lejana(coef, th, radio)
    peso = np.sin(th)
    sigma_sca = float(2.0 * np.pi * np.trapezoid(np.abs(amp) ** 2 * peso, th))
    amp_0 = amplitud_lejana(coef, np.array([0.0]), radio)[0]
    sigma_ext = float(4.0 * np.pi / abs(k0) * np.imag(amp_0))
    return {
        "sigma_sca": sigma_sca,
        "sigma_ext": sigma_ext,
        "desvio_relativo": abs(sigma_ext - sigma_sca) / abs(sigma_sca) if sigma_sca else float("inf"),
    }


def residuos(
    frequency: float, medio: Material, solido: ElasticSolid, radio: float, theta: float = 1.1,
) -> dict[str, float]:
    """Residuos relativos de las tres condiciones de borde, después de resolver.

    Se evalúan en un ángulo genérico: con simetría axial todos los términos son ``∝ Pₙ`` o
    ``∝ Pₙ′``, así que el sistema queda determinado por un solo punto.
    """
    lam, mu = _propiedades(solido)
    k0 = complex(medio.wavenumber(frequency))
    omega = 2.0 * np.pi * frequency
    datos = coeficientes(frequency, medio, solido, radio)
    k_alpha, k_beta = datos["k_alpha"], datos["k_beta"]
    A, B, C = datos["A"], datos["B"], datos["C"]

    r = radio
    x = np.cos(theta)
    res = np.zeros(3, dtype=complex)
    esc = np.zeros(3)
    # escala de σ_rθ: NO se normaliza por sí misma (esa razón da ~1 aunque el residuo sea 1e-17,
    # que fue lo que me confundió la primera vez). Se usa la escala de la presión y la tensión
    # normal, que son las que el sistema iguala.
    escala_normal = 0.0
    for n in range(datos["n_terms"] + 1):
        coef_leg = np.zeros(n + 1)
        coef_leg[n] = 1.0
        P = float(_leg.legval(x, coef_leg))
        dP = float(-np.sqrt(1.0 - x * x) * _leg.legval(x, _leg.legder(coef_leg)))

        xa, xb, x0 = k_alpha * r, k_beta * r, k0 * r
        f = complex(spherical_jn(n, xa))
        f1 = k_alpha * _j_deriv(n, xa)
        f2 = k_alpha**2 * _j_second(n, xa)
        g = complex(spherical_jn(n, xb))
        g1 = k_beta * _j_deriv(n, xb)
        g2 = k_beta**2 * _j_second(n, xb)
        jn = complex(spherical_jn(n, x0))
        jn1 = k0 * _j_deriv(n, x0)
        hn = _h1(n, x0)
        hn1 = k0 * _h1_deriv(n, x0)
        nn = n * (n + 1)

        u_r_sol = (A[n] * f1 + B[n] * nn * g / r**2) * P
        u_r_flu = (jn1 + C[n] * hn1) / (float(medio.rho) * omega**2) * P
        sig_rr = (A[n] * (-lam * nn * f / r**2 + 2.0 * lam * f1 / r + (lam + 2.0 * mu) * f2)
                  + B[n] * (-4.0 * mu * nn * g / r**3 + 2.0 * mu * nn * g1 / r**2)) * P
        p = (jn + C[n] * hn) * P
        sig_rt = (A[n] * (-2.0 * mu * f / r**2 + 2.0 * mu * f1 / r)
                  + B[n] * (mu * nn * g / r**3 - 2.0 * mu * g1 / r**2 + mu * g2 / r)) * dP

        res[0] += u_r_sol - u_r_flu
        res[1] += sig_rr + p
        res[2] += sig_rt
        esc[0] += abs(u_r_flu) + abs(u_r_sol)
        esc[1] += abs(p) + abs(sig_rr)
        escala_normal += abs(p) + abs(sig_rr)
    return {
        "u_r": float(abs(res[0]) / max(esc[0], 1e-300)),
        "sigma_rr": float(abs(res[1]) / max(esc[1], 1e-300)),
        "sigma_rt": float(abs(res[2]) / max(escala_normal, 1e-300)),
    }

def campo_dispersado(
    coef: dict[str, Any], points, radio: float, direccion=(0.0, 0.0, 1.0)
) -> np.ndarray:
    """Presión dispersada **exacta** en puntos arbitrarios (no sólo en el campo lejano).

        p_s(x) = Σₙ (2n+1)·iⁿ·Cₙ·hₙ⁽¹⁾(k₀r)·Pₙ(cosθ)

    Se evalúa la serie con las Hankel completas en vez de la asíntota ``e^{ikr}/r``: si se compara
    contra el FEM en ``r = 1.5a``, la asíntota introduce un error de orden ``1/(k₀r)`` que se
    confundiría con error de discretización.
    """
    C = np.asarray(coef["C"])
    k0 = coef["k0"]
    d = np.asarray(direccion, dtype=float)
    d = d / np.linalg.norm(d)
    x = np.atleast_2d(np.asarray(points, dtype=float))
    r = np.linalg.norm(x, axis=1)
    cos_t = np.clip(x @ d / np.where(r > 0, r, 1.0), -1.0, 1.0)
    salida = np.zeros(x.shape[0], dtype=complex)
    for n, cn in enumerate(C):
        if cn == 0.0:
            continue
        coef_leg = np.zeros(n + 1)
        coef_leg[n] = 1.0
        # hₙ⁽¹⁾ vectorizada: `_h1` devuelve un escalar de Python y no acepta arrays de r
        hankel = spherical_jn(n, k0 * r) + _I * spherical_yn(n, k0 * r)
        salida += (2 * n + 1) * (_I ** n) * cn * hankel * _leg.legval(cos_t, coef_leg)
    return salida
