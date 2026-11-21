"""Dispersión de una onda plana por una esfera fluida: solución por ondas parciales.

Es la referencia analítica de **Anderson (1950)**, *J. Acoust. Soc. Am.* 22:426–431 ("Sound
Scattering from a Fluid Sphere"), la versión acústica de la serie de Mie. Fue el dato que faltaba
para validar la dispersión penetrable en 3D (D40); la esfera **elástica** sumergida es Faran (1951)
y queda como la etapa siguiente.

Planteo
-------
Onda plana ``p_i = exp(i·k₀·z)`` incidiendo sobre una esfera de radio ``a``, densidad ``ρ₁`` y
velocidad ``c₁``, sumergida en un fluido ``ρ₀, c₀``. Con la convención del repo (``e^{−iωt}``,
onda viajando hacia ``+z``):

    p_i   = Σ (2n+1) iⁿ jₙ(k₀r) Pₙ(cosθ)
    p_s   = Σ (2n+1) iⁿ Aₙ hₙ⁽¹⁾(k₀r) Pₙ(cosθ)
    p_int = Σ (2n+1) iⁿ Bₙ jₙ(k₁r) Pₙ(cosθ)

Condiciones en ``r = a``: continuidad de presión total y de flujo normal ``(1/ρ)∂_r p``:

    Rₙ + AₙHₙ = BₙJₙ
    (1/ρ₀)(Rₙ′ + AₙHₙ′) = (1/ρ₁) BₙJₙ′

con ``Rₙ = jₙ(k₀a)``, ``Rₙ′ = k₀jₙ′(k₀a)``, ``Hₙ = hₙ⁽¹⁾(k₀a)``, ``Hₙ′ = k₀hₙ⁽¹⁾′(k₀a)``,
``Jₙ = jₙ(k₁a)``, ``Jₙ′ = k₁jₙ′(k₁a)``. Despejando:

    Aₙ = (D·Rₙ − Rₙ′) / (Hₙ′ − D·Hₙ),     D = (ρ₁/ρ₀)·(Jₙ′/Jₙ)
    Bₙ = (Rₙ + AₙHₙ) / Jₙ

La suma no está restringida a ``ka ≪ 1``: eso es la aproximación de Rayleigh, no esta solución.
Lo único que hay que cuidar es truncar la serie con criterio (``_auto_n_terms``).

Qué se puede usar como test
---------------------------
- **Esfera casada** (``ρ₁ = ρ₀``, ``c₁ = c₀``): ``Aₙ = 0`` para todo ``n``, no hay dispersión.
- **Límite rígido** (``ρ₁ → ∞``): ``Aₙ → −Rₙ′/Hₙ′`` (``∂_r p = 0`` en la superficie).
- **Límite blando** (``c₁ → 0``): ``Aₙ → −Rₙ/Hₙ`` (``p = 0``).
- **Continuidad** de ``p`` y ``(1/ρ)∂_r p`` en ``r = a`` con los tres campos evaluados.
- **Teorema óptico**: la sección eficaz total (extinción) tiene que coincidir con la parte
  imaginaria de la amplitud hacia adelante.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy.special import h1vp, jvp, spherical_jn, spherical_yn  # noqa: F401

from ..materials import Material

_I = 1j


def _hankel1(n: np.ndarray, x: complex) -> np.ndarray:
    """``hₙ⁽¹⁾(x) = jₙ(x) + i·yₙ(x)``, vectorizado en ``n``."""
    return spherical_jn(n, x) + _I * spherical_yn(n, x)


def _hankel1_deriv(n: np.ndarray, x: complex) -> np.ndarray:
    """Derivada radial de ``hₙ⁽¹⁾``: ``jₙ′ + i·yₙ′``."""
    return spherical_jn(n, x, derivative=True) + _I * spherical_yn(n, x, derivative=True)


def _check_material(material, nombre: str) -> Material:
    if not isinstance(material, Material):
        raise TypeError(f"{nombre} tiene que ser un Material; llegó {type(material).__name__}")
    return material


def _check_arguments(frequency: float, radius: float) -> None:
    if not np.isfinite(frequency) or frequency <= 0.0:
        raise ValueError(f"frequency tiene que ser > 0; llegó {frequency!r}")
    if not np.isfinite(radius) or radius <= 0.0:
        raise ValueError(f"radius tiene que ser > 0; llegó {radius!r}")


def _auto_n_terms(k: complex, radius: float, n_max_extra: int = 15) -> int:
    """Trunca cuando el término asintótico de la serie ya es despreciable.

    La regla práctica es ``n_max ≈ ka + 4·(ka)^(1/3) + 2``; se agregan términos hasta que el
    coeficiente de mayor orden cae por debajo de ``1e-12`` del máximo, con tope duro para no
    quedarse colgado si ``ka`` es enorme.
    """
    ka = abs(k) * radius
    n_estimado = int(ka + 4.0 * ka ** (1.0 / 3.0) + 2.0)
    return int(min(max(n_estimado, 8) + n_max_extra, 4000))


def _coeficientes(k0: complex, k1: complex, rho0: float, rho1: float, radius: float,
                  n_max: int) -> tuple[np.ndarray, np.ndarray]:
    """``(Aₙ, Bₙ)`` de la derivación del docstring."""
    n = np.arange(n_max + 1, dtype=float)
    x0 = k0 * radius
    x1 = k1 * radius

    r_n = spherical_jn(n, x0)
    rp_n = k0 * spherical_jn(n, x0, derivative=True)
    h_n = _hankel1(n, x0)
    hp_n = k0 * _hankel1_deriv(n, x0)
    j_n = spherical_jn(n, x1)
    jp_n = k1 * spherical_jn(n, x1, derivative=True)

    # G = (ρ₀/ρ₁)·(Jₙ′/Jₙ). El cociente va así y no al revés: con ρ₀/ρ₁, el límite ρ₁→∞ da
    # G → 0 y Aₙ → −Rₙ′/Hₙ′, que es la condición RÍGIDA (∂_r p = 0). Invertido daba el límite
    # BLANDO (p = 0), y como además la presión sigue siendo continua por construcción, el error
    # sólo aparecía en el salto de flujo (medido: 73.7 % constante en ka).
    with np.errstate(divide="ignore", invalid="ignore"):
        g = (rho0 / rho1) * (jp_n / j_n)
        a_n = (g * r_n - rp_n) / (hp_n - g * h_n)
        b_n = (r_n + a_n * h_n) / j_n
    return a_n, b_n


def _as_points(points) -> np.ndarray:
    pts = np.asarray(points, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 3:
        raise ValueError(f"points tiene que ser (M, 3); llegó {pts.shape}")
    return pts


def _evaluar(k: complex, n: np.ndarray, coef: np.ndarray, r: np.ndarray,
             cos_th: np.ndarray, radial: str, derivada: bool = False) -> np.ndarray:
    """Suma ``Σ (2n+1) iⁿ coefₙ · fₙ(kr) · Pₙ(cosθ)`` con ``f`` = ``j`` o ``h⁽¹⁾``.

    ``Pₙ`` se evalúa con la recurrencia de Legendre a mano: `scipy.special.lpn` cambió de firma
    entre versiones y no vale la pena atarse a eso.
    """
    x = k * r
    if radial == "j":
        fn = (spherical_jn(n[:, None], x[None, :], derivative=derivada) * (k if derivada else 1.0))
    elif radial == "h":
        fn = (_hankel1_deriv(n[:, None], x[None, :]) * (k if derivada else 1.0)
              if derivada else _hankel1(n[:, None], x[None, :]))
    else:
        raise ValueError(f"radial desconocida: {radial!r}")

    orden = n.astype(int)
    pols = np.empty((orden.size, cos_th.size), dtype=float)
    pols[0] = 1.0
    if orden.size > 1:
        pols[1] = cos_th
    for i in range(2, orden.size):
        pn = orden[i]
        pols[i] = ((2.0 * pn - 1.0) * cos_th * pols[i - 1] - (pn - 1.0) * pols[i - 2]) / pn

    pesos = (2.0 * n + 1.0) * (_I ** n) * coef
    return pesos @ (fn * pols)


def fluid_sphere_scattering(
    frequency: float,
    medium: Material,
    sphere: Material,
    radius: float,
    points,
    *,
    boundary: str = "fluid",
) -> dict[str, Any]:
    """Campo dispersado y total de una esfera fluida, evaluado en ``points``.

    Parameters
    ----------
    frequency : frecuencia [Hz].
    medium : fluido exterior (de ahí salen ``k₀`` y ``ρ₀``).
    sphere : fluido de la esfera (``k₁``, ``ρ₁``). Con ``boundary='rigid'`` se ignora.
    radius : radio de la esfera [m].
    points : ``(M, 3)`` puntos del exterior (``|x| > radius``).
    boundary : ``'fluid'`` (Anderson, esfera penetrable) o ``'rigid'`` (pared rígida, límite
        ``ρ₁ → ∞``). ``'soft'`` no se expone: la esfera blanda no tiene sentido físico acá.

    Returns
    -------
    dict con ``scattered``, ``incident``, ``total`` (arrays complejos), ``ka``, ``n_terms`` y
    ``coefficients`` (``Aₙ``).
    """
    _check_arguments(frequency, radius)
    _check_material(medium, "medium")
    pts = _as_points(points)
    r = np.linalg.norm(pts, axis=1)
    if np.any(r < radius * (1.0 - 1e-12)):
        raise ValueError(
            "hay puntos dentro de la esfera: esta función evalúa el campo EXTERIOR "
            f"(r mínimo {r.min():.6g} m contra radio {radius:.6g} m)"
        )

    if boundary not in ("fluid", "rigid"):
        raise ValueError(f"boundary desconocida: {boundary!r}; opciones ['fluid', 'rigid']")

    k0 = complex(medium.wavenumber(frequency))
    rho0 = float(medium.rho)
    n_max = _auto_n_terms(k0, radius)
    n = np.arange(n_max + 1, dtype=float)

    if boundary == "rigid":
        # ∂_r p_total = 0 en la superficie ⇒ Aₙ = −Rₙ′/Hₙ′
        x0 = k0 * radius
        rp_n = k0 * spherical_jn(n, x0, derivative=True)
        hp_n = k0 * _hankel1_deriv(n, x0)
        a_n = -rp_n / hp_n
        b_n = np.zeros_like(a_n)
    else:
        _check_material(sphere, "sphere")
        k1 = complex(sphere.wavenumber(frequency))
        a_n, b_n = _coeficientes(k0, k1, rho0, float(sphere.rho), radius, n_max)

    cos_th = pts[:, 2] / np.maximum(r, 1e-300)

    # incidente (coef = 1, jₙ) y dispersado (Aₙ, hₙ⁽¹⁾) en el exterior
    inc = _evaluar(k0, n, np.ones_like(n), r, cos_th, "j")
    sca = _evaluar(k0, n, a_n, r, cos_th, "h")

    # derivadas radiales: hacen falta para verificar la continuidad de flujo en r = a sin
    # diferencias finitas (que a ka grande introducen error propio)
    d_inc = _evaluar(k0, n, np.ones_like(n), r, cos_th, "j", derivada=True)
    d_sca = _evaluar(k0, n, a_n, r, cos_th, "h", derivada=True)

    salida: dict[str, Any] = {
        "scattered": sca,
        "incident": inc,
        "total": inc + sca,
        "d_scattered": d_sca,
        "d_total": d_inc + d_sca,
        "ka": float(abs(k0) * radius),
        "n_terms": int(n_max),
        "coefficients": a_n,
        "boundary": boundary,
        "k_out": k0,
        "rho_out": rho0,
    }
    if boundary == "fluid":
        salida["interior"] = _evaluar(k1, n, b_n, r, cos_th, "j")
        salida["d_interior"] = _evaluar(k1, n, b_n, r, cos_th, "j", derivada=True)
        salida["k_in"] = k1
        salida["rho_in"] = float(sphere.rho)
    return salida
