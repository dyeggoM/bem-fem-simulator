"""Métricas de focalización sobre el campo complejo de presión (Fase 4B).

Todas trabajan sobre arreglos 1D: ``coord`` es la coordenada lateral (``y``, en metros) e
``intensity`` la intensidad ``|p|²`` sobre una línea de evaluación. Se usan para comparar un
transmisor (una apertura) contra un arreglo sintetizado:

* :func:`peak_location`      — dónde está el foco.
* :func:`lateral_fwhm`       — ancho del lóbulo principal a ``−level`` (0.5 ⇒ FWHM) [m].
* :func:`sidelobe_level_db`  — máximo fuera del lóbulo principal relativo al pico [dB].
* :func:`intensity_gain_db`  — ganancia de intensidad de un campo respecto de una referencia [dB]
  (``10·log10`` del cociente de intensidades = ``20·log10`` del cociente de presiones).
* :func:`field_correlation`  — correlación compleja normalizada, invariante a fase global.

Las interpolaciones son lineales (en ``|p|²`` para los cruces de nivel); con grillas finas el sesgo
es despreciable frente a los milímetros que se reportan, y el test lo verifica por refinamiento.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "field_correlation",
    "intensity_gain_db",
    "lateral_fwhm",
    "peak_location",
    "sidelobe_level_db",
]


def _as_line(coord, intensity) -> tuple[np.ndarray, np.ndarray]:
    c = np.asarray(coord, dtype=float).reshape(-1)
    i = np.asarray(intensity, dtype=float).reshape(-1)
    if c.shape != i.shape:
        raise ValueError(f"coord e intensity deben tener igual longitud; {c.shape} vs {i.shape}")
    if c.size < 3:
        raise ValueError("hacen falta al menos 3 muestras")
    if np.any(np.diff(c) <= 0.0):
        raise ValueError("coord debe ser estrictamente creciente")
    if np.any(~np.isfinite(i)):
        raise ValueError("intensity contiene valores no finitos")
    if np.any(i < 0.0):
        raise ValueError("intensity debe ser >= 0 (es |p|²)")
    return c, i


def peak_location(points, field) -> np.ndarray:
    """Punto de ``points`` donde ``|field|²`` es máximo (empate ⇒ el primero).

    ``points``: ``(M, 2)`` [m] (o ``(M,)`` para una línea 1D, y entonces devuelve un escalar).
    """
    pts = np.asarray(points, dtype=float)
    f = np.asarray(field, dtype=complex).reshape(-1)
    if pts.ndim == 2:
        if pts.shape[0] != f.size:
            raise ValueError(f"points {pts.shape} y field {f.shape} no coinciden")
    elif pts.ndim == 1:
        if pts.shape[0] != f.size:
            raise ValueError(f"points {pts.shape} y field {f.shape} no coinciden")
    else:
        raise ValueError(f"points debe ser (M, 2) o (M,); recibido {pts.shape}")
    idx = int(np.argmax(np.abs(f) ** 2))
    return pts[idx].copy()


def lateral_fwhm(coord, intensity, level: float = 0.5) -> float:
    """Ancho a ``level`` del pico por interpolación lineal de los cruces [m].

    ``level = 0.5`` ⇒ FWHM. Busca el primer cruce descendente a cada lado del pico y interpola
    linealmente en ``|p|²``. Falla si el nivel no se cruza dentro de la grilla (lóbulo recortado):
    es preferible un error a un número inventado.
    """
    if not (0.0 < float(level) < 1.0):
        raise ValueError(f"level debe estar en (0, 1); recibido {level!r}")
    c, i = _as_line(coord, intensity)
    peak = float(i.max())
    if peak <= 0.0:
        raise ValueError("intensity es idénticamente nula")
    target = float(level) * peak
    imax = int(np.argmax(i))

    def crossing(indices) -> float:
        prev = imax
        for k in indices:
            if i[k] < target:
                if i[k] == i[prev]:  # pragma: no cover - grilla degenerada
                    return float(c[k])
                frac = (target - i[prev]) / (i[k] - i[prev])
                return float(c[prev] + frac * (c[k] - c[prev]))
            prev = k
        raise ValueError(
            f"el nivel {level:g} del pico no se cruza dentro de la grilla "
            f"[{c[0]:.4g}, {c[-1]:.4g}] m: el lóbulo está recortado"
        )

    left = crossing(range(imax - 1, -1, -1))
    right = crossing(range(imax + 1, c.size))
    return right - left


def _main_lobe_bounds(i: np.ndarray, imax: int) -> tuple[int, int]:
    """Índices (inclusive) que delimitan el lóbulo principal: los mínimos locales que lo flanquean.

    Si no hay mínimo local de un lado (el patrón decae monótonamente hasta el borde), el límite es
    el borde de la grilla.
    """
    left = 0
    for k in range(imax - 1, 0, -1):
        if i[k] <= i[k - 1] and i[k] <= i[k + 1]:
            left = k
            break
    right = i.size - 1
    for k in range(imax + 1, i.size - 1):
        if i[k] <= i[k - 1] and i[k] <= i[k + 1]:
            right = k
            break
    return left, right


def sidelobe_level_db(coord, intensity) -> float:
    """Máximo fuera del lóbulo principal relativo al pico [dB] (negativo o ``-inf``).

    El lóbulo principal se delimita con los mínimos locales (ceros) que flanquean al pico. Si no
    hay nada fuera, devuelve ``-inf``.
    """
    c, i = _as_line(coord, intensity)
    peak = float(i.max())
    if peak <= 0.0:
        raise ValueError("intensity es idénticamente nula")
    imax = int(np.argmax(i))
    left, right = _main_lobe_bounds(i, imax)
    outside = np.concatenate([i[:left], i[right + 1:]])
    if outside.size == 0:
        return -np.inf
    return float(10.0 * np.log10(float(outside.max()) / peak))


def intensity_gain_db(field_ref, field_test) -> float:
    """``10·log10(max|test|² / max|ref|²)`` [dB] (``= 20·log10`` del cociente de presiones).

    Referencia típica: el campo de **un** elemento. Con ``N`` elementos iguales sumando en fase la
    ganancia es ``20·log10(N)`` dB en presión, que es ``10·log10(N²)`` dB de intensidad.
    """
    ref = np.asarray(field_ref, dtype=complex).reshape(-1)
    test = np.asarray(field_test, dtype=complex).reshape(-1)
    p_ref = float(np.max(np.abs(ref) ** 2))
    p_test = float(np.max(np.abs(test) ** 2))
    if p_ref <= 0.0:
        raise ValueError("field_ref es idénticamente nulo")
    if p_test <= 0.0:
        return -np.inf
    return float(10.0 * np.log10(p_test / p_ref))


def field_correlation(a, b) -> float:
    """``|<a,b>| / (‖a‖·‖b‖)`` con ``<a,b> = Σ conj(a)·b``; invariante a fase global.

    ``1`` ⇒ idénticos salvo fase/amplitud global. Falla con campos nulos.
    """
    va = np.asarray(a, dtype=complex).reshape(-1)
    vb = np.asarray(b, dtype=complex).reshape(-1)
    if va.shape != vb.shape:
        raise ValueError(f"a y b deben tener la misma forma; {va.shape} vs {vb.shape}")
    na, nb = float(np.linalg.norm(va)), float(np.linalg.norm(vb))
    if na <= 0.0 or nb <= 0.0:
        raise ValueError("alguno de los campos es nulo")
    return float(abs(np.vdot(va, vb)) / (na * nb))
