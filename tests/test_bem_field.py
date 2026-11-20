"""Contraste del modelo de paneles finitos en espacio libre (BEM de capa simple) contra la
integral de Rayleigh con baffle que usa el dashboard.

Un test por criterio de aceptación, y cada uno imprime el número medido:

1. ``test_panel_angosto_coincide_con_el_cilindro_pulsante`` — panel de ``1e-4·λ`` contra la solución
   exacta del cilindro pulsante 2D (``p = i·ρ·c·v·H₀⁽¹⁾(k·r)/H₀⁽¹⁾′(k·a)``) y su límite de fuente
   lineal; error relativo < 2 %.
2. ``test_campo_cumple_la_ecuacion_de_helmholtz`` — ``|∇²p + k²p| / (|k|²·|p|) < 1 %`` con
   diferencias finitas centradas, y ``test_campo_lejano_decae_como_uno_sobre_raiz_de_r`` —
   ``|p|·√r`` constante al 3 %.
3. ``test_contraste_con_rayleigh_con_baffle`` — arreglo enfocado 250 kHz, ``D = 40 mm``,
   ``F = 90 mm`` en agua: diferencia en dB sobre el eje y en el foco contra
   ``dashboard.field.evaluate``. Se esperaba ~6 dB "por el factor 2 del baffle"; lo que sale es
   0.00 dB y el test documenta por qué (el factor 2 del baffle y el de la lámina bilateral son el
   mismo 2; el precio del baffle está en el semicampo de atrás y en la potencia).
4. ``test_reciprocidad`` y ``test_potencia_radiada_coincide_con_la_fuente_lineal`` — reciprocidad
   de la función de Green al 1 % y balance de energía (potencia en campo lejano contra
   ``ω·ρ·Q²/8``) al 2 %.
5. ``test_validacion_de_entradas`` — toda entrada inválida falla con ``ValueError`` y mensaje.

Los demás tests cubren la API (dict / ``Source`` / ``Panel``), la radiación a ambos lados, la
convergencia de la cuadratura y el costo de la malla 80×80.
"""

from __future__ import annotations

import math
import time

import numpy as np
import pytest
from scipy.special import hankel1

from ultrasound_sim.dashboard import field as campo
from ultrasound_sim.dashboard.bem_field import (
    Panel,
    free_space_field,
    free_space_pressure,
    line_source_pressure,
    panels_from_sources,
)
from ultrasound_sim.materials import WATER, Material

#: Frecuencia de trabajo de todos los tests [Hz] (λ = 6.0 mm en agua).
F = 250e3
WAVELENGTH = WATER.wavelength(F)
K = WATER.wavenumber(F)
OMEGA = 2.0 * math.pi * F

#: Tolerancias de aceptación.
TOL_ANALITICO = 0.02
TOL_HELMHOLTZ = 0.01
TOL_DECAY = 0.03
TOL_RECIPROCIDAD = 0.01
TOL_ENERGIA = 0.02


def relative_error(calculated, reference) -> float:
    """Error relativo L2 entre dos campos."""
    return float(np.linalg.norm(np.ravel(calculated) - np.ravel(reference))
                 / np.linalg.norm(np.ravel(reference)))


def focused_sources(n_elements: int = 64, aperture: float = 40e-3,
                    focus: float = 90e-3) -> list:
    """Fuentes del arreglo enfocado del dashboard (misma ley de fase que la UI)."""
    return campo.array_sources(n_elements, aperture, focus, 0.0, F, WATER)


# --------------------------------------------------------------------------- #
# criterio 1 — límite analítico exacto
# --------------------------------------------------------------------------- #


def test_panel_angosto_coincide_con_el_cilindro_pulsante():
    r"""Un panel de ``1e-4·λ`` coincide con el cilindro pulsante exacto dentro del 2 %.

    La lámina bilateral de ancho ``w`` y velocidad ``v`` equivale a una fuente de línea de caudal
    ``Q = 2·w·v`` [m²/s] (las dos caras). El cilindro de radio ``a`` con velocidad radial
    ``v_c = Q/(2·π·a)`` tiene, con ``p = A·H₀⁽¹⁾(k·r)`` y ``u_r = A·k·H₀⁽¹⁾′/(i·ω·ρ)``,

        ``A = i·ρ·c·v_c / H₀⁽¹⁾′(k·a)``   ⇒   ``p = i·ρ·c·v_c·H₀⁽¹⁾(k·r)/H₀⁽¹⁾′(k·a)``.

    Ojo: la forma ``−ρ·c·v·H₀/H₀′`` que circula en algunos apuntes tiene un ``−i`` de menos —no
    reproduce el límite ``k·a → 0`` ni ``Q = 2·π·a·v_c``—; acá se verifica la forma derivada y,
    además, que a ``k·a = 0.031`` coincide con la fuente lineal (el error residual es
    ``O((k·a)²·ln k·a)``).
    """
    width = 1e-4 * WAVELENGTH
    radius = WAVELENGTH / 200.0
    Q = 2.0 * width
    v_cylinder = Q / (2.0 * math.pi * radius)
    H0_prime = -hankel1(1, K * radius)

    points = np.array([[0.010, 0.000], [0.050, 0.010], [0.300, -0.050]])
    p_panel = free_space_pressure(points, [Panel(0.0, width, 1.0)], F, WATER, n_per_panel=8)

    r = np.hypot(points[:, 0], points[:, 1])
    p_cylinder = 1j * WATER.rho * WATER.c * v_cylinder * hankel1(0, K * r) / H0_prime
    p_line = line_source_pressure(points, Q, F, WATER)

    err_cylinder = relative_error(p_panel, p_cylinder)
    err_line = relative_error(p_cylinder, p_line)

    print(f"\n  panel w={width:.3e} m (1e-4·λ) vs cilindro a={radius:.3e} m (ka={abs(K * radius):.4f}):"
          f" error_rel={err_cylinder:.3e} (criterio {TOL_ANALITICO})"
          f" | cilindro vs fuente lineal: {err_line:.3e}")
    assert err_cylinder < TOL_ANALITICO
    # el cilindro finito ya es el límite de fuente lineal dentro de la misma tolerancia
    assert err_line < TOL_ANALITICO
    assert np.all(np.isfinite(p_panel))


# --------------------------------------------------------------------------- #
# criterio 2 — consistencia del solver
# --------------------------------------------------------------------------- #


def test_campo_cumple_la_ecuacion_de_helmholtz():
    """``∇²p + k²p = 0`` en el interior del dominio, con error relativo < 1 %.

    Laplaciano de 5 puntos con ``h = λ/40``; el error de truncado es ``O((k·h)²) ≈ 0.001``, así que
    lo que se mide es el error del campo, no el del stencil.
    """
    panels = panels_from_sources(focused_sources(16, 20e-3, 60e-3))
    h = WAVELENGTH / 40.0
    puntos = [(0.030, 0.000), (0.045, 0.009), (0.060, 0.000), (0.035, -0.012)]

    peor = 0.0
    medidas = []
    for x0, y0 in puntos:
        xs = np.array([x0, x0 + h, x0 - h])
        ys = np.array([y0, y0 + h, y0 - h])
        p = free_space_field(xs, ys, panels, F, WATER)
        laplaciano = (p[0, 1] + p[0, 2] + p[1, 0] + p[2, 0] - 4.0 * p[0, 0]) / h**2
        residual = abs(laplaciano + K**2 * p[0, 0]) / (abs(K) ** 2 * abs(p[0, 0]))
        medidas.append(residual)
        peor = max(peor, residual)

    print(f"\n  residual |∇²p+k²p|/(|k|²|p|) con h=λ/40 en {len(puntos)} puntos: "
          f"{', '.join(f'{v:.2e}' for v in medidas)} (peor {peor:.2e}, criterio {TOL_HELMHOLTZ})")
    assert peor < TOL_HELMHOLTZ


def test_campo_lejano_decae_como_uno_sobre_raiz_de_r():
    """En campo lejano ``|p|·√r`` es constante al 3 % (fuente lineal dominante).

    Panel de 1 mm (0.17·λ) sin fase: su caudal neto no es nulo, así que lejos domina el término
    monopolar ``∝ H₀⁽¹⁾(k·r) ∝ r^(−1/2)``.
    """
    radios = np.linspace(0.3, 3.0, 12)
    p = free_space_field(radios, np.zeros_like(radios), [Panel(0.0, 1e-3, 1.0)], F, WATER)
    amplitud = np.abs(p[0]) * np.sqrt(radios)
    dispersion = float((amplitud.max() - amplitud.min()) / amplitud.mean())

    print(f"\n  |p|·√r entre r=0.3 m y r=3 m: min={amplitud.min():.6e} max={amplitud.max():.6e}"
          f" dispersión={dispersion:.2e} (criterio {TOL_DECAY})")
    assert dispersion < TOL_DECAY


# --------------------------------------------------------------------------- #
# criterio 3 — contraste con la integral de Rayleigh con baffle
# --------------------------------------------------------------------------- #


def test_contraste_con_rayleigh_con_baffle():
    """Diferencia en dB, sobre el eje y en el foco, entre espacio libre y Rayleigh con baffle.

    Arreglo enfocado realista: 250 kHz, ``D = 40 mm``, ``F = 90 mm``, agua, 64 elementos. Se
    esperaba ~6 dB por el factor 2 del baffle; **lo medido es 0.00 dB** y el test lo documenta: el
    factor 2 de la imagen del baffle y el factor 2 de la lámina que radia por las dos caras son el
    mismo 2, así que las dos integrales son idénticas (fase opuesta). Los 6 dB aparecen sólo si se
    compara contra una densidad unilateral, que no satisface el salto de velocidad.
    """
    aperture, focus = 40e-3, 90e-3
    sources = focused_sources(64, aperture, focus)
    xs = np.linspace(1e-3, 0.150, 80)
    ys = np.array([0.0])

    p_free = free_space_field(xs, ys, sources, F, WATER)
    p_baffle = campo.evaluate(xs, ys, sources, F, WATER)

    # el foco exacto (no tiene por qué ser un nodo de la malla)
    p_free_focus = free_space_field(np.array([focus]), ys, sources, F, WATER)[0, 0]
    p_baffle_focus = campo.evaluate(np.array([focus]), ys, sources, F, WATER)[0, 0]

    db_focus = 20.0 * math.log10(abs(p_free_focus) / abs(p_baffle_focus))
    db_eje = 20.0 * np.log10(np.abs(p_free[0]) / np.abs(p_baffle[0]))
    fase = np.angle(p_free_focus / p_baffle_focus)
    # lo que daría comparar contra una densidad unilateral (mitad de σ): el "6 dB" esperado
    db_unilateral = 20.0 * math.log10(abs(p_free_focus) / (0.5 * abs(p_baffle_focus)))

    print(f"\n  foco (x={focus * 1e3:.0f} mm, y=0): |p_free|={abs(p_free_focus):.6e} Pa"
          f" |p_baffle|={abs(p_baffle_focus):.6e} Pa"
          f" -> ΔdB={db_focus:+.4f} (criterio ±0.01)")
    print(f"  eje completo ({xs.size} puntos): |ΔdB| máximo={np.abs(db_eje).max():.2e};"
          f" desfase free−baffle={fase:+.6f} rad (|π|={math.pi:.6f})")
    print(f"  lectura con densidad unilateral (no es este modelo): {db_unilateral:+.4f} dB")

    assert abs(db_focus) < 0.01, "el modelo libre y el de baffle deben coincidir en módulo"
    assert np.abs(db_eje).max() < 0.01
    # fase opuesta: el convenio global de aperture.py (p = −(ωρ/2)vH₀) es −p_free
    assert abs(abs(fase) - math.pi) < 1e-6
    # los 6 dB esperados son 20·log10(2); quedan como referencia de la lectura unilateral
    assert abs(db_unilateral - 20.0 * math.log10(2.0)) < 0.01


def test_radiacion_a_ambos_lados_y_simetria():
    """La lámina radia también hacia ``x < 0`` (con ``p(−x) = p(x)``), cosa que el baffle no hace.

    Es la diferencia física que el dashboard no puede representar: ``dashboard.field.evaluate``
    rechaza ``x ≤ 0``.
    """
    sources = focused_sources(16, 20e-3, 60e-3)
    xs = np.array([2e-3, 10e-3, 30e-3])
    p_delante = free_space_field(xs, np.zeros(3), sources, F, WATER)
    p_detras = free_space_field(-xs, np.zeros(3), sources, F, WATER)
    diferencia = float(np.abs(p_delante - p_detras).max() / np.abs(p_delante).max())

    print(f"\n  |p(−x) − p(x)| / |p|máx = {diferencia:.2e} (radiación simétrica a los dos lados);"
          f" el dashboard rechaza x ≤ 0")
    assert diferencia < 1e-12

    with pytest.raises(ValueError, match="baffle"):
        campo.evaluate(-xs, np.zeros(3), sources, F, WATER)


# --------------------------------------------------------------------------- #
# criterio 4 — reciprocidad y balance de energía
# --------------------------------------------------------------------------- #


def test_reciprocidad():
    """Intercambiar fuente y receptor da el mismo valor (función de Green simétrica) al 1 %.

    Se usan paneles puntuales (``1e-4·λ``) sobre la lámina: el campo en ``B`` debido a una fuente en
    ``A`` debe ser igual al campo en ``A`` debido a una fuente en ``B``.
    """
    width = 1e-4 * WAVELENGTH
    pares = [((0.0, -5e-3), (0.0, 12e-3)), ((0.0, -18e-3), (0.0, 18e-3)),
             ((0.0, 3e-3), (0.0, 40e-3))]
    errores = []
    for (ya, yb) in [(a[1], b[1]) for a, b in pares]:
        p_ab = free_space_pressure([[0.0, yb]], [Panel(ya, width, 1.0)], F, WATER)[0]
        p_ba = free_space_pressure([[0.0, ya]], [Panel(yb, width, 1.0)], F, WATER)[0]
        errores.append(abs(p_ab - p_ba) / abs(p_ab))

    print(f"\n  reciprocidad en {len(errores)} pares: "
          f"{', '.join(f'{e:.2e}' for e in errores)} (criterio {TOL_RECIPROCIDAD})")
    assert max(errores) < TOL_RECIPROCIDAD


def test_potencia_radiada_coincide_con_la_fuente_lineal():
    """Balance de energía: la potencia en campo lejano es ``ω·ρ·Q²/8`` [W/m] al 2 %.

    Integrando ``|p|²/(2·ρ·c)`` sobre un círculo de radio ``R = 2 m`` (campo lejano). El resultado
    depende del **cuadrado** de la densidad de capa simple, así que este test es el que detecta un
    factor 2 perdido en ``σ = −2·i·ω·ρ·v``: con la mitad de densidad la potencia caería 4×.
    """
    width = 1e-3 * WAVELENGTH
    radius = 2.0
    n = 2000
    theta = (np.arange(n) + 0.5) * 2.0 * np.pi / n
    points = np.column_stack([radius * np.cos(theta), radius * np.sin(theta)])

    p = free_space_pressure(points, [Panel(0.0, width, 1.0)], F, WATER)
    potencia = float(np.sum(np.abs(p) ** 2 / (2.0 * WATER.rho * WATER.c))
                     * (2.0 * np.pi * radius / n))
    Q = 2.0 * width
    esperada = OMEGA * WATER.rho * Q**2 / 8.0
    error = abs(potencia - esperada) / esperada

    print(f"\n  potencia en r={radius} m: {potencia:.6e} W/m vs ωρQ²/8 = {esperada:.6e} W/m"
          f" -> error_rel={error:.2e} (criterio {TOL_ENERGIA})")
    assert error < TOL_ENERGIA


def test_medio_con_atenuacion_coincide_con_la_fuente_lineal():
    """Con ``k`` complejo (α > 0) el panel angosto sigue coincidiendo con la fuente lineal.

    La atenuación entra por ``k = ω/c + i·α(f)`` en ``H₀⁽¹⁾``, sin ningún término extra: es el mismo
    camino que usa :func:`line_source_pressure`, así que la comparación mide que no se perdió la
    parte imaginaria. ``α = 100 Np/m`` a 250 kHz es un parámetro de test (≈0.43 dB/cm), no un valor
    de tejido: ``SOFT_TISSUE`` está marcado VERIFY en ``materials.py``.
    """
    lossy = Material(name="agua_con_pérdidas_test", c=WATER.c, rho=WATER.rho,
                     alpha0=400.0, alpha_power=1.0, f_ref=1.0e6,
                     provenance="parámetro de test, no es un dato físico")
    k = lossy.wavenumber(F)
    width = 1e-4 * WAVELENGTH
    points = np.array([[0.010, 0.000], [0.050, 0.010], [0.300, -0.050]])

    p_panel = free_space_pressure(points, [Panel(0.0, width, 1.0)], F, lossy, n_per_panel=8)
    p_line = line_source_pressure(points, 2.0 * width, F, lossy)
    error = relative_error(p_panel, p_line)

    print(f"\n  k={k:.3f} 1/m (Im/Re={k.imag / k.real:.1%}): panel vs fuente lineal "
          f"error_rel={error:.2e}")
    assert k.imag > 1.0, "el test no estaría ejercitando k complejo"
    assert error < 1e-3


# --------------------------------------------------------------------------- #
# criterio 5 — validación de entradas
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "cambio, fragmento",
    [
        ({"frequency": 0.0}, "frequency debe ser > 0"),
        ({"frequency": -1.0}, "frequency debe ser > 0"),
        ({"frequency": float("nan")}, "frequency debe ser finita"),
        ({"panels": []}, "panels está vacío"),
        ({"panels": None}, "panels es None"),
        ({"panels": [Panel(0.0, 0.0, 1.0)]}, "width debe ser > 0"),
        ({"panels": [Panel(0.0, 1e-3, float("nan"))]}, "'amp' debe ser finito"),
        ({"panels": [Panel(float("inf"), 1e-3, 1.0)]}, "'y' debe ser finito"),
        ({"panels": [{"width": 1e-3}]}, "clave 'y'"),
        ({"panels": ["panel"]}, "cada panel debe ser"),
        ({"n_per_panel": 1}, "n_per_panel debe ser >= 2"),
        ({"xs": [[0.05]]}, "xs debe ser un array 1-D"),
        ({"xs": []}, "xs está vacío"),
        ({"xs": [float("nan")]}, "xs contiene NaN o inf"),
        ({"ys": [float("inf")]}, "ys contiene NaN o inf"),
        ({"material": "agua"}, "material debe ser un materials.Material"),
        ({"equation": "burton_miller"}, "equation debe ser"),
        ({"equation": "chief"}, "equation debe ser"),
        # punto sobre la lámina: singularidad logarítmica de H₀
        ({"xs": [0.0]}, "sobre la lámina"),
    ],
)
def test_validacion_de_entradas(cambio, fragmento):
    """Toda entrada inválida falla con ``ValueError`` y un mensaje que dice qué está mal."""
    base: dict = dict(xs=[0.05], ys=[0.0], panels=[Panel(0.0, 1e-3, 1.0)],
                      frequency=F, material=WATER)
    base.update(cambio)
    with pytest.raises(ValueError, match=fragmento):
        free_space_field(**base)


def test_validacion_de_puntos_y_salida_finita():
    """``free_space_pressure`` valida ``points`` y el resultado es finito y de la forma pedida."""
    with pytest.raises(ValueError, match="points debe tener forma"):
        free_space_pressure([[0.05, 0.0, 0.0]], [Panel(0.0, 1e-3)], F, WATER)
    with pytest.raises(ValueError, match="points está vacío"):
        free_space_pressure(np.zeros((0, 2)), [Panel(0.0, 1e-3)], F, WATER)
    with pytest.raises(ValueError, match="points contiene NaN"):
        free_space_pressure([[np.nan, 0.0]], [Panel(0.0, 1e-3)], F, WATER)

    out = free_space_field([0.01, 0.02], [0.0, 0.01, -0.01], [Panel(0.0, 1e-3, 1.0)], F, WATER)
    print(f"\n  salida: dtype={out.dtype} forma={out.shape} finito={np.all(np.isfinite(out))}")
    assert out.shape == (3, 2)
    assert out.dtype == complex
    assert np.all(np.isfinite(out))
    assert np.all(np.abs(out) > 0.0)


# --------------------------------------------------------------------------- #
# API, cuadratura y costo
# --------------------------------------------------------------------------- #


def test_panel_dict_source_y_panel_dan_exactamente_lo_mismo():
    """La convención de datos del dashboard (``dict`` o ``Source``) da el mismo resultado que
    :class:`Panel`."""
    punto = [0.03, 0.01]
    p_source = free_space_pressure([punto], [campo.Source(y=0.002, width=1e-3, amp=0.7, phase=0.4)],
                                   F, WATER)[0]
    p_dict = free_space_pressure([punto], [{"y": 0.002, "width": 1e-3, "amp": 0.7, "phase": 0.4}],
                                 F, WATER)[0]
    p_panel = free_space_pressure([punto], [Panel(0.002, 1e-3, 0.7, 0.4)], F, WATER)[0]

    print(f"\n  Source={p_source:.6e} dict={p_dict:.6e} Panel={p_panel:.6e}")
    assert p_source == p_dict == p_panel
    # también se acepta un único panel en forma de dict (sin lista)
    assert free_space_pressure([punto], {"y": 0.002, "width": 1e-3, "amp": 0.7, "phase": 0.4},
                               F, WATER)[0] == p_panel


def test_cuadratura_converge_con_n_per_panel():
    """Con ``n_per_panel = 6`` ya está en el redondeo para este integrando suave.

    Documenta que el valor por defecto (16) es holgado y que el error del campo, cuando la malla no
    toca la lámina, no depende de la cuadratura.
    """
    panels = panels_from_sources(focused_sources(16, 20e-3, 60e-3))
    puntos = np.array([[0.03, 0.0], [0.045, 0.009], [0.06, -0.01]])
    fino = free_space_pressure(puntos, panels, F, WATER, n_per_panel=48)
    diferencias = {n: float(np.abs(free_space_pressure(puntos, panels, F, WATER, n_per_panel=n)
                                   - fino).max() / np.abs(fino).max())
                   for n in (6, 8, 16)}
    print(f"\n  |p(n) − p(48)|/|p|máx: " + ", ".join(f"n={n}: {d:.2e}" for n, d in diferencias.items()))
    assert max(diferencias.values()) < 1e-10


def test_costo_malla_80x80():
    """Costo del cálculo para una malla 80×80 con el arreglo de 64 elementos (se reporta el tiempo).

    No se compara contra un umbral ajustado (depende de la máquina); sólo se exige que termine y que
    el resultado sea finito.
    """
    sources = focused_sources(64, 40e-3, 90e-3)
    xs = np.linspace(0.5e-3, 40e-3, 80)
    ys = np.linspace(-20e-3, 20e-3, 80)

    inicio = time.perf_counter()
    p = free_space_field(xs, ys, sources, F, WATER)
    segundos = time.perf_counter() - inicio

    print(f"\n  malla 80×80, {len(sources)} paneles, n_per_panel=16: {segundos:.2f} s"
          f" ({len(sources) * 16} nodos por punto), forma={p.shape}, "
          f"finito={np.all(np.isfinite(p))}")
    assert p.shape == (80, 80)
    assert np.all(np.isfinite(p))
    assert segundos < 120.0
