"""BEM de contorno para Helmholtz 2D en dominio **exterior** con condición de Neumann.

Convención del repo (``docs/decisions.md`` D4): tiempo ``exp(−i·ω·t)``, ``∇²p + k²p = 0``,
``k = ω/c + i·α`` con ``Im(k) ≥ 0``, ondas salientes ``∝ exp(+i·k·r)`` y

    G(x, y) = (i/4) · H₀⁽¹⁾(k·|x − y|)          (Green 2D saliente)

Formulación (directa, incógnita = presión total sobre la frontera)
------------------------------------------------------------------
Sea ``Ω`` el fluido exterior y ``S`` el contorno del cuerpo. La normal ``n̂`` de la malla es la
**saliente del cuerpo**, que apunta hacia el fluido (``n̂ = r̂`` en un círculo; la malla es CCW y
``n̂ = (t_y, −t_x)``, coherente con ``meshing.outward_normals``). Con ``p`` total, ``p_inc`` la onda
plana incidente y ``q = ∂p/∂n`` el dato de Neumann, la representación integral del campo dispersado
``p_s = p − p_inc`` (que satisface la condición de radiación) llevada al contorno da

    C(x)·p(x) = p_inc(x) + ∫_S [ p(y)·∂G/∂n_y − G(x, y)·q(y) ] ds_y ,      C = 1/2 en S lisa

El término ``p_inc`` es el "término libre" del campo incidente: aplicando el teorema de extinción a
``p_inc`` (que es entero, satisface Helmholtz dentro del cuerpo) su integral de frontera se anula
fuera del cuerpo, y por eso el lado derecho de la ecuación de la presión **total** es simplemente
``p_inc − ∫G·q``. Colocando ``x`` en el punto medio de cada panel:

    (½·I − D)·p = p_inc − S·q,        S_ij = ∫_{Γ_j} G ds,   D_ij = ∫_{Γ_j} ∂G/∂n_y ds

que es un sistema de segunda especie (bien condicionado, ver más abajo). Con pared rígida
``q = 0`` y el sistema se reduce a ``(½·I − D)·p = p_inc``.

El campo dispersado en un punto del fluido se recupera con la misma representación:

    p_s(x) = ∫_S [ p(y)·∂G/∂n_y − G(x, y)·q(y) ] ds_y ,   x ∉ S

y ``p_total = p_inc + p_s``. (Equivalentemente, resolviendo para ``p_s`` con el dato ``q_s = q − q_inc``
se obtiene el mismo resultado; el test ``test_scattered_formulation_matches_total_formulation`` lo
usa como verificación cruzada de las matrices de capa simple.)

Elementos y cuadratura
----------------------
Elementos **constantes** por panel (presión y ``q`` constantes en cada segmento), colocación en el
punto medio. Se elige constante y no lineal por dos razones concretas:

1. el término diagonal de ``D`` es exactamente nulo y el de ``S`` se resuelve por extracción exacta
   de la singularidad logarítmica (ver abajo), así que no queda ninguna integral singular tratada a
   ciegas;
2. para contornos suaves la convergencia observada es O(h²) en el campo (medidas abajo), más que
   suficiente para el 2 % de la Fase 3. Los elementos lineales duplicarían las incógnitas sin que
   este problema los necesite; quedan como mejora si alguna vez hace falta más precisión.

Integrales singulares:

* ``S_ii`` (panel que contiene al punto de colocación): es la única integral con una singularidad
  logarítmica real, y se resuelve por **extracción de la singularidad**, no por cuadratura a
  ciegas. Con ``h = L/2``, ``z = k·h`` y ``γ`` la constante de Euler::

      S_ii = (i/4)·2·∫₀^h H₀⁽¹⁾(k·r) dr
           = (i/4)·(2/k)·∫₀^z H₀⁽¹⁾(t) dt

  y dentro de la última integral se separa la parte logarítmica con la asíntota de ``Y₀``::

      H₀⁽¹⁾(t) = 1 + (2i/π)(ln(t/2) + γ) + R(t),      R(t) = O(t²·ln t) suave

  ``∫₀^z [1 + (2i/π)(ln(t/2)+γ)] dt = z + (2i/π)·z·(ln(z/2) + γ − 1)`` es cerrada, y ``R`` se
  integra con Gauss-Legendre sobre ``t = z·τ²`` (la sustitución cuadrática aplana el ``t²·ln t``
  del origen; con 12 nodos ya converge al redondeo). **Ojo**: ``∫₀^z t·H₀⁽¹⁾(t) dt = z·H₁⁽¹⁾(z)`` es
  cierta, pero para ``∫₀^z H₀⁽¹⁾ dt`` no hay forma cerrada elemental (involucra funciones de
  Struve); usar esa identidad acá da un ``S_ii`` con el orden de magnitud equivocado.
* ``D_ii``: para un panel **recto** el núcleo de doble capa es idénticamente nulo sobre él, porque
  ``∂G/∂n_y ∝ (y − x)·n̂`` y todos los ``y`` del panel cumplen ``(y − x_i)·n̂ = 0`` cuando ``x_i`` es
  su punto medio. El valor principal de la integral es entonces 0 y el salto ``1/2`` ya está
  explícito en la ecuación ``(½·I − D)``: el elemento diagonal es ``0`` sin aproximación. (Es el
  mismo hecho que en el caso clásico "la integral del núcleo de doble capa sobre el elemento que
  contiene al punto vale ``∓½``": acá el ``½`` es el término de salto y la parte de valor principal
  sobre el panel plano es nula. Si el panel estuviera curvado, habría que integrar el resto.)
* Pares **cercanos** (``dist(x_i, Γ_j) < L_j/2``): la cuadratura de Gauss-Legendre se aplica sobre
  una subdivisión geométrica del panel (razón 2) que se refina hacia el punto de ``Γ_j`` más
  próximo a ``x_i``, hasta que cada subsegmento es del tamaño de esa distancia. Sin esto la
  cuadratura de Gauss fija pierde varios dígitos en los paneles vecinos.
* Pares lejanos: Gauss-Legendre de orden ``quad_order`` (16 por defecto) sobre el panel completo,
  vectorizado sobre todos los pares.

Condicionamiento y discretización
---------------------------------
El sistema es de segunda especie, así que ``cond(½·I − D)`` se mantiene bajo y crece despacio con
``ka``. Medido con ``panels_for_frequency`` (16 paneles por longitud de onda, mínimo 64):

===========  =========  ===========  ==================
``ka``       paneles    paneles/λ    ``cond(A)`` (medido)
===========  =========  ===========  ==================
1            64         64           2.4
2            64         32           2.5
5            128        25.6         7.4
===========  =========  ===========  ==================

La discretización usa ``N = 2^m`` paneles uniformes (``≥ 16`` por longitud de onda, ``≥ 64`` en
total). Con esos números el error relativo L2 del campo dispersado contra la serie exacta es
``1.3e-3`` (ka=1), ``2.7e-3`` (ka=2) y ``2.4e-3`` (ka=5), y baja como O(h²) al duplicar ``N``.

Elección de ecuación (``equation``)
-----------------------------------
La ecuación de arriba, ``(½I − D)·p = p_inc − S·q``, es de segunda especie y por eso tiene
**resonancias espurias** en las frecuencias propias del problema interior de Dirichlet
(``J_n(ka) = 0``): en esos ``ka`` el operador ``½I − D`` es (casi) singular y la solución explota.
Se ofrecen tres ecuaciones; las tres resuelven la misma incógnita ``p``, cambia sólo la ecuación.

* ``"standard"``: la de arriba. Barata y exacta fuera de las resonancias, inservible en ellas.
* ``"burton_miller"``: combinación lineal de la ecuación de doble capa y su derivada normal. Con
  la incógnita dispersada ``u = p_s`` (``q_s = q − q_inc``)::

      [(½I − D) + η·W]·u = −[S + η·(½I − K')]·q_s

  y con la incógnita total ``p`` se suma el mismo operador aplicado a ``p_inc`` (que es el
  ``p_inc − S·q_inc + η·W·p_inc`` de la extinción del campo incidente, ver más abajo). ``η = i/k``
  (complejo, de módulo ``1/k``). ``W`` es el operador hipersingular (derivada
  normal de la doble capa) y ``K'`` la doble capa **adjunta** (derivada normal en el punto de
  colocación). Como la ecuación de doble capa falla en los ``ka`` con ``J_n = 0`` y la
  hipersingular falla en los ``ka`` con ``J_n' = 0``, y los dos conjuntos de autovalores del
  interior (Dirichlet y Neumann) están intercalados y nunca coinciden, la combinación con
  ``Im η ≠ 0`` es inyectiva en **todas** las frecuencias. El sistema es cuadrado (``N×N``) y algo
  peor condicionado que el estándar (entradas ``~η/h``), pero resoluble directo.
* ``"chief"``: ecuación estándar más ``M`` filas extra imponiendo que la representación integral
  del campo dispersado se anule en puntos **interiores** (teorema del campo nulo). Por defecto
  ``M = 4`` puntos sobre la circunferencia ``r = a/2``, uniformes en ángulo. El sistema es
  sobredeterminado (``N+M`` filas, ``N`` incógnitas) y se resuelve por mínimos cuadrados
  (``numpy.linalg.lstsq``); ``residual`` reporta el residuo relativo de esa solución.

Los dos términos nuevos se discretizan por colocación en la misma malla de paneles constantes.
``W`` es hipersingular (``O(1/r²)``) y se regulariza por **sustracción de la singularidad
estática**: se resta el núcleo de Laplace en forma cerrada y se reemplaza su integral sobre el
panel propio por el valor principal finito de Hadamard (``2/(π·L)``, con ``L`` la longitud del
panel), integrando el resto con Gauss-Legendre; en el resto se resta analíticamente el término
``1/z`` de ``H₁`` para que la cancelación no sea numérica (ver ``_hypersingular_kernels``). El
núcleo de ``K'`` es ``O(1/r)`` y su elemento diagonal es exactamente nulo por el mismo argumento
que ``D_ii``. ``K'`` incluye sólo el valor principal: el salto ``½`` está explícito en la ecuación.

Identidades usadas para escribir el miembro derecho de Burton-Miller con la incógnita total
(verificadas numéricamente contra campos radiantes fabricados, convergentes como el error de
discretización):

    (½I − D)·p_inc = p_inc − S·q_inc            (extinción del campo incidente)
    W·u = (−½I + K')·q,   (D − S)·u = 0 en el interior   (u radiante)

``ka = 1, 2, 5`` no coinciden con ninguna resonancia: ahí las tres ecuaciones coinciden dentro
del error de discretización (test ``tests/test_bem_formulations.py``).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import hankel1

from ..materials import Material

__all__ = [
    "EQUATIONS",
    "CircleBoundary",
    "ExteriorSolution",
    "assemble_operators",
    "circle_boundary",
    "panels_for_frequency",
    "solve_exterior_neumann",
]

#: Paneles por longitud de onda usados por ``panels_for_frequency``.
PANELS_PER_WAVELENGTH = 16.0

#: Mínimo de paneles (los ``ka`` chicos no necesitan más, pero el costo es despreciable).
MIN_PANELS = 64

#: Orden de la cuadratura de Gauss-Legendre por subsegmento.
DEFAULT_QUAD_ORDER = 16

#: Máxima subdivisión geométrica de un panel para los pares cercanos.
DEFAULT_MAX_SUBDIVISIONS = 12

#: Ecuaciones integrales disponibles (ver "Elección de ecuación" en el docstring del módulo).
EQUATIONS = ("standard", "burton_miller", "chief")

#: Puntos interiores de CHIEF por defecto.
DEFAULT_CHIEF_POINTS = 4

#: Radio (relativo a ``radius``) de la circunferencia donde se ubican los puntos de CHIEF.
CHIEF_RADIUS_FRACTION = 0.5


# --------------------------------------------------------------------------- #
# malla de contorno
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class CircleBoundary:
    """Contorno circular discretizado en ``n_panels`` segmentos rectos (CCW).

    Atributos
    ---------
    radius : radio [m].
    n_panels : cantidad de paneles.
    p0, p1 : ``(N, 2)`` extremos de cada panel; ``p0 → p1`` recorre el círculo en sentido CCW.
    midpoints : ``(N, 2)`` puntos de colocación (punto medio de cada panel).
    normals : ``(N, 2)`` normales unitarias **salientes** del cuerpo (apuntan al fluido).
    lengths : ``(N,)`` longitud de cada panel [m].
    """

    radius: float
    n_panels: int
    p0: np.ndarray
    p1: np.ndarray
    midpoints: np.ndarray
    normals: np.ndarray
    lengths: np.ndarray

    @property
    def perimeter(self) -> float:
        return float(np.sum(self.lengths))

    @property
    def panel_length(self) -> float:
        return float(self.lengths[0])


def circle_boundary(radius: float, n_panels: int) -> CircleBoundary:
    """Discretiza el círculo ``r = radius`` en ``n_panels`` segmentos iguales orientados CCW.

    Los vértices se ubican en ``θ_j = 2π·j/N`` y cada panel es la cuerda ``[θ_j, θ_{j+1}]``. La
    normal saliente se calcula como ``(t_y, −t_x)`` con ``t = p1 − p0`` (misma convención que
    ``meshing.outward_normals`` para mallas CCW) y se verifica que apunte hacia afuera.
    """
    if not (radius > 0.0):
        raise ValueError(f"radius debe ser > 0; recibido {radius!r}")
    if n_panels < 8:
        raise ValueError(f"n_panels debe ser >= 8; recibido {n_panels!r}")

    theta = 2.0 * np.pi * np.arange(n_panels) / n_panels
    vertices = radius * np.column_stack([np.cos(theta), np.sin(theta)])

    p0 = vertices
    p1 = np.roll(vertices, -1, axis=0)
    midpoints = 0.5 * (p0 + p1)
    tangent = p1 - p0
    lengths = np.hypot(tangent[:, 0], tangent[:, 1])
    normals = np.column_stack([tangent[:, 1], -tangent[:, 0]]) / lengths[:, None]

    if not np.all(np.sum(midpoints * normals, axis=1) > 0.0):  # pragma: no cover - salvaguarda
        raise RuntimeError("las normales calculadas no son salientes; revisar la orientación")

    return CircleBoundary(
        radius=float(radius),
        n_panels=int(n_panels),
        p0=p0,
        p1=p1,
        midpoints=midpoints,
        normals=normals,
        lengths=lengths,
    )


def panels_for_frequency(radius: float, wavenumber: complex,
                         panels_per_wavelength: float = PANELS_PER_WAVELENGTH,
                         minimum: int = MIN_PANELS) -> int:
    """Cantidad de paneles ``2^m`` recomendada para un ``ka`` dado.

    La circunferencia mide ``2π·a`` y hay ``ka`` longitudes de onda en ella, así que con
    ``panels_per_wavelength`` paneles por longitud de onda se necesitan ``≈ ppw·ka`` paneles.
    Se redondea a la potencia de dos siguiente y se respeta el mínimo ``minimum``.
    """
    if not (panels_per_wavelength > 0.0):
        raise ValueError("panels_per_wavelength debe ser > 0")
    if minimum < 8:
        raise ValueError("minimum debe ser >= 8")
    ka = abs(complex(wavenumber)) * radius
    needed = max(float(minimum), panels_per_wavelength * ka)
    return int(2 ** int(np.ceil(np.log2(needed))))


# --------------------------------------------------------------------------- #
# cuadratura sobre un panel
# --------------------------------------------------------------------------- #


def _segment_quadrature(p0: np.ndarray, p1: np.ndarray, x: np.ndarray,
                        order: int = DEFAULT_QUAD_ORDER,
                        max_subdivisions: int = DEFAULT_MAX_SUBDIVISIONS):
    """Nodos y pesos de Gauss-Legendre sobre el segmento ``p0-p1``, graduados hacia ``x``.

    El panel se parte en subsegmentos de tamaño geométricamente creciente (razón 2) a partir del
    punto del panel más cercano a ``x``; la cantidad de subsegmentos se elige para que el más
    pequeño sea del orden de ``dist(x, panel)``. Así el núcleo logarítmico (``G``) o casi singular
    (``∂G/∂n``) queda resuelto por la cuadratura incluso cuando ``x`` está sobre el panel vecino.

    Devuelve ``(t, w)`` con ``t ∈ [0, 1]`` parámetros sobre el panel y ``w`` los pesos en metros
    (``ds = L·dt`` ya incluido).
    """
    d = p1 - p0
    length = float(np.hypot(d[0], d[1]))
    t_star = float(np.clip(np.dot(x - p0, d) / (length * length), 0.0, 1.0))
    dist = float(np.hypot(*(p0 + t_star * d - x)))

    nodes, weights = leggauss(order)
    unit_nodes = 0.5 * (nodes + 1.0)
    unit_weights = 0.5 * weights

    params, wts = [], []
    for direction, end in ((1.0, 1.0), (-1.0, 0.0)):
        span = abs(end - t_star)
        if span <= 0.0:
            continue
        if dist >= span:
            n_sub = 1
        else:
            n_sub = int(np.clip(np.ceil(np.log2(span / max(dist, 1e-15))), 1, max_subdivisions))
        base = span / (2.0 ** n_sub - 1.0)
        cuts = [t_star]
        acc = 0.0
        for k in range(n_sub):
            acc += base * 2.0 ** k
            cuts.append(t_star + direction * acc)
        cuts[-1] = end
        cuts = np.sort(np.asarray(cuts))
        for ta, tb in zip(cuts[:-1], cuts[1:]):
            params.append(ta + unit_nodes * (tb - ta))
            wts.append(unit_weights * (tb - ta) * length)

    return np.concatenate(params), np.concatenate(wts)


def _distance_to_panels(x: np.ndarray, boundary: CircleBoundary) -> np.ndarray:
    """Distancia de cada punto de ``x`` (M) a cada panel (N); forma ``(M, N)``."""
    d = (boundary.p1 - boundary.p0)[None, :, :]                     # (1,N,2)
    rel = x[:, None, :] - boundary.p0[None, :, :]                   # (M,N,2)
    length2 = np.sum(d * d, axis=2, keepdims=True)                  # (1,N,1)
    t = np.clip(np.sum(rel * d, axis=2, keepdims=True) / length2, 0.0, 1.0)
    closest = boundary.p0[None, :, :] + t * d
    return np.hypot(*(x[:, None, :] - closest).transpose(2, 0, 1))


def _self_panel_single_layer(k: complex, length: float, order: int = 24) -> complex:
    """``∫_{Γ_i} G ds`` sobre el panel propio, por extracción de la singularidad logarítmica.

    ``S_ii = (i/4)·2·∫₀^h H₀⁽¹⁾(k·r) dr`` con ``h = length/2``. Separando

        H₀⁽¹⁾(t) = 1 + (2i/π)(ln(t/2) + γ) + R(t)

    la parte logarítmica integra en forma cerrada y ``R = O(t² ln t)`` se resuelve con
    Gauss-Legendre sobre ``t = k·h·τ²`` (ver docstring del módulo: **no** existe forma cerrada
    elemental para ``∫₀^z H₀ dt``).
    """
    h = 0.5 * length
    z = k * h
    tau, w = leggauss(order)
    tau = 0.5 * (tau + 1.0)
    w = 0.5 * w
    t = z * tau * tau
    jac = 2.0 * z * tau
    remainder = hankel1(0, t) - (1.0 + (2j / np.pi) * (np.log(t / 2.0) + np.euler_gamma))
    smooth = np.sum(w * remainder * jac)
    log_part = z + (2j / np.pi) * z * (np.log(z / 2.0) + np.euler_gamma - 1.0)
    return 0.25j * (2.0 / k) * (smooth + log_part)


def _layer_matrices(points: np.ndarray, boundary: CircleBoundary, k: complex,
                    order: int = DEFAULT_QUAD_ORDER,
                    max_subdivisions: int = DEFAULT_MAX_SUBDIVISIONS,
                    handle_self_panel: bool = False):
    """Matrices de capa simple ``S`` y doble ``D``, forma ``(M, N)``.

    ``S[i, j] = ∫_{Γ_j} G(x_i, y) ds``  y  ``D[i, j] = ∫_{Γ_j} ∂G/∂n_y(x_i, y) ds``,
    con ``n_y`` la normal saliente del cuerpo en ``Γ_j``.

    Si ``handle_self_panel`` es True se asume que ``points`` son los puntos medios de los paneles
    (mismo orden) y se usan las expresiones exactas de ``S_ii`` y ``D_ii``.
    """
    N = boundary.n_panels
    M = points.shape[0]
    S = np.zeros((M, N), dtype=complex)
    D = np.zeros((M, N), dtype=complex)

    nodes, weights = leggauss(order)
    unit_nodes = 0.5 * (nodes + 1.0)
    unit_weights = 0.5 * weights
    dvec = boundary.p1 - boundary.p0
    # nodos de cuadratura de cada panel en el plano: (N, nq, 2)
    quad_points = boundary.p0[:, None, :] + unit_nodes[None, :, None] * dvec[:, None, :]
    quad_weights = unit_weights[None, :] * boundary.lengths[:, None]   # (N, nq)

    # bloque de filas para no materializar (M, N, nq, 2) cuando M y N son grandes
    block = max(1, int(4.0e6 / max(N * order, 1)))
    for start in range(0, M, block):
        stop = min(start + block, M)
        xb = points[start:stop]                                        # (m, 2)
        diff = quad_points[None, :, :, :] - xb[:, None, None, :]       # (m, N, nq, 2)
        r = np.hypot(diff[..., 0], diff[..., 1])
        with np.errstate(invalid="ignore", divide="ignore"):
            G = 0.25j * hankel1(0, k * r)
            dot = np.einsum("mnqd,nd->mnq", diff, boundary.normals)
            dG_dn = -(1j * k / 4.0) * hankel1(1, k * r) * dot / r
        S[start:stop] = np.einsum("mnq,nq->mn", G, quad_weights)
        D[start:stop] = np.einsum("mnq,nq->mn", dG_dn, quad_weights)

    # pares cercanos: cuadratura graduada panel por panel
    dist = _distance_to_panels(points, boundary)
    near = dist < 0.5 * boundary.lengths[None, :]
    if handle_self_panel:
        np.fill_diagonal(near, False)
    for i, j in zip(*np.nonzero(near)):
        t, w = _segment_quadrature(boundary.p0[j], boundary.p1[j], points[i],
                                   order, max_subdivisions)
        qp = boundary.p0[j][None, :] + t[:, None] * dvec[j][None, :]
        diff = qp - points[i]
        r = np.hypot(diff[:, 0], diff[:, 1])
        S[i, j] = np.sum(w * 0.25j * hankel1(0, k * r))
        D[i, j] = np.sum(w * (-(1j * k / 4.0) * hankel1(1, k * r)
                              * (diff @ boundary.normals[j]) / r))

    if handle_self_panel:
        # S_ii por extracción de la singularidad logarítmica y D_ii exactamente nulo
        for i in range(N):
            S[i, i] = _self_panel_single_layer(k, boundary.lengths[i])
        np.fill_diagonal(D, 0.0)

    return S, D


# --------------------------------------------------------------------------- #
# operadores de la derivada normal (Burton-Miller)
# --------------------------------------------------------------------------- #


def _hypersingular_kernels(k: complex, r: np.ndarray, rs: np.ndarray,
                           wx: np.ndarray, wy: np.ndarray, c):
    r"""Núcleos ``W`` (hipersingular) y ``K'`` (doble capa adjunta).

    Con ``w = x − y``, ``r = |w|`` y las normales ``n_x`` (en el punto de colocación) y ``n_y``
    (en el panel):

    * ``K'(x, y) = ∂G/∂n_x = −(i·k/4)·H₁(k·r)·(w·n_x)/r``   (débilmente singular, ``O(1/r)``);
    * ``W(x, y) = ∂²G/∂n_x∂n_y = −(i·k/4)·[ k·H₀(k·r)·a·b + H₁(k·r)·(c − 2·a·b)/r ]``,
      con ``a = (w·n_x)/r``, ``b = (w·n_y)/r``, ``c = n_x·n_y``.

    Se devuelve la descomposición ``W = W_static + W_regular`` con

        W_static = −(1/2π)·(c − 2·a·b)/r²

    el núcleo de Laplace (``O(1/r²)``, integrable sólo como valor principal finito) y ``W_regular``
    suave. Para que la resta no sea numérica se resta **analíticamente** el término dominante de
    ``H₁``::

        H₁(z) = −2i/(π·z) + H₁ˢ(z),   H₁ˢ(z) = H₁(z) + 2i/(π·z) = O(z)

    con lo que el término ``1/r²`` del núcleo de Helmholtz resulta exactamente ``W_static`` y el
    resto queda ``−(i·k²/4)·H₀·a·b − (i·k/4)·H₁ˢ·(c − 2·a·b)/r``, que tiende a ``0`` y es
    integrable con Gauss-Legendre sin cancelación catastrófica. (Restar directamente ``H₁/r`` de
    ``1/r²`` sí la tendría: ambos términos son ``~10¹²`` para ``r ~ 10⁻⁶``.)

    ``rs`` es ``r`` con los ceros reemplazados por ``1`` (evita ``0/0`` en los nodos que caen sobre
    el punto de colocación; ese caso lo resuelve el valor principal finito de :func:`_W_self_panel`).
    """
    r2 = r * r
    with np.errstate(divide="ignore", invalid="ignore"):
        ab = wx * wy / np.where(r2 > 0.0, r2, 1.0)
    W_static = -(1.0 / (2.0 * np.pi)) * (c - 2.0 * ab) / np.where(r2 > 0.0, r2, 1.0)
    hankel1_smooth = hankel1(1, k * rs) + 2j / (np.pi * (k * rs))
    W_regular = (-(1j * k * k / 4.0) * hankel1(0, k * r) * ab
                 - (1j * k / 4.0) * hankel1_smooth * (c - 2.0 * ab) / rs)
    K_adjoint = (1j * k / 4.0) * hankel1(1, k * r) * wx / rs
    return W_static, W_regular, K_adjoint


def _W_self_panel(k: complex, length: float, order: int = 48) -> complex:
    r"""``∫_{Γ_i} W ds`` sobre el panel propio, por el valor principal finito de Hadamard.

    Para el panel recto que contiene a ``x_i`` (punto medio) es ``a = b = 0`` y ``c = 1`` a lo
    largo de todo el panel, así que el núcleo estático vale ``W_static = −1/(2π·s²)`` con ``s`` la
    coordenada a lo largo del panel, y su integral

        f.p.∫_{−h}^{h} −ds/(2π·s²) = −(1/2π)·(−2/h) = 1/(π·h) = 2/(π·L)

    con ``h = L/2``. El resto (``W_regular``, suave y par en ``s``) se integra con Gauss-Legendre.
    Que la suma de esta diagonal con las integrales de ``W_static`` sobre los demás paneles dé
    ``0`` para densidad constante (``W·1 = 0`` es exacto en cualquier curva cerrada) es la
    verificación que se usó para fijar el valor ``2/(π·L)``.
    """
    h = 0.5 * length
    tau, weights = leggauss(order)
    s = h * tau
    rs = np.where(np.abs(s) > 0.0, np.abs(s), 1.0)
    hankel1_smooth = hankel1(1, k * rs) + 2j / (np.pi * (k * rs))
    regular = -(1j * k / 4.0) * hankel1_smooth / rs
    return 2.0 / (np.pi * length) + h * np.sum(weights * regular)


def _normal_derivative_matrices(points: np.ndarray, boundary: CircleBoundary, k: complex,
                                order: int = DEFAULT_QUAD_ORDER,
                                max_subdivisions: int = DEFAULT_MAX_SUBDIVISIONS,
                                handle_self_panel: bool = False):
    """Matrices ``W`` (hipersingular) y ``K'`` (doble capa adjunta), forma ``(M, N)``.

    ``W[i, j] = ∫_{Γ_j} ∂²G/∂n_x∂n_y ds`` y ``K'[i, j] = ∫_{Γ_j} ∂G/∂n_x ds``, con ``n_x`` la
    normal en ``points[i]`` y ``n_y`` la normal saliente del cuerpo en ``Γ_j``. Ambas incluyen
    sólo el valor principal (los saltos ``∓½`` son explícitos en las ecuaciones, igual que el
    ``½`` de ``(½I − D)``).

    Si ``handle_self_panel`` es True se asume que ``points`` son los puntos medios de los paneles
    (mismo orden) y se usan :func:`_W_self_panel` y ``K'_ii = 0`` (exactos, ver docstring del
    módulo). Misma estrategia de cuadratura que :func:`_layer_matrices`: Gauss-Legendre
    vectorizado para los pares lejanos y cuadratura graduada hacia el punto para los cercanos.
    """
    N = boundary.n_panels
    M = points.shape[0]
    W = np.zeros((M, N), dtype=complex)
    K_adjoint = np.zeros((M, N), dtype=complex)

    nodes, weights = leggauss(order)
    unit_nodes = 0.5 * (nodes + 1.0)
    unit_weights = 0.5 * weights
    dvec = boundary.p1 - boundary.p0
    quad_points = boundary.p0[:, None, :] + unit_nodes[None, :, None] * dvec[:, None, :]
    quad_weights = unit_weights[None, :] * boundary.lengths[:, None]      # (N, nq)

    block = max(1, int(4.0e6 / max(N * order, 1)))
    for start in range(0, M, block):
        stop = min(start + block, M)
        xb = points[start:stop]                                        # (m, 2)
        diff = xb[:, None, None, :] - quad_points[None, :, :, :]        # (m, N, nq, 2)
        r2 = np.einsum("mnqd,mnqd->mnq", diff, diff)
        r = np.sqrt(r2)
        rs = np.where(r > 0.0, r, 1.0)
        wx = np.einsum("mnqd,md->mnq", diff, boundary.normals[start:stop])
        wy = np.einsum("mnqd,nd->mnq", diff, boundary.normals)
        c = (boundary.normals[start:stop] @ boundary.normals.T)[:, :, None]
        W_static, W_regular, Kq = _hypersingular_kernels(k, r, rs, wx, wy, c)
        W[start:stop] = np.einsum("mnq,nq->mn", W_static + W_regular, quad_weights)
        K_adjoint[start:stop] = np.einsum("mnq,nq->mn", Kq, quad_weights)

    # pares cercanos: cuadratura graduada panel por panel
    dist = _distance_to_panels(points, boundary)
    near = dist < 0.5 * boundary.lengths[None, :]
    if handle_self_panel:
        np.fill_diagonal(near, False)
    for i, j in zip(*np.nonzero(near)):
        t, w = _segment_quadrature(boundary.p0[j], boundary.p1[j], points[i],
                                   order, max_subdivisions)
        qp = boundary.p0[j][None, :] + t[:, None] * dvec[j][None, :]
        diff = points[i][None, :] - qp
        r = np.hypot(diff[:, 0], diff[:, 1])
        rs = np.where(r > 0.0, r, 1.0)
        nx = boundary.normals[i]
        ny = boundary.normals[j]
        W_static, W_regular, Kq = _hypersingular_kernels(
            k, r, rs, diff @ nx, diff @ ny, complex(nx @ ny))
        W[i, j] = np.sum(w * (W_static + W_regular))
        K_adjoint[i, j] = np.sum(w * Kq)

    if handle_self_panel:
        for i in range(N):
            W[i, i] = _W_self_panel(k, boundary.lengths[i])
        np.fill_diagonal(K_adjoint, 0.0)

    return W, K_adjoint


# --------------------------------------------------------------------------- #
# ensamblado y resolución
# --------------------------------------------------------------------------- #


def assemble_operators(boundary: CircleBoundary, k: complex,
                       quad_order: int = DEFAULT_QUAD_ORDER,
                       max_subdivisions: int = DEFAULT_MAX_SUBDIVISIONS):
    """Ensambla ``S`` y ``D`` sobre el contorno (``(N, N)``), sin invertir nada.

    Útil para inspeccionar condicionamiento o para reutilizar las matrices con varios datos de
    Neumann.
    """
    return _layer_matrices(boundary.midpoints, boundary, k, quad_order,
                           max_subdivisions, handle_self_panel=True)


@dataclass(frozen=True)
class ExteriorSolution:
    """Resultado de :func:`solve_exterior_neumann`.

    Atributos
    ---------
    p_total, p_scattered, p_incident : ``(M,)`` complejo, campos en ``points_eval`` [Pa].
    p_boundary : ``(N,)`` presión total resuelta sobre los paneles [Pa].
    k : número de onda usado [1/m].
    n_panels : incógnitas del sistema (una por panel).
    condition_number : ``cond`` de la matriz del sistema resuelto (``½·I − D`` en la ecuación
        estándar, su combinación con ``η·W`` en Burton-Miller y la matriz aumentada en CHIEF).
    boundary : la malla usada.
    formulation : ``"total"`` o ``"scattered"`` (la incógnita).
    equation : ``"standard"``, ``"burton_miller"`` o ``"chief"`` (la ecuación).
    n_equations : filas del sistema (``N``; ``N + M`` en CHIEF).
    residual : ``‖A·x − b‖ / ‖b‖`` de la solución devuelta (``~1e-16`` en los sistemas cuadrados y
        resueltos directo; el residuo de mínimos cuadrados en CHIEF).
    """

    p_total: np.ndarray
    p_scattered: np.ndarray
    p_incident: np.ndarray
    p_boundary: np.ndarray
    k: complex
    n_panels: int
    condition_number: float
    boundary: CircleBoundary
    formulation: str
    equation: str = "standard"
    n_equations: int | None = None
    residual: float = 0.0


def _check_equation(equation) -> str:
    """Valida la elección de ecuación integral."""
    if equation not in EQUATIONS:
        opciones = ", ".join(repr(e) for e in EQUATIONS)
        raise ValueError(
            f"equation debe ser uno de {opciones}; recibido {equation!r}"
        )
    return equation


def _check_chief_points(chief_points) -> int:
    """Valida la cantidad de puntos interiores de CHIEF."""
    if isinstance(chief_points, bool) or not isinstance(chief_points, (int, np.integer)):
        raise ValueError(
            f"chief_points debe ser un entero >= 1; recibido {chief_points!r}"
        )
    if int(chief_points) < 1:
        raise ValueError(f"chief_points debe ser >= 1; recibido {chief_points!r}")
    return int(chief_points)


def _check_bm_eta(bm_eta, k: complex) -> complex:
    """Valida (o elige) el ``η`` de Burton-Miller; por defecto ``η = i/k``."""
    if bm_eta is None:
        if abs(k) == 0.0:
            raise ValueError(
                "burton_miller con k = 0 no está definido (η = i/k); "
                "usar frequency > 0 o pasar bm_eta"
            )
        return 1j / k
    try:
        eta = complex(bm_eta)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"bm_eta debe ser un escalar complejo no nulo; recibido {bm_eta!r}"
        ) from exc
    if not np.isfinite(eta) or eta == 0:
        raise ValueError(
            f"bm_eta debe ser un escalar complejo finito y no nulo; recibido {bm_eta!r}"
        )
    return eta


def _chief_interior_points(boundary: CircleBoundary, m: int) -> np.ndarray:
    """``(M, 2)`` puntos interiores uniformes sobre ``r = 0.5·a`` para las filas de CHIEF."""
    theta = 2.0 * np.pi * np.arange(m) / m
    radius = CHIEF_RADIUS_FRACTION * boundary.radius
    return radius * np.column_stack([np.cos(theta), np.sin(theta)])


def _relative_residual(matrix: np.ndarray, solution: np.ndarray, rhs: np.ndarray) -> float:
    """``‖A·x − b‖ / ‖b‖`` (evita dividir por cero si ``b`` es nulo)."""
    scale = float(np.linalg.norm(rhs))
    if scale == 0.0:
        scale = 1.0
    return float(np.linalg.norm(matrix @ solution - rhs) / scale)


def solve_exterior_neumann(points_eval, frequency: float, material: Material, radius: float,
                           n_panels: int | None = None, angle: float = 0.0,
                           p_inc: float = 1.0, neumann_data=None,
                           boundary: CircleBoundary | None = None,
                           formulation: str = "total",
                           quad_order: int = DEFAULT_QUAD_ORDER,
                           equation: str = "standard",
                           chief_points: int | None = None,
                           bm_eta: complex | None = None) -> ExteriorSolution:
    """Resuelve la dispersión por un cilindro rígido (Neumann exterior) con BEM.

    Parámetros
    ----------
    points_eval : ``(M, 2)`` puntos del fluido donde se quiere el campo (``r > radius``).
    frequency : frecuencia [Hz].
    material : medio exterior (``materials.Material``), de ahí salen ``c`` y ``α(f)``.
    radius : radio del cilindro [m] (se ignora si se pasa ``boundary``).
    n_panels : paneles; ``None`` ⇒ :func:`panels_for_frequency`.
    angle : ángulo de incidencia [rad].
    p_inc : amplitud de la onda plana incidente [Pa].
    neumann_data : dato ``q = ∂p/∂n`` sobre los paneles. ``None`` o ``0`` ⇒ pared rígida. Puede ser
        un escalar o un array ``(N,)`` alineado con ``boundary.midpoints``.
    boundary : contorno ya construido (opcional, para reutilizar la malla).
    formulation : ``"total"`` (incógnita ``p``, por defecto) o ``"scattered"`` (incógnita ``p_s``,
        con ``q_s = q − q_inc``); las dos deben dar el mismo campo y se usan como verificación
        cruzada.
    quad_order : orden de Gauss-Legendre por subsegmento.
    equation : ecuación integral a resolver, ver "Elección de ecuación" en el docstring del módulo.
        ``"standard"`` (por defecto, ``(½I − D)·p = p_inc − S·q``) sufre resonancias espurias en
        los ``ka`` con ``J_n(ka) = 0``; ``"burton_miller"`` y ``"chief"`` las eliminan a costa de
        un poco más de trabajo por frecuencia. En un barrido de ``ka`` hay que usar una de las dos.
    chief_points : cantidad de puntos interiores de CHIEF (``equation="chief"``). ``None`` ⇒
        :data:`DEFAULT_CHIEF_POINTS` (4). Sólo válido con ``equation="chief"``.
    bm_eta : ``η`` de Burton-Miller (``equation="burton_miller"``). ``None`` ⇒ ``η = i/k``, que es
        el valor que evita el fallo en *todas* las resonancias. Sólo válido con
        ``equation="burton_miller"``.

    Devuelve
    --------
    :class:`ExteriorSolution`.

    Ejemplos
    --------
    En un barrido de frecuencia conviene Burton-Miller (mismo costo que CHIEF, sistema cuadrado)::

        sol = solve_exterior_neumann(pts, f, WATER, a, equation="burton_miller")

    y para auditar el resultado, ``sol.n_equations``, ``sol.residual`` y
    ``sol.condition_number`` dicen cuántas filas tenía el sistema, cuánto lo satisface la solución
    y cuánto está condicionado.
    """
    if formulation not in ("total", "scattered"):
        raise ValueError("formulation debe ser 'total' o 'scattered'")
    equation = _check_equation(equation)
    if chief_points is not None and equation != "chief":
        raise ValueError(
            f"chief_points solo se usa con equation='chief'; equation={equation!r}"
        )
    if bm_eta is not None and equation != "burton_miller":
        raise ValueError(
            f"bm_eta solo se usa con equation='burton_miller'; equation={equation!r}"
        )
    if not isinstance(material, Material):
        raise TypeError("material debe ser un materials.Material")

    pts = np.asarray(points_eval, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 2:
        raise ValueError(f"points_eval debe tener forma (M, 2); recibido {pts.shape}")

    k = material.wavenumber(frequency)

    if boundary is None:
        if n_panels is None:
            n_panels = panels_for_frequency(radius, k)
        boundary = circle_boundary(radius, n_panels)
    else:
        radius = boundary.radius

    r = np.hypot(pts[:, 0], pts[:, 1])
    if np.any(r <= radius * (1.0 + 1e-12)):
        raise ValueError(
            "points_eval debe quedar fuera del cilindro (r > radius); "
            f"r mínimo = {r.min():.6g} m, radius = {radius:.6g} m"
        )

    n_unknowns = boundary.n_panels
    q = np.zeros(n_unknowns, dtype=complex) if neumann_data is None else np.asarray(
        neumann_data, dtype=complex)
    if q.ndim == 0:
        q = np.full(n_unknowns, complex(q))
    if q.shape != (n_unknowns,):
        raise ValueError(f"neumann_data debe ser escalar o de forma ({n_unknowns},)")

    S, D = assemble_operators(boundary, k, quad_order)

    cos_a, sin_a = np.cos(angle), np.sin(angle)
    mid = boundary.midpoints
    p_inc_b = p_inc * np.exp(1j * k * (mid[:, 0] * cos_a + mid[:, 1] * sin_a))
    # ∂p_inc/∂n sobre el contorno (normal saliente del cuerpo)
    q_inc = p_inc * 1j * k * (boundary.normals @ np.array([cos_a, sin_a])) * np.exp(
        1j * k * (mid[:, 0] * cos_a + mid[:, 1] * sin_a))

    A = 0.5 * np.eye(n_unknowns) - D

    # -- elección de ECUACIÓN -------------------------------------------------
    #
    # En las tres variantes se resuelve la incógnita ``u`` que pide ``formulation`` (``p`` total o
    # ``p_s``); en la estándar el camino es exactamente el de siempre.
    if equation == "standard":
        if formulation == "total":
            rhs = p_inc_b - S @ q
            matrix = A
            n_equations = n_unknowns
            u = np.linalg.solve(matrix, rhs)
            S_eval, D_eval = _layer_matrices(pts, boundary, k, quad_order)
            p_boundary = u
            p_scattered = D_eval @ p_boundary - S_eval @ q
        else:
            q_s = q - q_inc
            rhs = -(S @ q_s)
            matrix = A
            n_equations = n_unknowns
            u = np.linalg.solve(matrix, rhs)
            S_eval, D_eval = _layer_matrices(pts, boundary, k, quad_order)
            p_boundary = u + p_inc_b
            p_scattered = D_eval @ u - S_eval @ q_s
        residual = _relative_residual(matrix, u, rhs)
    else:
        # Burton-Miller y CHIEF comparten los operadores de la derivada normal:
        #   W  = ∂²G/∂n_x∂n_y  (hipersingular, regularizado por sustracción estática)
        #   K' = ∂G/∂n_x       (doble capa adjunta, valor principal; K'_ii = 0)
        # con las identidades  W·u = (−½I + K')·q_s  y  (D − S)·u = 0  en el interior, válidas para
        # todo campo radiante u (verificadas contra campos fabricados en el docstring del módulo).
        W, K_adjoint = _normal_derivative_matrices(
            mid, boundary, k, quad_order, handle_self_panel=True)
        half_neg = -0.5 * np.eye(n_unknowns)
        q_shift = q - q_inc                      # dato de Neumann del campo dispersado

        if formulation == "total":
            rhs_e = p_inc_b - S @ q
            incident_bm = W @ p_inc_b            # [(½I − D) + ηW]·p_inc = (p_inc − S·q_inc) + η·W·p_inc
            q_eval = q
        else:
            rhs_e = -(S @ q_shift)
            incident_bm = np.zeros(n_unknowns, dtype=complex)
            q_eval = q_shift

        if equation == "burton_miller":
            eta = _check_bm_eta(bm_eta, k)
            matrix = A + eta * W
            rhs = rhs_e + eta * ((half_neg + K_adjoint) @ q_shift + incident_bm)
            u = np.linalg.solve(matrix, rhs)
            n_equations = n_unknowns
        else:  # chief
            m_points = (_check_chief_points(chief_points) if chief_points is not None
                        else DEFAULT_CHIEF_POINTS)
            x_int = _chief_interior_points(boundary, m_points)
            S_int, D_int = _layer_matrices(x_int, boundary, k, quad_order)
            # filas extra: (D − S)·u = 0 en el interior (teorema del campo nulo); con la incógnita
            # total el término incidente D_int·p_inc se suma al miembro derecho.
            matrix = np.vstack([A, D_int])
            rhs_int = S_int @ q_shift
            if formulation == "total":
                rhs_int = rhs_int + D_int @ p_inc_b
            rhs = np.concatenate([rhs_e, rhs_int])
            u, *_ = np.linalg.lstsq(matrix, rhs, rcond=None)
            n_equations = n_unknowns + m_points

        S_eval, D_eval = _layer_matrices(pts, boundary, k, quad_order)
        if formulation == "total":
            p_boundary = u
        else:
            p_boundary = u + p_inc_b
        p_scattered = D_eval @ u - S_eval @ q_eval
        residual = _relative_residual(matrix, u, rhs)

    p_inc_eval = p_inc * np.exp(1j * k * (pts[:, 0] * cos_a + pts[:, 1] * sin_a))

    return ExteriorSolution(
        p_total=p_inc_eval + p_scattered,
        p_scattered=p_scattered,
        p_incident=p_inc_eval,
        p_boundary=p_boundary,
        k=k,
        n_panels=n_unknowns,
        condition_number=float(np.linalg.cond(matrix)),
        boundary=boundary,
        formulation=formulation,
        equation=equation,
        n_equations=n_equations,
        residual=residual,
    )
