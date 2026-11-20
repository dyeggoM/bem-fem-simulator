"""Campo 2D de paneles radiantes **finitos** en espacio libre (BEM de capa simple).

Este módulo calcula lo mismo que el motor del dashboard (:mod:`ultrasound_sim.dashboard.field`)
pero **sin baffle**: los paneles son una lámina delgada de ancho finito, con el fluido a ambos
lados, que radia hacia adelante y hacia atrás y difracta en sus bordes. Sirve para medir cuánto del
campo que muestra el dashboard se debe al supuesto de baffle rígido infinito.

Modelo y derivación
-------------------
Convención del repo (``docs/decisions.md`` D4): tiempo ``exp(−i·ω·t)``, ``∇²p + k²p = 0``,
``k = ω/c + i·α`` (:meth:`materials.Material.wavenumber`), Green saliente 2D

    G(x, y) = (i/4)·H₀⁽¹⁾(k·|x − y|),        (∇² + k²)·G = −δ(x − y)

y Euler linealizado ``ρ·∂u/∂t = −∇p`` con ``u = ∇p/(i·ω·ρ)``.

Un panel es un segmento sobre ``x = 0`` con velocidad normal prescrita ``v`` (misma convención que
``dashboard.field.Source``: ``y``, ``width``, ``amp``, ``phase``). La lámina es de espesor nulo y
**radia por las dos caras**: la cara ``+x`` empuja el fluido hacia ``+x`` y la cara ``−x`` lo empuja
hacia ``−x``, de modo que

    ∂p/∂x|₀₊ = i·ω·ρ·v ,        ∂p/∂x|₀₋ = −i·ω·ρ·v .

Con la representación de capa simple ``p(x) = ∫_S σ(y)·G(x, y) ds`` (que cumple
``(∇² + k²)p = −σ·δ_S``) el salto de derivada normal a través de la lámina es
``∂p/∂x|₀₊ − ∂p/∂x|₀₋ = −σ``, entonces

    σ(y) = −2·i·ω·ρ·v(y)        y        p(x) = −2·i·ω·ρ · Σ_j ∫_{panel j} v_j · G(x, y) dy .

Es la única integral que hay que calcular: **este problema no tiene sistema lineal que resolver**, y
eso es un resultado físico, no una simplificación. La simetría de la lámina (misma ``v`` en las dos
caras) implica ``p(−x, y) = p(x, y)``, o sea que la presión es continua a través de la lámina y su
salto ``[p] = 0``. No hay cuerpo cerrado, no hay dominio interior y por lo tanto no aparece la
incógnita de presión superficial que resuelve un BEM de contorno cerrado como
:mod:`ultrasound_sim.bem.helmholtz2d`: acá la densidad de capa simple queda determinada por el salto
de velocidad y no queda nada que invertir. Todo el trabajo numérico es la cuadratura de esa capa
simple, con la singularidad logarítmica de ``H₀`` resuelta por refinamiento graduado hacia el punto de
evaluación (ver "Cuadratura").

Por eso ``equation`` sólo acepta ``"single_layer"``: Burton-Miller y CHIEF existen para eliminar las
resonancias espurias de los problemas **interiores** de un cuerpo cerrado (``J_n(ka) = 0`` y
``J_n′(ka) = 0``); una lámina abierta no encierra dominio interior y no tiene esas resonancias. Pasar
``"burton_miller"`` o ``"chief"`` levanta ``ValueError`` con esa explicación en vez de fingir que se
eligió una ecuación distinta.

Contraste con la integral de Rayleigh con baffle (el resultado que se busca)
--------------------------------------------------------------------------
El dashboard usa (``aperture.py``, ``dashboard/field.py``)

    p_baffle(x) = +2·i·ω·ρ · Σ_j ∫ v_j · G(x, y) dy ,     x > 0

donde el factor 2 es la imagen del baffle rígido. Comparando con la expresión de arriba:

    p_free(x) = −p_baffle(x)      para todo x > 0

es decir **el mismo módulo y fase opuesta**, no 6 dB menos. La razón es que el factor 2 del baffle y
el factor 2 de la lámina bilateral son *el mismo* 2: en el baffle lo aporta la fuente imagen (que para
una apertura plana está en la misma posición que la fuente real, por eso duplica la presión); en la
lámina libre lo aportan las dos caras, que también radian en fase hacia ``x > 0``. Consecuencias
medibles (test ``tests/test_bem_field.py``):

* ``|p_free|/|p_baffle| = 1`` sobre el eje y en el foco (0.00 dB) para la misma ``v`` prescrita;
* el precio real del baffle es el **semicampo de atrás** y la **potencia**: la lámina libre radia
  hacia ``x > 0`` la misma potencia que el baffle *y* otra igual hacia ``x < 0`` (2× la total,
  +3 dB), y el dashboard ni siquiera puede representar ``x < 0``. Normalizando por potencia total
  radiada en vez de por velocidad, la presión al frente queda 3 dB abajo;
* los 6 dB esperados aparecen si se compara el baffle contra una densidad **unilateral**
  (``−i·ω·ρ·∫v·G``, la mitad de ``σ``); esa densidad no satisface el salto de velocidad de una lámina
  que radia por las dos caras (daría ``v/2`` en cada una), así que no es el modelo de este módulo.

Cuidado con el signo global: ``aperture.py`` documenta ``p = −(ω·ρ/2)·v·H₀`` para una apertura
puntual, que es exactamente ``−p_free``. Es un convenio global (equivalente a ``v → −v``), no una
diferencia física; las comparaciones en dB usan ``|p|`` y el desfase de π se reporta aparte.

Variante con dorso rígido (explícitamente **no** implementada)
--------------------------------------------------------------
Si el transductor tuviera el dorso rígido (una sola cara radiante, ``∂p/∂x|₀₋ = 0`` en vez de
``−i·ω·ρ·v``), el problema deja de ser simétrico: el salto de presión ``[p]`` ya no es nulo, hay que
resolver una ecuación integral de primera especie (hipersingular) sobre una superficie abierta y el
campo al frente queda ~6 dB por debajo del de baffle para la misma velocidad, porque la cara trasera
no radia y el baffle duplica la fuente con su imagen. Ese caso no es "radiación a ambos lados" y este
módulo no lo implementa ni lo reporta como si lo hiciera: los 6 dB que aparecen en la literatura
corresponden a esa geometría (o a una densidad unilateral), no a la lámina bilateral.

Límites analíticos (validación)
-------------------------------
* **Fuente lineal.** Para ``ancho → 0`` la lámina es una fuente de línea con caudal por unidad de
  longitud ``Q = 2·width·v`` [m²/s] y

      p(r) = −i·ω·ρ·Q·G = (ω·ρ·Q/4)·H₀⁽¹⁾(k·r) .

  Sale de ``(∇² + k²)p = i·ω·ρ·Q·δ`` (continuidad + Euler) o, equivalente, del límite ``k·r → 0`` de
  la solución del cilindro.
* **Cilindro pulsante.** Cilindro de radio ``a`` con velocidad radial ``v_c``: con ``p = A·H₀⁽¹⁾(k·r)``
  y ``u_r = A·k·H₀⁽¹⁾′(k·r)/(i·ω·ρ)``, la condición ``u_r(a) = v_c`` da

      A = i·ρ·c·v_c / H₀⁽¹⁾′(k·a)   ⇒   p(r) = i·ρ·c·v_c · H₀⁽¹⁾(k·r) / H₀⁽¹⁾′(k·a) ,   r > a.

  El factor ``+i`` es correcto con la convención de este repo; la forma ``−ρ·c·v·H₀/H₀′`` que circula
  en algunos apuntes tiene un ``−i`` de menos (no reproduce ni el límite ``k·a → 0`` ni
  ``Q = 2·π·a·v_c``). Con ``Q`` fijo y ``a → 0`` las dos expresiones coinciden como ``O((k·a)²·ln k·a)``,
  y el test lo verifica: a ``k·a = 0.031`` la diferencia es ``2.2e-3``.

Cuadratura
----------
Gauss-Legendre de ``n_per_panel`` nodos por panel. El integrando es suave en ``x > 0`` salvo que el
punto caiga cerca de la lámina: ahí ``H₀(k·r)`` tiene el pico logarítmico de ``r → 0``. Cada panel se
subdivide en ``n_sub`` tramos iguales con

    n_sub = ceil(log₂( (width/2) / d_min )),        1 ≤ n_sub ≤ max_subdivisions

con ``d_min`` la distancia mínima del retículo de evaluación a ese panel (``0`` si lo cruza). Los
puntos que caen *sobre* un panel (``|x| < 1e-4·λ`` dentro de su ancho) se rechazan con ``ValueError``:
ahí el valor depende del detalle de la cuadratura y no es un dato útil.

Costo medido: 80×80 puntos, 64 paneles de ``0.625 mm`` (λ/10 a 250 kHz), ``n_per_panel = 16``:
≈ 2–5 s según la carga de la máquina (1024 nodos de cuadratura por punto; ver
``test_costo_malla_80x80``). El costo es ``O(n_puntos · n_paneles · n_per_panel)``, dominado por las
``H₀⁽¹⁾``; con ``n_per_panel = 8`` el resultado ya no cambia y el tiempo baja a la mitad.

Supuestos declarados
--------------------
* Panel de espesor nulo, rígido en su propio plano y sin interacción fluido–estructura: la velocidad
  normal es un dato prescrito, igual que en el dashboard. La mecánica del transductor no está en el
  modelo.
* Medio homogéneo e isótropo con ``c``, ``ρ`` y ``α(f)`` de ``materials.Material`` (no se introduce
  ninguna constante nueva; con ``k`` complejo la atenuación entra sola).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import hankel1

from ..materials import Material

__all__ = [
    "EQUATIONS",
    "Panel",
    "free_space_field",
    "free_space_pressure",
    "line_source_pressure",
    "panels_from_sources",
]

#: Ecuaciones disponibles. Sólo hay una capa simple porque el salto de presión es nulo (ver el
#: docstring del módulo); Burton-Miller y CHIEF no aplican a una superficie abierta.
EQUATIONS = ("single_layer",)

#: Nodos de Gauss-Legendre por panel (y por tramo, si el panel se subdivide).
DEFAULT_NODES_PER_PANEL = 16

#: Máxima subdivisión de un panel para los puntos cercanos.
DEFAULT_MAX_SUBDIVISIONS = 12

#: Distancia mínima admitida entre un punto de evaluación y la lámina, en longitudes de onda.
MIN_DISTANCE_IN_WAVELENGTHS = 1e-4

#: Elementos por bloque en la evaluación vectorizada: acota la memoria de ``(puntos, nodos)`` y
#: mantiene los temporales en caché (bloques de ~10⁵ elementos son los más rápidos, medido).
_MAX_BLOCK_ELEMENTS = 100_000


# --------------------------------------------------------------------------- #
# paneles
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Panel:
    """Panel radiante plano sobre ``x = 0`` (misma convención que ``dashboard.field.Source``).

    Atributos
    ---------
    y : centro del panel sobre el eje ``y`` [m].
    width : ancho del panel [m]; cubre ``y ± width/2``.
    amp : módulo de la velocidad normal prescrita ``|v_n|`` [m/s].
    phase : fase de la velocidad normal [rad]; ``v_n = amp·exp(i·phase)``.
    """

    y: float
    width: float
    amp: float = 1.0
    phase: float = 0.0

    @property
    def velocity(self) -> complex:
        """Velocidad normal compleja ``v_n = amp·exp(i·phase)`` [m/s]."""
        return complex(self.amp * math.cos(self.phase), self.amp * math.sin(self.phase))

    @classmethod
    def from_dict(cls, d: Mapping) -> "Panel":
        """Construye un panel desde un diccionario con las claves del dashboard.

        ``y`` y ``width`` son obligatorios; ``amp`` y ``phase`` toman los mismos valores por defecto
        que ``dashboard.field.Source.from_dict`` (``1.0`` y ``0.0``).
        """
        if "y" not in d:
            raise ValueError("cada panel necesita la clave 'y' (centro sobre el eje y, en metros)")
        if "width" not in d:
            raise ValueError("cada panel necesita la clave 'width' (ancho en metros)")
        return cls(
            y=float(d["y"]),
            width=float(d["width"]),
            amp=float(d.get("amp", 1.0)),
            phase=float(d.get("phase", 0.0)),
        )


def _as_panel(item) -> Panel:
    """Convierte un diccionario, un :class:`Panel` o un ``dashboard.field.Source`` en ``Panel``."""
    if isinstance(item, Panel):
        return item
    if isinstance(item, Mapping):
        return Panel.from_dict(item)
    # duck typing: acepta dashboard.field.Source y cualquier objeto con y/width/amp/phase
    attrs = ("y", "width", "amp", "phase")
    if all(hasattr(item, a) for a in attrs):
        try:
            return Panel(y=float(item.y), width=float(item.width),
                         amp=float(item.amp), phase=float(item.phase))
        except (TypeError, ValueError) as exc:
            raise ValueError(f"panel con atributos no numéricos: {item!r}") from exc
    raise ValueError(
        "cada panel debe ser un dict con 'y'/'width' (y opcionalmente 'amp'/'phase'), un Panel "
        f"o un objeto con esos atributos; recibido {type(item).__name__}"
    )


def panels_from_sources(sources: Iterable) -> list[Panel]:
    """Convierte una lista de fuentes del dashboard (``Source`` o dicts) en :class:`Panel`."""
    return [_as_panel(s) for s in sources]


# --------------------------------------------------------------------------- #
# validación
# --------------------------------------------------------------------------- #


def _check_material(material) -> Material:
    if not isinstance(material, Material):
        raise ValueError(
            "material debe ser un materials.Material (por ejemplo materials.WATER); "
            f"recibido {type(material).__name__}"
        )
    return material


def _check_frequency(frequency) -> float:
    try:
        value = float(frequency)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"frequency debe ser un número en Hz; recibido {frequency!r}") from exc
    if not math.isfinite(value):
        raise ValueError(f"frequency debe ser finita; recibido {frequency!r}")
    if value <= 0.0:
        raise ValueError(
            f"frequency debe ser > 0 Hz (k = ω/c depende de ella); recibido {frequency!r}"
        )
    return value


def _check_axis(values, name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 1:
        raise ValueError(f"{name} debe ser un array 1-D de coordenadas [m]; recibido forma {arr.shape}")
    if arr.size == 0:
        raise ValueError(f"{name} está vacío: no hay puntos de evaluación")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contiene NaN o inf; todas las coordenadas deben ser finitas [m]")
    return arr


def _check_panels(panels: Sequence) -> list[Panel]:
    if panels is None:
        raise ValueError("panels es None: hace falta al menos un panel radiante")
    if isinstance(panels, Mapping):
        panels = [panels]          # se acepta un único panel en forma de dict
    items = list(panels)
    if not items:
        raise ValueError("panels está vacío: hace falta al menos un panel radiante")
    out = []
    for j, item in enumerate(items):
        panel = _as_panel(item)
        for name in ("y", "width", "amp", "phase"):
            value = getattr(panel, name)
            if not math.isfinite(value):
                raise ValueError(
                    f"panel {j}: '{name}' debe ser finito; recibido {value!r} "
                    "(y y width en metros, phase en radianes, amp en m/s)"
                )
        if panel.width <= 0.0:
            raise ValueError(f"panel {j}: width debe ser > 0 m; recibido {panel.width!r}")
        out.append(panel)
    return out


def _check_equation(equation) -> str:
    if equation not in EQUATIONS:
        raise ValueError(
            f"equation debe ser {EQUATIONS[0]!r}; recibido {equation!r}. "
            "Una lámina abierta con velocidad prescrita no tiene sistema lineal que resolver "
            "(el salto de presión es nulo por simetría), así que Burton-Miller y CHIEF —que "
            "eliminan resonancias interiores de cuerpos cerrados— no aplican: no hay dominio "
            "interior ni resonancias espurias que regularizar."
        )
    return equation


def _check_options(n_per_panel, max_subdivisions) -> tuple[int, int]:
    if isinstance(n_per_panel, bool) or not isinstance(n_per_panel, (int, np.integer)):
        raise ValueError(f"n_per_panel debe ser un entero >= 2; recibido {n_per_panel!r}")
    if int(n_per_panel) < 2:
        raise ValueError(f"n_per_panel debe ser >= 2; recibido {n_per_panel!r}")
    if isinstance(max_subdivisions, bool) or not isinstance(max_subdivisions, (int, np.integer)):
        raise ValueError(
            f"max_subdivisions debe ser un entero >= 1; recibido {max_subdivisions!r}"
        )
    if int(max_subdivisions) < 1:
        raise ValueError(f"max_subdivisions debe ser >= 1; recibido {max_subdivisions!r}")
    return int(n_per_panel), int(max_subdivisions)


def _reject_points_on_sheet(points: np.ndarray, panels: list[Panel], wavelength: float) -> None:
    """Rechaza puntos demasiado cerca de la lámina (singularidad logarítmica de ``H₀``)."""
    tol = MIN_DISTANCE_IN_WAVELENGTHS * wavelength
    near = np.abs(points[:, 0]) < tol
    if not np.any(near):
        return
    for j, panel in enumerate(panels):
        inside = np.abs(points[:, 1] - panel.y) <= 0.5 * panel.width
        bad = np.nonzero(near & inside)[0]
        if bad.size:
            ejemplos = ", ".join(
                f"({points[q, 0]:.3e}, {points[q, 1]:.3e}) m" for q in bad[:3])
            raise ValueError(
                "hay puntos de evaluación sobre la lámina (|x| < "
                f"{MIN_DISTANCE_IN_WAVELENGTHS:g}·λ = {tol:.3e} m dentro del ancho del panel "
                f"{j}): {ejemplos}. Ahí H₀⁽¹⁾(k·r) tiene una singularidad logarítmica y el "
                "valor depende de la cuadratura; evaluar a |x| mayor."
            )


# --------------------------------------------------------------------------- #
# cuadratura
# --------------------------------------------------------------------------- #


def _panel_quadrature(panels: list[Panel], points: np.ndarray,
                      n_per_panel: int, max_subdivisions: int):
    """Nodos, pesos [m] y velocidades [m/s] de todos los paneles, concatenados.

    Cada panel se subdivide en ``n_sub`` tramos iguales según la distancia mínima del retículo al
    panel, para resolver el pico logarítmico de ``H₀`` cuando hay puntos muy cerca de la lámina.
    """
    base_nodes, base_weights = leggauss(n_per_panel)
    unit_nodes = 0.5 * (base_nodes + 1.0)          # en [0, 1]
    unit_weights = 0.5 * base_weights

    y_nodes, weights, velocities = [], [], []
    abs_x = np.abs(points[:, 0])                                 # (M,)
    for panel in panels:
        half = 0.5 * panel.width
        # distancia de cada punto al segmento (0 si la proyección en y cae dentro del panel)
        dy_out = np.maximum(np.abs(points[:, 1] - panel.y) - half, 0.0)
        d_min = float(np.min(np.hypot(abs_x, dy_out)))

        if d_min >= half:
            n_sub = 1
        else:
            ratio = half / max(d_min, 1e-300)
            n_sub = int(np.clip(np.ceil(np.log2(ratio)), 1, max_subdivisions))

        edges = np.linspace(panel.y - half, panel.y + half, n_sub + 1)
        velocity = panel.velocity
        for a, b in zip(edges[:-1], edges[1:]):
            length = b - a
            y_nodes.append(a + length * unit_nodes)
            weights.append(length * unit_weights)
            velocities.append(np.full(n_per_panel, velocity, dtype=complex))

    return (np.concatenate(y_nodes), np.concatenate(weights), np.concatenate(velocities))


# --------------------------------------------------------------------------- #
# API pública
# --------------------------------------------------------------------------- #


def line_source_pressure(points, source_strength: float, frequency: float,
                          material: Material) -> np.ndarray:
    """Presión de una fuente lineal 2D de caudal ``Q`` en el origen [Pa].

    ``p(r) = −i·ω·ρ·Q·G = (ω·ρ·Q/4)·H₀⁽¹⁾(k·r)``, que es el límite ``ancho → 0`` de
    :func:`free_space_field` y la referencia de la validación.

    Parámetros
    ----------
    points : ``(M, 2)`` posiciones de evaluación [m].
    source_strength : caudal por unidad de longitud ``Q`` [m²/s]. Para un panel de ancho ``w`` y
        velocidad normal ``v`` vale ``Q = 2·w·v`` (las dos caras).
    frequency : frecuencia [Hz] (> 0).
    material : medio (``materials.Material``); de ahí salen ``c``, ``ρ`` y ``α(f)``.

    Devuelve
    --------
    ``(M,)`` complejo [Pa].
    """
    frequency = _check_frequency(frequency)
    material = _check_material(material)
    pts = np.asarray(points, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 2:
        raise ValueError(f"points debe tener forma (M, 2) en metros; recibido {pts.shape}")
    if not np.all(np.isfinite(pts)):
        raise ValueError("points contiene NaN o inf")
    try:
        strength = float(source_strength)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"source_strength debe ser un número en m²/s; recibido {source_strength!r}"
        ) from exc
    if not math.isfinite(strength):
        raise ValueError(f"source_strength debe ser finito; recibido {source_strength!r}")

    k = complex(material.wavenumber(frequency))
    omega = 2.0 * math.pi * frequency
    r = np.hypot(pts[:, 0], pts[:, 1])
    return -1j * omega * material.rho * strength * 0.25j * hankel1(0, k * r)


def _check_points(points) -> np.ndarray:
    pts = np.asarray(points, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 2:
        raise ValueError(f"points debe tener forma (M, 2) en metros; recibido {pts.shape}")
    if pts.shape[0] == 0:
        raise ValueError("points está vacío: no hay puntos de evaluación")
    if not np.all(np.isfinite(pts)):
        raise ValueError("points contiene NaN o inf; todas las coordenadas deben ser finitas [m]")
    return pts


def free_space_pressure(points, panels, frequency: float, material: Material, *,
                        n_per_panel: int = DEFAULT_NODES_PER_PANEL,
                        max_subdivisions: int = DEFAULT_MAX_SUBDIVISIONS,
                        equation: str = "single_layer") -> np.ndarray:
    """Presión [Pa] en puntos arbitrarios: versión punto a punto de :func:`free_space_field`.

    Misma física y mismos parámetros que :func:`free_space_field`; acá los puntos de evaluación son
    ``(M, 2)`` y pueden estar a cualquier lado de la lámina (``x > 0`` o ``x < 0``). Es la forma
    natural para integrar sobre una curva (potencia en campo lejano) o para verificar reciprocidad.

    Parámetros
    ----------
    points : ``(M, 2)`` posiciones de evaluación ``(x, y)`` [m].
    panels, frequency, material, n_per_panel, max_subdivisions, equation : ver
        :func:`free_space_field`.

    Devuelve
    --------
    ``(M,)`` complejo [Pa].
    """
    frequency = _check_frequency(frequency)
    material = _check_material(material)
    _check_equation(equation)
    n_per_panel, max_subdivisions = _check_options(n_per_panel, max_subdivisions)
    panel_list = _check_panels(panels)
    pts = _check_points(points)

    wavelength = material.wavelength(frequency)
    _reject_points_on_sheet(pts, panel_list, wavelength)

    k = complex(material.wavenumber(frequency))
    omega = 2.0 * math.pi * frequency

    y_nodes, weights, velocities = _panel_quadrature(
        panel_list, pts, n_per_panel, max_subdivisions)
    # amplitud y peso juntos: Σ_j v_j·w_j·(i/4)·H₀(k·r), con el prefactor −2iωρ de la lámina
    amplitude = velocities * weights
    prefactor = 0.25j * (-2j * omega * material.rho)

    out = np.empty(pts.shape[0], dtype=complex)
    abs_x = np.abs(pts[:, 0])
    # |x| para que la simetría p(−x) = p(x) sea exacta también en punto flotante.
    # Bloques de puntos para no materializar (M, nodos) entero cuando la malla es grande.
    block = max(1, int(_MAX_BLOCK_ELEMENTS // max(y_nodes.size, 1)))
    for start in range(0, pts.shape[0], block):
        stop = min(start + block, pts.shape[0])
        r = np.hypot(abs_x[start:stop, None], pts[start:stop, 1][:, None] - y_nodes[None, :])
        out[start:stop] = (amplitude[None, :] * hankel1(0, k * r)).sum(axis=1)

    out *= prefactor
    if not np.all(np.isfinite(out)):
        raise ValueError(
            "el campo calculado contiene NaN o inf: revisar que ningún punto esté sobre la lámina "
            "ni que las entradas tengan magnitudes absurdas"
        )
    return out


def free_space_field(xs, ys, panels, frequency: float, material: Material, *,
                     n_per_panel: int = DEFAULT_NODES_PER_PANEL,
                     max_subdivisions: int = DEFAULT_MAX_SUBDIVISIONS,
                     equation: str = "single_layer") -> np.ndarray:
    """Presión de un conjunto de paneles finitos radiantes en espacio libre 2D (sin baffle).

    Calcula ``p(x) = −2·i·ω·ρ · Σ_j ∫_{panel j} v_j·G(x, y) dy`` con ``G = (i/4)·H₀⁽¹⁾``, que es la
    solución exacta (capa simple) de una lámina delgada con velocidad normal prescrita en las dos
    caras: los paneles radian hacia ``x > 0`` **y** hacia ``x < 0``, y difractan en sus bordes. Ver
    el docstring del módulo para la derivación y para el contraste con la integral de Rayleigh con
    baffle (mismo módulo, fase opuesta).

    Parámetros
    ----------
    xs : ``(nx,)`` coordenadas ``x`` de la malla [m]. Pueden ser negativas (detrás de la lámina);
        no pueden caer sobre un panel (``|x| < 1e-4·λ`` dentro de su ancho).
    ys : ``(ny,)`` coordenadas ``y`` de la malla [m].
    panels : secuencia de paneles. Cada uno puede ser un :class:`Panel`, un ``dict`` con
        ``{"y", "width", "amp", "phase"}`` (misma convención que ``dashboard.field.Source``) o un
        objeto con esos atributos. ``y`` y ``width`` en metros, ``amp`` en m/s, ``phase`` en rad.
    frequency : frecuencia [Hz] (> 0).
    material : medio (``materials.Material``); de ahí salen ``c``, ``ρ`` y ``α(f)``.
    n_per_panel : nodos de Gauss-Legendre por panel (o por tramo si el panel se subdivide).
    max_subdivisions : tope de subdivisión de un panel para puntos cercanos a la lámina.
    equation : sólo ``"single_layer"``; ver el docstring del módulo (Burton-Miller y CHIEF no
        aplican a una superficie abierta y levantan ``ValueError``).

    Devuelve
    --------
    ``(ny, nx)`` complejo [Pa], con ``out[i, j] = p(x = xs[j], y = ys[i])`` (filas = ``y``).

    Ejemplos
    --------
    El arreglo enfocado del dashboard, en espacio libre::

        from ultrasound_sim.dashboard import field
        fuentes = field.array_sources(64, 0.040, 0.090, 0.0, 250e3, WATER)
        p = free_space_field(xs, ys, fuentes, 250e3, WATER)
    """
    frequency = _check_frequency(frequency)
    material = _check_material(material)
    panel_list = _check_panels(panels)
    xs = _check_axis(xs, "xs")
    ys = _check_axis(ys, "ys")

    grid_x, grid_y = np.meshgrid(xs, ys)                       # (ny, nx)
    pts = np.column_stack([grid_x.ravel(), grid_y.ravel()])
    values = free_space_pressure(
        pts, panel_list, frequency, material,
        n_per_panel=n_per_panel, max_subdivisions=max_subdivisions, equation=equation)
    return values.reshape(ys.size, xs.size)
