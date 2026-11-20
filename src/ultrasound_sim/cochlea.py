"""Fase 5a — corte transversal simplificado de la cóclea (sin membrana basilar).

Geometría (lo que el módulo construye)
--------------------------------------
Un corte 2D en plano ``x-y`` con simetría de revolución alrededor del origen::

    r < a              fluido interior   ("perilinfa"):  materials.WATER
    a <= r <= a + t    cápsula ótica     (sólido elástico):  materials.ElasticSolid
    a + t < r < R_out  fluido exterior   (medio incidente):  materials.WATER

``a`` es el radio interior de la cápsula (cavidad perilinfática), ``t`` su espesor y ``R_out`` el
radio del contorno exterior, donde se impone la ABC de primer orden. La onda plana incidente viaja
en ``+x`` (entra por ``x < 0``), igual que en el caso acoplado de la Fase 2.

Qué es aproximación y qué no
----------------------------
- **Aproximación geométrica**: la cóclea real es un conducto espiral con dos scalae separadas por
  la membrana basilar y las de Reissner. Acá se modela **sólo** la cápsula ótica como corona
  circular y un único compartimento fluido; la membrana basilar y las dos scalae son la Fase 5b.
  La corona es el corte transversal de la cápsula, no su forma espiral.
- **Aproximación de la ABC**: la ABC de primer orden sobre el **campo total** es exacta en
  incidencia normal sobre una cara plana alineada (Fase 1) y tiene un piso de reflexión espuria en
  contornos curvos (D12). Acá el contorno es un círculo concéntrico, así que la onda saliente es casi
  cilíndrica y el error cae con ``1/(k·R_out)``; por eso ``R_out = a + t + 2·λ`` como mínimo. El
  módulo **no** corrige ese error, lo mide: con un dispersor fuerte y ``ka ~ 8–17`` el campo total
  contra la serie exacta de un cilindro rígido da **21–28 % en L2** (medias de amplitud dentro de
  1–2 %), y no mejora al agrandar el dominio (``pad`` 2.5/4/6 λ → 21/24/27 %), o sea que el piso es
  la ABC-1 en contorno curvo más la dispersión P1, no la distancia al contorno.
- **Aproximación de material**: ``ElasticSolid`` isótropo y sin pérdidas; los valores de hueso son
  el ballpark del plan (D5), ``c_s`` **no** está en el plan y va marcado ``VERIFY``.

Transmisión al fluido interior: qué se puede medir y qué no
-----------------------------------------------------------
``t_interior`` (promedio de ``|p|`` sobre ``r = a/2``, normalizado por la incidente) es la métrica
que pide el encargo. **No es la transmisión de una placa** y no converge a ella: la cavidad es
cerrada y sin pérdidas, así que (i) su campo es una onda estacionaria y (ii) la pared recibe energía
por toda la circunferencia, no sólo por la cara iluminada. Medido a 100 kHz con ``t = 2 mm``
(``|T|`` de la placa equivalente = 0.7959): ``a/t = 5, 10, 20`` dan ``-16 %, +29 %, +0.3 %`` y
``a/t = 40`` da ``+19 %`` — sin tendencia a cero, con ±15 % de dispersión entre puntos vecinos (el
nivel se mueve ±30 % en un barrido de frecuencia de ±5 %). Con una cáscara **casi fluida**
(impedancia casada, ``|T| = 1``) el mismo ``t_interior`` cae dentro de **0.5–1.5 %** de la placa, que
es la validación de la cadena completa (malla + acoplamiento + ABC + término incidente).

Para comparar contra una placa hace falta (a) un lado transmitido **abierto** (ducto o semiespacio,
lo que hace la Fase 2) o (b) amortiguar la cavidad. Ninguna de las dos cosas está en la geometría
cerrada que pide la Fase 5a, y por eso el test de límite plano reporta el número medido en vez de
forzar el 5 %.

Ventana redonda (Fase 5b-1): el parche de presión nula
-----------------------------------------------------
``window_half_angle`` abre un **arco de la cara interna de la cápsula** (``r = a``, centrado en
``window_center_angle``, por defecto ``0`` = cara de sombra) donde se impone ``p = 0`` (Dirichlet).
Es la idealización declarada de la **ventana redonda**: membrana flexible, masa e impedancia
despreciables ⇒ nodo de presión. Supuestos, escritos porque no se heredan de ninguna fuente (D20):

- Es un parche de **presión nula**, no una membrana con ``E``, espesor ni pretensión: por D20 no hay
  conjunto 2D validado para eso, y la rigidez por unidad de longitud exigiría declarar el ancho
  fuera del plano. Tampoco es la **platina del estribo** (ventana oval): ésa es entrada, va sobre la
  cara iluminada y su idealización sería una presión impuesta, no un nodo.
- El parche se aplica sobre los nodos del anillo ``r = a`` dentro del arco. El área efectiva queda
  cuantizada por ``n_theta`` (a 20 mm y ``h = 0.75 mm``, ``0.72°`` por sector): con 5° de semiángulo
  el arco tiene 13 nodos a 95 kHz y 15 a 105 kHz si se remalla, o sea que el tamaño del parche salta
  ~15 % entre mallas. Para barrer frecuencia sin ese ruido hay que **reusar la malla**
  (``dataclasses.replace(section, frequency=f)``), que es lo que hace
  ``tests/test_cochlea_window.py``.

**Lo que el parche p=0 NO es: una salida de energía.** Una pared blanda es un reflector perfecto
(``R = −1``): el flujo de Poynting real que la cruza es idénticamente nulo, porque ``p = 0`` en la
superficie. Medido: ``P_ventana/P_inc ≈ 3e−3`` real contra ``−2.1e−2`` reactiva (a 10°, en el anillo
justo por dentro de ``a``; en ``r = a`` la real es exactamente 0) y ``P_cavidad/P_inc ≈ 1e−5`` con el
parche abierto o cerrado — el balance de energía se cierra, pero **todos los términos son cero**. Lo
que el parche hace es *desintonizar* los modos de la cavidad (baja la presión interior y mezcla
armónicos), no amortiguarla: la dispersión de ``t_interior`` en un barrido de ±5 % apenas cambia
(17.8 % cerrada → 12.9 % a 10°, 16.5 % a 5° y 17.2 % a 20°; ver el test de dispersión).

Por eso ``window_model="open"`` (opcional, mismo arco) **saca los elementos sólidos del arco** y deja
que el fluido pase por la ventana hasta el dominio exterior, donde la ABC sí absorbe potencia real.
Es el límite de membrana transparente (impedancia nula, no sólo presión nula) y es el único de los dos
que puede actuar como salida. Medido: con 10° de semiángulo la dispersión baja a 12.9 % y con 20° a
8.6 % (std/media; 11.6 % en semirrango), y la potencia que cruza la ventana crece monótona con el
área (0.21 % → 1.12 % → 8.84 % de ``P_inc``). Aun así, con un arco fisiológicamente razonable el
criterio de 10 % no se alcanza: la cavidad sigue siendo un resonador de Q alto y el amortiguador de
banda ancha de la cóclea real es la **membrana basilar** (5b-2), no la ventana. Los dos modelos quedan
porque el negativo es el dato: dice qué falta.

Ventana redonda (Fase 5b-1b): la impedancia resistiva
----------------------------------------------------
``window_model='impedance'`` reemplaza el parche por la condición **física** de la ventana,
``p = Z_r·(v·n)`` sobre el mismo arco (``r = a``). Con ``v_n = (1/(iωρ))·∂p/∂n`` sale
``∂p/∂n = (i·ω·ρ/Z_r)·p``, o sea **el mismo término de borde que la ABC-1 con ``k → k_w``**,
``k_w = ω·ρ_f/Z_r``: la forma débil agrega ``−∫ (i·k_w/ρ_f)·p·q ds`` sobre las aristas del arco. No
hace falta un modelo de membrana: ``Z_r`` es un dato de entrada y la condición es exactamente la de
un contorno **resistivo** (absorbente), que sí transporta potencia real
``P = ½·∫ |p|²·Re(1/Z_r) ds > 0`` hacia afuera (con ``Re Z_r > 0``) — al revés que ``p = 0``, donde
``P = p·(v·n) = 0`` por construcción (D22).

Dos identidades que no son aproximaciones y por eso sirven de test:

- con ``Z_r = ρ_f·c_f`` el coeficiente ``k_w = ω/c_f`` **es** ``k_f``: la condición coincide
  término a término con la ABC-1 aplicada sobre ese arco, así que el resultado tiene que reproducir
  el caso ABC-1 sobre el mismo arco (y caer a máquina sobre la matriz, ver
  ``tests/test_cochlea_impedance.py``);
- con ``Z_r`` real y finito la potencia real por el arco **no es nula** y **sale** (``P > 0``),
  medible de dos maneras independientes: el flujo por diferencia radial un anillo por dentro
  (:func:`arc_power_flux`, lo mismo que se mide en los otros modelos) y la integral de borde
  ``½∫|p|²·Re(1/Z_r) ds`` (:func:`window_boundary_power`). Además ``P_inc − P_disp = P_ventana``,
  que es el balance de energía del caso.

``Z_r`` es un **supuesto declarado** (D20 no tiene impedancia citable de la ventana redonda): el
valor por defecto es ``ρ_f·c_f`` (agua casada, el que hace la identidad con la ABC) y la sensibilidad
se barre 0.5 / 1 / 2 × ρc en vez de fijar un número con aire de medido.

*Implementación*: ``solve_coupled`` no admite términos de borde en aristas **interiores** (el arco
``r = a`` es la interfaz fluido/sólido, no frontera de la malla), así que el caso ``'impedance'``
ensambla el mismo sistema acoplado acá, con los bloques públicos de ``fem``
(``p1_elastic_local``, ``p1_tri_local``, ``abc_boundary_block``, ``abc_incident_rhs``) más el
término del arco. La comprobación de que el ensamblado es el mismo: con el término del arco apagado
tiene que reproducir ``solve_cross_section`` de la Fase 5a dígito a dígito.

Convención de tiempo y atenuación: D4 (``e^{-iωt}``), igual que todo el repo. Todos los materiales
que entran son de ``materials.py`` o construidos explícitamente por quien llama, con ``provenance``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from . import meshing
from .analytic.helmholtz1d import StackResponse, stack_response
from .fem.assembly import p1_tri_local, weights_to_local
from .fem.coupled2d import (
    CoupledDofs,
    CoupledSolution,
    abc_boundary_block,
    abc_incident_rhs,
    interface_edges,
    solve_coupled,
)
from .fem.elastic2d import p1_elastic_local
from .fem.helmholtz2d import PlaneWave, plane_wave
from .materials import WATER, ElasticSolid, Material

__all__ = [
    "BONE_BALLPARK",
    "CochleaCrossSection",
    "CrossSectionSolution",
    "arc_power_flux",
    "build_cross_section",
    "intercepted_power",
    "layer_from_solid",
    "plate_reference",
    "plate_speed",
    "radial_power_flux",
    "scattered_power",
    "solve_arc_condition",
    "solve_cross_section",
    "wall_injection_power",
    "window_arc_edges",
    "window_boundary_power",
]

_TWO_PI = 2.0 * np.pi

#: Fracciones de ``a`` donde la malla fuerza anillos de nodos en la cavidad. ``a/2`` es la
#: circunferencia de la métrica de presión media; los demás anillos existen para que el flujo radial
#: (:func:`radial_power_flux`) tenga anillos consecutivos en cualquier radio de la cavidad.
_PROBE_FRACTIONS = (0.20, 0.30, 0.40, 0.50, 0.60, 0.75, 0.90)

#: Hueso por defecto. ``c_p`` y ``ρ`` son el punto medio del rango del plan; ``c_s`` **no** está en
#: el plan (D5): queda marcado y no debe usarse en una validación contra literatura.
BONE_BALLPARK = ElasticSolid(
    name="bone_cortical_ballpark",
    rho=1950.0,
    c_p=3500.0,
    c_s=1900.0,
    provenance=(
        "VERIFY — c_p=3500 y ρ=1950 son el punto medio del rango del plan "
        "(3000–4000 m/s, 1900–2000 kg/m³, 'ballpark, verify'); c_s=1900 m/s no figura en el plan: "
        "es un ballpark de hueso cortical, no verificado contra literatura (D5)."
    ),
)


# ---------------------------------------------------------------------------
# Geometría
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CochleaCrossSection:
    """Malla y metadatos del corte transversal (cápsula elástica + fluido interior).

    Atributos
    ---------
    nodes : ``(N, 2)`` coordenadas [m].
    triangles : ``(M, 3)`` triángulos CCW.
    regions : ``(M,)`` ``'fluid'`` o ``'solid'`` por elemento (lo que consume
        :func:`fem.coupled2d.solve_coupled`).
    abc_edges : ``(B, 2)`` aristas de la circunferencia ``r = R_out``, orientadas por el triángulo
        CCW (normal saliente ``+r``).
    frequency : frecuencia [Hz].
    a, thickness, outer_radius : geometría [m].
    ppw : puntos por longitud de onda pedidos (en el fluido).
    h : paso efectivo de malla [m] (el mínimo entre ``λ/ppw`` y ``t/min_elements_across_wall``).
    n_theta : sectores angulares (par y múltiplo de 4: los ejes son rayos de nodos).
    radii : radios de los anillos de nodos (sin el centro), crecientes.
    ring_start : índice del primer nodo de cada anillo en ``nodes`` (``ring_start[k]`` es el nodo
        ``j = 0`` del anillo ``k``).
    band_region : región de la banda radial entre anillos consecutivos.
    probe_radius : ``a/2``, radio de la circunferencia de la métrica ``t_interior``.
    window_half_angle : semiángulo [**grados**] del arco de ventana redonda, o ``None`` (cavidad
        cerrada, comportamiento de la Fase 5a).
    window_center_angle : centro del arco [**grados**]; ``0`` = ``+x``, la cara de **sombra** (la onda
        incidente viaja en ``+x``).
    window_model : ``'pressure_release'`` (parche ``p = 0``, idealización declarada), ``'open'``
        (se quitan los elementos sólidos del arco: el fluido pasa y la ABC exterior absorbe) o
        ``'impedance'`` (``p = Z_r·(v·n)`` sobre el arco: contorno resistivo, salida de energía real).
    window_impedance : ``Z_r`` [Pa·s/m] del modelo ``'impedance'``; complejo (``Re Z_r > 0`` para
        que absorba). ``None`` ⇒ ``ρ_f·c_f`` del fluido, que es el valor que hace la condición
        idéntica a la ABC-1 sobre el arco. Ver el docstring del módulo (D20/D22).
    """

    nodes: np.ndarray
    triangles: np.ndarray
    regions: np.ndarray
    abc_edges: np.ndarray
    frequency: float
    a: float
    thickness: float
    outer_radius: float
    ppw: float
    h: float
    n_theta: int
    radii: np.ndarray
    ring_start: np.ndarray
    band_region: np.ndarray
    probe_radius: float
    fluid: Material
    bone: ElasticSolid
    window_half_angle: float | None = None
    window_center_angle: float = 0.0
    window_model: str = "pressure_release"
    window_impedance: complex | None = None

    # -- consultas ----------------------------------------------------------

    @property
    def has_window(self) -> bool:
        return self.window_half_angle is not None and self.window_half_angle > 0.0

    def window_nodes(self) -> np.ndarray:
        """Nodos del anillo ``r = a`` dentro del arco de la ventana (vacío si no hay ventana).

        Son los nodos donde :func:`solve_cross_section` impone ``p = 0`` en el modelo
        ``'pressure_release'``. El arco se mide en el plano ``x-y`` desde ``window_center_angle``
        (grados) y su área efectiva está cuantizada por ``n_theta``.
        """
        if not self.has_window:
            return np.zeros(0, dtype=int)
        return _ring_arc(self, self.a, self.window_half_angle, self.window_center_angle)

    def window_edges(self, half_angle: float | None = None, center_angle: float | None = None):
        """Aristas de la **interfaz** ``r = a`` dentro del arco (las que usa el modelo de impedancia).

        Son aristas interiores fluido/sólido (el arco es la cara interna de la cápsula), así que no
        son frontera de la malla y por eso el término ``−∫(i·k_w/ρ)·p·q ds`` se ensambla en
        :func:`solve_arc_condition`, no en ``solve_coupled``.
        """
        return window_arc_edges(self, half_angle, center_angle)


    @property
    def n_nodes(self) -> int:
        return int(self.nodes.shape[0])

    @property
    def n_elements(self) -> int:
        return int(self.triangles.shape[0])

    @property
    def ppw_effective(self) -> float:
        return float(self.fluid.wavelength(self.frequency) / self.h)

    @property
    def outer_wall_radius(self) -> float:
        return self.a + self.thickness

    def describe(self) -> str:
        return (
            f"corte: a={self.a * 1e3:.2f} mm, t={self.thickness * 1e3:.2f} mm, "
            f"R_out={self.outer_radius * 1e3:.2f} mm, a/t={self.a / self.thickness:.1f}, "
            f"f={self.frequency / 1e3:.1f} kHz, ppw={self.ppw_effective:.1f} "
            f"(h={self.h * 1e3:.3f} mm), {self.n_nodes} nodos, {self.n_elements} elementos "
            f"({int(np.count_nonzero(self.regions == 'solid'))} sólidos), n_theta={self.n_theta}"
        )

    def ring_nodes(self, radius: float) -> np.ndarray:
        """Índices de los nodos del anillo de radio ``radius`` (exacto en la malla)."""
        r = np.linalg.norm(self.nodes, axis=1)
        sel = np.flatnonzero(np.isclose(r, float(radius), rtol=0.0, atol=1e-9 * self.h))
        if sel.size == 0:
            raise ValueError(
                f"no hay anillo de nodos en r={radius!r}; los anillos son "
                f"{np.array2string(self.radii, precision=6)}"
            )
        return sel

    def ray_nodes(self, theta: float) -> np.ndarray:
        """Nodos sobre el rayo ``theta`` (ordenados por radio; el centro va primero)."""
        x, y = self.nodes[:, 0], self.nodes[:, 1]
        r = np.hypot(x, y)
        ang = np.arctan2(y, x)
        tol = 0.25 * _TWO_PI / self.n_theta
        d = np.abs(np.arctan2(np.sin(ang - theta), np.cos(ang - theta)))
        sel = np.flatnonzero(d < tol)
        return sel[np.argsort(r[sel])]


def _ring_radii(a, thickness, outer_radius, h, min_across_wall, fractions):
    """Radios de los anillos de nodos, con anillos exactos en ``a/2``, ``a`` y ``a+t``.

    Los anillos de la cavidad incluyen las fracciones fijas ``fractions`` de ``a`` (sondas que
    existen en la malla, sin interpolar) y se subdividen con paso ``<= h``.
    """
    edges = [0.0] + [f * a for f in sorted(set(fractions) | {0.5}) if 0.0 < f < 1.0] + [a]
    out = []
    for r0, r1 in zip(edges[:-1], edges[1:]):
        n = max(1, int(np.ceil((r1 - r0) / h)))
        out.extend(np.linspace(r0, r1, n + 1)[1:].tolist())
    n_shell = max(int(min_across_wall), int(np.ceil(thickness / h)))
    out.extend(np.linspace(a, a + thickness, n_shell + 1)[1:].tolist())
    n_par = max(1, int(np.ceil((outer_radius - a - thickness) / h)))
    out.extend(np.linspace(a + thickness, outer_radius, n_par + 1)[1:].tolist())
    return np.asarray(out, dtype=float)


def _ring_arc(section, radius, half_angle_deg, center_deg=0.0) -> np.ndarray:
    """Nodos del anillo ``r = radius`` con ``|θ − center| <= half_angle`` (ángulos en grados).

    El arco se mide en el plano ``x-y`` con ``θ = atan2(y, x)`` y se cierra por el camino corto
    (``atan2(sin, cos)``), así que ``center_deg = 180`` funciona igual que ``0``.
    """
    ring = section.ring_nodes(radius)
    ang = np.arctan2(section.nodes[ring, 1], section.nodes[ring, 0])
    c = np.radians(float(center_deg))
    d = np.arctan2(np.sin(ang - c), np.cos(ang - c))
    return ring[np.abs(d) <= np.radians(float(half_angle_deg)) + 1e-12]


def _validate_window(half_angle, center_angle, model):
    """Normaliza y valida los parámetros de la ventana redonda."""
    model = str(model)
    if model not in ("pressure_release", "open", "impedance"):
        raise ValueError(
            "window_model debe ser 'pressure_release', 'open' o 'impedance'; "
            f"recibido {model!r}"
        )
    center_angle = float(center_angle)
    if not np.isfinite(center_angle):
        raise ValueError(f"window_center_angle debe ser finito; recibido {center_angle!r}")
    if half_angle is None:
        return None, center_angle, model
    half_angle = float(half_angle)
    if not np.isfinite(half_angle) or not (0.0 < half_angle <= 90.0):
        raise ValueError(
            "window_half_angle debe estar en (0, 90] grados (0 = sin ventana, use None; > 90 "
            f"invadiría la cara iluminada); recibido {half_angle!r}"
        )
    return half_angle, center_angle, model


def _validate_impedance(impedance, fluid: Material) -> complex:
    """``Z_r`` normalizado: ``None`` ⇒ ``ρ_f·c_f`` (agua casada); ``Re Z_r > 0`` obligatorio.

    El caso ``Z_r = ρ_f·c_f`` es el que hace ``k_w = ω·ρ_f/Z_r = k_f``, o sea la condición
    **idéntica** a la ABC-1 sobre el arco. Por eso es el default: el único valor que no inventa
    nada nuevo y que además sirve de test interno. ``Re Z_r > 0`` es lo que garantiza que el arco
    **absorba** (``P = ½∫|p|²·Re(1/Z_r) ds > 0``); una parte reactiva pura no transporta potencia
    real y volvería al caso ``p = 0``/pared blanda, así que se rechaza en vez de devolver un cero
    silencioso.
    """
    z = complex(fluid.impedance if impedance is None else complex(impedance))
    if not (np.isfinite(z.real) and np.isfinite(z.imag)):
        raise ValueError(f"window_impedance debe ser finito; recibido {z!r}")
    if z.real <= 0.0:
        raise ValueError(
            "window_impedance necesita Re(Z_r) > 0 para que la ventana transporte potencia real "
            f"hacia afuera; recibido {z!r} (Re = {z.real:g})"
        )
    return z


def _window_hole_mask(nodes, triangles, a, thickness, half_angle, center_angle) -> np.ndarray:
    """Elementos de la banda ``[a, a+t]`` cuyo centroide cae dentro del arco de la ventana."""
    cent = nodes[triangles].mean(axis=1)
    r = np.linalg.norm(cent, axis=1)
    ang = np.arctan2(cent[:, 1], cent[:, 0])
    c = np.radians(float(center_angle))
    d = np.arctan2(np.sin(ang - c), np.cos(ang - c))
    in_wall = (r > a) & (r < a + thickness)
    return in_wall & (np.abs(d) <= np.radians(float(half_angle)))


def build_cross_section(
    frequency: float,
    a: float,
    thickness: float,
    outer_radius: float | None = None,
    ppw: float = 20.0,
    fluid: Material = WATER,
    bone: ElasticSolid = BONE_BALLPARK,
    min_elements_across_wall: int = 3,
    pad_wavelengths: float = 2.0,
    *,
    window_half_angle: float | None = None,
    window_center_angle: float = 0.0,
    window_model: str = "pressure_release",
    window_impedance: complex | None = None,
) -> CochleaCrossSection:
    """Malla del corte transversal: fluido interior, corona elástica y fluido exterior.

    Parameters
    ----------
    frequency : frecuencia [Hz] (``> 0``).
    a : radio interior de la cápsula [m] (pared interior de la perilinfa).
    thickness : espesor de la cápsula [m].
    outer_radius : radio del contorno exterior [m]. Si es ``None``,
        ``a + thickness + pad_wavelengths·λ_fluido``.
    ppw : puntos por longitud de onda en el fluido; el paso se limita además a
        ``thickness / min_elements_across_wall`` para no perder el espesor de la pared. **Ojo**: eso
        hace que subir ``ppw`` sólo refine mientras ``λ/ppw < t/min_elements_across_wall``; para
        refinar más hay que subir ``min_elements_across_wall`` (o bajar ``t``).
    fluid, bone : materiales (``materials.WATER`` y :data:`BONE_BALLPARK` por defecto).
    min_elements_across_wall : elementos radiales mínimos a través de la pared.
    pad_wavelengths : margen exterior en longitudes de onda cuando ``outer_radius`` es ``None``.
    window_half_angle : semiángulo [**grados**] del arco de la ventana redonda sobre la cara interna
        de la cápsula (``r = a``); ``None`` (default) = cavidad cerrada, el problema de la Fase 5a.
    window_center_angle : centro del arco [grados]; ``0`` es ``+x``, la cara de sombra.
    window_model : ``'pressure_release'`` (default) deja la malla intacta y el parche ``p = 0`` se
        impone en :func:`solve_cross_section`; ``'open'`` **quita** los elementos sólidos del arco
        (el fluido pasa por la ventana hasta el dominio exterior); ``'impedance'`` deja la malla
        intacta y el arco pasa a ser un contorno resistivo ``p = Z_r·(v·n)``. Ver el docstring del
        módulo.
    window_impedance : ``Z_r`` [Pa·s/m] del modelo ``'impedance'``; si es ``None`` se guarda el
        default ``ρ_f·c_f`` (agua casada) resuelto con el fluido de la sección.

    Returns
    -------
    :class:`CochleaCrossSection`.

    Notes
    -----
    El mallado es polar estructurado: ``n_theta`` sectores iguales y anillos concéntricos, con
    anillos forzados en ``r = a/2`` (circunferencia de la métrica), ``r = a`` y ``r = a + t``
    (las dos caras de la cápsula). Así la interfaz fluido/sólido es exactamente conforme y todos
    los nodos de los anillos usados como sondas existen en la malla, sin interpolación.

    Con ``window_model='open'`` los triángulos cuyo centroide cae en la banda ``[a, a+t]`` y dentro
    del arco pasan a ``'fluid'``: la cápsula queda con una **abertura** de lados rectos (paredes
    rígidas, ``∂p/∂n = 0``), que es el límite de membrana transparente. La malla no cambia de
    topología y sigue siendo conforme.
    """
    frequency = float(frequency)
    a = float(a)
    thickness = float(thickness)
    window_half_angle, window_center_angle, window_model = _validate_window(
        window_half_angle, window_center_angle, window_model
    )
    if not frequency > 0.0:
        raise ValueError(f"frequency debe ser > 0; recibido {frequency!r}")
    if not a > 0.0:
        raise ValueError(f"a debe ser > 0; recibido {a!r}")
    if not thickness > 0.0:
        raise ValueError(f"thickness debe ser > 0; recibido {thickness!r}")
    if not isinstance(fluid, Material):
        raise TypeError(f"fluid debe ser materials.Material; recibido {type(fluid).__name__}")
    if not isinstance(bone, ElasticSolid):
        raise TypeError(
            f"bone debe ser materials.ElasticSolid; recibido {type(bone).__name__}"
        )
    window_impedance = (
        _validate_impedance(window_impedance, fluid) if window_model == "impedance" else None
    )
    if ppw < 6.0:
        raise ValueError(f"ppw debe ser >= 6; recibido {ppw!r}")
    if min_elements_across_wall < 1:
        raise ValueError("min_elements_across_wall debe ser >= 1")

    lam = fluid.wavelength(frequency)
    if outer_radius is None:
        outer_radius = a + thickness + float(pad_wavelengths) * lam
    outer_radius = float(outer_radius)
    if not outer_radius > a + thickness:
        raise ValueError(
            f"outer_radius ({outer_radius!r}) debe ser > a + thickness ({a + thickness!r})"
        )

    h = min(lam / float(ppw), thickness / float(min_elements_across_wall))
    radii = _ring_radii(
        a, thickness, outer_radius, h, min_elements_across_wall, _PROBE_FRACTIONS
    )
    radii = np.unique(np.round(radii, 12))

    n_theta = int(np.ceil(_TWO_PI * outer_radius / h))
    n_theta = max(24, n_theta + (n_theta % 4))

    theta = _TWO_PI * np.arange(n_theta) / n_theta
    ct, st = np.cos(theta), np.sin(theta)

    # nodos: 0 = centro; después un bloque (n_theta, 2) por anillo
    coords = [np.zeros((1, 2))]
    ring_start = np.empty(radii.size, dtype=int)
    nid = 1
    for k, r in enumerate(radii):
        ring_start[k] = nid
        coords.append(np.column_stack([r * ct, r * st]))
        nid += n_theta
    nodes = np.vstack(coords)

    def node(k: int, j: int) -> int:
        """Nodo ``j`` del anillo ``k`` (``k = -1`` ⇒ centro)."""
        if k < 0:
            return 0
        return int(ring_start[k] + (j % n_theta))

    tris: list[tuple[int, int, int]] = []
    band_region: list[str] = []

    # abanico central
    mid_r = 0.5 * radii[0]
    for j in range(n_theta):
        tris.append((0, node(0, j), node(0, j + 1)))
        band_region.append("solid" if a < mid_r < a + thickness else "fluid")

    for k in range(radii.size - 1):
        mid = 0.5 * (radii[k] + radii[k + 1])
        label = "solid" if a < mid < a + thickness else "fluid"
        for j in range(n_theta):
            p, q = node(k, j), node(k, j + 1)
            s, t = node(k + 1, j), node(k + 1, j + 1)
            tris.append((p, s, t))
            tris.append((p, t, q))
            band_region.extend([label, label])

    triangles = np.asarray(tris, dtype=int)
    neg = meshing.triangle_areas(nodes, triangles) < 0.0
    if np.any(neg):
        triangles[neg] = triangles[neg][:, [0, 2, 1]]

    regions = np.asarray([str(v) for v in band_region])
    # control: la banda radial y el centroide de cada elemento tienen que coincidir
    cent_r = np.linalg.norm(nodes[triangles].mean(axis=1), axis=1)
    solid_geom = (cent_r > a) & (cent_r < a + thickness)
    hole = np.zeros(triangles.shape[0], dtype=bool)
    if window_model == "open" and window_half_angle is not None:
        hole = _window_hole_mask(
            nodes, triangles, a, thickness, window_half_angle, window_center_angle
        )
        regions = regions.copy()
        regions[hole] = "fluid"
    if not np.array_equal(regions == "solid", solid_geom & ~hole):
        bad = int(np.flatnonzero((regions == "solid") != (solid_geom & ~hole))[0])
        raise RuntimeError(
            f"marcación de regiones inconsistente en el elemento {bad}: la banda dice "
            f"{regions[bad]!r} y el centroide está en r={cent_r[bad]:.6g}"
        )

    edges = meshing.boundary_edges(triangles)
    r_edges = np.linalg.norm(nodes[edges].mean(axis=1), axis=1)
    abc = edges[np.isclose(r_edges, outer_radius, rtol=0.0, atol=0.05 * h)]
    if abc.shape[0] != n_theta:
        raise RuntimeError(
            f"la frontera exterior tiene {abc.shape[0]} aristas y se esperaban {n_theta}"
        )

    return CochleaCrossSection(
        nodes=nodes,
        triangles=triangles,
        regions=regions,
        abc_edges=abc,
        frequency=frequency,
        a=a,
        thickness=thickness,
        outer_radius=outer_radius,
        ppw=float(ppw),
        h=float(h),
        n_theta=n_theta,
        radii=radii,
        ring_start=ring_start,
        band_region=np.asarray(band_region),
        probe_radius=0.5 * a,
        fluid=fluid,
        bone=bone,
        window_half_angle=window_half_angle,
        window_center_angle=window_center_angle,
        window_model=window_model,
        window_impedance=window_impedance,
    )


# ---------------------------------------------------------------------------
# Referencia analítica: placa plana del mismo espesor
# ---------------------------------------------------------------------------


def layer_from_solid(solid: ElasticSolid) -> Material:
    """Capa fluida equivalente a un sólido en incidencia normal (``Z = ρ·c_p``).

    Es la convención con la que el plan (Fase 2) y `tests/test_coupled_fluid_solid.py` comparan
    una capa elástica contra `stack_response`: en deformación uniaxial la onda longitudinal viaja a
    ``c_p``. Una placa **libre** (sin rodillo) en plano deformación transmite a la velocidad de
    placa, algo menor; la diferencia se mide y se reporta en el test de límite plano.
    """
    return Material(
        name=f"{solid.name}_capa_equivalente",
        c=solid.c_p,
        rho=solid.rho,
        provenance="capa elástica equivalente en incidencia normal (Z = ρ·c_p), Fase 2",
    )


def plate_speed(solid: ElasticSolid) -> float:
    """Velocidad de onda longitudinal de **placa libre** en plano deformación.

    ``c_placa = sqrt(4·μ·(λ+μ) / (ρ·(λ+2μ)))`` con ``λ, μ`` los módulos de Lamé del sólido. Es la
    velocidad de la onda extensional de una placa delgada con caras libres (plano deformación);
    ``c_p = sqrt((λ+2μ)/ρ)`` es la de deformación uniaxial (la que impone un rodillo).
    """
    lam, mu = solid.lam, solid.mu
    return float(np.sqrt(4.0 * mu * (lam + mu) / (solid.rho * (lam + 2.0 * mu))))


def plate_reference(
    frequency: float,
    thickness: float,
    fluid: Material = WATER,
    bone: ElasticSolid = BONE_BALLPARK,
    *,
    speed: float | None = None,
) -> StackResponse:
    """``R``/``T`` analíticos de una placa plana del mismo espesor (fluido / placa / fluido).

    ``speed=None`` usa ``c_p`` del sólido (convención de la Fase 2, la que pide el plan);
    ``speed=plate_speed(bone)`` usa la velocidad de placa libre.
    """
    if speed is None:
        layer = layer_from_solid(bone)
    else:
        layer = Material(
            name=f"{bone.name}_velocidad_placa",
            c=float(speed),
            rho=bone.rho,
            provenance="onda extensional de placa libre en plano deformación",
        )
    return stack_response(frequency, [fluid, layer, fluid], [float(thickness)])


# ---------------------------------------------------------------------------
# Medición sobre el campo
# ---------------------------------------------------------------------------


def intercepted_power(section: CochleaCrossSection, amplitude: complex = 1.0) -> float:
    """Potencia [W/m] que la onda plana intercepta en la silueta de la cápsula.

    ``I·2·(a+t)`` con ``I = |p_inc|²/(2·ρ·c)``. Es la normalización de una transmisión "por unidad
    de ancho proyectado": en el límite de placa plana, la potencia transmitida dividida por este
    valor tiende a ``|T|²``.
    """
    amp = abs(complex(amplitude))
    return float(
        amp**2 * 2.0 * (section.a + section.thickness) / (2.0 * section.fluid.rho * section.fluid.c)
    )


def _flux_from_rings(section, values, radii) -> float:
    """Flujo radial neto [W/m] a partir de los valores en 2 o 3 anillos consecutivos.

    Con tres anillos usa diferencia **central** (error ``O(h²)`` y factor 4 menor que la lateral);
    con dos, lateral. ``I_r = Im(p*·∂_r p)/(2·ω·ρ)``, ``P = 2π·r·⟨I_r⟩`` y Parseval para la media
    angular.
    """
    vals = [np.asarray(v, dtype=complex) for v in values]
    radii = [float(r) for r in radii]
    n_theta = section.n_theta
    spectra = [np.fft.fft(v) / n_theta for v in vals]
    for x in spectra:
        if not np.all(np.isfinite(x)):
            raise ValueError("los anillos del flujo tienen nodos sin gdl de presión")
    if len(radii) == 3:
        r0, r1, r2 = radii
        # diferencia central no uniforme (segundo orden)
        h1, h2 = r1 - r0, r2 - r1
        a = -h2 / (h1 * (h1 + h2))
        b = (h2 - h1) / (h1 * h2)
        c = h1 / (h2 * (h1 + h2))
        dx = a * spectra[0] + b * spectra[1] + c * spectra[2]
        integral = float(np.sum(np.imag(np.conj(spectra[1]) * dx)))
        return float(np.pi * r1 * integral / (_TWO_PI * section.frequency * section.fluid.rho))
    r1, r2 = radii
    dx = (spectra[1] - spectra[0]) / (r2 - r1)
    r_mid = 0.5 * (r1 + r2)
    integral = float(np.sum(np.imag(np.conj(spectra[0]) * dx)))
    return float(np.pi * r_mid * integral / (_TWO_PI * section.frequency * section.fluid.rho))


def _rings_around(section, radius) -> tuple[list[float], list[int]]:
    """Anillos de la malla alrededor de ``radius``: 3 (central) o 2 si está en el borde."""
    radii = np.asarray(section.radii, dtype=float)
    i = int(np.argmin(np.abs(radii - float(radius))))
    if 0 < i < radii.size - 1:
        sel = [i - 1, i, i + 1]
    elif i == 0:
        sel = [0, 1]
    else:
        sel = [radii.size - 2, radii.size - 1]
    rs = [float(radii[j]) for j in sel]
    nodes = [section.ring_nodes(r) for r in rs]
    return rs, nodes


def radial_power_flux(section: CochleaCrossSection, p: np.ndarray, radius: float) -> float:
    """Flujo neto de potencia radial [W/m] a través del círculo ``r = radius`` (positivo afuera).

    Cuadratura exacta sobre los nodos del anillo y ``∂_r p`` por diferencia finita radial (central
    si hay tres anillos), con la convención ``e^{−iωt}``: ``I_r = Im(p*·∂_r p)/(2·ω·ρ)``,
    ``P = 2π·r·⟨I_r⟩``, y ``⟨p*·q⟩ = Σ_m conj(X_m)·Y_m`` por Parseval.
    """
    rs, nodes = _rings_around(section, radius)
    return _flux_from_rings(section, [p[n] for n in nodes], rs)


def _cavity_rings(section, radius=None, count=3):
    """Anillos **estrictamente dentro** de la cavidad (``r < a``) y su índice posicional.

    El anillo ``r = a`` es la interfaz: sus vecinos radiales exteriores son nodos sólo sólidos, sin
    gdl de presión, así que un flujo derivado con diferencia radial no puede usar ``r = a`` como
    anillo central. Se evalúa en el anillo de malla más cercano por dentro.
    """
    radii = np.asarray(section.radii, dtype=float)
    inside = radii[radii < section.a - 1e-12]
    if inside.size < 2:
        raise ValueError(
            "la cavidad no tiene dos anillos de nodos por dentro de r = a; no se puede derivar el "
            f"flujo radial (radios: {np.array2string(radii, precision=6)})"
        )
    count = int(min(max(int(count), 2), inside.size))
    if radius is None:
        i = max(1, inside.size - 2)  # último trío: se evalúa en inside[-2]
    else:
        i = int(np.argmin(np.abs(inside - float(radius))))
    sel = [j for j in range(i - 1, i + 2) if 0 <= j < inside.size][:count]
    rings = [np.asarray(section.ring_nodes(r), dtype=int) for r in inside[sel]]
    return [float(inside[j]) for j in sel], rings


def _sector_flux(section, values, radii, mask):
    """``(real, reactiva)`` del flujo radial por el sector ``mask`` de los anillos dados.

    Cuadratura directa (no Parseval) sobre los nodos del sector: ``P = Σ_j I_r,j · ds`` con
    ``ds = 2π·r/n_theta``, ``I_r = Im(p*∂_r p)/(2ωρ)`` y el mismo esquema radial que
    :func:`_flux_from_rings`. Positivo = hacia afuera. La parte reactiva es ``Re(p*∂_r p)/(2ωρ)``:
    no transporta potencia neta, es el intercambio oscilante (campo cercano / pared blanda).
    """
    vals = [np.asarray(v, dtype=complex)[mask] for v in values]
    n_theta = section.n_theta
    omega_rho = _TWO_PI * section.frequency * section.fluid.rho
    if len(radii) == 3:
        r0, r1, r2 = radii
        h1, h2 = r1 - r0, r2 - r1
        a, b, c = -h2 / (h1 * (h1 + h2)), (h2 - h1) / (h1 * h2), h1 / (h2 * (h1 + h2))
        dx = a * vals[0] + b * vals[1] + c * vals[2]
        r_eval, ref = r1, vals[1]
    else:
        dx = (vals[1] - vals[0]) / (radii[1] - radii[0])
        r_eval, ref = 0.5 * (radii[0] + radii[1]), vals[0]
    ds = _TWO_PI * r_eval / n_theta
    prod = np.conj(ref) * dx
    real = float(np.sum(np.imag(prod)) * ds / (2.0 * omega_rho))
    reactive = float(np.sum(np.real(prod)) * ds / (2.0 * omega_rho))
    return real, reactive


def arc_power_flux(
    section: CochleaCrossSection,
    p: np.ndarray,
    half_angle: float,
    *,
    center_angle: float = 0.0,
    radius: float | None = None,
) -> tuple[float, float]:
    """``(real, reactiva)`` [W/m] del flujo radial por un **arco** de la cavidad.

    Se evalúa sobre los anillos de malla por dentro de ``r = a`` (por defecto el último trío) y el
    sector ``|θ − center| <= half_angle`` (grados). Con ``half_angle = 180`` y ``radius`` igual al
    anillo central de :func:`radial_power_flux` reproduce su valor: es la validación de la
    cuadratura directa contra la de Parseval.

    Positivo = potencia saliendo de la cavidad por ese arco. **Ojo con la interpretación cerca de
    una pared blanda**: en ``r = a`` exacto la parte real se anula porque ``p = 0`` en el parche; el
    valor a un anillo por dentro incluye circulación reactiva del campo cercano que no sale del
    dominio (ver el docstring del módulo).
    """
    radii, rings = _cavity_rings(section, radius)
    mid = rings[1] if len(rings) == 3 else rings[0]
    ang = np.arctan2(section.nodes[mid, 1], section.nodes[mid, 0])
    c = np.radians(float(center_angle))
    d = np.arctan2(np.sin(ang - c), np.cos(ang - c))
    mask = np.abs(d) <= np.radians(float(half_angle)) + 1e-12
    if not np.any(mask):
        return 0.0, 0.0
    return _sector_flux(section, [p[r] for r in rings], radii, mask)


def window_arc_edges(
    section: CochleaCrossSection,
    half_angle: float | None = None,
    center_angle: float | None = None,
) -> np.ndarray:
    """Aristas de la interfaz ``r = a`` dentro del arco de la ventana (``(E, 2)``, vacío si no hay).

    El arco es la **cara interna de la cápsula**: ahí las aristas son interiores (un elemento de
    fluido de un lado y uno sólido del otro), no frontera de la malla. Se toman las aristas cuyos
    **dos** nodos caen dentro de ``|θ − center| <= half_angle``, que es el mismo conjunto que usa
    ``p = 0`` para el parche (:meth:`CochleaCrossSection.window_nodes`), así los dos modelos cubren
    el mismo arco y los números son comparables.

    El largo del arco queda cuantizado por ``n_theta`` (a 20 mm y ``h = 0.75 mm``, ``0.72°`` por
    sector), igual que en ``pressure_release``.
    """
    half = section.window_half_angle if half_angle is None else half_angle
    center = section.window_center_angle if center_angle is None else center_angle
    if half is None or float(half) <= 0.0:
        return np.zeros((0, 2), dtype=int)
    x = section.nodes
    iface, _ = interface_edges(x, section.triangles, section.regions)
    if iface.shape[0] == 0:
        return np.zeros((0, 2), dtype=int)
    mid = 0.5 * (x[iface[:, 0]] + x[iface[:, 1]])
    on_inner = np.abs(np.linalg.norm(mid, axis=1) - section.a) <= 0.25 * section.h
    c = np.radians(float(center))
    tol = np.radians(float(half)) + 1e-12

    def _inside(idx):
        ang = np.arctan2(x[idx, 1], x[idx, 0])
        d = np.arctan2(np.sin(ang - c), np.cos(ang - c))
        return np.abs(d) <= tol

    keep = on_inner & _inside(iface[:, 0]) & _inside(iface[:, 1])
    return np.asarray(iface[keep], dtype=int).reshape(-1, 2)


def window_boundary_power(
    section: CochleaCrossSection,
    p: np.ndarray,
    *,
    impedance: complex | None = None,
    half_angle: float | None = None,
    center_angle: float | None = None,
) -> float:
    """Potencia real [W/m] que cruza **el arco**, integrada sobre la superficie (``½∫|p|²·Re(1/Z_r) ds``).

    Es la medida **exacta** del modelo de impedancia, no una diferencia finita: con ``p = Z_r·v_n``
    sobre el arco, ``P = ½·∫ Re(p·v_n*) ds = ½·∫ |p|²·Re(1/Z_r) ds``. Positivo = **saliendo** de la
    cavidad (la normal saliente del fluido en ``r = a`` es ``+r``). Con ``Re Z_r > 0`` es positivo
    por construcción, y es justamente lo que ``p = 0`` no puede dar (allí ``p = 0`` ⇒ ``P = 0``,
    D22).

    La integral se hace con la masa P1 de arista (``L/6·[[2,1],[1,2]]``), así que el valor es
    cuadrático en ``p`` y no depende de diferencias radiales ni de qué anillo se elija. Se compara
    contra :func:`arc_power_flux` (flujo radial un anillo por dentro) en los tests: son dos rutas
    independientes y tienen que coincidir.
    """
    z = _validate_impedance(impedance, section.fluid)
    edges = window_arc_edges(section, half_angle, center_angle)
    if edges.shape[0] == 0:
        return 0.0
    x = section.nodes
    lengths = np.linalg.norm(x[edges[:, 1]] - x[edges[:, 0]], axis=1)
    m_loc = weights_to_local(lengths)
    pe = np.asarray(p, dtype=complex)[edges]
    dens = np.einsum("ei,ej,eij->e", pe, np.conj(pe), m_loc).real
    return float(0.5 * (1.0 / z).real * np.sum(dens))


def _assemble_and_solve_arc(
    section: CochleaCrossSection,
    *,
    arc_edges: np.ndarray,
    k_arc: complex,
    arc_order: int = 1,
    abc_order: int = 1,
    amplitude: complex = 1.0,
    direction=(1.0, 0.0),
) -> CoupledSolution:
    """Sistema acoplado del corte + término de borde en el arco (impedancia o ABC).

    Es el mismo ensamblado que ``fem.coupled2d.solve_coupled`` —mismos bloques, mismos signos,
    misma numeración de gdl— con **un** término extra: ``−∫ (i·k_arc/ρ_f)·p·q ds`` sobre
    ``arc_edges``. No lleva ``periodic_pairs`` ni ``dirichlet``: el modelo de impedancia no los
    necesita, y así el camino del parche ``p = 0`` sigue siendo el de ``solve_coupled`` intacto. La
    razón de no llamar a ``solve_coupled`` es que sus ``abc_edges`` tienen que ser **frontera** de la
    malla y el arco de la ventana es interfaz fluido/sólido (interior).

    Si ``k_arc`` vale ``k_f`` y ``arc_order=1`` el resultado tiene que ser idéntico al de la ABC-1
    sobre ese arco; con ``arc_edges`` vacío, idéntico a ``solve_cross_section`` sin ventana. Las dos
    cosas se verifican en ``tests/test_cochlea_impedance.py`` (es la prueba de que este ensamblado
    no introduce un factor ni un signo propios).
    """
    x = np.asarray(section.nodes, dtype=float)
    t = np.asarray(section.triangles, dtype=int)
    reg = np.asarray(section.regions)
    n = x.shape[0]
    solid_e = reg == "solid"
    if not np.any(solid_e) or np.all(solid_e):
        raise ValueError("el corte necesita elementos 'solid' y 'fluid'")
    solid_tris = t[solid_e]
    fluid_tris = t[~solid_e]

    solid_node = np.zeros(n, dtype=bool)
    solid_node[solid_tris.ravel()] = True
    fluid_node = np.zeros(n, dtype=bool)
    fluid_node[fluid_tris.ravel()] = True

    reps_u = np.flatnonzero(solid_node)
    reps_p = np.flatnonzero(fluid_node)
    n_u = 2 * reps_u.size
    u_dof = -np.ones((n, 2), dtype=int)
    u_dof[reps_u, 0] = 2 * np.arange(reps_u.size)
    u_dof[reps_u, 1] = 2 * np.arange(reps_u.size) + 1
    p_dof = -np.ones(n, dtype=int)
    p_dof[reps_p] = n_u + np.arange(reps_p.size)
    size = n_u + reps_p.size

    fluid = section.fluid
    frequency = section.frequency
    omega = _TWO_PI * frequency
    k_f = fluid.wavenumber(frequency)

    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    vals: list[np.ndarray] = []

    def _add(rdofs, cdofs, local):
        nr = rdofs.shape[1]
        nc = cdofs.shape[1]
        rows.append(np.repeat(rdofs, nc, axis=1).ravel())
        cols.append(np.tile(cdofs, (1, nr)).ravel())
        vals.append(np.asarray(local).reshape(-1))

    k_s, m_s, _ = p1_elastic_local(x, solid_tris, section.bone)
    _add(
        u_dof[solid_tris].reshape(-1, 6),
        u_dof[solid_tris].reshape(-1, 6),
        k_s - section.bone.rho * omega**2 * m_s,
    )

    m_f, kf, _ = p1_tri_local(x, fluid_tris)
    _add(p_dof[fluid_tris], p_dof[fluid_tris], (kf - k_f**2 * m_f) / fluid.rho)

    iface, iface_n = interface_edges(x, t, reg)
    a_i, b_i = iface[:, 0], iface[:, 1]
    lengths = np.linalg.norm(x[b_i] - x[a_i], axis=1)
    ud = u_dof[iface].reshape(-1, 4)
    pd = p_dof[iface]
    if np.any(ud < 0) or np.any(pd < 0):
        raise RuntimeError("un nodo de interfaz no tiene los tres gdl; revisá la marcación de regions")
    m_loc = weights_to_local(lengths)
    c_local = np.empty((iface.shape[0], 4, 2), dtype=float)
    nx_e = iface_n[:, 0][:, None, None]
    ny_e = iface_n[:, 1][:, None, None]
    c_local[:, 0, :] = (nx_e * m_loc)[:, 0, :]
    c_local[:, 1, :] = (ny_e * m_loc)[:, 0, :]
    c_local[:, 2, :] = (nx_e * m_loc)[:, 1, :]
    c_local[:, 3, :] = (ny_e * m_loc)[:, 1, :]
    _add(ud, pd, c_local)
    _add(pd, ud, omega**2 * np.transpose(c_local, (0, 2, 1)))

    b_full = np.zeros(size, dtype=complex)

    inc = plane_wave(frequency, fluid, direction, amplitude=amplitude)
    abc = np.asarray(section.abc_edges, dtype=int).reshape(-1, 2)
    center_abc = (0.0, 0.0) if abc_order == 2 else None
    blk = abc_boundary_block(
        x, abc, k_f, rho=fluid.rho, abc_order=abc_order, abc_center=center_abc
    )
    _add_blocks(rows, cols, vals, blk, p_dof)
    rhs = abc_incident_rhs(
        x, abc, inc, k_f, rho=fluid.rho, abc_order=abc_order, abc_center=center_abc
    )
    has = p_dof >= 0
    np.add.at(b_full, p_dof[has], rhs[has])

    arc = np.asarray(arc_edges, dtype=int).reshape(-1, 2)
    if arc.shape[0]:
        arc_center = (0.0, 0.0) if arc_order == 2 else None
        a_blk = abc_boundary_block(
            x, arc, complex(k_arc), rho=fluid.rho, abc_order=arc_order, abc_center=arc_center
        )
        _add_blocks(rows, cols, vals, a_blk, p_dof)

    a_full = sp.coo_matrix(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
        shape=(size, size),
        dtype=complex,
    ).tocsr()

    scale = float(np.linalg.norm(b_full))
    full = spla.splu(a_full.tocsc()).solve(b_full)
    resid = a_full @ full - b_full
    residual = float(np.linalg.norm(resid) / (scale if scale > 0.0 else 1.0))

    p = np.full(n, np.nan + 0.0j, dtype=complex)
    u = np.full((n, 2), np.nan + 0.0j, dtype=complex)
    p[has] = full[p_dof[has]]
    has_u = u_dof[:, 0] >= 0
    u[has_u, 0] = full[u_dof[has_u, 0]]
    u[has_u, 1] = full[u_dof[has_u, 1]]

    return CoupledSolution(
        p=p,
        u=u,
        dof_map=CoupledDofs(
            u=u_dof, p=p_dof, n_solid=n_u, n_pressure=reps_p.size, size=size
        ),
        size=size,
        frequency=frequency,
        wavenumber=k_f,
        residual=residual,
        interface_edges=iface,
        interface_normals=iface_n,
    )


def _add_blocks(rows, cols, vals, block: sp.spmatrix, p_dof: np.ndarray) -> None:
    """Agrega un bloque de presión (``N × N``, numeración de nodos) al sistema, vía ``p_dof``."""
    coo = block.tocoo()
    row_dof = p_dof[coo.row]
    col_dof = p_dof[coo.col]
    if np.any(row_dof < 0) or np.any(col_dof < 0):
        raise RuntimeError(
            "el término de borde toca nodos sin gdl de presión; revisá las aristas del arco"
        )
    rows.append(np.asarray(row_dof, dtype=int).ravel())
    cols.append(np.asarray(col_dof, dtype=int).ravel())
    vals.append(np.asarray(coo.data, dtype=complex).ravel())


def wall_injection_power(
    section: CochleaCrossSection,
    solution: CoupledSolution,
) -> tuple[float, float]:
    """``(real, reactiva)`` [W/m] del trabajo ``p·(u·n)`` sobre la cara interna de la cápsula.

    **Diagnóstico, no medida de energía**: se calcula del lado sólido sobre las aristas de interfaz
    en ``r = a``, con ``P = −(ω/2)·Im ∫ p·(u·n)* ds`` (``n`` saliente del sólido, ``v_s = −iω·u``) y
    la masa P1 de arista ``L/6·[[2,1],[1,2]]``. La tracción discreta P1 tiene error ``O(h)``, así que
    la parte real —que debería ser exactamente 0, porque la potencia neta que entra a una cavidad sin
    pérdidas es 0— sale del orden del error: medido ``−0.575`` sin ventana (``ppw=20``), ``−0.538``
    (``ppw=16``) y ``−0.610`` (``ppw=24``) en unidades de ``P_inc``, mientras ``t_interior`` cambia
    sólo 1.3 % entre esas mallas; y el trabajo sobre las dos caras del anillo (``a`` y ``a+t``) no se
    cancela (``−0.575`` contra ``+0.425``) como debería. Por eso las métricas de potencia de
    :class:`CrossSectionSolution` salen **del lado fluido** (:func:`radial_power_flux`,
    :func:`arc_power_flux`) y esta función queda sólo como auditoría de que esa ruta no sirve.
    """
    sol: CoupledSolution = solution
    x = section.nodes
    iface = np.asarray(sol.interface_edges, dtype=int).reshape(-1, 2)
    normals = np.asarray(sol.interface_normals, dtype=float).reshape(-1, 2)
    mid = 0.5 * (x[iface[:, 0]] + x[iface[:, 1]])
    r_mid = np.linalg.norm(mid, axis=1)
    on_inner = np.abs(r_mid - section.a) <= 0.5 * section.h
    iface, normals = iface[on_inner], normals[on_inner]
    if iface.shape[0] == 0:
        return 0.0, 0.0
    length = np.linalg.norm(x[iface[:, 1]] - x[iface[:, 0]], axis=1)
    p_a = np.asarray(sol.p[iface[:, 0]], dtype=complex)
    p_b = np.asarray(sol.p[iface[:, 1]], dtype=complex)
    if not (np.all(np.isfinite(p_a)) and np.all(np.isfinite(p_b))):
        raise ValueError("la cara interna tiene aristas de interfaz sin gdl de presión")
    u_a = np.asarray(sol.u[iface[:, 0]], dtype=complex)
    u_b = np.asarray(sol.u[iface[:, 1]], dtype=complex)
    un_a = np.einsum("ij,ij->i", u_a, normals)
    un_b = np.einsum("ij,ij->i", u_b, normals)
    # ∫ p·(u·n) ds = L/6·[2·p_a·un_a + p_a·un_b + p_b·un_a + 2·p_b·un_b]
    integ = length / 6.0 * (2.0 * p_a * un_a + p_a * un_b + p_b * un_a + 2.0 * p_b * un_b)
    omega = _TWO_PI * section.frequency
    total = complex(np.sum(integ))
    return float(-0.5 * omega * total.imag), float(-0.5 * omega * total.real)


def scattered_power(
    section: CochleaCrossSection,
    p: np.ndarray,
    amplitude: complex = 1.0,
    direction=(1.0, 0.0),
    offset_wavelengths: float = 0.6,
) -> float:
    """Potencia dispersada [W/m] medida en el fluido exterior.

    Es el flujo saliente del campo **dispersado** ``p_s = p − p_i`` (la incidente se resta
    analíticamente nodo por nodo) a través de un círculo en ``r = a + t + offset·λ``. Para un blanco
    **cerrado y sin pérdidas** toda la potencia interceptada vuelve a salir y esto vale
    ``P_inc = 2(a+t)·|p_inc|²/(2ρc)``: la diferencia mide absorción espuria de la ABC y error de
    discretización, no una pérdida física.
    """
    lam = section.fluid.wavelength(section.frequency)
    inc = plane_wave(section.frequency, section.fluid, direction, amplitude=amplitude)
    r_nom = section.a + section.thickness + float(offset_wavelengths) * lam
    rs, nodes = _rings_around(section, r_nom)
    ps = [p[idx] - inc.value(section.nodes[idx]) for idx in nodes]
    return _flux_from_rings(section, ps, rs)


@dataclass(frozen=True)
class CrossSectionSolution:
    """Campo resuelto y las métricas de transmisión del corte.

    Atributos
    ---------
    p : ``(N,)`` presión nodal compleja (``NaN`` donde el nodo no toca fluido).
    u : ``(N, 2)`` desplazamiento nodal complejo (``NaN`` donde no toca sólido).
    solution : :class:`fem.coupled2d.CoupledSolution` completa (gdl, residuo, interfaz).
    section : la geometría.
    incident_amplitude : amplitud compleja de la onda plana incidente [Pa].
    t_interior : ``mean(|p|)`` sobre ``r = a/2`` normalizado por ``|p_inc|``. Es la métrica pedida
        ("presión media en ``r = a/2``"). Incluye la onda estacionaria de la cavidad cerrada.
    t_interior_complex : ``|mean(p)|`` sobre ``r = a/2`` normalizado (responde al armónico ``n = 0``).
    t_interior_rms : ``sqrt(mean(|p|²))`` sobre ``r = a/2`` normalizado.
    t_center : ``|p|`` en el origen normalizado.
    power_scattered_fraction : potencia dispersada ``p_s = p − p_i`` sobre la interceptada. En un
        blanco **cerrado y sin pérdidas** toda la energía vuelve a salir: vale 1. Es el balance de
        energía que esta topología admite (no hay lado transmitido abierto).
    power_cavity_fraction : flujo neto de potencia dentro de la cavidad (``r = 0.75·a``) sobre la
        interceptada. En un sistema sin pérdidas vale 0 (la cavidad no absorbe).
    window_half_angle : semiángulo [grados] efectivamente impuesto, o ``None`` (cavidad cerrada).
    window_model : ``'pressure_release'``, ``'open'`` o ``'impedance'`` (``'abc1_arc'``/
        ``'abc2_arc'`` son las referencias de auditoría de :func:`solve_arc_condition`).
    window_node_count : nodos con ``p = 0`` (modelo ``'pressure_release'``); 0 si no hay ventana o
        si el modelo es otro.
    window_edge_count : aristas del arco con el término de borde (modelo ``'impedance'`` y sus
        referencias ``'abc1_arc'``/``'abc2_arc'``); 0 en los otros modelos.
    window_impedance : ``Z_r`` usado en el modelo ``'impedance'`` (``None`` en los otros).
    power_window_real_fraction : potencia real saliendo de la cavidad por el arco de la ventana,
        sobre la interceptada, medida en el anillo de malla justo por dentro de ``r = a``. Con el
        parche ``p = 0`` es ``≈ 0`` **por construcción**: una pared blanda no transporta potencia
        real (el valor exacto en ``r = a`` es cero; a un anillo por dentro queda la circulación
        reactiva, que está en el campo ``reactive``). Con ``'impedance'`` es **positiva** y del
        orden de ``|p|²/(2·Z_r)``.
    power_window_reactive_fraction : parte reactiva del mismo flujo (intercambio oscilante, sin
        transporte neto de energía).
    power_window_boundary_fraction : potencia real del modelo ``'impedance'`` integrada **sobre el
        arco** (``½∫|p|²·Re(1/Z_r) ds``, :func:`window_boundary_power`), sobre la interceptada. Es
        la medida exacta de la condición, independiente de la diferencia radial; 0 en los otros
        modelos.
    """

    p: np.ndarray
    u: np.ndarray
    solution: CoupledSolution
    section: CochleaCrossSection
    incident_amplitude: complex
    t_interior: float
    t_interior_complex: float
    t_interior_rms: float
    t_center: float
    power_scattered_fraction: float
    power_cavity_fraction: float
    window_half_angle: float | None = None
    window_model: str = "pressure_release"
    window_node_count: int = 0
    window_edge_count: int = 0
    power_window_real_fraction: float = 0.0
    power_window_reactive_fraction: float = 0.0
    window_impedance: complex | None = None
    power_window_boundary_fraction: float = 0.0

    def as_dict(self) -> dict[str, object]:
        return {
            "t_interior": self.t_interior,
            "t_interior_complex": self.t_interior_complex,
            "t_interior_rms": self.t_interior_rms,
            "t_center": self.t_center,
            "power_scattered_fraction": self.power_scattered_fraction,
            "power_cavity_fraction": self.power_cavity_fraction,
            "window_half_angle": self.window_half_angle,
            "window_model": self.window_model,
            "window_node_count": self.window_node_count,
            "window_edge_count": self.window_edge_count,
            "window_impedance": self.window_impedance,
            "power_window_real_fraction": self.power_window_real_fraction,
            "power_window_reactive_fraction": self.power_window_reactive_fraction,
            "power_window_boundary_fraction": self.power_window_boundary_fraction,
            "residual": self.solution.residual,
            "size": self.solution.size,
        }

    def describe(self) -> str:
        if self.window_half_angle is None:
            win = "sin ventana (cavidad cerrada)"
        elif self.window_model == "impedance":
            z = complex(self.window_impedance) if self.window_impedance is not None else complex(0)
            win = (
                f"ventana impedance ±{self.window_half_angle:g}° "
                f"({self.window_edge_count} aristas, "
                f"Z_r = {z.real / 1e6:.4g}{z.imag / 1e6:+.3g}j MRayl)"
            )
        elif self.window_model in ("abc1_arc", "abc2_arc"):
            win = (
                f"referencia {self.window_model} ±{self.window_half_angle:g}° "
                f"({self.window_edge_count} aristas)"
            )
        else:
            win = (
                f"ventana {self.window_model} ±{self.window_half_angle:g}° "
                f"({self.window_node_count} nodos p=0)"
            )
        boundary = (
            f"   P_ventana(borde)/P_inc = {self.power_window_boundary_fraction:+.3e}"
            if self.window_model == "impedance"
            else ""
        )
        return (
            f"{self.section.describe()}\n"
            f"{win}\n"
            f"t_interior(|p| en r=a/2) = {self.t_interior:.5f}   "
            f"t_interior_complex = {self.t_interior_complex:.5f}   "
            f"t_center = {self.t_center:.5f}\n"
            f"P_disp/P_inc = {self.power_scattered_fraction:.5f}   "
            f"P_cavidad/P_inc = {self.power_cavity_fraction:.5f}\n"
            f"P_ventana/P_inc = {self.power_window_real_fraction:+.3e} "
            f"(reactiva {self.power_window_reactive_fraction:+.3e}){boundary}   "
            f"residuo = {self.solution.residual:.1e}"
        )


def _measure_solution(
    section: CochleaCrossSection,
    sol: CoupledSolution,
    *,
    half: float | None,
    center_angle: float,
    model: str,
    amplitude: complex,
    direction,
    impedance: complex | None = None,
    node_count: int = 0,
    edge_count: int = 0,
) -> CrossSectionSolution:
    """Métricas de transmisión/potencia sobre un campo ya resuelto (común a todos los modelos)."""
    amp = abs(complex(amplitude))
    if amp == 0.0:
        raise ValueError("la amplitud incidente no puede ser nula")

    # --- circunferencia r = a/2 (métrica pedida) ----------------------------
    ring = section.ring_nodes(section.probe_radius)
    vals = sol.p[ring]
    if not np.all(np.isfinite(vals)):
        raise RuntimeError("la circunferencia de la sonda tiene nodos sin gdl de presión")
    t_interior = float(np.mean(np.abs(vals))) / amp
    t_interior_complex = float(abs(np.mean(vals))) / amp
    t_interior_rms = float(np.sqrt(np.mean(np.abs(vals) ** 2))) / amp

    origin = int(np.argmin(np.linalg.norm(section.nodes, axis=1)))
    t_center = float(abs(sol.p[origin])) / amp

    # --- potencias ----------------------------------------------------------
    p_inc_total = intercepted_power(section, amplitude)
    power_scattered = scattered_power(section, sol.p, amplitude, direction)
    power_cavity = radial_power_flux(section, sol.p, 0.75 * section.a) / p_inc_total

    win_real = win_react = win_boundary = 0.0
    if half is not None:
        win_real, win_react = arc_power_flux(section, sol.p, half, center_angle=center_angle)
        win_real /= p_inc_total
        win_react /= p_inc_total
        if model == "impedance":
            win_boundary = (
                window_boundary_power(
                    section, sol.p, impedance=impedance,
                    half_angle=half, center_angle=center_angle,
                )
                / p_inc_total
            )

    return CrossSectionSolution(
        p=sol.p,
        u=sol.u,
        solution=sol,
        section=section,
        incident_amplitude=complex(amplitude),
        t_interior=t_interior,
        t_interior_complex=t_interior_complex,
        t_interior_rms=t_interior_rms,
        t_center=t_center,
        power_scattered_fraction=float(power_scattered / p_inc_total),
        power_cavity_fraction=float(power_cavity),
        window_half_angle=half,
        window_model=model,
        window_node_count=int(node_count),
        window_edge_count=int(edge_count),
        power_window_real_fraction=float(win_real),
        power_window_reactive_fraction=float(win_react),
        window_impedance=None if impedance is None else complex(impedance),
        power_window_boundary_fraction=float(win_boundary),
    )


def solve_arc_condition(
    section: CochleaCrossSection,
    *,
    half_angle: float | None = None,
    center_angle: float | None = None,
    condition: str = "impedance",
    impedance: complex | None = None,
    abc_order: int = 1,
    amplitude: complex = 1.0,
    direction=(1.0, 0.0),
) -> CrossSectionSolution:
    """Resuelve el corte con una condición de borde sobre el **arco** ``r = a``.

    Es el motor del modelo ``window_model='impedance'`` y, con las condiciones de referencia, la
    forma de auditar el ensamblado en ``tests/test_cochlea_impedance.py``.

    Parameters
    ----------
    half_angle, center_angle : arco en grados; ``None`` usa los de la sección.
    condition :

        - ``'impedance'``: ``p = Z_r·(v·n)``, o sea ``−∫(i·k_w/ρ_f)·p·q ds`` con
          ``k_w = ω·ρ_f/Z_r`` (el término de la ABC-1 con ``k → k_w``). Es la física del modelo.
        - ``'abc1'``: **referencia**, la ABC-1 de siempre aplicada sobre el mismo arco
          (``k_w = k_f``). Con ``Z_r = ρ_f·c_f`` tiene que coincidir con ``'impedance'``.
        - ``'abc2'``: **referencia**, BGT-2 sobre el mismo arco (``abc_center = 0``, ``r_e = a``).
          Es otra condición (agrega los términos de curvatura y angular), no equivalente a la
          impedancia: sirve para cuantificar cuánto aportan esos términos en un arco de ``ka ≈ 8``.
    impedance : ``Z_r`` [Pa·s/m]; ``None`` ⇒ ``ρ_f·c_f``.
    abc_order : orden de la ABC del contorno **exterior** (``r = R_out``), igual que en
        :func:`solve_cross_section`.

    Returns
    -------
    :class:`CrossSectionSolution` (``window_model`` = ``condition``, con ``'abc1_arc'``/
    ``'abc2_arc'`` para las referencias).
    """
    if not isinstance(section, CochleaCrossSection):
        raise TypeError(f"section debe ser CochleaCrossSection; recibido {type(section).__name__}")
    half, center, _ = _validate_window(
        section.window_half_angle if half_angle is None else half_angle,
        section.window_center_angle if center_angle is None else center_angle,
        "impedance",
    )
    if half is None:
        raise ValueError("solve_arc_condition necesita un arco: window_half_angle no puede ser None")
    if abc_order not in (1, 2):
        raise ValueError(f"abc_order debe ser 1 o 2; recibido {abc_order!r}")
    condition = str(condition)
    if condition not in ("impedance", "abc1", "abc2"):
        raise ValueError(
            f"condition debe ser 'impedance', 'abc1' o 'abc2'; recibido {condition!r}"
        )

    k_f = section.fluid.wavenumber(section.frequency)
    z = _validate_impedance(impedance, section.fluid)
    edges = window_arc_edges(section, half, center)
    if edges.shape[0] == 0:
        raise ValueError(
            f"el arco de ±{half:g}° no contiene aristas de interfaz en r = a "
            f"(n_theta={section.n_theta}: cada sector son {360.0 / section.n_theta:.3f}°)"
        )
    if condition == "impedance":
        # p = Z_r·v_n  con  v_n = ∂p/∂n/(iωρ)  ⇒  ∂p/∂n = (i·ω·ρ/Z_r)·p  ⇒  k_w = ω·ρ/Z_r
        k_arc = _TWO_PI * section.frequency * section.fluid.rho / z
        arc_order = 1
        model = "impedance"
        z_used = z
    else:
        k_arc = k_f
        arc_order = 1 if condition == "abc1" else 2
        model = condition + "_arc"
        z_used = None

    sol = _assemble_and_solve_arc(
        section,
        arc_edges=edges,
        k_arc=k_arc,
        arc_order=arc_order,
        abc_order=abc_order,
        amplitude=amplitude,
        direction=direction,
    )
    return _measure_solution(
        section,
        sol,
        half=half,
        center_angle=center,
        model=model,
        amplitude=amplitude,
        direction=direction,
        impedance=z_used,
        node_count=0,
        edge_count=edges.shape[0],
    )


def solve_cross_section(
    section: CochleaCrossSection,
    *,
    amplitude: complex = 1.0,
    direction=(1.0, 0.0),
    window_half_angle: float | None = None,
    window_center_angle: float | None = None,
    window_model: str | None = None,
    window_impedance: complex | None = None,
    abc_order: int = 1,
) -> CrossSectionSolution:
    """Resuelve el corte con el solver acoplado de la Fase 2 y mide la transmisión.

    Parameters
    ----------
    section : :class:`CochleaCrossSection` de :func:`build_cross_section`.
    amplitude, direction : onda plana incidente en el fluido exterior (amplitud y dirección).
    window_half_angle : semiángulo [**grados**] del arco de ventana redonda sobre ``r = a``.
        ``None`` (default) usa lo que traiga la sección — que en el caso de la Fase 5a también es
        ``None``, o sea **cavidad cerrada, comportamiento idéntico a la Fase 5a**. Un valor acá
        sobreescribe (y permite abrir la ventana sin remallar) el de la sección.
    window_center_angle : centro del arco [grados]; ``None`` usa el de la sección (``0`` = ``+x``,
        la cara de sombra).
    window_model : ``None`` usa el de la sección. ``'pressure_release'`` impone ``p = 0`` en los
        nodos del arco (idealización declarada de la ventana redonda); ``'open'`` sólo tiene sentido
        si la malla se construyó con ``window_model='open'`` (no queda pared que presurizar);
        ``'impedance'`` impone ``p = Z_r·(v·n)`` sobre las aristas del arco y **transporta potencia
        real** hacia afuera (delega en :func:`solve_arc_condition`).
    window_impedance : ``Z_r`` [Pa·s/m] del modelo ``'impedance'``; ``None`` ⇒ ``ρ_f·c_f``.
    abc_order : orden de la ABC del contorno exterior. ``1`` (default) es el de la Fase 5a: los
        números de 5a se reproducen exactos. ``2`` (BGT-2) mejora la cara de sombra (D21) y usa
        ``abc_center=(0, 0)``.

    Returns
    -------
    :class:`CrossSectionSolution`.

    Notes
    -----
    La ABC se aplica sobre ``r = R_out`` con el término incidente sólo en las caras iluminadas
    (``n·d < 0``), que es lo que hace `fem.coupled2d.solve_coupled`. El fluido es el mismo dentro y
    fuera (``section.fluid``), así que las potencias son directamente comparables.

    La ventana ``pressure_release`` entra como ``dirichlet={'p': nodos del arco}``: gdl esenciales
    nulos, sin tocar el solver. **No es una salida de energía** (ver el docstring del módulo): es una
    pared blanda, y el balance de energía se cierra con todos los términos en cero. El modelo
    ``'impedance'`` sí es una salida real y su término de borde se ensambla en
    :func:`solve_arc_condition`.
    """
    if not isinstance(section, CochleaCrossSection):
        raise TypeError(
            f"section debe ser CochleaCrossSection; recibido {type(section).__name__}"
        )
    half, center, model = _validate_window(
        section.window_half_angle if window_half_angle is None else window_half_angle,
        section.window_center_angle if window_center_angle is None else window_center_angle,
        section.window_model if window_model is None else window_model,
    )
    if model == "open" and section.window_model != "open":
        raise ValueError(
            "window_model='open' requiere que la malla se haya construido con "
            "build_cross_section(..., window_model='open'): los elementos sólidos del arco tienen "
            "que estar fuera, si no la ventana queda tapada por la pared"
        )
    if abc_order not in (1, 2):
        raise ValueError(f"abc_order debe ser 1 o 2; recibido {abc_order!r}")

    if model == "impedance":
        z = _validate_impedance(
            section.window_impedance if window_impedance is None else window_impedance,
            section.fluid,
        )
        if half is None:
            raise ValueError(
                "window_model='impedance' necesita un arco (window_half_angle > 0): sin ventana "
                "no hay superficie donde imponer p = Z_r·(v·n)"
            )
        return solve_arc_condition(
            section,
            half_angle=half,
            center_angle=center,
            condition="impedance",
            impedance=z,
            abc_order=abc_order,
            amplitude=amplitude,
            direction=direction,
        )

    f = section.frequency
    inc: PlaneWave = plane_wave(f, section.fluid, direction, amplitude=amplitude)
    wnodes = _ring_arc(section, section.a, half, center) if half is not None else np.zeros(0, int)
    if half is not None and wnodes.size < 2:
        raise ValueError(
            f"el arco de ±{half:g}° no contiene nodos de la malla en r = a "
            f"(n_theta={section.n_theta}: cada sector son {360.0 / section.n_theta:.3f}°)"
        )
    extra = {}
    if abc_order == 2:
        extra["abc_order"] = 2
        extra["abc_center"] = (0.0, 0.0)
    dirichlet = {"p": wnodes} if (half is not None and model == "pressure_release") else None
    sol = solve_coupled(
        section.nodes,
        section.triangles,
        section.regions,
        section.fluid,
        section.bone,
        abc_edges=section.abc_edges,
        p_inc=inc,
        frequency=f,
        dirichlet=dirichlet,
        **extra,
    )
    return _measure_solution(
        section,
        sol,
        half=half,
        center_angle=center,
        model=model,
        amplitude=amplitude,
        direction=direction,
        node_count=int(wnodes.size) if dirichlet is not None else 0,
    )
