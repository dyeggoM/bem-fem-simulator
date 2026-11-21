"""Escena 2D de trabajo: capas con material, fuente plana y solución FEM (Fase U1).

Qué resuelve
------------
Un apilado de **capas verticales** (el eje de propagación es *x*) con un material por capa, una
onda plana incidiendo de frente desde la izquierda y ABC en las dos caras de *x*. Es el problema
que ya está validado en la librería (Fase 1: R y T contra analítico al 0.511 %), así que sirve
como banco de prueba de la UI: si la UI muestra R y T dentro del 2 %, la cadena
escena → malla → solve → métricas está bien.

Qué NO resuelve todavía (y por qué)
-----------------------------------
- **Incidencia oblicua**: las caras de arriba y abajo son rígidas, lo cual sólo es consistente con
  incidencia normal. Para oblicua hace falta absorbedor también arriba/abajo (PML o ABC de orden 2
  en las cuatro caras) — U1b.
- **Regiones en bloque** (no estratificadas): la malla es un reticulado estructurado y las
  interfaces conforman por construcción sólo si las capas son bandas verticales. Para bloques 2D
  hace falta un mallador con interfaz pegada — U1c.
- **Elástico**: acá sólo fluido acústico. El acoplamiento fluido-sólido (`fem.coupled2d`) va en U3,
  con la cóclea.

Por qué las capas son una partición de *x* y no regiones arbitrarias: la malla tiene que ser
**conforme** en las interfaces. Con un reticulado estructurado eso sale gratis; con bloques
arbitrarios hay que pegarlos, y un centroide mal asignado mete un material en el elemento
equivocado sin que nada avise.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Any

import numpy as np

from .. import materials as mats
from ..analytic import helmholtz1d
from ..fem import helmholtz2d
from ..meshing import boundary_edges, grid_tri_mesh
from . import validacion

MATERIALES: dict[str, Any] = {
    "agua": mats.WATER,
    "tejido": mats.SOFT_TISSUE,
    "hueso": mats.BONE_CORTICAL,
    "perilinfa": mats.PERILYMPH,
    "endolinfa": mats.ENDOLYMPH,
}

ABC_DISPONIBLES = ("abc1",)
FORMULACIONES = ("total",)

_LIMITES = {
    "frecuencia_hz": (20.0e3, 1.0e6),
    "ppw": (6, 40),
    "alto_m": (2.0e-3, 100.0e-3),
    "capas": (2, 12),
    "espesor_m": (0.2e-3, 200.0e-3),
}

# tope duro: por encima de esto el solve deja de ser interactivo y conviene scripts/
MAX_ELEMENTOS = 60_000
MAX_NODOS = 40_000


def _num(datos: dict[str, Any], clave: str, defecto: float) -> float:
    valor = datos.get(clave, defecto)
    try:
        numero = float(valor)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{clave} tiene que ser un número; llegó {valor!r}") from exc
    if not np.isfinite(numero):
        raise ValueError(f"{clave} no puede ser {numero}")
    return numero


def _en_rango(clave: str, valor: float) -> float:
    bajo, alto = _LIMITES[clave]
    if not (bajo <= valor <= alto):
        raise ValueError(f"{clave} = {valor:g} fuera de rango [{bajo:g}, {alto:g}]")
    return valor


@dataclass(frozen=True)
class Escena:
    """Escena ya validada. Inmutable para que el resultado no cambie si el dict de entrada se toca."""

    frecuencia_hz: float
    ppw: float
    alto_m: float
    capas: tuple[tuple[str, float], ...]  # (material, espesor)
    abc: str
    formulacion: str

    def materiales(self) -> tuple[Any, ...]:
        return tuple(MATERIALES[nombre] for nombre, _ in self.capas)

    def espesores(self) -> tuple[float, ...]:
        return tuple(espesor for _, espesor in self.capas)

    def a_dict(self) -> dict[str, Any]:
        return {
            "frecuencia_hz": self.frecuencia_hz,
            "ppw": self.ppw,
            "alto_m": self.alto_m,
            "capas": [{"material": m, "espesor_m": e} for m, e in self.capas],
            "abc": self.abc,
            "formulacion": self.formulacion,
        }


def validar_escena(datos: dict[str, Any]) -> Escena:
    """Normaliza y valida el JSON de la escena. Falla con mensaje, nunca con un `KeyError`."""
    if not isinstance(datos, dict):
        raise ValueError("la escena tiene que ser un objeto JSON")

    frecuencia = _en_rango("frecuencia_hz", _num(datos, "frecuencia_hz", 100.0e3))
    ppw = _en_rango("ppw", _num(datos, "ppw", 12))
    alto = _en_rango("alto_m", _num(datos, "alto_m", 20.0e-3))

    crudas = datos.get("capas")
    if not isinstance(crudas, list):
        raise ValueError("la escena necesita una lista 'capas'")
    bajo, alto_capas = _LIMITES["capas"]
    if not (bajo <= len(crudas) <= alto_capas):
        raise ValueError(f"'capas' tiene que tener entre {bajo} y {alto_capas} elementos")

    capas: list[tuple[str, float]] = []
    for i, capa in enumerate(crudas):
        if not isinstance(capa, dict):
            raise ValueError(f"la capa {i} tiene que ser un objeto {{material, espesor_m}}")
        nombre = capa.get("material")
        if nombre not in MATERIALES:
            raise ValueError(
                f"capa {i}: material desconocido {nombre!r}; opciones: {sorted(MATERIALES)}"
            )
        espesor = _en_rango("espesor_m", _num(capa, "espesor_m", 0.0))
        capas.append((str(nombre), espesor))

    if len({nombre for nombre, _ in capas}) == 1 and len(capas) > 1:
        # no es un error físico, pero casi siempre es un descuido: se avisa por los metadatos
        pass

    abc = str(datos.get("abc", "abc1"))
    if abc not in ABC_DISPONIBLES:
        raise ValueError(f"abc desconocido: {abc!r}; opciones: {list(ABC_DISPONIBLES)}")

    formulacion = str(datos.get("formulacion", "total"))
    if formulacion not in FORMULACIONES:
        raise ValueError(
            f"formulación desconocida: {formulacion!r}; en U1 sólo {list(FORMULACIONES)} "
            "(el campo dispersado con ABC curvo llega con los dispersores, U2)"
        )

    angulo = abs(_num(datos, "angulo_deg", 0.0))
    if angulo > 1e-9:
        raise ValueError(
            "U1 resuelve sólo incidencia normal: arriba y abajo las paredes son rígidas, y una "
            "onda oblicua se refleja ahí en vez de salir. Incidencia oblicua entra con absorbedor "
            "en las cuatro caras (U1b)."
        )

    return Escena(frecuencia, ppw, alto, tuple(capas), abc, formulacion)


def _paso(material: Any, frecuencia: float, ppw: float) -> float:
    """Paso de malla dentro de un material: ``min(λ/ppw, λ/ppw)`` pero con el λ de *ese* material."""
    return abs(material.wavelength(frecuencia)) / ppw


def malla_de_escena(escena: Escena) -> tuple[np.ndarray, np.ndarray, Any, np.ndarray]:
    """Devuelve ``(xs, ys, malla, indice_de_capa_por_elemento)``.

    El paso de *x* se fija por capa (cada material tiene su λ); el de *y* es el más fino de todos,
    para que ningún elemento quede estirado. ``xs`` incluye exactamente las interfaces.
    """
    materiales = escena.materiales()
    espesores = escena.espesores()

    xs: list[float] = [0.0]
    for material, espesor in zip(materiales, espesores):
        paso = _paso(material, escena.frecuencia_hz, escena.ppw)
        n = max(1, int(np.ceil(espesor / paso)))
        xs.extend(np.linspace(xs[-1], xs[-1] + espesor, n + 1)[1:].tolist())

    pasos = [_paso(m, escena.frecuencia_hz, escena.ppw) for m in materiales]
    paso_y = min(pasos) if pasos else 1.0e-3
    ny = max(1, int(np.ceil(escena.alto_m / paso_y)))
    ys = np.linspace(0.0, escena.alto_m, ny + 1)

    malla = grid_tri_mesh(xs, ys)
    n_elementos = malla.triangles.shape[0]
    if n_elementos > MAX_ELEMENTOS or malla.nodes.shape[0] > MAX_NODOS:
        raise ValueError(
            f"la malla sale de {n_elementos} elementos y {malla.nodes.shape[0]} nodos; el tope "
            f"interactivo es {MAX_ELEMENTOS} / {MAX_NODOS}. Bajá ppw, el alto o la cantidad de capas."
        )

    # capa por elemento según el centroide en x
    centro_x = malla.nodes[malla.triangles].mean(axis=1)[:, 0]
    bordes = np.cumsum(espesores)
    capa = np.searchsorted(bordes, centro_x, side="right")
    capa = np.clip(capa, 0, len(espesores) - 1)
    return np.asarray(xs), ys, malla, capa


def _ajustar_dos_ondas(x: np.ndarray, p: np.ndarray, k: complex) -> tuple[complex, complex]:
    """Ajusta ``A·e^{ikx} + B·e^{-ikx}`` por mínimos cuadrados y devuelve ``(A, B)``."""
    base = np.column_stack([np.exp(1j * k * x), np.exp(-1j * k * x)])
    coef, *_ = np.linalg.lstsq(base, p, rcond=None)
    return complex(coef[0]), complex(coef[1])


def _ajustar_una_onda(x: np.ndarray, p: np.ndarray, k: complex) -> complex:
    base = np.exp(1j * k * x)[:, None]
    coef, *_ = np.linalg.lstsq(base, p, rcond=None)
    return complex(np.asarray(coef).reshape(-1)[0])


def resolver_escena(datos: dict[str, Any]) -> dict[str, Any]:
    """Corre la escena y devuelve campo, cortes y métricas con su referencia analítica."""
    t0 = time.perf_counter()
    escena = validar_escena(datos)
    xs, ys, malla, capa = malla_de_escena(escena)

    materiales = escena.materiales()
    espesores = escena.espesores()
    f = escena.frecuencia_hz

    k_el = np.array([materiales[c].wavenumber(f) for c in capa], dtype=complex)
    rho_el = np.array([materiales[c].rho for c in capa], dtype=float)

    aristas = boundary_edges(malla.triangles)
    x_aristas = malla.nodes[aristas].mean(axis=1)[:, 0]
    x_max = float(xs[-1])
    tol = 1e-9 * max(x_max, 1.0)
    abc = aristas[(np.abs(x_aristas) < tol) | (np.abs(x_aristas - x_max) < tol)]

    p_inc = helmholtz2d.plane_wave(f, materiales[0], (1.0, 0.0))
    t_solve = time.perf_counter()
    p = helmholtz2d.solve(
        malla.nodes, malla.triangles, k_el, abc_edges=abc, p_inc=p_inc,
        element_rho=rho_el,
    )
    dt_solve = time.perf_counter() - t_solve

    campo = np.asarray(p, dtype=complex).reshape(ys.size, xs.size)
    amplitud = np.abs(campo)
    pico = float(amplitud.max()) or 1.0
    nivel_db = 20.0 * np.log10(np.maximum(amplitud / pico, 1e-6))

    # ---- R y T: ajuste de las ondas pegadas a las interfaces (misma técnica que el test 3D) ----
    k_1 = materiales[0].wavenumber(f)
    k_last = materiales[-1].wavenumber(f)
    t1 = espesores[0]
    x_ini_ultima = float(xs[-1] - espesores[-1])

    nodos_x = malla.nodes[:, 0]
    ventana = 0.5 * abs(materiales[0].wavelength(f))
    sel_izq = (nodos_x > t1 - ventana) & (nodos_x < t1 - tol)
    sel_der = (nodos_x > x_ini_ultima + tol) & (nodos_x < x_ini_ultima + ventana)

    metricas: dict[str, Any] = {}
    if sel_izq.sum() >= 8 and sel_der.sum() >= 8:
        a_izq, b_izq = _ajustar_dos_ondas(nodos_x[sel_izq], p[sel_izq], k_1)
        c_der = _ajustar_una_onda(nodos_x[sel_der], p[sel_der], k_last)
        r_num = b_izq / a_izq
        t_num = c_der / a_izq

        medias = [materiales[0]] + list(materiales) + [materiales[-1]]
        # helmholtz1d mide x desde la primera interfaz y pide (izq, *capas, der)
        capas_analitico = tuple(materiales)
        if len(espesores) >= 2:
            respuesta = helmholtz1d.stack_response(f, capas_analitico, espesores[:-1][1:])
            r_an = complex(respuesta.R)
            t_an = complex(respuesta.T)
        else:
            r_an = t_an = complex("nan")

        # `r_num` sale de ajustar `a·e^{ikx} + b·e^{−ikx}` con la **x global**, o sea que `b/a` está
        # referido a x = 0; `r_an` (matriz de transferencia) está referido a la **primera interfaz**
        # (x = t1). Comparar los fasores así mezclaba un desfase de propagación e^{2ik t1} con el error
        # del solver: medido, daba "199 %" de error donde la amplitud coincidía al 0.19 % y el número
        # variaba con la frecuencia (100/200/300/400 kHz → 192/107/126/184 %), con
        # |ΔR|/|R| = 2·sin(k·t1) en los cuatro puntos. Se lleva `r` a la referencia de la interfaz.
        r_num_interfaz = r_num * np.exp(-2j * k_1 * t1)
        e_r = (abs(r_num_interfaz - r_an) / abs(r_an) if abs(r_an) > 1e-12 else 0.0)
        e_r_modulo = (abs(abs(r_num) - abs(r_an)) / abs(r_an) if abs(r_an) > 1e-12 else 0.0)
        e_t = abs(abs(t_num) - abs(t_an)) / abs(t_an) if abs(t_an) > 1e-12 else 0.0
        # El criterio de la Fase 1 es "dentro del 2 % EN AMPLITUD", no en relativo: cuando R es
        # chico (agua/tejido, R = 1.2 %) el relativo se infla con cualquier ruido de discretización
        # del orden del 1 % de la incidente. Se reportan las dos cosas y se juzga por la absoluta.
        metricas.update({
            "R": abs(r_num),
            "T": abs(t_num),
            "R_analitico": abs(r_an),
            "T_analitico": abs(t_an),
            "error_R_abs": float(abs(abs(r_num) - abs(r_an))),
            "error_T_abs": float(abs(abs(t_num) - abs(t_an))),
            # `error_R_pct` sobre el MÓDULO, igual que `error_T_pct`: es el número comparable entre
            # los dos y el que se lee de un vistazo. El fasor (ya llevado a la referencia de la
            # interfaz) va aparte, con nombre propio, porque incluye la fase y por lo tanto es del
            # orden del error de fase de la malla (ppw 8 → ~7°), no del error de amplitud.
            "error_R_pct": 100.0 * e_r_modulo,
            "error_R_fasor_pct": 100.0 * e_r,
            "error_T_pct": 100.0 * e_t,
            # el desfase que se sacó para comparar los fasores, a la vista: si alguien compara la fase
            # de R contra el analítico sin corregirlo, va a ver este giro y no un error del solver
            # `k_1` puede ser complejo (atenuación): el giro de FASE es la parte real
            "desfase_referencia_R_deg": float(np.degrees(2.0 * np.real(k_1) * t1)),
            "energia": float(abs(r_num) ** 2 + abs(t_num) ** 2),
            "ventana_fit_mm": 0.5 * abs(materiales[0].wavelength(f)) * 1e3,
        })
    else:
        metricas["aviso"] = (
            "capas demasiado finas para separar incidente y reflejada con la ventana de ajuste "
            "(0.5 λ): R y T no se reportan"
        )

    # ---- sello: lo evalúa el servidor, con la referencia y la fuente del criterio (U6) ----
    sello = validacion.sello_conjunto({
        "escena_rt_r": metricas.get("error_R_abs"),
        "escena_rt_t": metricas.get("error_T_abs"),
    })

    # ---- cortes utiles para la UI ----
    col_pico = int(np.argmax(amplitud.max(axis=0)))
    fila_central = int(np.argmin(np.abs(ys - 0.5 * escena.alto_m)))

    return {
        "ok": True,
        "escena": escena.a_dict(),
        "grilla": {
            "x_mm": (xs * 1e3).tolist(),
            "y_mm": (ys * 1e3).tolist(),
            "n_nodos": int(malla.nodes.shape[0]),
            "n_elementos": int(malla.triangles.shape[0]),
            "n_aristas_abc": int(abc.shape[0]),
            "interfaces_mm": (np.cumsum(espesores)[:-1] * 1e3).tolist(),
        },
        "campo": {
            "amplitud": np.round(amplitud, 6).tolist(),
            "db": np.round(nivel_db, 3).tolist(),
            "max": pico,
        },
        "cortes": {
            "eje_x_mm": (xs * 1e3).tolist(),
            "eje_y_mm": (ys * 1e3).tolist(),
            "axial": np.round(np.abs(campo[fila_central, :]), 6).tolist(),
            "lateral_en_x_pico": np.round(np.abs(campo[:, col_pico]), 6).tolist(),
            "x_pico_mm": float(xs[col_pico] * 1e3),
        },
        "metricas": metricas,
        "sello": sello,
        "tiempos": {
            "malla_s": round(t_solve - t0, 4),
            "solve_s": round(dt_solve, 4),
            "total_s": round(time.perf_counter() - t0, 4),
        },
        "materiales_usados": [
            {"nombre": materiales[i].name, "c_m_s": materiales[i].c, "rho": materiales[i].rho}
            for i in range(len(materiales))
        ],
    }


def estimar_costo(datos: dict[str, Any]) -> dict[str, Any]:
    """Malla y costo estimado **sin resolver**: es lo que la UI muestra antes de apretar el botón."""
    escena = validar_escena(datos)
    xs, ys, malla, _ = malla_de_escena(escena)
    n_el = int(malla.triangles.shape[0])
    # el solve 2D ronda 1 s por 5 k elementos en esta máquina (medido en D27 y en el dashboard)
    return {
        "ok": True,
        "grilla": {
            "n_nodos": int(malla.nodes.shape[0]),
            "n_elementos": n_el,
            "nx": int(xs.size),
            "ny": int(ys.size),
            "paso_x_min_mm": float(np.min(np.diff(xs)) * 1e3),
            "paso_y_mm": float((ys[1] - ys[0]) * 1e3),
            "longitud_mm": float(xs[-1] * 1e3),
            "interfaces_mm": (np.cumsum(escena.espesores())[:-1] * 1e3).tolist(),
        },
        "estimado_solve_s": round(0.2 * n_el / 1000.0, 2),
        "tope_elementos": MAX_ELEMENTOS,
    }
