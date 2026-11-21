"""Corte transversal de cóclea desde la UI (Fase U3): campo interior, no sólo los chips.

Qué resuelve
------------
La escena de la Fase 5 completa: fluido interior (perilinfa), cápsula elástica con espesor y
``c_s`` configurables, ventana redonda y membrana basilar. El solver es `cochlea.solve_cross_section`
—el acoplado fluido-sólido de la Fase 2—, exactamente el que midió los números de D29/D30, sin
reimplementar nada.

Lo que la UI agrega sobre el panel de D31: el **mapa del campo interior** (|p| interpolado sobre el
disco) y el balance de energía completo, además de ``t_interior`` y el armónico ``n = 0``.

Qué NO hace
-----------
- **Q (ancho de banda de la resonancia)**: no se mide acá. Necesita un barrido de frecuencia alrededor
  del pico (≈ 9 solves a 2 s cada uno) y eso es el trabajo de U5, no de un botón de esta pantalla.
  La UI lo dice en vez de mostrar un número que no midió.
- **Variación a lo largo de la espira**: esto es un corte transversal, no la cóclea desenrollada.
  Un solo ``a``, un solo espesor, una sola membrana.
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np
from scipy.interpolate import LinearNDInterpolator

from .. import cochlea, materials as mats
from . import validacion

FRECUENCIA_HZ = 1.0e5
RADIO_M = 20.0e-3

_LIMITES = {
    "thickness_mm": (0.5, 8.0),
    "c_s": (1000.0, 2500.0),
    "alpha_np_m": (0.0, 40.0),
    "window_half_angle": (0.0, 40.0),
    "ppw": (12.0, 28.0),
}

#: c_p y ρ del hueso: ballpark VERIFY de D20, no un valor de la literatura. Se fijan acá y se
#: declaran en la respuesta para que la UI pueda decirlo.
HUESO_CP = 3500.0
HUESO_RHO = 1950.0

PUNTOS_MAPA = 121
MAX_NODOS = 120_000


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


def validar_cochlea(datos: dict[str, Any]) -> dict[str, Any]:
    """Normaliza y valida. Falla con mensaje, nunca con un `KeyError`."""
    if not isinstance(datos, dict):
        raise ValueError("el pedido tiene que ser un objeto JSON")

    thickness_mm = _en_rango("thickness_mm", _num(datos, "thickness_mm", 2.0))
    c_s = _en_rango("c_s", _num(datos, "c_s", 1900.0))
    if c_s >= HUESO_CP:
        raise ValueError(f"c_s ({c_s:g}) tiene que ser menor que c_p ({HUESO_CP:g} m/s)")
    alpha = _en_rango("alpha_np_m", _num(datos, "alpha_np_m", 0.0))
    ppw = _en_rango("ppw", _num(datos, "ppw", 20.0))

    membrana = bool(datos.get("membrana", True))
    ventana = bool(datos.get("ventana", True))
    semiangulo = _en_rango("window_half_angle", _num(datos, "window_half_angle", 10.0))
    if ventana and semiangulo <= 0.0:
        raise ValueError("con ventana abierta, window_half_angle tiene que ser > 0")

    return {
        "thickness_mm": thickness_mm,
        "c_s": c_s,
        "alpha_np_m": alpha,
        "ppw": ppw,
        "membrana": membrana,
        "ventana": ventana,
        "window_half_angle": semiangulo,
    }


def _seccion(p: dict[str, Any]):
    hueso = mats.ElasticSolid(
        name="hueso configurado en la UI", rho=HUESO_RHO, c_p=HUESO_CP, c_s=p["c_s"],
        provenance=f"editado en el dashboard; c_p={HUESO_CP:g} y rho={HUESO_RHO:g} son el "
                   f"ballpark VERIFY de D20",
    )
    kwargs: dict[str, Any] = dict(
        ppw=p["ppw"], pad_wavelengths=2.5, min_elements_across_wall=2, bone=hueso,
    )
    if p["ventana"]:
        kwargs.update(
            window_half_angle=p["window_half_angle"], window_center_angle=0.0,
            window_model="impedance",
            window_impedance=complex(mats.WATER.rho * mats.WATER.c),
        )
    if p["membrana"]:
        kwargs["membrane_angle"] = 90.0
    return cochlea.build_cross_section(
        FRECUENCIA_HZ, RADIO_M, p["thickness_mm"] * 1e-3, **kwargs
    )


def _json_limpio(valor: Any) -> Any:
    """Deja el diccionario serializable: complejos a ``{re, im}`` y numpy a tipos nativos.

    `sol.as_dict()` trae ``window_impedance`` complejo: sin esto, `json.dumps` revienta con
    "Object of type complex is not JSON serializable" y el navegador no recibe nada.
    """
    if isinstance(valor, dict):
        return {str(k): _json_limpio(v) for k, v in valor.items()}
    if isinstance(valor, (list, tuple)):
        return [_json_limpio(v) for v in valor]
    if isinstance(valor, complex):
        return {"re": valor.real, "im": valor.imag}
    if isinstance(valor, np.integer):
        return int(valor)
    if isinstance(valor, np.floating):
        return None if not np.isfinite(float(valor)) else float(valor)
    if isinstance(valor, np.bool_):
        return bool(valor)
    if isinstance(valor, np.ndarray):
        return _json_limpio(valor.tolist())
    if isinstance(valor, float) and not np.isfinite(valor):
        return None
    return valor


def estimar_cochlea(datos: dict[str, Any]) -> dict[str, Any]:
    """Malla del corte sin resolver: es lo que la UI muestra antes de apretar el botón."""
    p = validar_cochlea(datos)
    t0 = time.perf_counter()
    seccion = _seccion(p)
    dt = time.perf_counter() - t0
    if seccion.n_nodes > MAX_NODOS:
        raise ValueError(
            f"la malla sale de {seccion.n_nodes} nodos; el tope interactivo es {MAX_NODOS}. "
            f"Bajá ppw (hoy {p['ppw']:g})."
        )
    return {
        "ok": True,
        "malla": {
            "n_nodos": int(seccion.n_nodes),
            "n_elementos": int(seccion.n_elements),
            "radio_mm": RADIO_M * 1e3,
            "cascara_mm": p["thickness_mm"],
        },
        "estimado_solve_s": round(0.5 + seccion.n_nodes / 20_000.0, 2),
        "tiempo_malla_s": round(dt, 3),
        "supuestos": {
            "c_p_hueso": HUESO_CP,
            "rho_hueso": HUESO_RHO,
            "frecuencia_hz": FRECUENCIA_HZ,
            "verificacion": "c_p y rho son ballpark VERIFY (D20); c_s, espesor y alpha los fija el "
                            "usuario",
        },
    }


def resolver_cochlea(datos: dict[str, Any]) -> dict[str, Any]:
    """Resuelve el corte y devuelve campo interior, balance de energía y las métricas de D29/D30."""
    t_inicio = time.perf_counter()
    p = validar_cochlea(datos)

    t0 = time.perf_counter()
    seccion = _seccion(p)
    t_malla = time.perf_counter() - t0
    if seccion.n_nodes > MAX_NODOS:
        raise ValueError(
            f"la malla sale de {seccion.n_nodes} nodos; el tope interactivo es {MAX_NODOS}."
        )

    kwargs: dict[str, Any] = {"alpha": p["alpha_np_m"]}
    if p["membrana"]:
        kwargs["membrane_impedance"] = cochlea.BASILAR_MEMBRANE_BALLPARK.impedance(FRECUENCIA_HZ)
    t0 = time.perf_counter()
    sol = cochlea.solve_cross_section(seccion, **kwargs)
    t_solve = time.perf_counter() - t0

    # ---- mapa del campo interior: se interpola sólo con los nodos de fluido ----
    nodos = np.asarray(seccion.nodes, dtype=float)
    p_nodal = np.asarray(sol.p, dtype=complex)
    radios = np.linalg.norm(nodos, axis=1)
    finito = np.isfinite(p_nodal)
    interior = finito & (radios < RADIO_M * (1.0 - 1e-9))

    a_mm = RADIO_M * 1e3
    eje = np.linspace(-a_mm, a_mm, PUNTOS_MAPA)
    gx, gy = np.meshgrid(eje, eje, indexing="xy")
    consulta = np.column_stack([gx.ravel(), gy.ravel()]) * 1e-3
    fuera = np.linalg.norm(consulta, axis=1) >= RADIO_M

    amplitud = np.full(consulta.shape[0], np.nan)
    if interior.sum() >= 10:
        interp = LinearNDInterpolator(
            nodos[interior],
            np.abs(p_nodal[interior]),
            fill_value=np.nan,
        )
        amplitud[~fuera] = interp(consulta[~fuera])
    amplitud = amplitud.reshape(PUNTOS_MAPA, PUNTOS_MAPA)

    pico = float(np.nanmax(amplitud)) if np.isfinite(amplitud).any() else 1.0
    pico = pico or 1.0
    with np.errstate(invalid="ignore", divide="ignore"):
        db = 20.0 * np.log10(np.maximum(amplitud / pico, 1e-6))
    # huecos como None (null): un NaN de Python no es JSON válido (D35)
    db_mapa = [
        [None if not np.isfinite(db[i, j]) else round(float(db[i, j]), 3)
         for j in range(PUNTOS_MAPA)]
        for i in range(PUNTOS_MAPA)
    ]

    metricas = _json_limpio(sol.as_dict()) if hasattr(sol, "as_dict") else {}
    # las cifras que la UI muestra primero, con nombre explícito
    resumen = {
        "t_interior": float(sol.t_interior),
        "t_interior_n0": float(sol.t_interior_complex),
        "absorcion_membrana": float(getattr(sol, "power_membrane_fraction", 0.0)),
        "ventana_real": float(getattr(sol, "power_window_real_fraction", 0.0)),
        "dispersada": float(getattr(sol, "power_scattered_fraction", 0.0)),
        "balance_total": float(getattr(sol, "power_total_fraction", 0.0)),
        "cavidad_neta": float(getattr(sol, "power_cavity_fraction", 0.0)),
    }
    # La identidad que documenta `cochlea.CrossSectionSolution` es `P_total + P_membrana = 0`
    # cuando la membrana es el único sumidero. NO es "dispersada + absorbida = 1": `dispersada`
    # (P_scattered/P_inc) puede pasar de 1 porque se mide en otro radio y con otra normalización.
    # Mezclarlas daría un "balance" de 1.85 que no significa nada; este sí cierra.
    resumen["balance_energia"] = float(
        abs(resumen["balance_total"] + resumen["absorcion_membrana"])
    )
    resumen["extincion"] = float(getattr(sol, "power_extinction_fraction", 0.0))

    # El criterio de t_interior sólo aplica al caso base de D30: con otros parámetros no hay
    # referencia contra la cual comparar, y eso se declara (sin_referencia), no se omite.
    es_caso_base = (
        abs(p["thickness_mm"] - 2.0) < 1e-9
        and abs(p["c_s"] - 1900.0) < 1e-9
        and abs(p["alpha_np_m"]) < 1e-9
        and p["membrana"] and p["ventana"]
    )
    sello = validacion.sello_conjunto({
        "cochlea_t_interior": resumen["t_interior"] if es_caso_base else None,
        "cochlea_balance": resumen["balance_energia"],
    })

    return {
        "ok": True,
        "sello": sello,
        "parametros": p,
        "mapa": {
            "x_mm": np.round(eje, 4).tolist(),
            "y_mm": np.round(eje, 4).tolist(),
            "db": db_mapa,
            "radio_mm": a_mm,
            "cascara_mm": p["thickness_mm"],
            "membrana_angulo_deg": 90.0 if p["membrana"] else None,
        },
        "metricas": resumen,
        "detalle": metricas,
        "malla": {
            "n_nodos": int(seccion.n_nodes),
            "n_elementos": int(seccion.n_elements),
            "nodos_fluido": int(interior.sum()),
        },
        "tiempos": {
            "malla_s": round(t_malla, 3),
            "solve_s": round(t_solve, 3),
            "total_s": round(time.perf_counter() - t_inicio, 3),
        },
        "referencia": {
            "t_interior_d30": 0.53187,
            "nota": "t_interior = 0.53187 es el valor de D30 con espesor 2 mm, c_s 1900 m/s, "
                    "alpha 0 y la membrana del ballpark",
        },
    }
