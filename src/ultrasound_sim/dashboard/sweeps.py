"""Barridos de sensibilidad de la cóclea como trabajo de cola (Fase U5).

Qué mide
--------
El ranking one-at-a-time del corte transversal: cuánto se mueve ``t_interior`` al variar **un**
parámetro por vez, con el resto en el caso base. Es el barrido de D24/D30 —los mismos rangos, los
mismos valores, la misma definición de "rango relativo" = ``(máx − mín) / t_interior(base)``—, así
que el resultado se puede contrastar contra números ya medidos:

    espesor 53.11 % · c_p 16.06 % · atenuación 9.60 % · c_s 8.22 % · ρ 1.03 % · Z_r 0.28 %

Esos son los pesos de D30 (con membrana). Si el barrido de la UI no los reproduce, hay un error acá
—no en el solver, que ya está validado— y el test lo fija.

Qué NO hace
-----------
- No corre la grilla 2D espesor × atenuación completa: eso son 52 puntos más y no cambia el ranking.
  Pedirla sigue siendo trabajo de `scripts/sensitivity_cochlea.py`.
- No mide Q. El ancho de banda de la resonancia necesita una definición fijada (¿−3 dB en |p| o en
  |p|²?, ¿sobre qué pico?) y una resolución de barrido que un Q de 287 exige por debajo del 0.35 %
  de la frecuencia. Queda explícito como pendiente en vez de inventar el número (D37).
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np

from .. import attenuation as att
from .. import cochlea, materials as mats

FRECUENCIA_HZ = 1.0e5
RADIO_M = 20.0e-3
PPW = 20.0
PAD_WAVELENGTHS = 2.5
MIN_ELEMENTS_ACROSS_WALL = 2
WINDOW_HALF_ANGLE_DEG = 10.0
WINDOW_CENTER_ANGLE_DEG = 0.0

BASE_CP = 3500.0
BASE_CS = 1900.0
BASE_RHO = 1950.0
BASE_ESPESOR_M = 2.0e-3

#: Puntos por parámetro. Los extremos tienen que ser los del barrido de D30: el "rango relativo" es
#: (máx − mín), así que recortar los extremos cambia el resultado.
PUNTOS = {
    # 13 puntos 1.0 … 7.0 mm cada 0.5: los del script. Con menos puntos se pierde un extremo
    # interior (t_interior NO es monótona en el espesor, D24/D29) y el rango sale 1.1 % corto.
    "espesor": [1.0e-3 + 0.5e-3 * i for i in range(13)],
    "c_p": [1752.0, 2500.0, 3285.0, 3500.0],
    "c_s": [1360.0, 1500.0, 1640.0, 1900.0],
    "rho": [1868.0, 1900.0, 1950.0],
    "atenuacion": ["sin_atenuacion", "white_2006", "pichardo_2011", "pinton_2012"],
    "Z_r": [0.5, 1.0, 2.0],   # los del script: sin impedancia citable de la ventana (D20)
}

#: Las tres leyes de D19 más el caso sin pérdidas. Los valores NO se escriben acá: se piden a
#: `attenuation.py` en la frecuencia del barrido (a 100 kHz: 3.777 / 6.243 / 15.312 Np/m). Escribir
#: un α "aproximado" en el barrido habría cambiado el peso de la atenuación sin que nada lo note.
LEYES_ALPHA = {
    "sin_atenuacion": None,
    "white_2006": att.SKULL_WHITE2006,
    "pichardo_2011": att.SKULL_PICHARDO_2011_POINTS,
    "pinton_2012": att.SKULL_PINTON_2012,
}


def alpha_de_ley(nombre: str) -> float:
    """α [Np/m] de la ley a la frecuencia del barrido; 0 si es el caso sin pérdidas."""
    ley = LEYES_ALPHA[nombre]
    return 0.0 if ley is None else float(ley.alpha(FRECUENCIA_HZ))


def _hueso(c_p: float, c_s: float, rho: float) -> mats.ElasticSolid:
    return mats.ElasticSolid(
        name=f"hueso_cp{c_p:.0f}_cs{c_s:.0f}_rho{rho:.0f}",
        rho=float(rho), c_p=float(c_p), c_s=float(c_s),
        provenance="punto del barrido de sensibilidad (D20/D30); c_s=1900 y c_p=3500 son el "
                   "ballpark VERIFY",
    )


def _t_interior(*, c_p: float, c_s: float, rho: float, espesor_m: float,
                alpha: float, z_factor: float) -> float:
    """Un punto del barrido: malla + solve con la membrana del ballpark puesta."""
    hueso = _hueso(c_p, c_s, rho)
    seccion = cochlea.build_cross_section(
        FRECUENCIA_HZ, RADIO_M, espesor_m,
        ppw=PPW, pad_wavelengths=PAD_WAVELENGTHS,
        min_elements_across_wall=MIN_ELEMENTS_ACROSS_WALL,
        bone=hueso,
        window_half_angle=WINDOW_HALF_ANGLE_DEG,
        window_center_angle=WINDOW_CENTER_ANGLE_DEG,
        window_model="impedance",
        window_impedance=complex(z_factor * mats.WATER.impedance),
        membrane_angle=90.0,
    )
    sol = cochlea.solve_cross_section(
        seccion, alpha=alpha,
        membrane_impedance=cochlea.BASILAR_MEMBRANE_BALLPARK.impedance(FRECUENCIA_HZ),
    )
    return float(sol.t_interior)


def total_de_puntos() -> int:
    return sum(len(v) for v in PUNTOS.values()) + 1  # +1 = el caso base


def barrido_sensibilidad(trabajo, **_ignorado: Any) -> dict[str, Any]:
    """Corre el OAT completo y devuelve el ranking con los rangos relativos.

    ``trabajo`` es el :class:`dashboard.jobs.Trabajo`: de ahí salen ``avanzar`` y ``cancelado``.
    """
    base = _t_interior(c_p=BASE_CP, c_s=BASE_CS, rho=BASE_RHO, espesor_m=BASE_ESPESOR_M,
                       alpha=0.0, z_factor=1.0)
    total = total_de_puntos()
    hechos = 1
    trabajo.avanzar(hechos / total, f"caso base: t_interior = {base:.5f}")

    por_parametro: dict[str, list[float]] = {}
    filas: list[dict[str, Any]] = []

    for nombre, valores in PUNTOS.items():
        for valor in valores:
            if trabajo.cancelado():
                raise RuntimeError("cancelado por el usuario")
            c_p, c_s, rho = BASE_CP, BASE_CS, BASE_RHO
            espesor_m, alpha, z_factor = BASE_ESPESOR_M, 0.0, 1.0
            etiqueta = str(valor)

            if nombre == "espesor":
                espesor_m = float(valor)  # type: ignore[arg-type]
                etiqueta = f"{espesor_m * 1e3:g} mm"
            elif nombre == "c_p":
                c_p = float(valor)  # type: ignore[arg-type]
                # con c_p = 1752 no es admisible c_s = 1900 (ν > 0.5): se baja a 1500
                c_s = 1500.0
                etiqueta = f"{c_p:g} m/s"
            elif nombre == "c_s":
                c_s = float(valor)  # type: ignore[arg-type]
                etiqueta = f"{c_s:g} m/s"
            elif nombre == "rho":
                rho = float(valor)  # type: ignore[arg-type]
                etiqueta = f"{rho:g} kg/m³"
            elif nombre == "atenuacion":
                alpha = alpha_de_ley(str(valor))
                etiqueta = f"{valor} (α = {alpha:.2f} Np/m)"
            elif nombre == "Z_r":
                z_factor = float(valor)  # type: ignore[arg-type]
                etiqueta = f"{z_factor:g} · ρc"

            valor_t = _t_interior(c_p=c_p, c_s=c_s, rho=rho, espesor_m=espesor_m,
                                  alpha=alpha, z_factor=z_factor)
            por_parametro.setdefault(nombre, []).append(valor_t)
            filas.append({
                "parametro": nombre,
                "valor": etiqueta,
                "t_interior": valor_t,
                "delta_pct": 100.0 * (valor_t / base - 1.0),
            })
            hechos += 1
            trabajo.avanzar(hechos / total, f"{nombre} = {etiqueta}: {valor_t:.5f}")

    ranking = []
    for nombre, valores in por_parametro.items():
        lo, hi = min(valores), max(valores)
        ranking.append({
            "parametro": nombre,
            "min": lo,
            "max": hi,
            "rango_pct": 100.0 * (hi - lo) / base,
        })
    ranking.sort(key=lambda r: r["rango_pct"], reverse=True)

    return resumen_para_ui({
        "ok": True,
        "base": {"t_interior": base, "espesor_mm": BASE_ESPESOR_M * 1e3, "c_p": BASE_CP,
                 "c_s": BASE_CS, "rho": BASE_RHO},
        "puntos": len(filas),
        "filas": filas,
        "ranking": ranking,
        # referencia ya medida (D30, con membrana): el contraste que la UI muestra
        "referencia_d30": {"espesor": 53.11, "c_p": 16.06, "atenuacion": 9.60,
                           "c_s": 8.22, "rho": 1.03, "Z_r": 0.28},
    })


def enviar_barrido(cola, etiqueta: str = "sensibilidad de la cóclea"):
    """Encola el barrido y devuelve el job (la UI hace polling de ``GET /api/jobs/<id>``)."""
    return cola.enviar("sensibilidad", etiqueta, barrido_sensibilidad)


def resumen_para_ui(resultado: dict[str, Any]) -> dict[str, Any]:
    """Compara los rangos medidos contra la referencia de D30, parámetro por parámetro."""
    ref = resultado["referencia_d30"]
    comparacion = []
    for fila in resultado["ranking"]:
        esperado = ref.get(fila["parametro"])
        desvio = None if esperado is None else 100.0 * (fila["rango_pct"] - esperado) / esperado
        comparacion.append({**fila, "d30_pct": esperado, "desvio_pct": desvio})
    return {**resultado, "comparacion": comparacion}


def estimar_demora_s() -> float:
    """Estimación honesta: cada punto es una malla (~0.5 s) y un solve (~1.5–2.5 s)."""
    return round(total_de_puntos() * 2.4, 0)
