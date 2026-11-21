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
from . import validacion
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
#: Pesos publicados en D30 (con membrana), en % del `t_interior` base. Es la referencia contra la
#: que se sella el barrido; vive acá y no en el JS (U6).
REFERENCIA_D30 = {"espesor": 53.11, "c_p": 16.06, "atenuacion": 9.60,
                  "c_s": 8.22, "rho": 1.03, "Z_r": 0.28}

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


#: Materiales con su `c` para predecir las resonancias de la losa.
MATERIALES_COCLEA = {
    "agua": mats.WATER, "tejido": mats.SOFT_TISSUE, "hueso": mats.BONE_CORTICAL,
    "perilinfa": mats.PERILYMPH, "endolinfa": mats.ENDOLYMPH,
}


def _hueso(c_p: float, c_s: float, rho: float) -> mats.ElasticSolid:
    return mats.ElasticSolid(
        name=f"hueso_cp{c_p:.0f}_cs{c_s:.0f}_rho{rho:.0f}",
        rho=float(rho), c_p=float(c_p), c_s=float(c_s),
        provenance="punto del barrido de sensibilidad (D20/D30); c_s=1900 y c_p=3500 son el "
                   "ballpark VERIFY",
    )


def _t_interior(*, c_p: float, c_s: float, rho: float, espesor_m: float,
                alpha: float, z_factor: float, frecuencia: float | None = None) -> float:
    """Un punto del barrido: malla + solve con la membrana del ballpark puesta."""
    f = FRECUENCIA_HZ if frecuencia is None else float(frecuencia)
    hueso = _hueso(c_p, c_s, rho)
    seccion = cochlea.build_cross_section(
        f, RADIO_M, espesor_m,
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
        membrane_impedance=cochlea.BASILAR_MEMBRANE_BALLPARK.impedance(f),
    )
    return float(sol.t_interior)


# --------------------------------------------------------------------------- barrido de frecuencia y Q
#
# DEFINICIÓN DE Q — fijada por Sim-Datos en la consulta del 2026-11-20, y es **operativa**, no un dato
# biológico: no existe un Q citable de la resonancia de la cápsula ótica humana en 0.1–1 kHz ni a
# 100 kHz–1 MHz. Se define, para un modo aislado y débilmente amortiguado:
#
#     Q_3dB = f_peak / (f2 − f1),  con |H|²(f) = |H|²_max / 2,  o sea |H| = |H|_max / √2
#
# Son los mismos puntos se use 20·log10|H| o 10·log10|H|²; **no** es "media amplitud". Y en un campo
# estacionario |p|² no representa por sí solo la intensidad neta, así que lo que se mide es el ancho de
# banda del pico de `t_interior`, no una potencia absorbida.
#
# Lo que sí está documentado y NO es esto: Q_ERB = CF/ERB ≈ 12.7 ± 1 a 1 kHz (Shera, Guinan y Oxenham
# 2002) es selectividad de los filtros cocleares, no la resonancia global de la cápsula, y Q_10dB se
# usa para respuestas de membrana basilar dependientes del nivel. Nada de eso se traslada acá.

BANDA_Q_POR_DEFECTO = (85.0e3, 115.0e3)   # Hz
PUNTOS_Q_POR_DEFECTO = 13


def _campana_decaimiento(valores: np.ndarray, frecuencias: np.ndarray) -> dict[str, Any]:
    """Q por ancho de banda a −3 dB en potencia (|H| = máx/√2) sobre el pico más alto."""
    valores = np.asarray(valores, dtype=float)
    if valores.size < 5 or not np.all(np.isfinite(valores)):
        return {"Q": None, "motivo": "curva incompleta para medir el ancho de banda"}

    i_pico = int(np.argmax(valores))
    pico = float(valores[i_pico])
    nivel = pico / np.sqrt(2.0)

    def _cruce(inicio: int, paso: int) -> float | None:
        """Camina desde el pico hacia un extremo e interpola dónde corta el nivel de −3 dB.

        Ojo con el rango: `range(i, -1)` con paso +1 es vacío (bug que tuve en la primera versión y
        dejaba el Q siempre en `None`). Por eso se avanza a mano con el paso explícito.
        """
        i = inicio
        while 0 <= i + paso < valores.size:
            j = i + paso
            a, b = valores[i], valores[j]
            if (a - nivel) * (b - nivel) <= 0 and a != b:
                peso = (nivel - a) / (b - a)
                return float(frecuencias[i] + peso * (frecuencias[j] - frecuencias[i]))
            i = j
        return None

    f_abajo = _cruce(i_pico, -1)
    f_arriba = _cruce(i_pico, +1)
    if f_abajo is None or f_arriba is None:
        # Si la curva no baja 3 dB dentro de la banda, el ancho a media potencia es AL MENOS el de la
        # banda: eso da una cota superior de Q, que es un resultado (la respuesta es ancha), no un
        # fracaso. Es el caso del corte con membrana: no hay un pico estrecho que medir.
        ancho_banda = float(frecuencias[-1] - frecuencias[0])
        return {
            "Q": None,
            "Q_cota_superior": float(frecuencias[i_pico] / ancho_banda),
            "f_pico_hz": float(frecuencias[i_pico]),
            "motivo": ("la respuesta no cae 3 dB dentro de la banda: el ancho a media potencia es "
                       "al menos el de la banda, así que Q no supera la cota superior"),
        }

    ancho = f_arriba - f_abajo
    if ancho <= 0:
        return {"Q": None, "motivo": "los cruces de −3 dB salieron invertidos"}
    return {
        "Q": float(frecuencias[i_pico] / ancho),
        "f_pico_hz": float(frecuencias[i_pico]),
        "f_abajo_hz": f_abajo,
        "f_arriba_hz": f_arriba,
        "ancho_hz": float(ancho),
        "nivel_3db": float(nivel),
        "resolucion_relativa": float(abs(frecuencias[1] - frecuencias[0]) / frecuencias[i_pico]),
    }


def barrido_frecuencia(trabajo, **datos: Any) -> dict[str, Any]:
    """t_interior(f) en una banda y el Q que sale de su pico, con la definición de arriba."""
    desde = float(datos.get("desde_hz", BANDA_Q_POR_DEFECTO[0]))
    hasta = float(datos.get("hasta_hz", BANDA_Q_POR_DEFECTO[1]))
    puntos = int(datos.get("puntos", PUNTOS_Q_POR_DEFECTO))
    if not (20.0e3 <= desde < hasta <= 1.0e6):
        raise ValueError(f"banda inválida: {desde:g}–{hasta:g} Hz")
    if not (5 <= puntos <= 41):
        raise ValueError(f"puntos fuera de rango [5, 41]: {puntos}")

    frecuencias = np.linspace(desde, hasta, puntos)
    valores: list[float] = []
    for i, f in enumerate(frecuencias):
        if trabajo.cancelado():
            raise RuntimeError("cancelado por el usuario")
        trabajo.avanzar((i + 1) / (puntos + 1), f"f = {f / 1e3:.1f} kHz")
        valores.append(_t_interior(c_p=BASE_CP, c_s=BASE_CS, rho=BASE_RHO,
                                   espesor_m=BASE_ESPESOR_M, alpha=0.0, z_factor=1.0,
                                   frecuencia=float(f)))

    arreglo = np.asarray(valores)
    q = _campana_decaimiento(arreglo, frecuencias)
    return {
        "ok": True,
        "banda_hz": [desde, hasta],
        "puntos": puntos,
        "frecuencias_hz": frecuencias.tolist(),
        "t_interior": arreglo.tolist(),
        "q": q,
        "definicion": "Q = f_pico/(f2 − f1) con |p|² = |p|²_max/2 (o sea |p| = |p|_max/√2)",
        "advertencia": (
            "Definición OPERATIVA: no existe un Q citable de la resonancia de la cápsula ótica. "
            "Q_ERB ≈ 12.7 a 1 kHz es selectividad de los filtros cocleares, no esta resonancia."
        ),
        "resolucion": (
            f"con {puntos} puntos el paso relativo es "
            f"{100 * abs(frecuencias[1] - frecuencias[0]) / frecuencias[int(np.argmax(arreglo))]:.2f} % "
            f"del pico; un Q mayor que ~1/paso no se resuelve, hay que refinar la banda"
        ),
    }


# --------------------------------------------------------------------------- barrido de la escena 2D
#
# La escena 2D resuelve en ~10–80 ms, así que un barrido de 40 frecuencias son un par de segundos:
# entra en un job por consistencia, no por necesidad.
#
# Criterio físico verificable: una losa de espesor t en un medio tiene MÁXIMOS DE TRANSMISIÓN en las
# resonancias de media onda, f_n = n·c_losa/(2t). Con hueso como fluido (c = 3500 m/s) y t = 10 mm
# eso da 175, 350, 525 kHz… Si el barrido no los encuentra, hay un error de espesor, de c o de fase.

def barrido_apertura(datos: dict, trabajo=None) -> dict[str, Any]:
    """FWHM lateral contra el paso entre elementos: la pregunta de diseño del proyecto.

    Con la apertura total fija (`d`), variar la cantidad de elementos `n` es variar el **pitch**
    (`d/n`). Es la variable que D25/D26 identifican como la que manda: para empatar con una apertura
    continua dentro del ~2 % hace falta `pitch ≤ λ/4`, y lo que pesa es la cantidad de elementos, no
    el ancho de cada uno.

    Se resuelve un campo por punto (no un barrido de frecuencia): cada punto es un mapa de verdad, así
    que 10 puntos son ~10 s con la grilla de la vista.
    """
    # import local: `server` importa este módulo, así que a nivel de módulo sería un ciclo
    from . import server

    material = str(datos.get("material", "agua"))
    f_khz = float(datos.get("frecuencia_khz", 250.0))
    apertura_mm = float(datos.get("apertura_mm", 60.0))
    foco_mm = float(datos.get("foco_mm", 90.0))
    modelo = str(datos.get("modelo", "baffle"))
    desde = int(datos.get("desde", 8))
    hasta = int(datos.get("hasta", 64))
    puntos = int(datos.get("puntos", 7))
    if not 2 <= desde < hasta <= 512:
        raise ValueError(f"pide 2 <= desde < hasta <= 512; llegó {desde}, {hasta}")
    if not 2 <= puntos <= 15:
        raise ValueError(f"puntos tiene que estar en [2, 15]: {puntos}")
    if not 5.0 <= apertura_mm <= 200.0:
        raise ValueError(f"apertura_mm fuera de [5, 200]: {apertura_mm}")

    valores = sorted({int(round(v)) for v in np.linspace(desde, hasta, puntos)})

    def un_punto(n: int) -> dict[str, Any]:
        d_m = apertura_mm * 1e-3
        # los nombres son los del endpoint (`n_elements`, `aperture_mm`, `focus_x_mm`): usar `n` y `d`
        # no fallaba, simplemente se ignoraba y el barrido devolvía el MISMO campo en todos los
        # puntos (FWHM 9.12 mm con 8 y con 64 elementos, sidelobe -9.24 dB idéntico). Un barrido que
        # no barre es peor que no tenerlo, porque parece un resultado.
        r = server.calcular({
            "material": material, "modelo": modelo, "frequency_khz": f_khz,
            "n_elements": n, "aperture_mm": apertura_mm,
            "focus_x_mm": foco_mm, "focus_y_mm": 0.0,
            "grid": {"x0": 0.005, "x1": 0.001 * (foco_mm * 2.2), "y0": -0.001 * apertura_mm * 0.9,
                     "y1": 0.001 * apertura_mm * 0.9, "nx": 120, "ny": 120},
        })
        m = r["metrics"]
        pitch = d_m / n
        # λ del MATERIAL elegido, no del agua: en hueso (c = 3500) la λ es más de dos veces más corta
        # y el pitch que cumple λ/4 baja con ella. Fijarla en 1500 sería medir el criterio de otro medio.
        lam = float(r["meta"]["lambda_mm"]) * 1e-3
        return {
            "n": n, "pitch_mm": pitch * 1e3, "pitch_over_lambda": pitch / lam,
            "fwhm_lateral_mm": m.get("fwhm_lateral_mm"),
            "fwhm_over_lambda": (m.get("fwhm_lateral_mm") / 1e3 / lam) if m.get("fwhm_lateral_mm") else None,
            "sidelobe_db": m.get("sidelobe_db"), "grating_deg": m.get("grating_deg"),
            "peak_x_mm": m.get("peak_x_mm"), "segundos": r["meta"]["seconds"],
        }

    filas = []
    for i, n in enumerate(valores):
        if trabajo is not None and trabajo.cancelado():
            raise RuntimeError("cancelado por el usuario")
        if trabajo is not None:
            trabajo.avanzar(i / (len(valores) + 1), f"n = {n} elementos (pitch {apertura_mm/n:.2f} mm)")
        filas.append(un_punto(n))

    # referencia cuasi-continua: mismo ancho, muchos elementos. Es el "empata con continuo" de D25.
    if trabajo is not None:
        trabajo.avanzar(len(valores) / (len(valores) + 1), "apertura cuasi-continua de referencia")
    ref = un_punto(400)

    peor = None
    if ref.get("fwhm_lateral_mm"):
        for f in filas:
            if f.get("fwhm_lateral_mm"):
                f["desvio_vs_continuo"] = abs(f["fwhm_lateral_mm"] - ref["fwhm_lateral_mm"]) / ref["fwhm_lateral_mm"]
                if peor is None or f["desvio_vs_continuo"] > peor["desvio_vs_continuo"]:
                    peor = f

    return {
        "ok": True,
        "barrido": {
            "material": material, "frecuencia_khz": f_khz, "apertura_mm": apertura_mm,
            "foco_mm": foco_mm, "modelo": modelo,
            "n": [f["n"] for f in filas],
            "pitch_mm": [f["pitch_mm"] for f in filas],
            "pitch_over_lambda": [f["pitch_over_lambda"] for f in filas],
            "fwhm_lateral_mm": [f["fwhm_lateral_mm"] for f in filas],
            "desvio_vs_continuo": [f.get("desvio_vs_continuo") for f in filas],
            "sidelobe_db": [f["sidelobe_db"] for f in filas],
            "grating_deg": [f["grating_deg"] for f in filas],
        },
        "filas": filas,
        "referencia_continua": ref,
        "resumen": {
            "n_puntos": len(filas),
            "peor_desvio": None if peor is None else float(peor["desvio_vs_continuo"]),
            "peor_n": None if peor is None else int(peor["n"]),
            "pitch_con_2pct": next((f["pitch_over_lambda"] for f in filas
                                    if f.get("desvio_vs_continuo") is not None
                                    and f["desvio_vs_continuo"] <= 0.02), None),
            "segundos_totales": float(sum(f["segundos"] for f in filas) + ref["segundos"]),
        },
    }


def _frecuencias(desde: float, hasta: float, puntos: int) -> np.ndarray:
    if not (20.0e3 <= desde < hasta <= 1.0e6):
        raise ValueError(f"banda inválida: {desde:g}–{hasta:g} Hz")
    if not (3 <= puntos <= 121):
        raise ValueError(f"puntos fuera de rango [3, 121]: {puntos}")
    return np.linspace(desde, hasta, puntos)


def _picos(valores: np.ndarray, frecuencias: np.ndarray) -> list[dict[str, float]]:
    """Máximos locales (interpolados parabólicamente) con su frecuencia."""
    picos: list[dict[str, float]] = []
    for i in range(1, valores.size - 1):
        if not (valores[i] > valores[i - 1] and valores[i] >= valores[i + 1]):
            continue
        a, b, c = valores[i - 1], valores[i], valores[i + 1]
        denom = a - 2.0 * b + c
        delta = 0.5 * (a - c) / denom if denom != 0.0 else 0.0
        f_pico = frecuencias[i] + delta * (frecuencias[i + 1] - frecuencias[i])
        picos.append({"f_hz": float(f_pico), "valor": float(b - 0.25 * (a - c) * delta)})
    return picos


def barrido_escena(trabajo, **datos: Any) -> dict[str, Any]:
    """R y T de la escena 2D en función de la frecuencia, con los picos de transmisión."""
    from . import scene as escena_mod

    base = datos.get("escena") or {}
    if not isinstance(base, dict):
        raise ValueError("'escena' tiene que ser el objeto de la escena 2D")
    desde = float(datos.get("desde_hz", 100.0e3))
    hasta = float(datos.get("hasta_hz", 400.0e3))
    puntos = int(datos.get("puntos", 30))
    frecuencias = _frecuencias(desde, hasta, puntos)

    filas: list[dict[str, Any]] = []
    for i, f in enumerate(frecuencias):
        if trabajo.cancelado():
            raise RuntimeError("cancelado por el usuario")
        trabajo.avanzar((i + 0.5) / puntos, f"f = {f / 1e3:.0f} kHz")
        try:
            r = escena_mod.resolver_escena(dict(base, frecuencia_hz=float(f)))
        except ValueError as exc:
            filas.append({"f_hz": float(f), "R": None, "T": None, "error": str(exc)})
            continue
        m = r["metricas"]
        filas.append({
            "f_hz": float(f),
            "R": m.get("R"),
            "T": m.get("T"),
            "R_analitico": m.get("R_analitico"),
            "T_analitico": m.get("T_analitico"),
            "error_R_abs": m.get("error_R_abs"),
            "error_T_abs": m.get("error_T_abs"),
            "elementos": r["grilla"]["n_elementos"],
        })
        trabajo.avanzar((i + 1) / puntos, f"f = {f / 1e3:.0f} kHz: T = {m.get('T') or 0:.3f}")

    validas = [fila for fila in filas if fila["T"] is not None]
    if not validas:
        raise ValueError("ninguna frecuencia del barrido dio R y T (¿banda fuera de rango?)")

    fs = np.array([fila["f_hz"] for fila in validas])
    ts = np.array([float(fila["T"]) for fila in validas])
    picos = _picos(ts, fs)

    # Si la escena es una losa en un medio, los picos deberían caer en n·c/(2t)
    capas = base.get("capas") or []
    prediccion: dict[str, Any] | None = None
    if len(capas) == 3:
        material_medio = str(capas[1].get("material", ""))
        espesor = float(capas[1].get("espesor_m", 0.0))
        if material_medio in MATERIALES_COCLEA and espesor > 0:
            c_losa = MATERIALES_COCLEA[material_medio].c
            f1 = c_losa / (2.0 * espesor)
            prediccion = {
                "material": material_medio, "c_m_s": c_losa, "espesor_mm": espesor * 1e3,
                "resonancias_hz": [n * f1 for n in range(1, 6) if n * f1 <= fs[-1] * 1.05],
                "formula": "f_n = n·c_losa/(2·t) — resonancias de media onda de la losa",
            }

    errores = [abs(fila["error_T_abs"]) for fila in validas if fila["error_T_abs"] is not None]
    return {
        "ok": True,
        "banda_hz": [desde, hasta],
        "puntos": len(filas),
        "filas": filas,
        "picos": picos,
        "prediccion": prediccion,
        "max_error_T": max(errores) if errores else None,
        "sello": validacion.sello_conjunto(
            # el peor error de T del barrido contra la matriz de transferencia
            {"escena_rt_t": max(errores) if errores else None}
        ),
    }


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

    desvios = [
        abs(f["rango_pct"] - REFERENCIA_D30[f["parametro"]]) / REFERENCIA_D30[f["parametro"]]
        for f in ranking if f["parametro"] in REFERENCIA_D30
    ]
    sello = validacion.sello_conjunto(
        {"barrido_d30": max(desvios) if desvios else None}
    )

    return resumen_para_ui({
        "ok": True,
        "base": {"t_interior": base, "espesor_mm": BASE_ESPESOR_M * 1e3, "c_p": BASE_CP,
                 "c_s": BASE_CS, "rho": BASE_RHO},
        "puntos": len(filas),
        "filas": filas,
        "ranking": ranking,
        # referencia ya medida (D30, con membrana): el contraste que la UI muestra
        "referencia_d30": dict(REFERENCIA_D30),
        "sello": sello,
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
