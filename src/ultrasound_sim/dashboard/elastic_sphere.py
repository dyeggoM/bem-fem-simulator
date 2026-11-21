"""Esfera **elástica** sumergida resuelta por el FEM acoplado, contra la serie de Faran (U7).

Qué agrega sobre la esfera penetrable 3D (D42)
----------------------------------------------
Aquella es un dispersor de **contraste de impedancia**: adentro sigue siendo un fluido. Ésta resuelve
el sólido de verdad, con dos potenciales (P y S) y las tres condiciones de interfaz, así que aparecen
las **ondas de corte** y la respuesta depende de `c_s` — que es la variable que importa para el hueso.

Las dos piezas ya estaban y cada una tiene su validación:

- `analytic/faran.py`: la serie exacta (D45), validada por unitariedad a 1e-16 y teorema óptico.
- `fem/coupled3d.py`: el FEM acoplado (D46), validado por límites y por el test estático exacto.

Acá sólo se juntan y se sellan.

Los dos modos de malla, y por qué el rápido no está validado
-----------------------------------------------------------
El costo manda: con dominio de `4a` y `h = a/4.5` son 21 398 nodos y ~32 s; con `2.5a` y `h = a/4` es
~1 s. El dominio chico deja ~9 % de error porque el ABC-2 no alcanza con `k₀b = 2.5` (D46, 3.er
apéndice). Por eso el modo rápido existe **etiquetado**, y su sello sale `fuera` — no se disimula: es
una vista para explorar, no un número para reportar.
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np

from ..analytic import faran
from ..fem import coupled3d
from ..materials import WATER, ElasticSolid
from ..meshing3d_sphere import sphere_shell_tet_mesh
from . import validacion

__all__ = ["MODOS", "MATERIALES", "costo", "resolver"]

#: Modos de malla: ``(b/a, h/a)``. El "precisa" es el que se reporta.
MODOS: dict[str, dict[str, Any]] = {
    "rapida": {"b_rel": 2.5, "h_rel": 4.0, "etiqueta": "rápida (dominio 2.5a, ~1 s)"},
    "precisa": {"b_rel": 4.0, "h_rel": 4.5, "etiqueta": "precisa (dominio 4a, ~30 s)"},
}

#: Materiales con `c_s`, que es lo que distingue este problema del de contraste de impedancia.
MATERIALES = {
    "acero": {"rho": 7800.0, "c_p": 5900.0, "c_s": 3200.0},
    "hueso": {"rho": 1950.0, "c_p": 3500.0, "c_s": 1900.0},
    "aluminio": {"rho": 2700.0, "c_p": 6320.0, "c_s": 3130.0},
    "tungsteno": {"rho": 19300.0, "c_p": 5200.0, "c_s": 2800.0},
    "tejido (casi fluido)": {"rho": 1050.0, "c_p": 1540.0, "c_s": 30.0},
}

_LIMITES = {"frecuencia_hz": (20e3, 400e3), "ka": (0.3, 3.0), "puntos_angulo": (12, 361)}

#: **Interruptor empírico**, no un criterio físico. Debajo de este contraste de impedancia (`ρ·c_p`
#: del sólido sobre `ρ·c` del agua) la restricción de traslación media nula rompe el resultado: medido
#: 4.55 (hueso) da 13 %, 1.08 (tejido) da **5264 %** con la restricción y 23 % sin ella. El 2.5 sale
#: de esos dos puntos.
#:
#: Sim-Datos (D48) confirmó que **no existe** un umbral citable: la magnitud que gobierna el error es
#: la respuesta dipolar `3/(2g+1)` con `g = ρ_s/ρ_0` (Bruus 2012, ecs. 22 y 26), no la impedancia. Se
#: deja como interruptor de seguridad mientras se reemplaza la formulación.
CONTRASTE_MINIMO = 2.5

#: Las dos formulaciones del sólido, con el número medido que las distingue. No es una preferencia:
#: la dispersada es la correcta en el continuo (da la traslación `3/(2g+1)` al 1-2 % y el caso
#: transparente exacto) pero es una resta de términos grandes y **contra Faran es peor** (19.49 % en
#: acero contra 4.29 %). La total+restricción es incorrecta pero numéricamente cómoda. D49.
FORMULACIONES: dict[str, dict[str, Any]] = {
    "total": {
        "etiqueta": "total + traslación nula (la sellada: 4.29 % en acero)",
        "nota": "Suprime el dipolo físico (D48); el 4 % está en parte cancelado, no ganado.",
    },
    "dispersada": {
        "etiqueta": "dispersada (la correcta: traslación exacta, 19.49 % en acero)",
        "nota": "Correcta en el continuo y mal condicionada en esta malla (D49).",
    },
}

#: Debajo de este `|p_s|` máximo (relativo a la incidente, que vale 1) el error *relativo* deja de ser
#: un indicador útil: el denominador se achica pero el error absoluto no. Ahí el sello pasa a
#: `sin_referencia` y se informa el error absoluto, en vez de reportar "23 %" de algo que vale 1/50.
#:
#: Sim-Datos (D48) pide reportar **siempre las dos** métricas (`E_inc` y `E_sca`) y avisa que no hay
#: umbral citable para descartar la relativa: cerca de un nulo del campo es muy sensible (Lavia et al.,
#: TetraScatt, fig. 4). Por eso las dos van en la salida y este corte sólo decide **el sello**.
AMPLITUD_MINIMA = 0.05


def _frecuencia(datos: dict) -> float:
    f = float(datos.get("frecuencia_hz", 100e3))
    lo, hi = _LIMITES["frecuencia_hz"]
    if not (lo <= f <= hi):
        raise ValueError(f"frecuencia_hz fuera de [{lo:g}, {hi:g}] Hz: {f:g}")
    return f


def _modo(datos: dict) -> tuple[str, dict]:
    nombre = str(datos.get("modo", "precisa"))
    if nombre not in MODOS:
        raise ValueError(f"modo desconocido {nombre!r}; hay {sorted(MODOS)}")
    return nombre, MODOS[nombre]


def _ka(datos: dict) -> float:
    ka = float(datos.get("ka", 1.0))
    lo, hi = _LIMITES["ka"]
    if not (lo <= ka <= hi):
        raise ValueError(f"ka fuera de [{lo}, {hi}]: {ka}")
    return ka


def _formulacion(datos: dict) -> str:
    nombre = str(datos.get("formulacion", "total"))
    if nombre not in FORMULACIONES:
        raise ValueError(f"formulacion desconocida {nombre!r}; hay {sorted(FORMULACIONES)}")
    return nombre


def _material(datos: dict) -> dict[str, float]:
    nombre = str(datos.get("material", "acero"))
    if nombre not in MATERIALES:
        raise ValueError(f"material desconocido {nombre!r}; hay {sorted(MATERIALES)}")
    return MATERIALES[nombre]


def costo(datos: dict | None = None) -> dict[str, Any]:
    """Costo estimado **antes** de resolver, para que la UI no lance a ciegas."""
    datos = datos or {}
    _, modo = _modo(datos)
    ka = _ka(datos)
    k0 = abs(complex(WATER.wavenumber(_frecuencia({"frecuencia_hz": datos.get("frecuencia_hz", 100e3)}))))
    radio = ka / k0
    b = modo["b_rel"] * radio
    h = radio / modo["h_rel"]
    # la cáscara domina el conteo: volumen ≈ (4/3)π(b³ − a³), unos 6 tets por cubo de lado h
    volumen = 4.0 / 3.0 * np.pi * (b**3 - radio**3)
    tets = 6.0 * volumen / h**3
    nodos = tets / 5.0
    segundos = 32.0 if modo["b_rel"] >= 4.0 else 1.3
    return {
        "nodos_estimados": int(nodos), "tets_estimados": int(tets), "segundos_estimados": segundos,
        "radio_m": radio, "dominio_m": b, "elemento_m": h, "modo": datos.get("modo", "precisa"),
        "aviso": (
            None if modo["b_rel"] >= 3.0 else
            "dominio chico: el ABC-2 deja ~9 % de error (D46). Sirve para explorar, no para reportar."
        ),
    }


def resolver(datos: dict, trabajo=None) -> dict[str, Any]:
    """Resuelve la esfera elástica y la contrasta contra Faran en la misma corrida.

    ``trabajo`` es el job de la cola (el worker lo pasa como primer argumento): se le informa el
    progreso con ``avanzar``, igual que el resto de los barridos.
    """
    def paso(fraccion: float, mensaje: str) -> None:
        if trabajo is not None:
            trabajo.avanzar(fraccion, mensaje)

    t0 = time.time()
    f = _frecuencia(datos)
    ka = _ka(datos)
    # OJO con los nombres: `_modo` devuelve el modo de MALLA ("rapida"/"precisa"), no el material.
    # Tenerlos cruzados hacía que la ficha dijera "rapida" donde va el material, y que el aviso de
    # bajo contraste no se disparara nunca.
    nombre_modo, modo = _modo(datos)
    nombre_material = str(datos.get("material", "acero"))
    prop = _material(datos)
    nombre_formulacion = _formulacion(datos)
    # en la dispersada no hay nada que restringir: el desplazamiento medio de `u_s` ya es el de la
    # física, y forzarlo suprimiría justo el dipolo que la formulación recupera
    fijar = bool(datos.get("fijar_rigidos", True)) if nombre_formulacion == "total" else False
    automatico = False

    agua = WATER
    k0 = abs(complex(agua.wavenumber(f)))
    radio = ka / k0
    solido = ElasticSolid(
        name=nombre_material, rho=prop["rho"], c_p=prop["c_p"], c_s=prop["c_s"],
        provenance=(
            f"benchmark de la UI (U7): {nombre_material}"
            + ("" if prop["c_s"] > 500 else
               " — c_s SIN FUENTE CITABLE en la banda 100 kHz-1 MHz (D48): supuesto ballpark,"
               " no dato")
        ),
    )

    # La restricción de traslación media nula supone que la onda plana no ejerce fuerza neta sobre
    # el cuerpo, y eso vale para un cuerpo que dispersa de verdad. Con contraste bajo el cuerpo se
    # mueve parcialmente CON el fluido y la suposición se rompe: medido, con impedancia 1.08× el agua
    # el campo sale **53×** el exacto (5264 % de error) con la restricción, y 1.13× sin ella (23 %).
    # Se desactiva sola y se informa, en vez de dejar que el usuario pida un número que no existe.
    contraste = (prop["rho"] * prop["c_p"]) / (agua.rho * float(agua.c))
    if nombre_formulacion == "total" and fijar and contraste < CONTRASTE_MINIMO:
        fijar = False
        automatico = True

    paso(0.05, f"malla: esfera de {radio*1e3:.2f} mm, dominio {modo['b_rel']:g}·a")
    mall = sphere_shell_tet_mesh(radio, modo["b_rel"] * radio, radio / modo["h_rel"])
    paso(0.35, f"resolviendo el acoplado ({mall.n_nodes} nodos)")

    try:
        res = coupled3d.resolver(
            mall, agua, solido, f, abc_order=2,
            fijar_rigidos=(fijar if nombre_formulacion == "total" else None),
            formulacion=nombre_formulacion,
        )
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(f"el solve acoplado falló: {exc}") from exc

    paso(0.8, "comparando contra Faran")
    n_ang = int(np.clip(datos.get("puntos_angulo", 121), *_LIMITES["puntos_angulo"]))
    th = np.linspace(0.0, np.pi, n_ang)
    r_muestra = 1.5 * radio
    puntos = r_muestra * np.stack([np.sin(th), np.zeros_like(th), np.cos(th)], axis=1)

    coef = faran.coeficientes(f, agua, solido, radio)
    exacto = faran.campo_dispersado(coef, puntos, radio)
    fem = coupled3d._interpolar(mall, res["p_s"], puntos)

    util = np.abs(exacto) > 1e-3 * np.abs(exacto).max()
    # Dos métricas del mismo campo, y no son intercambiables. El MÁXIMO es el que sella (criterio),
    # pero es frágil: la esfera del FEM es facetada y eso corre los lóbulos unos grados, lo que infla
    # el máximo sin que la amplitud esté mal (D49, apéndice 3: refinar 2.2× sólo bajó 19.49 → 18.84 %).
    # La L2 pesa toda la curva y no depende de un punto. Se reportan las dos.
    error = float(np.abs(fem[util] - exacto[util]).max() / np.abs(exacto[util]).max())
    delta = np.abs(fem[util] - exacto[util])
    error_l2 = float(np.sqrt(np.sum(delta**2) / np.sum(np.abs(exacto[util]) ** 2)))
    razon = float(np.abs(fem).max() / np.abs(exacto).max())
    # error absoluto, con la incidente (1 Pa) como escala: es el que sigue teniendo sentido cuando el
    # campo dispersado es muy débil
    error_abs = float(np.abs(fem[util] - exacto[util]).max())

    # Resolución de la onda de CORTE: la malla se elige con la longitud de onda del fluido, pero un
    # sólido blando (c_s chico) tiene λ_s = c_s/f corta y puede quedar sin resolver. Medido con grasa
    # (c_s = 100): λ_s = 1.0 mm contra un elemento de 0.53 mm, o sea 1.9 elementos por λ, y el error
    # sin restricción se queda en 353 % aunque el contraste sea casi 1. No es el mismo problema que
    # el de la restricción y hay que decirlo aparte.
    elemento = radio / modo["h_rel"]
    lambda_s = solido.c_s / f
    elementos_por_corte = lambda_s / elemento

    paso(0.9, "campo en el plano de propagación")
    mapa = _mapa_plano(mall, res["p_s"], modo["b_rel"] * radio, radio)

    # `sello_conjunto` y no `sello`: la UI espera la forma con la lista de criterios, igual que el
    # resto de las vistas (una sola forma para todos los sellos, D39).
    # Si el campo dispersado es tan débil que el error relativo no dice nada, se pasa None: el
    # registro lo convierte en `sin_referencia`, que es el estado honesto (no un aprobado por omisión).
    medido = error if np.abs(exacto).max() >= AMPLITUD_MINIMA else None
    sello = validacion.sello_conjunto({"volumen_esfera_elastica": medido})

    return {
        "ok": True,
        "metrica": {
            "error_relativo": error, "error_l2_relativo": error_l2,
            "error_absoluto": error_abs, "razon_amplitudes": razon,
            "ka": float(coef["ka"]), "contraste_impedancia": contraste,
            "restriccion_automatica": automatico,
            "elementos_por_lambda_corte": float(elementos_por_corte),
            "lambda_corte_mm": lambda_s * 1e3, "elemento_mm": elemento * 1e3,
            "radio_mm": radio * 1e3, "frecuencia_hz": f, "material": nombre_material,
            "n_puntos": int(util.sum()), "nodos": mall.n_nodes,
            "elementos_solidos": int(len(mall.tets_solid)),
            "elementos_fluidos": int(len(mall.tets_fluid)),
            "segundos": time.time() - t0, "modo": nombre_modo,
            "fijar_rigidos": fijar, "formulacion": nombre_formulacion,
            "residuo_rigido": None if res["residuo_rigido"] is None
            else float(np.abs(res["residuo_rigido"]).max()),
        },
        "angulos": {
            "theta_deg": np.degrees(th).tolist(),
            "fem_modulo": np.abs(fem).tolist(),
            "faran_modulo": np.abs(exacto).tolist(),
            "fem_fase_deg": np.degrees(np.angle(fem)).tolist(),
            "faran_fase_deg": np.degrees(np.angle(exacto)).tolist(),
        },
        "mapa": mapa,
        "sello": sello,
    }


def _mapa_plano(mall, p_s, dominio: float, radio: float, n: int = 81) -> dict[str, Any]:
    """``|p_s|`` en el plano x-z, interpolado punto por punto sobre los tetraedros del fluido.

    Se grafica el **dispersado** y no el total: en el total manda la incidente (que vale ~1) y la
    perturbación de pocos por ciento es invisible en una escala de dB. Es la misma lección de D42.

    No se usa `griddata`: interpolar en 3D sobre un plano es degenerado para Qhull (los puntos de
    consulta son coplanares y el símplex no tiene volumen). Se localiza el tetraedro contenedor, que
    además devuelve `nan` fuera del fluido — y ahí dentro está la esfera, que se enmascara en vez de
    rellenarse.
    """
    eje = np.linspace(-dominio, dominio, n)
    xx, zz = np.meshgrid(eje, eje)
    consulta = np.stack([xx.ravel(), np.zeros(xx.size), zz.ravel()], axis=1)
    z = coupled3d._interpolar(mall, p_s, consulta).reshape(n, n)
    z = np.where(np.isnan(z), np.nan, np.abs(z))
    return {
        "x_mm": (eje * 1e3).tolist(), "z_mm": (eje * 1e3).tolist(),
        "modulo": np.where(np.isnan(z), None, z).tolist(),
        "radio_mm": radio * 1e3, "dominio_mm": dominio * 1e3,
        "maximo": float(np.nanmax(z)) if np.any(~np.isnan(z)) else None,
    }
