"""Barrido en `ka` de la esfera elástica: dónde vale cada formulación (U7b).

Por qué existe
--------------
D50 midió que la formulación dispersada da la traslación `3/(2g+1)` al 1-2 % en el límite de Rayleigh
y, sin embargo, contra Faran da 19 % en acero a `ka = 1`. Un número suelto no dice **dónde** empieza a
fallar. Este barrido resuelve el mismo material en varios `ka` y grafica el error de cada formulación
en función de `ka`, que es la variable que gobierna la dispersión.

Lo que el barrido NO hace
-------------------------
No elige por vos: devuelve el error de las dos formulaciones en cada punto, con su sello. Si la
dispersada es peor en todo el rango, eso es un resultado, no un fallo del barrido, y la UI lo dice.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from . import elastic_sphere, validacion

__all__ = ["barrido_ka", "puntos_ka"]

#: Mismo rango que la vista: por debajo de 0.3 no hay dispersión que medir y por encima de 3 el
#: patrón ya tiene muchos lóbulos y el máximo angular se vuelve más frágil todavía.
_KA_MIN, _KA_MAX = 0.3, 3.0


def puntos_ka(desde: float, hasta: float, n: int) -> list[float]:
    desde, hasta = float(desde), float(hasta)
    if not (_KA_MIN <= desde <= _KA_MAX and _KA_MIN <= hasta <= _KA_MAX):
        raise ValueError(f"ka fuera de [{_KA_MIN}, {_KA_MAX}]: {desde}, {hasta}")
    if hasta <= desde:
        raise ValueError(f"hasta ({hasta}) tiene que ser mayor que desde ({desde})")
    n = int(n)
    if not 2 <= n <= 25:
        raise ValueError(f"puntos tiene que estar en [2, 25]: {n}")
    return [float(x) for x in np.linspace(desde, hasta, n)]


def barrido_ka(datos: dict, trabajo=None) -> dict[str, Any]:
    """Resuelve la esfera en cada `ka` para la formulación pedida y devuelve error contra Faran.

    ``trabajo`` es el job de la cola: se le informa el avance por punto, que es lo que mira la barra.
    """
    def paso(fraccion: float, mensaje: str) -> None:
        if trabajo is not None:
            trabajo.avanzar(fraccion, mensaje)

    valores = puntos_ka(
        datos.get("desde", 0.5), datos.get("hasta", 2.0), datos.get("puntos", 7),
    )
    formulacion = str(datos.get("formulacion", "total"))
    if formulacion not in elastic_sphere.FORMULACIONES:
        raise ValueError(f"formulacion desconocida {formulacion!r}")
    material = str(datos.get("material", "acero"))
    modo = str(datos.get("modo", "rapida"))

    # Comparar las dos formulaciones en la misma corrida: es lo que muestra el cruce (D51) sin tener
    # que lanzar dos barridos y acordarse de los números. Cuesta el doble, y la UI lo dice.
    comparar = bool(datos.get("comparar", False))
    formulaciones = ["total", "dispersada"] if comparar else [formulacion]

    total_pasos = len(valores) * len(formulaciones)
    hecho = 0
    por_formulacion: dict[str, list[dict[str, Any]]] = {}
    for form in formulaciones:
        filas: list[dict[str, Any]] = []
        for ka in valores:
            # cancelación cooperativa entre puntos (desde dentro de un solve no se puede interrumpir)
            if trabajo is not None and trabajo.cancelado():
                raise RuntimeError("cancelado por el usuario")
            paso(hecho / total_pasos, f"{form}: ka = {ka:.2f} ({hecho + 1}/{total_pasos})")
            r = elastic_sphere.resolver({
                "ka": ka, "material": material, "modo": modo, "formulacion": form,
                "frecuencia_hz": datos.get("frecuencia_hz", 100e3),
                "puntos_angulo": int(datos.get("puntos_angulo", 61)),
                "con_mapa": False,
            })
            m = r["metrica"]
            filas.append({
                "ka": float(m["ka"]),
                "error_relativo": float(m["error_relativo"]),
                "error_l2_relativo": float(m["error_l2_relativo"]),
                "razon_amplitudes": float(m["razon_amplitudes"]),
                "elementos_por_lambda_corte": float(m["elementos_por_lambda_corte"]),
                "sello": r["sello"]["estado"],
                "nodos": int(m["nodos"]),
                "segundos": float(m["segundos"]),
            })
            hecho += 1
        por_formulacion[form] = filas

    filas = por_formulacion[formulacion]

    paso(0.99, "listo")
    peor = max(filas, key=lambda f: f["error_relativo"])
    dentro = [f for f in filas if f["sello"] == "validado"]

    # cuando se comparan las dos, se manda la serie de cada una y además el cruce: el primer `ka`
    # donde deja de ganar la total. Ese punto es el resultado que interesa, no una curva más.
    series = {
        form: {
            "ka": [f["ka"] for f in fs],
            "error_relativo": [f["error_relativo"] for f in fs],
            "error_l2_relativo": [f["error_l2_relativo"] for f in fs],
            "sello": [f["sello"] for f in fs],
        }
        for form, fs in por_formulacion.items()
    }
    cruce = None
    if comparar:
        for i, ka in enumerate(valores):
            if por_formulacion["dispersada"][i]["error_relativo"] < por_formulacion["total"][i]["error_relativo"]:
                cruce = float(ka)
                break

    return {
        "ok": True,
        "barrido": {
            "material": material, "modo": modo, "formulacion": formulacion,
            "n_puntos": len(filas), "ka": [f["ka"] for f in filas],
            "error_relativo": [f["error_relativo"] for f in filas],
            "error_l2_relativo": [f["error_l2_relativo"] for f in filas],
            "razon_amplitudes": [f["razon_amplitudes"] for f in filas],
            "sello": [f["sello"] for f in filas],
            "comparar": comparar, "series": series if comparar else None,
            "cruce_desde_ka": cruce,
        },
        "filas": filas,
        "resumen": {
            "peor_error": float(peor["error_relativo"]), "peor_ka": float(peor["ka"]),
            "mejor_error": float(min(f["error_relativo"] for f in filas)),
            "n_dentro": len(dentro),
            "dentro_de": [f["ka"] for f in dentro],
            "tolerancia": validacion.CRITERIOS["volumen_esfera_elastica"].tolerancia,
            "segundos_totales": float(sum(f["segundos"] for f in filas)),
        },
    }
