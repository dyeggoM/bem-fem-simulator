"""Cargar y exportar un arreglo de fuentes planeado (U10).

Formato
-------
Dos, los dos con las mismas columnas y unidades explícitas en la cabecera:

* **CSV** con cabecera ``y_mm,ancho_mm,amp,fase_rad`` (una fila por elemento).
* **JSON**: ``{"fuentes": [{"y_mm": …, "ancho_mm": …, "amp": …, "fase_rad": …}, …]}``.

Se aceptan también los nombres en inglés (``y,width,amp,phase``) porque el que exporta el sistema escribe
en inglés y no tiene sentido obligar a renombrar para volver a cargar.

Lo que se valida, y por qué
---------------------------
1. **Unidades**: las columnas se llaman ``*_mm`` a propósito. Se convierte a metros y se dice en el
   resultado, en vez de asumir que el número está en metros (32 filas en metros serían 32 km).
2. **Elementos que no se solapan**: dos elementos con el mismo ``y`` o pisados no son un arreglo físico.
   Se compara con ``|Δy| >= (w_i + w_j)/2``; si no, se dice cuáles.
3. **``amp`` > 0 y ``ancho_mm`` > 0**: un elemento de ancho cero no radia y una amplitud negativa es una
   fase de π disfrazada.
"""

from __future__ import annotations

import csv
import io
import json
from typing import Any

__all__ = ["ErrorArreglo", "cargar_arreglo", "exportar_arreglo", "COLUMNAS"]


class ErrorArreglo(ValueError):
    """El archivo no describe un arreglo de fuentes usable."""


COLUMNAS = ("y_mm", "ancho_mm", "amp", "fase_rad")

_ALIAS = {
    "y_mm": ("y_mm", "y", "posicion_mm", "posición_mm"),
    "ancho_mm": ("ancho_mm", "width_mm", "width", "ancho"),
    "amp": ("amp", "amplitud", "amplitude", "amp_rel"),
    "fase_rad": ("fase_rad", "phase", "phase_rad", "fase"),
}


def _clave(fila: dict, columna: str) -> Any:
    for nombre in _ALIAS[columna]:
        if nombre in fila:
            return fila[nombre]
    return None


def _a_float(valor, columna: str, indice: int) -> float:
    if valor is None or (isinstance(valor, str) and not valor.strip()):
        raise ErrorArreglo(f"fila {indice}: falta la columna '{columna}'")
    try:
        return float(valor)
    except (TypeError, ValueError) as exc:
        raise ErrorArreglo(f"fila {indice}: '{columna}' = {valor!r} no es un número") from exc


def _validar_solapes(fuentes: list[dict]) -> None:
    orden = sorted(range(len(fuentes)), key=lambda i: fuentes[i]["y"])
    for a, b in zip(orden, orden[1:]):
        fa, fb = fuentes[a], fuentes[b]
        separacion = fb["y"] - fa["y"]
        minimo = 0.5 * (fa["width"] + fb["width"])
        if separacion < minimo - 1e-12:
            raise ErrorArreglo(
                f"los elementos {a + 1} y {b + 1} se solapan: están a {separacion * 1e3:.3f} mm y "
                f"necesitan al menos {minimo * 1e3:.3f} mm ((ancho_1 + ancho_2)/2)"
            )


def cargar_arreglo(texto: str, formato: str | None = None) -> list[dict]:
    """Devuelve ``[{y, width, amp, phase}]`` en unidades SI, listo para la vista de campo."""
    texto = (texto or "").strip()
    if not texto:
        raise ErrorArreglo("el archivo está vacío")

    if formato is None:
        formato = "json" if texto.lstrip().startswith(("{", "[")) else "csv"

    if formato == "json":
        try:
            crudo = json.loads(texto)
        except json.JSONDecodeError as exc:
            raise ErrorArreglo(f"JSON inválido: {exc.msg} (línea {exc.lineno})") from exc
        if isinstance(crudo, dict):
            crudo = crudo.get("fuentes") or crudo.get("sources") or crudo.get("arreglo")
        if not isinstance(crudo, list):
            raise ErrorArreglo(
                "el JSON tiene que ser una lista de fuentes o un objeto con la clave 'fuentes'"
            )
        filas = crudo
    elif formato == "csv":
        lector = csv.DictReader(io.StringIO(texto))
        if not lector.fieldnames:
            raise ErrorArreglo("el CSV no tiene cabecera")
        normalizada = {c.strip().lower(): c for c in lector.fieldnames}
        tiene_y = any(a in normalizada for a in _ALIAS["y_mm"])
        filas = []
        for fila in lector:
            filas.append({k.strip().lower(): v for k, v in fila.items() if k})
        if not tiene_y:
            raise ErrorArreglo(
                f"el CSV no tiene columna de posición; cabecera leída: {sorted(normalizada)} "
                f"(se espera {', '.join(COLUMNAS)})"
            )
        if not filas:
            raise ErrorArreglo("el CSV tiene cabecera pero ninguna fila de fuentes")
    else:
        raise ErrorArreglo(f"formato {formato!r} no soportado: 'csv' o 'json'")

    if not isinstance(filas, list) or not filas:
        raise ErrorArreglo("el archivo no tiene ninguna fuente")

    fuentes: list[dict] = []
    for i, fila in enumerate(filas):
        if not isinstance(fila, dict):
            raise ErrorArreglo(f"fila {i + 1}: se esperaba un objeto")
        y_mm = _a_float(_clave(fila, "y_mm"), "y_mm", i + 1)
        ancho_mm = _a_float(_clave(fila, "ancho_mm"), "ancho_mm", i + 1)
        amp = _a_float(_clave(fila, "amp"), "amp", i + 1)
        fase = _a_float(_clave(fila, "fase_rad"), "fase_rad", i + 1)
        if not (ancho_mm > 0.0):
            raise ErrorArreglo(f"fila {i + 1}: ancho_mm = {ancho_mm} tiene que ser > 0")
        if not (amp > 0.0):
            raise ErrorArreglo(
                f"fila {i + 1}: amp = {amp} tiene que ser > 0 (una amplitud negativa es una fase de π)"
            )
        fuentes.append({"y": y_mm * 1e-3, "width": ancho_mm * 1e-3, "amp": amp, "phase": fase})

    _validar_solapes(fuentes)
    return fuentes


def exportar_arreglo(fuentes) -> str:
    """CSV con cabecera, en mm y radianes: lo mismo que `cargar_arreglo` acepta."""
    salida = io.StringIO()
    escritor = csv.writer(salida, lineterminator="\n")
    escritor.writerow(list(COLUMNAS))
    for f in fuentes:
        y = f["y"] if isinstance(f, dict) else getattr(f, "y")
        ancho = f["width"] if isinstance(f, dict) else getattr(f, "width")
        amp = f["amp"] if isinstance(f, dict) else getattr(f, "amp")
        fase = f["phase"] if isinstance(f, dict) else getattr(f, "phase")
        escritor.writerow([f"{y * 1e3:.4f}", f"{ancho * 1e3:.4f}", f"{amp:.6g}", f"{fase:.6g}"])
    return salida.getvalue()
