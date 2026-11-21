"""Materiales de la cadena piel → hueso → cóclea, con su fuente (U9).

Cada valor sale de `docs/materiales-cadena.md`, que es la respuesta de Sim-Datos con cita por valor. **No
se agrega ningún número acá que no esté en ese documento**, y lo que no tiene fuente citable va como
``NC`` y se usa como valor declarado con ``verificado=False``, no como medición.

Tres advertencias que la UI tiene que arrastrar:

1. No existe un conjunto experimental completo validado de 0,1 a 1 MHz para todas las capas: la mayoría
   son **parámetros de modelado publicados**, no mediciones de cada tejido.
2. **Piel, grasa, músculo y los líquidos se modelan como fluidos** (sólo `c_L`): Sim-Datos no da `c_S`
   citable para tejido blando en esta banda, y un `c_S` inventado cambiaría la física sin avisar. El hueso
   sí entra como sólido elástico porque hay `c_S` de la literatura.
3. Donde `α` es ``NC`` se usa 0 y queda dicho: no se reemplaza por la atenuación del agua.
"""

from __future__ import annotations

from dataclasses import dataclass

from .materials import ElasticSolid, Material

__all__ = ["MedioCadena", "MEDIOS_CADENA", "medios_de", "resumen_medios"]

_DB_A_NP = 1.0 / 8.685889638065035


def _db_cm_a_np_m(valor_db_cm: float) -> float:
    """Convierte dB/cm a Np/m: ``1 dB/cm = 100 dB/m`` y ``1 Np/m = 8.6859 dB/m``.

    Escrito como función y no como constante en línea porque la primera versión lo dejó en
    ``valor*10/8.686`` (un factor 10 de menos) y el número seguía pareciendo razonable.
    """
    return valor_db_cm * 100.0 * _DB_A_NP


@dataclass(frozen=True)
class MedioCadena:
    nombre: str
    etiqueta: str
    tipo: str                     # 'fluido' | 'solido'
    material: Material | ElasticSolid
    fuente: str
    verificado: bool
    nota: str = ""


def _fluido(nombre, etiqueta, rho, c, alpha0=None, power=1.0, fuente="", verificado=False, nota=""):
    return MedioCadena(
        nombre=nombre, etiqueta=etiqueta, tipo="fluido",
        material=Material(name=nombre, rho=rho, c=c, alpha0=alpha0, alpha_power=power,
                          f_ref=1e6, provenance=fuente),
        fuente=fuente, verificado=verificado, nota=nota,
    )


def _solido(nombre, etiqueta, rho, c_p, c_s, fuente, verificado=False, nota=""):
    return MedioCadena(
        nombre=nombre, etiqueta=etiqueta, tipo="solido",
        material=ElasticSolid(name=nombre, rho=rho, c_p=c_p, c_s=c_s, provenance=fuente),
        fuente=fuente, verificado=verificado, nota=nota,
    )


MEDIOS_CADENA: dict[str, MedioCadena] = {
    "agua": _fluido(
        "agua", "agua", 1000.0, 1500.0, alpha0=None,
        fuente="valor de referencia del proyecto (sin pérdidas)", verificado=True,
        nota="se usa como medio de acople a los dos lados de la cadena",
    ),
    "piel": _fluido(
        "piel", "piel", 1116.0, 1537.0, alpha0=4.6, power=1.0,
        fuente="α = 4,6·f Np/m (Pichardo 2023, BabelBrain tabla 1, vía Duck)",
        verificado=False,
        nota="parámetro de MODELADO, sin validación experimental en 0,1–1 MHz; sin c_S citable",
    ),
    "grasa": _fluido(
        "grasa", "grasa", 950.0, 1478.0, alpha0=_db_cm_a_np_m(0.48),
        power=1.1,
        fuente="α = 0,48·f^1,1 dB/cm (Gu y Jing 2019, tabla I)",
        verificado=False,
        nota="grasa genérica, no subcutánea craneal; parámetros de simulación",
    ),
    "musculo": _fluido(
        "musculo", "músculo", 1050.0, 1547.0, alpha0=_db_cm_a_np_m(1.09),
        power=1.1,
        fuente="α = 1,09·f^1,1 dB/cm (Gu y Jing 2019, tabla I)",
        verificado=False,
        nota="sin orientación de fibras; no asumir equivalencia paralela/perpendicular",
    ),
    "hueso_cortical": _solido(
        "hueso_cortical", "hueso cortical", 1896.5, 2471.0, 1600.0,
        fuente="ρ 1896,5 (Pichardo 2023) · c_L 2471 ± 90 m/s a 0,836 MHz (Pichardo, Sin y Hynynen 2011) "
               "· c_S 1545–1758 m/s (Pichardo 2017)",
        verificado=False,
        nota="c_S se toma en el centro del rango citado (1600 m/s): es una elección, no un valor medido",
    ),
    "diploe": _solido(
        "diploe", "diploe", 1738.0, 2300.0, 1400.0,
        fuente="ρ 1738 (Pichardo 2023) · c_L 2300 ± 100 m/s a 0,836 MHz · c_S 1252–1574 m/s",
        verificado=False,
        nota="c_S en el centro del rango citado (1400 m/s): elección, no medición",
    ),
    "perilinfa": _fluido(
        "perilinfa", "perilinfa", 1000.0, 1500.0, alpha0=None,
        fuente="aproximación publicada a agua (Obrist 2019 §1.2); α = NC en la banda",
        verificado=False,
        nota="α no citable: se usa 0 y NO se reemplaza por la del agua",
    ),
    "endolinfa": _fluido(
        "endolinfa", "endolinfa", 1000.0, 1500.0, alpha0=None,
        fuente="parámetro de modelo (Modelling of magnetic vestibular stimulation, 2025, tabla 1)",
        verificado=False,
        nota="α no citable: se usa 0",
    ),
}


def medios_de(nombres) -> tuple[dict[str, Material], dict[str, ElasticSolid]]:
    """Separa los nombres pedidos en ``(fluidos, solidos)`` para `fem.multicapa`."""
    fluidos: dict[str, Material] = {}
    solidos: dict[str, ElasticSolid] = {}
    for nombre in nombres:
        if nombre not in MEDIOS_CADENA:
            raise ValueError(
                f"medio {nombre!r} no está en la cadena: disponibles {sorted(MEDIOS_CADENA)}"
            )
        medio = MEDIOS_CADENA[nombre]
        if isinstance(medio.material, ElasticSolid):
            solidos[nombre] = medio.material
        else:
            fluidos[nombre] = medio.material
    return fluidos, solidos


def t_capas_analitico(capas: list[tuple[str, float]], frecuencia: float,
                      medio_exterior: str = "agua") -> float | None:
    """``|t|`` exacto para capas **fluidas** por matriz de transferencia (incidencia normal).

    Devuelve ``None`` si hay alguna capa sólida: con cizalla la fórmula de capas fluidas no vale y no se
    inventa una referencia.

    Cada capa aporta ``M_j = [[cos φ, i·Z·sin φ], [i·sin φ/Z, cos φ]]`` con ``φ = k·d`` y ``Z = ρ·c``,
    que relaciona ``[p; u]`` entre sus dos caras. Con ``[p0; u0] = [1+R; (1-R)/Z1]`` y ``p_N/u_N = Z_N``:

        R = (β - α)/(α + β),  α = M11 - Z_N·M21,  β = (Z_N·M22 - M12)/Z1

    y la transmitida es ``p_N = M11·p0 + M12·u0`` — **no** ``1+R``, que es la presión en la primera
    interfaz y en una capa de cuarto de onda pasa de 1 (el transformador de impedancia).
    """
    import numpy as np

    z1 = None
    zinc = np.eye(2, dtype=complex)
    for nombre, espesor in capas:
        medio = MEDIOS_CADENA[nombre]
        mat = medio.material
        if isinstance(mat, ElasticSolid):
            return None                 # con cizalla la fórmula de capas fluidas no vale
        z = mat.rho * mat.c
        if z1 is None:
            z1 = z
        k = complex(mat.wavenumber(frecuencia))
        # con atenuación k es complejo: la matriz sigue valiendo, cos y sin de argumento complejo
        phi = k * espesor
        capa = np.array([[np.cos(phi), 1j * z * np.sin(phi)],
                         [1j * np.sin(phi) / z, np.cos(phi)]], dtype=complex)
        zinc = zinc @ capa

    ext = MEDIOS_CADENA[medio_exterior].material
    z_ext = ext.rho * ext.c
    if z1 is None:
        return 1.0
    alfa = zinc[0, 0] - z_ext * zinc[1, 0]
    beta = (z_ext * zinc[1, 1] - zinc[0, 1]) / z1
    reflejado = (beta - alfa) / (alfa + beta)
    p0 = 1.0 + reflejado
    u0 = (1.0 - reflejado) / z1
    return float(abs(zinc[0, 0] * p0 + zinc[0, 1] * u0))


def resumen_medios() -> list[dict]:
    """Lo que la UI muestra para que el usuario sepa con qué está simulando."""
    salida = []
    for medio in MEDIOS_CADENA.values():
        mat = medio.material
        item = {
            "nombre": medio.nombre,
            "etiqueta": medio.etiqueta,
            "tipo": medio.tipo,
            "rho": mat.rho,
            "c": getattr(mat, "c", None) or getattr(mat, "c_p", None),
            "c_s": getattr(mat, "c_s", None),
            "fuente": medio.fuente,
            "verificado": medio.verificado,
            "nota": medio.nota,
        }
        salida.append(item)
    return salida
