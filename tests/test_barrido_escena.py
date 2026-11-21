"""Tests del barrido de frecuencia de la escena 2D.

Criterio físico verificable: una losa de espesor `t` en un medio transmite al máximo en las
**resonancias de media onda**, `f_n = n·c_losa/(2t)`. Con hueso como fluido (c = 3500 m/s) y
t = 10 mm: 175, 350 kHz. Si el barrido no las encuentra —o las encuentra corridas— hay un error de
espesor, de velocidad o de fase, y ninguna de las tres se ve mirando un solo punto de frecuencia.

El barrido usa el solver de la escena (10–80 ms por punto), así que 41 puntos son un par de segundos:
los tests no necesitan marca `slow`.
"""

from __future__ import annotations

import numpy as np
import pytest

from ultrasound_sim.dashboard import sweeps

LOSA = {
    "ppw": 12,
    "alto_m": 20.0e-3,
    "capas": [
        {"material": "agua", "espesor_m": 20.0e-3},
        {"material": "hueso", "espesor_m": 10.0e-3},
        {"material": "agua", "espesor_m": 20.0e-3},
    ],
}


class TrabajoFalso:
    def __init__(self, cancelar_desde: int | None = None) -> None:
        self.progreso = 0.0
        self.mensajes: list[str] = []
        self._cancelar_desde = cancelar_desde

    def avanzar(self, fraccion: float, mensaje: str) -> None:
        self.progreso = fraccion
        self.mensajes.append(mensaje)

    def cancelado(self) -> bool:
        return self._cancelar_desde is not None and len(self.mensajes) >= self._cancelar_desde


# --------------------------------------------------------------------------- validación


def test_banda_invalida():
    with pytest.raises(ValueError):
        sweeps.barrido_escena(TrabajoFalso(), escena=LOSA, desde_hz=400e3, hasta_hz=100e3)


def test_puntos_fuera_de_rango():
    with pytest.raises(ValueError):
        sweeps.barrido_escena(TrabajoFalso(), escena=LOSA, puntos=2)


def test_escena_invalida():
    with pytest.raises(ValueError):
        sweeps.barrido_escena(TrabajoFalso(), escena="no soy una escena")


def test_se_puede_cancelar():
    trabajo = TrabajoFalso(cancelar_desde=3)
    with pytest.raises(RuntimeError) as exc:
        sweeps.barrido_escena(trabajo, escena=LOSA, puntos=21)
    assert "cancelado" in str(exc.value)


# --------------------------------------------------------------------------- física


def test_las_resonancias_de_media_onda_caen_donde_corresponde():
    """Los picos de T tienen que estar dentro del 3 % de n·c/(2t) = 175 y 350 kHz."""
    r = sweeps.barrido_escena(TrabajoFalso(), escena=LOSA, desde_hz=100e3, hasta_hz=400e3, puntos=41)
    predichas = r["prediccion"]["resonancias_hz"]
    encontradas = [p["f_hz"] for p in r["picos"]]
    print(f"\n  predichas {[round(f/1e3, 1) for f in predichas]} kHz · "
          f"encontradas {[round(f/1e3, 1) for f in encontradas]} kHz")

    assert len(encontradas) >= 2, f"no encontró las dos primeras resonancias: {encontradas}"
    for f_pred in predichas[:2]:
        cercana = min(encontradas, key=lambda f: abs(f - f_pred))
        assert abs(cercana - f_pred) / f_pred < 0.03, (
            f"la resonancia de {f_pred/1e3:.0f} kHz aparece en {cercana/1e3:.1f} kHz"
        )


def test_en_resonancia_T_sube_y_R_baja():
    """Física de la losa: en la resonancia de media onda la transmisión es máxima y la reflexión
    mínima. Si T y R subieran juntas, el barrido estaría mezclando magnitudes."""
    r = sweeps.barrido_escena(TrabajoFalso(), escena=LOSA, desde_hz=150e3, hasta_hz=200e3, puntos=21)
    filas = [f for f in r["filas"] if f["T"] is not None]
    i_T = max(range(len(filas)), key=lambda i: filas[i]["T"])
    i_R = min(range(len(filas)), key=lambda i: filas[i]["R"])
    print(f"\n  T máxima en {filas[i_T]['f_hz']/1e3:.1f} kHz ({filas[i_T]['T']:.3f}) · "
          f"R mínima en {filas[i_R]['f_hz']/1e3:.1f} kHz ({filas[i_R]['R']:.3f})")
    assert abs(filas[i_T]["f_hz"] - filas[i_R]["f_hz"]) <= 5e3
    assert filas[i_T]["T"] > 0.9
    assert filas[i_R]["R"] < 0.15


def test_el_barrido_trae_su_sello_y_los_datos_por_punto():
    r = sweeps.barrido_escena(TrabajoFalso(), escena=LOSA, puntos=9)
    assert r["puntos"] == 9
    assert len(r["filas"]) == 9
    assert all("f_hz" in fila for fila in r["filas"])
    assert r["sello"]["estado"] in ("validado", "fuera", "sin_referencia")
    assert r["sello"]["criterios"][0]["id"] == "escena_rt_t"
    assert r["max_error_T"] is not None


def test_el_error_crece_al_subir_la_frecuencia_con_la_misma_malla():
    """Con ppw fijo, a más frecuencia entran menos elementos por longitud de onda: el error de T
    tiene que ser peor arriba. Es lo que hace que el sello marque 'fuera' en 100–400 kHz a ppw 12, y
    por eso el barrido informa el peor error en vez de un promedio."""
    r = sweeps.barrido_escena(TrabajoFalso(), escena=LOSA, desde_hz=100e3, hasta_hz=400e3, puntos=31)
    filas = [f for f in r["filas"] if f["T"] is not None and f["error_T_abs"] is not None]
    abajo = np.mean([f["error_T_abs"] for f in filas[:10]])
    arriba = np.mean([f["error_T_abs"] for f in filas[-10:]])
    print(f"\n  error medio abajo {100*abajo:.2f} % · arriba {100*arriba:.2f} %")
    assert arriba > abajo, "el error no crece con la frecuencia: revisar la malla"
