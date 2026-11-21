"""Tests del barrido de frecuencia y del Q (consulta a Sim-Datos del 2026-11-20).

La definición de Q **no es biológica**: Sim-Datos confirmó que no existe un Q citable de la resonancia
de la cápsula ótica (ni en 0.1–1 kHz ni a 100 kHz–1 MHz) y fijó la definición operativa:

    Q = f_pico/(f2 − f1)  con  |H|² = |H|²_max/2,  o sea  |H| = |H|_max/√2

Son los mismos cruces con 20·log10|H| o con 10·log10|H|²; **no** es "media amplitud". Lo que sí está
documentado, y no es esto: Q_ERB ≈ 12.7 ± 1 a 1 kHz (Shera, Guinan y Oxenham 2002) es selectividad de
los filtros cocleares, no la resonancia global de la cápsula.

Los tests del ancho de banda usan una función sintética para no pagar 13 solves; el test de
integración usa el solver real y está marcado `slow`.
"""

from __future__ import annotations

import numpy as np
import pytest

from ultrasound_sim.dashboard import sweeps


class TrabajoFalso:
    def __init__(self) -> None:
        self.progreso = 0.0
        self.mensajes: list[str] = []

    def avanzar(self, fraccion: float, mensaje: str) -> None:
        self.progreso = fraccion
        self.mensajes.append(mensaje)

    def cancelado(self) -> bool:
        return False


def _campana(f0: float, q: float, x: np.ndarray) -> np.ndarray:
    """Lorentziana con ancho a media potencia f0/q, para tener un Q conocido de antemano."""
    gamma = f0 / q
    return 1.0 / np.sqrt(1.0 + ((x - f0) / (gamma / 2.0)) ** 2)


# --------------------------------------------------------------------------- ancho de banda


@pytest.mark.parametrize("q_esperado", [1.0, 2.0, 5.0])
def test_el_q_recupera_el_ancho_de_media_potencia(q_esperado):
    """Con una campana de Q conocido, la medición tiene que devolver ese Q."""
    f0 = 100e3
    x = np.linspace(f0 * 0.5, f0 * 1.5, 401)
    y = _campana(f0, q_esperado, x)
    r = sweeps._campana_decaimiento(y, x)
    print(f"\n  Q esperado {q_esperado}: medido {r['Q']:.3f} · f_pico {r['f_pico_hz'] / 1e3:.2f} kHz")
    assert r["Q"] == pytest.approx(q_esperado, rel=0.05)


def test_una_respuesta_plana_da_cota_superior_y_no_un_q():
    """Si no baja 3 dB dentro de la banda, el ancho es al menos el de la banda: Q queda acotado.

    Es el caso real del corte con membrana: no hay pico estrecho que medir, y eso **es** el resultado
    (la membrana rompió la resonancia), no un fallo del barrido.
    """
    f0 = 100e3
    x = np.linspace(85e3, 115e3, 13)
    y = 1.0 + 0.02 * np.cos(2 * np.pi * (x - f0) / 20e3)   # ondulación suave, nunca cae 30 %
    r = sweeps._campana_decaimiento(y, x)
    assert r["Q"] is None
    assert r["Q_cota_superior"] == pytest.approx(100e3 / 30e3, rel=0.1)
    assert "ancho" in r["motivo"]


def test_un_pico_en_el_borde_no_da_q_medible():
    x = np.linspace(80e3, 120e3, 21)
    y = np.exp(-((x - 120e3) ** 2) / (2 * (3e3) ** 2))   # el máximo cae en el último punto
    r = sweeps._campana_decaimiento(y, x)
    assert r["Q"] is None
    assert "cota" in r["motivo"] or "ancho" in r["motivo"]


def test_curva_corta_no_inventa_nada():
    r = sweeps._campana_decaimiento(np.array([1.0, 2.0]), np.array([1.0, 2.0]))
    assert r["Q"] is None
    assert "incompleta" in r["motivo"]


def test_la_resolucion_se_declara():
    """Un Q mayor que 1/paso relativo no se resuelve: la respuesta lo tiene que decir."""
    r = sweeps.barrido_frecuencia(TrabajoFalso(), desde_hz=98e3, hasta_hz=102e3, puntos=5)
    assert "paso relativo" in r["resolucion"]
    assert r["puntos"] == 5


# --------------------------------------------------------------------------- validación del pedido


def test_banda_invalida():
    with pytest.raises(ValueError):
        sweeps.barrido_frecuencia(TrabajoFalso(), desde_hz=120e3, hasta_hz=100e3)


def test_puntos_fuera_de_rango():
    with pytest.raises(ValueError):
        sweeps.barrido_frecuencia(TrabajoFalso(), puntos=100)


# --------------------------------------------------------------------------- integración (slow)


@pytest.mark.slow
def test_el_corte_con_membrana_no_tiene_un_q_resoluble():
    """Integración real: la membrana ensancha tanto la respuesta que Q queda acotado.

    Coherente con D29 (`Q` baja de 287 a < 5 cuando entra la membrana): acá se mide con la definición
    operativa y sale una cota, no un valor. El número exacto depende de la banda, así que se exige la
    coherencia física (cota finita y menor que la del caso sin membrana), no una cifra.
    """
    r = sweeps.barrido_frecuencia(TrabajoFalso(), desde_hz=85e3, hasta_hz=115e3, puntos=13)
    print(f"\n  Q {r['q']['Q']} · cota {r['q'].get('Q_cota_superior')} · "
          f"pico {r['q']['f_pico_hz'] / 1e3:.1f} kHz")
    assert r["q"]["Q"] is None
    assert r["q"]["Q_cota_superior"] is not None
    assert r["q"]["Q_cota_superior"] < 10.0
    assert len(r["t_interior"]) == 13
    assert all(0.0 < v < 2.0 for v in r["t_interior"])
