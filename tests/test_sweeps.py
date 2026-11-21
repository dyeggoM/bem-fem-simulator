"""Tests del barrido de sensibilidad (U5).

El criterio de aceptación de U5 es duro y verificable: el barrido de la UI tiene que reproducir el
ranking de D30 —los pesos ya publicados— porque si no coincide, el que está mal es el barrido de la
UI, no el solver (que tiene su propia validación en D29/D30).

Los tests con ``slow`` corren los 31 puntos (~90 s). El rápido corre una versión reducida para
fijar la forma y la lógica de cancelación.
"""

from __future__ import annotations

import json

import pytest

from ultrasound_sim.dashboard import sweeps

#: Pesos publicados en D30 (con membrana), en % del `t_interior` base.
D30 = {"espesor": 53.11, "c_p": 16.06, "atenuacion": 9.60, "c_s": 8.22, "rho": 1.03, "Z_r": 0.28}


class TrabajoFalso:
    """Lo mínimo que el barrido usa del job real."""

    def __init__(self) -> None:
        self.progreso = 0.0
        self.mensajes: list[str] = []
        self._cancelar_desde: int | None = None

    def avanzar(self, fraccion: float, mensaje: str) -> None:
        self.progreso = fraccion
        self.mensajes.append(mensaje)

    def cancelado(self) -> bool:
        return self._cancelar_desde is not None and len(self.mensajes) >= self._cancelar_desde


@pytest.fixture()
def puntos_reducidos(monkeypatch):
    """La lógica del barrido con un solver **sintético**: nada de resolver FEM.

    Los tests de estructura, cancelación y cálculo del rango no tienen por qué pagar 30 solves; lo
    que verifican es la aritmética del ranking y el manejo del job. El único test que ejercita el
    solver de verdad es el `slow` de aceptación.
    """
    monkeypatch.setattr(sweeps, "PUNTOS", {
        "espesor": [1.0e-3, 7.0e-3],
        "c_s": [1360.0, 1900.0],
        "Z_r": [0.5, 2.0],
    })

    def falso(*, c_p, c_s, rho, espesor_m, alpha, z_factor):
        # dependencia suave y monótona con cada cosa, para que los rangos sean predecibles
        return 0.5 + 0.2 * (espesor_m / 0.007) - 0.05 * (c_s / 1900.0) + 0.01 * z_factor

    monkeypatch.setattr(sweeps, "_t_interior", falso)


def test_las_leyes_de_atenuacion_salen_de_attenuation_py():
    """No se escriben α a mano: a 100 kHz son 3.777 / 6.243 / 15.312 Np/m, no números redondos."""
    assert sweeps.alpha_de_ley("sin_atenuacion") == 0.0
    assert sweeps.alpha_de_ley("white_2006") == pytest.approx(3.777, abs=0.01)
    assert sweeps.alpha_de_ley("pichardo_2011") == pytest.approx(6.243, abs=0.01)
    assert sweeps.alpha_de_ley("pinton_2012") == pytest.approx(15.312, abs=0.01)


def test_total_de_puntos_cuenta_el_base():
    assert sweeps.total_de_puntos() == sum(len(v) for v in sweeps.PUNTOS.values()) + 1


def test_barrido_reducido_arma_ranking_y_comparacion(puntos_reducidos):
    trabajo = TrabajoFalso()
    r = sweeps.barrido_sensibilidad(trabajo)

    assert r["ok"] is True
    assert r["puntos"] == 6
    assert trabajo.progreso == 1.0
    assert len(trabajo.mensajes) == 7   # 6 puntos + el mensaje del caso base
    assert {f["parametro"] for f in r["ranking"]} == {"espesor", "c_s", "Z_r"}
    # el ranking viene ordenado de mayor a menor rango
    rangos = [f["rango_pct"] for f in r["ranking"]]
    assert rangos == sorted(rangos, reverse=True)
    # la comparación trae el valor de D30 y su desvío
    assert all("d30_pct" in f and "desvio_pct" in f for f in r["comparacion"])
    json.dumps(r, allow_nan=False)


def test_el_barrido_se_puede_cancelar(puntos_reducidos):
    trabajo = TrabajoFalso()
    trabajo._cancelar_desde = 2   # a los dos mensajes, cancela
    with pytest.raises(RuntimeError) as exc:
        sweeps.barrido_sensibilidad(trabajo)
    assert "cancelado" in str(exc.value)


def test_el_rango_se_calcula_contra_el_base_y_no_contra_el_minimo(puntos_reducidos):
    """`rango_pct = (máx − mín) / t_interior(base)`: si se dividiera por el mínimo, todos los pesos
    quedarían inflados y el ranking cambiaría de escala sin cambiar de orden — el tipo de error que
    pasa desapercibido mirando sólo el orden."""
    trabajo = TrabajoFalso()
    r = sweeps.barrido_sensibilidad(trabajo)
    base = r["base"]["t_interior"]
    for fila in r["ranking"]:
        assert fila["rango_pct"] == pytest.approx(100.0 * (fila["max"] - fila["min"]) / base, rel=1e-9)


@pytest.mark.slow
def test_el_barrido_completo_reproduce_el_ranking_de_d30():
    """Aceptación de U5: los seis pesos, dentro del 5 % de lo publicado en D30."""
    trabajo = TrabajoFalso()
    r = sweeps.barrido_sensibilidad(trabajo)
    medidos = {f["parametro"]: f["rango_pct"] for f in r["ranking"]}
    print("\n  " + " · ".join(f"{k} {v:.2f} % (D30 {D30[k]:.2f})" for k, v in medidos.items()))
    for nombre, esperado in D30.items():
        assert medidos[nombre] == pytest.approx(esperado, rel=0.05), (
            f"{nombre}: {medidos[nombre]:.2f} % contra {esperado:.2f} % de D30"
        )
    # y el orden
    assert [f["parametro"] for f in r["ranking"]] == ["espesor", "c_p", "atenuacion", "c_s",
                                                      "rho", "Z_r"]
    assert r["base"]["t_interior"] == pytest.approx(0.53187, rel=0.01)
