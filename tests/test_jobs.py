"""Tests de la cola de trabajos (U5).

La cola es la pieza que evita que un barrido de 90 s deje el navegador colgado. Lo que hay que
fijar acá no es "corre rápido" sino: que el estado avance, que el error de un job no mate al worker,
que cancelar funcione en los dos momentos (en cola y corriendo) y que la retención no crezca sin
límite.
"""

from __future__ import annotations

import threading
import time

import pytest

from ultrasound_sim.dashboard.jobs import MAX_JOBS, Cola


def _esperar(cola: Cola, job_id: str, estados: tuple[str, ...], timeout: float = 20.0):
    limite = time.time() + timeout
    while time.time() < limite:
        estado = cola.estado(job_id)
        assert estado is not None
        if estado["estado"] in estados:
            return estado
        time.sleep(0.02)
    raise AssertionError(f"el job no llegó a {estados} en {timeout} s: {cola.estado(job_id)}")


def test_un_job_corre_y_devuelve_el_resultado():
    cola = Cola()
    trabajo = cola.enviar("prueba", "suma", lambda t, a, b: a + b, 2, 3)
    assert trabajo.estado == "en_cola"
    estado = _esperar(cola, trabajo.id, ("listo",))
    assert estado["resultado"] == 5
    assert estado["progreso"] == 1.0
    assert estado["segundos"] is not None and estado["segundos"] >= 0.0


def test_el_progreso_y_el_mensaje_se_actualizan():
    cola = Cola()

    def trabajo_lento(t):
        for i in range(4):
            t.avanzar((i + 1) / 4, f"paso {i + 1}")
            time.sleep(0.05)
        return "ok"

    trabajo = cola.enviar("prueba", "progreso", trabajo_lento)
    _esperar(cola, trabajo.id, ("running", "corriendo"), timeout=5.0)
    time.sleep(0.06)
    intermedio = cola.estado(trabajo.id)
    assert 0.0 < intermedio["progreso"] <= 1.0
    assert intermedio["mensaje"].startswith("paso")
    final = _esperar(cola, trabajo.id, ("listo",))
    assert final["resultado"] == "ok"


def test_un_error_no_mata_al_worker():
    cola = Cola()

    def falla(t):
        raise ValueError("esto se rompe a propósito")

    roto = cola.enviar("prueba", "falla", falla)
    estado = _esperar(cola, roto.id, ("error",))
    assert "ValueError" in estado["error"]
    assert "a propósito" in estado["error"]

    # el worker sigue vivo: el siguiente job corre
    bueno = cola.enviar("prueba", "sigue", lambda t: 42)
    assert _esperar(cola, bueno.id, ("listo",))["resultado"] == 42


def test_cancelar_en_cola_no_llega_a_correr():
    cola = Cola()
    barrera = threading.Event()

    def bloquea(t):
        barrera.wait(timeout=10)
        return "primero"

    primero = cola.enviar("prueba", "bloquea", bloquea)
    _esperar(cola, primero.id, ("corriendo",), timeout=5.0)

    segundo = cola.enviar("prueba", "cancelame", lambda t: "no debería correr")
    assert segundo.estado == "en_cola"
    cancelado = cola.cancelar(segundo.id)
    assert cancelado["estado"] == "cancelado"

    barrera.set()
    _esperar(cola, primero.id, ("listo",))
    time.sleep(0.2)
    assert cola.estado(segundo.id)["estado"] == "cancelado"
    assert cola.estado(segundo.id)["resultado"] is None


def test_cancelar_corriendo_corta_entre_pasos():
    cola = Cola()

    def largo(t):
        for _ in range(200):
            if t.cancelado():
                raise RuntimeError("cancelado por el usuario")
            time.sleep(0.01)
        return "terminó"

    trabajo = cola.enviar("prueba", "largo", largo)
    _esperar(cola, trabajo.id, ("corriendo",), timeout=5.0)
    time.sleep(0.05)
    cola.cancelar(trabajo.id)
    estado = _esperar(cola, trabajo.id, ("cancelado", "listo"), timeout=10.0)
    assert estado["estado"] == "cancelado"
    assert estado["resultado"] is None


def test_job_desconocido_devuelve_none():
    cola = Cola()
    assert cola.estado("noexiste") is None
    assert cola.cancelar("noexiste") is None


def test_la_retencion_no_crece_sin_limite():
    """Al pasarse del tope se descartan los más viejos **ya terminados**.

    Se encola en serie (esperando cada uno) a propósito: los jobs viejos se podan, así que
    preguntar por uno de los primeros devuelve `None` — y eso es lo correcto, no un bug. Lo que hay
    que garantizar es que la retención no crezca y que los últimos estén.
    """
    cola = Cola()
    for i in range(MAX_JOBS + 3):
        trabajo = cola.enviar("prueba", f"job {i}", lambda t, i=i: i)
        _esperar(cola, trabajo.id, ("listo",), timeout=15.0)
    vivos = cola.listar()
    assert len(vivos) <= MAX_JOBS
    assert vivos[-1]["etiqueta"] == f"job {MAX_JOBS + 2}"
    # el primero ya no está: se podó
    assert not any(j["etiqueta"] == "job 0" for j in vivos)


def test_listar_devuelve_todos_los_campos_que_usa_la_ui():
    cola = Cola()
    trabajo = cola.enviar("sensibilidad", "x", lambda t: {"a": 1})
    _esperar(cola, trabajo.id, ("listo",))
    estado = cola.estado(trabajo.id)
    for clave in ("id", "tipo", "etiqueta", "estado", "progreso", "mensaje", "error",
                  "segundos", "resultado"):
        assert clave in estado
    assert estado["resultado"] == {"a": 1}
