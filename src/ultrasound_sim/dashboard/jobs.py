"""Cola de trabajos en proceso para las corridas largas de la UI (Fase U5).

Por qué existe
--------------
Un barrido de sensibilidad son ~25 solves de ~2 s, y un 3D son minutos. Un `POST` síncrono deja al
navegador colgado sin señal de vida, y si el usuario cierra la pestaña el trabajo sigue ocupando la
máquina sin que nadie lo sepa. Acá cada corrida larga es un **job** con estado, progreso, mensaje y
cancelación.

Decisiones
----------
- **Un solo worker.** La máquina es compartida y cada solve ya usa varios núcleos por BLAS: correr
  dos barridos en paralelo no acelera, sólo enturbia los tiempos. Los jobs se encolan.
- **Cancelación cooperativa.** El worker consulta ``cancelado()`` entre pasos (no entre operaciones
  de álgebra lineal): cancelar un solve de 2 s tarda hasta 2 s. Interrumpir dentro del solve no es
  posible sin threads dentro de SuperLU, y no vale la pena.
- **El resultado se guarda tal cual** en memoria del proceso, con tope de jobs vivos para que un
  script mal escrito no se coma la RAM.
"""

from __future__ import annotations

import threading
import time
import uuid
from collections import OrderedDict
from typing import Any, Callable

#: Tope de jobs retenidos (los viejos se descartan primero). 20 resultados de barrido son ~2 MB.
MAX_JOBS = 20

_ESTADOS = ("en_cola", "corriendo", "listo", "error", "cancelado")


class Trabajo:
    """Un job: estado, progreso y resultado. Todos los accesos van bajo ``_lock``."""

    def __init__(self, tipo: str, etiqueta: str) -> None:
        self.id = uuid.uuid4().hex[:12]
        self.tipo = tipo
        self.etiqueta = etiqueta
        self.estado = "en_cola"
        self.progreso = 0.0
        self.mensaje = "en cola"
        self.resultado: Any = None
        self.error: str | None = None
        self.t_creado = time.time()
        self.t_inicio: float | None = None
        self.t_fin: float | None = None
        self._cancelar = threading.Event()

    # ---- API para el worker ----
    def cancelado(self) -> bool:
        return self._cancelar.is_set()

    def avanzar(self, fraccion: float, mensaje: str) -> None:
        self.progreso = min(max(float(fraccion), 0.0), 1.0)
        self.mensaje = mensaje

    def pedir_cancelacion(self) -> None:
        self._cancelar.set()

    # ---- API para la UI ----
    def a_dict(self) -> dict[str, Any]:
        segundos = None
        if self.t_inicio is not None:
            fin = self.t_fin if self.t_fin is not None else time.time()
            segundos = round(fin - self.t_inicio, 2)
        return {
            "id": self.id,
            "tipo": self.tipo,
            "etiqueta": self.etiqueta,
            "estado": self.estado,
            "progreso": round(self.progreso, 4),
            "mensaje": self.mensaje,
            "error": self.error,
            "segundos": segundos,
            "resultado": self.resultado if self.estado == "listo" else None,
        }


class Cola:
    """Cola FIFO con un worker. ``enviar(fn, ...)`` devuelve el id del job."""

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._jobs: OrderedDict[str, Trabajo] = OrderedDict()
        self._pendientes: list[tuple[Trabajo, Callable[..., Any], tuple, dict]] = []
        self._despertar = threading.Condition(self._lock)
        self._worker: threading.Thread | None = None

    # ---- público ----
    def enviar(self, tipo: str, etiqueta: str, fn: Callable[..., Any],
               *args: Any, **kwargs: Any) -> Trabajo:
        trabajo = Trabajo(tipo, etiqueta)
        with self._despertar:
            self._jobs[trabajo.id] = trabajo
            self._podar()
            self._pendientes.append((trabajo, fn, args, kwargs))
            if self._worker is None or not self._worker.is_alive():
                self._worker = threading.Thread(target=self._loop, daemon=True,
                                                name="dashboard-jobs")
                self._worker.start()
            self._despertar.notify()
        return trabajo

    def estado(self, job_id: str) -> dict[str, Any] | None:
        with self._lock:
            trabajo = self._jobs.get(job_id)
            return trabajo.a_dict() if trabajo else None

    def cancelar(self, job_id: str) -> dict[str, Any] | None:
        with self._lock:
            trabajo = self._jobs.get(job_id)
            if trabajo is None:
                return None
            if trabajo.estado == "en_cola":
                # todavía no arrancó: se saca de la fila y se marca
                self._pendientes = [p for p in self._pendientes if p[0] is not trabajo]
                trabajo.estado = "cancelado"
                trabajo.mensaje = "cancelado antes de arrancar"
                trabajo.t_fin = time.time()
            elif trabajo.estado == "corriendo":
                trabajo.pedir_cancelacion()
                trabajo.mensaje = "cancelando… (termina el paso en curso)"
            return trabajo.a_dict()

    def _podar(self) -> None:
        """Descarta los jobs más viejos que ya terminaron, hasta volver al tope.

        Se llama al enviar **y al terminar**: si sólo se podara al enviar, una ráfaga de jobs que
        todavía están en cola no deja podar nada y la retención queda por encima del tope para
        siempre (lo detectó el test de retención). Los que corren o esperan no se tocan.
        """
        for viejo_id in list(self._jobs):
            if len(self._jobs) <= MAX_JOBS:
                break
            if self._jobs[viejo_id].estado in ("corriendo", "en_cola"):
                continue
            self._jobs.pop(viejo_id, None)

    def listar(self) -> list[dict[str, Any]]:
        with self._lock:
            return [t.a_dict() for t in self._jobs.values()]

    # ---- worker ----
    def _loop(self) -> None:
        while True:
            with self._despertar:
                while not self._pendientes:
                    self._despertar.wait(timeout=30.0)
                trabajo, fn, args, kwargs = self._pendientes.pop(0)
                if trabajo.estado == "cancelado":
                    continue
                trabajo.estado = "corriendo"
                trabajo.t_inicio = time.time()
                trabajo.mensaje = "arrancando"
            try:
                resultado = fn(trabajo, *args, **kwargs)
                with self._lock:
                    if trabajo.cancelado():
                        trabajo.estado = "cancelado"
                        trabajo.mensaje = "cancelado"
                    else:
                        trabajo.resultado = resultado
                        trabajo.estado = "listo"
                        trabajo.progreso = 1.0
                        trabajo.mensaje = "listo"
            except Exception as exc:  # el worker no puede morir por un job
                with self._lock:
                    if trabajo.cancelado():
                        # el worker levanta al ver la bandera: eso es una cancelación, no un fallo
                        trabajo.estado = "cancelado"
                        trabajo.mensaje = "cancelado"
                    else:
                        trabajo.estado = "error"
                        trabajo.error = f"{type(exc).__name__}: {exc}"
                        trabajo.mensaje = "error"
            finally:
                with self._lock:
                    trabajo.t_fin = time.time()
                    self._podar()


COLA = Cola()
