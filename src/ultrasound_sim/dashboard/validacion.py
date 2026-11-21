"""Trazabilidad: el sello de cada resultado, calculado en el servidor (Fase U6).

El problema que resuelve
------------------------
Hasta U5 cada vista se fabricaba su propio sello en JavaScript con umbrales escritos a mano:

    const CRITERIO_AMPLITUD = 0.02;      // escena
    const REF_T_INTERIOR = 0.53187;      // cóclea
    const CRITERIO_GREEN = 10.0;         // 3D

Tres copias del mismo tipo de decisión, en el lugar donde nadie las audita y donde además nadie
puede saber de dónde salió el 0.53187. Si mañana cambia un criterio hay que buscar el número en el
JS, y si el criterio no aplica a la configuración corrida, la UI lo declara "validado" igual.

Acá cada criterio es **un dato con su fuente**: qué se compara, contra qué referencia, con qué
tolerancia y de qué decisión (ADR) o fase sale. El resultado de cada corrida trae su sello ya
evaluado, y la UI no decide nada — sólo pinta el estado que recibe.

Los tres estados posibles, y ninguno más:

- ``validado``: el medido entra en la tolerancia alrededor de la referencia.
- ``fuera``: no entra. El número se puede mostrar, pero no como validado.
- ``sin_referencia``: la configuración no tiene contra qué compararse (por ejemplo, parámetros
  fuera del caso base de D30). **No es un aprobado por omisión** — es el estado explícito de "esto
  no está verificado", que es lo que evita el malentendido más caro.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Any


@dataclass(frozen=True)
class Criterio:
    """Un criterio de aceptación con su procedencia."""

    id: str
    descripcion: str
    referencia: float
    tolerancia: float
    modo: str          # "abs" (|medido − ref| ≤ tol) o "rel" (|medido − ref|/|ref| ≤ tol)
    unidad: str
    fuente: str

    def evaluar(self, medido: float | None) -> dict[str, Any]:
        if medido is None:
            return {
                **asdict(self), "medido": None, "estado": "sin_referencia", "margen": None,
                # "no se midió" era engañoso: el número existe, lo que falta es contra qué compararlo
                "texto": (f"SIN REFERENCIA — {self.descripcion}: no hay contra qué comparar esta "
                          f"configuración. La referencia cubre sólo su caso base ({self.fuente})."),
            }

        desvio = abs(float(medido) - self.referencia)
        margen = desvio if self.modo == "abs" else (
            desvio / abs(self.referencia) if self.referencia else float("inf")
        )
        estado = "validado" if margen <= self.tolerancia else "fuera"

        if self.modo == "abs":
            comparacion = f"|{medido:.4g} − {self.referencia:.4g}| = {desvio:.4g} {self.unidad}"
            limite = f"≤ {self.tolerancia:.4g} {self.unidad}"
        else:
            comparacion = f"desvío {100.0 * margen:.2f} %"
            limite = f"≤ {100.0 * self.tolerancia:.2f} %"

        veredicto = "VALIDADO" if estado == "validado" else "FUERA DEL CRITERIO"
        texto = (
            f"{veredicto} — {self.descripcion}: {comparacion} ({limite}). Fuente: {self.fuente}."
        )
        return {
            **asdict(self),
            "medido": float(medido),
            "estado": estado,
            "margen": float(margen),
            "texto": texto,
        }


#: Registro de criterios. Cada entrada es un número que ya está medido y documentado en el repo: no
#: se agrega uno nuevo sin la corrida que lo respalde.
CRITERIOS: dict[str, Criterio] = {
    "escena_rt_r": Criterio(
        id="escena_rt_r",
        descripcion="error de R de la escena 2D contra la matriz de transferencia",
        referencia=0.0, tolerancia=0.02, modo="abs", unidad="amplitud",
        fuente="Fase 1 (R y T < 2 % en amplitud) / D34",
    ),
    "escena_rt_t": Criterio(
        id="escena_rt_t",
        descripcion="error de T de la escena 2D contra la matriz de transferencia",
        referencia=0.0, tolerancia=0.02, modo="abs", unidad="amplitud",
        fuente="Fase 1 (R y T < 2 % en amplitud) / D34",
    ),
    "dispersor_fem_serie": Criterio(
        id="dispersor_fem_serie",
        descripcion="FEM contra la serie de Bessel (cilindro rígido)",
        referencia=0.0, tolerancia=0.02, modo="abs", unidad="L2 relativo",
        fuente="Fase 4 (discrepancia < 2 %) / D35: la librería da 0.94 % a ka = 1",
    ),
    "dispersor_bem_serie": Criterio(
        id="dispersor_bem_serie",
        descripcion="BEM contra la serie de Bessel (cilindro rígido)",
        referencia=0.0, tolerancia=0.02, modo="abs", unidad="L2 relativo",
        fuente="Fase 4 / D35: la librería da 0.05 % a ka = 1",
    ),
    "cochlea_t_interior": Criterio(
        id="cochlea_t_interior",
        descripcion="t_interior del corte de cóclea contra el caso base de D30",
        referencia=0.53187, tolerancia=0.01, modo="rel", unidad="adimensional",
        fuente="D30 (espesor 2 mm, c_s 1900 m/s, α = 0, membrana del ballpark)",
    ),
    "cochlea_balance": Criterio(
        id="cochlea_balance",
        descripcion="cierre del balance de energía del corte (P_total + P_membrana)",
        referencia=0.0, tolerancia=0.01, modo="abs", unidad="de la potencia interceptada",
        fuente="D29 / D36: la identidad documentada es P_total + P_membrana = 0",
    ),
    "volumen_green": Criterio(
        id="volumen_green",
        descripcion="error contra el Green 3D de una fuente puntual",
        referencia=0.0, tolerancia=0.10, modo="abs", unidad="L2 relativo",
        fuente="D32 (8.4 % a ppw 8 en caja de 1.1 λ) / D38: la UI da 7.0 % a ppw 10",
    ),
    "volumen_rt_r": Criterio(
        id="volumen_rt_r",
        descripcion="error de R de la escena 3D contra la matriz de transferencia",
        referencia=0.0, tolerancia=0.02, modo="abs", unidad="amplitud",
        fuente="Fase 1 / D38",
    ),
    "volumen_rt_t": Criterio(
        id="volumen_rt_t",
        descripcion="error de T de la escena 3D contra la matriz de transferencia",
        referencia=0.0, tolerancia=0.02, modo="abs", unidad="amplitud",
        fuente="Fase 1 / D38",
    ),
    "barrido_d30": Criterio(
        id="barrido_d30",
        descripcion="peor desvío del ranking de sensibilidad contra los pesos de D30",
        referencia=0.0, tolerancia=0.05, modo="abs", unidad="% relativo",
        fuente="D30 (53.11 / 16.06 / 9.60 / 8.22 / 1.03 / 0.28 %) / D37",
    ),
}


def sello(criterio_id: str, medido: float | None) -> dict[str, Any]:
    """Evalúa un criterio contra el medido y devuelve el sello listo para la UI."""
    if criterio_id not in CRITERIOS:
        raise KeyError(f"criterio desconocido: {criterio_id!r}; hay {sorted(CRITERIOS)}")
    return CRITERIOS[criterio_id].evaluar(medido)


def sello_conjunto(pares: dict[str, float | None]) -> dict[str, Any]:
    """Sella varios criterios y resume el peor estado.

    El resumen **no** promedia: si uno falla, el conjunto falla. Un promedio de criterios dejaría
    pasar un error grande compensado por otro chico, que en validación es exactamente lo que no se
    quiere.
    """
    sellos = [sello(criterio_id, medido) for criterio_id, medido in pares.items()]
    estados = {s["estado"] for s in sellos}
    if "fuera" in estados:
        resumen = "fuera"
    elif "sin_referencia" in estados:
        resumen = "sin_referencia"
    else:
        resumen = "validado"
    return {
        "estado": resumen,
        "criterios": sellos,
        "texto": " ".join(s["texto"] for s in sellos),
    }


def catalogo() -> list[dict[str, Any]]:
    """Todos los criterios con su procedencia: es lo que la UI puede listar sin inventar nada."""
    return [asdict(c) for c in CRITERIOS.values()]
