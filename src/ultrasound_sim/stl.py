"""Lector de STL con control de que la superficie sirva para mallar un volumen (U9).

Por qué un lector propio
------------------------
`gmsh` lee STL, pero no dice si la malla está rota: si el STL tiene un agujero, gmsh mallará *algo* y el
solve dará un número. En un simulador eso es peor que un error, porque el número parece válido.

Lo que se controla acá, antes de mallar:

1. **Formato**: binario o ASCII, detectado por tamaño y no por el texto (hay STL binarios que empiezan
   con la palabra ``solid`` y rompen los lectores que miran sólo la primera línea).
2. **Cerrada**: cada arista compartida por **exactamente dos** triángulos. Se reporta cuántas aristas
   tienen uno (agujero) y cuántas más de dos (superficie que se toca a sí misma).
3. **Normales hacia afuera**: el volumen por divergencia tiene que ser **positivo**. Un STL con las
   caras invertidas malla igual y da un sólido con masa negativa; el acople sale con el signo cambiado.

El volumen se calcula con la fórmula de la divergencia sobre las caras (suma de tetraedros desde el
origen), que no necesita que el sólido contenga al origen.
"""

from __future__ import annotations

import struct
from dataclasses import dataclass
from pathlib import Path

import numpy as np

__all__ = [
    "MallaSTL", "leer_stl", "leer_stl_bytes", "verificar_stl", "escerrada", "volumen_con_signo",
    "cubo_stl", "seccion_por_plano", "area_contorno",
]


class _Nombre:
    """Sólo para que los mensajes de error digan el nombre del archivo."""

    def __init__(self, nombre: str):
        self.name = nombre


@dataclass
class MallaSTL:
    """Triángulos ``(M, 3, 3)`` en metros y las normales declaradas ``(M, 3)``."""

    triangulos: np.ndarray
    normales: np.ndarray
    origen: str = ""

    @property
    def n_triangulos(self) -> int:
        return int(self.triangulos.shape[0])

    def bbox(self) -> tuple[np.ndarray, np.ndarray]:
        v = self.triangulos.reshape(-1, 3)
        return v.min(axis=0), v.max(axis=0)

    def diagonal(self) -> float:
        lo, hi = self.bbox()
        return float(np.linalg.norm(hi - lo))

    def resumen(self) -> dict:
        lo, hi = self.bbox()
        return {
            "n_triangulos": self.n_triangulos,
            "bbox_min_mm": np.round(lo * 1e3, 3).tolist(),
            "bbox_max_mm": np.round(hi * 1e3, 3).tolist(),
            "diagonal_mm": round(self.diagonal() * 1e3, 3),
            "volumen_mm3": round(volumen_con_signo(self) * 1e9, 3),
        }


def _leer_binario(datos: bytes) -> tuple[np.ndarray, np.ndarray]:
    n = struct.unpack_from("<I", datos, 80)[0]
    esperado = 84 + 50 * n
    if len(datos) != esperado:
        raise ValueError(
            f"STL binario inconsistente: dice {n} triángulos ({esperado} bytes) y el archivo tiene "
            f"{len(datos)} bytes"
        )
    crudo = np.frombuffer(datos, dtype=np.uint8, count=50 * n, offset=84)
    crudo = crudo.reshape(n, 50)
    flotantes = crudo[:, :48].copy().view("<f4").reshape(n, 12)
    normales = flotantes[:, 0:3].astype(float)
    triangulos = flotantes[:, 3:12].reshape(n, 3, 3).astype(float)
    return triangulos, normales


def _leer_ascii(texto: str) -> tuple[np.ndarray, np.ndarray]:
    normales, triangulos = [], []
    actual_n = np.zeros(3)
    vertices: list[list[float]] = []
    for linea in texto.splitlines():
        partes = linea.split()
        if not partes:
            continue
        clave = partes[0].lower()
        if clave == "facet" and len(partes) >= 5 and partes[1].lower() == "normal":
            actual_n = np.array([float(v) for v in partes[2:5]])
            vertices = []
        elif clave == "vertex" and len(partes) >= 4:
            vertices.append([float(v) for v in partes[1:4]])
        elif clave == "endfacet" and len(vertices) == 3:
            normales.append(actual_n)
            triangulos.append(vertices)
    if not triangulos:
        raise ValueError("el STL ASCII no tiene ninguna cara ('facet ... endfacet')")
    return np.asarray(triangulos, dtype=float), np.asarray(normales, dtype=float)


def leer_stl(ruta: str | Path) -> MallaSTL:
    """Lee un STL binario o ASCII desde un archivo. Las coordenadas se toman como **metros**."""
    ruta = Path(ruta)
    return leer_stl_bytes(ruta.read_bytes(), ruta.name)


def leer_stl_bytes(datos: bytes, nombre: str = "STL") -> MallaSTL:
    """Igual que `leer_stl` pero desde memoria (lo que llega por HTTP)."""
    ruta = _Nombre(nombre)
    if len(datos) < 84:
        raise ValueError(f"{ruta.name}: el archivo es demasiado chico para ser un STL ({len(datos)} bytes)")
    # el tamaño decide: un binario mide exactamente 84 + 50*n, y un ASCII casi nunca coincide
    n_bin = struct.unpack_from("<I", datos, 80)[0] if len(datos) >= 84 else 0
    if len(datos) == 84 + 50 * n_bin and n_bin > 0:
        triangulos, normales = _leer_binario(datos)
        origen = "binario"
    else:
        try:
            texto = datos.decode("utf-8", errors="strict")
        except UnicodeDecodeError as exc:
            raise ValueError(
                f"{ruta.name}: no es un STL binario consistente (el tamaño no cierra) ni texto UTF-8 válido"
            ) from exc
        triangulos, normales = _leer_ascii(texto)
        origen = "ascii"
    return MallaSTL(triangulos=triangulos, normales=normales, origen=origen)


def _aristas(triangulos: np.ndarray) -> np.ndarray:
    """Aristas ``(3M, 2)`` con los vértices ordenados, para poder contarlas sin importar la orientación."""
    t = np.asarray(triangulos, dtype=float)
    aristas = np.concatenate([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]], axis=0)
    aristas = np.sort(aristas, axis=1)
    # clave por vértice redondeado: dos caras que comparten una arista escriben el mismo float
    return np.ascontiguousarray(np.round(aristas, 12))


def escerrada(triangulos: np.ndarray) -> tuple[bool, dict]:
    """``True`` si cada arista está en exactamente dos triángulos, con el detalle de lo que falla."""
    a = _aristas(triangulos)
    claves, cuentas = np.unique(a.reshape(-1, 6), axis=0, return_counts=True)
    uno = int(np.count_nonzero(cuentas == 1))
    dos = int(np.count_nonzero(cuentas == 2))
    mas = int(np.count_nonzero(cuentas > 2))
    detalle = {
        "aristas_una_cara": uno,      # agujeros
        "aristas_dos_caras": dos,     # las buenas
        "aristas_mas_de_dos": mas,    # superficie que se toca a sí misma
        "aristas_totales": int(claves.shape[0]),
    }
    return (uno == 0 and mas == 0), detalle


def volumen_con_signo(malla: MallaSTL | np.ndarray) -> float:
    """Volumen por divergencia: positivo si las normales apuntan hacia afuera."""
    t = malla.triangulos if isinstance(malla, MallaSTL) else np.asarray(malla, dtype=float)
    a, b, c = t[:, 0], t[:, 1], t[:, 2]
    return float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0)


def verificar_stl(malla: MallaSTL, nombre: str = "STL") -> dict:
    """Valida que la superficie sirva para mallar un volumen. Lanza con el motivo concreto si no."""
    cerrada, detalle = escerrada(malla.triangulos)
    if not cerrada:
        partes = []
        if detalle["aristas_una_cara"]:
            partes.append(
                f"{detalle['aristas_una_cara']} arista(s) con una sola cara (agujero en la superficie)"
            )
        if detalle["aristas_mas_de_dos"]:
            partes.append(
                f"{detalle['aristas_mas_de_dos']} arista(s) con más de dos caras (la superficie se toca "
                f"a sí misma)"
            )
        raise ValueError(
            f"{nombre}: la malla no es cerrada y no se puede mallar un volumen — " + "; ".join(partes)
            + f" (de {detalle['aristas_totales']} aristas)"
        )
    vol = volumen_con_signo(malla)
    if vol <= 0:
        raise ValueError(
            f"{nombre}: las normales apuntan hacia adentro (volumen por divergencia {vol * 1e9:.3f} mm³ "
            f"≤ 0). Hay que invertir el orden de los vértices de cada cara; mallar así daría un sólido "
            f"con masa negativa y el acople saldría con el signo cambiado."
        )
    return {"cerrada": True, "volumen_m3": vol, **detalle}


def seccion_por_plano(malla: MallaSTL, eje: str = "z", valor: float = 0.0,
                      tol: float = 1e-12) -> list[np.ndarray]:
    """Corta la superficie con un plano y devuelve los **contornos cerrados** ``(k, 2)`` en el plano.

    El simulador es 2D: un STL describe una superficie en 3D, así que lo que entra al modelo es **la
    sección** —el mismo criterio que ya está documentado para las aperturas, que se modelan como su corte
    lineal—. Esto corta y devuelve polígonos cerrados, listos para mallar.

    Cada triángulo que cruza el plano aporta un segmento; los segmentos se encadenan por extremos. Se
    devuelve una lista porque una sección puede tener varias piezas (o agujeros).
    """
    idx = {"x": 0, "y": 1, "z": 2}
    if eje not in idx:
        raise ValueError(f"eje {eje!r} no es 'x', 'y' ni 'z'")
    j = idx[eje]
    otros = [k for k in (0, 1, 2) if k != j]

    t = np.asarray(malla.triangulos, dtype=float)
    d = t[:, :, j] - valor

    segmentos: list[tuple[np.ndarray, np.ndarray]] = []
    for tri, dd in zip(t, d):
        # sin cruce estricto: los tres del mismo lado (o tocando el plano sin cruzarlo)
        if np.all(dd > 0) or np.all(dd < 0):
            continue
        cortes = []
        for a, b in ((0, 1), (1, 2), (2, 0)):
            if (dd[a] > 0) != (dd[b] > 0):
                s = dd[a] / (dd[a] - dd[b])
                cortes.append(tri[a] + s * (tri[b] - tri[a]))
            elif abs(dd[a]) <= tol:
                cortes.append(tri[a])
        if len(cortes) < 2:
            continue
        p0 = np.array([cortes[0][k] for k in otros])
        p1 = np.array([cortes[1][k] for k in otros])
        if np.linalg.norm(p1 - p0) > tol:
            segmentos.append((p0, p1))

    if not segmentos:
        raise ValueError(
            f"el plano {eje} = {valor} no corta la malla: la superficie queda entera de un lado"
        )

    # encadenar segmentos por extremos redondeados
    clave = lambda p: (round(p[0], 10), round(p[1], 10))
    pendientes: dict[tuple, list[tuple]] = {}
    for p0, p1 in segmentos:
        pendientes.setdefault(clave(p0), []).append((p0, p1))

    contornos: list[np.ndarray] = []
    usados: set[tuple] = set()
    for p0, p1 in segmentos:
        k0 = clave(p0)
        if (k0, clave(p1)) in usados:
            continue
        cadena = [p0]
        actual = p1
        usados.add((k0, clave(p1)))
        while True:
            k = clave(actual)
            cadena.append(actual)
            if k == clave(cadena[0]) and len(cadena) > 2:
                break
            siguiente = None
            for par in pendientes.get(k, []):
                marca = (k, clave(par[1]))
                if marca not in usados:
                    siguiente = par
                    usados.add(marca)
                    break
            if siguiente is None:
                break
            actual = siguiente[1]
        if len(cadena) > 3:
            contornos.append(np.array(cadena[:-1] if np.allclose(cadena[0], cadena[-1]) else cadena))

    if not contornos:
        raise ValueError("no se pudo cerrar ningún contorno: la superficie no es cerrada")
    return contornos


def area_contorno(contorno: np.ndarray) -> float:
    """Área con signo por la fórmula del polígono (positiva si el contorno va en sentido antihorario)."""
    c = np.asarray(contorno, dtype=float)
    x, y = c[:, 0], c[:, 1]
    return float(0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def cubo_stl(lado: float = 1.0, centro=(0.0, 0.0, 0.0)) -> MallaSTL:
    """Cubo cerrado con normales hacia afuera, para tests y para el botón de ejemplo."""
    h = lado / 2.0
    cx, cy, cz = centro
    v = np.array([
        [-h, -h, -h], [h, -h, -h], [h, h, -h], [-h, h, -h],
        [-h, -h, h], [h, -h, h], [h, h, h], [-h, h, h],
    ]) + np.array([cx, cy, cz])
    caras = [
        (0, 3, 2), (0, 2, 1),   # z = -h
        (4, 5, 6), (4, 6, 7),   # z = +h
        (0, 1, 5), (0, 5, 4),   # y = -h
        (3, 7, 6), (3, 6, 2),   # y = +h
        (0, 4, 7), (0, 7, 3),   # x = -h
        (1, 2, 6), (1, 6, 5),   # x = +h
    ]
    tri = np.array([[v[i], v[j], v[k]] for i, j, k in caras], dtype=float)
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    normales = n / np.linalg.norm(n, axis=1, keepdims=True)
    return MallaSTL(triangulos=tri, normales=normales, origen="generado")
