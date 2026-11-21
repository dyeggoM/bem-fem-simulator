"""La cadena piel → hueso → cóclea resuelta por FEM y contrastada contra la referencia (U9).

Qué hace
--------
Arma un dominio 2D rectangular con agua a la entrada, las capas pedidas en el medio y agua a la salida,
incide una onda plana normal y resuelve con `fem.multicapa.solve_multicapa`. Devuelve el mapa de presión,
el corte axial y **el coeficiente de transmisión medido**, comparado contra la matriz de transferencia
cuando la comparación es válida.

Cuándo NO hay referencia
------------------------
Con alguna capa sólida elástica (hueso) hay cizalla y la fórmula de capas fluidas no vale. En ese caso
`t_teorico` sale `None` y la UI dice «sin referencia analítica», que es la verdad: no se compara contra un
número que no corresponde al mismo problema.
"""

from __future__ import annotations

import time
from typing import cast

import numpy as np

from .. import cadena as materiales
from ..materials import Material
from ..fem import helmholtz2d, multicapa
from ..meshing import boundary_edges, edge_midpoints, grid_tri_mesh

__all__ = ["resolver_cadena", "opciones_cadena"]

#: Presupuesto de nodos para la malla de la cadena. El solve corre en la cola (background), así que el
#: tope no es por bloquear la UI sino por no lanzar una corrida de minutos sin avisar. Medido: ~93 000
#: nodos es ~10 s de solve; 150 000 queda en el orden de 30 s.
MAX_NODOS = 150_000


def opciones_cadena() -> dict:
    return {
        "medios": materiales.resumen_medios(),
        "aviso": (
            "La mayoría de estos valores son parámetros de MODELADO publicados, no mediciones de cada "
            "tejido: no existe un conjunto experimental completo validado de 0,1 a 1 MHz. Piel, grasa, "
            "músculo y líquidos entran como FLUIDOS porque no hay c_S citable para tejido blando; el "
            "hueso entra como sólido elástico."
        ),
    }


def _capas(datos: dict) -> list[tuple[str, float]]:
    crudas = datos.get("capas") or []
    if not crudas:
        raise ValueError("hacen falta 'capas': [{'medio': ..., 'espesor_mm': ...}, ...]")
    salida = []
    for i, capa in enumerate(crudas):
        medio = str(capa.get("medio", "")).strip()
        if not medio:
            raise ValueError(f"la capa {i + 1} no dice el medio")
        espesor = float(capa.get("espesor_mm", 0.0)) * 1e-3
        if not (espesor > 0):
            raise ValueError(f"la capa {i + 1} ({medio}) tiene espesor {capa.get('espesor_mm')!r} mm")
        salida.append((medio, espesor))
    return salida


def resolver_cadena(datos: dict) -> dict:
    t0 = time.perf_counter()
    capas = _capas(datos)
    f = float(datos.get("frecuencia_hz", 500e3))
    ppw = float(datos.get("ppw", 14.0))
    agua_mm = float(datos.get("agua_mm", 12.0))
    alto_mm = float(datos.get("alto_mm", 24.0))
    nx = int(datos.get("nx", 140))
    ny = int(datos.get("ny", 90))

    fluidos, solidos = materiales.medios_de([m for m, _ in capas] + ["agua"])
    # el agua tiene que estar como fluido; si el usuario la pone como capa, ya está
    agua = fluidos.get("agua")
    if agua is None:
        agua = cast(Material, materiales.MEDIOS_CADENA["agua"].material)
        fluidos["agua"] = agua

    lam = agua.wavelength(f)
    h = lam / ppw
    entrada = agua_mm * 1e-3
    espesor_total = sum(e for _, e in capas)
    x0_dominio = float(datos.get("estructura_x_mm", agua_mm)) * 1e-3
    salida = max(entrada, 8e-3)
    # el dominio se arma de modo que la estructura entre completa, con agua a los dos lados
    largo = max(x0_dominio + espesor_total + salida, entrada + espesor_total + salida)

    n_largo = int(round(largo / h))
    n_alto = int(round(2 * (alto_mm * 1e-3) / h))
    # En vez de rechazar la corrida, se baja el ppw hasta que entre: el usuario pidió una frecuencia y un
    # dominio, no un número de elementos por longitud de onda. Se avisa qué se hizo y con cuánto quedó,
    # porque un mapa más grueso puede no resolver la capa más fina.
    aviso_malla = None
    if n_largo * n_alto > MAX_NODOS:
        ppw_original = ppw
        while n_largo * n_alto > MAX_NODOS and ppw > 4.0:
            ppw = ppw * 0.9
            h = lam / ppw
            n_largo = int(round(largo / h))
            n_alto = int(round(2 * (alto_mm * 1e-3) / h))
        if n_largo * n_alto > MAX_NODOS:
            raise ValueError(
                f"a {f / 1e3:.0f} kHz el dominio necesita ~{n_largo * n_alto} nodos incluso con "
                f"4 elementos por longitud de onda (tope {MAX_NODOS}). Bajá el tamaño del dominio "
                f"(agua o alto) o la frecuencia."
            )
        aviso_malla = (
            f"se bajó de {ppw_original:.1f} a {ppw:.1f} elementos por longitud de onda para entrar en "
            f"{MAX_NODOS} nodos; la capa más fina puede quedar subresuelta"
        )
    xs = np.linspace(0.0, largo, n_largo + 1)
    ys = np.linspace(-alto_mm * 1e-3, alto_mm * 1e-3, n_alto + 1)
    malla = grid_tri_mesh(xs, ys)
    nodos, tris = malla.nodes, malla.triangles
    t_malla = time.perf_counter()

    centro_x = nodos[tris][:, :, 0].mean(axis=1)
    centro_y = nodos[tris][:, :, 1].mean(axis=1)
    medio_elem = np.full(centro_x.shape, "agua", dtype=object)
    corte = []

    # La estructura puede ser una **franja finita** dentro del agua (se puede mover) o atravesar todo el
    # dominio (capas infinitas, el caso de una interfaz plana). Con `estructura_alto_mm` en 0 o mayor que
    # el dominio, se comporta como antes.
    alto_est = float(datos.get("estructura_alto_mm", 0.0)) * 1e-3
    x0_est = float(datos.get("estructura_x_mm", agua_mm))
    y0_est = float(datos.get("estructura_y_mm", 0.0)) * 1e-3
    finita = 0.0 < alto_est < 2.0 * alto_mm * 1e-3

    x = x0_est * 1e-3
    for nombre, espesor in capas:
        dentro = (centro_x >= x) & (centro_x < x + espesor)
        if finita:
            dentro &= np.abs(centro_y - y0_est) <= alto_est / 2.0
        medio_elem[dentro] = nombre
        corte.append((nombre, x, x + espesor))
        x += espesor
    fin_estructura = x

    ar = boundary_edges(tris)
    mid = edge_midpoints(nodos, ar)

    # ---- excitación: onda plana o un arreglo de fuentes --------------------
    crudas = datos.get("fuentes") or []
    if crudas:
        from ..arreglos import cargar_arreglo  # noqa: F401  (documenta el formato aceptado)

        # Las fuentes son aperturas sobre el baffle en x = 0: el dominio deja de tener ABC en la entrada
        # (ahí hay baffle rígido con las aperturas) y pasa a tenerlas en la salida y en las laterales,
        # porque una fuente localizada radia en todas direcciones.
        en_baffle = np.abs(mid[:, 0]) < 1e-12
        abc = ar[~en_baffle]
        neumann = ar[en_baffle]
        medios_baffle = edge_midpoints(nodos, neumann)
        v = np.zeros(neumann.shape[0], dtype=complex)
        fuera_de_apertura = 0
        for j, (cx, cy) in enumerate(medios_baffle):
            for s in crudas:
                y_s = float(s["y_mm"]) * 1e-3
                ancho = float(s["ancho_mm"]) * 1e-3
                if abs(cy - y_s) <= 0.5 * ancho + 1e-12:
                    v[j] = float(s.get("amp", 1.0)) * np.exp(1j * float(s.get("fase_rad", 0.0)))
                    break
        if not np.any(v):
            fuera_de_apertura = int(neumann.shape[0])
            raise ValueError(
                "ninguna fuente cae sobre el baffle de la malla: revisá las posiciones (y en mm) y que "
                "estén dentro del alto del dominio"
            )
        # La forma débil suma ∫ (1/ρ)·(∂p/∂n)·q ds, y ∂p/∂n = +i·ω·ρ·v (n saliente = −x̂, la fuente empuja
        # hacia +x). El ρ se cancela: se pasa i·ω·v y no depende del medio.
        dp_dn = 1j * 2.0 * np.pi * f * v
        r = multicapa.solve_multicapa(
            nodos, tris, medio_elem, fluidos, solidos, f, abc_edges=abc,
            neumann_edges=neumann, dp_dn=dp_dn,
        )
        excitacion = "fuentes"
        n_fuentes = len(crudas)
    else:
        # ABC sólo en entrada y salida: en las paredes laterales la onda plana tiene ∂p/∂n = 0 y una ABC
        # ahí refleja (D56)
        abc = ar[(np.abs(mid[:, 0]) < 1e-12) | (np.abs(mid[:, 0] - largo) < 1e-12)]
        onda = helmholtz2d.plane_wave(f, agua, (1.0, 0.0))
        r = multicapa.solve_multicapa(
            nodos, tris, medio_elem, fluidos, solidos, f, abc_edges=abc, onda=onda
        )
        excitacion = "onda plana"
        n_fuentes = 0
    t_solve = time.perf_counter()

    # ---- referencia analítica (sólo con onda plana y sin sólidos) ----------
    # La matriz de transferencia es para una onda plana. Con un arreglo de fuentes el campo no es una
    # onda plana, así que no hay contra qué comparar y no se inventa una.
    t_teorico = (materiales.t_capas_analitico(capas, f)
                 if (not solidos and excitacion == "onda plana") else None)
    hay_solido = bool(solidos)

    # ---- campo transmitido medido ------------------------------------------
    # Banda de una celda de ancho, no la línea exacta: en una malla no estructurada |y| < 1e-12 toma
    # 2 o 3 nodos y la media no dice nada (el mismo error que costó 45 % en D60).
    # Banda de UNA celda a cada lado. Con media celda quedaba vacía cuando el eje caía justo entre dos
    # filas de nodos (ymin = h/2 exacto), y el corte axial salía sin puntos.
    eje = np.abs(nodos[:, 1]) <= h
    # con sólidos, los nodos del sólido no tienen gdl de presión: se toman sólo los que sí
    con_p = np.asarray(r.dof_p, dtype=int)[eje] >= 0
    x_m = nodos[eje, 0][con_p]
    p_m = np.abs(r.p_todos[eje])[con_p]
    # la zona de medición arranca 4 mm después de la última interfaz, para no leer el campo cercano
    zona = x_m > (fin_estructura + 0.004)
    # Con onda plana el campo detrás vale |t|·|p_incidente| y la media ES el coeficiente. Con fuentes el
    # campo no es una onda plana y ese cociente no significa nada: se reporta el pico detrás, en la misma
    # escala relativa que el mapa.
    if excitacion == "fuentes":
        t_medido = None
        pico_detras = float(p_m[zona].max()) if zona.any() else None
    else:
        t_medido = float(p_m[zona].mean()) if zona.any() else None
        pico_detras = None

    # ---- mapa --------------------------------------------------------------
    X, Y = np.meshgrid(np.linspace(0, largo, nx), np.linspace(-alto_mm * 1e-3, alto_mm * 1e-3, ny))
    from matplotlib.tri import LinearTriInterpolator, Triangulation

    tri = Triangulation(nodos[:, 0], nodos[:, 1], tris)
    dentro = tri.get_trifinder()(X, Y) >= 0
    # El campo es un fasor: hay que interpolar PARTE REAL Y PARTE IMAGINARIA por separado.
    # `LinearTriInterpolator(tri, z_complejo)` tira la parte imaginaria con un ComplexWarning y el mapa
    # terminaba mostrando |Re(p)|, que no es la amplitud.
    re = np.asarray(LinearTriInterpolator(tri, np.real(r.p_todos))(X, Y), dtype=float)
    im = np.asarray(LinearTriInterpolator(tri, np.imag(r.p_todos))(X, Y), dtype=float)
    campo = np.hypot(re, im)
    campo = np.where(dentro, campo, np.nan)
    pico = np.nanmax(campo) if np.isfinite(campo).any() else 1.0
    db = 20.0 * np.log10(np.maximum(campo / (pico or 1.0), 1e-6))

    t_mapa = time.perf_counter()

    return {
        "ok": True,
        "capas": [
            {"medio": nombre, "etiqueta": materiales.MEDIOS_CADENA[nombre].etiqueta,
             "desde_mm": round(a * 1e3, 3), "hasta_mm": round(b * 1e3, 3),
             "espesor_mm": round((b - a) * 1e3, 3)}
            for nombre, a, b in corte
        ],
        "estructura": {
            "x_mm": round(corte[0][1] * 1e3, 3) if corte else 0.0,
            "y_mm": round(float(datos.get("estructura_y_mm", 0.0)), 3),
            "ancho_mm": round((fin_estructura - (corte[0][1] if corte else 0.0)) * 1e3, 3),
            "alto_mm": (round(alto_est * 1e3, 3) if finita else None),
            "finita": finita,
        },
        "metricas": {
            "t_medido": t_medido,
            "pico_detras": pico_detras,
            "t_teorico": t_teorico,
            "error_pct": (None if (t_teorico is None or t_medido is None or t_teorico == 0)
                          else round(100 * abs(t_medido - t_teorico) / t_teorico, 2)),
            "hay_fuentes": excitacion == "fuentes",
            "hay_solido": hay_solido,
            "referencia": (
                "sin referencia: se excita con un arreglo de fuentes, no con una onda plana"
                if excitacion == "fuentes"
                else (
                    "matriz de transferencia de capas fluidas (la estructura es finita: la referencia "
                    "vale para capas infinitas, así que el error incluye la difracción de los bordes)"
                    if t_teorico is not None and finita
                    else "matriz de transferencia de capas fluidas" if t_teorico is not None
                    else "sin referencia analítica: hay capa sólida y la fórmula de capas fluidas no "
                         "incluye cizalla"
                )
            ),
            "frecuencia_hz": f,
            "ppw": ppw,
            "nodos": int(nodos.shape[0]),
            "elementos": int(tris.shape[0]),
            "espesor_total_mm": round(espesor_total * 1e3, 3),
            "lambda_agua_mm": round(lam * 1e3, 3),
            "ppw_efectivo": round(ppw, 2),
            "aviso_malla": aviso_malla,
            "excitacion": excitacion,
            "n_fuentes": n_fuentes,
        },
        "grilla": {
            "x_mm": np.round(np.linspace(0, largo, nx) * 1e3, 4).tolist(),
            "y_mm": np.round(np.linspace(-alto_mm, alto_mm, ny), 4).tolist(),
            "db": np.round(db, 3).tolist(),
        },
        "corte_axial": {
            "x_mm": np.round(x_m * 1e3, 4).tolist(),
            "abs_p": np.round(p_m, 6).tolist(),
        },
        "medios": {n: {"etiqueta": m.etiqueta, "tipo": m.tipo, "fuente": m.fuente,
                       "verificado": m.verificado, "nota": m.nota}
                   for n, m in materiales.MEDIOS_CADENA.items()},
        "tiempos": {
            "malla_s": round(t_malla - t0, 3),
            "solve_s": round(t_solve - t_malla, 3),
            "pos_s": round(t_mapa - t_solve, 3),
        },
    }
