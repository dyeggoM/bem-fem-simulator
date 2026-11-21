"""Escena 3D por capas con rebanadas (Fase U4).

Qué resuelve
------------
Un apilado de capas en una caja 3D (eje de propagación *x*), con el núcleo P1 en tetraedros de
`fem.helmholtz3d`, ABC en las dos caras de *x* y paredes rígidas en *y*/*z* — o sea incidencia
normal, igual que la escena 2D de U1 pero en 3D.

Por qué rebanadas y no volumen: a 18 081 nodos (D32) el campo es un sólido y el trazado completo no
se lee. La UI muestra el plano *x-y* a mitad de *z* y el *x-z* a mitad de *y*, que es lo que deja ver
la onda estacionaria en las capas y el frente plano.

El techo está medido (D32) y se declara antes de resolver: 5 λ × 1 λ × 1 λ a ppw 20 son 18 081 nodos
y ~5 s; un cubo de 3 λ a ppw 12 no termina en 10 min con este ensamblado. No es memoria, es el costo
del solve denso por BLAS. Por eso hay tope de nodos con mensaje, no un cuelgue.
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np

from .. import materials as mats
from ..analytic import helmholtz1d
from ..fem import helmholtz3d as h3
from ..meshing3d import boundary_faces, faces_on_plane, grid_tet_mesh
from . import validacion

MATERIALES: dict[str, Any] = {
    "agua": mats.WATER,
    "tejido": mats.SOFT_TISSUE,
    "hueso": mats.BONE_CORTICAL,
    "perilinfa": mats.PERILYMPH,
    "endolinfa": mats.ENDOLYMPH,
}

ABC_DISPONIBLES = ("abc1", "bgt2")
MAX_NODOS = 90_000

_LIMITES = {
    "frecuencia_hz": (20.0e3, 400.0e3),
    "ppw": (6.0, 24.0),
    "alto_m": (0.5e-3, 45.0e-3),   # 3 λ a 100 kHz: arriba de eso el tope de nodos dice que no
    "capas": (2, 6),
    "espesor_m": (0.2e-3, 60.0e-3),
}


def _num(datos: dict[str, Any], clave: str, defecto: float) -> float:
    valor = datos.get(clave, defecto)
    try:
        numero = float(valor)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{clave} tiene que ser un número; llegó {valor!r}") from exc
    if not np.isfinite(numero):
        raise ValueError(f"{clave} no puede ser {numero}")
    return numero


def _en_rango(clave: str, valor: float) -> float:
    bajo, alto = _LIMITES[clave]
    if not (bajo <= valor <= alto):
        raise ValueError(f"{clave} = {valor:g} fuera de rango [{bajo:g}, {alto:g}]")
    return valor


def validar_volumen(datos: dict[str, Any]) -> dict[str, Any]:
    """Normaliza y valida el pedido 3D. Falla con mensaje, nunca con un `KeyError`."""
    if not isinstance(datos, dict):
        raise ValueError("el pedido tiene que ser un objeto JSON")

    frecuencia = _en_rango("frecuencia_hz", _num(datos, "frecuencia_hz", 100.0e3))
    ppw = _en_rango("ppw", _num(datos, "ppw", 12.0))
    alto = _en_rango("alto_m", _num(datos, "alto_m", 1.5e-3))

    fuente = str(datos.get("fuente", "puntual"))
    if fuente not in ("plano", "puntual", "esfera"):
        raise ValueError(
            f"fuente desconocida: {fuente!r}; opciones ['plano', 'puntual', 'esfera']"
        )

    crudas = datos.get("capas")
    if not isinstance(crudas, list):
        raise ValueError("el pedido necesita una lista 'capas'")
    _, alto_capas = _LIMITES["capas"]
    # puntual y esfera van en un medio homogéneo: una sola capa es válida
    bajo = 1 if fuente in ("puntual", "esfera") else _LIMITES["capas"][0]
    if not (bajo <= len(crudas) <= alto_capas):
        raise ValueError(f"'capas' tiene que tener entre {bajo} y {alto_capas} elementos")

    capas: list[tuple[str, float]] = []
    for i, capa in enumerate(crudas):
        if not isinstance(capa, dict):
            raise ValueError(f"la capa {i} tiene que ser {{material, espesor_m}}")
        nombre = capa.get("material")
        if nombre not in MATERIALES:
            raise ValueError(
                f"capa {i}: material desconocido {nombre!r}; opciones: {sorted(MATERIALES)}"
            )
        espesor = _en_rango("espesor_m", _num(capa, "espesor_m", 0.0))
        capas.append((str(nombre), espesor))

    abc = str(datos.get("abc", "abc1"))
    if abc not in ABC_DISPONIBLES:
        raise ValueError(f"abc desconocido: {abc!r}; opciones {list(ABC_DISPONIBLES)}")

    if fuente == "esfera":
        ka = _num(datos, "ka", 1.0)
        if not (0.2 <= ka <= 6.0):
            raise ValueError(f"ka = {ka:g} fuera de [0.2, 6]")
        material_esfera = str(datos.get("material_esfera", "hueso"))
        if material_esfera not in MATERIALES:
            raise ValueError(
                f"material de la esfera desconocido: {material_esfera!r}; "
                f"opciones {sorted(MATERIALES)}"
            )
        return {
            "frecuencia_hz": frecuencia,
            "ppw": ppw,
            "alto_m": alto,
            "capas": tuple(capas),
            "abc": abc,
            "fuente": "esfera",
            "ka": ka,
            "material_esfera": material_esfera,
        }

    if fuente == "puntual" and len({c.get("material") for c in crudas}) > 1:
        raise ValueError(
            "la fuente puntual de U4 va en un medio homogéneo: con capas de distinto material el "
            "Green 3D deja de ser la referencia y el contraste de R/T se mide con 'plano'"
        )

    angulo = abs(_num(datos, "angulo_deg", 0.0))
    if angulo > 1e-9:
        raise ValueError(
            "la escena 3D de U4 resuelve sólo incidencia normal: en y/z las paredes son rígidas y "
            "una onda oblicua se reflejaría ahí en vez de salir."
        )

    return {
        "frecuencia_hz": frecuencia,
        "ppw": ppw,
        "alto_m": alto,
        "capas": tuple(capas),
        "abc": abc,
        "fuente": fuente,
    }


def _malla(p: dict[str, Any]):
    materiales = [MATERIALES[nombre] for nombre, _ in p["capas"]]
    espesores = [espesor for _, espesor in p["capas"]]
    pasos = [abs(m.wavelength(p["frecuencia_hz"])) / p["ppw"] for m in materiales]

    xs: list[float] = [0.0]
    for espesor, paso in zip(espesores, pasos):
        n = max(1, int(np.ceil(espesor / paso)))
        xs.extend(np.linspace(xs[-1], xs[-1] + espesor, n + 1)[1:].tolist())
    if p["fuente"] in ("puntual", "esfera"):
        # La caja va CENTRADA en el origen. El caso puntual lo necesita porque la fuente es el nodo
        # más cercano al origen y quedaría sobre la ABC (medido: 48 % de error que no bajaba al
        # refinar). La esfera lo necesita por simetría: su centro tiene que ser el origen.
        centro = xs[-1] / 2.0
        xs = [x - centro for x in xs]

    paso_transversal = min(pasos)
    n_t = max(1, int(np.ceil(p["alto_m"] / paso_transversal)))
    ys = np.linspace(-p["alto_m"] / 2.0, p["alto_m"] / 2.0, n_t + 1)
    zs = np.linspace(-p["alto_m"] / 2.0, p["alto_m"] / 2.0, n_t + 1)

    malla = grid_tet_mesh(xs, ys, zs)
    if malla.n_nodes > MAX_NODOS:
        raise ValueError(
            f"la malla sale de {malla.n_nodes} nodos y el tope interactivo es {MAX_NODOS}. "
            f"Bajá ppw (hoy {p['ppw']:g}), el alto (hoy {p['alto_m'] * 1e3:g} mm) o el espesor "
            f"total. Para más que esto hace falta FEniCSx+MPI o resolver el exterior por BEM (D32)."
        )
    centro_x = malla.nodes[malla.tets].mean(axis=1)[:, 0]
    bordes = np.cumsum(espesores)
    capa = np.clip(np.searchsorted(bordes, centro_x, side="right"), 0, len(espesores) - 1)
    return np.asarray(xs), ys, zs, malla, capa


def estimar_volumen(datos: dict[str, Any]) -> dict[str, Any]:
    """Nodos, elementos y costo estimado sin resolver: lo que la UI dice antes de apretar."""
    p = validar_volumen(datos)
    xs, ys, zs, malla, _ = _malla(p)
    n = int(malla.n_nodes)
    # medido en D32: 18 081 nodos en 5.3 s con este ensamblado
    return {
        "ok": True,
        "malla": {
            "n_nodos": n,
            "n_elementos": int(malla.tets.shape[0]),
            "nx": int(xs.size),
            "ny": int(ys.size),
            "nz": int(zs.size),
            # con fuente puntual la malla va centrada: xs[-1] es media longitud, no la longitud
            "longitud_mm": float((xs[-1] - xs[0]) * 1e3),
            "alto_mm": float(p["alto_m"] * 1e3),
            "paso_x_min_mm": float(np.min(np.diff(xs)) * 1e3),
            "interfaces_mm": (np.cumsum([e for _, e in p["capas"]])[:-1] * 1e3).tolist(),
        },
        "estimado_solve_s": round(5.3 * (n / 18_081.0) ** 1.5, 1),
        "tope_nodos": MAX_NODOS,
        "techo_conocido": "un cubo de 3 λ a ppw 12 no termina en 10 min (D32)",
    }


def _ajustar(x: np.ndarray, p: np.ndarray, k: complex, dos_ondas: bool) -> np.ndarray:
    base = (np.column_stack([np.exp(1j * k * x), np.exp(-1j * k * x)]) if dos_ondas
            else np.exp(1j * k * x)[:, None])
    coef, *_ = np.linalg.lstsq(base, p, rcond=None)
    return np.asarray(coef).reshape(-1)


def resolver_volumen(datos: dict[str, Any]) -> dict[str, Any]:
    """Resuelve la caja 3D y devuelve dos rebanadas, R/T y los tiempos."""
    t_inicio = time.perf_counter()
    p = validar_volumen(datos)
    xs, ys, zs, malla, capa = _malla(p)
    f = p["frecuencia_hz"]
    materiales = [MATERIALES[nombre] for nombre, _ in p["capas"]]

    k_el = np.array([materiales[c].wavenumber(f) for c in capa], dtype=complex)
    rho_el = np.array([materiales[c].rho for c in capa], dtype=float)

    caras = boundary_faces(malla.tets)
    x_max = float(xs[-1])
    t_solve = time.perf_counter()
    # sólo se llenan en el modo puntual; se inicializan acá para que la rama de métricas no dependa
    # de qué rama eligió el intérprete
    err = float("nan")
    amplitud = float("nan")
    anillo = np.zeros(malla.nodes.shape[0], dtype=bool)
    radio = 0.0
    dentro = np.zeros(malla.tets.shape[0], dtype=bool)
    ok = np.zeros(1, dtype=bool)

    if p["fuente"] == "esfera":
        # Dispersión penetrable: el contraste es una esfera de otro material. Necesita el operador
        # local y el término de contraste con el fondo (D41), no el campo total.
        from ..analytic import sphere as esfera_analitica
        from scipy.interpolate import LinearNDInterpolator

        material_esfera = MATERIALES[p["material_esfera"]]
        k0 = materiales[0].wavenumber(f)
        rho0 = materiales[0].rho
        radio = p["ka"] / abs(k0)
        centroides = malla.nodes[malla.tets].mean(axis=1)
        dentro = np.linalg.norm(centroides, axis=1) < radio
        if not dentro.any():
            raise ValueError(
                f"ningún elemento cayó dentro de la esfera de radio {radio * 1e3:.2f} mm: "
                f"subí ppw (hoy {p['ppw']:g}) para que la malla la resuelva"
            )
        k_el = np.where(dentro, material_esfera.wavenumber(f), k0)
        rho_el = np.where(dentro, material_esfera.rho, rho0)
        abc = caras
        p_inc = h3.plane_wave(f, materiales[0], (0.0, 0.0, 1.0))
        p_s = h3.solve_scattered(
            malla.nodes, malla.tets, k_el, abc_faces=abc, element_rho=rho_el,
            background_k=k0, background_rho=rho0, p_inc=p_inc,
        )
        total = h3.total_field(p_s, p_inc, malla.nodes)
        # Los mapas muestran el campo DISPERSADO, no el total: en el total manda la incidente
        # (|p| ≈ 1 en todas partes) y una perturbación del 5 % es invisible en una escala de 30 dB.
        # La métrica sí se mide sobre el total, que es lo que da la serie de Anderson.
        campo = p_s
        dt_solve = time.perf_counter() - t_solve

        # error contra la serie de Anderson sobre r = 1.5a
        th = np.linspace(0.0, np.pi, 48)
        ph = np.linspace(0.0, 2.0 * np.pi, 20, endpoint=False)
        TH, PH = np.meshgrid(th, ph, indexing="ij")
        rr = 1.5 * radio
        pts = np.column_stack([
            (rr * np.sin(TH) * np.cos(PH)).ravel(),
            (rr * np.sin(TH) * np.sin(PH)).ravel(),
            (rr * np.cos(TH)).ravel(),
        ])
        serie = esfera_analitica.fluid_sphere_scattering(
            f, materiales[0], material_esfera, radio, pts
        )
        interp = LinearNDInterpolator(
            malla.nodes, np.column_stack([np.real(total), np.imag(total)]), fill_value=np.nan
        )
        vals = interp(pts)
        ok = ~np.isnan(vals[:, 0])
        fem = vals[ok, 0] + 1j * vals[ok, 1]
        ref = np.asarray(serie["total"])[ok]
        err = (float(np.linalg.norm(fem - ref) / np.linalg.norm(ref))
               if ok.sum() > 8 else float("nan"))
        amplitud = float(np.abs(p_s).max())
        anillo = np.zeros(malla.nodes.shape[0], dtype=bool)
    elif p["fuente"] == "puntual":
        # Una fuente puntual radia en todas las direcciones: si las caras laterales fueran rígidas,
        # la onda se reflejaría ahí y el mapa sería un patrón de cavidad, no un Green 3D. Por eso
        # acá la ABC va en las SEIS caras.
        abc = caras
        idx = int(np.argmin(np.linalg.norm(malla.nodes, axis=1)))
        fuente = malla.nodes[idx]
        campo = h3.solve(
            malla.nodes, malla.tets, k_el, abc_faces=abc, element_rho=rho_el,
            point_sources=fuente[None, :], point_volume_velocities=np.array([1.0]), frequency=f,
        )
        dt_solve = time.perf_counter() - t_solve

        # validación contra el Green 3D en un anillo, como el test de D32
        radio = np.linalg.norm(malla.nodes - fuente, axis=1)
        lam = abs(materiales[0].wavelength(f))
        anillo = (radio > 0.25 * lam) & (radio < 0.45 * lam)
        omega = 2.0 * np.pi * f
        r_seguro = np.where(radio > 0.0, radio, 1.0)
        exacto = (-1j * omega * materiales[0].rho * np.exp(1j * materiales[0].wavenumber(f) * r_seguro)
                  / (4.0 * np.pi * r_seguro))
        err = (float(np.linalg.norm(np.asarray(campo)[anillo] - exacto[anillo])
                     / np.linalg.norm(exacto[anillo])) if anillo.sum() > 8 else float("nan"))
        amplitud = (float(np.mean(np.abs(np.asarray(campo)[anillo]) * radio[anillo])
                          / (omega * materiales[0].rho / (4 * np.pi))) if anillo.sum() > 8 else float("nan"))
    else:
        i_izq = faces_on_plane(malla.nodes, caras, 0, float(xs[0]))
        i_der = faces_on_plane(malla.nodes, caras, 0, x_max)
        abc = np.vstack([caras[i_izq], caras[i_der]])
        p_inc = h3.plane_wave(f, materiales[0], (1.0, 0.0, 0.0))
        campo = h3.solve(
            malla.nodes, malla.tets, k_el, abc_faces=abc, element_rho=rho_el, p_inc=p_inc,
            abc_order=1 if p["abc"] == "abc1" else 2,
        )
        dt_solve = time.perf_counter() - t_solve

    # ---- rebanadas ----
    forma = (xs.size, ys.size, zs.size)
    volumen = np.asarray(campo, dtype=complex).reshape(forma)
    k_z = zs.size // 2
    k_y = ys.size // 2
    rebanada_xy = np.abs(volumen[:, :, k_z])          # (nx, ny) a z = 0
    rebanada_xz = np.abs(volumen[:, k_y, :])          # (nx, nz) a y = 0
    pico = float(volumen.__abs__().max()) or 1.0

    def _db(matriz):
        with np.errstate(invalid="ignore", divide="ignore"):
            db = 20.0 * np.log10(np.maximum(matriz / pico, 1e-6))
        return np.round(db, 2).tolist()

    # ---- R y T por ajuste pegado a las interfaces (la receta del test 3D, D32) ----
    espesores = [e for _, e in p["capas"]]
    k1 = materiales[0].wavenumber(f)
    k_last = materiales[-1].wavenumber(f)
    t1 = espesores[0]
    x_ini_ultima = float(xs[-1] - espesores[-1])
    nodos_x = malla.nodes[:, 0]
    ventana = 0.5 * abs(materiales[0].wavelength(f))
    tol = 1e-9 * max(x_max, 1.0)
    sel_izq = (nodos_x > t1 - ventana) & (nodos_x < t1 - tol)
    sel_der = (nodos_x > x_ini_ultima + tol) & (nodos_x < x_ini_ultima + ventana)

    metricas: dict[str, Any] = {}
    if p["fuente"] == "esfera":
        metricas = {
            "error_serie_pct": 100.0 * err if np.isfinite(err) else None,
            "radio_mm": radio * 1e3,
            "ka": p["ka"],
            "material_esfera": p["material_esfera"],
            "elementos_en_esfera": int(dentro.sum()),
            "amplitud_dispersada": amplitud,
            "puntos_comparados": int(ok.sum()),
            "nota": "referencia: ondas parciales de Anderson (1950) para esfera fluida, evaluada "
                    "sobre r = 1.5a",
        }
    elif p["fuente"] == "puntual":
        metricas = {
            "error_green_pct": 100.0 * err if np.isfinite(err) else None,
            "amplitud_green": amplitud if np.isfinite(amplitud) else None,
            "nodos_anillo": int(anillo.sum()),
            "nota": "un monopolo de caudal Q tiene que dar p = −iωρQ·e^{ikr}/(4πr); se compara en "
                    "el anillo 0.25–0.45 λ, donde el Green no es singular",
        }
    elif sel_izq.sum() >= 8 and sel_der.sum() >= 8:
        a_izq, b_izq = _ajustar(nodos_x[sel_izq], campo[sel_izq], k1, True)
        c_der = _ajustar(nodos_x[sel_der], campo[sel_der], k_last, False)[0]
        r_num = b_izq / a_izq
        t_num = c_der / a_izq
        respuesta = helmholtz1d.stack_response(f, tuple(materiales), tuple(espesores[:-1][1:]))
        r_an, t_an = complex(respuesta.R), complex(respuesta.T)
        metricas.update({
            "R": abs(r_num), "R_analitico": abs(r_an),
            "T": abs(t_num), "T_analitico": abs(t_an),
            "error_R_abs": float(abs(abs(r_num) - abs(r_an))),
            "error_T_abs": float(abs(abs(t_num) - abs(t_an))),
            "energia": float(abs(r_num) ** 2 + abs(t_num) ** 2),
        })
    elif p["fuente"] == "plano":
        metricas["aviso"] = "capas demasiado finas para la ventana de ajuste (0.5 λ)"

    if p["fuente"] == "esfera":
        sello = validacion.sello_conjunto(
            {"volumen_esfera_serie": (None if metricas.get("error_serie_pct") is None
                                      else metricas["error_serie_pct"] / 100.0)}
        )
    elif p["fuente"] == "puntual":
        sello = validacion.sello_conjunto(
            {"volumen_green": (None if metricas.get("error_green_pct") is None
                               else metricas["error_green_pct"] / 100.0)}
        )
    else:
        sello = validacion.sello_conjunto({
            "volumen_rt_r": metricas.get("error_R_abs"),
            "volumen_rt_t": metricas.get("error_T_abs"),
        })

    return {
        "ok": True,
        "sello": sello,
        "escena": {**p, "capas": [{"material": m, "espesor_m": e} for m, e in p["capas"]]},
        "malla": {
            "n_nodos": int(malla.n_nodes),
            "n_elementos": int(malla.tets.shape[0]),
            "n_caras_abc": int(abc.shape[0]),
        },
        "rebanadas": {
            "xy": {
                "x_mm": np.round(xs * 1e3, 4).tolist(),
                "y_mm": np.round(ys * 1e3, 4).tolist(),
                # .T: Plotly lee z[fila][columna] como [y][x]; la matriz sale [x][y]
                "db": _db(rebanada_xy.T),
                "donde": f"z = {zs[k_z] * 1e3:.2f} mm",
            },
            "xz": {
                "x_mm": np.round(xs * 1e3, 4).tolist(),
                "z_mm": np.round(zs * 1e3, 4).tolist(),
                "db": _db(rebanada_xz.T),
                "donde": f"y = {ys[k_y] * 1e3:.2f} mm",
            },
            "interfaces_mm": (np.cumsum(espesores)[:-1] * 1e3).tolist(),
            "pico": pico,
        },
        # Nube completa para la isosuperficie 3D: la vista dibujaba sólo dos rebanadas porque el
        # bundle de Plotly era `cartesian` (sin gl3d). Con el bundle completo se puede mostrar el
        # volumen. Se manda el módulo normalizado y los nodos en mm: ~9 000 nodos son ~400 KB de JSON.
        "nube": {
            "x_mm": np.round(malla.nodes[:, 0] * 1e3, 4).tolist(),
            "y_mm": np.round(malla.nodes[:, 1] * 1e3, 4).tolist(),
            "z_mm": np.round(malla.nodes[:, 2] * 1e3, 4).tolist(),
            "abs": (np.abs(campo) / (np.abs(campo).max() or 1.0)).tolist(),
            "max_abs": float(np.abs(campo).max()),
        },
        "metricas": metricas,
        "tiempos": {
            "malla_s": round(t_solve - t_inicio, 3),
            "solve_s": round(dt_solve, 3),
            "total_s": round(time.perf_counter() - t_inicio, 3),
        },
    }
