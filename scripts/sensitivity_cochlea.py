#!/usr/bin/env python
"""Análisis de sensibilidad global del corte de cápsula coclear sobre ``t_interior``.

Qué mide
--------
``t_interior`` = ``mean(|p|)`` sobre la circunferencia ``r = a/2`` (la métrica de la Fase 5a),
normalizada por la amplitud de la onda plana incidente. Es **la presión media en el fluido
interior**, no un coeficiente de transmisión: la cavidad es cerrada y sin pérdidas, así que su
campo es una onda estacionaria que recibe energía por toda la circunferencia (D18). Eso manda
sobre cómo se lee todo lo de abajo, y está medido acá (ver «Hallazgo»).

Caso base (el de D22 / ``tests/test_cochlea_impedance.py``)
-----------------------------------------------------------
``a = 20 mm``, ``t = 2 mm``, ``f = 100 kHz``, ``ppw = 20``, ventana de **impedancia** de ``±10°``
centrada en la cara de sombra (``+x``), ``Z_r = ρ_f·c_f``, hueso ``c_p = 3500``,
``c_s = 1900``, ``ρ = 1950`` kg/m³, sin atenuación. Reproduce ``t_interior = 0.972180`` contra
``solve_cross_section`` (diferencia 0.0 a máquina; ver ``tests/test_sensitivity_cochlea.py``).

Método
------
**One-at-a-time (OAT)**: se mueve un parámetro por vez y los demás quedan en el caso base. Se
reporta mínimo, máximo, rango relativo ``(máx−mín)/base`` y el ranking por ese rango. Además se
corre la **grilla completa espesor × atenuación** (los dos que se espera que dominen), que es lo
que después hace falta para decidir si vale la pena agregar la membrana basilar (5b-2).

Rangos y procedencia (todos son *supuestos declarados*, D20)
-------------------------------------------------------------
===================  ==================================  ==============================
parámetro            rango                               fuente / estado
===================  ==================================  ==============================
``c_p`` hueso        1752 / 2500 / 3285 / 3500 m/s       Ammi 2008 (temporal humano,
                                                         1.03 MHz): 1752–3285. 3500 es
                                                         el base, extremo **alto**.
``c_s`` hueso        1360 / 1500 / 1640 / 1900 m/s       White 2006 (calvaria) 1500±140;
                                                         1900 es el base **VERIFY**, no
                                                         medido (D5).
``ρ`` hueso          1868 / 1900 / 1950 kg/m³            Peterson & Dechow 2003.
espesor ``t``        1.0 … 7.0 mm (paso 0.5)            cápsula ótica más fina que el
                                                         cráneo reportado (6–8 mm): el
                                                         extremo bajo es el supuesto.
atenuación           ninguna / White 2006 /               leyes de ``attenuation.py`` (D19).
                     Pichardo 2011 / Pinton 2012
``Z_r`` ventana       0.5 / 1 / 2 × ρ_f·c_f               sin valor citable (D20); se barre
                                                         como en D22 en vez de fijarlo.
===================  ==================================  ==============================

``c_p`` y ``c_s`` no son independientes: ``c_s < c_p`` siempre (si no, ν > 0.5). En el base eso
se cumple (1900 < 3500), pero **no** con ``c_p = 1752``. Por eso el barrido de ``c_p`` se corre
entero con ``c_s = 1500`` fijo (el sustituto de White 2006, dentro de su ±140) y así el rango de
``c_p`` no se contamina con ``c_s``; el efecto de ``c_s`` se mide aparte, a ``c_p = 3500``.

Motor de cálculo (y por qué no es una llamada a ``solve_cross_section``)
-----------------------------------------------------------------------
``ElasticSolid`` no lleva atenuación, así que para las tres leyes hace falta un sólido
viscoelástico: se usan módulos de Lamé complejos con
``k = ω/c + i·α`` (D4) ⇒ ``c̃ = c/(1 + i·α·c/ω)``, ``λ̃ = ρ(c̃_p² − 2c̃_s²)``, ``μ̃ = ρ·c̃_s²``.
El sistema acoplado se ensambla acá con los **bloques públicos** de ``fem`` (``p1_elastic_local``,
``p1_tri_local``, ``abc_boundary_block``, ``abc_incident_rhs``, ``interface_edges``) y el término
de borde del arco ``−∫(i·k_w/ρ_f)·p·q ds`` con ``k_w = ω·ρ_f/Z_r``, exactamente como
``cochlea._assemble_and_solve_arc``; la única diferencia es que la rigidez del sólido se arma con
``D`` complejo. Con ``α = 0`` el ensamblado es el mismo y da ``t_interior`` idéntico a
``solve_cross_section`` (comprobado en el test, tolerancia 1 %).

Dos decisiones de cálculo, medidas:

- **Factorización con ``permc_spec='MMD_AT_PLUS_A'``**: sobre el mismo sistema (45 181 gdl,
  345 377 nnz) SuperLU tarda **43.5 s con el default COLAMD** y **1.9 s con MMD_AT_PLUS_A**
  (misma solución; sólo cambia el pivoteo). Sin eso el barrido no es corrible.
- La atenuación se aplica igual a la onda P y a la S (``α`` único). Las leyes de D19 son
  coeficientes *efectivos* de cráneo (reflexión, dispersión y conversión de modo incluidos en el
  número), así que un sólido homogéneo con esa ``α`` **no** reproduce el mecanismo, sólo su
  magnitud. Es una cota, no un modelo de pérdidas.

Hallazgo (y es el resultado principal)
--------------------------------------
1. **El espesor domina, y por lejos**: rango relativo ~72 % del valor base entre 1 y 7 mm.
2. **``t_interior`` NO es monótona en el espesor**: 1.0996 / 0.9492 / 0.9722 / 1.0629 / 1.1839 /
   1.2691 / 1.3649 / 1.3167 / 1.1920 / 1.0405 / 0.8968 / 0.7683 / 0.6635 para 1.0 … 7.0 mm
   (paso 0.5). No es un artefacto de malla (persiste a ppw=30: 0.9805 / 1.2938 / 0.8533 en
   2/4/6 mm) ni un bug del ensamblado (α=0 reproduce el solver del repo al dígito): es la cavidad
   resonante de D18. El armónico que sí se comporta como una transmisión es el **n=0
   (respiración)**, ``|mean(p)|``: baja 5.1× de 0.3866 a 0.0761, punto a punto salvo un repunte
   final de 0.0757 a 0.0761 (+0.5 % del tramo, dentro del piso numérico de 1.4 %), igual que
   ``|T|`` de la placa equivalente (0.933 → 0.437, estrictamente decreciente).
3. Después del espesor, el que manda es ``c_p`` (~16 % al bajar a 1752 m/s), y ese es el dato que
   más falta: el rango medido de Ammi está **muy por debajo** de nuestro 3500.
4. Atenuación, ``Z_r``, ``c_s`` y ``ρ`` mueven ``t_interior`` menos de 6 % cada uno, del orden del
   piso numérico (cambiar ppw 20→30 mueve el caso base 0.9 % y los casos con pared gruesa hasta
   5 %). A 100 kHz la atenuación es **desconocida** (D20: no hay dato bajo 270 kHz), pero aun la
   extrapolación más agresiva (Pinton, lineal) mueve 5.5 %: no cambia ninguna conclusión.

Uso
---
    .venv/bin/python scripts/sensitivity_cochlea.py                 # barrido completo (~4 min)
    .venv/bin/python scripts/sensitivity_cochlea.py --step 1        # espesor cada 1 mm
    .venv/bin/python scripts/sensitivity_cochlea.py --ppw 12        # grilla reducida
    .venv/bin/python scripts/sensitivity_cochlea.py --abc-order 2   # BGT-2 en el contorno
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

_REPO = Path(__file__).resolve().parents[1]
if (_REPO / "src").is_dir() and str(_REPO / "src") not in sys.path:
    sys.path.insert(0, str(_REPO / "src"))

from ultrasound_sim.attenuation import (  # noqa: E402
    SKULL_PICHARDO_2011_POINTS,
    SKULL_PINTON_2012,
    SKULL_WHITE2006,
)
from ultrasound_sim.cochlea import (  # noqa: E402
    build_cross_section,
    plate_reference,
    window_arc_edges,
)
from ultrasound_sim.fem.assembly import (  # noqa: E402
    p1_tri_local,
    triangle_element_areas,
    triangle_gradients,
    weights_to_local,
)
from ultrasound_sim.fem.coupled2d import (  # noqa: E402
    abc_boundary_block,
    abc_incident_rhs,
    interface_edges,
)
from ultrasound_sim.fem.elastic2d import p1_elastic_local, strain_matrix  # noqa: E402
from ultrasound_sim.fem.helmholtz2d import plane_wave  # noqa: E402
from ultrasound_sim.materials import WATER, ElasticSolid  # noqa: E402

# --------------------------------------------------------------------------- #
# Configuración del caso base
# --------------------------------------------------------------------------- #

FREQUENCY_HZ = 1.0e5
RADIUS_M = 20.0e-3
BASE_THICKNESS_M = 2.0e-3
PPW = 20.0
PAD_WAVELENGTHS = 2.5  # el del caso 0.972180 (D22); 2.0 cambiaría el número base
MIN_ELEMENTS_ACROSS_WALL = 2
WINDOW_HALF_ANGLE_DEG = 10.0
WINDOW_CENTER_ANGLE_DEG = 0.0  # +x = cara de sombra
Z_FLUID = WATER.impedance

#: ``t_interior`` del caso base según ``solve_cross_section`` (D22 / test_cochlea_impedance).
T_INTERIOR_REFERENCE = 0.972180

OUT_CSV = _REPO / "out" / "sensitivity_cochlea.csv"

#: Las tres leyes de D19 más el caso sin pérdidas. ``None`` = material sin atenuación.
LAWS = {
    "sin_atenuacion": None,
    "white_2006": SKULL_WHITE2006,
    "pichardo_2011": SKULL_PICHARDO_2011_POINTS,
    "pinton_2012": SKULL_PINTON_2012,
}

#: Parámetros físicos que entran al ranking (los controles numéricos quedan afuera).
RANKED = ("espesor", "c_p", "atenuacion", "Z_r", "c_s", "rho")

#: Orden de presentación en consola.
ORDER = ("espesor", "c_p", "c_s", "rho", "atenuacion", "Z_r")


# --------------------------------------------------------------------------- #
# Escenario y motor
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Scenario:
    """Un punto del espacio de parámetros del modelo de cápsula."""

    c_p: float = 3500.0
    c_s: float = 1900.0
    rho: float = 1950.0
    thickness_m: float = BASE_THICKNESS_M
    law: str = "sin_atenuacion"
    z_factor: float = 1.0

    def bone(self) -> ElasticSolid:
        return ElasticSolid(
            name=f"bone_cp{self.c_p:.0f}_cs{self.c_s:.0f}_rho{self.rho:.0f}",
            rho=float(self.rho),
            c_p=float(self.c_p),
            c_s=float(self.c_s),
            provenance=(
                "punto del barrido de sensibilidad (D20): c_p Ammi 2008 / c_s White 2006 "
                "(1900 = VERIFY) / ρ Peterson & Dechow 2003"
            ),
        )

    def alpha(self, frequency: float = FREQUENCY_HZ) -> float:
        """Atenuación de amplitud [Np/m] de la ley elegida a ``frequency``."""
        law = LAWS[self.law]
        return 0.0 if law is None else float(law.alpha(frequency))

    def impedance(self) -> complex:
        return complex(self.z_factor * Z_FLUID)

    @property
    def label(self) -> str:
        return f"{self.law}·t={self.thickness_m * 1e3:g}mm·{self.c_p:.0f}/{self.c_s:.0f}/{self.rho:.0f}"


def solve_section(section, *, alpha: float = 0.0, abc_order: int = 1):
    """Resuelve el corte acoplado con la rigidez del sólido posiblemente compleja.

    Es el mismo ensamblado de ``cochlea._assemble_and_solve_arc`` (mismos bloques, mismos signos,
    misma numeración de gdl) con ``D`` complejo cuando ``alpha > 0``:

        c̃ = c / (1 + i·α·c/ω)   ⇒   k = ω/c + i·α        (D4, ``e^{−iωt}``)
        λ̃ = ρ·(c̃_p² − 2·c̃_s²),  μ̃ = ρ·c̃_s²

    Devuelve ``(p, u)`` nodales complejos (``NaN`` donde el nodo no toca ese medio).
    """
    x = np.asarray(section.nodes, dtype=float)
    t = np.asarray(section.triangles, dtype=int)
    reg = np.asarray(section.regions)
    n = x.shape[0]
    fluid = section.fluid
    solid = section.bone
    omega = 2.0 * np.pi * section.frequency
    k_f = fluid.wavenumber(section.frequency)

    solid_e = reg == "solid"
    if not np.any(solid_e) or np.all(solid_e):
        raise ValueError("el corte necesita elementos 'solid' y 'fluid'")
    solid_tris, fluid_tris = t[solid_e], t[~solid_e]

    solid_node = np.zeros(n, dtype=bool)
    solid_node[solid_tris.ravel()] = True
    fluid_node = np.zeros(n, dtype=bool)
    fluid_node[fluid_tris.ravel()] = True
    reps_u = np.flatnonzero(solid_node)
    reps_p = np.flatnonzero(fluid_node)
    n_u = 2 * reps_u.size
    u_dof = -np.ones((n, 2), dtype=int)
    u_dof[reps_u, 0] = 2 * np.arange(reps_u.size)
    u_dof[reps_u, 1] = 2 * np.arange(reps_u.size) + 1
    p_dof = -np.ones(n, dtype=int)
    p_dof[reps_p] = n_u + np.arange(reps_p.size)
    size = n_u + reps_p.size

    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    vals: list[np.ndarray] = []

    def _add(rd, cd, loc) -> None:
        rd, cd = np.asarray(rd), np.asarray(cd)
        if rd.ndim == 1:
            rd, cd = rd[:, None], cd[:, None]
        nr, nc = rd.shape[1], cd.shape[1]
        rows.append(np.repeat(rd, nc, axis=1).ravel())
        cols.append(np.tile(cd, (1, nr)).ravel())
        vals.append(np.asarray(loc).reshape(-1))

    # --- sólido: rigidez con módulos complejos, masa real ----------------------
    _, m_s, area = p1_elastic_local(x, solid_tris, solid)
    if alpha > 0.0:
        cp = solid.c_p / (1.0 + 1j * alpha * solid.c_p / omega)
        cs = solid.c_s / (1.0 + 1j * alpha * solid.c_s / omega)
        lam_c = solid.rho * (cp**2 - 2.0 * cs**2)
        mu_c = solid.rho * cs**2
    else:
        lam_c, mu_c = complex(solid.lam), complex(solid.mu)
    d = np.array(
        [
            [lam_c + 2.0 * mu_c, lam_c, 0.0],
            [lam_c, lam_c + 2.0 * mu_c, 0.0],
            [0.0, 0.0, mu_c],
        ]
    )
    b = strain_matrix(triangle_gradients(x, solid_tris))
    stiff = area[:, None, None] * np.einsum("mki,kl,mlj->mij", b, d, b)
    _add(
        u_dof[solid_tris].reshape(-1, 6),
        u_dof[solid_tris].reshape(-1, 6),
        stiff - solid.rho * omega**2 * m_s,
    )

    # --- fluido (dentro y fuera es el mismo material) -------------------------
    m_f, kf, _ = p1_tri_local(x, fluid_tris)
    _add(p_dof[fluid_tris], p_dof[fluid_tris], (kf - k_f**2 * m_f) / fluid.rho)

    # --- acoplamiento p–u en r = a y r = a + t --------------------------------
    iface, iface_n = interface_edges(x, t, reg)
    if iface.shape[0] == 0:
        raise ValueError("el corte no tiene interfaz fluido/sólido")
    lengths = np.linalg.norm(x[iface[:, 1]] - x[iface[:, 0]], axis=1)
    ud = u_dof[iface].reshape(-1, 4)
    pd = p_dof[iface]
    if np.any(ud < 0) or np.any(pd < 0):
        raise RuntimeError("un nodo de interfaz no tiene los tres gdl")
    m_loc = weights_to_local(lengths)
    c_local = np.empty((iface.shape[0], 4, 2), dtype=float)
    nx_e = iface_n[:, 0][:, None, None]
    ny_e = iface_n[:, 1][:, None, None]
    c_local[:, 0, :] = (nx_e * m_loc)[:, 0, :]
    c_local[:, 1, :] = (ny_e * m_loc)[:, 0, :]
    c_local[:, 2, :] = (nx_e * m_loc)[:, 1, :]
    c_local[:, 3, :] = (ny_e * m_loc)[:, 1, :]
    _add(ud, pd, c_local)
    _add(pd, ud, omega**2 * np.transpose(c_local, (0, 2, 1)))

    # --- ABC exterior + término incidente ------------------------------------
    b_full = np.zeros(size, dtype=complex)
    inc = plane_wave(section.frequency, fluid, (1.0, 0.0), amplitude=1.0)
    abc = np.asarray(section.abc_edges, dtype=int).reshape(-1, 2)
    center = (0.0, 0.0) if abc_order == 2 else None
    blk = abc_boundary_block(
        x, abc, k_f, rho=fluid.rho, abc_order=abc_order, abc_center=center
    )
    coo = blk.tocoo()
    _add(p_dof[coo.row], p_dof[coo.col], coo.data)
    rhs = abc_incident_rhs(
        x, abc, inc, k_f, rho=fluid.rho, abc_order=abc_order, abc_center=center
    )
    has_p = p_dof >= 0
    np.add.at(b_full, p_dof[has_p], rhs[has_p])

    # --- ventana resistiva en el arco r = a ----------------------------------
    if section.has_window and section.window_model == "impedance":
        arc = window_arc_edges(section)
        if arc.shape[0] == 0:
            raise ValueError(
                f"el arco de ±{section.window_half_angle:g}° no tiene aristas de interfaz"
            )
        z = complex(section.window_impedance)
        k_arc = 2.0 * np.pi * section.frequency * fluid.rho / z
        arc_blk = abc_boundary_block(x, arc, k_arc, rho=fluid.rho, abc_order=1)
        coo = arc_blk.tocoo()
        _add(p_dof[coo.row], p_dof[coo.col], coo.data)

    a_full = sp.coo_matrix(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
        shape=(size, size),
        dtype=complex,
    ).tocsr()

    # MMD_AT_PLUS_A: 1.9 s contra 43.5 s del default COLAMD en este sistema (medido).
    full = spla.splu(a_full.tocsc(), permc_spec="MMD_AT_PLUS_A").solve(b_full)
    p = np.full(n, np.nan + 0.0j, dtype=complex)
    p[has_p] = full[p_dof[has_p]]
    u = np.full((n, 2), np.nan + 0.0j, dtype=complex)
    has_u = u_dof[:, 0] >= 0
    u[has_u, 0] = full[u_dof[has_u, 0]]
    u[has_u, 1] = full[u_dof[has_u, 1]]
    return p, u, size


_CACHE: dict[tuple, dict] = {}


def evaluate(
    scn: Scenario,
    *,
    ppw: float = PPW,
    abc_order: int = 1,
) -> dict:
    """Construye, resuelve y mide un escenario. Cacheado (el barrido repite puntos)."""
    key = (scn, float(ppw), int(abc_order))
    if key in _CACHE:
        return _CACHE[key]

    t0 = time.time()
    section = build_cross_section(
        FREQUENCY_HZ,
        RADIUS_M,
        scn.thickness_m,
        ppw=ppw,
        pad_wavelengths=PAD_WAVELENGTHS,
        min_elements_across_wall=MIN_ELEMENTS_ACROSS_WALL,
        bone=scn.bone(),
        window_half_angle=WINDOW_HALF_ANGLE_DEG,
        window_center_angle=WINDOW_CENTER_ANGLE_DEG,
        window_model="impedance",
        window_impedance=scn.impedance(),
    )
    p, _, size = solve_section(
        section, alpha=scn.alpha(), abc_order=abc_order
    )
    ring = section.ring_nodes(section.probe_radius)
    vals = p[ring]
    if not np.all(np.isfinite(vals)):
        raise RuntimeError("la circunferencia de la sonda tiene nodos sin gdl de presión")

    out = {
        "scenario": scn,
        "t_interior": float(np.mean(np.abs(vals))),
        "t_interior_complex": float(abs(np.mean(vals))),
        "t_interior_rms": float(np.sqrt(np.mean(np.abs(vals) ** 2))),
        "plate_T": float(abs(plate_reference(FREQUENCY_HZ, scn.thickness_m).T)),
        "gdl": int(size),
        "nodes": int(section.n_nodes),
        "alpha_np_m": float(scn.alpha()),
        "seconds": time.time() - t0,
        "ppw": float(ppw),
        "abc_order": int(abc_order),
    }
    _CACHE[key] = out
    return out


# --------------------------------------------------------------------------- #
# Puntos del barrido
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class SweepPoint:
    """Un punto etiquetado de un barrido OAT."""

    parametro: str
    valor: str
    scenario: Scenario
    notas: str = ""


def thickness_points(values) -> list[SweepPoint]:
    return [
        SweepPoint(
            "espesor",
            f"{t * 1e3:.1f} mm",
            Scenario(thickness_m=t),
            "cápsula 1–7 mm; el cráneo reportado es 6–8 mm (la cápsula ótica es más fina)",
        )
        for t in values
    ]


def oat_points(thickness_values) -> list[SweepPoint]:
    """Todos los barridos one-at-a-time, con los demás parámetros en el caso base."""
    pts: list[SweepPoint] = []

    # c_p — con c_s = 1500 fijo: 1900 no es admisible con c_p = 1752 (ν > 0.5)
    for cp in (1752.0, 2500.0, 3285.0, 3500.0):
        pts.append(
            SweepPoint(
                "c_p",
                f"{cp:.0f} m/s",
                Scenario(c_p=cp, c_s=1500.0),
                "Ammi 2008 (temporal, 1.03 MHz): 1752–3285. Barrido con c_s=1500 fijo "
                "porque c_s=1900 no es admisible con c_p=1752",
            )
        )

    for cs in (1360.0, 1500.0, 1640.0, 1900.0):
        pts.append(
            SweepPoint(
                "c_s",
                f"{cs:.0f} m/s",
                Scenario(c_s=cs),
                "White 2006 (sustituto de calvaria) 1500±140; 1900 es el base VERIFY (D5)",
            )
        )

    for rho in (1868.0, 1900.0, 1950.0):
        pts.append(
            SweepPoint(
                "rho",
                f"{rho:.0f} kg/m³",
                Scenario(rho=rho),
                "Peterson & Dechow 2003 (cortical del temporal)",
            )
        )

    pts.extend(thickness_points(thickness_values))

    for name in ("sin_atenuacion", "white_2006", "pichardo_2011", "pinton_2012"):
        scn = Scenario(law=name)
        a = scn.alpha()
        pts.append(
            SweepPoint(
                "atenuacion",
                name,
                scn,
                f"α(100 kHz) = {a:.2f} Np/m; leyes de D19 (efectivas de cráneo, no absorción "
                f"local); a 100 kHz son extrapolación (D20: no hay dato bajo 270 kHz)",
            )
        )

    for zf in (0.5, 1.0, 2.0):
        pts.append(
            SweepPoint(
                "Z_r",
                f"{zf:g} × ρ_f·c_f",
                Scenario(z_factor=zf),
                "sin impedancia citable de la ventana redonda (D20): se barre, no se fija",
            )
        )
    return pts


def grid_points(thickness_values, laws) -> list[SweepPoint]:
    """Grilla completa espesor × atenuación."""
    return [
        SweepPoint(
            "espesor_x_atenuacion",
            f"{t * 1e3:.1f} mm + {law}",
            Scenario(thickness_m=t, law=law),
            "grilla completa: los dos parámetros que se espera que dominen",
        )
        for t in thickness_values
        for law in laws
    ]


# --------------------------------------------------------------------------- #
# Análisis
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Row:
    parametro: str
    valor: str
    t_interior: float
    notas: str


@dataclass(frozen=True)
class Report:
    base: dict
    rows: tuple[Row, ...]
    grid: tuple[Row, ...]
    ranking: tuple[tuple[str, float, float, float], ...]  # (parámetro, mín, máx, rango rel.)

    def as_csv_rows(self) -> list[Row]:
        return list(self.rows) + list(self.grid)


def analyse(
    *,
    thickness_values=None,
    ppw: float = PPW,
    abc_order: int = 1,
    include_grid: bool = True,
    verbose: bool = True,
) -> Report:
    """Corre el OAT completo + la grilla espesor × atenuación y arma el ranking."""
    if thickness_values is None:
        thickness_values = [0.001 + 0.0005 * i for i in range(13)]  # 1.0 … 7.0 mm

    base = evaluate(Scenario(), ppw=ppw, abc_order=abc_order)
    t_base = base["t_interior"]

    if verbose:
        print("=" * 96)
        print("SENSIBILIDAD GLOBAL DE LA TRANSMISIÓN AL FLUIDO INTERIOR (cápsula coclear)")
        print("=" * 96)
        print(
            f"  a = {RADIUS_M * 1e3:.1f} mm   f = {FREQUENCY_HZ / 1e3:.0f} kHz   "
            f"ppw = {ppw:g}   ABC orden {abc_order}   pad = {PAD_WAVELENGTHS:g}λ   "
            f"ventana impedance ±{WINDOW_HALF_ANGLE_DEG:g}°   Z_r = ρ_f·c_f"
        )
        print(
            f"  base: c_p = 3500 m/s, c_s = 1900 m/s, ρ = 1950 kg/m³, t = 2.0 mm, "
            f"sin atenuación"
        )
        print(
            f"  → t_interior(base) = {t_base:.6f}   (referencia solve_cross_section "
            f"{T_INTERIOR_REFERENCE:.6f}, Δ = {100 * (t_base / T_INTERIOR_REFERENCE - 1):+.2f} %)"
            f"   gdl = {base['gdl']}   nodos = {base['nodes']}"
        )

    # ---- OAT ---------------------------------------------------------------
    points = oat_points(thickness_values)
    rows: list[Row] = []
    per_param: dict[str, list[float]] = {}
    if verbose:
        print("\n" + "-" * 96)
        print("ONE-AT-A-TIME (un parámetro por vez; el resto en el caso base)")
        print("-" * 96)
        print(f"  {'parámetro':>12}  {'valor':>18}  {'t_interior':>10}  {'Δ vs base':>10}   notas")

    for pt in points:
        res = evaluate(pt.scenario, ppw=ppw, abc_order=abc_order)
        v = res["t_interior"]
        per_param.setdefault(pt.parametro, []).append(v)
        note = pt.notas + f" | α = {res['alpha_np_m']:.2f} Np/m"
        rows.append(Row(pt.parametro, pt.valor, v, note))
        if verbose:
            print(
                f"  {pt.parametro:>12}  {pt.valor:>18}  {v:10.6f}  "
                f"{100 * (v / t_base - 1):+9.2f} %   {pt.notas[:56]}"
            )

    # ---- controles numéricos (no son parámetros físicos) --------------------
    controls: list[Row] = []
    fine = evaluate(Scenario(), ppw=2.0 * ppw, abc_order=abc_order)
    bgt2 = evaluate(Scenario(), ppw=ppw, abc_order=2)
    controls.append(
        Row(
            "control_numerico_ppw",
            f"ppw {2 * ppw:g} vs {ppw:g}",
            fine["t_interior"],
            f"piso de discretización: {100 * (fine['t_interior'] / t_base - 1):+.2f} % "
            "(misma física, malla más fina)",
        )
    )
    controls.append(
        Row(
            "control_numerico_abc_order",
            "BGT-2 (abc_order=2)",
            bgt2["t_interior"],
            f"piso de frontera (D21): {100 * (bgt2['t_interior'] / t_base - 1):+.2f} %",
        )
    )
    if abc_order == 2:
        controls = controls[1:]  # el control de ABC ya es el caso base
    rows.extend(controls)

    if verbose:
        print("\n  controles numéricos (no son parámetros físicos):")
        for r in controls:
            print(
                f"  {r.parametro:>26}  {r.valor:>22}  {r.t_interior:10.6f}  "
                f"{100 * (r.t_interior / t_base - 1):+9.2f} %   {r.notas}"
            )

    # ---- ranking -----------------------------------------------------------
    ranking = []
    for name in RANKED:
        vals = per_param.get(name)
        if not vals:
            continue
        lo, hi = min(vals), max(vals)
        ranking.append((name, lo, hi, (hi - lo) / t_base))
    ranking.sort(key=lambda x: x[3], reverse=True)

    if verbose:
        print("\n" + "-" * 96)
        print("RANKING DE INFLUENCIA  (rango relativo = (máx − mín) / t_interior base)")
        print("-" * 96)
        print(f"  {'#':>2}  {'parámetro':>12}  {'mín':>10}  {'máx':>10}  {'rango rel.':>11}  barra")
        for i, (name, lo, hi, rng) in enumerate(ranking, start=1):
            bar = "#" * max(1, int(round(60 * rng / ranking[0][3])))
            print(f"  {i:>2}  {name:>12}  {lo:10.6f}  {hi:10.6f}  {100 * rng:10.2f} %  {bar}")
        noise = abs(fine["t_interior"] / t_base - 1.0)
        print(
            f"\n  piso numérico medido (ppw {ppw:g} → {2 * ppw:g}) = {100 * noise:.2f} %: "
            "los parámetros por debajo de eso están dentro del error de malla."
        )

    # ---- grilla espesor × atenuación ---------------------------------------
    grid: list[Row] = []
    laws = tuple(LAWS)
    if verbose and include_grid:
        print("\n" + "-" * 96)
        print("GRILLA COMPLETA  espesor × atenuación   (t_interior)")
        print("-" * 96)
        print(f"  {'t [mm]':>7} " + " ".join(f"{law:>15}" for law in laws) + "   |T| placa")
    by_t: dict[float, dict[str, float]] = {}
    for t in thickness_values:
        by_t.setdefault(float(t), {})
    if include_grid:
        for pt in grid_points(thickness_values, laws):
            res = evaluate(pt.scenario, ppw=ppw, abc_order=abc_order)
            grid.append(Row(pt.parametro, pt.valor, res["t_interior"], pt.notas))
            by_t[pt.scenario.thickness_m][pt.scenario.law] = res["t_interior"]
    else:
        # sin grilla, la fila "sin atenuación" sale del OAT (ya resuelta y cacheada)
        by_t = {
            float(t): {
                "sin_atenuacion": evaluate(
                    Scenario(thickness_m=t), ppw=ppw, abc_order=abc_order
                )["t_interior"]
            }
            for t in thickness_values
        }
    for t in sorted(by_t):
        by_t[t]["plate_T"] = evaluate(
            Scenario(thickness_m=t), ppw=ppw, abc_order=abc_order
        )["plate_T"]
    if verbose and include_grid:
        for t_mm in sorted(by_t):
            row = by_t[t_mm]
            print(
                f"  {t_mm * 1e3:7.1f} "
                + " ".join(f"{row.get(law, float('nan')):15.5f}" for law in laws)
                + f"   {row['plate_T']:.4f}"
            )

    # ---- hallazgo: monotonía ----------------------------------------------
    t_sorted = sorted(by_t)
    t_vals = [by_t[t]["sin_atenuacion"] for t in t_sorted]
    mono = all(b < a for a, b in zip(t_vals, t_vals[1:]))
    breath = [
        evaluate(Scenario(thickness_m=t), ppw=ppw, abc_order=abc_order)[
            "t_interior_complex"
        ]
        for t in t_sorted
    ]
    plate = [by_t[t]["plate_T"] for t in t_sorted]
    breath_mono = all(b < a for a, b in zip(breath, breath[1:]))
    plate_mono = all(b < a for a, b in zip(plate, plate[1:]))
    if verbose:
        print("\n" + "-" * 96)
        print("HALLAZGO — ¿más hueso, menos pasa?")
        print("-" * 96)
        print(
            "  t_interior(t) sin atenuación: "
            + " ".join(f"{v:.4f}" for v in t_vals)
            + f"\n  monótona decreciente: {mono}"
        )
        print(
            "  |mean(p)| (armónico n=0, 'respiración'): "
            + " ".join(f"{v:.4f}" for v in breath)
            + f"\n  monótona decreciente: {breath_mono}"
        )
        print(
            "  |T| de la placa equivalente          : "
            + " ".join(f"{v:.4f}" for v in plate)
            + f"\n  monótona decreciente: {plate_mono}"
        )
        if not mono:
            print(
                "  ⇒ mean(|p|) en la cavidad cerrada es una amplitud de onda ESTACIONARIA (D18), "
                "no una\n    transmisión: oscila con el espesor. Lo que se comporta como "
                "transmisión es el armónico n=0."
            )
        if not breath_mono:
            steps = [
                (t_sorted[i + 1], breath[i + 1] - breath[i]) for i in range(len(breath) - 1)
            ]
            t_up, d_up = max(steps, key=lambda s: s[1])
            i_up = t_sorted.index(t_up)
            prev = breath[i_up - 1]
            print(
                f"    (el n=0 baja {breath[0] / breath[-1]:.1f}× en total; el mayor tramo que sube "
                f"es {prev:.4f} → {breath[i_up]:.4f} a t = {t_up * 1e3:g} mm, "
                f"{100 * d_up / prev:+.1f} % del tramo, dentro del piso numérico)"
            )

    return Report(base=base, rows=tuple(rows), grid=tuple(grid), ranking=tuple(ranking))


def write_csv(report: Report, path: Path = OUT_CSV) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["parametro", "valor", "t_interior", "notas"])
        for row in report.as_csv_rows():
            w.writerow([row.parametro, row.valor, f"{row.t_interior:.6f}", row.notas])
        w.writerow(
            [
                "ranking",
                "; ".join(f"{n}={100 * r:.2f}%" for n, _, _, r in report.ranking),
                f"{report.base['t_interior']:.6f}",
                "ranking por rango relativo (máx−mín)/base, mayor = más influyente",
            ]
        )
    return path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--step", type=float, default=0.5, help="paso del barrido de espesor [mm]")
    ap.add_argument("--ppw", type=float, default=PPW, help="puntos por longitud de onda")
    ap.add_argument("--abc-order", type=int, default=1, choices=(1, 2))
    ap.add_argument("--out", type=Path, default=OUT_CSV)
    ap.add_argument("--no-csv", action="store_true")
    args = ap.parse_args(argv)

    step = args.step * 1e-3
    n = int(round((7.0e-3 - 1.0e-3) / step)) + 1
    thickness = [1.0e-3 + step * i for i in range(n)]

    t0 = time.time()
    report = analyse(thickness_values=thickness, ppw=args.ppw, abc_order=args.abc_order)
    if not args.no_csv:
        path = write_csv(report, args.out)
        print(f"\n  CSV: {path}  ({len(report.as_csv_rows())} filas)")
    print(f"  tiempo total: {time.time() - t0:.1f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
