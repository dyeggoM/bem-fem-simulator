"""Helmholtz 3D P1 con ABC: malla, onda plana, convergencia y costo.

Qué verifica (y qué no)
-----------------------
1. **La malla**: volúmenes positivos, volumen total igual al de la caja, y caras de frontera con
   normal saliente. Si la orientación de los tetraedros está mal, todo lo demás es ruido.
2. **Onda plana con ABC**: el campo total tiene que ser ``p_i``. El criterio fuerte no es un número
   absoluto —P1 en 3D con una caja de 2 λ no baja del 3 %— sino el **orden de convergencia**: al
   refinar, el error tiene que caer como ``h²`` (medido 1.88 / 1.94 / 1.97).
3. **La ABC hace su trabajo**: con la pared rígida (sin ABC) el error es enorme porque la onda se
   refleja; con ABC baja a un dígito. Es la comparación que dice que la ABC está bien puesta.
4. **Costo**: tabla de nodos/tets/tiempo para dos resoluciones, que es lo que decide qué geometría
   es viable en 3D.

Sobre la fuente puntual: el camino existe (``_add_point_sources``) pero su amplitud todavía no está
validada contra el Green 3D — está anotado en ``docs/decisions.md`` (D32) como pendiente, y no se
testea acá para no fijar un número que no se entiende.
"""

from __future__ import annotations

import time

import numpy as np
import pytest

from ultrasound_sim import materials
from ultrasound_sim.fem import helmholtz3d as h3
from ultrasound_sim.meshing3d import (
    boundary_faces,
    box_tet_mesh,
    box_tet_mesh_from_spacing,
    faces_on_plane,
    tet_signed_volumes,
)

WATER = materials.WATER
FREQUENCY = 100e3
LAMBDA = WATER.wavelength(FREQUENCY)
K = WATER.wavenumber(FREQUENCY)


def _caja_1d(ppw: float, largo_lambdas: float = 2.0) -> tuple:
    """Caja ``largo × 1 × 1`` λ con ABC sólo en las dos caras de x.

    Las caras laterales son rígidas y eso es **consistente** con una onda plana que viaja en x
    (``n·d = 0``): la ABC se aplica donde la onda cruza, no en las paredes paralelas.
    """
    m = box_tet_mesh_from_spacing(
        (0.0, 0.0, 0.0), (largo_lambdas * LAMBDA, LAMBDA, LAMBDA), LAMBDA / ppw
    )
    faces = boundary_faces(m.tets)
    i0 = faces_on_plane(m.nodes, faces, 0, 0.0)
    i1 = faces_on_plane(m.nodes, faces, 0, largo_lambdas * LAMBDA)
    return m, np.vstack([faces[i0], faces[i1]])


# --------------------------------------------------------------------------- malla
def test_la_malla_de_caja_es_valida():
    m = box_tet_mesh(0.0, 2e-3, 0.0, 1e-3, 0.0, 3e-3, 4, 2, 6)
    vol = tet_signed_volumes(m.nodes, m.tets)
    print(f"\n  {m.nodes.shape[0]} nodos, {m.tets.shape[0]} tets, "
          f"volumen mínimo {vol.min():.3e} m³")
    assert m.tets.shape == (6 * 4 * 2 * 6, 4), "cada celda da 6 tetraedros"
    assert m.nodes.shape[0] == 5 * 3 * 7
    assert np.all(vol > 0.0), "todos los volúmenes tienen que ser positivos (orientación)"
    assert vol.sum() == pytest.approx(2e-3 * 1e-3 * 3e-3, rel=1e-12)

    faces = boundary_faces(m.tets)
    # celdas de la superficie: 2·(nx·ny + ny·nz + nx·nz) = 2·(8+12+24) = 88; y cada una aporta
    # 2 triángulos (cada cara del cubo se parte en dos)
    celdas_superficie = 2 * (4 * 2 + 2 * 6 + 4 * 6)
    assert faces.shape[0] == 2 * celdas_superficie


# --------------------------------------------------------------------------- onda plana
def test_onda_plana_con_abc_baja_el_error():
    """Con ABC el error es de discretización; sin ABC (pared rígida) hay onda estacionaria."""
    m, abc = _caja_1d(12.0)
    pw = h3.plane_wave(FREQUENCY, WATER, (1.0, 0.0, 0.0))
    exacto = pw.value(m.nodes)

    con_abc = h3.solve(m.nodes, m.tets, K, abc_faces=abc, p_inc=pw)
    sin_abc = h3.solve(m.nodes, m.tets, K)
    err_con = np.linalg.norm(con_abc - exacto) / np.linalg.norm(exacto)
    err_sin = np.linalg.norm(sin_abc - exacto) / np.linalg.norm(exacto)

    print(f"\n  error L2 con ABC = {100 * err_con:.2f} %   |   sin ABC (rígida) = {100 * err_sin:.1f} %")
    assert err_con < 0.15, f"la ABC no está funcionando: {100 * err_con:.1f} %"
    assert err_sin > 5.0 * err_con, "sin ABC tiene que reflejar mucho más que con ABC"


def test_convergencia_de_orden_dos():
    """P1 en 3D: el error tiene que caer como h² al refinar (se acepta 1.6–2.2)."""
    errores = {}
    for ppw in (8.0, 12.0, 16.0):
        m, abc = _caja_1d(ppw)
        pw = h3.plane_wave(FREQUENCY, WATER, (1.0, 0.0, 0.0))
        p = h3.solve(m.nodes, m.tets, K, abc_faces=abc, p_inc=pw)
        exacto = pw.value(m.nodes)
        errores[ppw] = np.linalg.norm(p - exacto) / np.linalg.norm(exacto)

    print("\n  ppw  error L2")
    for ppw, e in errores.items():
        print(f"  {ppw:4.0f}  {100 * e:6.2f} %")
    ordenes = [
        np.log(errores[a] / errores[b]) / np.log(b / a)
        for a, b in ((8.0, 12.0), (12.0, 16.0))
    ]
    print("  orden medido: " + " / ".join(f"{o:.2f}" for o in ordenes))
    assert all(1.6 < o < 2.2 for o in ordenes), f"orden fuera de rango: {ordenes}"
    assert errores[16.0] < 0.08, f"error a ppw=16 demasiado alto: {100 * errores[16.0]:.1f} %"


def test_el_pico_apunta_hacia_donde_va_la_onda():
    """La fase tiene que avanzar en +x: ``arg p`` crece con ``x`` a razón de ``k``."""
    m, abc = _caja_1d(12.0)
    pw = h3.plane_wave(FREQUENCY, WATER, (1.0, 0.0, 0.0))
    p = h3.solve(m.nodes, m.tets, K, abc_faces=abc, p_inc=pw)

    eje = np.isclose(m.nodes[:, 1], 0.0) & np.isclose(m.nodes[:, 2], 0.0) & (m.nodes[:, 0] > 0.3 * LAMBDA)
    x = m.nodes[eje, 0]
    orden = np.argsort(x)
    pendiente = np.polyfit(x[orden], np.unwrap(np.angle(p[eje][orden])), 1)[0]
    print(f"\n  d(arg p)/dx = {pendiente:.3e} 1/m   contra k = {K.real:.3e} 1/m")
    assert pendiente == pytest.approx(K.real, rel=0.05)


def test_fuente_puntual_reproduce_el_green_3d():
    """Un monopolo de caudal ``Q`` tiene que dar ``p = −i·ω·ρ·Q·e^{ikr}/(4πr)``.

    El signo importa: con la convención del repo (``e^{−iωt}``, D4) y un RHS ``−iωQ`` lumpado en el
    nodo, la solución es ``p = S·G`` con ``S = −iωρQ`` y ``G = e^{ikr}/(4πr)``. Equivocarlo da un
    desfase constante de 180° sin tocar la amplitud — que es exactamente lo que pasó la primera vez
    que se miró esto a mano.
    """
    f = FREQUENCY
    L = 1.1 * LAMBDA
    m = box_tet_mesh_from_spacing((-L / 2, -L / 2, -L / 2), (L / 2, L / 2, L / 2), LAMBDA / 8)
    faces = boundary_faces(m.tets)
    idx = int(np.argmin(np.linalg.norm(m.nodes, axis=1)))
    fuente = m.nodes[idx]

    p = h3.solve(m.nodes, m.tets, K, abc_faces=faces, element_rho=WATER.rho,
                 point_sources=fuente[None, :], point_volume_velocities=np.array([1.0]),
                 frequency=f)

    r = np.linalg.norm(m.nodes - fuente, axis=1)
    anillo = (r > 0.25 * LAMBDA) & (r < 0.45 * LAMBDA)
    omega = 2.0 * np.pi * f
    # sólo en el anillo: en el nodo de la fuente r = 0 y el Green diverge (ahí no se evalúa)
    r_seguro = np.where(r > 0.0, r, 1.0)
    exacto = -1j * omega * WATER.rho * np.exp(1j * K * r_seguro) / (4.0 * np.pi * r_seguro)
    err = np.linalg.norm(p[anillo] - exacto[anillo]) / np.linalg.norm(exacto[anillo])
    amplitud = float(np.mean(np.abs(p[anillo]) * r[anillo]) / (omega * WATER.rho / (4 * np.pi)))

    print(f"\n  anillo {0.25:.2f}–{0.45:.2f} λ ({int(anillo.sum())} nodos): "
          f"error L2 {100 * err:.2f} %   |p|·r / teórico = {amplitud:.4f}")
    # La caja es chica a propósito (1.1 λ, la ABC queda a 0.55 λ de la fuente) y ppw=8: el 8.4 %
    # medido es eso, no el Green. Lo que se exige con fuerza es la AMPLITUD (el factor ρ y el −i).
    assert err < 0.10, f"no reproduce el Green 3D: {100 * err:.1f} %"
    assert amplitud == pytest.approx(1.0, abs=0.02)


def test_estratificado_rt_contra_la_matriz_de_transferencia():
    """Una losa en un medio homogéneo: ``R`` y ``|T|`` contra la matriz de transferencia 1D.

    Por qué no se compara el campo punto a punto: la **dispersión de P1 se acumula** con la
    distancia recorrida (medido: el error crece 13 % → 27 % → 39 % de izquierda a derecha, y es
    *idéntico* en el caso homogéneo). Eso es el mismo techo que D10 documentó en 2D, no un problema
    del estratificado. Por eso se ajustan las dos ondas en una ventana **pegada a cada interfaz**
    (0.5 λ a cada lado) y se comparan los coeficientes.

    ``R`` es la relación entre dos ondas del *mismo* medio, así que no depende de ninguna referencia
    de fase; ``|T|`` tampoco. La fase de ``T`` sí depende de si se refiere a ``x = 0`` o a la salida
    de la losa, y por eso no se exige.
    """
    from ultrasound_sim.analytic import helmholtz1d

    # material de prueba declarado (no es un valor de la literatura): contraste fuerte para que R
    # sea medible con ppw razonable
    fluido = materials.Material(name="fluido de prueba", c=2000.0, rho=1200.0,
                                provenance="TEST: sintético, 2000 m/s y 1200 kg/m³")
    f = FREQUENCY
    lam = LAMBDA
    ppw, pad, t = 12.0, 1.5 * lam, 0.5 * lam
    L = pad + t + pad
    m = box_tet_mesh(-pad, t + pad, -0.4 * lam, 0.4 * lam, -0.4 * lam, 0.4 * lam,
                     int(round(L / (lam / ppw))), int(round(0.8 * lam / (lam / ppw))),
                     int(round(0.8 * lam / (lam / ppw))))
    cen = m.nodes[m.tets].mean(axis=1)[:, 0]
    dentro = (cen > 0.0) & (cen < t)
    k = np.where(dentro, fluido.wavenumber(f), K)
    rho = np.where(dentro, fluido.rho, WATER.rho)

    faces = boundary_faces(m.tets)
    i0 = faces_on_plane(m.nodes, faces, 0, -pad)
    i1 = faces_on_plane(m.nodes, faces, 0, t + pad)
    p = h3.solve(m.nodes, m.tets, k, abc_faces=np.vstack([faces[i0], faces[i1]]), element_rho=rho,
                 p_inc=h3.plane_wave(f, WATER, (1.0, 0.0, 0.0)))

    x = m.nodes[:, 0]

    def ajustar(sel, kk, dos_ondas):
        xs = x[sel]
        a = np.column_stack([np.exp(1j * kk * xs), np.exp(-1j * kk * xs)]) if dos_ondas \
            else np.exp(1j * kk * xs)[:, None]
        coef, *_ = np.linalg.lstsq(a, p[sel], rcond=None)
        return coef

    izq = ajustar((x > -0.5 * lam) & (x < -1e-9), K, True)
    der = ajustar((x > t + 1e-9) & (x < t + 0.5 * lam), fluido.wavenumber(f), False)
    r_num = izq[1] / izq[0]
    t_num = der[0] / izq[0]

    analitico = helmholtz1d.stack_response(f, (WATER, fluido, WATER), (t,))
    e_r = abs(r_num - analitico.R) / abs(analitico.R)
    e_t = abs(abs(t_num) - abs(analitico.T)) / abs(analitico.T)

    print(f"\n  R: numérico {r_num:.5f}  analítico {analitico.R:.5f}  (error {100 * e_r:.2f} %)")
    print(f"  |T|: numérico {abs(t_num):.5f}  analítico {abs(analitico.T):.5f}  "
          f"(error {100 * e_t:.2f} %)")
    assert e_r < 0.05, f"R fuera de 5 %: {100 * e_r:.1f} %"
    assert e_t < 0.05, f"|T| fuera de 5 %: {100 * e_t:.1f} %"


def test_campo_dispersado_nulo_en_un_medio_homogeneo():
    """Sin dispersor, ``p_s = 0`` **exacto**: la onda incidente ya satisface el problema dispersado.

    Es el sanity check del ensamblado dispersado: si el lado derecho tuviera algo de más (un término
    incidente mal puesto, un signo en la ABC), acá aparecería un campo espurio en vez de cero.

    Lo que este test **no** dice: que el camino sirva para un dispersor penetrable. `assemble_scattered`
    no arma el término de contraste voluminoso `−(k² − k_i²)·p_i/ρ`, así que con un hueco blando
    también da ``p_s = 0`` (medido). Para un obstáculo **rígido** sí está el camino: `neumann_faces`
    con `dp_dn_incident = −∂p_i/∂n`. Queda como pendiente en D32.
    """
    pad = 2.0 * LAMBDA
    m = box_tet_mesh_from_spacing((-pad, 0.0, 0.0), (pad, LAMBDA, LAMBDA), LAMBDA / 12)
    faces = boundary_faces(m.tets)
    abc = np.vstack([faces[faces_on_plane(m.nodes, faces, 0, -pad)],
                     faces[faces_on_plane(m.nodes, faces, 0, pad)]])
    p_s = h3.solve_scattered(m.nodes, m.tets, K, abc_faces=abc, element_rho=WATER.rho)

    print(f"\n  |p_s| máximo / |p_i| medio = {np.abs(p_s).max() / np.abs(K).real:.3e}  "
          f"(|p_s| máximo crudo = {np.abs(p_s).max():.3e})")
    assert np.abs(p_s).max() < 1e-9, f"apareció campo dispersado espurio: {np.abs(p_s).max():.3e}"


# --------------------------------------------------------------------------- costo
@pytest.mark.slow
def test_presupuesto_de_costo():
    """La tabla que decide qué es viable: nodos, tets y segundos por resolución."""
    print("\n  ppw   nodos    tets   tiempo")
    for ppw in (8.0, 12.0):
        m, abc = _caja_1d(ppw)
        t0 = time.perf_counter()
        h3.solve(m.nodes, m.tets, K, abc_faces=abc,
                 p_inc=h3.plane_wave(FREQUENCY, WATER, (1.0, 0.0, 0.0)))
        dt = time.perf_counter() - t0
        print(f"  {ppw:4.0f} {m.nodes.shape[0]:7d} {m.tets.shape[0]:7d}  {dt:7.2f} s")
        assert dt < 120.0
