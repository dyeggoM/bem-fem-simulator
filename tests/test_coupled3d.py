"""Fase 6b: malla de esfera sumergida, elasticidad 3D y acoplamiento fluido–sólido.

Qué está validado acá y qué no
------------------------------
**Validado** (tests de esta archivo):

- La malla: volumen que converge al de la esfera al refinar y **interfaz conforme** (cada cara
  interior pertenece a exactamente dos tetraedros).
- El ensamblado elástico: los seis modos de cuerpo rígido caen en el núcleo de ``K`` a precisión de
  máquina y la energía de una deformación uniforme da exacta.
- El acoplamiento, por sus **límites**: sin contraste la dispersión tiende a cero, y con un sólido
  muy rígido el resultado coincide con la esfera rígida de `sphere.py` dentro del error del dominio.

**No validado todavía**: la comparación con **Faran** para un sólido real (acero) da **22 %** a
``ka = 1``. El patrón angular coincide en forma y fase, pero la amplitud sale 18 % baja, y el
desplazamiento interior del FEM queda casi uniforme (rígido) mientras Faran predice deformación. Los
tests de acoplamiento de acá usan tolerancias holgadas y **medidas**, no una meta.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from ultrasound_sim import materials as mats
from ultrasound_sim.analytic import faran
from ultrasound_sim.fem import coupled3d, elastic3d
from ultrasound_sim.meshing3d import box_tet_mesh
from ultrasound_sim.meshing3d_sphere import sphere_shell_tet_mesh

F = 100.0e3
AGUA = mats.WATER
ACERO = mats.ElasticSolid(
    name="acero (benchmark)", rho=7800.0, c_p=5900.0, c_s=3200.0, provenance="benchmark (D5)"
)
RIGIDO = mats.ElasticSolid(
    name="casi rigido", rho=1.0e6, c_p=1.0e5, c_s=6.0e4, provenance="límite numérico"
)


def _radio() -> float:
    return 1.0 / abs(complex(AGUA.wavenumber(F)))


# --------------------------------------------------------------------------- malla


def test_el_volumen_de_las_esferas_converge():
    a = _radio()
    errores = []
    for div in (3.0, 4.5):
        malla = sphere_shell_tet_mesh(a, 4.0 * a, a / div)
        exacto = 4.0 / 3.0 * np.pi * a**3
        errores.append(abs(malla.solid_volume() - exacto) / exacto)
    print(f"\n  error de volumen: {[f'{100*e:.2f} %' for e in errores]}")
    assert errores[0] < 0.05
    assert errores[1] < errores[0], "refinar no mejora el volumen"


def test_la_interfaz_es_conforme_y_esta_sobre_la_esfera():
    """Cada cara de interfaz tiene que ser frontera del sólido **y** del fluido, y estar en r = a."""
    a = _radio()
    malla = sphere_shell_tet_mesh(a, 4.0 * a, a / 3.5)

    def caras_una_vez(tets):
        from ultrasound_sim.meshing3d_sphere import _faces_once, _tet_faces

        return {tuple(sorted(f)) for f in _faces_once(_tet_faces(tets))[0]}

    por_solido = caras_una_vez(malla.tets_solid)
    por_fluido = caras_una_vez(malla.tets_fluid)
    interfaz = {tuple(sorted(f)) for f in malla.interface_faces}

    assert interfaz <= por_solido, "hay caras de interfaz que no son frontera del sólido"
    assert interfaz <= por_fluido, "hay caras de interfaz que no son frontera del fluido"

    r = np.linalg.norm(malla.nodes[malla.interface_faces].mean(axis=1), axis=1)
    print(f"\n  {len(interfaz)} caras · radio {r.min()/a:.4f}–{r.max()/a:.4f} a "
          f"· área {malla.interface_areas().sum()/(4*np.pi*a**2):.4f}")
    # 0.96 y no 0.97: los centroides de una superficie facetada caen por dentro. Con h = a/3.5 el
    # mínimo medido es 0.9617a; el área, en cambio, converge al 1.6 %.
    assert r.min() > 0.95 * a and r.max() <= a * (1 + 1e-9)
    assert abs(malla.interface_areas().sum() / (4 * np.pi * a**2) - 1.0) < 0.03


def test_la_normal_de_la_interfaz_apunta_hacia_afuera_del_solido():
    a = _radio()
    malla = sphere_shell_tet_mesh(a, 4.0 * a, a / 3.5)
    n = malla.interface_normals()
    centros = malla.nodes[malla.interface_faces].mean(axis=1)
    radial = centros / np.linalg.norm(centros, axis=1)[:, None]
    cosenos = np.einsum("ij,ij->i", n, radial)
    print(f"\n  n·r̂ mínimo {cosenos.min():.4f} (positivo = sale del sólido)")
    # No se exige que sea ~1: en una superficie facetada la normal de la cara se aparta de la radial
    # tanto como el ángulo del facetado (medido: 0.86 con h = a/3.5). Lo que importa es el SIGNO:
    # una cara con la normal invertida daría ~−1 y rompería el signo de todo el acoplamiento.
    assert cosenos.min() > 0.5, "alguna normal de interfaz no sale del sólido"


# --------------------------------------------------------------------------- elasticidad


def test_los_modos_de_cuerpo_rigido_estan_en_el_nucleo():
    malla = box_tet_mesh(0, 0.010, 0, 0.012, 0, 0.009, 4, 5, 3)
    k, _ = elastic3d.p1_elastic3d_matrices(malla.nodes, malla.tets, ACERO)
    modos = elastic3d.rigid_body_modes(malla.nodes)
    peor = max(
        np.abs(k @ modos[:, j]).max() / (np.abs(k).max() * np.abs(modos[:, j]).max())
        for j in range(6)
    )
    print(f"\n  peor residuo relativo en los 6 modos rígidos: {peor:.2e}")
    assert peor < 1e-12


def test_la_energia_de_una_deformacion_uniforme_es_exacta():
    """P1 reproduce campos lineales sin error: el test no tiene tolerancia de malla."""
    malla = box_tet_mesh(0, 0.010, 0, 0.012, 0, 0.009, 5, 6, 4)
    k, _ = elastic3d.p1_elastic3d_matrices(malla.nodes, malla.tets, ACERO)
    rng = np.random.default_rng(7)
    a = rng.normal(size=(3, 3))
    e = 0.5 * (a + a.T)
    u = (malla.nodes @ e).reshape(-1)
    energia_fe = 0.5 * u @ (k @ u)

    d = elastic3d.elasticity_matrix(ACERO)
    eps = np.array([e[0, 0], e[1, 1], e[2, 2], 2 * e[1, 2], 2 * e[0, 2], 2 * e[0, 1]])
    energia_ex = 0.5 * malla.volumes().sum() * eps @ d @ eps
    print(f"\n  FE {energia_fe:.8e} · exacta {energia_ex:.8e}")
    assert abs(energia_fe - energia_ex) / abs(energia_ex) < 1e-12


def test_la_masa_total_es_rho_por_volumen():
    malla = box_tet_mesh(0, 0.010, 0, 0.012, 0, 0.009, 4, 5, 3)
    _, m = elastic3d.p1_elastic3d_matrices(malla.nodes, malla.tets, ACERO)
    uno = np.ones(3 * malla.n_nodes)
    total = uno @ (m @ uno)
    exacto = 3.0 * ACERO.rho * malla.volumes().sum()
    print(f"\n  masa {total:.6e} · 3ρV {exacto:.6e}")
    assert abs(total - exacto) / exacto < 1e-12


def test_el_lame_reconstruye_las_velocidades():
    lam, mu = elastic3d.lame_parameters(ACERO)
    assert np.isclose(np.sqrt((lam + 2 * mu) / ACERO.rho), ACERO.c_p)
    assert np.isclose(np.sqrt(mu / ACERO.rho), ACERO.c_s)


# --------------------------------------------------------------------------- acoplamiento


@pytest.fixture(scope="module")
def malla_esfera():
    """Malla **rápida** para los tests de límites: dominio de 2.5a, 1.3 s por solve.

    No es la malla de producción: a 2.5a el ABC-2 deja ~9 % de error. La comparación con Faran en el
    dominio recomendado (4a) va marcada `slow`, porque son 32 s.
    """
    a = _radio()
    with pytest.warns(UserWarning, match="ABC-2 domina"):
        return sphere_shell_tet_mesh(a, 2.5 * a, a / 3.0)


def _puntos(a, n=13):
    th = np.linspace(0.0, np.pi, n)
    return 1.5 * a * np.stack([np.sin(th), np.zeros_like(th), np.cos(th)], axis=1)


def test_sin_contraste_no_dispersa(malla_esfera):
    """Una esfera indistinguible del fluido no puede dispersar: si aparece dispersión, el término
    incidente del acoplamiento está mal (es el caso que delató el `k₀` de más en Faran, D45)."""
    casi_agua = mats.ElasticSolid(
        name="casi agua", rho=AGUA.rho, c_p=float(AGUA.c), c_s=30.0, provenance="límite"
    )
    # con el cuerpo transparente hay que APAGAR la restricción de rígidos (ver el otro test)
    r = coupled3d.resolver(malla_esfera, AGUA, casi_agua, F, abc_order=2, fijar_rigidos=False)
    dispersion = np.abs(r["p_s"]).max()
    print(f"\n  |p_s| máximo con contraste nulo: {dispersion:.3e} (incidente = 1 Pa)")
    # 3 % de la incidente en esta malla (h = a/3): es el error de discretización del contraste nulo,
    # no dispersión física. Con h = a/5 baja a 0.3 %, así que la cota está holgada a propósito para
    # no fijar un umbral que dependa del mallador.
    assert dispersion < 5e-2


def test_el_limite_rigido_coincide_con_la_esfera_rigida(malla_esfera):
    """Con un sólido muy rígido, el acoplado tiene que dar lo mismo que la esfera rígida de
    `sphere.py` — con la tolerancia del **dominio**, no una ideal: el mismo dominio resuelto con el
    Helmholtz ya validado (D41) da 10.35 % en esta malla, así que exigir 2 % sería exigirle al
    acoplado más que al fluido solo."""
    from ultrasound_sim.analytic import sphere

    a = _radio()
    r = coupled3d.resolver(malla_esfera, AGUA, RIGIDO, F, abc_order=2)
    pts = _puntos(a)
    fem = coupled3d._interpolar(malla_esfera, r["p_s"], pts)
    rig = sphere.fluid_sphere_scattering(F, AGUA, AGUA, a, pts, boundary="rigid")["scattered"]
    util = np.abs(rig) > 1e-3 * np.abs(rig).max()
    error = np.abs(fem[util] - rig[util]).max() / np.abs(rig[util]).max()
    print(f"\n  error contra esfera rígida: {100*error:.2f} %")
    assert error < 0.20


def test_el_desplazamiento_del_solido_es_del_orden_del_de_una_presion_de_1_pa(malla_esfera):
    """Escala física: con 1 Pa de incidente y ``ρcω`` del orden de 2.4e6, ``u ~ p/(ρcω)``. Un
    desplazamiento 10³ veces más grande significaría una rigidez que no entró."""
    r = coupled3d.resolver(malla_esfera, AGUA, ACERO, F, abc_order=2)
    escala = 1.0 / (AGUA.rho * float(AGUA.c) * 2 * np.pi * F)
    print(f"\n  |u|max {np.abs(r['u']).max():.3e} m · escala p/(ρcω) {escala:.3e} m")
    assert 1e-3 * escala < np.abs(r["u"]).max() < 1e3 * escala

# --------------------------------------------------------------------------- contra Faran


def test_contra_faran_con_los_rigidos_restringidos(malla_esfera):
    """La comparación con Faran con los modos rígidos restringidos.

    Por qué hace falta la restricción: el sólido libre tiene seis modos de cuerpo rígido **sin
    rigidez elástica**. En una onda plana la fuerza neta exacta sobre un cuerpo real es cero, y la
    ecuación de movimiento da entonces desplazamiento medio nulo; pero la malla no cancela esa fuerza
    exactamente y el residuo excita la traslación sin freno. Medido: el desplazamiento normal en la
    interfaz salía 20× el exacto y la amplitud dispersada 16 % baja.

    Con la restricción el error cae de 18.15 % a **9.70 %**, que es el piso del dominio: el mismo
    malla resuelta sólo en fluido con el desplazamiento *exacto* de Faran en la interfaz da 9.73 %.
    Es decir, el sólido deja de aportar error y lo que queda es malla + ABC.
    """
    a = _radio()
    coef = faran.coeficientes(F, AGUA, ACERO, a)
    pts = _puntos(a)
    exacto = faran.campo_dispersado(coef, pts, a)
    util = np.abs(exacto) > 1e-3 * np.abs(exacto).max()

    r = coupled3d.resolver(malla_esfera, AGUA, ACERO, F, abc_order=2, fijar_rigidos=True)
    fem = coupled3d._interpolar(malla_esfera, r["p_s"], pts)
    error = np.abs(fem[util] - exacto[util]).max() / np.abs(exacto[util]).max()
    razon = np.abs(fem).max() / np.abs(exacto).max()
    print(f"\n  error contra Faran: {100*error:.2f} % · razón de amplitudes {razon:.3f}")
    assert error < 0.13
    assert abs(razon - 1.0) < 0.08


def test_la_fuerza_neta_que_la_malla_no_cancela_es_chica(malla_esfera):
    """El multiplicador de la restricción **es** la fuerza y el torque netos que la malla no supo
    cancelar. Se devuelve para poder medirlo en vez de que quede escondido en un desplazamiento
    espurio."""
    r = coupled3d.resolver(malla_esfera, AGUA, ACERO, F, abc_order=2, fijar_rigidos=True)
    residuo = r["residuo_rigido"]
    assert residuo is not None and residuo.size == 6
    print(f"\n  |residuo rígido| = {np.abs(residuo).max():.2e}")
    assert np.abs(residuo).max() < 1e-5


def test_sin_la_restriccion_el_error_es_peor_y_se_ve(malla_esfera):
    """Sin restricción el error es ~2× peor. Se deja medido: si algún cambio de formulación lo
    arregla, este test avisa (y habría que quitar la restricción, que no es gratis)."""
    a = _radio()
    coef = faran.coeficientes(F, AGUA, ACERO, a)
    pts = _puntos(a)
    exacto = faran.campo_dispersado(coef, pts, a)
    util = np.abs(exacto) > 1e-3 * np.abs(exacto).max()
    r = coupled3d.resolver(malla_esfera, AGUA, ACERO, F, abc_order=2, fijar_rigidos=False)
    fem = coupled3d._interpolar(malla_esfera, r["p_s"], pts)
    error = np.abs(fem[util] - exacto[util]).max() / np.abs(exacto[util]).max()
    print(f"\n  error SIN restricción: {100*error:.2f} % (con restricción: ~9.7 %)")
    assert error < 0.25


def test_la_restriccion_no_se_usa_cuando_el_cuerpo_es_transparente(malla_esfera):
    """Con la esfera indistinguible del medio, el cuerpo **sí** tiene desplazamiento medio distinto de
    cero (la onda incidente pasa de largo y ejerce fuerza neta), así que restringir los rígidos mete
    una dispersión falsa de 1.6 Pa sobre una incidente de 1 Pa. Es la contraprueba de que la
    restricción es una condición física con dominio de validez, no un truco que siempre ayuda."""
    casi_agua = mats.ElasticSolid(
        name="casi agua", rho=AGUA.rho, c_p=float(AGUA.c), c_s=30.0, provenance="límite"
    )
    sin = coupled3d.resolver(malla_esfera, AGUA, casi_agua, F, abc_order=2, fijar_rigidos=False)
    con = coupled3d.resolver(malla_esfera, AGUA, casi_agua, F, abc_order=2, fijar_rigidos=True)
    print(f"\n  |p_s| sin restricción {np.abs(sin['p_s']).max():.3e} · con restricción "
          f"{np.abs(con['p_s']).max():.3e}")
    assert np.abs(sin["p_s"]).max() < 5e-2
    assert np.abs(con["p_s"]).max() > 0.5

# --------------------------------------------------------------------------- dominio y piso


def test_el_dominio_chico_avisa():
    """El ABC-2 necesita dominio: a 2.5a deja ~9 % de error. Se avisa en vez de devolver un número
    malo sin decir nada."""
    a = _radio()
    with pytest.warns(UserWarning, match="ABC-2 domina"):
        sphere_shell_tet_mesh(a, 2.5 * a, a / 4.0)
    # el dominio recomendado no avisa
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        sphere_shell_tet_mesh(a, 4.0 * a, a / 4.0)


@pytest.mark.slow
def test_contra_faran_en_el_dominio_de_produccion():
    """El número que se reporta: dominio de 4a y h = a/4.5.

    Medido: **4.21 %** (con el fluido solo y datos exactos, 4.06 % — o sea que el sólido aporta 0.15
    puntos y el resto es malla+ABC). Costo: 32 s y 21 398 nodos, por eso la marca `slow`.
    """
    a = _radio()
    mall = sphere_shell_tet_mesh(a, 4.0 * a, a / 4.5)
    coef = faran.coeficientes(F, AGUA, ACERO, a)
    pts = _puntos(a)
    exacto = faran.campo_dispersado(coef, pts, a)
    util = np.abs(exacto) > 1e-3 * np.abs(exacto).max()

    r = coupled3d.resolver(mall, AGUA, ACERO, F, abc_order=2)
    fem = coupled3d._interpolar(mall, r["p_s"], pts)
    error = np.abs(fem[util] - exacto[util]).max() / np.abs(exacto[util]).max()
    print(f"\n  4a: error {100*error:.2f} % · {mall.n_nodes} nodos (vs 10.20 % a 2.5a)")
    assert error < 0.06
