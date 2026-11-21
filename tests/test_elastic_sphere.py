"""Tests de la vista "esfera elástica" (U7): el cliente de `coupled3d` + `faran` en la UI.

Lo que se verifica acá es el **contrato de la vista**, no la física: la física está validada en
`test_coupled3d.py` y `test_faran.py`. Concretamente:

- que el costo estimado sea coherente con el modo elegido y avise cuando el dominio es chico;
- que el modo "rápida" **no** se declare validado (su error es del ABC, y la UI tiene que decirlo);
- que el sello del modo "precisa" salga del registro de criterios y no de una constante;
- que el mapa venga enmascarado dentro de la esfera y no con valores inventados.
"""

from __future__ import annotations

import numpy as np
import pytest

from ultrasound_sim.dashboard import elastic_sphere, validacion

F = 100.0e3


def test_los_dos_modos_declaran_su_costo_y_su_aviso():
    rapida = elastic_sphere.costo({"modo": "rapida", "ka": 1.0})
    precisa = elastic_sphere.costo({"modo": "precisa", "ka": 1.0})
    print(f"\n  rápida: {rapida['nodos_estimados']} nodos · {rapida['segundos_estimados']} s")
    print(f"  precisa: {precisa['nodos_estimados']} nodos · {precisa['segundos_estimados']} s")
    assert rapida["nodos_estimados"] < precisa["nodos_estimados"]
    assert rapida["aviso"] and "ABC" in rapida["aviso"]
    assert precisa["aviso"] is None


def test_el_radio_sale_del_ka_y_la_frecuencia():
    """`ka = k₀·a`: el radio no es un parámetro libre, es `ka/k₀`. Si la UI lo tratara como dato
    independiente, dos configuraciones con el mismo `ka` darían esferas distintas."""
    agua_k0 = abs(complex(__import__("ultrasound_sim.materials", fromlist=["WATER"]).WATER
                          .wavenumber(F)))
    costo = elastic_sphere.costo({"modo": "rapida", "ka": 2.0})
    assert np.isclose(costo["radio_m"], 2.0 / agua_k0)


def test_parametros_fuera_de_rango_se_rechazan():
    for datos in ({"ka": 0.1}, {"ka": 5.0}, {"frecuencia_hz": 5e3}, {"modo": "otro"}):
        with pytest.raises(ValueError):
            elastic_sphere.costo(datos)
    with pytest.raises(ValueError):
        elastic_sphere.resolver({"material": "plastilina"})


@pytest.fixture(scope="module")
def corrida_rapida():
    """Una sola corrida del modo rápido (~4 s) para los tests que miran la salida."""
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return elastic_sphere.resolver(
            {"modo": "rapida", "ka": 1.0, "material": "acero", "puntos_angulo": 61,
             "fijar_rigidos": True}
        )


def test_el_modo_rapido_no_se_declara_validado(corrida_rapida):
    """Su error (~11 %) es del ABC, no del solver: el sello tiene que decir FUERA y la UI agregar
    por qué. Si algún día el modo rápido saliera validado, sería porque el criterio dejó de medir."""
    sello = corrida_rapida["sello"]
    error = corrida_rapida["metrica"]["error_relativo"]
    print(f"\n  error {100*error:.2f} % · sello {sello['estado']}")
    assert sello["estado"] == "fuera"
    assert sello["criterios"][0]["id"] == "volumen_esfera_elastica"
    assert error > sello["criterios"][0]["tolerancia"]


def test_el_sello_sale_del_registro_y_no_de_una_constante(corrida_rapida):
    criterio = validacion.CRITERIOS["volumen_esfera_elastica"]
    assert criterio.referencia == 0.0
    assert "Faran" in criterio.fuente and "D46" in criterio.fuente
    # el mismo error, sellado a mano contra el registro, tiene que dar el mismo estado
    medido = corrida_rapida["metrica"]["error_relativo"]
    assert validacion.sello("volumen_esfera_elastica", medido)["estado"] == \
        corrida_rapida["sello"]["estado"]


def test_la_curva_angular_viene_comparada_punto_a_punto(corrida_rapida):
    a = corrida_rapida["angulos"]
    n = len(a["theta_deg"])
    assert n == len(a["fem_modulo"]) == len(a["faran_modulo"]) == 61
    assert a["theta_deg"][0] == pytest.approx(0.0) and a["theta_deg"][-1] == pytest.approx(180.0)
    razon = corrida_rapida["metrica"]["razon_amplitudes"]
    print(f"\n  razón de amplitudes {razon:.3f} (1.0 = misma amplitud que Faran)")
    assert 0.8 < razon < 1.2


def test_el_mapa_enmascara_la_esfera_en_vez_de_rellenarla(corrida_rapida):
    """Dentro de la esfera no hay fluido: si aparecieran valores ahí, la interpolación estaría
    inventando un campo donde el problema no lo define."""
    mapa = corrida_rapida["mapa"]
    z = np.array([[np.nan if v is None else v for v in fila] for fila in mapa["modulo"]])
    x = np.array(mapa["x_mm"])[None, :]
    zz = np.array(mapa["z_mm"])[:, None]
    dentro = x**2 + zz**2 < (0.95 * mapa["radio_mm"]) ** 2
    print(f"\n  puntos dentro de la esfera: {int(dentro.sum())} · con valor: "
          f"{int(np.sum(~np.isnan(z[dentro])))}")
    assert np.all(np.isnan(z[dentro]))


def test_el_residuo_de_fuerza_se_reporta(corrida_rapida):
    """La fuerza neta que la malla no cancela es un diagnóstico, no un detalle: si creciera, el
    resultado dejaría de tener sentido aunque el error contra Faran todavía parezca razonable."""
    residuo = corrida_rapida["metrica"]["residuo_rigido"]
    print(f"\n  |residuo rígido| = {residuo:.2e} N")
    assert residuo is not None and residuo < 1e-4

# --------------------------------------------------------------------------- contraste bajo


def test_el_contraste_bajo_desactiva_la_restriccion_solo():
    """Con contraste bajo la suposición de la restricción (fuerza neta nula) deja de valer.

    Medido: con impedancia 1.08× el agua, el campo sale **53×** el exacto con la restricción puesta y
    1.27× sin ella. Dejarla activa porque "es lo recomendado" sería entregar un número que no existe.
    """
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = elastic_sphere.resolver(
            {"modo": "rapida", "ka": 1.0, "material": "tejido (casi fluido)", "puntos_angulo": 31,
             "fijar_rigidos": True}
        )
    m = r["metrica"]
    print(f"\n  contraste {m['contraste_impedancia']:.2f} · desactivada sola "
          f"{m['restriccion_automatica']} · razón de amplitudes {m['razon_amplitudes']:.2f}")
    assert m["contraste_impedancia"] < elastic_sphere.CONTRASTE_MINIMO
    assert m["restriccion_automatica"] is True
    assert m["fijar_rigidos"] is False
    # sin la restricción el campo queda del orden del exacto, no 53×
    assert m["razon_amplitudes"] < 3.0


def test_un_campo_dispersado_debil_no_se_sella_por_su_error_relativo():
    """Cuando `|p_s|` es una fracción chica de la incidente, el error *relativo* explota por el
    denominador aunque el error absoluto sea bueno. Ahí el estado honesto es `sin_referencia`, con el
    error absoluto a la vista — no un `fuera` sobre un 23 % de algo que vale 1/50."""
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = elastic_sphere.resolver(
            {"modo": "rapida", "ka": 1.0, "material": "tejido (casi fluido)", "puntos_angulo": 31}
        )
    m = r["metrica"]
    print(f"\n  error relativo {100*m['error_relativo']:.1f} % · absoluto {m['error_absoluto']:.2e} Pa "
          f"· sello {r['sello']['estado']}")
    assert r["sello"]["estado"] == "sin_referencia"
    assert m["error_absoluto"] < 0.05          # menos del 5 % de la incidente (1 Pa)

def test_la_resolucion_de_la_onda_de_corte_se_reporta():
    """La malla se elige con la longitud de onda del FLUIDO, pero un sólido blando tiene `λ_s = c_s/f`
    corta. Medido: acero 53.6 elementos por `λ_s`, hueso 31.8, tejido casi fluido 0.5.

    Importa distinguirlo: el 13 % del hueso **no** es falta de resolución de corte (31.8 elementos por
    λ_s), es el modelo de la restricción. Y a la grasa (`c_s = 100`, 1.9 elementos por λ_s) no le
    alcanza con apagar la restricción.
    """
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rapida = elastic_sphere.resolver(
            {"modo": "rapida", "ka": 1.0, "material": "acero", "puntos_angulo": 12}
        )["metrica"]
        blanda = elastic_sphere.resolver(
            {"modo": "rapida", "ka": 1.0, "material": "tejido (casi fluido)", "puntos_angulo": 12}
        )["metrica"]
    print(f"\n  acero: {rapida['elementos_por_lambda_corte']:.1f} elems por λ_s · "
          f"tejido: {blanda['elementos_por_lambda_corte']:.1f}")
    assert rapida["elementos_por_lambda_corte"] > 20
    assert blanda["elementos_por_lambda_corte"] < 1
    assert blanda["lambda_corte_mm"] < blanda["elemento_mm"]
