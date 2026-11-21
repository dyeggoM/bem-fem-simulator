"""La cadena piel→cóclea por el dashboard: medida contra la referencia, y lo que dice cuando no la hay.

Es el primer lugar donde el dashboard usa el acople de N regiones, así que lo que se fija acá es que el
número que muestra la UI sea el que corresponde: cuando todas las capas son fluidas se compara contra la
matriz de transferencia y el error tiene que ser chico; cuando hay un sólido elástico **no** hay
referencia y se dice, en vez de comparar contra una fórmula que no incluye cizalla.
"""

from __future__ import annotations

import numpy as np
import pytest

from ultrasound_sim.dashboard import cadena_field


def _datos(capas, **extra):
    base = {
        "capas": [{"medio": m, "espesor_mm": e} for m, e in capas],
        "frecuencia_hz": 500e3, "ppw": 14, "nx": 80, "ny": 60,
    }
    base.update(extra)
    return base


def test_cadena_fluida_coincide_con_la_matriz_de_transferencia():
    r = cadena_field.resolver_cadena(_datos([("piel", 2.0), ("perilinfa", 5.0)]))
    m = r["metricas"]
    print(f"\n  t medido {m['t_medido']:.4f} · teórico {m['t_teorico']:.4f} · error {m['error_pct']:.2f} %")
    assert m["t_teorico"] is not None
    assert m["error_pct"] is not None and m["error_pct"] < 5.0
    assert not m["hay_solido"]


def test_con_hueso_no_se_inventa_una_referencia():
    """El hueso entra como sólido elástico: la fórmula de capas fluidas no aplica."""
    r = cadena_field.resolver_cadena(_datos([("piel", 2.0), ("hueso_cortical", 5.0)]))
    m = r["metricas"]
    print(f"\n  con hueso: medido {m['t_medido']:.4f} · teórico {m['t_teorico']} · "
          f"referencia: {m['referencia']}")
    assert m["hay_solido"] is True
    assert m["t_teorico"] is None
    assert m["error_pct"] is None
    assert "sin referencia" in m["referencia"]


def test_el_caso_transparente_transmite_uno():
    """Piel de agua: no hay contraste y la transmisión es 1."""
    r = cadena_field.resolver_cadena(_datos([("agua", 5.0)]))
    m = r["metricas"]
    print(f"\n  capa de agua: medido {m['t_medido']:.4f}")
    assert abs(m["t_medido"] - 1.0) < 0.02


def test_capas_invalidas_fallan_con_el_motivo():
    with pytest.raises(ValueError, match="espesor"):
        cadena_field.resolver_cadena(_datos([("piel", 0.0)]))
    with pytest.raises(ValueError, match="medio"):
        cadena_field.resolver_cadena(_datos([("no_existe", 1.0)]))
    with pytest.raises(ValueError, match="capas"):
        cadena_field.resolver_cadena({"frecuencia_hz": 500e3})


def test_el_mapa_usa_el_modulo_del_fasor_y_no_solo_la_parte_real():
    """`LinearTriInterpolator` con un fasor complejo descarta la parte imaginaria con un
    ComplexWarning: el mapa mostraba |Re(p)|, que no es la amplitud. Se detecta porque en las capas el
    campo tiene fase y |Re(p)| != |p|."""
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("error", np.exceptions.ComplexWarning)
        r = cadena_field.resolver_cadena(_datos([("piel", 2.0), ("perilinfa", 5.0)]))
    db = np.array(r["grilla"]["db"])
    print(f"\n  mapa: rango {db.min():.1f} … {db.max():.1f} dB")
    assert db.max() <= 0.0 + 1e-6
    # una onda plana atravesando capas casi no se modula: el mapa es suave y no tiene nulos profundos.
    # Lo que se está comprobando es que no salió ninguna ComplexWarning (arriba) y que la escala es
    # coherente; un mapa de sólo la parte real tendría un rango mucho mayor y sin sentido físico.
    assert db.min() > -20.0


def test_la_estructura_finita_se_puede_mover_y_el_campo_cambia():
    """Con la estructura finita, moverla en x tiene que cambiar el campo; con capas infinitas no."""
    base = _datos([("piel", 2.0), ("perilinfa", 5.0)], estructura_alto_mm=14.0)
    a = cadena_field.resolver_cadena({**base, "estructura_x_mm": 10.0})
    b = cadena_field.resolver_cadena({**base, "estructura_x_mm": 18.0})
    print(f"\n  x=10 mm: t {a['metricas']['t_medido']:.4f} · x=18 mm: t {b['metricas']['t_medido']:.4f}")

    assert a["estructura"]["finita"] is True
    assert a["estructura"]["alto_mm"] == 14.0
    assert a["estructura"]["x_mm"] == 10.0 and b["estructura"]["x_mm"] == 18.0
    assert abs(a["metricas"]["t_medido"] - b["metricas"]["t_medido"]) > 1e-3, (
        "mover la estructura no cambió el campo: se está resolviendo siempre lo mismo"
    )


def test_mover_en_y_saca_la_estructura_del_eje():
    """La pieza también se mueve transversalmente: en el eje deja de estar enfrente."""
    base = _datos([("piel", 2.0), ("perilinfa", 5.0)], estructura_alto_mm=12.0)
    centrada = cadena_field.resolver_cadena({**base, "estructura_y_mm": 0.0})
    corrida = cadena_field.resolver_cadena({**base, "estructura_y_mm": 10.0})
    print(f"\n  centrada: t {centrada['metricas']['t_medido']:.4f} · "
          f"desplazada 10 mm: t {corrida['metricas']['t_medido']:.4f}")
    assert centrada["estructura"]["y_mm"] == 0.0
    assert corrida["estructura"]["y_mm"] == 10.0
    assert corrida["metricas"]["t_medido"] > centrada["metricas"]["t_medido"], (
        "con la pieza fuera del eje la onda pasa casi sin obstáculo por el eje"
    )


def test_las_capas_infinitas_siguen_siendo_el_default():
    r = cadena_field.resolver_cadena(_datos([("piel", 2.0), ("perilinfa", 5.0)]))
    assert r["estructura"]["finita"] is False
    assert r["estructura"]["alto_mm"] is None
    assert r["metricas"]["error_pct"] is not None


def _fuentes_foco(n=16, apertura_mm=24.0, foco_mm=30.0, f=500e3):
    k = 2 * np.pi * f / 1500.0
    salida = []
    for i in range(n):
        y_mm = -apertura_mm / 2 + apertura_mm * (i + 0.5) / n
        salida.append({"y_mm": y_mm, "ancho_mm": apertura_mm / n, "amp": 1.0,
                       "fase_rad": -k * float(np.hypot(foco_mm, y_mm)) * 1e-3})
    return salida


def test_se_puede_excitar_cada_la_cadena_con_un_arreglo_de_fuentes():
    """La cadena no tiene que estar limitada a una onda plana: con un arreglo, el campo enfoca y
    atraviesa las capas. Y ahí NO hay coeficiente exacto contra el cual comparar, así que no se inventa."""
    r = cadena_field.resolver_cadena({**_datos([("piel", 2.0), ("perilinfa", 5.0)]),
                                      "fuentes": _fuentes_foco(), "alto_mm": 40.0})
    m = r["metricas"]
    print(f"\n  excitación {m['excitacion']} · {m['n_fuentes']} fuentes · "
          f"pico detrás {m['pico_detras']:.3e} · t_medido {m['t_medido']}")
    assert m["excitacion"] == "fuentes"
    assert m["n_fuentes"] == 16 and m["hay_fuentes"] is True
    assert m["t_medido"] is None, "con fuentes el 'transmitido' no es un coeficiente"
    assert m["pico_detras"] is not None and m["pico_detras"] > 0
    assert "arreglo de fuentes" in m["referencia"]
    assert np.all(np.isfinite(np.array(r["grilla"]["db"])))


def test_una_fuente_fuera_del_dominio_lo_dice():
    base = _datos([("piel", 2.0)], alto_mm=20.0)
    with pytest.raises(ValueError, match="ninguna fuente cae"):
        cadena_field.resolver_cadena({
            **base,
            "fuentes": [{"y_mm": 500.0, "ancho_mm": 2.0, "amp": 1.0, "fase_rad": 0.0}],
        })


def test_a_alta_frecuencia_se_ajusta_la_malla_en_vez_de_fallar():
    """A 600 kHz con el dominio por defecto la malla pasaba el tope y la corrida se caía con un
    ValueError. Ahora se baja el ppw para entrar y se avisa, en vez de bloquear."""
    r = cadena_field.resolver_cadena({**_datos([("piel", 2.0), ("perilinfa", 5.0)]),
                                      "frecuencia_hz": 600e3, "ppw": 14})
    m = r["metricas"]
    print(f"\n  600 kHz: {m['nodos']} nodos · ppw {m['ppw_efectivo']} · aviso {m['aviso_malla']}")
    assert m["nodos"] <= cadena_field.MAX_NODOS
    assert m["ppw_efectivo"] <= 14


def test_el_mapa_cubre_agua_capas_y_agua():
    """La grilla tiene que incluir el agua de entrada y de salida, no sólo las capas."""
    r = cadena_field.resolver_cadena(_datos([("piel", 2.0), ("perilinfa", 5.0)]))
    x = np.array(r["grilla"]["x_mm"])
    capas = r["capas"]
    assert x[0] == 0.0
    assert x[-1] > capas[-1]["hasta_mm"], "el mapa no llega a la salida"
    assert capas[0]["desde_mm"] > 0, "las capas tendrían que arrancar después del agua de entrada"
    assert np.all(np.isfinite(np.array(r["grilla"]["db"])))
