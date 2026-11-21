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


def test_el_mapa_cubre_agua_capas_y_agua():
    """La grilla tiene que incluir el agua de entrada y de salida, no sólo las capas."""
    r = cadena_field.resolver_cadena(_datos([("piel", 2.0), ("perilinfa", 5.0)]))
    x = np.array(r["grilla"]["x_mm"])
    capas = r["capas"]
    assert x[0] == 0.0
    assert x[-1] > capas[-1]["hasta_mm"], "el mapa no llega a la salida"
    assert capas[0]["desde_mm"] > 0, "las capas tendrían que arrancar después del agua de entrada"
    assert np.all(np.isfinite(np.array(r["grilla"]["db"])))
