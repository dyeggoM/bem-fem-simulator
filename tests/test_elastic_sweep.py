"""El barrido en `ka` de la esfera elástica (U7b): rango, contrato y coherencia con el criterio."""

from __future__ import annotations

import pytest

from ultrasound_sim.dashboard import elastic_sweep


def test_los_puntos_salen_equiespaciados_y_en_orden():
    valores = elastic_sweep.puntos_ka(0.5, 2.0, 7)
    assert len(valores) == 7
    assert valores[0] == pytest.approx(0.5)
    assert valores[-1] == pytest.approx(2.0)
    assert valores == sorted(valores)
    print(f"\n  ka: {[round(v, 3) for v in valores]}")


def test_el_rango_rechaza_lo_que_no_tiene_sentido():
    """Los tres casos que un usuario puede escribir de verdad: invertido, fuera del rango de la
    vista, y un número de puntos absurdo. Se rechazan con un mensaje, no se recortan en silencio."""
    with pytest.raises(ValueError):
        elastic_sweep.puntos_ka(2.0, 0.5, 5)      # invertido
    with pytest.raises(ValueError):
        elastic_sweep.puntos_ka(0.1, 2.0, 5)      # debajo de 0.3
    with pytest.raises(ValueError):
        elastic_sweep.puntos_ka(0.5, 4.0, 5)      # encima de 3
    with pytest.raises(ValueError):
        elastic_sweep.puntos_ka(0.5, 2.0, 1)      # un solo punto no es un barrido


def test_el_barrido_devuelve_un_punto_por_ka_y_comparable_con_el_criterio():
    """Tres puntos en modo rápido (~6 s). Se verifica el contrato: una fila por `ka`, el error y el
    sello de cada una, y que el resumen use la tolerancia del registro y no una constante."""
    r = elastic_sweep.barrido_ka({
        "desde": 0.6, "hasta": 1.4, "puntos": 3, "material": "acero",
        "modo": "rapida", "formulacion": "total", "puntos_angulo": 41,
    })
    assert len(r["filas"]) == 3
    assert len(r["barrido"]["ka"]) == r["barrido"]["error_relativo"].__len__() == 3
    for fila in r["filas"]:
        assert fila["error_relativo"] > 0
        assert fila["sello"] in ("validado", "fuera", "sin_referencia")
        assert fila["elementos_por_lambda_corte"] > 0

    from ultrasound_sim.dashboard import validacion

    tol = validacion.CRITERIOS["volumen_esfera_elastica"].tolerancia
    print(f"\n  tolerancia del registro {tol} · peor {100*r['resumen']['peor_error']:.2f} %"
          f" en ka {r['resumen']['peor_ka']:.2f} · dentro {r['resumen']['n_dentro']}/3")
    assert r["resumen"]["tolerancia"] == tol
    # el rápido usa dominio 2.5a: con el ABC de 2º orden el error queda por encima del criterio, y el
    # barrido tiene que decirlo en vez de esconderlo
    assert r["resumen"]["n_dentro"] == 0
    assert all(f["sello"] == "fuera" for f in r["filas"])


def test_la_formulacion_se_rechaza_si_no_existe():
    with pytest.raises(ValueError):
        elastic_sweep.barrido_ka({"formulacion": "magica", "puntos": 2})


def test_comparar_devuelve_las_dos_series_y_donde_cruzan():
    """D51: las dos formulaciones se cruzan, y el barrido lo dice con un número (`cruce_desde_ka`)
    en vez de dejar que el usuario compare dos corridas a ojo."""
    r = elastic_sweep.barrido_ka({
        "desde": 0.8, "hasta": 1.6, "puntos": 2, "material": "acero", "modo": "rapida",
        "comparar": True, "puntos_angulo": 41,
    })
    b = r["barrido"]
    assert b["comparar"] is True
    assert set(b["series"]) == {"total", "dispersada"}
    for form in ("total", "dispersada"):
        assert len(b["series"][form]["ka"]) == 2
        assert len(b["series"][form]["error_relativo"]) == 2
    # el cruce es None (no hay punto donde gane la dispersada) o un ka del rango, nunca un número
    # fuera de la lista
    cruce = b["cruce_desde_ka"]
    print(f"\n  cruce desde ka: {cruce} · total {[round(100*v,1) for v in b['series']['total']['error_relativo']]}"
          f" · dispersada {[round(100*v,1) for v in b['series']['dispersada']['error_relativo']]}")
    assert cruce is None or cruce in b["series"]["total"]["ka"]


def test_sin_comparar_no_se_paga_el_doble():
    r = elastic_sweep.barrido_ka({
        "desde": 0.8, "hasta": 1.6, "puntos": 2, "material": "acero", "modo": "rapida",
        "puntos_angulo": 41,
    })
    assert r["barrido"]["comparar"] is False
    assert r["barrido"]["series"] is None
    assert r["barrido"]["cruce_desde_ka"] is None


def test_el_barrido_se_cancela_entre_puntos():
    """La cancelación del dashboard es cooperativa: el worker consulta entre pasos. Si el barrido no
    la mirara, el botón de cancelar de la UI no haría nada en una corrida de 25 puntos."""

    class TrabajoQueCancela:
        def __init__(self):
            self.llamadas = 0

        def cancelado(self):
            self.llamadas += 1
            return self.llamadas > 1     # deja pasar el primer punto y corta en el segundo

        def avanzar(self, fraccion, mensaje):
            pass

    with pytest.raises(RuntimeError, match="cancelado"):
        elastic_sweep.barrido_ka(
            {"desde": 0.6, "hasta": 1.4, "puntos": 3, "material": "acero", "modo": "rapida",
             "formulacion": "total", "puntos_angulo": 41},
            TrabajoQueCancela(),
        )
