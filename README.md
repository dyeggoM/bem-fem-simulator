# bem-fem-simulator

Simulador de ultrasonido en tejido blando, hueso y fluidos cocleares, con **dos métodos numéricos
independientes** (FEM y BEM) que se contrastan entre sí: donde no coinciden hay un error o una
suposición mal puesta.

El objetivo del proyecto es estudiar **focalización con transductores reales y arreglos de fuentes**,
y decidir con números si un arreglo puede igualar a un transductor monolítico.

## Qué hay implementado

| Módulo | Qué hace |
|---|---|
| `materials.py` | Fluidos y sólidos con procedencia obligatoria; `VERIFY` donde la fuente no existe |
| `meshing.py` | Mallas 1D, rectangulares y coronas polares, con orientación CCW y normales salientes verificadas |
| `analytic/helmholtz1d.py` | Matriz de transferencia para medios estratificados (referencia de R y T) |
| `analytic/scattering.py` | Serie de Jacobi-Anger + Hankel, cilindro infinito 2D |
| `fem/assembly.py` | Ensamblado P1 (1D y 2D), vectorizado sobre mallas no uniformes |
| `fem/helmholtz1d.py`, `fem/helmholtz2d.py` | Helmholtz FEM con atenuación, ABC-1 y BGT-2, campo total y dispersado |
| `fem/elastic2d.py`, `fem/coupled2d.py` | Elastodinámica en plano deformación y acoplamiento fluido-sólido en bloque |
| `bem/helmholtz2d.py` | BEM exterior propio (Green `(i/4)H₀⁽¹⁾`): ecuación estándar, Burton-Miller y CHIEF |
| `aperture.py`, `arrays.py`, `focus.py` | Integral de Rayleigh sobre baffle rígido, arreglos con delays, métricas de foco |
| `transducers.py` + `data/transducers.toml` | Catálogo editable de transductores, con procedencia obligatoria |
| `attenuation.py` | Leyes de atenuación del hueso verificadas contra la fuente, y modelos rechazados |
| `cochlea.py` | Corte transversal: cápsula ótica elástica + fluido interior + ventana de impedancia |

## Verificación

Cada fase tiene un criterio medible y el test lo imprime. Los números no son "anda bien":

- **1D vs analítico**: error L2 `0.399 %`; R de una interfaz contra `(Z₂−Z₁)/(Z₂+Z₁)` a `1e-10`.
- **2D acústico**: onda plana `0.51 %`; reflexión espuria de la ABC `−73 dB`; R/T de una capa al `0.05 %`.
- **2D elástico acoplado**: |T| de placa agua–acero/hueso–agua al `0.079 %` del analítico; Snell con
  conversión de modo a `0.42°`; resonancia de espesor al `0.012 %`.
- **BEM**: cilindro rígido vs serie de Bessel al `0.13 %`; Burton-Miller y CHIEF sobreviven las
  resonancias donde la ecuación estándar se degrada.
- **Contraste FEM vs BEM** (lo que da sentido a tener dos métodos): `0.97 %` a ka=1 y `0.59 %` a ka=2.
- **Focalización**: equivalencia arreglo ↔ apertura continua medida en función del pitch.

```bash
python -m venv .venv && .venv/bin/pip install numpy scipy matplotlib pytest
.venv/bin/python -m pytest -m "not slow"   # ciclo rápido
.venv/bin/python -m pytest                 # incluye barridos
```

## Estructura

```
src/ultrasound_sim/{materials,meshing,attenuation,transducers,cochlea}.py
src/ultrasound_sim/{analytic,fem,bem}/    # un directorio por método
src/ultrasound_sim/dashboard/             # dashboard interactivo local (HTTP + Plotly)
scripts/    bench_frequency.py · cross_validate.py · sensitivity_cochlea.py · focusing_study.py · report.py
tests/      un archivo por módulo, con los criterios de la fase
docs/       plan.md (fases y criterios) · decisions.md (ADR: convenciones y lo descartado)
```

## Dashboard

```bash
.venv/bin/python -m ultrasound_sim.dashboard        # http://127.0.0.1:8765/
```

Sliders de frecuencia, cantidad de elementos, apertura y foco; mapa de intensidad con los cortes
axial y lateral; y un editor de fuentes donde cada elemento tiene posición, ancho, amplitud y fase.
Las fuentes se **arrastran** con el mouse sobre el mapa (doble clic agrega una en esa altura), y la
configuración se **descarga/carga** como JSON o se autoguarda en el navegador — recargar no pierde el
trabajo. El motor es vectorizado (1–2 s por mapa de 120×120, 0.4 s durante el arrastre) y coincide con
`aperture.pressure` a `1e-15`.

Dos cosas que el dashboard dice en voz alta en vez de esconder: cuando el pico real no coincide con
el foco pedido (corrimiento no paraxial, D13) y cuándo el foco pedido está más allá del campo cercano
de la apertura.

El panel **cóclea · corte transversal** resuelve aparte el corte de la cápsula (sólido elástico +
fluido + membrana basilar) y expone los **tres parámetros que la medición todavía no fija** —
espesor de la cápsula, `c_s` del hueso y atenuación α— con un tooltip cada uno que dice qué es, su
rango publicado y cuánto pesa en el resultado (D30, D31).

Hay dos modelos de radiación, y la diferencia entre ellos está medida: *baffle rígido* (integral de
Rayleigh, sólo semiespacio) y *espacio libre* (lámina finita, capa simple, radia a los dos lados).
En `x > 0` dan el mismo módulo con 180° de fase —una identidad, no una aproximación— así que el
supuesto de baffle no sesga el campo delantero; lo único que niega es la radiación trasera, que en el
modelo libre es tan fuerte como la delantera (D28).

### Lo que agregó U8 (para trabajar, no para mirar una vez)

- **Registro de criterios**: la vista `criterios` muestra los 12 criterios del servidor con su
  tolerancia y su fuente, y al lado el número medido en la sesión. Un criterio sin corrida dice
  **sin medir**, que no es un aprobado.
- **Exportar**: un botón transversal baja un PNG de cada gráfica visible y un JSON con sus trazas. En
  la esfera además: JSON completo, CSV de `|p_s|(θ)` y un historial de corridas que se **superponen**
  en la gráfica para comparar formulaciones o materiales.
- **No se repite el solve**: cada vista cara guarda su último resultado (localStorage, 12 h) y lo
  repinta al volver, diciendo de cuándo es.
- **Barrido en `ka` de la esfera elástica**: grafica el error de cada formulación contra Faran con la
  tolerancia del criterio como línea, para ver **dónde** vale cada una y no sólo cuánto se equivoca
  en un punto.

## Limitaciones conocidas (medidas, no estimadas)

- **El material manda sobre la malla.** Con el modelo completo (membrana basilar incluida) el
  espesor de la cápsula explica el `53 %` de la dispersión del resultado, `c_p` el `16 %`, la
  atenuación el `10 %` y `c_s` el `8 %`; por debajo del `5 %` (`ρ`, `Z_r`) no se distingue del piso
  numérico de `0.93 %`. Refinar la malla no cambia ninguna conclusión hasta que los datos de
  material mejoren (D24, D30).
- **No hay datos de atenuación de hueso humano por debajo de 270 kHz**, que es justo la banda de
  diseño. Dos estudios serios (Pichardo 2011, White 2006) difieren por `2.4×–3.3×`.
- **La cápsula ótica no tiene propiedades medidas citables**; los valores usados son un sustituto
  isotrópico declarado, conservador en transmisión.
- **La convergencia de malla no es monótona** en el modelo de cápsula: refinar la pared puede empeorar
  el resultado al cambiar la mecánica discreta.

## Política de datos

Ningún número entra sin fuente. Los materiales llevan `provenance`; lo que la especificación no fija
se marca `VERIFY`; los modelos que no reproducen su propia fuente quedan registrados como
**rechazados** en `attenuation.REJECTED_MODELS`, con el test que lo demuestra.

Licencia: **MIT** — ver `LICENSE`. Copyright (c) 2026 Juan Diego Morales del Corral.
