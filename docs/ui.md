# Guía del dashboard

```bash
PYTHONPATH=src .venv/bin/python -m ultrasound_sim.dashboard --port 8765   # http://127.0.0.1:8765
```

El servidor es local y no tiene autenticación: escucha en `127.0.0.1`. Las corridas largas van a una
**cola con progreso y cancelación** (la cancelación es cooperativa: se atiende entre pasos, no dentro
de un solve).

## Las ocho vistas

| vista | qué resuelve | costo | sello |
|---|---|---|---|
| **campo** | apertura 2D sobre baffle rígido o en espacio libre; fuentes arrastrables | 0.4–1.2 s por mapa | métricas de foco, sin criterio de validación |
| **escena FEM** | capas planas 2D contra la matriz de transferencia | < 0.1 s | `escena_rt_*` |
| **dispersor BEM** | cilindro rígido por BEM **y** FEM contra la serie de Bessel | ~2 s | `dispersor_bem_serie`, `dispersor_fem_serie` |
| **cóclea** | corte transversal de la cápsula (sólido + fluido + membrana) | ~3 s | `cochlea_t_interior`, `cochlea_balance` |
| **barridos** | sensibilidad de la cóclea, frecuencia de la escena, **apertura (FWHM vs pitch)** | 10 s – varios min | según el barrido |
| **3D** | FEM de volumen: **isosuperficie 3D real** (gl3d) + dos rebanadas; fuente puntual contra Green, esfera contra la serie | 2–30 s | `volumen_*` |
| **esfera elástica** | sólido elástico real (ondas P y S) contra Faran, dos formulaciones | 2 s (rápida) / ~30 s (precisa) | `volumen_esfera_elastica` |
| **criterios** | el registro completo de criterios y lo medido en la sesión | instantáneo | — |

## Primer uso

La primera vez que se abre, el dashboard lanza solo un **tour de 9 pasos** que resalta cada control y
explica para qué está. Se cierra con `Esc` o «saltear», y se repite con el botón **tour** de arriba.
El tour **no cambia de vista** para resaltar las otras pestañas: resalta el botón del nav, porque
cambiar de vista dispara solves y un tour que lanza simulaciones para explicarse es peor que no tenerlo.

## Mover las fuentes con el mouse

Dos modos, con botón explícito en el panel:

- **mover con el mouse** (por defecto): arrastrar sobre el mapa mueve la fuente más cercana; en modo
  «arreglo» no hay fuentes sueltas y lo que se mueve es el **foco en y** (el haz se inclina). Doble
  clic agrega una fuente en esa altura.
- **zoom**: el arrastre vuelve a ser el de Plotly.

El arrastre funciona desde cualquier punto del mapa, no sólo sobre las asas de 12 px, y hay asas en
los dos modos (32 en un arreglo de 32 elementos). El cursor cambia a mano abierta para que se note.

*Por qué existía la confusión*: con el `dragmode` de Plotly activo, el gesto lo tomaba el zoom y las
asas sólo existían en modo «fuentes»; arrastrar en modo arreglo no hacía nada y el mapa se acercaba.

## Layout y legibilidad

- **El menú de vistas es una grilla de 2 columnas.** Antes eran diez botones en una fila con `flex:1`
  y en pantallas angostas había que scrollear horizontalmente para llegar a los últimos. «tour» y
  «exportar vista» salieron del nav a su propia fila.
- **La acción principal está fija al pie del panel.** «resolver», «barrer» y «lanzar» se distinguen en
  azul (`.btn.accion`) y su fila se mueve a una barra `sticky` al pie: en una columna de controles
  larga, el botón quedaba fuera de pantalla. En «campo» la barra avisa que **se recalcula solo**.
- **Los tooltips son claros y propios.** Ni el `title` nativo ni los hoverlabels de Plotly: ambos eran
  oscuros sobre oscuro. El globo de los `?` va con fondo `#f2f5f9` y texto oscuro, y los hoverlabels
  de las gráficas usan el mismo esquema.
- **Cada vista dice qué está mirando** en una tarjeta al principio del panel. Las vistas resuelven
  problemas físicos **distintos** (una apertura, capas planas, un cilindro dispersor, un corte de
  cóclea): sin decirlo, parece que el dashboard muestra cosas incoherentes.
- **Los estáticos se sirven con `no-store`.** Con `no-cache` y sin ETag, el navegador podía seguir
  usando el JS viejo y un arreglo de UI "no surtía efecto".

## Flujos de trabajo

**Comparar dos configuraciones.** En `campo`: *fijar referencia* congela el campo actual; cambiá lo
que quieras y el mapa pasa a mostrar la **diferencia en dB** (±20 dB, escala divergente), con el
perfil de la referencia punteado en los cortes. Si cambiás frecuencia o malla, el vector de
coordenadas deja de coincidir y **no se resta**: se avisa en la línea de estado.

**Comparar dos formulaciones.** En `esfera elástica` hay un selector `total` / `dispersada` con el
error medido de cada una, y el barrido en `ka` puede correr **las dos en una sola corrida** y marcar
el `ka` donde cruzan.

**No repetir un solve.** Cada vista cara guarda su último resultado en `localStorage` (12 h) y lo
repinta al volver, diciendo de cuándo es. Restaurar no re-guarda: si lo hiciera, la fecha se renovaría
en cada recarga y un resultado de ayer parecería de ahora.

**Exportar.** `exportar vista` (arriba a la derecha) baja un PNG de cada gráfica visible y un JSON con
sus trazas; no sabe nada de cada vista, lee lo que está en pantalla. En la esfera además hay JSON
completo, CSV de `|p_s|(θ)` y un historial de corridas que se superponen en la gráfica.

**Informe de la sesión.** En `criterios`: un Markdown con las corridas guardadas, sus métricas y su
sello **con su fuente**, más la tabla de criterios. Un criterio sin corrida dice `sin medir`, que no
es un aprobado.

## Trampas conocidas (medidas, no supuestas)

- **El modo rápido de la esfera no sirve para comparar formulaciones.** Con dominio 2.5a el error del
  ABC domina y las dos formulaciones miden ese error, no su diferencia: con el rápido la dispersada
  parecía mejor en `ka` bajo y con el preciso pasa lo contrario (D51).
- **El máximo sobre la curva angular es frágil.** La esfera del FEM es facetada y eso corre los
  lóbulos unos grados; el máximo lo castiga sin que la amplitud esté mal. Por eso cada punto trae
  también la norma L2.
- **El sello de la esfera elástica vale en una ventana (`ka` ≈ 1.1–1.6)** con la formulación total;
  fuera de ahí da 6.8 % y 9.2 %, las dos fuera del criterio (D51).
- **`error_R_pct` de la escena es el error de amplitud**; el fasor va aparte (`error_R_fasor_pct`) con
  el desfase de referencia que se le sacó (D52). Si comparás fases a mano contra el analítico, tenés
  que aplicar ese giro.
