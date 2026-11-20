# Decisiones (ADR corto)

Una entrada por decisión: contexto → decisión → consecuencia. Si algo se revierte, se anota, no se
borra.

## D1 — Layout del paquete

`src/ultrasound_sim/` con `analytic/` (soluciones exactas usadas como referencia), `fem/`, `bem/`,
`materials.py`, `meshing.py`. `tests/` espeja esa estructura.

*Consecuencia*: cada método tiene su carpeta y puede desarrollarse en paralelo sin conflictos; las
referencias analíticas viven aparte porque también se usan en los tests.

## D2 — FEM propio (P1/P2 Galerkin) en vez de scikit-fem/dolfinx, por ahora

La parte difícil de las fases 1–2 no es el ensamblado, es la frontera absorbente y el acoplamiento
fluido-sólido. Con ~300 líneas de Galerkin P1 sobre mallas triangulares (numpy/scipy) se controla
cada matriz, cada término de borde y la convención de signos, y se gana independencia real frente al
BEM para la Fase 4. scikit-fem/dolfinx se reevalúan en Fase 5/6, cuando la malla y las estructuras
finas dominen el problema.

*Consecuencia*: menos dependencias, más control; se paga con ensamblado propio y sin mallador de
calidad (mallas estructuradas y, más adelante, Gmsh para geometrías reales).

## D3 — BEM acústico 2D propio primero; bempp-cl queda como opción

bempp-cl exige OpenCL y su instalación no es trivial; el plan permite "bempp-cl o el tuyo". El BEM
2D acústico (colocación/Galerkin con `G = (i/4)H₀⁽¹⁾`) es acotado y verificable contra la serie de
Bessel del cilindro. Con una implementación propia el contraste FEM/BEM de la Fase 4 es realmente
independiente.

*Consecuencia*: se implementa Green, matrices de capa simple/doble, condición de radiación y
resolución de sistema denso. Si el tamaño crece, se migra a bempp-cl sin cambiar la interfaz de
`bem/`.

## D4 — Convención de tiempo y de atenuación (única en todo el repo)

- Convención: **e^(−iωt)**, ondas salientes ∝ e^(+ikr). Helmholtz: ∇²p + k²p = 0.
- Atenuación por ley de potencias: α(f) = α₀·(f/f_ref)^n  [Np/m], y el número de onda complejo
  **k = ω/c + i·α(f)** (parte imaginaria positiva ⇒ e^(ikx) decae con x, coherente con e^(−iωt)).
- Green 2D: G = (i/4)·H₀⁽¹⁾(k·r).
- Unidades SI internas (m, kg, s, Pa); la API pública acepta Hz y devuelve Pa.

*Consecuencia*: cualquier módulo que se desvíe de esto rompe el contraste FEM/BEM. Los tests de
`test_conventions.py` fijan el signo.

## D5 — Valores de material con procedencia obligatoria

`Material` lleva `provenance` (de dónde salió el número). Del plan vienen: agua/perilinfa/endolinfa
c=1500 m/s, ρ=1000; hueso c≈3000–4000, ρ≈1900–2000 (**ballpark**, el plan pide verificar). Lo que
no esté en el plan o en una fuente citada **no se inventa**: se declara `VERIFY` y se usa sólo en
tests adimensionales.

*Consecuencia*: no hay "números plausibles" silenciosos en la librería; la atenuación de hueso se
parametriza (α₀, n) y sus valores hay que traerlos de literatura o del usuario.

## D6 — Tolerancias de aceptación

Se fijan por fase en `docs/plan.md` (Fase 0: <1 % L2; Fase 1: 2 % en R/T; Fase 3: 2 % vs serie;
Fase 4: 2 % FEM vs BEM). Los tests fallan por encima de eso, no "cerca".

## D7 — Python del repo

`venv` local del repo (`.venv`) creado con el Python 3.14.9 que trae Hermes, porque el sistema no
tiene `pip`. Dependencias: numpy, scipy, matplotlib; `pytest` como dev. Gmsh/PyVista/dolfinx no se
instalan hasta que una fase los pida explícitamente.

## D8 — Techo inicial 100 kHz, hoja de ruta a >1 MHz, focalización desde 250 kHz

Decisión de Diego: empezar limitado a **100 kHz** (los transductores de 40 kHz son además los más
económicos), con hoja de ruta explícita para soportar >1 MHz, optimizando en el camino.

Física que hay que tener presente al ampliar: con cóclea ≈ 9 mm, a 40 kHz λ=37.5 mm (λ/4 del blanco)
y el FWHM de un foco realista (D=10 mm, F=20 mm) es ~75 mm: *no existe focalización*, la comparación
mediría iluminación. A 250 kHz son 12 mm, a 1 MHz 3 mm. Por eso los tests de focalización (Fase 4B)
corren en 250 kHz – 1 MHz; a 40–100 kHz se verifica lo contrario: que el campo llega casi uniforme
(predicción refutable, no un caso de descarte).

*Consecuencia*: las fases 0–5 se cierran a ≤100 kHz, donde la malla es baratísima (20 mm de dominio
a 100 kHz con ppw=20 son 1 458 elementos). La Fase 4B y la escalera de optimización (D10) son las que
habilitan subir.


## D9 — Techo de 3D con el stack actual: ~250 kHz

Con ppw=10 sobre 20 mm de dominio: 3D a 250 kHz son 222 k elementos; a 1 MHz, 14.2 M. Con un FEM
P1 en numpy/scipy el techo práctico está en ~250 kHz. Por encima: FEniCSx+MPI, k-Wave
pseudoespectral, o BEM exterior (escala con superficie). En 2D, 1 MHz son 35.5 k elementos: entra
sin cambios.

*Consecuencia*: la Fase 6 (3D) se planifica a ≤250 kHz salvo que se cambie de stack, y esa decisión
se toma con un benchmark, no por preferencia.

## D10 — Escalera de optimización para subir de 100 kHz a >1 MHz (medida, no supuesta)

Base: `scripts/bench_frequency.py`, ducto con paredes Neumann, dominio de 20 mm, P1.

Medición que manda: con P1 el error se acumula con las longitudes de onda recorridas, no con el
tamaño de la malla. A 100 kHz (1.33 λ) ppw=20 da 2.10 %; a 1 MHz (13.3 λ) ppw=20 da 19.64 %. Bajar
de 2 % a 1 MHz pide ppw ≈ 60 → ~1.3 M elementos 2D. La conclusión es que *el límite es el orden del
elemento*, no la memoria.

Escalones, en orden de esfuerzo/beneficio:

1. **`splu` en vez de `spsolve`** (hecho en el bench): a 1 MHz, 11.56 s → 0.68 s. 17× gratis. Es sólo
   elegir la API correcta.
2. **Reutilizar la factorización entre excitaciones** (lo que pide un arreglo: misma matriz, N rhs):
   28.8× con 4 rhs a 1 MHz, creciente con N. Requisito de diseño de la Fase 4B.
3. **Elementos P2 (o p-FEM)**: el error de dispersión pasa de O((kh)²) a O((kh)⁴), así que el ppw
   necesario cae en raíz cuadrada y el costo con él. **Aún sin medir**: es el primer ítem a validar
   contra el mismo ducto, con el mismo criterio de 2 %.
4. **PML en vez de ABC de primer orden**: la ABC falla a incidencia rasante — medido: ABC en los
   cuatro lados de un cuadrado da 45 % de error a 8 λ, mientras en ducto (incidencia normal) da
   0.5–2 %. Para dominios abiertos en 2D hace falta PML.
5. **BEM exterior para el problema de foco**: escala con la superficie; ya verificado a 0.24 % de
   error en ka=5 con 128 paneles. Es el camino natural para 250 kHz–1 MHz en dominio no acotado.
6. **3D por encima de 250 kHz**: FEniCSx+MPI o k-Wave pseudoespectral. Decisión con benchmark.

*Consecuencia*: cada escalón se acepta sólo si mide < 2 % en el mismo problema y mejora el tiempo
medido. No se sube de frecuencia cambiando de métrica.

## D11 — Dos invariantes de la formulación FEM y una limitación del BEM

Descubiertos implementando las fases 1 y 3; quedan como invariantes porque rompen los resultados de
forma silenciosa.

1. **Peso 1/ρ obligatorio cuando el contraste es de densidad.** La forma `∫(∇p∇v − k²pv)` sólo es
   correcta con ρ uniforme. Con agua/acero el R sale 0.43 en vez de 0.937 y *converge* a ese valor
   equivocado al refinar: no hay error visible, sólo un resultado mal. La forma fuerte correcta es
   `∇·((1/ρ)∇p) + (k²/ρ)p = 0`, y el FEM la implementa con `element_rho`/`edge_rho`. El default
   `ρ=1` deja la ecuación en la forma simple, así que pasar materiales con ρ distinta sin pasar ρ es
   el error a evitar.
2. **Cara de sombra ≠ cara iluminada en la ABC.** En las caras de sombra el campo total ya es
   saliente y no hay término incidente; aplicar el término fuente ahí mete una fuente espuria (llega
   a R≈1 de reflexión falsa). El criterio es `n·d < 0` sobre la dirección de la onda incidente.

3. **Resonancias espurias del BEM de segunda especie** (limitación conocida, no bug): en los ceros de
   `J_n(ka)` la ecuación integral directa se vuelve singular y la solución explota. ka=1,2,5 no las
   tocan, pero un *barrido de frecuencia* —que es justo lo que pide la Fase 4B— sí las cruza. Antes
   del barrido hay que pasar a Burton-Miller o CHIEF, y el test de aceptación es que la solución no
   cambie entre formulaciones en un ka resonante.


## D12 — La ABC de primer orden sobre campo total no vale en contornos curvos

Encontrado por el contraste de la Fase 4 (no por un test unitario), con números:

- Formulación de campo total + ABC de 1er orden en la corona FEM: *18.2 %* de error a ka=1 y
  *17.8 %* a ka=2 contra la serie exacta, mientras el BEM daba 0.05 %.
- *No converge con la malla*: ppw=20 → 18.2 %, ppw=40 → 18.0 %.
- *Sí mejora al alejar el contorno*: R = 2a/3a/4a → 19.9/18.2/18.1 %, y R=6a → 13.7 %.
- Desglose por cara: *4.3 % en la iluminada y 17–24 % en la de sombra*, máximo en θ=0.

Causa: en la cara de sombra el código impone `∂p/∂n = i·k·p` sobre el campo TOTAL. Eso es exacto sólo
si la onda llega normal al borde (cierto en una caja alineada: por eso la Fase 1 medía 0.05 %) y falso
en un contorno curvo, donde `n·d ≠ 1`. La máscara `illuminated` evita el problema en ductos y lo
crea en geometría curva.

Solución adoptada: resolver el campo **dispersado** `p_s`, que es puramente saliente en todo el
contorno, con la onda incidente entrando como dato Neumann en la superficie rígida
(`∂p_s/∂n = −∂p_i/∂n`). Forma débil:

    ∫(∇p_s·∇v − k²p_s v)/ρ − ∫_ABC (i k/ρ) p_s v = *−* ∫_Γ (1/ρ)(∂p_i/∂n) v

El signo del lado derecho es *menos*. Verificación independiente: con `+` la solución 1-D da
`p_s = −exp(+ikx)` y error relativo exacto *2.000*; con `−`, 1.4e-05. En el anillo, el signo invertido
da 78 % de error.

Segundo hallazgo: con la formulación correcta, *la ABC de 1er orden tiene un piso de ~6.5 %* en
contorno curvo (ppw 20/40/80 → 6.65/6.43/6.39 %): es reflexión de la onda cilíndrica, no
discretización. Para el criterio de 2 % hace falta el término de curvatura (BGT-2, `abc_order=2`),
que es el que usan los tests de Fase 4 con `abc_center` en el eje del cilindro.

*Consecuencia*: `fem.solve()` (campo total) queda para problemas tipo caja/ducto, donde es válido y
está validado; todo problema con contorno curvo o dispersor usa `fem.solve_scattered()`. Esta
distinción tiene que estar en el docstring de los dos, y está.

## D13 — Fuentes, arreglos y focalización: tres cosas que el cálculo independiente corrigió

De la Fase 4B, verificadas con una cuadratura propia (trapecio denso) separada del código del paquete.

1. **Signo de la fase de apuntamiento.** La especificación que escribí decía `φ_j = −k(y_j·d)` y está
   *mal*: con `e^{−iωt}` el campo lejano lleva `exp(−i k y sinθ)`, así que para apuntar a `d` hace
   falta `φ_j = +k(y_j·d)` — que es además lo que sale de la ley de foco `φ_j = −k|F−y_j|` con
   `d = (F−y_j)/|F−y_j|`. Verificado: con `+k` el lóbulo principal cae en −29.38° para un pedido de
   −30°; con `−k` apunta al espejo. La ley de *foco* (`−k|F−y_j|`) estaba bien.

2. **El foco geométrico no es el máximo axial a 100 kHz.** D=60 mm (4λ), F=40 mm, f/#=0.67,
   semiángulo 37°: el máximo axial cae en *29.5 mm*, 10.5 mm hacia la apertura. Tres vías
   independientes coinciden (trapecio directo 29.5 mm; cuadratura del paquete 30.8 mm; la misma
   apertura a 500 kHz enfoca en 39.4 mm). Consecuencia: el criterio "pico a menos de 1 mm de F" no
   es físico a 100 kHz y se reemplazó por banda medida + dirección + "en F la intensidad a menos de
   1 dB del máximo" (meseta).

3. **La ganancia de N elementos no es 20log10(N) fuera del régimen paraxial.** Medido en f/#=0.67:
   15.95 dB contra 18.06 dB teóricos, porque los elementos de borde iluminan el foco a 37° y su
   propia directividad los atenúa (|p_j| cae a 0.54 en el borde, y las fases quedan dentro de 5.6°).
   La identidad 20log10(N) sólo se recupera cuando la apertura total es pequeña frente a F.

Y la equivalencia arreglo↔transmisor, medida con la métrica que importa (FWHM en el foco real):
λ/8 → +1.5 %, λ/4 → +6.3 %, λ/2 → +27 %, λ → +48 %. Por *correlación de campo en la región* la
exigencia se relaja (λ/2 ya da 0.998): **la métrica define la respuesta, así que cada número se
reporta con su métrica.**

## D14 — Catálogo de transductores en TOML, con procedencia obligatoria

El usuario carga sus equipos en `data/transducers.toml`; `transducers.load_catalog()` los convierte
en `TransducerSpec` y de ahí en aperturas con foco o apuntamiento. Decisiones:

- **Formato TOML** y no YAML/JSON: comentarios, legible a mano, y `tomllib` está en la stdlib
  (Python ≥3.11) — cero dependencias nuevas.
- **`source` obligatorio y no vacío**: un número sin procedencia no entra, igual que en D5. Un campo
  faltante falla nombrando el equipo y el campo; no hay defaults silenciosos.
- **`verified = false` por defecto**, y el test del repo exige que *todo* lo que viene en el archivo
  versionado esté sin verificar hasta que exista un datasheet. Así nadie confunde el baseline
  sintético con un producto.
- **Kerf no soportado**: si `element_width_m ≠ pitch_m` el importador tira error en vez de
  aproximar elementos contiguos. Cuando el modelo soporte separación, se agrega ahí.
- El único equipo versionado es `referencia_lambda8_100k`, marcado *sintético* y trazable a nuestro
  estudio de Fase 4B; reproduce su FWHM (8.65 mm vs 8.77 mm del estudio), así que sirve de baseline
  y de test de integración del importador.

## D15 — Resonancias del BEM: el criterio intuitivo estaba mal, y por qué importa igual

Yo esperaba que la ecuación estándar fallara >20 % en ka = cero de J_n. **No ocurre** en este
problema, y la razón es física: el modo espurio interior es `J_n(ka)·cos(nθ)`, y la proyección de una
onda plana sobre él vale `∝ J_n(ka) = 0` — la compatibilidad de Fredholm se satisface por simetría, así
que el modo no se excita. Medido: *0.45 %* en el cero de J₀ (2.4048) y *8.9 %* en el de J₂ (5.1356);
ese 8.9 % **no converge** al refinar (8.6 → 9.2e-2 de N=64 a N=512): no es discretización, es el modo
espurio contaminando sin explotar. El fallo O(1) sólo aparece a décimas de por mil del cero (28.6 % en
ka=5.1433) y depende de la malla, así que no sirve como aserción.

Lo que sí es robusto y quedó medido: el *condicionamiento* del sistema estándar explota al acercarse a
la resonancia (3.4e2 con 64 paneles → 2.1e4 con 512, ∝1/h) mientras BM y CHIEF se quedan en ~7. Esa
es la señal de la enfermedad, no el campo.

**Consecuencia para el proyecto**: para un círculo con onda plana BM/CHIEF no cambian gran cosa (peor
caso del barrido 0.5–8: estándar 1.67 %, BM 1.07 %, CHIEF 0.63 %). Pero la Fase 5 tiene geometría
*no* circular y fuentes *descentradas*: ahí la simetría que salva al estándar desaparece y la
resonancia se excita de verdad. Por eso la opción queda implementada y expuesta al usuario
(`equation="standard"|"burton_miller"|"chief"`), con la estándar como default por compatibilidad.

Invariante que se usa como control interno: `W·1 = 0` en toda curva cerrada, que se cumple a *1e-12*
para N=32…256 con la sustracción estática del término `1/z` de H₁ y valor principal de Hadamard en el
auto-panel. Si esa identidad se rompe, el hipersingular está mal regularizado.

## D16 — Dos niveles de verificación

La suite ya pasaba los 6 minutos por los barridos (BEM 60 puntos y resonancia acoplada). Se separó:

- *ciclo rápido*: `pytest -m "not slow"` → *130 tests en 79 s*.
- *completo*: `pytest` → agrega los barridos pesados (marcados `slow`).

*Consecuencia*: la verificación que se corre en cada cambio es la rápida; los barridos se corren antes
de cerrar una fase. Se acepta el riesgo de que el barrido lento avise tarde, a cambio de no dejar de
correr los tests.

## D17 — Tres parámetros del solver acoplado que cambian el resultado, no son cosmética

Implementando la Fase 2 aparecieron tres cosas que, si se omiten, dan un número *plausible y
equivocado*:

- **`dirichlet` (rodillo `u_y = 0`)**: sin él, la placa del ducto transmite con una velocidad efectiva
  ≈ 8 % menor y la comparación contra la fórmula con `c_p` falla sin que nada "explote". En un ducto
  con incidencia normal el confinamiento lateral es parte del problema físico; el solver no puede
  adivinarlo.
- **`periodic_pairs`** (incidencia oblicua): con `H = λ_f/sinθ` la condición periódica es exacta y
  elimina las reflexiones de pared que contaminaban el FFT de ángulos. Con paredes rígidas en su lugar,
  el ángulo medido se desplaza.
- **`frequency`**: ni `Material` ni `ElasticSolid` la llevan (sólo `wavenumber(f)`), así que hay que
  pasarla; si se omite se deduce de `p_inc.k`, y eso sólo es válido si hay onda incidente.

*Consecuencia*: los tres son keyword-only y documentados; los tests de Fase 2 los ejercitan, y el test
de estado de interfaz (p y u contra el analítico, 0.86 % y −0.49°) es el que fija la convención de
signos del acoplamiento.

## D18 — El criterio de aceptación tiene que coincidir con la topología del modelo

Cuarta vez que un criterio mío de escritorio cae al medirlo, y esta es la lección más útil de todas:

Pedí que la cápsula ótica, con curvatura despreciable, transmitiera como una placa plana. **No puede.**
La cápsula es una *cavidad cerrada sin pérdidas*: el campo interior es una onda estacionaria que
recibe energía por toda la circunferencia, así que la "transmisión al fluido interior" no converge al
`T` de una placa — oscila ±30 % con ±5 % de frecuencia. Medido: +19 % de desvío a a/t=40 y tendencia no
monótona con la curvatura. El código está bien: con la pared casi transparente da −0.5…−1.5 % contra la
placa, y el campo dispersado coincide con la serie exacta del cilindro rígido al 0.8 % en potencia.

Regla que queda: **un coeficiente de transmisión necesita un lado transmitido abierto** (o absorción
dentro del dominio). Si la topología es cerrada, lo que se mide son estados estacionarios, no
coeficientes.

*Consecuencia fisiológica, que es la buena noticia*: la cóclea real **no** es una cavidad cerrada.
Tiene ventana oval (entrada) y ventana redonda (salida de presión) y la membrana basilar como
partición. Es decir, la Fase 5b no es un adorno de realismo: **es lo que hace medible la transmisión**.
Y como la ABC de primer orden en contorno curvo aporta 17–21 % en la cara de sombra sin mejorar con el
tamaño del dominio, antes de 5b hay que subir la frontera a BGT-2 o PML (D10, escalón 4).

## D19 — Atenuación del hueso: tres fuentes verificadas, dos modelos rechazados

Diego pasó un documento con cuatro parametrizaciones. Se contrastó cada una contra la fuente primaria
antes de meterla:

- **Pinton et al. 2012** (IEEE TUFFC 39:299, PMID 22225300): *13.3 ± 0.97 dB/cm* medidos en muestras
  de ~8 mm a 1 MHz, dependencia lineal. **Verificada.** Ojo: de esos 13.3, la absorción pura es
  2.7 dB/cm; el resto es reflexión, dispersión y conversión de modo.
- **White et al. 2006 / Treeby & Cox 2014** (citado en Robertson 2017): *8.83 dB/cm/MHz^1.43*.
  **Verificada.**
- **Pichardo et al. 2011** (Phys Med Biol 56:219): en lugar de una ley, se cargan sus *tres puntos
  medidos* (33.8 / 231 / 334 Np/m a 270 / 836 / 1402 kHz). **Verificada** como dato.
- **Rechazado: `27.31·f^1.722`** que circulaba como "Pichardo". No reproduce sus datos: a 1.4 MHz da
  **48.9 dB/cm** contra los **29.0** medidos. Queda en `REJECTED_MODELS` con la razón y un test que lo
  documenta, para que no vuelva por copiar y pegar.
- **Rechazado: `20·f`**: el propio documento lo declara modelo ilustrativo.

**Hallazgo**: los tres puntos de Pichardo *no son una ley de potencia*. Entre 270 y 836 kHz exigen
n ≈ 1.68 y entre 836 y 1402 kHz n ≈ 0.72; un ajuste global (n = 1.44) se desvía *−20 % en el punto
medio*. Por eso la recomendación es usar los puntos medidos con interpolación por tramos, y el ajuste
queda marcado como aproximación. El solver sólo acepta leyes de potencia hoy: usar la tabla requiere
un `k(f)` tabulado, que es trabajo pendiente.

**La trampa del factor 2** (corregida con Sim-Datos): un coeficiente en dB/cm es una pérdida *en dB*,
y el dB se define con 20·log para amplitud y 10·log para intensidad, así que **el número en dB/cm es el
mismo para las dos** (13.3 dB/cm de amplitud son 13.3 dB/cm de intensidad). Lo que se duplica es el
*exponente en Np*: la intensidad cae como `exp(−2αz)` con el mismo α. Al código va α_amp = 153.1 Np/m.
Duplicar los dB (mi error inicial) infla el informe por 2 sin romper ningún test.

**Impacto medido** (placa de 6 mm, agua–cráneo 2800 m/s / 1900 kg/m³–agua, incidencia normal):

| f | \|T\| sin pérdidas | \|T\| con Pinton | \|T\| con White |
|---|---|---|---|
| 100 kHz | 0.5320 | 0.5091 (−0.38 dB) | 0.5264 (−0.09 dB) |
| 500 kHz | 0.8161 | 0.4729 (−4.74 dB) | 0.6206 (−2.38 dB) |
| 1 MHz   | 0.6168 | 0.2722 (−7.11 dB) | 0.3660 (−4.53 dB) |

Pérdida sola en 6 mm (Pinton, 100 / 250 / 500 kHz / 1 MHz): 0.8 / 2.0 / 4.0 / 8.0 dB.
*Consecuencia*: a 100 kHz el cráneo es casi transparente (0.4 dB de amplitud más allá de las
interfaces); a 1 MHz domina el problema. Eso respalda el techo de 100 kHz de D8, pero por una razón
distinta de la que suponíamos: no es que 1 MHz sea caro de calcular, es que el hueso no deja pasar.

## D20 — Consulta a Sim-Datos: qué hay de verdad, y el hueco justo en nuestra banda

Consulta única (2026-11-19) por las seis magnitudes que el repo tenía sin fuente. Resultado por bloque,
con la respuesta resumida y lo que implica:

1. **Atenuación ósea 40 kHz–1 MHz.** No existe un α intrínseco citable para cráneo humano ni cápsula
   ótica *por debajo de 270 kHz*. Puntos disponibles: Pichardo 2011 (33.8±7.6 / 231±15 / 334±42 Np/m a
   270 / 836 / 1402 kHz; **son coeficientes efectivos por inversión, no mediciones locales, y las
   muestras estaban fijadas en formol** pese al "freshly excised" del título) y White 2006 (14±17 /
   53±43 / 70±28 Np/m a 272 / 548 / 840 kHz). *Discrepan 2.4× y 3.3× — no se promedian.* No hay n
   defendible en 40–100 kHz. Único dato cercano: costilla porcina a 40 kHz, 3.4±0.9 dB/cm (Suchkova
   1998), que es pérdida efectiva con interfaces, no absorción volumétrica.
   → *Implicación*: **en nuestra banda de diseño (≤100 kHz) la atenuación del hueso es desconocida, no
   "pequeña"**. Cualquier número que usemos ahí es extrapolación y hay que decirlo en cada informe.

2. **Cápsula ótica.** No hay conjunto experimental citable de c_P, c_S, ρ y α_S. Nuestro
   3500/1900/1950 m/s/m/s/kg/m³ es un *sustituto isotrópico declarado*, no un dato medido. Referencias
   del hueso temporal: c_P = 1752–3285 m/s a 1.03 MHz (Ammi 2008), cortical ρ ≈ 1868–1922 kg/m³
   (Peterson & Dechow 2003); sustituto de calvaria 2820±40 m/s y c_S = 1500±140 m/s (White 2006).
   → 3500 m/s está en el extremo alto: nuestro modelo es *conservador* en transmisión (más impedancia,
   menos pasa). Queda anotado en `materials.BONE_CORTICAL`.

3. **Tejido blando.** Gu & Jing 2019 (tabla I): grasa 0.48, músculo 1.09, conectivo 1.57
   dB·cm⁻¹·MHz⁻ⁿ con n = 1.1 (pérdidas de amplitud). Piel: sin par α0/n verificado; usarlo en
   40–100 kHz es extrapolación. Mast 2000 para relaciones empíricas.
   → `SOFT_TISSUE` sigue en `VERIFY`: ahora sabemos *por qué* y con qué se reemplazaría.

4. **Fluidos cocleares.** No hay mediciones diferenciadas citables: perilinfa/endolinfa como agua es un
   *supuesto acuoso publicado* (Grieser 2016). Agua a 37 °C: c ≈ 1524 m/s, ρ ≈ 993.3 kg/m³,
   α ≈ 0.00139·f_MHz² dB/cm (n = 2). Error de usar 1500/1000: c −1.57 %, ρ +0.67 %, Z −0.91 %.
   → No corresponde fijar α = 0 como propiedad: es idealización invíscida, y las pérdidas por
   interacción con membranas y paredes no están en la absorción del agua.

5. **Membrana basilar y ventanas (Fase 5b).** Espesor 1.26–1.92 µm basal y 0.53–0.89 µm apical
   (Hrncirik 2023); E aparente 9.7 / 6.0 / 6.4 MPa (Ishii 1995, *un solo espécimen*); ventana redonda
   69.4±4.3 µm y E 9.8 MPa; ρ = 1100 kg/m³ es *asignación de modelo* (Tubelli 2025), no medición. La
   ventana oval no es una membrana libre equivalente (platina del estímulo + ligamento anular).
   → Para 2D no hay conjunto validado: una rigidez por unidad de longitud exige declarar ancho fuera
   del plano, apoyos y pretensión. **5b tiene que declarar sus supuestos, no heredar constantes.**

6. **Transductores.** Aire 40 kHz: emisor de 10 mm, arreglo 16×16 = 256 con pitch 10 mm (Marzo y
   Drinkwater 2019). Contacto 40 kHz: Langevin de 45 mm de diámetro radiante (ficha técnica). Terapéutico
   transcraneal: aperturas 200–300 mm y 256–1024 elementos, F = 150 mm a 250/500 kHz. Sin rango de pitch
   verificable para esos equipos (y en arreglos no periódicos el pitch no existe: hacen falta coordenadas).
   → El catálogo (`data/transducers.toml`) usa estos rangos como plantilla, con `verified = false`.

*Consecuencia de proyecto*: a 1 MHz, la incertidumbre entre estudios de atenuación (2.4–3.3×) **domina
cualquier error numérico del solver** (0.1–13 %). No tiene sentido refinar la malla para bajar de 2 % si
la física del material está acotada al 300 %. Y a ≤100 kHz falta el dato directamente.

## D21 — BGT-2 funciona, el piso de 8–9 % es el acoplamiento P1 (y refinar puede empeorar)

Implementado `abc_order=2` en el solver acoplado. Resultados medidos:

- **El port es exacto**: `‖B_acoplado − B_helmholtz2d‖/‖B‖ = 1.5e-16` sobre un círculo con kR = 3.
- **Fluido homogéneo, sin dispersor**: el error contra la onda plana pasa de 20.55 % (ABC-1) a **9.00 %**
  (BGT-2) con ppw=20, y a **4.14 %** con ppw=30. La ABC-1 *no converge* (20.55 → 18.75 %): es defecto de
  modelo, no de malla. BGT-2 sí converge.
- **Cápsula casi rígida**: el sesgo de sombra desaparece — 16.93 % → **9.38 %** y 21.87 % → **9.72 %**
  (+5.1 y +7.0 dB) y las caras quedan parejas (9.4 vs 8.3 %).
- **No-regresión**: la placa de Fase 2 da idénticos |T| con las dos ABC (acero 0.29345, hueso 0.79639;
  desvío +0.001 %), residuo 1.4e-12, mismo tamaño de sistema (45 181 gdl) y tiempo indistinguible.

**Pero el criterio <5 %/<7 % no se alcanza, y no es la frontera.** La prueba: con una cáscara de
impedancia *casada con el agua* —sin dispersor alguno— el mismo modelo da 9.00 % (ppw 20) → 6.27 % (24)
→ 4.14 % (30), el mismo valor y la misma tendencia que el residuo de la cápsula. El piso es la
**discretización P1 del acoplamiento**, no la condición de contorno.

**Y un aviso que hay que tener escrito**: refinar la pared de la cápsula *empeora*. A `h ≤ t/3`
(4 elementos en el espesor) el error salta a **32–36 % con las dos ABC y en todas las caras**, con la
presión de cavidad intacta (0.0124). No es un bug: al cambiar la mecánica discreta de la pared se cambia
el problema. *La convergencia de malla no es monótona en este modelo y no se puede asumir.*

*Consecuencia*: bajar de ~9 % exige **P2 (o mayor orden) en la interfaz**, no mejor frontera — es el
escalón 3 de D10. Y por D20 esa mejora *no es la prioridad*: con la atenuación acotada a 2.4–3.3× y sin
datos bajo 270 kHz, bajar el error de discretización de 9 % a 2 % no cambia ninguna conclusión.

## D22 — Una superficie de presión nula NO es una salida de energía (y la ventana correcta es resistiva)

Intenté modelar la ventana redonda como parche `p = 0` para abrir la cavidad (D18). **Físicamente no
puede funcionar, y la razón es básica**: la potencia que atraviesa una superficie es `p·(v·n)`, así que
si `p = 0` la potencia es exactamente cero. Medido:

- `P_ventana/P_inc` real = *+3.1e-3* contra reactiva *−2.1e-2* (×6.7), y en `r = a` la real es
  *exactamente 0* (27 nodos con `p = 0.0` literal).
- La dispersión de la transmisión con ±5 % de frecuencia apenas baja: 17.8 % (cerrada) → 12.9 % (10°).
- Y la transmisión *decrece* con el tamaño de la ventana (1.006 / 0.937 / 0.830 para 5°/10°/20°), al
  revés de lo esperado: `p = 0` es una **pared blanda** (`R = −1`) que desintoniza, no un sumidero.

Variante de diagnóstico `window_model="open"` (se quitan los elementos sólidos del arco): ésa sí
transporta potencia real y crece monótona con el área (0.21 % → 1.12 % → 8.84 % de `P_inc` para
5°/10°/20°), pero *no es la cóclea*: la ventana redonda es una membrana, no un agujero.

**Conclusión y corrección**: la ventana tiene que ser una **impedancia resistiva**,
`p = Z_r·(v·n)`, que sí transporta potencia real. La forma débil es el mismo término de la ABC con
`k → ω·ρ/Z_r`; con `Z_r = ρc` se reproduce la ABC exactamente, lo que sirve de test interno. `Z_r` es
un **supuesto declarado**: Sim-Datos no tiene impedancia citable de la ventana redonda (D20), así que se
barre (0.5 / 1 / 2 × ρc) y se reporta la sensibilidad en vez de fijar un valor con aire de medido.

*Consecuencia*: el criterio de <10 % de dispersión sigue `xfail` con el número (12.9 % con 10° en el
modelo `p=0`); lo que falta no es malla ni frontera, es **el modelo de la ventana**.

*Regla que queda de esta seguidilla* (es la quinta condición que escribí y la física desmiente): **todo
criterio de aceptación tiene que nombrar la magnitud física que mide y su unidad, y hay que comprobar
que esa magnitud pueda ser distinta de cero en la configuración elegida.** El balance de energía "entra =
sale" se cumplía con todos los términos en 1e-14 y no probaba nada, porque la salida era nula por
construcción. Un criterio que se satisface por vacío es peor que no tener criterio.

## D23 — La ventana resistiva funciona; el problema es dónde está, no cómo se modela

`window_model="impedance"` (`p = Z_r·v_n`, `Z_r` por defecto `ρ_f·c_f`). Medido (a=20 mm, t=2 mm,
ppw=20, 100 kHz, ventana 10°):

- **Transporta potencia real**: `P_ventana/P_inc = +2.101e−2` por el borde (contra +3.08e−3 del parche
  `p=0`: *5.5× más*), con real/reactiva = 18.4 y balance de energía al **1.55 %** contra el flujo por el
  contorno exterior. Ese criterio —el que antes se cumplía por vacío— ahora se cumple con contenido.
- **Identidad con la ABC**: con `Z_r = ρ_f·c_f` la matriz de borde es *bit-idéntica* a la de ABC-1
  (`max|B_imp − B_abc1|/max|B| = 0.00e+00`, `t_interior` 0.972180 vs 0.972180). La condición de
  impedancia contiene a la ABC como caso particular, y eso quedó como test.
- **Contra ABC-2 la diferencia es 2.20 %**, y no es un error de factor: BGT-2 agrega el término de
  curvatura `1/(2ka) ≈ 6 %` a `ka ≈ 8`. Son dos condiciones distintas, no dos implementaciones de la
  misma.
- **Amortigua de verdad**: el pico de resonancia a 98 kHz baja de 1.443 (cerrada) a **1.207**
  (−16.4 %) transportando el 2.1 % de `P_inc`.

**Pero no alcanza para el criterio de banda ancha, y la razón es geométrica**: la ventana está en la
cara de sombra, donde `|p|_rms` es sólo **0.377** de la incidente; absorbe 2–3 % de la potencia. La
dispersión con ±5 % de frecuencia queda en **15.5 %** a 10° (peor que el 12.9 % del parche `p=0`, que
*desintoniza sin absorber*: un número mejor por la razón equivocada). Se cumple con 20° (10.7 %) y con
20° + `Z/2` (9.6 %).

*Conclusión*: una membrana pequeña en la cara de sombra no puede ser el amortiguador de banda ancha. El
amortiguador natural es la **membrana basilar**, que cruza el conducto donde la presión sí es alta. Es
5b-2, y ahora está justificado por medición, no por fisiología de manual.

## D24 — Sensibilidad: manda la geometría, no la atenuación (y la métrica estaba contaminada)

Barrido one-at-a-time con los rangos de D20 (`scripts/sensitivity_cochlea.py`, 85 puntos + CSV).
Caso base: `t_interior` = 0.972180, idéntico al solver (Δ 0).

**Ranking por rango relativo `(máx−mín)/base`:**

1. **espesor del hueso** (1–7 mm): 0.6635 … 1.3649 → *72.15 %*
2. **`c_p`** (1752–3500 m/s): 0.9712 … 1.1276 → *16.09 %*
3. atenuación (sin pérdidas / White / Pichardo / Pinton): 0.9191 … 0.9722 → *5.46 %*
4. `Z_r` de la ventana (0.5–2 ρc): 3.18 %
5. `c_s` (1360–1640): 2.17 %
6. `ρ` (1868–1950): 1.19 %

Controles: ppw 20→40 mueve 1.36 %; ABC-1→BGT-2 mueve −0.73 %. O sea que `Z_r`, `c_s` y `ρ` están *al
nivel del piso numérico* y no vale calibralos antes de subir la discretización (D21).

**Dos correcciones que salen de acá:**

- **`t_interior` no es monótona en el espesor** (1.10 / 0.95 / 0.97 / 1.06 / 1.18 / 1.27 / **1.36** a
  4 mm / 1.32 / 1.19 / 1.04 / 0.90 / 0.77 / 0.66 de 1 a 7 mm). No es malla (a ppw=30 el patrón se repite:
  0.9805 / 1.2938 / 0.8533 a 2/4/6 mm) ni bug (con α=0 reproduce el solver al dígito): es la **cavidad
  resonante** de D18 otra vez. La magnitud que sí se comporta como transmisión es el **armónico n=0**
  (modo de respiración): 0.3866 → 0.0761, un factor *5.1* que acompaña al `|T|` de placa (0.933 → 0.437).
  → *Regla*: junto a `t_interior` hay que reportar siempre el armónico n=0. `t_interior` solo se puede
  "mejorar" desintonizando, y eso ya nos pasó con `p=0` (D23).
- **A 100 kHz la atenuación NO es la prioridad**: mueve 5.5 %, contra 72 % del espesor y 16 % de `c_p`.
  Mi afirmación de que "la atenuación domina" vale para 1 MHz con pared de 6 mm, no para esta geometría.
  A esta frecuencia lo que manda es *cuánto hueso hay y a qué velocidad*, no cuánto absorbe.

**`c_p` y `c_s` no son independientes**: `c_s = 1900` con `c_p = 1752` da ν > 0.5 (imposible). El barrido de
`c_p` corrió con `c_s = 1500` fijo y el de `c_s` a `c_p = 3500`; quedó declarado en el CSV.

**Qué pedir primero, en orden**: (1) `c_p` de la cápsula ótica (o del temporal) en nuestra banda — es el
único material que mueve >10 %, y nuestro 3500 está en el extremo alto del rango medido, así que el
modelo es *conservador por ~16 %*; (2) espesor real de la cápsula (domina con 72 %); (3) atenuación a
≤100 kHz, que sigue desconocida pero pesa poco acá.

*Advertencia de método*: el 72 % del espesor es en buena parte resonancia de cavidad, y la membrana
basilar (5b-2) es justamente lo que rompe ese Q. **Hay que re-correr este barrido después de 5b-2**; el
ranking no es definitivo. -> *Hecho*: el barrido con membrana esta en **D30** y este ranking queda superado (el espesor baja a 53 %, `c_s` sube a 8 % y `Z_r` se cae a 0.3 %).

---

## D25 — Focalización con transductores reales: manda la cantidad de elementos, no el `pitch/λ`

15 configuraciones medidas con `scripts/focusing_study.py` (apertura 2D sobre baffle rígido; el
supuesto 2D está en el docstring y en la columna `notas` de las 15 filas) → `out/focusing_study.csv`
(27 columnas) y 15 figuras. FWHM lateral en el plano del pico axial real de cada configuración;
correlación de campo complejo contra la apertura continua del mismo diámetro.

**1. Los cuatro arreglos documentados alcanzan la apertura continua dentro del 2 %:**

| configuración | f (kHz) | N | pitch/λ | FWHM (mm) | vs continua | corr |
|---|---|---|---|---|---|---|
| `aire_40k_linea` | 40 | 16 | 0.267 | 26.64 | +1.3 % | 0.9999 |
| `transcraneal_250k` | 250 | 256 | 0.195 | 3.37 | +0.4 % | 1.0000 |
| `transcraneal_500k` | 500 | 256 | 0.391 | 1.71 | +1.6 % | 0.9992 |
| `diseno_1m` | 1000 | 256 | 0.521 | 1.15 | +1.8 % | 0.9991 |

**2. El ancla de 4B** (D = 60 mm = 4 λ, F = 40 mm, f/# = 0.67, 100 kHz, agua): λ/4 alcanza el 2 % y
λ/2 el 7 % — no hace falta llegar a λ/8 (ver D26, la tabla anterior estaba mal).

**3. La que no llega**: `langevin_40k_foco100` (N = 1) → +13.4 %, corr 0.989. Un pistón único no
tiene perfil de fase: no puede enfocar, sólo iluminar.

**4. A 40 kHz no hay foco que comparar**: FWHM 26.6 mm contra una cóclea de ~9 mm. Ahí se mide
iluminación, no focalización (D8).

**5. El mecanismo** (16 elementos, diámetro fijo, pitch/λ de 0.267 a 6.667): la correlación cae
0.9999 → 0.1844 y aparecen lóbulos de grataje en −36.9° (250 kHz), −17.5° (500 kHz) y −8.6° (1 MHz),
los tres donde predice `arcsin(λ/pitch)`. Con 256 elementos y λ/1.9 la correlación sigue en 0.9991:
**el error de cuantización del perfil de fase depende de N, no de `pitch/λ` solo**.

**6. El caso literal del criterio de grataje no existe**: pitch = λ con steer = 0 predice el lóbulo
en ±90°, justo donde el factor de elemento `sinc(w·sinθ/λ)` tiene un cero exacto (medido −75.7 dB).
Se verificó donde el lóbulo es interior: pitch = λ, steer −30°, N = 32 → predicho +30.00°, medido
+29.96°.

---

## D26 — La tabla 4B tenía los elementos al doble de ancho (bug mío, no del solver)

Al contrastar D25 con mi verificación independiente aparecieron dos discrepancias: λ/2 daba 10.71 mm
contra 9.16 mm, y pitch = λ 23.4 mm contra 10.96 mm. El diagnóstico, en orden, porque el orden es la
lección:

1. **Cuadratura**: `quad_order` 16 → 128 y mi densidad 601 → 4801 puntos no movieron *nada* en
   ninguno de los dos. Descartada la integración.
2. **Un solo elemento**, mismo punto: `aperture.pressure()` y Gauss-200 coincidían a `1.000000` de
   cociente y `0.00°` de fase. Kernel, normalización y ancho del elemento, correctos.
3. **Mi trapecio** parametrizaba `y = y_c + width·sin(t)` en vez de `y_c + (width/2)·sin(t)`: cada
   elemento de 7.5 mm se integraba como uno de 15 mm, solapado con los vecinos. De ahí la amplitud
   2.43× y los lóbulos más anchos.

Corregido, mi trapecio da `1.045312e6` contra `1.045309e6` del paquete, y λ/2 = *9.157 mm* contra
*9.16*. **La tabla de 4B venía del trapecio, así que estaba mal**; se reemplazó por las cifras del
estudio y el criterio de aceptación cambia: `pitch ≤ λ/4` para el 2 %, no λ/8.

**Lo que esto cierra**: el estudio había detectado una incoherencia interna en la tabla vieja (λ/2
con corr 0.998 *y* +27 % de FWHM a la vez). Era un síntoma del bug, no de la métrica: con +6.1 % y
0.9977 la tabla es consistente.

**Lo que no cubre**: sigue siendo 2D. Un arreglo hemisférico real no es un corte lineal con la misma
ley de fase; esto vale como cota de comparación, no como predicción de un dispositivo 3D.

---

## D27 — El dashboard expuso dos errores de ejes (y un colapso de layout que no era del solver)

Pedido: una forma de ver lo que hace el simulador, y después poder agregar y mover fuentes. Se
implementó `ultrasound_sim.dashboard`: servidor HTTP de la biblioteca estándar (sin dependencias
nuevas), Plotly servido local, y un motor vectorizado. Decisión de arquitectura: **una fuente es un
dato, no un parámetro** (`{"y", "width", "amp", "phase"}`), así que mover una fuente es cambiar su
`y` y el arreglo enfocado no es más que un generador de esa lista.

**Motor**: los nodos de Gauss de todas las fuentes se concatenan y se evalúan por fila. Contra
`aperture.pressure` el error es *9.6e-16* (Gauss de 8 nodos alcanza porque el integrando es analítico
en la malla). Costo: 60×60 en 117–477 ms según N; 120×120 en 1.0–2.0 s. Eso es lo que hace que los
sliders se sientan vivos.

**Los errores que encontró, en orden de aparición**:

1. **Axial y lateral invertidos en `metrics`** (Python). `p` es `(len(ys), len(xs))`, o sea
   *filas = y*, así que `intensidad[iy, :]` varía **x** y es el perfil *axial*, no el lateral. Estaban
   cruzados, y no se notaba porque en el test las dos longitudes coincidían por casualidad (40 y 40).
2. **FWHM en metros publicado como mm**: `focus.lateral_fwhm` devuelve metros; la clave decía `_mm`
   sin multiplicar por 1e3. El valor "0.025" era 25 mm.
3. **El mismo cruce de ejes, otra vez, en el JavaScript**: `p2_db[i][ix]` etiquetado como axial. Lo
   encontró la revisión visual del render, no un test.
4. **El índice del pico por comparación de floats**: buscar `peak_x_mm` en `x_mm` con tolerancia 1e-6
   fallaba por redondeo y caía a la columna 0, dando −6.7 dB en vez de 0 dB. Se reemplazó por un
   `argmax` sobre la matriz.
5. **Lóbulo lateral medido contra el flanco del lóbulo principal** (daba ≈ −3 dB en vez de −13 dB):
   ahora es el máximo *después del primer nulo*, y da −12.1 dB contra −13.3 dB teórico para N=32.

**Lo que no era del solver**: Plotly abortaba con *"Something went wrong with axis scaling"*. La causa
era el CSS: en vista angosta el `div` del mapa quedaba en *741×2 px* (`flex:1` con `min-height:0`
dentro de un contenedor de alto automático), y Plotly no puede escalar un eje de 2 px. Se corrigió el
layout, no el cálculo. Relacionado: el aspecto 1:1 no se puede pedir con `scaleanchor` + rangos fijos
(es lo que dispara el error); se calculan los rangos a partir del área de trazado real, y como la barra
de color se come ~65 px, el segundo paso mide el `px/mm` efectivo y corrige.

**Verificación**: `tests/test_dashboard.py`, 12 tests: el motor contra `aperture.pressure`, la ley de
fase contra `arrays.LinearArray.focused`, pico/grataje/FWHM contra lo teórico (`0.886·λ·f/#` → 6.73 mm
contra 6.65 mm), y el servidor real en un puerto efímero por HTTP (assets, POST, GET, 400 y 404).

**Arrastrar fuentes y guardar configuraciones (segunda pasada del dashboard)**: las fuentes viven
sobre el baffle (`x = 0`), así que arrastrarlas es 1-D: el asa mueve `y`, no hay grado de libertad en
`x`. El overlay de asas es propio (no de Plotly) y se *reposiciona* en cada render en vez de recrearse —
si se recreara, el `pointer capture` moriría a mitad del jalón y el arrastre se cortaría en el primer
recálculo. `setPointerCapture` va en un `try` (con pointers sintéticos tira) y los listeners de
`pointermove`/`pointerup` viven en `window`, así el asa sigue al puntero aunque se salga del cuadrado.
Durante el jalón la malla baja a 80×80 (~0.4 s) y al soltar vuelve a la elegida: el campo sigue al
dedo sin pagar 2 s por paso.

Configuración: `instantanea()` produce
`{formato: 'ultrasound-sim/config', version, material, frequency_khz, modo, n_elements, aperture_mm,
focus_x_mm, focus_y_mm, grid, sources}`; se descarga como JSON, se carga por `<input type=file>` y se
autoguarda en `localStorage` cada 3 s, así recargar la página no pierde el trabajo. `aplicar()` valida
el campo `formato` y avisa en la UI si el archivo no es del simulador. El test
`test_la_configuracion_que_escribe_el_cliente_es_la_que_el_servidor_entiende` manda ese mismo formato
al endpoint, de modo que renombrar un campo en el JS rompe el test.

Verificado en un navegador real, no sólo con tests: arrastrar el asa 0 cuarenta píxeles llevó la fuente
de −29.06 a −37.34 mm, actualizó la tabla y recalculó el campo (pico 75.8 → 71.1 mm, FWHM 7.06 →
6.93 mm); el doble clic agregó la fuente 33 en la altura del clic; cargar un JSON por el input cambió
material, frecuencia, elementos y malla, y recalculó; y tras recargar la página la configuración
seguía ahí.

**Deuda declarada**: no hay todavía "borrar todas" ni deshacer un arrastre; el rango de frecuencia no
está limitado por el techo de P1 (D8), así que se puede pedir 1 MHz en 2D y el resultado es correcto
pero caro de interpretar.

---

## D28 — El baffle no cambia el campo delantero; lo que niega es la radiación trasera

El dashboard tenía un solo modelo de radiación: la integral de Rayleigh sobre baffle rígido infinito,
`p = 2·i·ω·ρ·Σ ∫ v·G`. Se agregó el mismo problema resuelto **sin baffle**
(`dashboard/bem_field.py`): la apertura es una lámina delgada de ancho finito con fluido a los dos
lados, representada por una capa simple. Como la lámina es simétrica, la densidad queda determinada
por el salto de velocidad —`σ = −2·i·ω·ρ·v`— y **no hay sistema lineal que invertir**: no hay cuerpo
cerrado ni dominio interior, así que Burton-Miller y CHIEF (que existen para cerrar la ecuación en
superficies cerradas) no aplican. Todo el trabajo numérico es la cuadratura con la singularidad
logarítmica de `H₀`, resuelta por refinamiento graduado.

**El resultado, que no era el esperado**: en `x > 0` los dos modelos dan **el mismo módulo con 180° de
fase**. Es una identidad, no una coincidencia numérica:

- `p_libre + p_baffle` sobre una malla en `x > 0`: máximo `6.5e-12` relativo → `p_libre = −p_baffle`
  exacto.
- Reparto de potencia con la lámina: atrás/adelante = *0.9871* (la lámina radia igual a los dos lados,
  como debe).

O sea que **el supuesto de baffle no sesga nada de lo que el dashboard venía mostrando delante**. Lo
que sí hace es negar el semiespacio trasero, donde la lámina libre radia con la misma fuerza: el
máximo global del modelo libre puede caer *detrás* del plano (medido: −75.6 mm en la configuración por
defecto). Eso es lo que la UI ahora deja ver y avisa en el panel.

**Costo**: 80×80 con 32 paneles en *1.30 s*, contra ~0.35 s del modelo de Rayleigh (la cuadratura
graduada es cara). A 120×120 son ~3 s, así que el modelo libre no baja a 80×80 durante el arrastre de
asas: el jalón se siente más lento con este modelo y está documentado.

**Verificación del módulo** (`tests/test_bem_field.py`, 31 tests):

| criterio | medido |
|---|---|
| panel de 1e-4·λ contra el cilindro pulsante exacto | **2.155e-3** relativo (fuente lineal: 2.151e-3) |
| residuo de `∇²p + k²p = 0` a h = λ/40 | **1.9e-3 … 2.7e-3** |
| decaimiento `|p|·√r` entre 0.3 y 3 m | dispersión **8.64e-7** |
| reciprocidad (3 pares) | **1.2e-16** |
| contraste con el modelo con baffle en el foco | **ΔdB = 0.0000**, desfase π |
| 19 entradas inválidas | `ValueError` con mensaje; salida `complex128` finita siempre |

**Dos correcciones de convención que salieron de acá**:

1. La fórmula que escribí en el enunciado de la tarea, `p = −ρc·v·H₀(kr)/H₀′(ka)`, **tenía un −i de
   menos**. La correcta, derivada de `u = ∇p/(i·ω·ρ)`, es `p = i·ρ·c·v·H₀(kr)/H₀′(ka)`, y es la que
   verifica el test contra el cilindro. El error era mío, en el pedido, no del módulo.
2. `aperture.py` y `dashboard/field.py` tienen un **−1 global** frente a `u = ∇p/(i·ω·ρ)`: implican
   `p = −ρc·v` para un pistón infinito con baffle. Es un cambio de fase benigno —todas las métricas
   del repo son de módulo— pero obliga a comparar módulos cuando se cruza con `bem_field`. Queda
   anotado y **sin arreglar**: cambiarlo rompería las comparaciones bit a bit con los números ya
   publicados, y no cambia ninguna conclusión.

**Lo que esto NO es**: no es un BEM de contorno cerrado resolviendo un sistema. Para contrastar contra
un obstáculo real (cápsula ótica) haría falta usar `bem/helmholtz2d.py` con el campo incidente de las
aperturas: eso queda como el siguiente paso natural, y es otra cosa que la capa simple de acá. El
docstring del módulo lo dice explícitamente y `equation` sólo acepta `single_layer`: pedir
`burton_miller` o `chief` levanta `ValueError` explicando que regularizan resonancias interiores de
cuerpos cerrados y acá `[p] = 0` por simetría.

---

## D29 — La membrana basilar rompe la resonancia de cavidad: Q de 287 a menos de 5

Fase 5b-2. La membrana es un **segmento con impedancia distribuida** sobre un diámetro del fluido
interior, `Z_m = R + i·ω·m + k/(i·ω)`, ensamblada con el mismo bloque de borde que la ventana
resistiva de D23 —no hay operador nuevo, y con `Z = ρ_f·c_f` el bloque es **bit-idéntico** a la ABC-1
(`max|ΔB|/max|B| = 0.00e+00`). La membrana divide el fluido en scala vestibuli y tympani y, sobre todo,
le da al problema el camino de pérdida que D18 pedía.

**Criterio 1 — amortiguamiento (pedía Q_antes/Q_después > 2)**:

- Sin membrana: pico en *98.00 kHz* = 1.6698 con FWHM *0.341 kHz* ⇒ **Q = 287**.
- Con membrana: máximo *0.6731* en 97.0 kHz; en [88, 108] kHz no baja a A/√2 (mínimo 0.4787 contra
  nivel 0.4760) ⇒ FWHM > 19.6 kHz ⇒ **Q < 5**. Relación **> 57×**.
- Control de malla a ppw=20: el pico da 1.735 contra 1.670 de esta malla (+3.9 %) ⇒ el pico es físico,
  no ruido de discretización. **Pero el Q no está convergido en malla**: D21 ya midió que este modelo
  no converge monótono, así que el 287 es un número de esta malla, no una propiedad del modelo.

**Criterio 2 — monotonía con el espesor** (`t_interior`, a = 20 mm, 100 kHz, ppw = 20). La columna
"sin membrana" reproduce D24 **dígito a dígito** (test con tolerancia 5e-5), que es la no-regresión
de todo lo anterior:

| t [mm] | D24 / sin membrana | con membrana | n=0 sin | n=0 con |
|---|---|---|---|---|
| 1.0 | 1.0996 | 0.6040 | 0.3866 | 0.1813 |
| 2.0 | 0.9722 | 0.5319 | 0.2547 | 0.1791 |
| 3.0 | 1.1839 | 0.4571 | 0.1936 | 0.1443 |
| 4.0 | 1.3649 | 0.3876 | 0.1657 | 0.1261 |
| 5.0 | 1.1920 | 0.3520 | 0.1220 | 0.1256 |
| 6.0 | 0.8968 | 0.3276 | 0.0831 | 0.1214 |
| 7.0 | 0.6635 | 0.3247 | 0.0761 | 0.1169 |

Sin membrana: rango 2.06× y **dos** cambios de signo en Δt (el lóbulo de resonancia de D24). Con
membrana: rango 1.88×, **un** cambio de signo, *12 de 13 tramos decrecientes* y el último **+0.97 %**.
Ese repunte es reproducible a ppw=30 (+1.1 %) y significa que **con la pared gruesa el nivel ya no lo
fija la pared sino la membrana**: la curva satura en ~0.32. No se fuerza la monotonía estricta; se
documenta el tramo.

**Criterio 5 — el armónico n=0 sigue acompañando** (regla de D24): sin membrana 0.3866 → 0.0761
(*5.08×*), con membrana 0.1813 → 0.1169 (*1.55×*). La membrana reduce la dependencia del espesor casi
a la mitad también en esa métrica.

**Criterios 3 y 4 — límites y energía**, con los números finos: `Z = 10⁶·ρc` contra el rígido impuesto
da Δ`t_interior` = *−1e-4 %* y `max|Δp|/max|p|` = *2e-4 %*; `Z = 10⁻⁶·ρc` contra `p = 0`, lo mismo; la
convergencia con el factor de impedancia es monótona (2.0 %, 0.02 %, 2e-4 %). El balance cierra al
*0.46 %* con la membrana (`−P_total(R_out)` = 0.7281 contra `P_memb` = 0.7327), 0.39 % con Z/2, 0.47 %
con 2Z y 0.11 % en el caso rígido con `P_abs = 0`. Con derivada radial de 2º orden el balance daba
1.97 %: por eso el flujo se mide con 5 anillos. La membrana absorbe *P_memb/P_inc = 0.733*. A 90°
(perpendicular a la incidencia) absorbe más que a 30°: la orientación no es cosmética. Se rechaza
`Re(Z) ≤ 0` en vez de devolver un cero silencioso (D22).

**Sensibilidad a las propiedades de la membrana** (que son VERIFY, pero conviene saber cuánto pesan):
a 100 kHz la reactancia total es el *0.91 %* de `R`, así que el comportamiento lo fija la resistencia;
multiplicar la masa por 10 mueve `t_interior` un *1.75 %* y multiplicar la rigidez por 10, un *0.02 %*.
Conclusión: mientras `|X| ≪ R`, el resultado no depende del ballpark de masa y rigidez —que es
justamente lo que no tenemos medido—, y eso es lo que hace que el criterio 1 sea defendible.

**Lo que sigue siendo VERIFY**: todas las propiedades de la membrana. `BASILAR_MEMBRANE_BALLPARK` usa
m = 0.022 kg/m² (1100 kg/m³ × 20 µm), R = ρ_f·c_f y k = 1e8 Pa/m, cada uno con su rango publicado y su
`provenance`, y **ninguno se usa como validación**: los criterios se sostienen por límites, balance de
energía y no-regresión contra D24, que son relativos.

**Bug reportado y no arreglado**: el redondeo `n_theta + n_theta % 4` no garantiza múltiplo de 4 (a
ppw=20 da 502). Se corrigió la descripción y se agregó redondeo a par cuando hay membrana; la fórmula
formula original quedó intacta para no mover los números publicados.

---

## D30 — Con la membrana puesta el ranking cambia: el espesor pesa menos, `c_s` deja de ser ruido

D24 lo pedía explícito y acá está: el mismo barrido one-at-a-time, con la membrana de D29
(`scripts/sensitivity_cochlea.py --membrane`, 85 puntos, *171.9 s*,
`out/sensitivity_membrane.csv`). El script reproduce el caso base sin membrana **al dígito**
(`t_interior = 0.97218`, la referencia de `solve_cross_section`), así que las dos columnas son
comparables.

| # | parámetro | rango rel. **con membrana** | D24 (sin membrana) | qué pasó |
|---|---|---|---|---|
| 1 | espesor | **53.11 %** | 72.15 % | sigue dominando, pero 19 puntos menos |
| 2 | `c_p` | **16.06 %** | 16.09 % | *idéntico*: no dependía de la resonancia |
| 3 | atenuación | **9.60 %** | 5.46 % | sube: con la membrana el camino importa más |
| 4 | `c_s` | **8.22 %** | 2.17 % | **estaba en el piso numérico y deja de estarlo** |
| 5 | `ρ` | 1.03 % | 1.19 % | sigue al nivel del piso (0.93 %) |
| 6 | `Z_r` | **0.28 %** | 3.18 % | colapsa: la pérdida ya no pasa por la ventana |

El piso numérico medido (ppw 20 → 40) es *0.93 %*: por debajo de eso el ranking no distingue.

**Las tres conclusiones que cambian de verdad**:

1. **El 72 % del espesor era, en buena parte, la resonancia.** Sin cavidad resonante baja a 53 %: sigue
   siendo el parámetro más influyente —y por lo tanto lo primero que hay que medir—, pero dejó de ser
   "casi todo".
2. **`c_s` sale del piso numérico.** Con 2.17 % estaba al nivel del error de malla, y D24 decía
   explícitamente "no vale calibrarlos antes de subir la discretización". Con la membrana da *8.22 %*,
   nueve veces el piso: ahora importa y ahora hay que reportarlo con su incertidumbre.
3. **`Z_r` se vuelve irrelevante (3.18 % → 0.28 %).** Es lo más limpio de todo: la ventana resistiva era
   el amortiguador en el modelo cerrado de D23, y con la membrana el camino de pérdida se mudó a ella.
   Confirma que la ventana **no** es la palanca del modelo.

**Lo que NO cambió, y hay que decirlo**: `t_interior` sigue sin ser estrictamente monótona (el último
tramo sube +0.97 %, D29) y el armónico `n = 0` tampoco (0.1813 → 0.1959 a 1.5 mm, +8.0 % de ese tramo,
dentro del piso). La membrana saca la resonancia grande, no toda la estructura fina.

**Qué le pide esto a la parte experimental**: el orden de D24 se mantiene en el primer puesto
(espesor) y sube en el tercero (atenuación), y *agrega* `c_s` a la lista de parámetros que hay que
medir — antes se podía ignorar por estar dentro del ruido numérico. `ρ` y `Z_r` se pueden seguir
dejando en el ballpark.














