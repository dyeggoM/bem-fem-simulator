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
ranking no es definitivo.













