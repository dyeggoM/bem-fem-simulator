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

**Ajuste de rendimiento que salió de acá**: al exponer el corte en la UI quedó en evidencia que
`solve_cross_section` usaba el ordenamiento **COLAMD** que trae scipy por defecto para el LU disperso.
Un A/B en el mismo proceso (mismo sistema, 45 240 gdl, 345 556 nnz) dio *MMD_AT_PLUS_A* ≈ 1.0 s contra
*COLAMD* 1.0–29.4 s según la corrida, con la solución **idéntica a 8 decimales** (`t_interior =
0.48069574` en los cuatro casos). La librería pasó a MMD; el barrido de D30 ya lo usaba, y ahora los
dos coinciden a *0.00e+00*.

---

## D31 — Los tres parámetros que faltan medir, editables en la UI

D30 dejó una lista corta de lo que la parte experimental tiene que medir: espesor de la cápsula (53 %),
atenuación bajo 270 kHz (10 %) y `c_s` (8 %, nuevo). Los tres estaban hardcodeados en el ballpark. Ahora
se editan en el panel **cóclea · corte transversal** del dashboard, cada uno con un tooltip que dice
qué es, su rango publicado, cuánto pesa y de qué ADR sale.

- **espesor** 0.5–8 mm (arranca en 2 mm): es el que más pesa; el rango medido es 1–7 mm.
- **c_s** 1000–2500 m/s (arranca en 1900): White 2006 da 1500 ± 140 m/s; 1900 es el ballpark del plan y
  está en el extremo alto. Tiene que ser menor que `c_p`.
- **α** 0–40 Np/m (arranca en 0): las tres leyes de D19 dan 3.78 / 6.24 / 15.31 Np/m a 100 kHz.

`solve_cross_section` no sabía de atenuación en el sólido (el script de sensibilidad tiene su propio
ensamblado con `D` compleja desde antes). Se agregó `alpha` [Np/m] a la librería por el camino de
módulos complejos (`p1_elastic_local_damped`, `c~ = c/(1 + i·α·c/ω)`), y el criterio de aceptación es
fuerte: **la librería tiene que reproducir el ensamblado independiente del script**. Ahora da
*0.00e+00* de diferencia relativa (antes 1.39e-14, antes de igualar el ordenamiento del LU), y
`alpha = 0` es **bit-idéntico** al camino sin el parámetro.

Números de referencia del panel (t = 2 mm, c_s = 1900 m/s): α = 0 → `t_interior = 0.5319` (el de D29),
α = 15.31 → `0.4808` (el de D30 con la ley de Pinton), y la membrana absorbe 72.9 % y 66.8 %
respectivamente. Con t = 6 mm y c_s = 1500: `0.2963`.

**Qué NO se expone y por qué**: `c_p = 3500 m/s` y `ρ = 1950 kg/m³` siguen fijos. `c_p` pesa 16 % y
también es VERIFY (D20), así que es un candidato obvio a ser el cuarto campo; queda anotado. La
membrana usa `BASILAR_MEMBRANE_BALLPARK` (VERIFY) y la app lo dice en el pie.

---

## D32 — Núcleo 3D: converge como P1 tiene que converger, y el techo es el costo

Fase 6a. Tres módulos nuevos: `meshing3d.py` (malla de cajas: cada celda en 6 tetraedros con
orientación positiva verificada), `fem/helmholtz3d.py` (P1 en tetraedros, ABC-1 y BGT-2 en caras
planas, campo total y dispersado, fuente puntual) y `aperture3d.py` (integral de Rayleigh en 3D con
pistón circular o rectangular y foco). No hay geometría de cóclea ni elasticidad 3D todavía.

**Lo medido** (`tests/test_helmholtz3d.py`, 5 tests; `tests/test_aperture3d.py`, 14):

| criterio | medido |
|---|---|
| malla: volúmenes positivos y suma = volumen de la caja | exacto a 1e-12; 176 caras de frontera en 4×2×6 |
| onda plana con ABC (caja 2 λ × 1 λ × 1 λ) | 20.50 % a ppw 8 → **3.53 % a ppw 20** |
| orden de convergencia | **1.88 / 1.94 / 1.97** (P1, como corresponde) |
| ABC contra pared rígida | **9.57 %** contra **100.0 %** |
| fase del campo: `d(arg p)/dx` contra `k` | 4.192e2 contra 4.189e2 (0.07 %) |
| fuente puntual contra el Green 3D | amplitud **0.9903** del teórico; L2 8.4 % en caja de 1.1 λ |
| pistón 3D en el eje contra Rayleigh-Sommerfeld | 0.0000 %; último máximo en `a²/λ − λ/4` exacto |
| `|p|` en el pico axial | `2·ρ·c·v₀` exacto |

**El costo es el techo, no la física**: 2 λ × 1 λ × 1 λ a ppw 20 son 18 081 nodos y 96 000 tets →
5.3 s por solve. Un cubo de 3 λ a ppw 12 (50 653 nodos) no terminó en 10 minutos. O sea: en 3D el
presupuesto se gasta en la LU dispersa, y por eso el dominio útil con este ensamblado es de pocas
longitudes de onda — exactamente lo que el plan anticipaba (D8/D10) cuando decía que a 1 MHz hay que
pasar a FEniCSx+MPI, a k-Wave o resolver el exterior por BEM.

**Dos errores que fueron míos, no del módulo**, y quedan escritos porque el primero es fácil de
repetir:

1. **Factor ρ**: la ecuación se ensambla dividida por `element_rho`, que por defecto vale **1**. Sin
   pasar `element_rho=W.rho` la fuente puntual da `1/ρ` de la amplitud (medido: 0.0010 del teórico).
   El módulo está bien; el que se olvidó la densidad fui yo.
2. **Signo de `i` en el Green**: con la convención del repo y un RHS `−i·ω·Q`, la solución es
   `p = −i·ω·ρ·Q·e^{ikr}/(4πr)`. Con el signo opuesto la amplitud coincide al 1 % y el desfase es un
   **+178° constante**: un error de fase perfectamente plausible si no se mira el número complejo.

**Estratificado en 3D — validado** (`test_estratificado_rt_contra_la_matriz_de_transferencia`): una
losa de un fluido de prueba (c = 2000 m/s, ρ = 1200) en agua, con contraste de impedancia real.
Contra la matriz de transferencia 1D: **R al 2.91 %** (complejo, 0.2493+0.2114j contra 0.2424+0.2179j)
y **|T| al 2.11 %** (0.9254 contra 0.9454).

Por qué se ajustan R y T en una ventana pegada a las interfaces y no se compara el campo punto a
punto: **la dispersión de P1 se acumula con la distancia**. Medido en la misma caja: el error crece
13 % → 27 % → 39 % de izquierda a derecha, y es *idéntico* en el caso homogéneo — es el mismo techo
que D10 documentó en 2D, no un problema del estratificado. La fase de `T` además depende de si se
refiere a `x = 0` o a la salida de la losa; R y |T| no, por eso son los que se exigen.

**Campo dispersado: sólo obstáculo rígido, y sin validar**. `assemble_scattered` no arma el término de
contraste voluminoso `−(k² − k_i²)·p_i/ρ`, así que con un dispersor **penetrable** devuelve `p_s = 0`
exacto (medido). Lo que sí está es el camino de obstáculo rígido (`neumann_faces` +
`dp_dn_incident = −∂p_i/∂n`), sin test todavía. Lo único verificado es el sanity: en medio homogéneo
`p_s = 0` **exacto**, que descarta términos incidentes espurios.

**Intento de cerrar esa brecha: no salió, y quedó revertido.** Se probaron dos formas del término de
contraste (`−(k²−k_b²)/ρ·p_i·v` y la exacta `∫(1/ρ)∇p_i·∇v − (k²/ρ)p_i·v`, que además cubre el
contraste de `ρ`), y después restringirlo a los elementos con contraste. Nada validó:

- Con el operador ya idéntico entre las dos formulaciones (verificado: `‖A_total − A_dispersado‖ = 0.0`
  en la misma malla), la relación `b_abc = A·p_i + b_contraste` **no se cumple**: difiere 92 % de
  `b_abc`, repartido por igual en los 2 624 nodos (no sólo en las caras ABC).
- Sobre el estratificado, donde R y T están validados al 2.91 %/2.11 %: el dispersado da
  `|R| = 0.96` sin restringir (vs 0.33 analítico, y con los dos signos igual) y `|R| = 76` restringido
  a la zona de contraste. `dispersado + incidente` contra campo total: 95 %–134 % de diferencia.

Diagnóstico que deja: el término **no** es el problema — la ecuación del dispersado tiene que usar el
operador del **fondo** (`k_b`, `ρ_b`) en todo el dominio, no el `k` local del elemento, y entonces el
contraste entra entero en el lado derecho. El código actual ensambla el operador con `element_k`
local (correcto para el campo total, no para la diferencia), y eso es incompatible con cualquier
forma del término que se agregue al lado derecho. Rehacerlo pide: operador con `k_b`, RHS con
contraste, y una referencia de validación externa (serie de esfera) antes de volver a mirar números.

**Lo que falta para llamar a esto "cóclea 3D"**: el término de contraste para dispersores penetrables
(y su validación contra algún caso con solución conocida), elasticidad 3D, y la geometría de CT. Cada
uno es una fase, no un detalle.

## D33 — La UI es cliente de los solvers, no un segundo simulador

**Decisión.** El dashboard pasa de "explorador de campo de una apertura" a **banco de simulación**
con tres piezas: una *escena* (regiones + materiales + fuentes + malla), un *job* (una corrida de un
solver concreto sobre esa escena) y un *resultado* sellado con el criterio contra el que se validó.
Todo número que muestre la UI sale de un módulo de la librería que tiene test. Si falta algo en la
librería, primero se agrega ahí con su criterio medible y después se enchufa.

**Por qué no seguir con el motor propio del dashboard.** El motor actual (D27) no es FEM ni BEM: es la
integral de Rayleigh sobre aperturas con baffle. Es *exacto dentro de su modelo* —coincide con
`aperture.pressure` a 9.6e-16— y por eso sirve para explorar fuentes rápido, pero no puede representar
heterogeneidades, interfaces, dispersores ni elasticidad. Todo lo que se midió en las fases 0–6 vive
en módulos que la UI no toca. Dejar que la UI calcule por su cuenta sería tener dos simuladores y un
solo conjunto de validaciones.

**Consecuencias medidas.**

- **Síncrono hasta 2D, asíncrono desde 3D.** Un solve 2D con 5 k elementos tarda ~1 s y el motor de
  campo actual 117 ms–2 s (60×60–120×120, D27): el `POST` síncrono alcanza para U1–U3. Un 3D de
  18 081 nodos tarda 5.3 s y no escala (un cubo de 3 λ a ppw 12 no termina en 10 min), así que U4
  obliga a cola de trabajos con progreso y cancelación, más un tope de elementos con mensaje antes
  de arrancar en vez de colgar la máquina.
- **3D se muestra por rebanadas**, no como volumen: a 18 k nodos el campo es un sólido y el trazado
  completo no aporta lectura. Los cortes axial/lateral que ya existen se reusan.
- **Escena en JSON, malla reproducible.** Guardar y cargar tiene que devolver la *misma* malla, no
  una equivalente: si cambia el conteo de elementos cambian los números y se pierde la comparación
  entre corridas.
- **Sin editor CAD.** Cajas, rectángulos, círculos y arcos cubren las fases 0–6. Un CAD es una
  promesa que hoy no se puede sostener.
- **Sello de trazabilidad.** El sello (criterio, número de referencia, margen) se calcula en la
  corrida. Si se escribe como constante en el JS, en tres meses nadie sabe qué era confiable: es
  exactamente el error que el test de `test_cochlea_damping.py` evita del lado de la librería.

**Costo asumido.** La UI no va a poder mostrar nada que la librería no sepa calcular, y eso incluye
huecos conocidos: dispersión penetrable 3D (D32), elasticidad 3D y geometría de CT. Van a aparecer
en la UI como *no disponible*, no como un resultado aproximado.

## D34 — U1: la escena 2D, y dos cosas que sólo se ven en el navegador

**Qué se entregó.** Una pestaña "escena FEM" donde se arma un apilado de capas (material + espesor,
del registro de 5 materiales con los dos fluidos cocleares), se estima el costo de malla **antes** de
resolver, y se corre `helmholtz2d.solve` de verdad. Devuelve el campo en dB, los cortes axial y
lateral, y R/T con su referencia analítica y el error. `grid_tri_mesh(xs, ys)` arma el reticulado con
las interfaces exactamente en los nodos: la malla es conforme por construcción, así que no hay que
pegar regiones ni adivinar de qué material es cada elemento.

**Criterio de aceptación cumplido** (el de la Fase 1, 2 % en amplitud). Agua/hueso/agua a ppw 12:
R 0.9018 contra 0.9037 (abs 0.0019), T 0.4323 contra 0.4281 (abs 0.0041), energía 1.0000, 11 ms.
Agua/tejido/agua: R abs 0.0060, T abs 0.0009.

**Hallazgo 1 — el criterio tiene que ser absoluto, no relativo.** Con agua/tejido R vale 0.012, y el
error *relativo* de R sale 218 % mientras el error *absoluto* es 0.006: casi todo es piso de
discretización (P1 a ppw 12 ya da ~0.6 % de error de amplitud). Juzgar por relativo habría rechazado
una escena que cumple. Por eso la UI muestra `error_R_abs` / `error_T_abs` y el sello se arma con
esos, no con los porcentajes.

**Hallazgo 2 — `[hidden]` no oculta si el CSS pone `display:flex`.** Las dos vistas (`main.vista`)
compartían el ancho y el mapa de la escena se quedaba en 0 px aunque Plotly ya tenía los datos
(17×45). Se veía en el navegador como un panel vacío sin ningún error en consola. Arreglado con
`[hidden] { display:none !important; }`. Es el mismo tipo de bug que D27 ya había encontrado con
`scaleanchor`: el render no falla, simplemente no se ve.

**Lo que parece un bug y es física.** Con agua/tejido en incidencia normal el mapa sale casi
uniforme: R = 1.2 %, así que la onda estacionaria varía 2.4 % y en una escala de 40 dB no se ve. Por
eso el caso por defecto es agua/hueso/agua (R = 0.90) y hay selector de rango (40/20/10 dB o lineal).
Con el rango en 10 dB aparecen las franjas.

**Límites declarados en la propia UI**: sólo incidencia normal (arriba y abajo las paredes son
rígidas; una onda oblicua se reflejaría ahí, y el pedido con `angulo_deg != 0` se rechaza con esa
razón en vez de devolver un campo que parece válido), sólo fluido acústico, y tope de 60 000
elementos con mensaje que dice qué bajar.

## D35 — U2: el dispersor en la UI, y tres formas de romperse que no se ven en los tests

**Qué se entregó.** Una pestaña "dispersor BEM" donde el mismo cilindro rígido en agua se resuelve
por tres caminos a la vez: BEM (exterior, ecuación seleccionable entre `standard`, Burton-Miller y
CHIEF), FEM (corona `a <= r <= R` en campo dispersado con dato Neumann y ABC-2) y la **serie de
Bessel** como referencia. Se ve el mapa del campo dispersado, las tres curvas `|p|` contra θ en
`r = 1.5a` sobrepuestas, y las tres discrepancias. La receta está copiada de `cross_validate.py` —el
script que ya tiene los números de D12— y no reinventada.

**Aceptación cumplida** (Fase 4, 2 %): ka=1 → FEM−BEM *0.97 %*, FEM−serie *0.94 %*, BEM−serie
*0.05 %*; ka=2 → *0.58 / 0.60 / 0.14 %*; ka=3 → *0.84 / 0.72 / 0.27 %*. Las tres ecuaciones dan el
mismo campo: BM sube el condicionamiento de 2.37 a 12.78 y su error contra la serie de 0.05 % a
0.24 %, que sigue siendo 8× más chico que el criterio — el precio conocido de regularizar para
poder cruzar resonancias en un barrido.

**Hallazgo 1 — `NaN` no es JSON.** El mapa tiene huecos (dentro del cilindro no hay campo) y el
`json.dumps` de Python escribe el literal `NaN`, que **no es JSON válido**: el navegador rechaza el
cuerpo entero con `Unexpected token 'N'` y la UI queda sin nada, sin un mensaje que explique por qué.
Ahora los huecos se mandan como `None` (→ `null`) y el servidor serializa con `allow_nan=False`, así
que un NaN futuro revienta en el servidor con un 500 y un mensaje en vez de morir en el navegador.

**Hallazgo 2 — `flex-basis` le gana a `height`.** La gráfica de comparación quedó en *12 px* de alto.
No era Plotly: `.corte.ancho { flex: 1 1 100% }` fija el `flex-basis` en 100 %, y en un contenedor
flex de columna con alto automático ese 100 % resuelve a 0 y pisotea el `height: 190px` de la regla
responsive. Se cambió a `flex: 1 1 auto` con `min-height`. Tercera variante del mismo tema en D27 y
D34: **el render no falla, simplemente no se ve**, y el único detector es mirar el navegador.

**Hallazgo 3 — un test con la premisa equivocada.** El primer test de "sanity física" pedía que el
campo decayera al alejarse. Falla por física correcta: el mapa muestra `|p_total|`, que lejos tiende
a la onda incidente (|p| → 1), porque el **total** es incidente + dispersado. Medido: esquina −2.9 dB
contra −3.8 dB cerca. Se reemplazó por la condición de radiación sobre `p_scattered`
(`|p_s(4a)| / |p_s(1.6a)| < 0.75`), que es lo que de verdad tiene que decaer.

**Costo y su optimización.** El mapa son 4 812 puntos (70×70) evaluados con el BEM: 11.3 s con la
cuadratura por defecto para las 64 curvas. Bajar `quad_order` a 4 para el mapa lo deja en *1.9 s* con
*0.00 %* de diferencia contra `quad_order=12` — los puntos del mapa están lejos de la superficie y el
integrando es suave. La cuadratura alta se reserva para las 64 curvas de comparación, donde 0.05 %
importa.

## D36 — U3: la cóclea en la UI, y dos cuentas que no cerraban

**Qué se entregó.** Una pestaña "cóclea" con el corte transversal completo (cápsula elástica,
perilinfa, membrana basilar, ventana redonda) resuelto por el mismo `solve_cross_section` de la
Fase 5, más de lo que tenía el panel de D31: el **mapa del campo interior** (|p| interpolado sobre
el disco a partir de los nodos de fluido), un corte por el centro, y el balance de energía. Los tres
sliders (espesor, `c_s`, α) siguen con sus pesos medidos en el tooltip (53 %, 8 %, 10 %), y se
sumaron dos interruptores: membrana y ventana.

**Aceptación cumplida**: `t_interior = 0.53187` contra 0.53187 de D30 (**desvío 0.00 %**), armónico
`n = 0` 0.17913, membrana absorbiendo 72.9 %, balance de energía 0.09 %. Sin membrana el mismo
corte da `t_interior` 0.97218 y absorción exactamente 0: el interruptor cambia la física.

**Hallazgo 1 — inventé un balance que no existe.** Sumé `power_scattered_fraction` +
`power_membrane_fraction` + `power_window_real_fraction` esperando 1 y me dio **1.85**. No estaba
mal el código: estaba mal la cuenta. `power_scattered_fraction` (1.12 en el caso base) es
`P_disp/P_inc` medido en otro radio y con otra normalización, y no se suma con la absorción. La
identidad que sí documenta `CrossSectionSolution` es **`P_total + P_membrana = 0`** cuando la
membrana es el único sumidero: da 0.00085, o sea cierra al **0.09 %**. La UI reporta esa y no la
inventada. Lección: cuando un "balance" no cierra, primero se lee la definición de cada término, no
se busca el bug en el solver.

**Hallazgo 2 — segundo tipo de objeto que revienta el JSON.** Después de `NaN` (D35) apareció
`window_impedance`, que es `complex`: `json.dumps` lanza `Object of type complex is not JSON
serializable` y el navegador no recibe nada. Ahora hay un saneador (`_json_limpio`) que convierte
complejos a `{re, im}`, numpy a tipos nativos y no finitos a `null`. Los dos casos comparten causa
—serializar un resultado científico tal cual— así que el saneo es un solo paso antes de responder,
y el test fija `json.dumps(..., allow_nan=False)` como contrato.

**Lo que se decidió no hacer.** `Q` (ancho de banda de la resonancia de la cavidad) es la métrica
que D29 usó para mostrar que la membrana rompe la resonancia (287 → <5), pero medirla exige un
barrido de frecuencia: ~9 solves de 4.3 s, con la malla reconstruida en cada punto. Eso es un
trabajo en cola (U5), no un botón de esta pantalla, así que la UI dice "necesita barrido de
frecuencia (U5)" en vez de mostrar un número que no midió.

**Nota sobre el aspecto circular.** Acá sí se usa `scaleanchor` (a diferencia de D27, donde
rompía el mapa): en D27 el `div` colapsaba a 2 px de alto y Plotly fallaba; con el alto resuelto
(271 px medidos) el anclaje funciona y el disco sale circular, que con dos circunferencias
concéntricas no es cosmético. Si en algún layout vuelve a colapsar, el síntoma va a ser
"Something went wrong with axis scaling" y la salida es sacar el anclaje, no tocar el solver.

## D37 — U5: la cola de trabajos, el barrido que reproduce D30, y el Q que no se mide

**Qué se entregó.** Una cola de trabajos en proceso (`dashboard/jobs.py`) con un solo worker,
progreso, cancelación cooperativa y retención acotada a 20 jobs; el barrido one-at-a-time de la
cóclea como trabajo de esa cola (`dashboard/sweeps.py`); y una pestaña "barridos" que lanza, sigue el
progreso, cancela y dibuja el ranking contra la referencia.

**Aceptación cumplida, y de la forma más exigente**: el barrido de la UI reproduce el ranking de D30
*casi dígito a dígito* —espesor 53.11 % (D30 53.11), c_p 16.06 (16.06), atenuación 9.60 (9.60), c_s
8.22 (8.22), ρ 1.03 (1.03), Z_r 0.28 (0.28)— con peor desvío 1.8 % y `t_interior` base 0.53187, en
89 s. La única diferencia de fondo con el script es la grilla 2D espesor × atenuación, que no cambia
el ranking.

**Tres cosas que se arreglaron mirando el resultado, no el código.**

- **Los puntos tienen que ser los del script.** Con 7 espesores (extremos 1 y 7 mm) el rango salía
  52.53 % en vez de 53.11: `t_interior` no es monótona en el espesor (D24/D29) y con 13 puntos
  aparece un extremo interior que con 7 se pierde. Y Z_r estaba en (0.6, 1.0, 1.4) en vez de
  (0.5, 1.0, 2.0): 0.17 % contra 0.28 %. Los puntos son parte del resultado.
- **Sin comparación no se puede declarar éxito.** El primer sello decía "reproduce D30, peor desvío
  0.0 %" con la tabla vacía: `Math.max(...[], 0)` da 0. Ahora, si no hay comparación, lo dice.
- **Un job cancelado no es un job fallido.** El worker levanta al ver la bandera y el `except`
  genérico lo marcaba como `error`. La cola ahora distingue: si la bandera estaba puesta, `cancelado`.
  Igual con la retención: podaba sólo al enviar, así que una ráfaga de jobs en cola dejaba el
  diccionario por encima del tope para siempre; ahora poda también al terminar.

**Por qué no se mide Q todavía.** Es la métrica con la que D29 mostró que la membrana rompe la
resonancia de la cavidad (287 → <5), y era el compromiso de esta fase. No se hace porque falta la
**definición**, no el código: no hay en el repo la función que la calculó, y "Q" puede ser el ancho de
−3 dB en `|p|` o en `|p|²`, sobre el pico de `t_interior` dentro de qué banda. Peor: un Q de 287
implica un ancho relativo de 0.35 %, así que un barrido de 13 puntos en ±8 % **no lo puede
resolver** — haría falta refinamiento adaptativo y del orden de 100 solves. Inventar el número
habría sido lo más fácil y lo más inútil. Queda como el primer trabajo de cuando se fije la
definición (probablemente con Diego: es una decisión de métrica, no de código).

**Decisión de diseño de la cola.** Un solo worker: cada solve ya usa varios núcleos vía BLAS, así
que dos barridos en paralelo no aceleran nada y ensucian los tiempos. La cancelación es cooperativa
—se consulta entre pasos— porque interrumpir dentro de SuperLU no es posible sin threads anidados, y
no vale la pena para esperar 2 s.

## D38 — U4: el 3D en la UI, y las tres cosas que mentían en pantalla

**Qué se entregó.** Una pestaña "3D" con dos modos. *Puntual*: monopolo en una caja de agua con ABC
en las seis caras, comparado contra el Green 3D en un anillo — error **7.0 %** y amplitud **0.969**
(9 261 nodos, 1.0 s). *Plano*: capas apiladas con ABC sólo en *x*, R y T contra la matriz de
transferencia — **R 0.9032 vs 0.9037** (abs 0.0005), **T 0.4343 vs 0.4281** (abs 0.0061), energía
1.0044. Rebanadas *x-y* y *x-z*. Tope de 90 000 nodos y el techo de D32 declarado antes de apretar.

**Hallazgo 1 — la fuente estaba pegada a la ABC (48 % de error que no bajaba al refinar).** La caja
arrancaba en `x = 0` y la fuente se elige como "el nodo más cercano al origen": caía sobre la cara.
El error era 48.7 % a 1.5 λ y 49.0 % a 2.5 λ — constante, que es la firma de un problema de *puesta
en escena*, no de discretización. Centrando la caja: **7.0 %**. Queda un test que verifica que el
nodo de la fuente no pertenece a ninguna cara de la frontera, porque el síntoma (error que no
converge) no dice nada por sí solo.

**Hallazgo 2 — la rebanada salía transpeusta.** Plotly lee `z[fila][columna]` como `[y][x]` y la
matriz sale `[x][y]`: sin `.T` los ejes quedan cruzados. No se ve como un error, se ve como un
patrón raro.

**Hallazgo 3 — una esfera que se veía como un lóbulo.** El campo *era* radialmente simétrico (el
Green cierra al 7 % en todo el anillo), pero el `div` es 740 × 270 px y Plotly estiraba 30 mm × 30 mm
dentro de eso: la circunferencia se veía como un elipse horizontal. Con `scaleanchor` + `constrain:
'domain'` el área de trazado queda cuadrada y la simetría se ve. Es la cuarta variante del mismo
tema (D27, D34, D35): **el render no falla, miente**. Acá costó dos vueltas porque el primer
diagnóstico (mirando el mapa) fue "el solver no da una onda esférica" y era falso.

**Nota de física que la UI dice explícitamente.** En el modo *plano* el campo es uniforme en *y/z*:
en incidencia normal un apilado 3D es un problema 1D. No es un defecto del render; es la física. El
modo puntual es el que muestra algo realmente 3D, y por eso es el que abre por defecto.

## D39 — U6: el sello vive en el servidor, y tres criterios que no son lo que parecían

**Qué se entregó.** `dashboard/validacion.py`: un registro de 10 criterios, cada uno con
descripción, referencia, tolerancia, unidad y **fuente** (fase o ADR). Los cinco endpoints de
simulación devuelven `sello` ya evaluado, con tres estados —`validado`, `fuera`, `sin_referencia`— y
la UI se limita a pintarlo. `GET /api/validacion` publica el catálogo completo para que la UI pueda
listar qué se valida y contra qué, sin inventar nada.

**Lo que se eliminó.** Hasta U5 cada vista tenía su propio umbral en JavaScript: `0.02` en la escena,
`0.53187` en la cóclea, `10.0` en el 3D, `2.0` en el dispersor. Cuatro copias del mismo tipo de
decisión, en el lugar donde nadie las audita y sin decir de dónde salió el número. Ahora hay un test
que falla si cualquiera de esos valores reaparece en el JS.

**Decisión 1 — `sin_referencia` es un estado, no un hueco.** El error clásico de este tipo de código
es que un `None` termine contando como "cumple" (nadie compara contra nada, entonces nada falla). El
registro tiene un tercer estado explícito y un test que verifica que **no** es `validado`. Caso real:
la cóclea con espesor 3 mm no tiene referencia de `t_interior` (D30 sólo cubre 2 mm), y la UI dice
"SIN REFERENCIA" en lugar de mostrar 0.45708 como si estuviera validado.

**Decisión 2 — el resumen no promedia.** Con varios criterios por corrida, el conjunto se queda con
el peor estado: un error grande no puede compensarse con uno chico. Promediar criterios de validación
es cómo se aprueba un simulador que no sirve.

**Hallazgo — el balance de energía de la cóclea empeora al refinar.** Con la membrana y espesor
4 mm, el cierre `P_total + P_membrana` da **1.42 %** a ppw 20, **2.76 %** a ppw 24 y **4.07 %** a
ppw 28. Al principio iba a escribir en el sello "subí ppw" (lo que dice el reflejo); medirlo antes
mostró que ese consejo era **falso** y que refinar empeora. El sello ahora reporta el estado sin
prescribir la solución, y el fenómeno queda como pregunta abierta: no es discretización (crece con
la malla), así que apunta a la medición del flujo —el radio y los anillos de la integral de potencia—
para esa geometría. Es el tipo de cosa que el sello destapó justamente porque ahora hay un lugar
único donde se declara el criterio.

**Consecuencia de diseño.** Un criterio nuevo no se agrega sin la corrida que lo respalde: el
registro es datos, y cada entrada nombra su ADR. Eso deja el repo en un estado donde "¿qué está
validado y contra qué?" se responde leyendo un archivo, no buscando en cinco JS.

## D40 — Consulta a Sim-Datos sobre los puntos abiertos: qué hay, qué no, y qué pasa a ser campo del usuario

Se consultó a `sim-datos` (2026-11-20) por los cuatro puntos abiertos. La respuesta ordena el trabajo:
**un punto tenía referencia analítica clásica y tres no tienen dato citable**, y esos tres pasan a ser
campos editables con tooltip en vez de valores que el código inventaría.

### 1. Q de la cavidad coclear — NO hay dato citable; la definición sí quedó fijada

No existe un Q medido de la resonancia de la cápsula ótica humana, ni en 0.1–1 kHz ni a
100 kHz–1 MHz. Lo que está documentado **no es esto**:

- `Q_ERB = CF/ERB ≈ 12.7 ± 1` a 1 kHz (Shera, Guinan y Oxenham, *PNAS* 2002) es **selectividad de los
  filtros cocleares**, no una resonancia global de la cápsula, y su método aplica desde 1 kHz.
- `Q_10dB` (pico/ancho 10 dB abajo) se usa para respuestas de membrana basilar y depende del nivel
  sonoro y del estado fisiológico.

Definición **operativa** adoptada (es la que implementa `sweeps._campana_decaimiento`):

    Q = f_pico / (f₂ − f₁)   con   |H|²(f) = |H|²_max / 2,  o sea   |H| = |H|_max / √2

Son los mismos cruces con `20·log₁₀|H|` o con `10·log₁₀|H|²`; **no es "media amplitud"**. En un campo
estacionario `|p|²` no representa por sí solo la intensidad neta, así que se mide el ancho de banda del
pico de `t_interior`, no una potencia absorbida. Equivalente temporal: `Q = ω₀τ_A/2` para
`A(t) ∝ exp(−t/τ_A)`.

**Medido con esa definición**: con la membrana puesta, en 85–115 kHz la respuesta **no cae 3 dB dentro
de la banda**, así que no hay un Q que reportar y sí una **cota superior: Q ≤ 3.25** (f_pico 97.5 kHz).
Es coherente con D29 (Q 287 → <5 al entrar la membrana) y es un resultado, no un fallo: el corte con
membrana no tiene un pico estrecho. Validado contra lorentzianas sintéticas de Q conocido: 1.000,
2.000 y 5.000 medidos exactos.

**Lo que se agregó a la UI**, siguiendo la instrucción de Diego: un bloque "Q del pico" con el barrido
como trabajo en cola (banda desde/hasta y cantidad de puntos, editables) y un **campo de Q de
referencia opcional con tooltip** que dice exactamente lo de arriba — que no existe un Q citable de
esta cavidad y que `Q_ERB ≈ 12.7` no es esto. Si el usuario carga un valor de su propia fuente, la UI
compara: con un valor medido reporta el desvío; con una cota dice si la referencia queda dentro o
fuera del alcance del modelo. **Comparar una cota con un número exacto y llamarlo "desvío" sería
inventar precisión**, y por eso son dos mensajes distintos.

### 2. Dispersión penetrable 3D — SÍ hay referencia analítica clásica

- **Esfera fluida**: Anderson (1950), *J. Acoust. Soc. Am.* 22:426–431 — solución por ondas parciales
  (la "serie de Mie acústica"). Esfera `a`, `ρ₁,c₁` en fluido `ρ₀,c₀`, con continuidad de `p` y de
  `(1/ρ)∂ₙp`. **No está restringida a `ka ≪ 1`**: la aproximación de Rayleigh es la que lo está.
- **Esfera sólida elástica**: Faran (1951), *J. Acoust. Soc. Am.* 23:405–418 — ondas longitudinales y
  de corte en un sólido isótropo sumergido; también solución completa por ondas parciales.

Esto cierra el hueco que D32 había dejado anotado: la fase de dispersión penetrable 3D tiene ahora
contra qué validarse, y el primer paso es la esfera fluida de Anderson (más simple y sin acoplamiento
elástico).

### 3. Cápsula ótica y CT — tensor sin dato específico; CT→c sí calibrado

- Peterson y Dechow (2003) midieron propiedades **direccionales 3D** por ultrasonido en cortical
  externa del temporal de 15 cráneos. Pero **cortical externa ≠ porción petrosa ≠ cápsula ótica**, y no
  hay tensor `E₁,E₂,E₃,ν_ij,G_ij` citable ni verificado para la cápsula ótica humana, ni validez
  dinámica hasta 1 MHz. No se traslada el tensor de la bóveda.
- **CT→velocidad sí tiene calibración experimental**: Webb et al. (2018), 91 fragmentos de dos cráneos,
  500 kHz, `c = m·HU + b` con `R²` entre 0.23 y 0.53 según adquisición. Ejemplo de la tabla I (GE,
  120 kVp, kernel BONEPLUS, reconstrucción convencional): **m = 0.75 m·s⁻¹·HU⁻¹, b = 1.32 km/s,
  R² = 0.46**. No es ley universal: depende de energía, reconstrucción y hueso.

### 4. α del hueso a 100 kHz — sigue siendo extrapolación del material

Hay datos publicados **por debajo** de 270 kHz, pero son **pérdida de inserción del sistema**, no α
intrínseca: Li et al. (2025), *Phys. Med. Biol.* 70:215024 (20 cráneos ex vivo, barrido 100–290 kHz en
pasos de 10 kHz) y Li et al. (2025), *JASA* 157:2336–2349 (a 220 kHz, 3.6 ± 3.4 dB en históricos y
2.9 ± 1.8 dB en conservados con Thiel). Incluyen interfaces, interferencia y geometría. Para **α del
material** a 100 kHz no hay dato citable que separe esas contribuciones: usar Pichardo/White sigue
siendo extrapolación, ahora contrastable contra mediciones de transmisión del sistema completo.

### Qué significa para la UI

Los tres puntos sin dato citable quedan como **campos que el usuario ingresa**, con tooltip que explica
qué hace el número y qué se sabe de él — nunca como constante en el código:

| campo | estado | dónde |
|---|---|---|
| Q de referencia | opcional, sin valor por defecto | pestaña cóclea, bloque "Q del pico" |
| propiedades de la membrana basilar | `VERIFY` (ballpark), editables | pestaña cóclea |
| tensor elástico de la cápsula | no existe citable → isótropo `VERIFY` (`c_p`, `ρ` fijos y declarados) | pestaña cóclea |
| `m`, `b` de CT→c | calibración citable con ejemplo (Webb et al.) | fase de geometría de CT |

La regla que sale de acá: cuando `sim-datos` dice "no hay dato citable", el número **no** se rellena —
se expone como entrada del usuario con la explicación de su efecto, y el resultado lo declara.

## D41 — Dispersión penetrable 3D: la serie de Anderson, y el error de signo que traía D32

**Qué se entregó.** `analytic/sphere.py`: la solución por ondas parciales de **Anderson (1950)** para
una onda plana sobre una esfera fluida (la "Mie acústica"), con `boundary='fluid'|'rigid'`, campos
exterior/interior, derivadas radiales y truncación automática. Y el FEM 3D penetrable funcionando:
`assemble_scattered` acepta `background_k` + `background_rho` + `p_inc` e integra el término de
contraste correcto.

**Validación de la referencia primero.** Una serie mal implementada es peor que ninguna: el FEM
"coincidiría" con un error. Cuatro chequeos, todos independientes:

- **Esfera casada** (`ρ₁ = ρ₀`, `c₁ = c₀`): `Aₙ = 6.1e-24`. No hay dispersión.
- **Continuidad en `r = a`** de `p` y de `(1/ρ)∂ᵣp` con los tres campos (incidente, dispersado,
  interior), a ka = 0.5 / 1 / 3 / 5: saltos de **1e-14 % a 1e-12 %**.
- **Límite rígido**: `ρ₁ → 1e9` da la condición `∂ᵣp = 0` al **0.0001 %**, y la derivada normal del
  total rígido se anula al **4.8e-13**.
- Coherencia con la truncación: `n_max ≈ ka + 4·(ka)^(1/3) + 2` da 23–28 términos para ka ≤ 5.

**El bug que D32 no pudo cerrar.** En el coeficiente de reflexión de la esfera,

    Aₙ = (G·Rₙ − Rₙ′) / (Hₙ′ − G·Hₙ),    G = (ρ₀/ρ₁)·(Jₙ′/Jₙ)

el cociente de densidades va **así**: con `G = (ρ₀/ρ₁)(J′/J)`, el límite `ρ₁ → ∞` da `G → 0` y
`Aₙ → −Rₙ′/Hₙ′`, que es la condición **rígida** (`∂ᵣp = 0`). Invertido (`ρ₁/ρ₀`) daba el límite
**blando** (`p = 0`). Y como la presión queda continua *por construcción* (la `Bₙ` sale de esa
ecuación), el error no aparecía en `p`: aparecía en el **salto de flujo, 73.7 %, constante en ka** —
la firma de un factor sistemático, no de discretización. Fue el chequeo del límite rígido el que lo
delató; mirando sólo el campo dispersado contra el FEM hubiera sido indistinguible de un error de
malla.

**El FEM contra la serie** (misma caja, esfera por centroides de elemento, ABC en las seis caras):

- ka = 1, ppw 12: **3.34 %** (2 197 nodos, 228 elementos en la esfera), 0.19 s.
- ka = 2, ppw 12: **2.73 %** (13 824 nodos, 1 446 elementos), 0.4 s.
- Convergencia a ka = 1: **5.72 %** (ppw 8) → 3.32 % (12) → 1.93 % (16) → 1.65 % (20, 1.1 s).
- Sin contraste (`k` y `ρ` iguales al fondo): `p_s = 0` **exacto** (1e-9).

**Por qué funciona ahora y no antes.** La forma débil correcta del dispersado es

    ∫(1/ρ)∇p_s·∇v − (k²/ρ)p_s v − ABC(p_s) = −C(p_i),
    C(p_i) = ∫(1/ρ − 1/ρ₀)∇p_i·∇v − (k²/ρ − k₀²/ρ₀)p_i v

El operador es el **local** (lo que ya armaba `assemble_scattered`) y el lado derecho lleva el
contraste **completo**: los dos términos restan el fondo (`1/ρ₀`, `k₀²/ρ₀`), no sólo `k²`. El intento
de D32 usaba `(∇p_i·∇v − k²p_i v)/ρ`, que no es `C(p_i)`; por eso con `ρ` variable la comparación se
iba a 128 % y con `ρ` constante (donde los términos de fondo se cancelaban dentro del volumen) la R
extraída salía 0.96 en vez de 0.33. Y se integra **sólo donde el medio difiere**: fuera del dispersor
el integrando es cero en el continuo, pero discretizado sobre todo el dominio deja el flujo espurio
de la frontera.

**Siguiente paso con referencia ya disponible**: Faran (1951) para la esfera **elástica** sumergida
(ondas P y S), que es el mismo tipo de serie con dos potenciales en el sólido.

## D42 — La esfera penetrable en la UI, y por qué el mapa mostraba casi nada

**Qué se entregó.** Un tercer modo en la pestaña 3D: *esfera penetrable*, con material de la esfera y
`ka` editables. Resuelve el campo dispersado con el término de contraste de D41 y sella contra la
serie de Anderson. Medido: **2.94 %** de error a ka = 1 (9 261 nodos, 180 elementos en la esfera,
960 puntos de comparación), *1.2 s*. Criterio nuevo en el registro: `volumen_esfera_serie`, 5 %.

**Hallazgo 1 — el mapa del campo total no muestra nada.** El primer intento graficó `|p_total|`: como
la onda incidente vale ≈ 1 en todo el dominio, la perturbación de la esfera (del orden del 5 %) es
invisible en una escala de 30 dB — el mapa salía **casi uniforme, con 0.3 dB de rango**. Mostrando el
**campo dispersado** el rango pasa a decenas de dB. Es la tercera vez que aparece el mismo patrón
(campo uniforme con agua/tejido en 2D, esfera en 3D): *cuando el campo total está dominado por la
incidente, lo que hay que graficar es el dispersado*.

**Hallazgo 2 — el plano equivocado.** Con la propagación en `z`, el corte `x-y` (a `z = 0`) es
perpendicular al haz y sale casi simétrico: una mancha oscura central. Los lóbulos viven en el plano
`x-z`. En el modo esfera se invirtió cuál es el mapa principal, y el título dice cuál es cuál.

**Lo que el mapa NO demuestra, y por qué igual sirve.** A ka = 1 con una esfera de contraste moderado
el campo dispersado está dominado por el modo monopolar: se ve un lóbulo central que decae
radialmente, no la asimetría dipolar de una esfera rígida. Eso es física, no un defecto — pero
tampoco es la validación: **la validación son los 960 puntos comparados contra la serie** (2.94 %),
que es un chequeo muchísimo más fuerte que mirar un mapa. La lección es que el sello numérico y la
imagen tienen que coincidir en lo que afirman, y cuando no, manda el número.

**Costo y tope.** 1.2 s con 9 261 nodos en una caja de 30 mm (2 λ). La malla tiene que resolver la
esfera: con `ppw = 6` y `ka = 0.2` no cae ningún elemento adentro y el pedido se **rechaza** con un
mensaje que dice qué subir, en vez de resolver una caja homogénea y devolver un error del 100 %.

## D43 — Barrido de frecuencia de la escena: la resonancia de media onda como criterio

**Qué se entregó.** Un barrido en frecuencia de la escena 2D como trabajo de la cola
(`sweeps.barrido_escena`), con banda y cantidad de puntos editables, que devuelve R y T por punto y
dibuja las **resonancias predichas** para que se vea si el modelo las encuentra.

**El criterio, que es físico y no una comparación consigo mismo.** Una losa de espesor `t` en un
medio transmite al máximo en las resonancias de media onda:

    f_n = n·c_losa / (2t)

Con hueso como fluido (c = 3500 m/s) y t = 10 mm: **175 y 350 kHz**. Medido en 100–400 kHz con 41
puntos: picos en **176.7 y 353.8 kHz** (≤ 1 % de la predicción, y el corrimiento es el paso del
barrido de 7.5 kHz). En la resonancia coinciden las tres cosas que tienen que coincidir: `T` sube a
**0.999**, `R` baja a **0.031** en la misma frecuencia, y la matriz de transferencia reproduce el
mismo pico. Un barrido que se equivocara de espesor, de velocidad o de fase se vería acá.

**El sello marca "fuera" a ppw 12, y está bien.** El error de T del barrido llega a **5.90 %** por
encima de 350 kHz, contra el 2 % de la Fase 1. La causa es la malla, no el barrido: con `ppw` fijo,
a 400 kHz la losa de 10 mm tiene menos elementos por longitud de onda. Subiendo a `ppw = 20` el peor
error baja a **2.08 %**. Es el caso opuesto al de D39 (donde refinar *empeoraba* el balance de la
cóclea): acá el consejo "subí ppw" sí es cierto, y por eso el barrido **informa el peor error en vez
de un promedio** — un promedio de R y T sobre la banda escondería justo el punto donde la malla deja
de alcanzar.

**Detalle de presentación.** El barrido usa la escena de la pestaña "escena FEM" si el usuario la
editó (mismo objeto `ESC` en el cliente) y si no una losa por defecto; y el gráfico de sensibilidad
tiene una guarda para no pisar el del barrido de frecuencia: sin ella, el último que corre borra al
otro, que es un bug de UI que no da error, sólo muestra el gráfico equivocado.

## D44 — Faran (1951): la derivación arrancada, y dónde se paró a propósito

**Estado: no entregado.** Se empezó la serie de la esfera elástica y se **borró** el módulo a medio
escribir en vez de dejarlo en el repo. Un coeficiente mal derivado no falla ruidosamente: da un campo
que parece razonable y no coincide con nada, y descubrirlo después cuesta más que escribirlo bien de
entrada. Esta ADR guarda lo que sí quedó verificado, para que la próxima sesión continúe en vez de
reempezar.

**Lo que quedó verificado con sympy** (no a mano). Descomposición
`u = ∇Φ + ∇×(∇×(Ψ r̂))`, con `Φ = f(r)P(θ)` y `Ψ = g(r)Q(θ)`:

    u_r = Φ_r + n(n+1)·Ψ/r²
    u_θ = (Φ_θ + Ψ_rθ)/r

Las dos salieron de sympy calculando el rotacional-rotacional en esféricas con simetría axial. El
chequeo vale: a mano había puesto en duda el segundo término y estaba equivocado en la duda, no en la
fórmula. Con `Ψ = 0` queda `u = ∇Φ` (`u_r = Φ_r`, `u_θ = Φ_θ/r`), que es el caso puramente
longitudinal.

**Lo que quedó derivado pero sin armar.** Las expresiones de `σ_rr` y `σ_rθ` en función de
`f, f′, f″, P, P″` y `g, g′, g″, Q, Q″` (con `ε` en esféricas y `σ = λ·tr(ε)I + 2μ·ε`) salieron de
sympy, pero están en la forma cruda: falta sustituir la ecuación de Legendre
(`P″ + cotθ·P′ = −n(n+1)P`), quedarse con `f = jₙ(k_α r)`, `g = jₙ(k_β r)` y armar el sistema **3×3**
por orden `n` para `(Aₙ, Bₙ, Cₙ)` con las tres condiciones en `r = a`:

    1. u_r(sólido) = u_r(fluido),  con u_r(fluido) = (1/(ρ₀ω²))·∂p/∂r
    2. σ_rr = −p
    3. σ_rθ = 0

**Cómo se valida cuando esté** — y esto sí está decidido, porque el problema de Faran es que no hay
tabla de referencia a mano:

1. **Teorema óptico**, que no necesita referencia: para un dispersor sin pérdidas,
   `σ_ext = σ_sca`, con `σ_ext = (4π/k₀²)·Im f(0)` y `σ_sca = ∫|f|²dΩ`. Es conservación de energía: si
   el sistema de coeficientes está mal, no cierra. Es el chequeo más fuerte disponible sin tabla.
2. **Residuos de las tres condiciones** en `r = a` después de resolver: tienen que bajar a ~1e-12.
3. **Contraste bajo**: con la esfera casi igual al fluido, `Cₙ → 0`.

**Lo que NO sirve como validación, y por qué.** Comparar contra `sphere.py` en el límite `μ → 0`
—que es el límite *físico* correcto, el sólido se vuelve fluido— es **numéricamente singular**:
`k_β = ω/c_s → ∞` y `jₙ(k_β a)` oscila sin límite. La comparación mediría el ruido del límite, no la
formulación.

**Nota de entorno.** La derivación usa `sympy`, que se instaló en el `.venv` para esto. **No** es
dependencia de ejecución: `analytic/faran.py` cuando exista sólo necesitará `scipy.special`. Si en
algún momento se agrega a `pyproject.toml`, va en las de desarrollo.

## D45 — Faran entregado: la esfera elástica, validada sin tabla de referencia

**Qué se entregó.** `analytic/faran.py`, la serie de ondas parciales de la esfera **elástica**
sumergida: dos potenciales adentro (P y S), tres condiciones en `r = a` y sistema 3×3 por orden `n`.
Es la contraparte de Anderson (esfera fluida, `sphere.py`).

**Cómo se valida, ya que no hay tabla a mano.** Cinco chequeos, y los cinco son independientes de la
tabla que no hay:

- **Condición unitaria por coeficiente**, `Re Cₙ = −|Cₙ|²`: se cumple a **1e-16–1e-19**, o sea
  precisión de máquina. Es el chequeo más fuerte porque es una relación cuadrática entre la parte
  real y el módulo de cada coeficiente: no se cumple de casualidad.
- **Teorema óptico**, `σ_ext = σ_sca`: desvío de **1.0e-3 %** (acero, hueso y aluminio), que es la
  cuadratura de `σ_sca`, no la serie.
- **Límite de fluido**: con `ρ` y `c_p` iguales al agua y `c_s` bajando, `|C₀|` cae
  **1.62e-3 → 4.04e-4 → 1.01e-4** para `c_s` 100/50/25, o sea como `c_s²`.
- **Aproximación al rígido**: contra la esfera rígida de `sphere.py` (otro camino, otra derivación),
  aluminio (Z/Z₀ = 11.4) **5.2e-3**, acero (30.7) **2.3e-3**, tungsteno (66.9) **1.2e-3** — monótona
  con la impedancia, como tiene que ser. A `ka = 0.2` el acero queda a **3.5e-5** del rígido.
- **Residuos** de las tres condiciones de borde: **1e-15 / 1e-17 / 1e-18**.

**Los tres errores propios que hubo que encontrar, porque los tres dan resultados plausibles.**

1. **`n = 0` en las recurrencias de Bessel.** `jₙ′ = jₙ₋₁ − (n+1)jₙ/x` pide `j₋₁` cuando `n = 0`, y
   `scipy` devuelve `inf`: el NaN se propagaba a todo el sistema y salía `sigma = nan`. Van aparte
   (`j₀′ = −j₁`, `h₀′ = −h₁`).
2. **Un `k₀` de más en la continuidad de desplazamiento.** La fila 1 decía
   `k₀/(ρ₀ω²)·∂h/∂r` cuando `∂h/∂r` **ya** trae el factor `k₀`: duplicado, el sistema entero salía
   como si cualquier esfera fuera rígida. El síntoma que lo delató fue el caso de control —un sólido
   con la misma `ρ` y `c_p` que el agua **tiene** que dar `C ≈ 0` y daba `0.2125`, el valor rígido—;
   el desajuste de tensiones era exactamente `k₀ = 419`. Que la condición unitaria y el teorema
   óptico se cumplieran igual es lo importante: **los dos son propiedades estructurales del sistema y
   los dos pasaban con la física mal**. No reemplazan a un caso límite con respuesta conocida.
3. **Un residuo normalizado contra sí mismo.** Para `σ_rθ` usé `|Σσ_rθ|/Σ|σ_rθ|`, que da ~1 aunque el
   residuo sea 1e-17 (los términos no se cancelan: cada uno ya es cero). Se normaliza contra la
   presión y la tensión normal, que son las que el sistema iguala.

**Lo que NO sirve como validación**, y quedó escrito antes de programar: comparar contra `sphere.py`
en el límite `μ → 0` exacto. Es físicamente correcto pero numéricamente singular (`k_β → ∞`). El
límite de fluido se testea con `c_s` chico pero finito, donde la convergencia se ve.

**Lo que todavía no hay.** Esto es la **referencia**, no una simulación: no hay solver FEM elástico 3D
que comparar contra ella (la fase 6b del plan). La serie sirve desde ya para dos cosas: dar la
respuesta exacta cuando esté el solver, y verificar los límites del propio solver (rígido y fluido)
antes de mirar un caso con contraste real.

## D46 — Fase 6b, primera mitad: malla de esfera, elasticidad 3D y acoplamiento (y lo que falta)

**Entregado y validado.**

- `meshing3d_sphere.sphere_shell_tet_mesh`: esfera sólida rodeada de una cáscara de fluido, con la
  interfaz **compartida nodo a nodo**, hecha con **Gmsh** (4.15.2, ahora dependencia). Volumen del
  sólido: error **3.51 % a h = a/3** y **1.57 % a h = a/4.5**; área de interfaz al **98.4 %** de
  `4πa²`; conformidad verificada (cada cara de interfaz es frontera de las dos regiones).
- `fem/elastic3d`: elasticidad vectorial P1 sobre tetraedros. Los seis **modos de cuerpo rígido**
  caen en el núcleo de `K` a **4.1e-16**, la **energía de una deformación uniforme** da exacta a
  **1e-16** (P1 reproduce campos lineales, así que no hay tolerancia de malla) y la masa total es
  `3ρV` exacta.
- `fem/coupled3d`: los dos bloques acoplados, con ABC-2 en la esfera exterior.

**Qué NO está validado, y hay que decirlo claro.** La comparación con **Faran** para acero a
`ka = 1` da **22 %**, y ese número **no** se reporta como validación. Lo que se sabe del error:

- El patrón angular coincide en **forma y fase** (razón de amplitudes 0.80–0.88 casi constante, fases
  dentro de unos pocos grados salvo en `θ = 0`), pero la amplitud sale **18 % baja**. Un déficit
  uniforme no es un error de forma: es algo que escala.
- El **desplazamiento interior** del FEM sale casi uniforme —el sólido se mueve como rígido— mientras
  Faran predice deformación. Para acero la respuesta exterior es casi rígida (los `Cₙ` de Faran
  difieren de la esfera rígida en 2.3e-3), así que la diferencia interior debería pesar poco… y sin
  embargo pesa.
- **A/B hecho**: el mismo dominio resuelto con el Helmholtz **ya validado** (D41) más ABC-2 da
  **10.35 %** contra la esfera rígida. O sea que el dominio y el ABC aportan ~10 % por sí solos en
  esta malla, y en el **límite rígido** el acoplado da **10.66 %** — coincidiendo con el camino
  validado. El acoplamiento está bien; lo que falta aislar es la respuesta **elástica**.

**Los cuatro errores propios de esta mitad** (los cuatro daban resultados plausibles):

1. **Malla de capas hecha a mano, no conforme.** Icósfera + capas de prismas: el split de 3
   tetraedros por prisma no respeta las diagonales de las caras compartidas. Medido: **2 240 caras
   internas** apareciendo una sola vez en vez de 320. Se abandonó por Gmsh en vez de pelear el split.
2. **`∇N₀` sumando el eje equivocado.** `-inv.sum(axis=1)` sobre un stack `(M,3,3)` no da la suma de
   columnas. Efecto: los gradientes no sumaban cero, una **traslación producía deformación** y la
   energía salía **5× más grande**. Se reescribió sumando los tres gradientes ya calculados.
3. **`∇Nᵢ` como fila de `J⁻¹` en vez de columna.** Es `J·∇N = eᵢ`. Con la transpuesta mal, los modos
   rígidos no se anulaban (6.6e-4) y la energía quedaba 4.6 % baja — un error del tamaño justo para
   confundirse con convergencia.
4. **Espacios reducidos.** Armar `K` sobre todos los nodos deja filas vacías en los nodos que sólo
   existen dentro de la esfera (44 medidas) y el sistema global sale singular. Lo mismo del lado del
   fluido. El error de SuperLU es "Factor is exactly singular", que no dice nada de la causa.

**Lo siguiente, concreto.** Aislar la respuesta elástica: comparar el desplazamiento interior del FEM
contra `Φ`/`Ψ` de Faran punto por punto (no sólo en el eje, donde `Pₙ(1) = 1` esconde los órdenes
`n ≥ 1`), y con eso decidir de qué lado está la diferencia. Recién después tiene sentido refinar:
la convergencia medida (22.4 → 19.6 → 17.5 % al pasar de h = a/3 a a/5) es **O(h)**, y un error de
primer orden no es el que deja una formulación correcta sobre una superficie facetada.

### D46 — apéndice: el 22 % queda **localizado**, no resuelto

Después de la primera mitad se hizo el diagnóstico que faltaba. Los tres resultados que importan:

**1. El campo interior de Faran está bien.** Se agregó `faran.campo_interior` (desplazamiento dentro
del sólido en puntos arbitrarios) y se validó contra lo que ya estaba validado: en ``r = a`` la
componente radial del sólido tiene que igualar a la del fluido, ``(1/ρ₀ω²)·∂p/∂r``. Da
**1.0000 en todos los ángulos** (±1e-4). Hacía falta medirlo fuera del eje: en el eje
``Pₙ(1) = 1`` y ``dPₙ/dθ = 0``, así que mirar sólo ahí **esconde todos los órdenes ``n ≥ 1```,
que es donde aparece la deformación.

**2. La rigidez elástica y la fuerza de acoplamiento están bien.** Test estático con respuesta
exacta: presión uniforme sobre la esfera, ``u_r = −p·r/(3K)``. Medido **0.9994** (desvío 1.5 %), y
0.9998 sobre la superficie. Ese test prueba a la vez ``K``, la matriz constitutiva y el término de
fuerza ``−∫p(n·v)``. (Detalle: el problema estático es singular en los seis modos rígidos, así que
hay que proyectarlos antes de comparar; sin eso los números salen con desvíos de 200 % que parecen
un bug de física y son una traslación arbitraria del solver.)

**3. El error está en el balance DINÁMICO, y se ve separando traslación de deformación.** En el
interior, para acero a ``ka = 1``:

- Faran: traslación **4.79e-15** m, deformación **3.86e-15** m (55 % del total).
- FEM: traslación **5.62e-14** m (12× la de Faran), deformación **6.63e-16** m (**0.17×**).

O sea: el sólido del FEM **casi no deforma** y **traslada 12 veces de más**. Ni la rigidez ni la
estática pueden producir eso — las dos están validadas. Lo que falla es cómo se equilibran las dos
contribuciones dinámicas en la interfaz: la fuerza del incidente ``−∫p_i(n·v)`` y la reacción
``∫p_s(n·v)``. En una esfera rígida esas dos se cancelan exactamente (no hay traslación en una onda
plana); si la cancelación queda desbalanceada, el residuo es una **fuerza neta** que traslada al
sólido sin deformarlo, que es exactamente el patrón medido.

**Siguiente paso, concreto y acotado.** Medir la fuerza neta en la solución del FEM:
``F = ∫_Γ (σ·n + p·n) dS`` sobre la interfaz. Si ``F ≠ 0``, el desbalance está en los dos términos
del acoplamiento y se arregla ahí; si ``F ≈ 0``, el problema es la *distribución* de la tracción
(el sólido recibe una presión casi uniforme en vez de la que corresponde), y el sospechoso pasa a ser
la proyección de ``n`` sobre los grados de libertad del sólido en ``C``.

Lo que **no** hay que hacer todavía: refinar. La convergencia medida es O(h) y el error no es de
discretización; refinar una formulación con un desbalance dinámico sólo hace más caro el mismo error.

### D46 — segundo apéndice: el 22 % **resuelto**, y el piso que queda (9.7 %)

**La aislación que faltaba.** Se resolvió **sólo el fluido** con el desplazamiento impuesto en la
interfaz, dos veces:

- con el desplazamiento **exacto** de Faran: error **9.73 %**, razón de amplitudes 0.980.
- con el desplazamiento que resuelve el **FEM acoplado**: error **21.52 %**, y el desplazamiento
  normal del FEM es **20.3×** el exacto.

O sea: el fluido y el ABC están bien (ese 9.73 % es el piso de la malla y el dominio, coherente con el
10.35 % que da el camino validado de D41 en la misma malla). **Todo el error estaba en el sólido.**

**La causa: los seis modos de cuerpo rígido.** El sólido libre no tiene rigidez elástica en esos
modos; la única que los contiene es la reacción del fluido, que es chica. Para una onda plana sobre un
cuerpo real la fuerza neta es **cero**, y la ecuación de movimiento `−ω²Mu = F` implica entonces
**desplazamiento medio nulo**. Pero la malla no cancela esa fuerza exactamente, el residuo excita la
traslación sin freno y esa traslación contamina todo el campo dispersado.

**La prueba y el arreglo.** Primero con una penalización (experimento): el error bajó de 19.36 % a
13.11 % y la razón de amplitudes pasó de 0.843 a **1.012**. Después, con la versión principiada —
restricción de Lagrange de traslación media nula, que además **devuelve los multiplicadores** (la
fuerza y el torque que la malla no canceló):

- `h = a/3.5`: 19.36 % → **10.20 %**, razón 0.843 → 0.970
- `h = a/4.5`: 18.15 % → **9.70 %**, razón 0.858 → 0.988
- `h = a/5.5`: 17.68 % → **9.21 %**, razón 0.865 → **0.997**
- residuo de fuerza: **4e-7 → 7.6e-8 N**

El límite rígido no cambia (9.54 % → 9.60 %): ahí el sólido ya no se mueve y la restricción no tiene
nada que atar.

**La restricción tiene dominio de validez, y se midió dónde se rompe.** Con la esfera
**indistinguible del medio**, apagarla es obligatorio: sin restricción la dispersión es 3.0e-2 (bien),
con restricción es **1.6 Pa sobre una incidente de 1 Pa**. Es coherente: si el cuerpo es transparente,
la onda pasa de largo y **sí** ejerce fuerza neta, así que el desplazamiento medio no es cero y
forzarlo mete dispersión falsa. `fijar_rigidos` queda como opción explícita de `resolver`, con default
`True` y el caso de excepción documentado y testeado.

**Qué queda, y es honesto llamarlo piso.** El 9.2-10.7 % restante **no es del sólido**: es el mismo
número que da el fluido solo con datos exactos (9.73 %). Descomposición de ese piso: superficie
facetada (el poliedro inscripto tiene radio efectivo 0.96-0.99 a), ABC-2 sobre una esfera de radio
2.5a (`k₀b = 2.5`, chico para el ABC) y elementos P1 en el fluido. Para bajarlo hay tres caminos
medibles y en ese orden: agrandar el dominio manteniendo el tamaño de elemento (el ABC pesa menos),
afinar la malla (el facetado pesa menos) y pasar a elementos de frontera curva (deal.II/FEniCSx,
repos de `docs/referencias.md`).

**Un barrido que no sirvió y por qué.** Se barrió `ka = 0.2 … 1.0` a radio variable: 41 % a `ka = 0.2`,
15 % a `ka = 0.5`, 19 % a `ka = 1`. El 41 % de `ka = 0.2` **no dice nada del acoplamiento**: con
`ka = 0.2` el dominio exterior queda en `k₀b = 0.5`, o sea media longitud de onda para el ABC. Es el
mismo error de antes (mirar un número sin preguntar qué lo domina).

### D46 — tercer apéndice: de dónde salen los 9 puntos, y el dominio pasa a 4a

**El experimento.** Fluido solo con el desplazamiento exacto de Faran en la interfaz (9.73 % de base),
moviendo **una cosa por vez**. Ese camino aísla dominio y malla del sólido, y cuesta 2 s en vez de 30.

**Dominio (`h = a/4.5` fijo):**

- `b = 2.5a` (`k₀b = 2.5`): **9.30 %** — 5 934 nodos, 2.0 s
- `b = 4.0a` (`k₀b = 4.0`): **4.06 %** — 21 398 nodos, 24 s
- `b = 6.0a` (`k₀b = 6.0`): **3.96 %** — 67 113 nodos, 303 s

O sea: **el ABC se lleva ~5 puntos de los 9**, y a 4a ya está por debajo del piso de malla. Ir a 6a
mejora 0.1 puntos por **12× el costo**. El ABC local no es el problema: el problema era ponerle un
dominio de `k₀b = 2.5`, que es media longitud de onda.

**Malla (en `b = 4a`):**

- `h = a/3`: 10.24 % (área de interfaz 0.9807 de la esfera)
- `h = a/4.5`: **4.06 %** (0.9913) — orden aparente 2.28
- `h = a/6`: **3.03 %** (0.9945) — orden 1.02

El orden **baja a 1** al refinar, y el área de la interfaz converge despacio: a partir de acá el que
manda es el **facetado** (el poliedro inscripto tiene área y radio menores que la esfera). Refinar más
no arregla eso — lo arreglan elementos de frontera curva, y por eso `docs/referencias.md` apunta a
deal.II/FEniCSx para ese paso.

**Acoplado completo, con el sólido resuelto:**

| configuración | error | nodos | tiempo |
|---|---|---|---|
| 2.5a, h=a/3.5 | 10.20 % | 2 858 | 1.3 s |
| 3.0a, h=a/4.5 | 7.62 % | 9 755 | 13.9 s |
| **4.0a, h=a/4.5** | **4.21 %** | 21 398 | 32 s |
| 4.0a, h=a/6 | 3.01 % | 47 574 | 211 s |

Comparado con el fluido solo y datos exactos a 4a (4.06 %), el **sólido aporta 0.15 puntos**: la
resolución elástica dejó de ser el cuello de botella.

**Decisión de API.** `sphere_shell_tet_mesh` toma `outer_radius` y `spacing` por defecto
(`4·radius` y `radius/4.5`) y **avisa** si el radio exterior baja de `3·radius`, con el número medido
en el mensaje. Un dominio chico sigue permitido —para una vista rápida sirve— pero ya no puede pasar
silenciosamente a un número que se reporte.

**Lo que queda sobre la mesa.** Si se quiere el 4 % **al costo de 1.3 s**, el camino es una condición
**Dirichlet-to-Neumann exacta** sobre la esfera exterior (expansión en armónicos esféricos con
`hₗ′(k b)/hₗ(k b)` por modo, en vez del ABC local): es exacta para una frontera esférica y permitiría
`b = 2.5a`. No está implementada; queda como el siguiente paso si el costo de 32 s molesta en la UI.

## D47 — U7: la esfera elástica en la UI, y el hallazgo del sólido de bajo contraste

**Qué se entregó.** Pestaña **"esfera elástica"**: material, `ka`, frecuencia, dos mallas (rápida y
precisa) y el interruptor de la restricción de modos rígidos, todo con tooltip. El servidor resuelve
con `coupled3d` y contrasta contra Faran en la misma corrida; la UI pinta el campo dispersado en el
plano x-z, el módulo `|p_s|` contra θ y la **diferencia de fase** (no las dos fases: sus saltos de
±180° se leen como discontinuidades del modelo y no lo son). El sello sale del registro
(`volumen_esfera_elastica`, 6 %), no de una constante en el JS.

**Error contra Faran por material**, malla precisa (dominio 4a, `h = a/4.5`, 21 398 nodos, ~40 s):

- acero: **4.29 %** (razón de amplitudes 0.969)
- aluminio: **3.99 %** (0.972)
- tungsteno: **4.71 %** (0.966)
- **hueso: 13.16 %** (1.058)

**El hueso no es un problema de malla: es de modelo.** Se probó apagar la restricción de rígidos y
empeora a **54.61 %** (contra 19.52 % del acero), así que la restricción no es la culpable — sigue
siendo imprescindible. Lo que pasa es que el hueso tiene un contraste de impedancia de **4.6** contra
**30.7** del acero: con contraste bajo, el cuerpo **se mueve parcialmente con el fluido**, y la
condición de traslación media nula deja de ser exacta. Es la misma frontera que apareció con la
esfera transparente (donde la restricción mete dispersión falsa de 1.6 Pa) pero en un régimen
intermedio, y por eso hay que decirlo en vez de ajustar el umbral hasta que entre.

**Lo que arreglaría eso, y no es una malla más fina.** Formular el sólido también en **campo
dispersado** (`u_s = u − u_i`, con `u_i` la incidente del fluido prolongada al interior): así el caso
transparente da `u_s = 0` **exactamente**, sin restricción, y los modos rígidos sólo se excitan por el
contraste real. Arrastra el mismo tipo de término que D41 (hay que aplicar el **operador completo** a
la incidente restringido al sólido, no `−(k²−k_b²)`), así que se hace con la misma disciplina: primero
la derivación, después el test de que el caso transparente da exactamente cero.

**Mientras tanto, la UI dice la verdad.** El sello marca `fuera` para el hueso con el número medido y
la fuente; no se subió el umbral a 13 % para que la tabla se viera verde. Un criterio que se ajusta
hasta que pasa deja de ser un criterio.

**Dos arreglos de plataforma que salieron de esto:**

- **Gmsh desde un hilo.** `gmsh.initialize()` reemplaza el manejador de SIGINT, y eso revienta en el
  hilo del worker de la cola: *"signal only works in main thread"*. Se pasa `interruptible=False`.
  Sin eso la vista no funcionaba desde la UI y sí desde la consola — el peor tipo de bug.
- **Interpolación para el mapa.** Probar los 8 672 tetraedros por cada uno de los 6 561 puntos del
  mapa tardaba **15.1 s** (el solve son 1.8 s). Con un KD-tree y 32 vecinos baja a **1.08 s**, y se
  verificó que da exactamente lo mismo que la búsqueda exhaustiva.

### D47 — apéndice: el tejido casi fluido y la restricción que no se puede imponer

El barrido por material que se había lanzado dejó un número que no se podía publicar así: con
**tejido casi fluido** (`c_s = 30 m/s`, contraste de impedancia **1.08×** el agua) el campo dispersado
sale **53×** el exacto — **5264 %** de error. No es un problema de malla: es la restricción de
traslación media nula, que **supone que la onda plana no ejerce fuerza neta**, y con contraste bajo el
cuerpo se mueve parcialmente **con** el fluido y la suposición se rompe. Sin la restricción, el mismo
caso da razón **1.13×** y error absoluto **7.2e-3 Pa** sobre una incidente de 1 Pa.

**Lo que se hizo, en vez de subir el umbral:**

1. **La restricción se desactiva sola** por debajo de un contraste de impedancia de **2.5**, y se
   informa (`restriccion_automatica`). El corte sale de dos puntos medidos: 4.55 (hueso) da 15 % con
   la restricción y 1.08 (tejido) da 5264 %; no de una fórmula.
2. **El sello deja de usar el error relativo cuando el campo es débil.** Si `|p_s|` máximo no llega al
   5 % de la incidente, el denominador se achica pero el error absoluto no, y "23 % de error" pasa a
   ser un número sin sentido que además asusta. En ese régimen el sello va a **`sin_referencia`** y se
   muestra el **error absoluto** contra la incidente.
3. **La UI lo dice en la ficha**: contraste, si la restricción se desactivó sola y por qué, y el error
   absoluto con su escala.

**Tabla completa, `ka = 1`, dominio 4a, `h = a/4.5`, 21 398 nodos:**

- acero (contraste 30.7): **4.29 %** — validado, es el caso del criterio
- aluminio (11.4): **3.99 %**
- tungsteno (66.9): **4.71 %**
- hueso (4.6): **13.16 %** — error de modelo por contraste moderado
- tejido casi fluido (1.08): con restricción **5264 %**; sin restricción, error absoluto *7.2e-3 Pa*

**La lección, que es la misma de siempre en este proyecto:** el criterio no se ajusta para que la
tabla quede verde. Un material que la formulación no sabe representar se marca, se explica y se
desactiva lo que corresponde — y el arreglo de fondo sigue siendo el mismo: formular el sólido en
campo dispersado, donde el caso transparente da exactamente cero **sin** necesidad de restringir nada.

**La tabla completa, y de dónde sale cada error.** Con `ka = 1`, dominio 4a, `h = a/4.5`:

- acero (contraste 30.7): **4.29 %** con restricción · 19.52 % sin ella
- aluminio (11.4): **3.99 %**
- tungsteno (66.9): **4.71 %**
- hueso (4.6): **13.16 %** con restricción · 54.61 % sin ella
- tejido casi fluido (1.08, `c_s = 30`): **5264 %** con restricción · 23.19 % sin ella
- tejido (`c_s = 200`): 1049 % con restricción · 111 % sin ella
- grasa (`c_s = 100`, contraste 0.92): 589 % con restricción · **353 % sin ella**

**La grasa agrega un tercer mecanismo, distinto de los otros dos.** Con `c_s = 100` la longitud de
onda de corte es `λ_s = c_s/f = 1.0 mm`, y el elemento mide 0.53 mm: **1.9 elementos por λ_s**. La
malla se elige con la longitud de onda del **fluido** (15 mm), y un sólido blando necesita su propia
malla más fina. Por eso a la grasa no le alcanza con apagar la restricción: el interior del sólido
está sub-resuelto.

Se agregó la métrica `elementos_por_lambda_corte` a la salida, y la UI avisa cuando baja de 4. Medido:
acero **53.6**, hueso **31.8** (o sea que el 13 % del hueso **no** es falta de resolución de corte, es
el modelo), tejido casi fluido **0.5**. Tres mecanismos distintos, tres diagnósticos distintos, y
ninguno se arregla subiendo el `ppw` del fluido.

## D48 — Consulta a Sim-Datos: la restricción de traslación media nula está mal planteada

Consulta del 2026-11-20 sobre (1) `c_s` de tejido blando en banda, (2) el parámetro correcto para
decidir "dispersor débil", (3) cómo reportar el error de un campo dispersado débil y (4) a partir de
qué contraste vale la suposición de desplazamiento medio nulo. La respuesta obliga a corregir tres
cosas de D47, y una de ellas es conceptual.

### 1. La premisa física era falsa

Se impuso `U_CM = 0` argumentando que "una onda plana no ejerce fuerza neta sobre un cuerpo real".
**Sim-Datos: eso confunde la fuerza armónica de primer orden con la fuerza media de radiación de
segundo orden.** Que la primera promedie cero en un ciclo no significa que su amplitud compleja sea
cero ni que el cuerpo no oscile. Para un sólido libre el balance es

    −ω²·M·U_CM = −∮ p_total·n dS

e imponer `U_CM = 0` **elimina la respuesta traslacional** salvo que la fuerza resultante sea también
cero.

Referencia directamente pertinente: **Feuillade y Clay (1999), *Anderson (1950) revisited*, JASA
106:553–564, §IV**, que demuestra que Faran y Anderson **contienen** el *rebound* traslacional y que
expandir alrededor de un origen fijo no lo elimina.

Criterio cuantitativo en el límite `k₀a ≪ 1` (esfera pequeña, fluido inviscido, respuesta dipolar),
de las ecs. 22 y 26 de Bruus (2012):

    v_partícula / v_incidente = 3 / (2g + 1),   con g = ρ_sólido / ρ_fluido

Con esa fórmula se entiende la tabla de D47 sin inventar nada: cuando `g ≈ 1` la traslación **acompaña
al fluido** y no debe anularse (tejido, g = 1.05 → 0.97); con acero (g = 7.8 → 0.18) la contribución
suprimida es chica. La inmovilidad aparece como límite de **gran inercia**, no de gran impedancia.

Consecuencia: **no existe un umbral en `Z_s/Z_0`**. El corte en 2.5 que se puso en D47 no es un
criterio físico; queda como **interruptor empírico de dos puntos medidos**, rotulado como tal, y el
arreglo es dejar de imponer la restricción resolviendo el sólido en campo dispersado (donde el
*rebound* aparece por construcción).

### 2. El parámetro de "dispersor débil" no es la impedancia

Para una esfera elástica los parámetros independientes son `ρ_s/ρ_0`, `c_p,s/c_0` y `c_s/c_0`; dos
materiales con el mismo producto `ρ·c_p` **no** dan la misma dispersión. El criterio publicado de Born
para dispersores fluidos (Sharma y Saha 2004, ecs. 4a–4b) es

    γ_ρ  = (ρ_s − ρ_0)/ρ_s ;  γ_κ = (κ_s − κ_0)/κ_0
    |γ_ρ| ≪ 1 ;  |γ_κ| ≪ 1 ;  x·|γ_ρ| < 1 ;  x·|γ_κ| < 1      (x = k₀a)

y **no** fija un porcentaje universal de contraste que garantice una precisión dada, ni se traslada
sin más a una esfera sólida con resonancias de corte. Rayleigh (`k₀a ≪ 1`) y Born (perturbación débil
del campo interno) son hipótesis distintas: un objeto chico puede tener contraste grande y dispersión
resonante.

### 3. Cómo reportar el error de un campo dispersado débil

Se reportan **las dos métricas, siempre, sobre los mismos puntos**, y no se elige una:

    E_inc = ‖p_s,num − p_s,ref‖∞ / |p_inc|        (precisión respecto de la excitación)
    E_sca = ‖p_s,num − p_s,ref‖∞ / ‖p_s,ref‖∞     (precisión respecto de la señal a resolver)

Es una recomendación de validación, **no una norma acústica**. No hay umbral citable para descartar
`E_sca`; cuando la referencia se anula el error relativo queda indefinido y cerca de un nulo es muy
sensible (Lavia et al., TetraScatt, fig. 4, atribuyen picos del error relativo a nulos de la
respuesta). O sea: el estado `sin_referencia` que se puso en D47 por debajo del 5 % de amplitud es
defendible, pero **no reemplaza** a `E_inc`: se muestran las dos.

Sobre secciones eficaces: son observables estándar (`dσ/dΩ = |F|²`, `σ_sca = ∫|F|²dΩ`), pero pierden
fase y detalle angular, y un patrón normalizado por su propio máximo puede ocultar un error
multiplicativo enorme. Para validar Faran: **campo complejo + σ_sca + patrón angular con escala
absoluta**, no `σ_sca` en lugar del campo.

### 4. `c_s` de tejido blando: no hay dato citable en la banda

- Las `c_s` de elastografía (1–5 m/s, μ ~ 1–20 kPa) son de **cientos de Hz**, no de la portadora
  ultrasónica; Chen et al. (2009) trabajan con ondas de corte de cientos de Hz.
- Hay mediciones ultrasónicas de tejido bovino (Madsen et al. 1983, JASA 74:1346, 2–14 MHz) pero
  **fuera de la banda** y con fuerte dispersión.
- En 100 kHz–1 MHz hay un dato biológico parcial: Feng, Li y Yun (2023, Nat. Commun. 14:4949), ondas
  de Rayleigh en piel humana *in vivo* hasta 500 kHz; a 500 kHz implican `c_R ≈ 200 m/s`, pero es una
  **onda superficial en piel estratificada**, no `c_s` volumétrica de músculo, grasa o hígado.
- Para tejido viscoelástico hace falta `μ*(ω)`, no sólo `μ` real: `k_s² = ρω²/μ*(ω)`.

**Conclusión para el simulador**: los materiales "tejido (casi fluido, c_s = 30)" y "grasa
(c_s = 100)" que están en la UI **no tienen fuente citable** en la banda. Quedan como *supuesto
explícito* rotulado en el tooltip, no como dato. Y no hay umbral universal de `c_s` para decidir
cuándo un modelo fluido es suficiente: hay que mirar la contribución de corte y su amortiguamiento.

### Lo que cambia en el repo

1. `docs/decisions.md` (este ADR): se corrige la justificación de D47. El umbral de contraste 2.5 deja
   de presentarse como criterio físico y pasa a interruptor empírico.
2. Los materiales blandos se rotulan en la UI como **ballpark sin fuente en banda**.
3. La métrica doble (`E_inc` y `E_sca`) ya está en la salida; el sello usa `E_sca` y el texto muestra
   las dos.
4. El paso de fondo —sólido en campo dispersado— deja de ser "una mejora" y pasa a ser **la
   corrección de un error conceptual**: la formulación actual suprime un término dipolar que Faran
   incluye.

## D49 — Sólido en campo dispersado: la formulación correcta es numéricamente peor

**Qué se hizo.** Se agregó `formulacion="dispersada"` a `coupled3d.resolver`: se resuelve
`u_s = u − u_i` en vez del desplazamiento total, y la incidente entra como **término de contraste**.
Es la respuesta a D48 (la restricción de traslación media nula suprime un dipolo físico).

**El detalle que hace que funcione: contraste, no residuo.** Restar el operador **completo**
`(K−ω²M)·u_i − b_s` deja un residuo O(h²) —medido: **0.35 %** de la tracción incidente, con los dos
términos valiendo 2.558e-7 y 2.549e-7— y ese residuo empuja un campo espurio de **~0.06 Pa**, del orden
del 18 % del campo real de un metal. La forma correcta es la lección de D41 trasladada a elasticidad:
aplicar sólo el **contraste del operador**, con el operador del fluido armado sobre la **misma malla y
la misma cuadratura**. Así, un cuerpo indistinguible del fluido da matrices **idénticas** y el término
es exactamente cero, sin depender de `h`.

**Validado en el límite, y esto sí es exacto.** Con `ρ` y `c_p` del agua y `c_s` bajando:

    c_s = 100 → |p_s| = 6.7e-2 Pa
    c_s =  30 → |p_s| = 2.6e-3 Pa
    c_s =  10 → |p_s| = 5.3e-4 Pa
    c_s =   3 → |p_s| = 1.2e-4 Pa      (incidente = 1 Pa)

Converge a cero como `c_s²`, **sin restringir nada**. Y el test de que las dos matrices coinciden da
diferencia exactamente 0.0. La formulación está bien.

**Y sin embargo, contra Faran es peor:**

- acero: **19.49 %** (razón 0.847) contra **4.29 %** de la total+restricción
- aluminio: **44.14 %** (razón 0.673) contra **3.99 %**

**Por qué (hipótesis con fundamento).** Para un metal el cuerpo casi no se mueve, así que
`u_s ≈ −u_i` con la magnitud de la incidente, y el término elástico `(K−ω²M)u_s` con la rigidez del
acero es enorme: la ecuación es una **resta de términos grandes**. Cualquier error de interpolación de
`u_i` en el interior se amplifica por `K` (~1e11). La formulación dispersada es correcta en el
continuo y **mal condicionada en esta malla**; la total+restricción es incorrecta en el continuo pero
bien condicionada.

**Consecuencia honesta, y hay que decirla: el 4.29 % estaba en parte de suerte.** No es que la
restricción acierte: suprime el término dipolar que Faran incluye, y en esta malla el error que
introduce se cancela parcialmente con el error de discretización del resto. Un número bueno por
cancelación no es una validación, y por eso el criterio `volumen_esfera_elastica` queda **calificado**
en D48/D49: vale para metales con la formulación total y **no** para materiales livianos.

**Decisión de default.** Se deja `formulacion="total"` con `fijar_rigidos=True` como default —es la
que reproduce 4.29 % y la que está sellada— y la dispersada se expone como opción documentada. Cambiar
el default a la formulación correcta pero 5-10× peor contra la referencia habría empeorado la UI sin
explicar nada; lo que corresponde es entender la amplificación antes de mover el default.

**Siguiente paso, concreto.** Antes de refinar a ciegas: (a) comparar la **traslación** que resuelve la
formulación dispersada contra la ley analítica `v_partícula/v_inc = 3/(2g+1)` en el límite `k₀a ≪ 1`
(Bruus 2012, ecs. 22 y 26) — si la traslación sale bien y el campo mal, el problema es el resto del
campo; si sale mal, es el acoplamiento con el fluido en un contorno que se mueve; (b) recién después
decidir si se refina el sólido, el fluido o se cambia de formulación.

#### D49 — apéndice: la corrida vieja, que confirma el diagnóstico

Quedó registrada una corrida del barrido anterior (con el operador **completo** restado, no el
contraste). Su límite transparente **no converge**:

    c_s = 100 → 1.22e-1 Pa
    c_s =  30 → 6.42e-2 Pa
    c_s =  10 → 5.72e-2 Pa      ← se estanca

Contra el operador de contraste, a `c_s = 10`: **5.3e-4 Pa**. Misma malla, mismo caso, **100× menos**.
O sea que el residuo O(h²) del operador completo no era un detalle de implementación: era el piso del
error, y por eso el límite no bajaba al bajar `c_s`. Es la evidencia de que el cambio de "restar el
operador completo" a "aplicar el contraste" era el paso correcto, más allá de que el número contra
Faran haya empeorado por otra razón.

Y el `fijar_rigidos` con default booleano era una trampa: pedir `formulacion="dispersada"` sin tocar
nada reventaba con un error que hablaba de la restricción. Ahora su default es `None` y se resuelve
solo según la formulación (`True` en total, `False` en dispersada); el error queda sólo para quien
pida explícitamente las dos cosas contradictorias.

#### D49 — apéndice 2: la tabla completa de la formulación dispersada, y el patrón

Barrido completo a `ka = 1`, dominio 4a, `h = a/4.5`. Se compara la formulación dispersada contra la
total+restricción, que es la que hoy está sellada:

| material | contraste `ρc_p` | dispersada | total+restricción |
|---|---|---|---|
| tungsteno | 66.9 | **10.70 %** (razón 0.914) | 4.71 % |
| acero | 30.7 | **19.49 %** (0.847) | 4.29 % |
| aluminio | 11.4 | **44.14 %** (0.673) | 3.99 % |
| hueso | 4.6 | **54.53 %** (0.628) | 13.16 % |

**El patrón: el error de la dispersada crece al bajar el contraste**, en el orden exacto de la
impedancia. Eso es coherente con la explicación de D49 y la precisa: el término de contraste
`(K_sólido − K_fluido)u_i` es una **diferencia de operadores grandes**, y el error de interpolación de
`u_i` en el interior no se hace pequeño al bajar el contraste. El campo dispersado verdadero sí se
achica con el contraste, así que el error espurio —que no se achica— domina cada vez más. De ahí que
tungsteno (contraste 66.9) sea el mejor caso de la dispersada y hueso (4.6) el peor.

**Lo que esto NO dice.** No dice que la dispersada esté mal implementada: su límite transparente
converge como `c_s²` y las dos matrices coinciden exactamente cuando los materiales son iguales. Dice
que **esta malla no resuelve el término de contraste** para contrastes moderados, y que el número
"bueno" de la total+restricción se sostiene por una cancelación entre el dipolo suprimido y el error
de discretización. Ninguna de las dos es una validación limpia del caso elástico, y así queda escrito.

**Nota de método.** La corrida de tests que se lanzó en paralelo con este barrido **se cortó por
timeout** (900 s, código 124): los dos jobs se peleaban la CPU, no era un fallo del código. Se repitió
en solitario para tener el veredicto.

**Insumo para el siguiente paso**: la ley `3/(2g+1)` (Bruus 2012) permite predecir la traslación
dipolar y separarla del resto del campo. Con la tabla de arriba, la hipótesis a testear es que la
traslación sale bien y el error está en el resto del campo (que es lo que el contraste moderado
amplifica). Eso decide si hay que refinar el sólido, el fluido, o cambiar de formulación.

## D50 — La traslación dipolar valida `3/(2g+1)` al 1-2 %, y la restricción la anula

**El test.** En el límite de Rayleigh (`k₀a ≪ 1`) la amplitud de oscilación del cuerpo tiene forma
cerrada: `|u_cm|/|u_i| = 3/(2g+1)` con `g = ρ_sólido/ρ_fluido` (Bruus 2012, ecs. 22 y 26). Es una
predicción **analítica, exacta en ese límite e independiente de la malla** — no hay nada que ajustar.
Se corre a `k₀a = 0.1` (`a = 0.239 mm`, dominio 4a, `h = a/4.5`) y se extrae el desplazamiento de
cuerpo rígido proyectando sobre el modo de traslación **con la matriz de masa**
(`a_d = t_dᵀMu / t_dᵀMt_d`), no promediando nodos: la malla no está uniformemente repartida.

| material | g | `3/(2g+1)` | total+restricción | total | dispersada |
|---|---|---|---|---|---|
| acero | 7.80 | 0.1807 | **0.0000** | 0.1772 | 0.1772 |
| aluminio | 2.70 | 0.4688 | **0.0000** | 0.4627 | 0.4627 |
| hueso | 1.95 | 0.6122 | **0.0001** | 0.6064 | 0.6064 |
| tejido `c_s=30` | 1.05 | 0.9677 | **10.09** | 0.9670 | 0.9670 |

**Tres lecturas, las tres importantes.**

1. **Las dos formulaciones correctas dan la traslación al 1-2 %.** El acoplamiento fluido-sólido y la
   masa añadida están bien. El error del 19-44 % a `ka = 1` **no viene de acá**.

2. **La restricción de traslación media nula no "aproxima" la traslación: la pone en cero.** Donde la
   física pide 0.18 (acero), 0.47 (aluminio) o 0.61 (hueso), devuelve 0.0000 / 0.0000 / 0.0001. Es la
   confirmación cuantitativa de D48: no es una condición de gauge inofensiva, es suprimir un grado de
   libertad que oscila de verdad. Y en el caso casi transparente **amplifica**: 10.09 contra 0.9677,
   un factor 10 de dispersión falsa.

3. **La formulación dispersada reproduce la traslación exactamente** (idéntica a la total sin
   restricción, hasta el dígito que se imprime). El operador de contraste no introduce sesgo en la
   respuesta dipolar, que es justo lo que se le pedía.

**Consecuencia para el siguiente paso.** Como la traslación está bien, el 19-44 % a `ka = 1` es un
problema de **resolución del resto del campo**: la estructura angular/radial del dispersado para un
objeto de tamaño comparable a la longitud de onda. Eso se ataca refinando (y midiendo si converge), no
cambiando la formulación. Se lanza la convergencia en `h` para el acero con la dispersada.

**Y una consecuencia para el sello.** El criterio `volumen_esfera_elastica` sigue apoyado en la
formulación total+restricción, cuya traslación es 0 por construcción. A `ka = 0.1` eso se ve como un
error del 100 % en la traslación y a `ka = 1` se disimula en el 4.29 %. Queda calificado: **el sello no
distingue "acertó" de "canceló"**, y eso es una limitación del criterio, no del método.

#### D49 — apéndice 3: **no es resolución de malla**, y la métrica de máximo es frágil

Convergencia en `h` de la formulación dispersada, acero a `ka = 1`, dominio 4a:

| `h` | nodos | error |
|---|---|---|
| a/4.5 | 21 398 | 19.49 % |
| a/6 | 47 574 | **18.84 %** |

Doblar el número de nodos baja el error **0.65 puntos** (orden aparente ≈ 0.14). Eso **descarta** la
hipótesis que había escrito en el apéndice anterior ("es resolución del resto del campo"): refinar no
lo arregla. El error es **sistemático de la formulación o de la métrica**, no de discretización.

Dos sospechosos, y el segundo es el que ahora me parece más probable:

1. **El término de contraste está mal en algo que no depende de `h`.** La identidad
   `(K_f − ω²M_f)u_i = −∮p_i(n·v)` es exacta en el continuo para cualquier `v` (se verificó con la
   integración por partes y el caso transparente converge como `c_s²`), así que si hay un error acá
   tiene que ser en el acoplamiento de superficie, no en el volumen.

2. **La métrica es frágil.** El error se mide como el **máximo sobre la curva angular** de la
   diferencia contra Faran. La esfera del FEM es **facetada**: con `h = a/6` hay ~6 facetas por radio,
   o sea del orden de 10° de paso angular, y eso **desplaza los lóbulos** en lugar de atenuarlos. Un
   lóbulo corrido unos grados da un error grande en el máximo aunque la magnitud global esté bien — y
   ese error no baja al refinar `h` tan rápido como baja el error de amplitud. Es el mismo tipo de
   trampa que el promedio de D43, al revés: el máximo premia el detalle más frágil.

**Lo que hay que hacer antes de tocar el código**: medir el error en **puntos fijos** (o en norma L2
sobre la curva, o en el coeficiente `C_0` de la serie) además del máximo. Si el error de los puntos
fijos baja con `h` y el máximo no, el problema es la métrica y la implementación está bien. Si ninguno
baja, es la formulación. Esa distinción decide el trabajo siguiente y todavía no la tengo.

**Mientras tanto**, la UI muestra las dos formulaciones con su número medido y la advertencia de que
la correcta da peor: no se esconde el resultado incómodo ni se presenta el 4.29 % como si fuera una
validación limpia.

## D51 — La ventana de validez de la esfera elástica es estrecha, y el barrido en `ka` lo muestra

**Qué se agregó.** `dashboard/elastic_sweep.py` y su vista: resuelve la esfera en varios `ka` con una
formulación y grafica el error contra Faran con la tolerancia del criterio como línea horizontal. Un
número suelto dice cuánto se equivoca; el barrido dice **en qué régimen vale**.

**Barrido preciso** (acero, dominio 4a, `h = a/4.5`, 61 puntos angulares):

| `ka` | total: error máx. | total: L2 | sello total | dispersada: máx. | dispersada: L2 |
|---|---|---|---|---|---|
| 0.75 | 6.84 % | 5.89 % | fuera | 19.50 % | 20.14 % |
| 1.17 | **4.88 %** | 3.09 % | **validado** | 20.24 % | 17.88 % |
| 1.58 | **4.78 %** | 3.12 % | **validado** | 14.54 % | 13.66 % |
| 2.00 | 9.20 % | 5.45 % | fuera | **7.43 %** | 7.69 % |

**Tres consecuencias, y la primera es la que hay que decir en voz alta.**

1. **El criterio del 6 % no vale en todo el rango: vale en `ka` ≈ 1.1–1.6.** En `ka = 0.75` la total da
   6.84 % y en `ka = 2` da 9.20 %, las dos **fuera**. O sea que decir "la esfera elástica está validada
   al 4.29 %" sin decir en qué ventana es exactamente el tipo de afirmación que este proyecto no hace.
   Se corrige la `fuente` del criterio para que nombre la ventana medida.

2. **Las dos formulaciones se cruzan.** La total gana en la ventana y en `ka` bajo; a `ka = 2` la
   dispersada es mejor (7.43 % contra 9.20 %). No hay una "mejor" en todo el rango: la formulación
   total es un instrumento con dominio de validez, no una aproximación universal. Es coherente con
   D48 —suprime un término cuyo peso cambia con `ka`— y con D49 —la dispersada es una resta de
   términos grandes que pesa más donde el campo dispersado es chico—.

3. **El modo rápido NO sirve para comparar formulaciones, y el barrido lo delató.** Con dominio 2.5a
   (5 puntos, acero): total 17.7 / 10.8 / 7.9 / 10.0 / 17.1 % y dispersada 15.2 / 19.6 / 20.1 / 11.3 /
   12.9 %, y ahí la dispersada parecía **mejor en `ka` bajo**. En el preciso ocurre al revés (6.84 %
   contra 19.50 % en `ka = 0.75`). La explicación es que con dominio 2.5a el error del ABC domina y
   las dos formulaciones miden ese error, no su diferencia. Es la lección de D43 otra vez, en otra
   escala: **no se comparan dos modelos donde manda un tercer error**. Por eso la UI muestra el modo y
   el aviso de dominio junto al barrido.

**Lo que queda anotado como límite del criterio, no como bug.** El criterio `volumen_esfera_elastica`
es un criterio de una ventana; la fuente ahora lo dice. Lo que falta para tener un criterio de todo el
rango es entender el error de la total fuera de la ventana —y ahí la hipótesis sigue siendo el dipolo
suprimido (D48), que pesa más justo donde la total falla: `ka` bajo y `ka` alto.

## D52 — El "error relativo de R" de la escena comparaba fasores con orígenes distintos

**Qué pasaba.** La escena reportaba `error_R_pct` en el rango de **100 % a 200 %** mientras el error
de amplitud era del **0.2 %** y el sello decía VALIDADO. Un número que se lee como "el modelo se
equivoca al doble" al lado de un "validado" es exactamente la clase de dato que no se puede dejar
pasar.

**La causa, medida.** Dos cosas a la vez:

1. **`R` se comparaba como fasor y `T` como módulo.** La línea era
   `e_r = |r_num − r_an| / |r_an|` para R y `e_t = | |t_num| − |t_an| | / |t_an|` para T. No era una
   decisión: era un descuido, y hacía que los dos números no fueran comparables entre sí.

2. **Los dos fasores tenían orígenes distintos.** `r_num` sale de ajustar `a·e^{ikx} + b·e^{−ikx}` con
   la `x` **global**, así que `b/a` está referido a `x = 0`; `r_an` (matriz de transferencia) está
   referido a la **primera interfaz** (`x = t₁`). Entre los dos hay un giro `e^{2ik·t₁}` que es
   propagación, no error.

**La verificación.** Si el desvío era sólo ese giro, tenía que cumplirse
`|ΔR|/|R| = 2·sin(k·t₁)`. Medido contra predicho, con `t₁ = 20 mm` de tejido:

| f | `|ΔR|/|R|` medido | `2·sin(k·t₁)` |
|---|---|---|
| 100 kHz | 1.917 | 1.907 |
| 200 kHz | 1.072 | 1.149 |
| 300 kHz | 1.259 | 1.215 |
| 400 kHz | 1.839 | 1.881 |

Cuatro de cuatro, dentro del 7 %. Es el giro de propagación, no un error del solver.

**La corrección.** Se lleva `r` a la referencia de la interfaz (`r·e^{−2ik t₁}`) antes de comparar
fasores, y se reportan **dos** números con nombre propio:

- `error_R_pct`: sobre el **módulo**, igual que `error_T_pct` — es el número comparable y el que se lee
  de un vistazo. Ahora da **0.4–2.0 %** (antes 107–199 %).
- `error_R_fasor_pct`: el fasor ya alineado. Da **6–8 %**, que es del orden del error de **fase** de
  la malla con `ppw = 8` (~7°), no del error de amplitud. Se reporta porque la fase es información
  real, pero no se mezcla con la amplitud.
- `desfase_referencia_R_deg`: el giro que se sacó (935° a 100 kHz, 3740° a 400 kHz), a la vista, para
  que quien compare fases sepa qué se corrigió.

**Por qué importa para el proyecto.** El criterio de la Fase 1 juzga "dentro del 2 % **en amplitud**" y
usa `error_R_abs`, que siempre estuvo bien: el sello nunca estuvo mal. Lo que estaba mal era un número
**al lado** del sello que lo contradecía. Un reporte que se contradice a sí mismo se deja de mirar, y
ahí se pierde el valor de tener el registro.

## D53 — El criterio de pitch `λ/4` medido en la UI: 0.28 λ para el 2 %

**Qué se agregó.** Un barrido de apertura en la pestaña `barridos`: con el ancho de apertura fijo
(60 mm), variar la cantidad de elementos es variar el **pitch**; cada punto es un campo de verdad y se
compara el FWHM lateral contra una apertura cuasi-continua (`n = 400`) del mismo ancho. Es D25/D26
llevado de script a instrumento.

**Medido** (60 mm, foco 90 mm, 250 kHz, agua, grilla 120×120):

| pitch | FWHM lateral | desvío vs continuo | grataje |
|---|---|---|---|
| 1.25 λ | 7.95 mm | **10.8 %** | 53.3° |
| 0.45 λ | 7.34 mm | 2.3 % | — |
| 0.28 λ | 7.30 mm | **1.7 %** | — |
| 0.20 λ | 7.29 mm | 1.6 % | — |
| 0.16 λ | 7.28 mm | 1.5 % | — |
| continuo (n=400) | 7.18 mm | — | — |

El pitch más grueso que cumple el 2 % cae en **0.28 λ**, o sea **1.1× el criterio `λ/4`** de D25/D26.
Coincide, y con el matiz que corresponde: el criterio es un orden, no un borde exacto; entre 0.25 λ y
0.28 λ el desvío está en 1.7–1.8 %.

**Dos errores propios que el barrido delató, y los dos del tipo que "parece un resultado":**

1. **No barría.** La primera versión pasaba `n` y `d` al payload de `calcular`, que espera `n_elements`
   y `aperture_mm`. No fallaba: las ignoraba y devolvía **el mismo campo cinco veces** (FWHM 9.12 mm y
   sidelobe −9.24 dB idénticos con 8 y con 64 elementos). Lo peor de ese modo de fallar es que el
   gráfico sale prolijo y plano, y sin un test que exija que el número **cambie** no se nota.
2. **La firma de la cola.** `jobs` llama `fn(trabajo, datos)` y el barrido esperaba `(datos, trabajo)`:
   moría con `'Trabajo' object has no attribute 'get'`, que no dice nada sobre la causa.

Ambos quedan fijados en `tests/test_barrido_apertura.py`: uno exige que el FWHM **cambie** con el
pitch, el otro que el punto del 2 % caiga entre 0.15 λ y 0.60 λ (o sea, ni un orden de magnitud
arriba ni abajo).

**Por qué va en el repo.** El criterio de pitch estaba en `docs/` y en un script; ahora se puede
reproducir desde la UI y el resultado trae su número al lado. Es la diferencia entre "D25 dice λ/4" y
"con estos parámetros, 0.28 λ da 1.7 %".

#### D53 — apéndice: el pitch que cumple el 2 % se mide en el medio elegido

El barrido usaba `λ = 1500/f` fija (agua) para normalizar el pitch. Eso mide el criterio de otro
medio: en hueso `λ` es 2.33× más corta y el pitch que cumple `λ/4` baja con ella. Ahora se toma
`λ` del material resuelto (el `lambda_mm` que devuelve el propio campo), y el barrido deja elegir
medio.

Medido, 60 mm y 250 kHz, 4 puntos:

| medio | `pitch/λ` para n = 8, 22, 36, 64 | pitch más grueso dentro del 2 % |
|---|---|---|
| agua | 1.25, 0.37, 0.22, 0.16 | **0.370 λ** |
| hueso | 0.536, 0.159, 0.095, 0.067 | **0.159 λ** |

El número en `λ` baja en hueso porque `λ` baja: en milímetros el pitch que cumple es el mismo
(60 mm / n), lo que cambia es contra qué se compara. Es la comprobación de que normalizar por `λ` del
medio es lo correcto y no un detalle de presentación.

## D54 — El tablero pasa a un solo problema; el FEM de ese problema todavía no cierra

**Qué se reordenó, y por qué.** El tablero tenía ocho pestañas que resolvían **problemas distintos**
(una apertura con baffle, ondas planas sobre capas, un cilindro dispersor, un corte de cóclea, una caja
3D). Comparar métodos así es imposible: no hay dos pestañas que resuelvan lo mismo. La estructura pasa
a: **BEM · fuentes** (el problema), **FEM · fuentes**, **escena capas**, **dispersor** y **3D**, con
cóclea, barridos, esfera y criterios marcados **para el futuro**.

Además la **frecuencia se sincroniza** entre campo, escena y 3D (antes campo arrancaba en 250 kHz y las
otras dos en 100 kHz: comparar obligaba a igualarlas a mano). Sincronizar **no** relanza corridas.

**Lo que está verificado.** El BEM de referencia (`dashboard.field.evaluate`) coincide con la solución
exacta de un monopolo 2D, `H₀(k·r)`:

| r (mm) | 30 | 60 | 100 | 140 |
|---|---|---|---|---|
| razón de módulos | 1.000 | 1.000 | 1.000 | 1.000 |
| diferencia de fase | 0.0° | 0.0° | 0.0° | 0.0° |

La referencia es buena. Eso convierte cualquier diferencia en un problema del FEM, no en una
ambigüedad sobre quién tiene razón.

**Lo que NO está verificado.** `dashboard/fuentes_fem.py` resuelve el mismo problema con P1 + ABC-1
(Neumann en el baffle con `∂p/∂n = −iωρv` en las aristas de cada fuente). Contra el BEM **no
reproduce el campo**: con una fuente de 2 mm en un dominio de 200×200 mm², la razón de amplitudes va de
**0.31 a 0.92** según el punto y las fases difieren entre **45° y 165°**. No es un factor global: un
desfase constante se corregiría con una constante, y acá el cociente complejo cambia punto a punto.

Sospechosos, en orden de probabilidad:

1. **La forma débil de `solve_scattered` no es la de un problema radiante.** Está escrita para un
   *dispersor*: impone `∂p_s/∂n = −∂p_i/∂n` e integra el dato con una convención de signo que se probó
   con `±` y ninguna de las dos alinea el campo. Probablemente haya un término de frontera de más (el
   `p_i` que se resta) que en un problema sin incidente no corresponde.
2. **El ABC-1 en un contorno rectangular abierto por el baffle.** El orden 2 exige contorno cerrado y
   no aplica; el orden 1 mete un error que no se midió por separado.
3. **La malla no resuelve el ancho de cada fuente**: quedaron 3 aristas para 2 mm.

**Decisión.** El módulo queda en el repo (la malla, la clasificación de aristas del baffle, el recorte
de la fuente por arista y la comparación contra el BEM son correctos y reutilizables) pero **la UI no
lo expone**: la pestaña se muestra **deshabilitada y con la etiqueta «en desarrollo»**, y el tooltip
dice exactamente qué falta. Mostrarla como si funcionara sería el peor de los tres caminos (esconderla
hace creer que no existe; mostrarla hace creer que funciona).

### Intento con el problema radiante — tampoco cierra

Se reescribió el ensamblado como problema **radiante**, que era el sospechoso nº1: `assemble(...,
p_inc=None)` da `A` y `b = 0` y deja el baffle como Neumann natural, y el dato de las fuentes entra sólo
en el lado derecho con

    b[nodo] += (1/ρ)·g·L/2 ,   g = ∂p/∂n = +i·ω·ρ·v

de modo que `(1/ρ)·g = i·ω·v` y el dato no depende de ρ. Nada de `solve_scattered`.

Medido con el arreglo de 16 fuentes (60 mm, foco 90 mm, 250 kHz, ppw 8, 32 148 nodos):
el error máximo es **170.8 %**, el L2 **123.3 %** y la razón de amplitudes **0.986**. La amplitud ya
está bien, la fase no. Pero el cociente complejo punto a punto **no es constante**: va de **0.11 a 9.61**
en módulo con **61° de dispersión** de fase. Eso descarta la hipótesis de convención —una convención da
un cociente constante— y deja el problema en el campo mismo.

Quedan dos sospechosos, y los dos son de discretización, no de formulación:

1. El **ABC de orden 1** en un contorno rectangular abierto por el baffle. No se aisló su error: habría
   que medir el mismo problema con dominio creciente y ver si el error baja como 1/(kR).
2. La **malla no resuelve el ancho de cada fuente** en el sentido correcto: con ppw 8 cada elemento de
   3.75 mm son 5 aristas, pero el dato se aplica con regla del punto medio sobre aristas, así que la
   velocidad queda muestreada en 5 saltos y eso inyecta armónicos espaciales cerca del baffle.

El test que hay que hacer antes de tocar nada más es el del **monopolo en dominio circular con ABC de
orden 2**, donde el contorno sí es cerrado y el ABC está validado en otros tests del repo. Si ahí el FEM
reproduce `H₀(k·r)`, el ensamblado radiante está bien y el problema es el contorno abierto. Si no,
el problema está en cómo se arma el dato.

**Estado**: `fuentes_fem.py` tiene el ensamblado radiante (mejor que el anterior, que reusaba el
dispersor) y la UI sigue sin exponerlo. D54 abre con esto, no con una promesa.

## D55 — STL entra de a una pieza; la cadena de capas se malla como geometría analítica

**El pedido**: cargar archivos STL para montar las estructuras a simular, y acoplar espacios para que el
ultrasonido atraviese piel, hueso y cóclea.

**La parte que funciona.** `stl.py` lee binario y ASCII y **verifica antes de mallar**: que cada arista
esté en exactamente dos triángulos y que el volumen por divergencia sea positivo. Un STL con un agujero
o con las caras dadas vuelta no malla a ciegas — se rechaza con el motivo escrito. Sin eso, un STL roto
produce un número y el número parece válido.

`meshing_volume.mallar_stl` convierte una pieza cerrada en tetraedros: la placa de 1 cm da **1000.0 mm³
exactos** con 702 tetraedros.

**La parte que no.** Dos STL que se tocan no se pueden mallar juntos, y no es un detalle de
implementación:

- `occ.importShapes` no acepta STL (sólo STEP/IGES/BREP). Hay que `merge` + `classifySurfaces` +
  `createGeometry`, y eso deja la geometría en el kernel **geo**, no en **occ**.
- Con las dos piezas en el kernel geo, `occ.fragment` falla («Unknown OpenCASCADE entity»).
- Sin fragment, las caras coincidentes están duplicadas y el mallador 3D se cae:
  «PLC Error: A segment and a facet intersect».
- Fusionar los nodos después de mallar la frontera (`removeDuplicateNodes` + `generate(3)`) rompe la
  correspondencia con la geometría: «Could not recover boundary mesh: error 2».

**La decisión**: para la cadena de capas —que es el caso piel → hueso → cóclea— se malla con
`mallar_capas`, que construye cada capa como **caja analítica** y las fragmenta con occ. Ahí la interfaz
es **una sola superficie** y la malla sale conforme por construcción: medido en piel 2 mm + hueso 5 mm +
líquido 4 mm sobre 20×20 mm, **4400.0 mm³** de volumen total, **3028 nodos y 0 duplicados**, y cada capa
cayendo en su rango de `x` (piel 0.21–1.79, hueso 2.17–6.88, líquido 7.17–10.83 mm).

`mallar_stl` acepta entonces **una pieza por malla** y se niega a más con un mensaje que dice por qué y
qué usar en su lugar. Es más útil una limitación explícita que una malla que parece buena y tiene la
interfaz rota: el acople entre regiones depende de que los nodos de la interfaz sean **los mismos**, y
una malla no conforme lo rompe en silencio.

**Lo que falta para STL multicapa**: mallar cada pieza con interfaz compartida requiere o bien geometry
kernel occ nativo desde el STL (no soportado por gmsh para STL) o bien acople con nodos no coincidentes
(mortar). Lo segundo es lo correcto numéricamente y es un desarrollo propio, no un ajuste. Queda como el
siguiente paso de U9.

## D56 — El acople de N regiones valida contra la teoría de capas (0.13 % en piel + hueso)

**Qué se hizo.** `fem/multicapa.py` resuelve **N regiones** (fluidos y sólidos elásticos mezclados) en una
malla conforme, sin tocar `solve_coupled` ni ninguna fórmula: reusa sus primitivas y agrega sólo el
ensamblado por grupo.

La observación que hace que esto sea chico: con malla conforme y P1, las interfaces entre dos regiones del
**mismo tipo** no necesitan bloque. Los dos lados comparten los mismos gdl, así que la continuidad de `u`
(sólido-sólido) y de `p` (fluido-fluido) es automática, el equilibrio de tracciones sale de la forma
débil al ensamblar cada región con **su** rigidez, y la continuidad de `(1/ρ)·∂p/∂n` sale sola porque cada
elemento lleva su `1/ρ`. La única interfaz que necesita bloque explícito es **fluido-sólido**, donde las
incógnitas son distintas.

**Lo medido** (agua | placa | agua, incidencia normal, ABC sólo en entrada y salida):

| caso | FEM | exacto | error |
|---|---|---|---|
| sin capa (onda plana sola) | 1.00184 | 1 | 0.18 % |
| capa con las propiedades del agua | 1.00187 | 1 | 0.19 % |
| placa de hueso 2.5 mm (λ/20) | 0.5892 | 0.5807 | **1.45 %** |
| piel 2 mm + hueso 5 mm (λ/18) | 0.9960 | 0.9947 | **0.13 %** |

Y el error **baja al refinar**: 5.33 % (λ/10) → **1.45 % (λ/20)** en la placa. Eso es lo que distingue
discretización de un error de formulación disfrazado.

**Dos errores propios, que quedan anotados porque los dos eran silenciosos.**

1. **ABC en las paredes laterales.** La primera versión la aplicaba a **toda** la frontera de fluido. Ahí
   la onda plana tiene `∂p/∂n = 0` y la ABC impone `∂p/∂n = i·k·p`: una condición espuria que refleja. El
   caso base daba |p| entre 0.82 y 1.37 y **no mejoraba al refinar** (λ/12: 0.819; λ/20: 0.848) — esa
   insensibilidad a `h` era la pista de que no era discretización. Con ABC sólo en entrada y salida da
   1.002. En un ducto con onda normal, las paredes son rígidas, no absorbentes.

2. **`1+R` no es la presión transmitida.** En la fórmula de matriz de transferencia del test devolvía
   `1+R`, que es la presión total en la **primera** interfaz. Para una capa de cuarto de onda eso da
   **1.814** donde el FEM da 0.581: es el transformador de impedancia, la presión en la cara de entrada
   sube, pero no es lo que sale del otro lado. Lo correcto es `p_N = M11·p0 + M12·u0`. La señal fue que
   `|t| > 1` sin pérdidas es imposible; el test tenía el error, no el solver.

**Lo que esto habilita**: la cadena piel → hueso → cóclea ya se puede resolver con materiales distintos por
región. Lo que sigue es conectarla al dashboard y contrastarla contra la transmisión craneal **medida**
(31 ± 7 % parietal a 500 kHz, Riis 2022), que es más exigente que la fórmula de capas porque incluye
geometría, curvatura y pérdidas.

## D57 — La cadena de tejidos entra al dashboard, con cuatro errores propios en el camino

**Qué se agregó.** Vista **«cadena piel→cóclea»**: lista de capas editable (medio + espesor), frecuencia,
elementos por longitud de onda, y un botón que resuelve. Devuelve el mapa de `|p|` en dB con las
**interfaces dibujadas y rotuladas**, el corte axial, un esquema de capas a escala (con el agua de entrada
y salida, que también es parte del problema) y cinco números: transmitido medido, exacto, error, espesor
total y λ en agua.

Medido en la UI: piel 2 mm + perilinfa 5 mm a 500 kHz → **0.9456 medido contra 0.9239 exacto, 2.35 %**.
Con hueso cortical la vista dice **«sin referencia»** y no muestra error, porque el hueso entra como
sólido elástico con cizalla y la fórmula de capas fluidas no lo incluye. Los materiales salen de
`docs/materiales-cadena.md` con su `fuente` y `verificado=false` a la vista (`ultrasound_sim/cadena.py`),
y la advertencia de que son parámetros de modelado va en el panel.

**Los cuatro errores, todos silenciosos, que vale la pena no repetir.**

1. **Clave doble en la persistencia.** `vistaGuardar(vista, …)` ya antepone `'usim.vista.'` y envuelve
   todo en `{fecha, resultado}`. Le pasé `'usim.vista.cadena'` y un objeto con su propio campo
   `resultado`: guardaba en `usim.vista.usim.vista.cadena` y la restauración leía otro nivel, así que la
   pestaña se abría **vacía y sin ningún error en consola**. Peor: como `init` es `async`, el fallo era un
   rechazo de promesa que no aparecía en `window.onerror`. Ahora hay un `.catch` que escribe el motivo en
   la línea de estado, y un test que fija que la clave se pasa pelada.

2. **`LinearTriInterpolator` con un fasor descarta la parte imaginaria.** matplotlib avisa con un
   `ComplexWarning` que se pierde entre los warnings de la suite; el mapa se estaba calculando como
   `|Re(p)|`, que no es la amplitud. Ahora se interpola parte real y parte imaginaria por separado
   (también en `fuentes_fem`, donde estaba el mismo defecto) y el test corre con
   `ComplexWarning` como **error**, no como aviso.

3. **El aviso del explorador de campo se colaba en la vista.** «El pico real está en 75.8 mm, no en
   F = 90 mm» es del foco de una apertura y en una cadena de capas no significa nada. Se oculta por CSS
   con `body.vista-cadena`, igual que se hizo con la escena, en vez de por `hidden` (que `app.js`
   re-muestra cuando llega su respuesta).

4. **`scaleanchor` sacaba las anotaciones del rango visible.** Con el div ancho, Plotly estiraba el eje x
   a −40…80 y recortaba el y; las etiquetas de capa y las líneas de interfaz quedaban **fuera** del área
   dibujada. Se veían «no aparecer» aunque estaban en el layout. Con rangos fijos en x e y se ven.

**Lo que queda de U9**: mover la estructura desde el dashboard (hoy las capas se editan como lista, no se
arrastran) y el STL multicapa, que necesita acople con nodos no coincidentes.

## D58 — La estructura se mueve con el mouse, y para eso tiene que ser finita

**Qué se agregó.** En la vista de cadena, la estructura puede ser una **pieza finita** (alto configurable)
en vez de capas que atraviesan el dominio. Con alto > 0 se dibuja un rectángulo con borde azul punteado
**sobre el mapa**, arrastrable en x e y; al soltar se vuelve a resolver. La posición también se fija con
dos deslizadores, y las interfaces de capa se **recortan al alto de la pieza** (si no, atraviesan todo el
dominio y parece una capa infinita, que es justo lo que no es).

**Por qué finita y no infinita.** Con capas infinitas, mover la estructura en x **no cambia la física**:
sólo cambia la fase de referencia. Un control que no hace nada es peor que no tenerlo. Con la pieza
finita, el campo la rodea y difracta en los bordes, y mover sí cambia lo que se mide.

**Medido en el navegador** (piel 2 mm + hueso 5 mm + perilinfa 5 mm, 500 kHz, alto 14 mm): arrastrar de
(12, 0) a (15.3, 9.2) mm cambia el transmitido de **1.1915 a 0.9282** y la malla pasa de 38 025 a 41 400
nodos. Por backend: x = 10 mm da t = 0.9436 y x = 18 mm da t = 0.9682; centrada en el eje t = 0.9443,
corrida 10 mm t = 0.9690 (la onda pasa casi sin obstáculo por el eje).

**Cuatro cosas que hubo que corregir, todas de presentación o de enganche.**

1. **El listener del arrastre no se registraba.** `$('mapaCadena').on('plotly_relayout', …)` estaba en
   `init()`, cuando el div **todavía no era un gráfico de Plotly** y no tenía `.on`; el arrastre movía el
   rectángulo y no recalculaba nada. Ahora se registra en `pintar()`, después del `Plotly.react`, y se
   quita el handler anterior para no acumularlos.
2. **La forma del evento no era la que esperaba.** Plotly emite `shapes[N].x0/x1/y0/y1`, con **N igual al
   índice del rectángulo** (después de las 2 líneas por interfaz, así que con 3 capas es 6, no un índice
   fijo). Filtrar por un prefijo armado a mano dejaba el arrastre mudo; ahora se acepta cualquier clave
   que empiece con `shapes` y se lee el shape **por nombre**.
3. **El rectángulo se confundía con el campo.** El sólido no tiene gdl de presión y salía negro, igual que
   un mínimo del colormap, así que parecía "zona sin campo" y no un objeto. Ahora va con relleno oscuro
   semitransparente y borde azul, `layer: above`.
4. **No había ninguna señal de que se arrastra.** Se agrega `cursor: grab` sobre el mapa cuando la pieza
   es finita, y el esquema de capas dice «pieza finita de N mm».

**Lo que esto cierra**: el control de posicionamiento de la estructura que pedía el objetivo, para
estructuras de capas. Mover una **geometría STL** con la misma interacción es el paso que sigue y depende
de tener el acople de STL multicapa (nodos no coincidentes), que es donde está el trabajo de fondo.

## D59 — El arreglo de fuentes y el transductor entran desde archivo (U10)

**Qué se agregó.** En el panel de la vista de campo:

* **Transductor del catálogo** (`data/transducers.toml`, que existía desde hace tiempo y la UI no usaba):
  al elegirlo aplica frecuencia, cantidad de elementos, apertura y foco, y **dice con qué está
  simulando**: `100 kHz · 32 elementos · apertura 60 mm · paso 0.125 λ · f/0.67 · sin lóbulo de rejilla`.
  El paso en λ se recalcula con el **medio elegido**, no con el agua: el mismo transductor da 0.125 λ en
  agua y **0.054 λ en hueso**, y ese factor es justo el que decide si hay lóbulo de rejilla (D53).
* **Cargar arreglo (CSV/JSON)** y **exportar**. El CSV lleva cabecera `y_mm,ancho_mm,amp,fase_rad`; las
  unidades están en el nombre de la columna a propósito.

**Lo que hace el servidor antes de aceptar un archivo** (`ultrasound_sim/arreglos.py`):

1. Convierte de mm a metros **diciéndolo**, en vez de asumir metros (32 filas en metros serían 32 km).
2. Rechaza elementos que se solapan: `|Δy| ≥ (ancho_i + ancho_j)/2`, nombrando **cuáles**.
3. Rechaza `ancho ≤ 0` y `amp ≤ 0` (una amplitud negativa es una fase de π disfrazada).
4. Si el archivo no pasa, **no se toca lo que está en pantalla**: el aviso sale en rojo con el motivo.

**Medido en el navegador**: cargar un CSV de 5 elementos aplica `−7.50 | 3.75 | 1.00 | 0.0` en la primera
fila y el contador pasa a «5 fuentes · apertura 18.8 mm». Cargar un archivo con dos elementos a 1 mm con
anchos de 4 mm da `los elementos 1 y 2 se solapan: están a 1.000 mm y necesitan al menos 4.000 mm` y
deja las 5 fuentes anteriores intactas. Exportar devuelve «exportadas 5 fuentes a CSV» y el archivo que
se exporta se vuelve a cargar sin pérdida (ida y vuelta en el test).

**Decisión sobre el nombre de las columnas**: el exportador escribe `y_mm,ancho_mm,amp,fase_rad` y el
cargador acepta además `y, width, amp, phase` en inglés, porque no tiene sentido obligar a renombrar un
archivo para volver a cargarlo.

**Lo que falta de U10**: no hay todavía un botón para **editar el catálogo** desde la UI (se edita el TOML
a mano, respetando que `source` no puede estar vacío), y el formato no lleva la geometría del transductor
curvo: el modelo sigue siendo un arreglo lineal, que es su corte 2D.
