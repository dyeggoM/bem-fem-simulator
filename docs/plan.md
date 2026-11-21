# Plan — simulador de ultrasonido en cóclea (BEM + FEM)

Estado: activo. Este documento manda. Lo que cambie se anota también en `docs/decisions.md`.

## Objetivo

Simulador de propagación de ultrasonido en tejido blando, hueso y fluidos cocleares, con **dos
métodos numéricos independientes** (FEM y BEM) para poder contrastarlos: la discrepancia entre ambos
es la señal de un error o de una suposición mal puesta.


## Estado entregado (al 2026-11-20)

Resumen de lo que está *cerrado con número medido*, no de lo que está escrito:

- **Fase 0 (1D)** — campo FEM vs analítico *0.399 %* L2; R/T de interfaz única contra impedancias.
- **Fase 1 (FEM 2D + ABC)** — R y T contra analítico *0.511 %*; |R| 0.93678 vs 0.93726; reflexión
  espuria de la ABC *−73 dB* (criterio: < −40 dB).
- **Fase 2 (hueso elástico acoplado)** — |T| *0.079 %*; balance |R|²+|T|² 0.999017 / 1.000075;
  Snell 52.53° vs 52.94°; resonancia de espesor +0.012 % (acero), +0.041 % (hueso).
- **Fase 3 (BEM 2D)** — cilindro rígido vs serie de Bessel *0.05 %* (ka=1) y *0.14 %* (ka=2);
  formulaciones `standard` / `burton_miller` / `chief` seleccionables.
- **Fase 4 (contraste BEM vs FEM)** — |FEM−BEM| *0.97 %* (ka=1), *0.59 %* (ka=2); |BEM−serie|
  0.05 %/0.14 %. Encontró el error de formulación D12, que ninguna prueba unitaria veía.
- **Fase 4B (fuentes, arreglos, focalización)** — 15 configuraciones (`out/focusing_study.csv`);
  los cuatro arreglos quedan dentro del 2 % de la apertura continua (aire 40k +1.3 %, 250k +0.4 %,
  500k +1.6 %, 1M +1.8 %); el Langevin monolítico no (+13.4 %). Criterio corregido: lo que manda es
  la cantidad de elementos, no el `pitch/λ` (λ/4 basta para el 2 %).
- **Fase 5a (corte de cóclea)** — malla converge +0.6 %, sin flujo espurio (1e-15), cáscara vs placa
  −1.52 %/−0.52 %. Dos criterios *no alcanzables y diagnosticados* (xfail con razón): la cápsula es
  una cavidad cerrada sin pérdidas, no hay un `T` que converja al de una placa.
- **Fase 5b-2 (membrana basilar)** — `Q` de la cavidad *287 → < 5*; absorción `P_memb/P_inc` 0.733;
  balance de energía *0.46 %*. Re-barrido de sensibilidad (D30): espesor 53 %, `c_p` 16 %,
  atenuación 9.6 %, `c_s` 8.2 %, `Z_r` 0.3 %.
- **Fase 6a (núcleo 3D)** — onda plana *3.53 %* a ppw 20, orden medido 1.88/1.94/1.97; ABC vs pared
  rígida 9.57 % vs 100 %; fuente puntual reproduce el Green 3D (amplitud 0.9903); pistón vs
  Rayleigh-Sommerfeld *0.0000 %*; R/T estratificado 3D *2.91 %* / *2.11 %* contra la matriz de
  transferencia. Techo: el costo, no la física.
- **UI (D27/D31)** — dashboard local con motor de campo, editor de fuentes con arrastre, panel de
  cóclea y métricas; 250 tests en verde (1 xfail documentado).

**Pendiente y por qué** (actualizado con la consulta a `sim-datos`, D40):

- **Dispersión penetrable 3D — CERRADA (D41).** `analytic/sphere.py` implementa la serie de **Anderson
  (1950)** validada por continuidad (salto de `p` y de flujo ≤ 1e-12 a ka = 0.5–5), por el límite
  rígido (`ρ→∞` da `∂ᵣp = 0` al 0.0001 %) y por la esfera casada (`Aₙ = 6e-24`). El FEM 3D con el
  término de contraste corregido la reproduce: **3.34 % (ka = 1)** y **2.73 % (ka = 2)** a ppw 12, con
  convergencia **5.72 % → 1.93 %** al ir de ppw 8 a 16. El bug de D32 era el cociente de densidades
  invertido en el contraste, más el fondo ausente en el término de rigidez.
  **En la UI**: la pestaña 3D tiene un tercer modo *esfera penetrable* (material y `ka` editables) que
  muestra el **campo dispersado** en el plano de propagación y sella contra Anderson
  (*2.94 % a ka = 1*, 9 261 nodos, 1.2 s). Se muestra el dispersado y no el total porque en el total
  manda la incidente y la perturbación desaparece en la escala de dB.
  **Faran (1951) — ENTREGADO (D45).** `analytic/faran.py`: la serie de la esfera **elástica**
  sumergida (dos potenciales P/S, tres condiciones de borde, sistema 3×3 por orden), validada sin
  tabla de referencia: condición unitaria a **1e-16**, teorema óptico al **1.0e-3 %**, límite de
  fluido cayendo como `c_s²` (**1.62e-3 → 4.04e-4 → 1.01e-4**) y aproximación monótona al rígido
  contra `sphere.py` (aluminio 5.2e-3 · acero 2.3e-3 · tungsteno 1.2e-3). 14 tests.
  **Fase 6b — primera mitad ENTREGADA (D46).** Malla esfera+cáscara de fluido con Gmsh (interfaz
  compartida nodo a nodo; volumen del sólido 3.51 % → 1.57 % al refinar) y elasticidad vectorial P1
  con acoplamiento fluido–sólido. Elasticidad **exacta**: modos rígidos en el núcleo de `K` a *4.1e-16*,
  energía de deformación uniforme exacta, masa `3ρV` exacta. El acoplamiento está bien por sus
  límites: sin contraste no dispersa y el límite rígido da *10.66 %* contra `sphere.py`, que es el
  *10.35 %* que da en el mismo dominio el Helmholtz ya validado (D41) — o sea, el dominio y el ABC.
  **Contra Faran: 9.2 % a `ka = 1`** (era 22 %), con razón de amplitudes *0.997*. El 22 % inicial se
  debía a los **seis modos de cuerpo rígido** del sólido libre: no tienen rigidez elástica y su
  excitación espuria por una fuerza neta que la malla no cancela contaminaba todo el campo. Se
  arregló con una restricción de Lagrange de traslación media nula, cuyos multiplicadores **son** la
  fuerza espuria y se devuelven como diagnóstico (*7.6e-8 N*). Convergencia al refinar: *10.20 %*
  (h=a/3.5) → *9.70 %* (a/4.5) → *9.21 %* (a/5.5).
  **El piso, descompuesto y bajado (D46, 3.er apéndice).** Fluido solo con datos exactos, una cosa
  por vez: el **dominio** se llevaba ~5 de los 9 puntos (`2.5a` → *9.30 %*, `4a` → *4.06 %*, `6a` →
  *3.96 %* por 12× el costo), y el resto es malla (`h=a/3` *10.24 %* → `a/4.5` *4.06 %* → `a/6`
  *3.03 %*, con el orden cayendo a 1 por el **facetado**). El acoplado completo en el dominio
  recomendado da **4.21 %** a `h=a/4.5` (32 s) y *3.01 %* a `a/6` (211 s): **el sólido aporta 0.15
  puntos**. `sphere_shell_tet_mesh` ahora toma `4·radius` y `radius/4.5` por defecto y **avisa** si el
  radio exterior baja de `3·radius`.
  Lo que queda sobre la mesa: una condición **Dirichlet-to-Neumann exacta** sobre la esfera exterior
  permitiría el 4 % al costo de 1.3 s (dominio chico), en vez de 32 s.
  **La restricción tiene dominio de validez medido**: con una esfera indistinguible del medio hay que
  apagarla, porque ahí la onda incidente **sí** ejerce fuerza neta y el desplazamiento medio no es
  cero (con restricción mete dispersión falsa de 1.6 Pa sobre 1 Pa). Está como opción explícita y
  testeada en los dos sentidos.

  **Faran (1951) — la derivación, guardada (D44).** La serie de la esfera elástica (dos
  potenciales P/S, tres condiciones de borde en `r = a`, sistema 3×3 por orden) quedó a mitad de
  derivación: `u_r = Φ_r + n(n+1)Ψ/r²` y `u_θ = (Φ_θ + Ψ_rθ)/r` **verificados con sympy**, y las
  expresiones de `σ_rr`/`σ_rθ` derivadas pero sin armar. Se borró el módulo a medio escribir a
  propósito: un coeficiente mal derivado no falla ruidosamente. La validación ya está decidida —el
  **teorema óptico** (`σ_ext = σ_sca`), que no necesita tabla de referencia, más los residuos de las
  tres condiciones y el límite de contraste bajo—. Lo que NO sirve como validación es comparar
  contra `sphere.py` con `μ → 0`: ese límite es físicamente correcto pero numéricamente singular
  (`k_β → ∞`).
- **Tensor elástico de la cápsula ótica — no hay dato citable.** Peterson y Dechow (2003) dan
  propiedades direccionales de cortical externa del temporal, no de la cápsula. Se sigue con isótropo
  `VERIFY` y los parámetros expuestos en la UI.
- **CT→velocidad — hay calibración.** Webb et al. (2018): `c = m·HU + b`, ejemplo `m = 0.75
  m·s⁻¹·HU⁻¹`, `b = 1.32 km/s`, `R² = 0.46` (depende de energía y reconstrucción). Entra como valores
  por defecto citados y editables en la fase de geometría de CT.
- **α del hueso a 100 kHz — sigue siendo extrapolación del material.** Li et al. (2025) tienen
  *pérdida de inserción* a 100–290 kHz, no α intrínseca: sirven para contrastar el modelo completo,
  no para fijar el coeficiente.
- **Q de la cavidad — no hay dato citable.** Se mide con la definición operativa de D40
  (`Q = f_pico/Δf` con `|p|²` a la mitad) y la UI permite cargar un Q de referencia propio.
- Referencia externa Elmer/openCFS, y el resto de la Fase 8 (ya completa en U1–U6).

## Física y decisiones de modelado (del plan de Diego)

1. **Los fluidos mandan primero.** Perilinfa/endolinfa se comportan como agua: acústica escalar,
   c ≈ 1500 m/s, ρ ≈ 1000 kg/m³.
2. **Tejido blando**: también fluido acústico, con atenuación.
3. **Hueso temporal (cápsula ótica)**: sólido elástico, denso y rígido, transporta ondas de corte.
   c longitudinal ≈ 3000–4000 m/s, ρ ≈ 1900–2000 kg/m³ (cifras de referencia del plan, a verificar
   contra literatura para el caso de uso). Tratarlo como fluido ignora la conversión de modo
   longitudinal ↔ corte en la interfaz tejido-hueso, que pesa mucho en incidencia oblicua.
   → El modelo realista es **acoplamiento fluido-sólido**: presión en fluidos, desplazamiento en
   hueso, continuidad de velocidad normal y de tracción
   (fluido: ∂p/∂n = ρ_f ω² u_n; sólido: σ·n = −p n).
4. **FEM es el método de trabajo** (heterogéneo, elástico, acoplado). **BEM encaja en las regiones
   homogéneas de fluido** y en dominios exteriores no acotados (cabeza, canal auditivo), donde evita
   fronteras absorbentes. **BEM elástico es mucho más difícil que BEM acústico: no se intenta
   temprano.**
5. **FEM-BEM híbrido** es el destino natural, pero es investigación: primero cada método por
   separado y contrastados en problemas sólo-acústicos.
6. **VIBRA** es referencia conceptual del acoplamiento vibro-acústico, no dependencia (alcance y
   licencia sin verificar).
7. Frecuencia decide todo: λ en agua ≈ 75 mm @ 20 kHz, 15 mm @ 100 kHz, 1.5 mm @ 1 MHz.
   Cóclea ≈ 9 mm de ancho, scalae ≈ 1 mm.
8. Geometría: segmentación CT/micro-CT (3D Slicer) → malla de superficie/volumen (Gmsh).
   Hounsfield → densidad y velocidad en hueso. Estructuras finas (membrana basilar, de Reissner,
   ventanas oval y redonda) piden cáscara delgada / condiciones efectivas / simplificación.
9. Atenuación en hueso alta y dependiente de frecuencia → **módulo complejo desde el principio**.
10. Dominio del tiempo para pulsos: referencia k-Wave (verificar soporte elástico actual en Python).

## Roadmap por fases, con criterio de aceptación

- **Fase 0 — Setup + 1D.** Helmholtz 1D a mano (FEM P1/P2 propio) y solución analítica de referencia.
  *Aceptación*: campo 1D FEM vs analítico, error relativo L2 < 1 % con ≥ 200 elementos; y R/T de una
  interfaz única contra la fórmula de impedancias.
- **Fase 1 — FEM acústico 2D.** Medios estratificados fluido/tejido, atenuación, frontera absorbente
  (ABC de primer orden / DtN). *Aceptación*: R y T de una capa contra el analítico dentro de **2 %**
  en amplitud; reflexión espuria de la ABC por debajo de **−40 dB** en incidencia normal.
- **Fase 2 — Hueso elástico 2D acoplado a fluido.** *Aceptación*: transmisión por placa
  fluido-sólido-fluido vs analítico, y ángulos de Snell con conversión de modo.
  **CERRADA.** Medido: incidencia normal, |T| FEM vs analítico con *0.079 %* de error en la tabla
  (acero 2 mm 0.29345 vs 0.29326; acero 5 mm 0.12678 vs 0.12670; hueso 2 mm 0.79639 vs 0.79587;
  hueso 5 mm 0.50903 vs 0.50863); balance |R|²+|T|² = 0.999017 (acero) y 1.000075 (hueso).
  Snell a 20° en hueso: θ_p 52.53° vs 52.94° (−0.42°), θ_s 25.75° vs 25.67° (+0.08°).
  Resonancia de espesor: acero 1.4752 MHz vs 1.4750 (+0.012 %), hueso 0.8754 vs 0.8750 (+0.041 %).
  A 20° en acero *no hay onda P transmitida* (sinθ_p = (5900/1500)·sin20° = 1.345 > 1, crítico a
  14.73°): queda como test de reflexión total, que es el caso que la Fase 5 necesita para incidencia
  oblicua sobre hueso.
  **Prerrequisito de la Fase 5: satisfecho** (ya hay elástico acoplado; la aproximación "hueso como
  fluido" queda sólo como caso de comparación, no como camino).
  *Formulación* (plano deformación, isotrópico):
  ``∇·σ + ρ ω² u = 0`` con ``σ = λ tr(ε) I + 2μ ε``;
  acoplamiento en la interfaz Γ con normal ``n`` saliente del sólido:
  ``∂p/∂n = ρ_f ω² (u·n)`` y ``σ·n = −p n``.
  Bloque lineal (convención **implementada**, con ``C_ij = ∫_Γ (n·φ_i) ψ_j ds`` y ``n`` saliente del
  sólido):
  ``[K_s − ρ_s ω² M_s, +C; +ω² Cᵀ, K_f − k² M_f/ρ_f] [u; p] = [f_s; f_f]``
  (equivalente al bloque con ``−C`` y ``−ω²Cᵀ`` que aparecía en el planteo: cambiar el signo global de
  ``C`` deja ``p`` idéntico y ``u`` con un desfase π, así que ``|T|``, ``R``, ángulos y resonancias no
  cambian. La pareja física es la de arriba y el test de estado de interfaz la fija.)
  (los signos salen de que la tracción del fluido sobre el sólido es ``−p n`` y de que
  ``v_n = −iω u_n`` con ``e^{−iωt}``).
  *Criterios*: (a) incidencia normal agua–placa de acero–agua: |T| vs la fórmula analítica de tres
  capas dentro de 2 %; (b) incidencia oblicua: ángulos longitudinal y de corte transmitidos medidos
  del campo FEM contra ``sinθ_p/c_p = sinθ_s/c_s = sinθ_i/c_f`` dentro de 1°; (c) resonancia de
  espesor de la placa dentro de 2 % de la frecuencia analítica. Materiales: agua (plan) y acero sólo
  como *benchmark paramétrico* con c_p, c_s, ρ explícitos en el test (no como material de librería,
  D5).
- **Fase 3 — BEM acústico 2D.** Dispersión exterior (implementación propia; bempp-cl como opción).
  *Aceptación*: cilindro rígido vs serie de Bessel, error < 2 % en el campo dispersado.
- **Fase 4 — Contraste FEM vs BEM** en problemas sólo-acústicos. *Aceptación*: discrepancia < 2 % en
  campo y en coeficientes, con mallas comparables; se reporta error y tiempo de cada uno.
  **Cerrada** (cilindro rígido, agua, 100 kHz, R=3a, ppw=20, ABC-2, 64 puntos en r=1.5a):
  |FEM−BEM| = *0.97 %* (ka=1) y *0.59 %* (ka=2); |FEM−serie| = 0.94 % / 0.60 %;
  |BEM−serie| = 0.05 % / 0.14 %. Error por cara: 0.58 % iluminada / 0.75 % sombra (ka=1).
  Tiempos: FEM 12–41 ms, BEM 203–242 ms (mismo problema, 854–3 146 elementos vs 64 paneles).
  Costo de la fase: encontró el error de formulación D12, que ninguna prueba unitaria veía.
- **Fase 4B — Fuentes, arreglos y focalización** (objetivo declarado por Diego: comparar
  transmisores comerciales contra un arreglo sintetizado).
  *Modelo de fuente*: apertura plana sobre baffle rígido infinito con velocidad normal prescrita
  ``v_n``. En 2D la integral de Rayleigh es
  ``p(x) = 2·i·ω·ρ · ∫_apertura v_n(y) · G(|x−y|) dy`` con ``G = (i/4)H₀⁽¹⁾(k r)``; el factor 2 es
  la imagen del baffle rígido. Un transductor = una apertura; un arreglo = N aperturas con amplitud
  y fase (delay) por elemento. Esto es exacto dentro del modelo, así que no hace falta malla ni
  datasheet para validarlo.
  *Métricas*: posición del pico, FWHM lateral en el foco, nivel de lóbulo lateral (dB), ganancia de
  intensidad respecto de un elemento, correlación de campo complejo en la región blanco.
  *Aceptación*: (a) ceros del pistón en ``sin θ = m·λ/D`` dentro de 0.5°; (b) lóbulo lateral del
  arreglo uniforme dentro de 2 dB de la fórmula de array factor
  ``|sin(Nψ/2)/(N sin(ψ/2))|``; (c) y (d) **equivalencia arreglo vs apertura continua** — *la*
  pregunta del usuario. Números de referencia (D=60 mm = 4λ, F=40 mm, f/#=0.67, 100 kHz, agua),
  medidos en el plano del máximo axial real, no en el geométrico:

  | pitch | N | FWHM | desvío vs continua | correlación 2D |
  |---|---|---|---|---|
  | λ   | 4  | 10.96 mm | +27.0 % | 0.9564 |
  | λ/2 | 8  |  9.16 mm |  +6.1 % | 0.9977 |
  | λ/4 | 16 |  8.75 mm |  +1.4 % | 0.9999 |
  | λ/8 | 32 |  8.66 mm |  +0.3 % | 1.0000 |

  Continua de referencia: FWHM *8.66 mm*, con el máximo axial en *29.6 mm* (no en F=40 mm: el
  corrimiento focal no paraxial de un foco de semiángulo 37° es −10.5 mm, confirmado por trapecio
  directo, por la cuadratura del paquete y por la misma apertura a 500 kHz, que enfoca en 39.4 mm).

  *Criterio corregido*: "equivalente" depende de la métrica y hay que decir cuál. Por *FWHM en el
  foco*, λ/4 ya da +1.4 % (cumple el 2 %) y λ/2 +6.1 % (cumple el 7 %); λ/8 lo baja a +0.3 %. Y el
  pitch solo no decide: con 256 elementos y D = 50–133 λ, λ/1.9 da +1.8 % (D25).
  Ganancia en F: 15.95 dB contra
  ``20log10(N)=18.06`` porque los elementos de borde iluminan a 37° y su propia directividad los
  atenúa (|p_j| cae a 0.54 en el borde): *la ganancia de N sólo vale en régimen paraxial*.
  (e) campo lejano decae como ``1/√r`` dentro de 3 %.


  *Prerrequisito del barrido en frecuencia*: pasar el BEM a Burton-Miller o CHIEF (D11.3), porque la
  formulación directa tiene resonancias espurias en los ceros de J_n(ka) y un barrido las cruza.
- **Fase 5 — Cóclea simplificada** (axisimétrica o corte 2D: cápsula ósea + espiral con fluido).
  *Aceptación*: contraste contra Elmer u openCFS como referencia externa.
  **5a (corte transversal: cápsula elástica + fluido interior) ENTREGADA, con dos criterios no
  alcanzables y bien diagnosticados**: `cochlea.py` construye la malla polar con anillos forzados en
  a/2, a, a+t y R_out, y resuelve sobre `coupled2d`. Medido (100 kHz, t=2 mm, hueso 3500/1900/1950):
  convergencia de malla +0.6 % (ppw 24→32) ✓; sin absorción ni flujo neto espurio (1e-15) ✓; cáscara
  casi transparente vs placa *−1.52 %* / *−0.52 %* ✓; campo dispersado vs serie de cilindro rígido
  *12.0 % L2* y potencia dispersada *0.8 %*.
  **No alcanzados (xfail con razón)**: (a) límite plano contra la placa da *+19 %* a a/t=40 y (b) la
  tendencia con curvatura no es monótona ([15.9, 28.5, 0.3] %). Causa: la cápsula es una **cavidad
  cerrada sin pérdidas** — recibe energía por toda la circunferencia y el campo interior es
  estacionario (oscila ±30 % con ±5 % de frecuencia), así que no hay un `T` que converja al de una
  placa. No es un error del código: es topología.
  **5b** = ventanas oval y redonda (la redonda como salida de energía) + membrana basilar separando
  las dos scalae. Eso no es sólo fisiología: es lo que le da al problema un camino de salida y hace
  medible la transmisión. Antes de 5b hace falta mejorar la frontera: con ABC-1 en contorno curvo el
  error de sombra es *17–21 %* y no baja al agrandar el dominio (pad 2.5/4/6 λ → 21/24/27 % L2).
  *Actualizado*: BGT-2 ya está (D21) y elimina el sesgo de sombra (16.9 → 9.4 %), pero el piso quedó en
  **8–9 % por la discretización P1 del acoplamiento** — comprobado con una cáscara de impedancia casada
  sin dispersor, que da el mismo 9.00 %. Bajar de ahí requiere P2 en la interfaz (escalón 3 de D10) y
  **hoy no es prioritario**: la incertidumbre de material (2.4–3.3×, D20) es mucho mayor.
  *Referencia externa (Elmer/openCFS): pendiente, no están instalados; los límites de arriba son el
  sustituto hasta entonces.*
  **5b-2 ENTREGADA (D29)**: la membrana basilar es un segmento con impedancia
  `Z_m = R + iωm + k/(iω)` sobre un diámetro del fluido interior. Rompe la resonancia de cavidad:
  el `Q` medido baja de **287** a **< 5** (factor > 57), y `t_interior` deja de tener el lóbulo que
  lo hacía no monótono (2 cambios de signo → 1, 12 de 13 tramos decrecientes y el último +0.97 %).
  El `n = 0` acompaña: 0.3866 → 0.0761 sin membrana (5.08×) contra 0.1813 → 0.1169 con membrana
  (1.55×). Propiedades de la membrana: **todas VERIFY** (no hay conjunto citable en nuestra banda);
  los criterios se validan por límites (Z→∞ rígido, Z→0 presión nula) y balance de energía, no por
  esos números.
  *Re-barrido de sensibilidad con la membrana puesta* (D30, 85 puntos): el espesor baja de 72 % a
  **53 %**, la atenuación sube de 5.5 % a **9.6 %**, `c_s` sale del piso numérico (2.2 % → **8.2 %**)
  y `Z_r` se cae a **0.3 %**. `c_p` no se mueve (16.1 %).
- **Fase 6 — 3D con geometría de CT** y mapas heterogéneos. *Aceptación*: resultados de literatura
  para skull / hueso temporal.
  **6a ENTREGADA (D32)**: núcleo acústico 3D — `meshing3d.py`, `fem/helmholtz3d.py` (P1, ABC-1 y
  BGT-2) y `aperture3d.py`. Onda plana al **3.53 % a ppw 20** con orden medido **1.88 / 1.94 / 1.97**;
  la ABC baja el error de **100 %** (pared rígida) a **9.57 %**; la fuente puntual reproduce el Green
  3D con amplitud **0.9903**; el pistón coincide con Rayleigh-Sommerfeld al **0.0000 %** y su último
  máximo axial cae en `a²/λ − λ/4`. El techo es el costo: 2 λ × 1 λ × 1 λ a ppw 20 son 18 081 nodos y
  **5.3 s** por solve; un cubo de 3 λ a ppw 12 no termina en 10 min. Falta para "cóclea 3D": R/T
  estratificado en 3D, dispersión contra obstáculo, elasticidad 3D y geometría de CT.
- **Fase 7 — Opcional**: híbrido FEM-BEM, pulsos en tiempo, GUI (frontend C# si se quiere).


## Fase 8 — UI de trabajo: del explorador de campo al banco de simulación

### Qué hay hoy, y dónde está el hueco

La UI actual (D27/D31) responde *"¿cómo se ve el campo de esta apertura?"*. No responde *"corré este
problema FEM/BEM y mostrame qué salió"*. Concretamente:

| lo que hay | lo que falta para trabajar |
|---|---|
| motor de campo tipo Rayleigh sobre aperturas (2D) | correr `fem.helmholtz2d`, `fem.elastic2d`, `fem.coupled2d`, `bem.helmholtz2d`, `fem.helmholtz3d` |
| una sola región implícita (baffle o espacio libre) | escena: regiones con material, interfaces, dispersores |
| sliders de apertura y foco | geometría y malla editables, con vista previa de la malla |
| campo en un plano fijo | campo de cualquier solver, cortes, y 3D con rebanadas |
| un solve por click, síncrono | trabajos con progreso y cancelación (un 3D tarda minutos) |
| panel de cóclea con 3 parámetros | la escena de cóclea completa (cápsula + scalae + membrana + ventanas) |
| sin barridos | barrido de frecuencia/geometría en cola, con CSV y comparación entre corridas |

El punto delicado no es pintar más cosas: es que **la UI no puede inventar física**. Cada resultado
tiene que decir contra qué número ya validado se comprobó.

### Principio rector

La UI es un **cliente de los solvers que ya existen**, no un segundo simulador. Todo número que
muestre sale de un módulo con test. Si la UI necesita algo que la librería no tiene, primero se
agrega a la librería con su criterio, después se enchufa.

### Las fases

- **U1 — Escena 2D + FEM acústico. ENTREGADA.** `dashboard/scene.py` (validación, malla, solve,
  métricas) + `meshing.grid_tri_mesh` (reticulado conforma por construcción) + pestaña "escena FEM"
  en el dashboard. Medido: agua/hueso/agua a ppw 12 — *R 0.9018 vs 0.9037* (abs 0.0019),
  *T 0.4323 vs 0.4281* (abs 0.0041), energía *1.0000*, solve *11 ms* con 1 152 elementos;
  agua/tejido/agua — R abs 0.0060, T abs 0.0009. 18 tests en `test_scene.py` + 5 de HTTP.
  Editor de regiones rectangulares (material por región),
  vista previa de malla con el conteo y el costo estimado, selector de formulación
  (total/dispersado), ABC-1/BGT-2, frecuencia, y mapa de |p| + cortes.
  *Aceptación*: desde la UI, la capa agua/tejido/agua reproduce R y T de la Fase 1 dentro de **2 %**
  (la librería da 0.511 % — la UI no debe empeorarlo más de lo que explica su propia malla), y la
  escena guarda/carga en JSON con la malla reproducible bit a bit.
  *Bloqueante técnico*: el servidor es síncrono hoy; un solve 2D de 5 k elementos tarda ~1 s, así que
  U1 puede quedar síncrona. Recién 3D obliga a trabajos asíncronos.
- **U2 — Dispersores y BEM. ENTREGADA (cilindro rígido).** `dashboard/scatterer.py` + pestaña
  "dispersor BEM": el mismo cilindro rígido en agua resuelto por BEM (exterior), FEM (corona en
  campo dispersado) y la serie de Bessel, con el mapa del campo dispersado, las tres curvas `|p|`
  contra θ en `r = 1.5a` y las tres discrepancias. Medido: ka=1 → *FEM−BEM 0.97 %*, *FEM−serie
  0.94 %*, *BEM−serie 0.05 %*; ka=2 → 0.58/0.60/0.14 %; ka=3 → 0.84/0.72/0.27 %. Las tres
  ecuaciones (standard/Burton-Miller/CHIEF) dan el mismo campo dentro del criterio. Costo 2.2 s
  por corrida (de los cuales 1.9 s es el mapa BEM).
  Falta de esta fase: elipse/polígono (contornos BEM no circulares) y dispersor **penetrable**
  (el mismo frente que D32 dejó abierto en 3D).
  Editor de regiones rectangulares (material por región),
  (`standard`/`burton_miller`/`chief`), y el **mismo problema por los dos métodos lado a lado**.
  *Aceptación*: el cilindro rígido contra la serie de Bessel da < 2 % mostrado en pantalla, y el panel
  de contraste FEM/BEM muestra la discrepancia y el tiempo de cada uno (referencia: 0.97 %/0.59 %,
  FEM 12–41 ms vs BEM 203–242 ms).
  Es la fase que más rinde: es el uso que Diego pidió desde el principio —comparar métodos y ver
  dónde no coinciden.
- **U3 — Cóclea de trabajo. ENTREGADA.** `dashboard/cochlea_field.py` + pestaña "cóclea": el corte
  completo de la Fase 5 (cápsula elástica, perilinfa, membrana basilar, ventana redonda) resuelto
  por `cochlea.solve_cross_section`, con el **mapa del campo interior** interpolado sobre el disco,
  corte por el centro, y las métricas de D29/D30 en vivo. Medido: `t_interior = 0.53187` (D30:
  0.53187, desvío *0.00 %*), armónico n=0 0.17913, membrana absorbe *72.9 %*, balance de energía
  *0.09 %*. Sin membrana: `t_interior` 0.97218 y absorción exactamente 0 (el interruptor mueve la
  física, no un parámetro cosmético). Costo *4.3 s* (malla 0.47 s + solve 4.3 s).
  **Q no se muestra**: medirlo pide un barrido de frecuencia (~9 solves) y eso es U5 — la UI lo dice
  en vez de inventar un número.
  La escena completa de 5a+5b: cápsula con espesor, dos scalae, membrana
  basilar con sus tres parámetros, ventanas. Vista del campo interior, `t_interior`, balance de
  energía y `Q` en vivo.
  *Aceptación*: reproduce `t_interior = 0.53187` dentro de 1 %, energía dentro de 0.46 % y `Q < 5`
  (los números de D29/D30), y los tres sliders ya expuestos siguen moviendo el resultado con los
  pesos medidos (53 %/8 %/10 %).
- **U4 — 3D. ENTREGADA (dos modos, con rebanadas).** `dashboard/volume3d.py` + pestaña "3D" con el
  plano *x-y* a mitad de *z* y el *x-z* a mitad de *y*. Dos modos: **puntual** (monopolo en agua, ABC
  en las seis caras, se compara con el Green 3D) y **plano** (capas, ABC sólo en *x*, R/T contra la
  matriz de transferencia). Medido: puntual → error Green *7.0 %* con amplitud *0.969* en 9 261 nodos
  y *1.0 s*; plano → *R 0.9032 vs 0.9037* (abs 0.0005), *T 0.4343 vs 0.4281* (abs 0.0061), energía
  *1.0044*. Tope de 90 000 nodos con mensaje, y el techo conocido de D32 declarado antes de resolver.
  Escena de cajas (o malla externa), solver P1 con ABC-1/BGT-2, y visualización por
  **rebanadas** (no volumen completo: con 18 k nodos el campo es un sólido y Plotly no lo pinta bien).
  *Aceptación*: reproduce D32 desde la UI —onda plana 3.53 % a ppw 20, pistón 0.0000 % contra
  Rayleigh-Sommerfeld, R/T estratificado 2.91 %/2.11 %— con el conteo de nodos y el tiempo a la vista.
  *Bloqueante*: acá sí hace falta cola de trabajos (`POST` devuelve job id, el front consulta estado y
  progreso, hay cancelar) y un tope duro de elementos con mensaje claro en vez de colgar la máquina.
- **U5 — Barridos con cola, progreso y cancelación. ENTREGADA (sensibilidad).** `dashboard/jobs.py`
  (cola FIFO, un worker, cancelación cooperativa, retención acotada) + `dashboard/sweeps.py` (barrido
  OAT de la cóclea) + pestaña "barridos" con barra de progreso, cancelar y ranking en vivo.
  *Aceptación cumplida*: el barrido de 31 puntos reproduce el ranking de D30 **exacto** — espesor
  *53.11 %*, c_p *16.06 %*, atenuación *9.60 %*, c_s *8.22 %*, ρ *1.03 %*, Z_r *0.28 %* (peor desvío
  *1.8 %*), en *89 s* con `t_interior` base 0.53187. Los α de las tres leyes se piden a
  `attenuation.py` (3.777 / 6.243 / 15.312 Np/m a 100 kHz), no se escriben a mano.
  **Barrido de frecuencia de la escena 2D — ENTREGADA (D43).** Barre la escena de la pestaña "escena
  FEM" en frecuencia desde la cola, mide R y T por punto y marca las **resonancias de media onda
  predichas** `f_n = n·c_losa/(2t)`. Medido con una losa de hueso de 10 mm: picos en *176.7 y
  353.8 kHz* contra *175 y 350* predichos (≤1 %), con *T = 0.999* y *R = 0.031* coincidiendo en la
  resonancia. El sello marca "fuera" a ppw 12 porque el error de T llega a *5.90 %* arriba de 350 kHz
  (a ppw 20 baja a *2.08 %*): es la malla, no el barrido, y el barrido lo dice en vez de promediar.
  Falta de esta fase: el barrido de las 15 configuraciones de D25 y la exportación a CSV desde la UI.
- **U7 — Esfera elástica. ENTREGADA (D47).** Pestaña nueva: material (acero, hueso, aluminio,
  tungsteno, tejido casi fluido), `ka`, frecuencia, dos mallas y el interruptor de la restricción de
  modos rígidos con su tooltip. El servidor resuelve con `coupled3d` y contrasta contra **Faran** en
  la misma corrida; la UI pinta el campo dispersado en el plano x-z, `|p_s|` contra θ y la
  **diferencia de fase**. Criterio nuevo `volumen_esfera_elastica` (6 %).
  Medido a `ka = 1`, dominio 4a, `h = a/4.5` (21 398 nodos, ~40 s): **acero 4.29 %**, aluminio
  *3.99 %*, tungsteno *4.71 %*, **hueso 13.16 %**.
  *El hueso es un error de MODELO, no de malla*: su contraste de impedancia es 4.6 contra 30.7 del
  acero, así que el cuerpo se mueve parcialmente con el fluido y la condición de traslación media
  nula deja de ser exacta. Apagar la restricción lo empeora a 54.6 %, o sea que sigue siendo
  imprescindible. El sello lo marca `fuera` con el número y la fuente: **no se subió el umbral para
  que la tabla se viera verde**.
  *Con tejido casi fluido (contraste 1.08) la restricción se rompe del todo*: el campo salía **53×**
  el exacto. Ahí la restricción **se desactiva sola** (umbral de contraste 2.5, medido) y el sello
  pasa a **`sin_referencia`**, porque con un campo dispersado del orden del 2 % de la incidente el
  error *relativo* no indica nada: se muestra el **error absoluto** contra la incidente
  (7.2e-3 Pa sobre 1 Pa).
  *Formulación dispersada IMPLEMENTADA (D49), y el resultado obliga a matizar.* `u_s = u − u_i` con
  el **operador de contraste** (no el completo: restar el completo deja un residuo O(h²) que empuja un
  campo espurio del 18 %). El límite transparente ahora **valida exacto**: `|p_s|` converge a cero
  como `c_s²` (6.7e-2 → 2.6e-3 → 5.3e-4 → 1.2e-4 para `c_s` 100/30/10/3) sin restringir nada.
  **Pero contra Faran es peor**: acero **19.49 %** y aluminio **44.14 %**, contra 4.29 %/3.99 % de la
  total+restricción. La razón: para un metal `u_s ≈ −u_i` y `(K−ω²M)u_s` con la rigidez del acero es
  una resta de términos grandes; el error de interpolación de la incidente se amplifica por `K`.
  **Consecuencia honesta**: el 4.29 % está en parte **cancelado** con el error de discretización, no
  ganado. El default sigue en la total (`fijar_rigidos=True`) porque es la sellada, y el criterio queda
  calificado a metales.
  *Barrido de apertura en la UI (U8, D53)*: FWHM lateral contra `pitch/λ` con la apertura continua
  como referencia. Medido a 60 mm, foco 90 mm, 250 kHz: 1.25 λ → **10.8 %** de desvío (y grataje a
  53°), 0.45 λ → 2.3 %, **0.28 λ → 1.7 %**, 0.16 λ → 1.5 %. El pitch más grueso que cumple el 2 % es
  **1.1× el `λ/4` de D25/D26**: coincide con el criterio y le pone el número al lado.
  *Barrido en `ka` en la UI (U7b)*: la vista resuelve la esfera en varios `ka` y grafica el error de
  cada formulación contra Faran, con la tolerancia del criterio como línea. Primer resultado con el
  modo rápido (dominio 2.5a, 5 puntos): la **dispersada gana en `ka = 0.5`** (15.2 % contra 17.7 %) y
  **pierde entre `ka = 1.0` y `1.5`** (19.6/20.1 % contra 10.8/7.9 %). O sea: la formulación correcta
  mejora justo en el límite de Rayleigh, que es donde la traslación `3/(2g+1)` pesa, y empeora donde
  el campo dispersado tiene estructura. Es la confirmación numérica de D50 en el rango, no sólo en el
  punto.
  *Traslación contrastada (D50), y el resultado es limpio*: a `k₀a = 0.1` las dos formulaciones
  correctas dan `|u_cm|/|u_i|` al **1-2 %** de `3/(2g+1)` (acero 0.1772 vs 0.1807; hueso 0.6064 vs
  0.6122). El acoplamiento fluido-sólido y la masa añadida están bien ⇒ **el 19-44 % a `ka = 1` es
  resolución del resto del campo, no de la física**. Y la restricción no "aproxima" la traslación: la
  pone en **0.0000** donde la ley pide 0.18, y en el caso casi transparente la **amplifica 10×**
  (10.09 contra 0.9677).
  *Siguiente paso*: convergencia en `h` de la dispersada en acero a `ka = 1` para confirmar que el
  19.49 % baja al refinar. Si no baja, el sospechoso es la estructura angular (elementos P1 rectos
  sobre una esfera facetada), no la formulación.
  **Y no es una mejora opcional (D48)**: Sim-Datos confirmó que la premisa de la restricción es
  **falsa** — una onda plana sí ejerce fuerza armónica de primer orden sobre el cuerpo, y Faran
  **incluye** ese *rebound* traslacional (Feuillade y Clay 1999, JASA 106:553, §IV). La amplitud
  suprimida escala como `3/(2g+1)` con `g = ρ_sólido/ρ_fluido` (Bruus 2012, ecs. 22 y 26): 0.18 en
  acero, 0.61 en hueso, 0.97 en tejido — exactamente el orden del error medido. La restricción es una
  **limitación de la implementación**, no del método, y el arreglo es la formulación dispersada.
  También de D48: no hay umbral citable de contraste para "dispersor débil" (los de Born son
  `|γ_ρ|, |γ_κ| ≪ 1` y `x·|γ| < 1`, sin porcentaje garantizado), y `c_s` de tejido blando **no tiene
  fuente en la banda** (las de elastografía son de cientos de Hz): esos materiales quedan rotulados
  como ballpark, no como dato.
- **U6 — Trazabilidad. ENTREGADA.** `dashboard/validacion.py` con un registro de *10 criterios*, cada
  uno con su descripción, referencia, tolerancia, unidad y **fuente** (fase o ADR). Cada endpoint
  devuelve `sello` ya evaluado —`validado` / `fuera` / `sin_referencia`— y la UI sólo lo pinta. El
  resumen **no promedia**: si un criterio falla, el conjunto falla. `sin_referencia` es un estado
  explícito, no un aprobado por omisión: la cóclea con espesor 3 mm lo muestra (el caso base de D30
  no la cubre). Test que fija el contrato: ningún JS del cliente puede contener `0.53187` ni
  `CRITERIO_*` — si aparecen, la validación volvió a partirse entre servidor y cliente.
  Cada resultado lleva su "sello": qué criterio, contra qué número de
  referencia y con cuánto margen. Un resultado sin sello se marca como *no verificado*.
  *Aceptación*: la UI rechaza mostrar un número como validado si el test que lo respalda no corrió
  para esa configuración —el sello sale de la corrida, no de una constante escrita a mano en el JS.
  Esta fase es la que evita que dentro de tres meses no sepamos qué era confiable.

- **U8 — Banco de simulación. ENTREGADA.** Lo que convierte el explorador en una herramienta con la
  que se trabaja horas sin perder el hilo. Cinco piezas, cada una con su decisión:
  1. **Vista `criterios`**: el registro de validación del servidor (12 criterios con tolerancia,
     unidad y fuente) cruzado con lo medido en la sesión. Un criterio sin corrida dice `sin medir`,
     no un aprobado por omisión.
  2. **Persistencia por vista** (`localStorage`, TTL 12 h): volver a una pestaña cara (esfera precisa
     ~30 s, 3D, cóclea, barridos) ya no cuesta otro solve; se repinta con su fecha. Restaurar no
     re-guarda, o el TTL nunca vencería.
  3. **Exportación transversal**: un botón que baja PNG de cada gráfica visible y un JSON con sus
     trazas, sin saber nada de cada vista. En la esfera, además, JSON completo, CSV de `|p_s|(θ)` y un
     historial de corridas que se **superponen** para comparar.
  4. **Comparar dos campos** en la vista de campo: referencia congelada y mapa de **diferencia en dB**
     (±20 dB, escala divergente), con el perfil de la referencia punteado en los cortes.
  5. **Barridos nuevos**: apertura (FWHM vs `pitch/λ`, D53) y `ka` de la esfera con "comparar", que
     resuelve las dos formulaciones en una corrida y marca el cruce (D51).
  *Aceptación cumplida*: las ocho vistas resuelven sin errores de consola y sellan; los barridos
  nuevos traen número medido (0.28 λ → 1.7 % de desvío; `ka` 1.17–1.58 validado y fuera de ahí 6.8 %
  y 9.2 %); y el informe de la sesión sale en Markdown con las corridas y sus fuentes.
  *Lo que la UI sigue sin hacer*: no mide nada propia (todo viene del servidor), no elige formulación
  por el usuario y no promedia criterios.

- **U9 — Estructuras propias: STL + acople multicapa. PENDIENTE.** Es el eje que falta para el objetivo
  declarado (ultrasonido atravesando piel → hueso → cóclea). Hoy `coupled2d`/`coupled3d` acoplan
  **fluido ↔ un sólido** con una malla etiquetada; no hay forma de meter una geometría propia ni de
  posicionarla.
  1. **Lector STL propio** (`stl.py`): binario y ASCII, sin dependencias, con control de que la malla
     sea **cerrada** (cada arista en exactamente 2 triángulos) y de que las normales apunten hacia
     afuera (volumen por divergencia > 0). Un STL abierto o invertido tiene que fallar con un mensaje
     que diga cuál de las dos cosas es, no mallar a ciegas.
  2. **STL → malla de volumen** (`meshing_volume.py`) con gmsh: un STL por región, `occ.fragment` para
     las que comparten superficie, etiqueta de región por volumen físico.
  3. **Acople multicapa**: N regiones sólidas contiguas + fluido, reusando la formulación de
     `coupled2d`/`coupled3d` (D41, sin tocar la física).
  4. **Posicionamiento desde el dashboard**: traslación y rotación de la estructura, y las fuentes
     alrededor de ella, con la misma interacción de arrastre que ya tienen las fuentes.
  *Aceptación*: (a) transmisión normal a través de una placa de hueso en agua dentro del **2 %** del
  coeficiente analítico de tres medios; (b) caso **transparente** (placa con las propiedades del agua)
  con reflexión **< 1e-6** relativa; (c) un STL cerrado de una esfera malla y da el mismo campo que la
  esfera analítica ya validada, dentro del error de esa malla; (d) un STL **abierto** falla con mensaje
  explícito; (e) rotar la placa 90° mueve la transmisión como `cos θ` en la ventana de validez;
  (f) la cadena piel + hueso modelada a **500 kHz** cae dentro de la **transmisión craneal medida de
  31 ± 7 % parietal** (Riis, Webb y Kubanek 2022) — es el único punto experimental citable en la banda y
  el que decide si el acople multicapa sirve.
  *Materiales*: `docs/materiales-cadena.md`, con cita por valor y **NC** donde no hay fuente. La mayoría
  son parámetros de modelado publicados, no mediciones de cada tejido, y la UI tiene que decirlo.

- **U10 — Arreglos de fuentes desde archivo. PENDIENTE.** `data/transducers.toml` existe y la UI **no lo
  usa**: cargar un arreglo planeado obliga a tipear 32 filas.
  1. Cargar el catálogo de transductores en la vista de campo (selector de equipo → `n`, `pitch`,
     frecuencia, foco).
  2. Cargar un arreglo propio desde CSV/JSON (posición, ancho, amplitud, fase por elemento), con
     validación de unidades y de que los elementos no se solapen.
  3. Exportar el arreglo editado con el mismo formato, para poder volver a cargarlo.
  *Aceptación*: cargar un arreglo de 16 elementos desde archivo da **el mismo campo** que armarlo a
  mano (diferencia relativa **< 1e-9**, es el mismo solve); un archivo con un elemento de ancho mayor
  que el pitch falla con mensaje; y el catálogo trae los campos `source`/`verified` a la vista, con
  `VERIFY` cuando no lo están.

- **U11 — Que el dashboard guíe. PARCIAL.** El tour de 9 pasos explica la pantalla una vez. Lo que
  falta es guía **por configuración**: cada vez que el usuario cambia algo que tiene consecuencias
  (material, frecuencia, número de fuentes, geometría), decir en el momento qué implica y qué se
  espera —por ejemplo «en hueso λ es 2.33× más corta: con este pitch vas a tener lóbulos de rejilla en
  ±X°», que ya se calcula—. No es un manual: es el sistema diciendo qué está pasando con lo que el
  usuario acaba de tocar.

### Orden y dependencias

U1 → U2 → U3 → U5, con U4 en paralelo cuando la cola de trabajos esté (que la introduce U4 misma).
U6 es transversal pero se hace al final de U3, cuando ya hay suficientes números para sellar.

Razón del orden: U1 da la infraestructura de escena que todo lo demás reusa; U2 entrega el uso
principal declarado (contrastar métodos) con el costo más bajo; U3 y U5 son lo que convierte esto en
herramienta de investigación; U4 es lo más caro y lo más expuesto al techo de cómputo.

### Lo que la UI *no* va a hacer

- No va a correr barridos paramétricos largos en el navegador: eso es `scripts/` + CSV.
- No va a prometer 3D a 1 MHz: el techo P1 está medido y la UI lo va a mostrar antes de arrancar.
- No va a tener editor CAD: cajas, rectángulos, círculos y arcos cubren lo que el modelo necesita hoy.

## Herramientas

- Núcleo FEM: FEniCSx o scikit-fem. BEM: bempp-cl. Validadores externos: openCFS / Elmer.
  Material de aprendizaje: tutoriales de BEM y AcousticBEM. **Descartado**: CalculiX, OpenFOAM.

## Régimen de frecuencia: barrido 40 kHz – 1 MHz, con la focalización como objetivo

Diego fijó el objetivo: **simular transmisores comerciales y comparar si un arreglo de fuentes logra
resultados equivalentes, sobre todo en focalización**. Eso obliga a un barrido, y el barrido cambia el
plan: la focalización no es un problema de "más frecuencia es mejor", es un problema de λ contra el
tamaño del blanco.

Números medidos (agua c=1500, cóclea ≈ 9 mm de ancho, apertura D=10 mm y foco F=20 mm):

- 40 kHz → λ = 37.5 mm; λ/ancho cóclea = 4.17; FWHM de foco ≈ *75 mm*. La cóclea es un cuarto de
  longitud de onda: **no hay focalización posible**, el campo que llega es casi uniforme.
- 100 kHz → λ = 15 mm; FWHM ≈ 30 mm. Sigue siendo más grande que la cóclea.
- 250 kHz → λ = 6 mm; FWHM ≈ 12 mm. Empieza a tener sentido hablar de foco.
- 500 kHz → λ = 3 mm; FWHM ≈ 6 mm.
- 1 MHz → λ = 1.5 mm; FWHM ≈ 3 mm. Foco comparable a la cóclea.

Consecuencia directa: **por debajo de ~250 kHz la comparación "transmisor vs arreglo" no mide
focalización, mide iluminación**. Si el objetivo es focalizar sobre la cóclea, el régimen útil es
250 kHz – 1 MHz, y 40 kHz sirve como referencia de "sin foco" (y para el caso de conducción ósea,
donde el objetivo no es un punto sino una región).

Costo de malla y precisión (medido con `scripts/bench_frequency.py`, ducto, dominio físico de 20 mm):

- 100 kHz (1.33 λ): ppw=10 → 392 elementos, error 7.66 %; ppw=20 → 1 458 elementos, *2.10 %*;
  ppw=40 → 5 832, *0.53 %*. Orden de convergencia medido: 1.86–2.00 (es P1, como corresponde).
- 500 kHz (6.67 λ): ppw=10 → 37.8 %; ppw=20 → 35 912 elementos, 9.83 %.
- 1 MHz (13.3 λ): ppw=10 → 72.7 %; ppw=20 → 142 578 elementos, 19.64 %.

*El cuello de botella no es la memoria: es la dispersión numérica de P1, que se acumula con la
cantidad de longitudes de onda recorridas.* Para bajar de 2 % a 1 MHz sobre 13 λ haría falta
ppw ≈ 60, o sea ~1.3 M elementos 2D: fuera de alcance con este ensamblado. Por eso en 2D el techo
cómodo de P1 es ~100–250 kHz, y subir de ahí es un problema de *orden de elemento*, no de RAM.
En 3D el techo práctico es aún más bajo; a 1 MHz hay que pasar a FEniCSx+MPI, a k-Wave o resolver
el exterior por BEM, que escala con la superficie y no con el volumen.

## Herramientas de inspección

- `scripts/bench_frequency.py` — costo de malla y error por ppw en un ducto.
- `scripts/cross_validate.py` — FEM contra BEM sobre el mismo problema.
- `scripts/sensitivity_cochlea.py` — barrido one-at-a-time del modelo de cóclea (CSV de 87 filas).
- `scripts/focusing_study.py` — las 15 configuraciones de focalización (CSV + figuras, 139 s).
- `scripts/report.py` — `out/report.html`, página de estado autocontenida (sin servidor).
- `python -m ultrasound_sim.dashboard` — *dashboard interactivo local* en `http://127.0.0.1:8765/`:
  sliders de frecuencia, elementos, apertura y foco; mapa del campo con los cortes axial y lateral;
  editor de fuentes (posición, ancho, amplitud y fase por fuente) sobre el mismo motor (D27).



