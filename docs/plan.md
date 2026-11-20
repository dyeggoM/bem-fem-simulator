# Plan — simulador de ultrasonido en cóclea (BEM + FEM)

Estado: activo. Este documento manda. Lo que cambie se anota también en `docs/decisions.md`.

## Objetivo

Simulador de propagación de ultrasonido en tejido blando, hueso y fluidos cocleares, con **dos
métodos numéricos independientes** (FEM y BEM) para poder contrastarlos: la discrepancia entre ambos
es la señal de un error o de una suposición mal puesta.

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



