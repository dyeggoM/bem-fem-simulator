# Propiedades de la cadena piel → hueso → cóclea

Todo lo de acá viene de una consulta a **Sim-Datos** (2026-11-21) con cita por valor. Nada está estimado
por el simulador: donde la fuente no da un número, dice **NC** (no citable) y el código tiene que tratarlo
como supuesto explícito, no rellenarlo.

Convención: `f` en **MHz**; `α` es atenuación de **amplitud**, `p(x) = p(0)·exp(−αx)` cuando `α` está en
Np/m. `α = α₀·fⁿ`.

**Advertencia que hay que arrastrar a toda la UI**: no existe un conjunto experimental completo validado
de 0,1 a 1 MHz para las ocho capas. La mayoría son **parámetros de modelado publicados**, no mediciones
directas de cada tejido. Un simulador que los presente como "propiedades medidas" está mintiendo.

## Capas

- **Piel** (homogeneizada, sin separar epidermis/dermis): ρ = **1116** kg/m³ · c_L = **1537** m/s ·
  c_S = NC · α = **4,6·f** Np/m (α₀ 4,6 Np·m⁻¹·MHz⁻¹, n = 1).
  *Fuente*: Pichardo (2023), BabelBrain, tabla 1, atribuido a Duck. Parámetros de modelado, sin
  validación experimental demostrada en toda la banda.

- **Grasa** (genérica, no subcutánea craneal): ρ = **950** kg/m³ · c_L = **1478** m/s · c_S = NC ·
  α = **0,48·f^1,1** dB/cm.
  *Fuente*: Gu y Jing (2019), tabla I. Parámetros de simulación, no ley experimental validada en
  0,1–1 MHz.

- **Músculo** (genérico, sin orientación de fibras): ρ = **1050** kg/m³ · c_L = **1547** m/s · c_S = NC ·
  α = **1,09·f^1,1** dB/cm.
  *Fuente*: Gu y Jing (2019), tabla I. No asumir equivalencia entre propagación paralela y perpendicular
  a las fibras.

- **Hueso cortical craneal**: ρ = **1896,5** kg/m³ (homogeneizado de modelado) ·
  c_L = **2384 ± 130** m/s a 0,270 MHz y **2471 ± 90** m/s a 0,836 MHz ·
  c_S = **1313–1577** m/s a 0,270 MHz y **1545–1758** m/s a 0,836 MHz (rango entre tres especímenes,
  modelo inverso) · α_L = **33 ± 9** y **240 ± 9** Np/m a esas dos frecuencias.
  Ley de modelado publicada: α_L = **163·f** Np/m (n = 1), que **no** es un ajuste exacto de esos puntos.
  *Fuente*: c y α de Pichardo, Sin y Hynynen (2011); c_S de Pichardo et al. (2017), tabla 6; ρ y ley de
  α de Pichardo (2023).

- **Diploe / hueso trabecular craneal**: ρ = **1738,0** kg/m³ · c_L = **2140 ± 130** m/s a 0,270 MHz y
  **2300 ± 100** m/s a 0,836 MHz · c_S = **1200–1365** y **1252–1574** m/s · α_L = **34 ± 13** y
  **216 ± 16** Np/m. Ley de modelado: α_L = **162·f** Np/m.
  *Fuente*: las mismas que el cortical.

- **Perilinfa**: ρ ≈ **1000** kg/m³ y c_L ≈ **1500** m/s, **aproximación publicada a agua**, no medición
  específica · α = **NC** (no se sustituye por la del agua).
  *Fuente*: Obrist (2019) §1.2; revisión de mecánica coclear (2014).

- **Endolinfa**: ρ ≈ **1000** kg/m³ · c_L = **1500** m/s como parámetro de modelo · α = **NC**.
  *Fuente*: "Modeling of magnetic vestibular stimulation…" (2025), tabla 1.

- **Membrana basilar**: espesor h = **7,5 µm** en la base → **2,5 µm** en el ápice · E = **50 MPa** en la
  base → **3 MPa** en el ápice. ρ, c_L, c_S y `α(f)` = **NC** en esta banda.
  *Fuente*: Wang, Keefe y Gan (2016) §III.D (parámetros del modelo adulto de Zhang y Gan 2011). Son
  valores **asignados a un modelo FEM**, no módulos dinámicos medidos a 100 kHz–1 MHz.

## Benchmarks medidos (esto es contra lo que hay que validar)

- **Atenuación craneal efectiva ≈ 20 dB/cm a 836 kHz** en cráneo humano. Es atenuación **longitudinal
  efectiva**, de transmisión ex vivo con ajuste inverso: no es absorción pura ni pérdida total universal
  de una trayectoria.
  *Fuente*: Pichardo, Sin y Hynynen (2011); citado también en la discusión de Pichardo et al. (2017).

- **Transmisión de presión a 500 kHz: 31 ± 7 % parietal · 13 ± 2 % frontal · 16 ± 4 % occipital.**
  Son transmisiones **regionales**, no coeficientes intrínsecos del material.
  *Fuente*: Riis, Webb y Kubanek (2022), *Acoustic properties across the human skull*.

## Precauciones que el código tiene que respetar

1. **Z no viene tabulada** en esas fuentes. Si el modelo necesita impedancia, se calcula `Z_L = ρ·c_L` y se
   reporta como derivada, no como valor medido.
2. `α` de cizalla ósea: 57 Np/m a 270 kHz y 373 Np/m a 836 kHz, como parámetros **globales ajustados**, no
   separados por cortical y diploe (Pichardo et al. 2017).
3. **No sumar** una pérdida total por transmisión a un modelo que ya reproduce reflexiones y dispersión:
   se duplicarían las pérdidas.
4. Los datos de **calvaria no validan** la porción petrosa del temporal que rodea la cóclea.
5. Tejidos blandos sin `c_S` **no** son sólidos sin cizalla: es que no hay valor citable.
