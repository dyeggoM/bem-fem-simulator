# Referencias externas — código y manuales que podemos consultar

Lista aportada por Diego (2026-11-20). No es una lista de dependencias: son **fuentes de consulta**
para cuando algo se traba — una formulación, una condición de borde, un benchmark, una constante.

Regla de uso: de acá se toma **cómo está escrito** algo (ecuación, convención, caso de verificación),
no se copia código sin revisarlo, y todo lo que entre al repo tiene que pasar por su propia
validación y quedar en `decisions.md` con la fuente. Un número que sale de uno de estos repos sin
cita citable se marca `VERIFY`, igual que cualquier otro.

## BEM acústico

- <https://github.com/lzhw1991/AcousticBEM> — BEM acústico en C++/MATLAB. Útil para contrastar la
  **ecuación integral de frontera**, el tratamiento de la singularidad del kernel y las variantes
  Burton–Miller/CHIEF (que ya tenemos: `bem/`, D21).
- <http://dehesa.freeshell.org/BEMLIB/> — BEMLIB (Pozrikidis): librería de BEM con los kernels
  escritos explícitamente. Bueno para verificar los coeficientes `c_ij` de la fórmula de
  representación y las integrales singulares.
- <https://www.boundary-element-method.com/acoustics/manual/abem.htm> — manual de ABEM. Referencia
  de formulación y de **casos de verificación** (esfera pulsante, esfera rígida, radiación de un
  pistón), que es justo el tipo de caso que usamos como puerta.
- <https://mhamilt.github.io/ca-knowledgebase.github.io/tutorials/bem/bem-tutorial-1.html> — tutorial
  de BEM acústico 2D/3D con las ecuaciones paso a paso. Útil para chequear convenciones de signo
  (`e^{-iωt}` vs `e^{+iωt}`) cuando dos implementaciones no coinciden.

## FEM y multipropósito

- <https://github.com/FEniCS> · <https://github.com/orgs/FEniCS/repositories> — FEniCS/FEniCSx. Es
  la salida que el plan ya contempla para pasar de MHz a 3D grande: formulación débil de alto nivel,
  mallado con Gmsh y paralelismo MPI. Ver D16.
- <https://github.com/FreeFem> — FreeFEM. Mismo uso; su lenguaje de EDP es cómodo para prototipar una
  formulación débil nueva antes de escribirla en P1 a mano.
- <https://github.com/ElmerCSC/elmerfem> — Elmer: multiphysics con módulo acústico y **acoplamiento
  fluido–estructura** ya resuelto. Es la referencia más directa para la fase 6b: cómo plantean la
  interfaz (continuidad de desplazamiento normal y de tracción) y qué elementos usan.
- <https://github.com/calculix> — CalculiX: solver de sólidos (elasticidad). Útil para contrastar la
  parte **elástica** sola, sin acoplar: modos propios de una esfera o de un tubo.
- <https://github.com/dealii/dealii> — deal.II: FEM adaptativo de alto orden. Referencia para
  estimación de error y mallas curvas (la superficie facetada es hoy una fuente de error, D46).
- <https://github.com/OpenFOAM> — OpenFOAM: volúmenes finitos, no FEM, pero su tratamiento de
  condiciones de contorno absorbentes y de interfaz es un contraste útil.
- <https://github.com/idaholab/moose> — MOOSE: framework multiphysics con módulo de acústica y de
  elasticidad. Bueno para ver cómo acoplan y qué tolerancias reportan.

## Acústica estructural / vibroacústica

- <https://github.com/MOPT-UFSC/VIBRA> — VIBRA (UFSC): análisis vibroacústico por FEM/BEM con
  **cáscaras y sólidos acoplados a fluido**. Es el más cercano a lo que hacemos en la fase 6b
  (esfera elástica sumergida) y a lo que viene con la cóclea.
- <https://gitlab.com/openCFS/cfs> — openCFS (ex openCFS): acústica y vibroacústica con acoplamiento,
  pensado para ingeniería. Referencia de formulación y de casos de validación publicados.

## Para qué los usaría, concretamente

1. **El 22 % de la fase 6b** (D46): el sospechoso es la interfaz fluido–sólido. Elmer y VIBRA son los
   que tienen esa formulación escrita y verificada; sirven para comparar los **dos términos** del
   acoplamiento (masa de superficie proyectada sobre la normal y su transpuesta) y el factor `ω²`.
2. **La superficie facetada**: deal.II y FEniCSx tienen elementos de frontera curva (isoparamétricos).
   Si el error resulta ser O(h) por el facetado, la salida no es refinar más fino sino curvar la
   frontera.
3. **Los benchmarks**: el manual de ABEM y VIBRA publican casos con respuesta conocida (esfera
   elástica en agua, cilindro, placa). Son candidatos a criterio nuevo en `validacion.py`, con la
   cita correspondiente.

## Lo que ya tenemos y no hace falta buscar afuera

- Esfera fluida y rígida: `analytic/sphere.py` (Anderson 1950), validada por continuidad y límites.
- Esfera elástica: `analytic/faran.py` (Faran 1951), validada por unitariedad a 1e-16 y teorema
  óptico (D45).
- Condiciones absorbentes 1 y 2: `fem/helmholtz3d.py` (ABC-1 y BGT-2).
