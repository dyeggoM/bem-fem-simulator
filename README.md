# ultrasound-simulator

Simulador de ultrasonido con dos métodos numéricos, para poder comparar y verificar resultados:

- **BEM** (boundary element method) — superficie/contorno, medio homogéneo por regiones.
- **FEM** (finite element method) — dominio completo, permite heterogeneidades e interfaces.

Objetivo: simular **múltiples fuentes de ultrasonido**, en distintos **materiales y combinaciones**,
y contrastar BEM contra FEM para validar la física (impedancia, reflexión/transmisión, atenuación,
campo de presión, tiempos de vuelo).

Python. La especificación completa llega en el documento de plan que Diego comparte en el chat
`Simulador`; hasta entonces, `docs/` es el sitio donde se anota lo acordado.
