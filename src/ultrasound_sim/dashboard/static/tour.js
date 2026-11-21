/* Tour guiado del dashboard.
 *
 * Por qué existe: la interfaz tiene ocho vistas, tres barridos y un registro de criterios, y casi todo
 * eso es invisible hasta que alguien te lo cuenta. El tour recorre los controles que importan, resalta
 * el elemento de turno y explica para qué está.
 *
 * Dos decisiones de diseño:
 *  - NO cambia de vista para resaltar las otras pestañas: resalta el botón del nav. Cambiar de vista
 *    dispara solves (la esfera precisa son 30 s) y un tour que lanza simulaciones para explicarse es
 *    peor que no tener tour.
 *  - Se muestra una sola vez (`localStorage`), y queda un botón para repetirlo cuando quieras.
 */

const TOUR_CLAVE = 'usim.tour.visto';

const TOUR_PASOS = [
  {
    sel: null,
    titulo: 'Qué es esto',
    texto: 'Un simulador de ultrasonido con dos métodos numéricos —BEM y FEM— sobre el mismo problema, '
      + 'para poder contrastarlos. Todo lo que se ve acá lo calcula el servidor: la página no hace '
      + 'física, sólo pide y dibuja. Este recorrido muestra los controles que importan.',
  },
  {
    sel: '#mapa',
    titulo: 'El campo',
    texto: 'El mapa de presión (dB) con los dos cortes: axial y lateral. El foco pedido y el pico real '
      + 'no siempre coinciden, y cuando no coinciden el dashboard lo dice arriba con el motivo '
      + '(corrimiento no paraxial, o foco más allá del campo cercano).',
  },
  {
    // `#arrastre` es un overlay con el mismo rect que el mapa: resaltarlo no se nota. Se resalta el
    // hint, que es el texto que explica el gesto.
    sel: '#hintMover',
    titulo: 'Mover las fuentes con el mouse',
    texto: 'Arrastrá sobre el mapa y se mueve la fuente más cercana. En modo «arreglo» no hay fuentes '
      + 'sueltas: lo que se mueve es el FOCO en y, o sea el haz se inclina. Las asas cuadradas marcan '
      + 'dónde está cada fuente o elemento. Doble clic sobre el mapa agrega una fuente.',
  },
  {
    sel: '#btnMover',
    titulo: 'Mover o hacer zoom',
    texto: 'Con «mover con el mouse» activo, arrastrar mueve y no hace zoom. Si preferís el zoom de '
      + 'Plotly, cambiá a «zoom»: el arrastre vuelve a acercar. Esta era la confusión más común: '
      + 'arrastrar hacía zoom y parecía que las fuentes no se podían mover.',
  },
  {
    sel: '#freq',
    titulo: 'Los parámetros',
    texto: 'Frecuencia, cantidad de elementos, apertura y foco. La cantidad de elementos es la que '
      + 'manda en la calidad del haz: el barrido de apertura (pestaña «barridos») mide que para no '
      + 'despegarse de una apertura continua hace falta un paso de λ/4.',
  },
  {
    sel: '#btnRefFijar',
    titulo: 'Comparar dos configuraciones',
    texto: 'Fijá el campo actual como referencia, cambiá lo que quieras y mirá la DIFERENCIA en dB en '
      + 'vez de dos mapas de a uno. Es el uso central: mover un transductor y ver qué cambió.',
  },
  {
    sel: 'nav.vistas',
    titulo: 'Las ocho vistas',
    texto: 'campo (aperturas 2D) · escena FEM (capas contra la matriz de transferencia) · dispersor '
      + 'BEM (el mismo cilindro por BEM, FEM y la serie exacta) · cóclea · barridos · 3D · esfera '
      + 'elástica (sólido con ondas de corte contra Faran) · criterios (el registro de validación).',
  },
  {
    sel: '#btnExportarVista',
    titulo: 'Exportar',
    texto: 'Baja un PNG de cada gráfica visible y un JSON con sus trazas. Funciona igual en las ocho '
      + 'vistas porque lee lo que está en pantalla.',
  },
  {
    sel: null,
    titulo: 'Listo',
    texto: 'Eso es todo. Cada vista guarda su último resultado: volver a una pestaña no cuesta otro '
      + 'solve. En «criterios» está el registro completo con la tolerancia y la fuente de cada '
      + 'criterio, y desde ahí se descarga el informe de la sesión. Podés repetir este tour con el '
      + 'botón «tour» de arriba.',
  },
];

const TOUR = { i: 0, activo: false };

function tour$ (sel) { return sel ? document.querySelector(sel) : null; }

function tourCaja() {
  let caja = document.getElementById('tourCaja');
  if (!caja) {
    caja = document.createElement('div');
    caja.id = 'tourCaja';
    caja.className = 'tour-caja';
    caja.innerHTML = `
      <div class="tour-paso"></div>
      <h3 class="tour-titulo"></h3>
      <p class="tour-texto"></p>
      <div class="tour-botones">
        <button class="btn mini" data-tour="saltear">saltear</button>
        <span class="tour-flex"></span>
        <button class="btn mini" data-tour="atras">atrás</button>
        <button class="btn mini" data-tour="siguiente">siguiente</button>
      </div>`;
    caja.querySelector('[data-tour="saltear"]').addEventListener('click', () => tourCerrar());
    caja.querySelector('[data-tour="atras"]').addEventListener('click', () => tourIr(TOUR.i - 1));
    caja.querySelector('[data-tour="siguiente"]').addEventListener('click', () => tourIr(TOUR.i + 1));
    document.body.appendChild(caja);
  }
  return caja;
}

function tourFoco() {
  let foco = document.getElementById('tourFoco');
  if (!foco) {
    foco = document.createElement('div');
    foco.id = 'tourFoco';
    foco.className = 'tour-foco';
    document.body.appendChild(foco);
  }
  return foco;
}

function tourIr(i) {
  if (i < 0) return;
  if (i >= TOUR_PASOS.length) { tourCerrar(); return; }
  TOUR.i = i;
  const paso = TOUR_PASOS[i];
  const caja = tourCaja();
  const foco = tourFoco();
  const el = tour$(paso.sel);

  caja.querySelector('.tour-paso').textContent = `${i + 1} / ${TOUR_PASOS.length}`;
  caja.querySelector('.tour-titulo').textContent = paso.titulo;
  caja.querySelector('.tour-texto').textContent = paso.texto;
  caja.querySelector('[data-tour="atras"]').disabled = i === 0;
  caja.querySelector('[data-tour="siguiente"]').textContent =
    i === TOUR_PASOS.length - 1 ? 'terminar' : 'siguiente';

  if (!el) {
    // paso sin elemento (bienvenida y cierre): se centra la caja y se apaga el resaltado
    foco.style.display = 'none';
    caja.classList.add('centrada');
    caja.style.left = '';
    caja.style.top = '';
    return;
  }
  caja.classList.remove('centrada');
  foco.style.display = 'block';
  el.scrollIntoView({ block: 'nearest', inline: 'nearest' });
  const r = el.getBoundingClientRect();
  const pad = 6;
  foco.style.left = `${r.left - pad}px`;
  foco.style.top = `${r.top - pad}px`;
  foco.style.width = `${r.width + 2 * pad}px`;
  foco.style.height = `${r.height + 2 * pad}px`;

  // la caja va debajo del elemento si entra; si no, arriba
  const ac = caja.getBoundingClientRect();
  let x = Math.min(Math.max(12, r.left), window.innerWidth - ac.width - 12);
  let y = r.bottom + 14;
  if (y + ac.height > window.innerHeight - 12) y = Math.max(12, r.top - ac.height - 14);
  caja.style.left = `${x}px`;
  caja.style.top = `${y}px`;
}

window.tourAbrir = function tourAbrir() {
  TOUR.activo = true;
  document.body.classList.add('tour-abierto');
  tourIr(0);
  window.addEventListener('resize', () => { if (TOUR.activo) tourIr(TOUR.i); });
};

function tourCerrar() {
  TOUR.activo = false;
  document.body.classList.remove('tour-abierto');
  const caja = document.getElementById('tourCaja');
  const foco = document.getElementById('tourFoco');
  if (caja) caja.remove();
  if (foco) foco.remove();
  try { localStorage.setItem(TOUR_CLAVE, '1'); } catch { /* nada */ }
}

document.addEventListener('keydown', (ev) => {
  if (!TOUR.activo) return;
  if (ev.key === 'Escape') tourCerrar();
  if (ev.key === 'ArrowRight') tourIr(TOUR.i + 1);
  if (ev.key === 'ArrowLeft') tourIr(TOUR.i - 1);
});

// primera visita: se abre solo. Después, sólo con el botón.
document.addEventListener('DOMContentLoaded', () => {
  const boton = document.getElementById('btnTour');
  if (boton) boton.addEventListener('click', () => window.tourAbrir());
  let visto = null;
  try { visto = localStorage.getItem(TOUR_CLAVE); } catch { /* nada */ }
  if (!visto) setTimeout(() => window.tourAbrir(), 1200);
});
