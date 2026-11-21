/* Dashboard BEM/FEM — cliente. Sin dependencias salvo Plotly, que se sirve local.
 *
 * Todo el estado vive en S; cada control llama a `pedir()`, que debouncea y aborta el pedido
 * anterior. El servidor devuelve el campo ya en dB y las métricas: acá no se calcula física.
 */

const S = {
  material: 'agua',
  modelo: 'baffle',
  freq: 250,
  modo: 'array',
  n: 32,
  d: 60,
  fx: 90,
  fy: 0,
  grid: 120,
  fuentes: null,          // null = el servidor arma el arreglo
  ultimo: null,
};

const $ = (id) => document.getElementById(id);
const fmt = (v, d = 2) => (v === null || v === undefined || Number.isNaN(v) ? '—' : v.toFixed(d));

const EJE = {
  gridcolor: '#232b36', zerolinecolor: '#2e3746', linecolor: '#232b36',
  tickfont: { size: 10, color: '#8b95a1' }, titlefont: { size: 11, color: '#8b95a1' },
};
const BASE = {
  paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
  font: { color: '#d7dde5', size: 11 },
  margin: { l: 52, r: 16, t: 26, b: 40 },
  hoverlabel: { bgcolor: '#1b2129', bordercolor: '#2e3746', font: { size: 11 } },
};
const CONFIG = { displayModeBar: true, displaylogo: false, responsive: true,
                 modeBarButtonsToRemove: ['select2d', 'lasso2d', 'autoScale2d'] };

let temporizador = null;
let control = null;
let arrastre = null;        // {i, y0} mientras se arrastra un asa
const asas = [];            // elementos .asa, uno por fuente (se reposicionan, no se recrean)
const BORDE = 1;            // borde de 1 px del div de Plotly: las coordenadas internas no lo cuentan

function pedir() {
  clearTimeout(temporizador);
  temporizador = setTimeout(enviar, 150);
}

/* ---------------------------------------------------------------------------------------------
 * Persistencia del último resultado de una vista.
 *
 * Las vistas caras (esfera precisa ~30 s, escena, volumen, cóclea) pierden el resultado al recargar
 * o al cambiar de pestaña y volver. Guardarlo en localStorage deja volver a mirarlo sin volver a
 * resolver — que es la diferencia entre usar el dashboard para explorar y usarlo para una sola foto.
 *
 * Se guarda por vista y con la fecha, así el usuario ve de cuándo es lo que está mirando. Si el
 * resultado no entra en la cuota, se descarta en silencio y se avisa por consola: nunca se rompe la
 * vista por no poder cachear.
 * ------------------------------------------------------------------------------------------- */

const VISTA_CLAVE = 'usim.vista.';
const VISTA_TTL_MS = 12 * 60 * 60 * 1000;   // 12 h: un resultado de ayer no debería parecer de ahora

window.vistaGuardar = function vistaGuardar(vista, resultado) {
  try {
    localStorage.setItem(VISTA_CLAVE + vista, JSON.stringify({
      fecha: new Date().toISOString(), resultado,
    }));
    return true;
  } catch (e) {
    try { localStorage.removeItem(VISTA_CLAVE + vista); } catch { /* nada */ }
    console.warn(`[${vista}] no se pudo guardar el resultado (cuota):`, e && e.message);
    return false;
  }
};

window.vistaCargar = function vistaCargar(vista) {
  try {
    const crudo = localStorage.getItem(VISTA_CLAVE + vista);
    if (!crudo) return null;
    const datos = JSON.parse(crudo);
    if (!datos || !datos.resultado) return null;
    const edad = Date.now() - Date.parse(datos.fecha);
    if (!Number.isFinite(edad) || edad > VISTA_TTL_MS) {
      localStorage.removeItem(VISTA_CLAVE + vista);
      return null;
    }
    return datos;
  } catch {
    return null;
  }
};

window.vistaBorrar = function vistaBorrar(vista) {
  try { localStorage.removeItem(VISTA_CLAVE + vista); } catch { /* nada */ }
};

/* ---------------------------------------------------------------------------------------------
 * Exportación de la vista activa.
 *
 * No sabe nada de la física ni del estado de cada vista: agarra las gráficas Plotly que están
 * VISIBLES en el DOM y las baja. Así funciona igual en las siete vistas y no hay que mantener siete
 * exportadores en paralelo (que es como se desincronizan).
 *
 * Qué baja:
 *   - un PNG por gráfica, con el nombre del contenedor para poder ordenarlos;
 *   - un JSON con las trazas completas (x, y, nombre, tipo) y el layout, que es lo que permite
 *     re-graficar o comparar sin volver a resolver.
 * ------------------------------------------------------------------------------------------- */

function contenedoresVisibles() {
  return [...document.querySelectorAll('.js-plotly-plot')].filter((gd) => {
    if (!gd || !gd.data || !gd.data.length) return false;
    const r = gd.getBoundingClientRect();
    return r.width > 20 && r.height > 20;   // una pestaña oculta mide 0
  });
}

function nombreDeGrafica(gd, i) {
  return gd.id || gd.parentElement.id || `grafica_${i + 1}`;
}

function descargarBlob(nombre, blob) {
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = nombre;
  document.body.appendChild(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 2000);
}

window.exportarVista = async function exportarVista() {
  const boton = document.getElementById('btnExportarVista');
  const graficas = contenedoresVisibles();
  if (!graficas.length) {
    if (boton) boton.textContent = 'no hay gráficas visibles';
    setTimeout(() => { if (boton) boton.textContent = 'exportar vista'; }, 2500);
    return;
  }
  if (boton) { boton.disabled = true; boton.textContent = 'exportando…'; }
  const vista = document.querySelector('.vista-btn.activo');
  const sello = (vista ? vista.dataset.vista : 'vista') + '_' + new Date()
    .toISOString().slice(0, 19).replace(/[:T]/g, '-');

  const trazas = [];
  for (let i = 0; i < graficas.length; i += 1) {
    const gd = graficas[i];
    const nombre = nombreDeGrafica(gd, i);
    trazas.push({ grafica: nombre, layout: { title: gd.layout && gd.layout.title }, datos: gd.data });
    try {
      const url = await Plotly.toImage(gd, { format: 'png', width: 1400, height: Math.max(500,
        Math.round(1400 * (gd.getBoundingClientRect().height / Math.max(1, gd.getBoundingClientRect().width)))) });
      const blob = await (await fetch(url)).blob();
      descargarBlob(`${sello}__${nombre}.png`, blob);
    } catch (e) {
      trazas.push({ grafica: nombre, error_png: String(e && e.message ? e.message : e) });
    }
  }
  descargarBlob(`${sello}__trazas.json`,
    new Blob([JSON.stringify({ formato: 'ultrasound-sim/dashboard/vista/1',
      fecha: new Date().toISOString(), vista: vista ? vista.dataset.vista : null, graficas: trazas },
    null, 2)], { type: 'application/json' }));

  if (boton) { boton.disabled = false; boton.textContent = 'exportar vista'; }
};

async function enviar() {
  if (control) control.abort();
  control = new AbortController();
  const xmax = Math.max(0.06, S.fx / 1000 * 1.6);
  const ym = Math.max(0.02, S.d / 1000 * 1.1);
  // con el modelo libre la malla se extiende al semiespacio trasero: ahí está lo que el baffle niega
  const x0 = S.modelo === 'libre' ? -xmax * 0.62 : 0.004;
  const cuerpo = {
    material: S.material,
    modelo: S.modelo,
    frequency_khz: S.freq,
    // mientras se arrastra un asa se baja la malla: 80×80 tarda ~0.4 s y el jalón se sigue
    grid: { x0, x1: xmax, y0: -ym, y1: ym,
            nx: arrastre ? 80 : S.grid, ny: arrastre ? 80 : S.grid },
  };
  if (S.modo === 'manual' && S.fuentes) {
    cuerpo.sources = S.fuentes;
  } else {
    cuerpo.n_elements = S.n;
    cuerpo.aperture_mm = S.d;
    cuerpo.focus_x_mm = S.fx;
    cuerpo.focus_y_mm = S.fy;
  }
  $('estado').className = 'pie';
  $('estado').textContent = 'calculando…';
  try {
    const r = await fetch('/api/field', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(cuerpo), signal: control.signal,
    });
    const datos = await r.json();
    if (!r.ok) throw new Error(datos.error || `HTTP ${r.status}`);
    S.ultimo = datos;
    pintar(datos);
    $('estado').textContent =
      `${datos.meta.n_sources} fuentes · ${datos.meta.modelo} · ${datos.meta.nx}×${datos.meta.ny} · ` +
      `λ = ${fmt(datos.meta.lambda_mm)} mm · ${fmt(datos.meta.seconds, 2)} s`;
  } catch (e) {
    if (e.name === 'AbortError') return;
    $('estado').className = 'pie error';
    $('estado').textContent = e.message;
  }
}

function pintar(r) {
  const m = r.metrics;
  $('mPico').textContent = `${fmt(m.peak_x_mm, 1)} · ${fmt(m.peak_y_mm, 1)} mm`;
  $('mFwhm').textContent = m.fwhm_lateral_mm ? `${fmt(m.fwhm_lateral_mm)} mm` : '—';
  $('mFwhmAx').textContent = m.fwhm_axial_mm ? `${fmt(m.fwhm_axial_mm)} mm` : '—';
  $('mLobe').textContent = m.sidelobe_db !== null ? `${fmt(m.sidelobe_db, 1)} dB` : '—';
  $('mGrat').textContent = m.grating_deg !== null ? `${fmt(m.grating_deg, 1)}°` : '—';
  $('mTime').textContent = `${fmt(r.meta.seconds, 2)} s`;

  // El foco geométrico no siempre es donde está el pico: si el pedido está más allá del campo
  // cercano (D²/4λ) la apertura no puede enfocar ahí. Mejor decirlo que dejar el mapa mintiendo.
  const aviso = $('aviso');
  const foco = r.meta.focus_x_mm;
  const campoCercano = r.meta.near_field_mm;
  if (foco && Math.abs(m.peak_x_mm - foco) > 0.15 * foco) {
    aviso.hidden = false;
    if (campoCercano && foco > campoCercano) {
      aviso.textContent =
        `El pico real está en ${fmt(m.peak_x_mm, 1)} mm, no en F = ${fmt(foco, 0)} mm: ` +
        `el campo cercano de esta apertura termina en ${fmt(campoCercano, 0)} mm y más allá no enfoca.`;
    } else {
      aviso.textContent =
        `El pico real está en ${fmt(m.peak_x_mm, 1)} mm, no en F = ${fmt(foco, 0)} mm: ` +
        'es el corrimiento focal no paraxial (D13), no un error del solver.';
    }
  } else {
    aviso.hidden = true;
  }

  const fuentes = (S.modo === 'manual' && S.fuentes)
    ? S.fuentes.map((f) => ({ x: 0, y: f.y * 1e3 }))
    : (() => {                       // arreglo: reconstruyo las posiciones para dibujarlas
        const paso = S.d / S.n, out = [];
        for (let j = 0; j < S.n; j += 1) out.push({ x: 0, y: -S.d / 2 + paso * (j + 0.5) });
        return out;
      })();

  const trazas = [{
    type: 'heatmap', x: r.x_mm, y: r.y_mm, z: r.p2_db,
    colorscale: 'Inferno', zmin: -40, zmax: 0, zhoverformat: '.1f',
    colorbar: { title: { text: 'dB', side: 'top' }, thickness: 11, len: 0.8,
                tickfont: { size: 9, color: '#8b95a1' }, outlinewidth: 0 },
    hovertemplate: 'x %{x:.1f} mm<br>y %{y:.1f} mm<br>%{z:.1f} dB<extra></extra>',
  }, {
    type: 'scatter', x: fuentes.map((f) => f.x), y: fuentes.map((f) => f.y),
    mode: 'markers', marker: { symbol: 'square', size: 6, color: '#58a6ff',
                               line: { width: 0.5, color: '#0d1117' } },
    hovertemplate: 'fuente y %{y:.2f} mm<extra></extra>', name: 'fuentes',
  }];
  // margen a la izquierda para que las fuentes sobre el baffle (x=0) no queden pisadas por el eje
  const xmax = r.x_mm[r.x_mm.length - 1];
  // Aspecto 1:1 sin `scaleanchor` (que con rangos fijos aborta con "Something went wrong with axis
  // scaling"): se elige la escala px/mm mayor de las dos y se centran ambos rangos en el dato.
  const caja = $('mapa');
  const xr = [-0.06 * xmax, xmax];
  const yr0 = [r.y_mm[0], r.y_mm[r.y_mm.length - 1]];
  const gleich = (wpx, hpx) => {                        // rangos con igual mm/px
    const escala = Math.max((xr[1] - xr[0]) / wpx, (yr0[1] - yr0[0]) / hpx);
    const cx = (xr[0] + xr[1]) / 2;
    const cy = (yr0[0] + yr0[1]) / 2;
    return { x: [cx - escala * wpx / 2, cx + escala * wpx / 2],
             y: [cy - escala * hpx / 2, cy + escala * hpx / 2] };
  };
  const rangoM = gleich(Math.max(120, caja.clientWidth - 62), Math.max(120, caja.clientHeight - 66));
  const xrango = rangoM.x;
  const yrango = rangoM.y;

  Plotly.react('mapa', trazas, {
    ...BASE, margin: { l: 54, r: 8, t: 26, b: 40 },
    title: { text: `${r.meta.material} · ${fmt(r.meta.frequency_khz, 0)} kHz`, font: { size: 12 },
             x: 0.01, y: 0.98, xanchor: 'left' },
    // Sin `scaleanchor`: con rangos fijos, Plotly no puede satisfacer el aspecto 1:1 dentro del
    // panel y aborta con "Something went wrong with axis scaling".
    xaxis: { ...EJE, title: 'x [mm]', range: xrango },
    yaxis: { ...EJE, title: 'y [mm]', range: yrango },
    shapes: [
      { type: 'line', x0: m.peak_x_mm, x1: m.peak_x_mm, yref: 'paper', y0: 0, y1: 1,
        line: { color: '#58a6ff', width: 1, dash: 'dot' } },
      { type: 'line', y0: m.peak_y_mm, y1: m.peak_y_mm, xref: 'paper', x0: 0, x1: 1,
        line: { color: '#58a6ff', width: 1, dash: 'dot' } },
    ],
    showlegend: false,
  }, CONFIG).then(() => {
    // Los márgenes reales no son los estimados (la barra de color se come ~65 px del área de
    // trazado), así que se compara el px/mm **medido** en cada eje y se corrige si difieren.
    const l = caja._fullLayout;
    if (!l || !l.xaxis || !l.yaxis) return;
    const ancho = l.xaxis._length;
    const alto = l.yaxis._length;
    if (!ancho || !alto || !l.xaxis.range || !l.yaxis.range) return;
    const pxmmX = ancho / (l.xaxis.range[1] - l.xaxis.range[0]);
    const pxmmY = alto / (l.yaxis.range[1] - l.yaxis.range[0]);
    if (Math.abs(pxmmX - pxmmY) > 0.02 * pxmmY) {
      const exacto = gleich(ancho, alto);
      Plotly.relayout('mapa', { xaxis: { range: exacto.x }, yaxis: { range: exacto.y } });
    }
    pintarAsas();
  });

  // p2_db es [fila=y][columna=x] (igual que en Python): axial fija y y varía x; lateral fija x.
  // El índice del pico se busca con un argmax sobre la matriz y no comparando coordenadas:
  // comparar floats redondeados fallaba y la línea caía en la columna 0 (−6.7 dB en vez de 0).
  let biy = 0;
  let bix = 0;
  let mejor = -Infinity;
  r.p2_db.forEach((fila, i) => fila.forEach((v, j) => {
    if (v > mejor) { mejor = v; biy = i; bix = j; }
  }));
  const axial = r.p2_db[biy];                                  // varía x
  const lateral = r.p2_db.map((fila) => fila[bix]);            // varía y

  // Rango automático por corte: en escala fija −40..0 dB un perfil axial que varía 5 dB se ve plano.
  const rango = (valores) => {
    const lo = Math.max(-45, Math.floor(Math.min(...valores)) - 3);
    return [lo, 2];
  };
  const lineaDb = (v) => ({ type: 'line', xref: 'paper', x0: 0, x1: 1, y0: -6, y1: -6,
                            line: { color: '#4a5566', width: 1, dash: 'dot' } });

  Plotly.react('axial', [{
    type: 'scatter', mode: 'lines', x: r.x_mm, y: axial, line: { color: '#58a6ff', width: 1.6 },
    hovertemplate: 'x %{x:.1f} mm<br>%{y:.1f} dB<extra></extra>',
  }], {
    ...BASE, margin: { l: 46, r: 12, t: 24, b: 34 },
    title: { text: 'corte axial (y = pico)', font: { size: 11 }, x: 0.01, y: 0.95, xanchor: 'left' },
    xaxis: { ...EJE, title: 'x [mm]' }, yaxis: { ...EJE, title: 'dB', range: rango(axial) },
    shapes: [lineaDb()],
  }, CONFIG);

  Plotly.react('lateral', [{
    type: 'scatter', mode: 'lines', x: r.y_mm, y: lateral, line: { color: '#3fb950', width: 1.6 },
    hovertemplate: 'y %{x:.1f} mm<br>%{y:.1f} dB<extra></extra>',
  }], {
    ...BASE, margin: { l: 46, r: 12, t: 24, b: 34 },
    title: { text: 'corte lateral en el pico', font: { size: 11 }, x: 0.01, y: 0.95, xanchor: 'left' },
    xaxis: { ...EJE, title: 'y [mm]' }, yaxis: { ...EJE, title: 'dB', range: rango(lateral) },
    shapes: [lineaDb()],
  }, CONFIG);
}

/* ---------------------------------------------------------------- controles */
function sincronizarEtiquetas() {
  $('freqVal').textContent = `${S.freq} kHz`;
  $('nVal').textContent = S.n;
  $('apVal').textContent = `${S.d} mm`;
  $('fxVal').textContent = `${S.fx} mm`;
  $('fyVal').textContent = `${S.fy} mm`;
  const c = { agua: 1500, tejido: 1540, hueso: 3500 }[S.material] || 1500;
  $('lamHint').textContent = `λ = ${(c / (S.freq * 1e3) * 1e3).toFixed(2)} mm (aprox.)`;
}

function filaFuente(src, i) {
  const tr = document.createElement('tr');
  const campos = [
    { k: 'y_mm', v: (src.y * 1e3).toFixed(2), paso: '0.1' },
    { k: 'width_mm', v: (src.width * 1e3).toFixed(2), paso: '0.1' },
    { k: 'amp', v: src.amp.toFixed(2), paso: '0.05' },
    { k: 'phase', v: src.phase.toFixed(1), paso: '0.1' },
  ];
  for (const c of campos) {
    const td = document.createElement('td');
    const inp = document.createElement('input');
    inp.type = 'number'; inp.step = c.paso; inp.value = c.v;
    inp.addEventListener('change', () => {
      const v = parseFloat(inp.value);
      if (Number.isNaN(v)) return;
      const f = S.fuentes[i];
      if (c.k === 'y_mm') f.y = v / 1e3;
      else if (c.k === 'width_mm') {
        if (v <= 0) { inp.value = (f.width * 1e3).toFixed(2); return; }
        f.width = v / 1e3;
      } else if (c.k === 'amp') f.amp = v;
      else f.phase = v;
      pedir();
    });
    td.appendChild(inp);
    tr.appendChild(td);
  }
  const td = document.createElement('td');
  const del = document.createElement('button');
  del.textContent = '×'; del.title = 'quitar';
  del.addEventListener('click', () => { S.fuentes.splice(i, 1); pintarFuentes(); pedir(); });
  td.appendChild(del);
  tr.appendChild(td);
  return tr;
}

function pintarFuentes() {
  const cuerpo = $('srcBody');
  cuerpo.innerHTML = '';
  (S.fuentes || []).forEach((f, i) => cuerpo.appendChild(filaFuente(f, i)));
}

function generarDesdeArreglo() {
  const paso = S.d / S.n, arr = [];
  const c = { agua: 1500, tejido: 1540, hueso: 3500 }[S.material] || 1500;
  const k = 2 * Math.PI * (S.freq * 1e3) / c;
  for (let j = 0; j < S.n; j += 1) {
    const y = -S.d / 2 + paso * (j + 0.5);
    arr.push({ y: y / 1e3, width: paso / 1e3, amp: 1,
               phase: -k * Math.hypot(S.fx, S.fy - y) / 1e3 });
  }
  S.fuentes = arr;
  pintarFuentes();
}

function cambiarModo(modo) {
  S.modo = modo;
  document.querySelectorAll('.tab').forEach((t) => t.classList.toggle('activo', t.dataset.modo === modo));
  $('ctrlArray').hidden = modo !== 'array';
  $('ctrlManual').hidden = modo !== 'manual';
  if (modo === 'manual' && !S.fuentes) generarDesdeArreglo();
  pintarFuentes();
  pedir();
}

document.querySelectorAll('.tab').forEach((t) => t.addEventListener('click', () => cambiarModo(t.dataset.modo)));

const enlazar = (id, clave, escala = 1) => {
  $(id).addEventListener('input', (e) => {
    S[clave] = clave === 'material' ? e.target.value : parseFloat(e.target.value) * escala;
    if (clave === 'n' || clave === 'd') S.fuentes = null;   // el arreglo manda de nuevo
    sincronizarEtiquetas(); pintarFuentes(); pedir();
  });
};
enlazar('freq', 'freq');
enlazar('nEl', 'n');
enlazar('ap', 'd');
enlazar('fx', 'fx');
enlazar('fy', 'fy');
$('mat').addEventListener('change', (e) => { S.material = e.target.value; sincronizarEtiquetas(); pedir(); });
$('model').addEventListener('change', (e) => {
  S.modelo = e.target.value;
  $('modeloHint').textContent = S.modelo === 'libre'
    ? 'Sin baffle: el campo también existe detrás. Los dos lóbulos son igual de fuertes '
      + '(180° de fase), así que el máximo global puede caer del lado de atrás.'
    : 'El baffle niega el semiespacio trasero: con este modelo no hay campo en x < 0.';
  pedir();
});
$('res').addEventListener('change', (e) => { S.grid = parseInt(e.target.value, 10); pedir(); });
$('btnManual').addEventListener('click', () => cambiarModo('manual'));
$('btnArray').addEventListener('click', () => { generarDesdeArreglo(); pedir(); });
$('btnAdd').addEventListener('click', () => {
  const ult = S.fuentes && S.fuentes.length ? S.fuentes[S.fuentes.length - 1] : null;
  const ancho = ult ? ult.width : (S.d / S.n) / 1e3;
  const y = ult ? ult.y + 2 * ancho : 0;
  S.fuentes.push({ y, width: ancho, amp: 1, phase: 0 });
  pintarFuentes(); pedir();
});


/* ---------------------------------------------------------------- asas arrastrables
 * Las fuentes viven sobre el baffle (x = 0), así que arrastrarlas es 1-D: sólo cambia su `y`.
 * El overlay es propio (no Plotly), y las asas se *reposicionan* en cada render en vez de
 * recrearse: si se recrearan, el pointer capture del arrastre se perdería a mitad del jalón.
 */
function pxDeY(l, y) {
  const a = l.yaxis;
  return BORDE + a._offset + (a.range[1] - y) / (a.range[1] - a.range[0]) * a._length;
}
function pxDeX(l, x) {
  const a = l.xaxis;
  return BORDE + a._offset + (x - a.range[0]) / (a.range[1] - a.range[0]) * a._length;
}
function yDePx(l, py) {
  const a = l.yaxis;
  const f = (py - BORDE - a._offset) / a._length;          // 0 = arriba
  return a.range[1] + (a.range[0] - a.range[1]) * f;
}

function pintarAsas() {
  const capa = $('arrastre');
  const l = $('mapa')._fullLayout;
  if (!l || !l.xaxis || !l.yaxis) return;
  const fuentes = S.modo === 'manual' && S.fuentes ? S.fuentes : null;
  const n = fuentes ? fuentes.length : 0;

  while (asas.length > n) capa.removeChild(asas.pop());
  for (let i = 0; i < n; i += 1) {
    let el = asas[i];
    if (!el) {
      el = document.createElement('div');
      el.className = 'asa';
      el.dataset.i = String(i);
      el.addEventListener('pointerdown', iniciarArrastre);
      capa.appendChild(el);
      asas[i] = el;
    }
    const x = pxDeX(l, 0);
    const y = pxDeY(l, fuentes[i].y * 1e3);
    el.style.left = `${x}px`;
    el.style.top = `${y}px`;
    el.title = `fuente ${i + 1}: y = ${(fuentes[i].y * 1e3).toFixed(2)} mm`;
  }
}

function iniciarArrastre(ev) {
  if (S.modo !== 'manual' || !S.fuentes) {
    $('estado').textContent = 'pasá a la pestaña "fuentes" para arrastrarlas';
    return;
  }
  ev.preventDefault();
  const i = parseInt(ev.currentTarget.dataset.i, 10);
  arrastre = { i, y0: S.fuentes[i].y };
  ev.currentTarget.classList.add('activa');
  // La captura es opcional: con pointers sintéticos (o navegadores viejos) tira y no hay que
  // romper el jalón por eso. Los listeners van en window para que el asa siga al puntero aunque
  // se salga del cuadrado.
  try { ev.currentTarget.setPointerCapture(ev.pointerId); } catch (e) { /* sin captura */ }
  window.addEventListener('pointermove', moverArrastre);
  window.addEventListener('pointerup', finArrastre);
  window.addEventListener('pointercancel', finArrastre);
}

function moverArrastre(ev) {
  if (!arrastre) return;
  const l = $('mapa')._fullLayout;
  if (!l || !l.yaxis) return;
  const caja = $('mapa').getBoundingClientRect();
  const py = ev.clientY - caja.top;
  let ymm = yDePx(l, py);
  const radio = (l.yaxis.range[1] - l.yaxis.range[0]) / 2;
  ymm = Math.max(-radio, Math.min(radio, ymm));          // no se sale del encuadre
  S.fuentes[arrastre.i].y = ymm / 1e3;
  const el = asas[arrastre.i];
  if (el) {
    el.style.top = `${pxDeY(l, ymm)}px`;
    el.title = `fuente ${arrastre.i + 1}: y = ${ymm.toFixed(2)} mm`;
  }
  pedir();                                               // recomputa con malla reducida
}

function finArrastre() {
  if (!arrastre) return;
  const el = asas[arrastre.i];
  if (el) el.classList.remove('activa');
  window.removeEventListener('pointermove', moverArrastre);
  window.removeEventListener('pointerup', finArrastre);
  window.removeEventListener('pointercancel', finArrastre);
  arrastre = null;
  pintarFuentes();
  enviar();                                              // vuelve a la malla elegida
}

/* ---------------------------------------------------------------- agregar con doble clic */
$('mapa').addEventListener('dblclick', (ev) => {
  const l = $('mapa')._fullLayout;
  if (!l || !l.yaxis) return;
  if (S.modo !== 'manual') cambiarModo('manual');
  const caja = $('mapa').getBoundingClientRect();
  const ymm = yDePx(l, ev.clientY - caja.top);
  const ancho = S.fuentes && S.fuentes.length ? S.fuentes[S.fuentes.length - 1].width : (S.d / S.n) / 1e3;
  S.fuentes.push({ y: ymm / 1e3, width: ancho, amp: 1, phase: 0 });
  pintarFuentes();
  pedir();
});

/* ---------------------------------------------------------------- guardar / cargar */
function instantanea() {
  return {
    formato: 'ultrasound-sim/config',
    version: 1,
    material: S.material,
    modelo: S.modelo,
    frequency_khz: S.freq,
    modo: S.modo,
    n_elements: S.n,
    aperture_mm: S.d,
    focus_x_mm: S.fx,
    focus_y_mm: S.fy,
    grid: S.grid,
    sources: S.modo === 'manual' && S.fuentes ? S.fuentes : null,
    // los tres parámetros de la cóclea viajan con la configuración: son parte de lo que el usuario fijó
    coclea: { thickness_mm: COC.espesor, c_s: COC.c_s, alpha_np_m: COC.alpha },
  };
}

function aplicar(cfg) {
  if (!cfg || cfg.formato !== 'ultrasound-sim/config') {
    throw new Error('no es un archivo de configuración del simulador');
  }
  S.material = cfg.material || S.material;
  S.modelo = cfg.modelo === 'libre' ? 'libre' : 'baffle';
  S.freq = Number(cfg.frequency_khz) || S.freq;
  S.n = Number(cfg.n_elements) || S.n;
  S.d = Number(cfg.aperture_mm) || S.d;
  S.fx = Number(cfg.focus_x_mm) || S.fx;
  S.fy = Number(cfg.focus_y_mm) || 0;
  S.grid = Number(cfg.grid) || S.grid;
  S.fuentes = Array.isArray(cfg.sources) && cfg.sources.length
    ? cfg.sources.map((s) => ({ y: Number(s.y), width: Number(s.width),
                                amp: Number(s.amp ?? 1), phase: Number(s.phase ?? 0) }))
    : null;
  if (cfg.coclea) {
    COC.espesor = Number(cfg.coclea.thickness_mm) || COC.espesor;
    COC.c_s = Number(cfg.coclea.c_s) || COC.c_s;
    COC.alpha = Number(cfg.coclea.alpha_np_m) || 0;
    $('cocEsp').value = COC.espesor;
    $('cocCs').value = COC.c_s;
    $('cocAlpha').value = COC.alpha;
    sincronizarCoclea();
  }
  $('mat').value = S.material;
  $('model').value = S.modelo;
  $('freq').value = S.freq;
  $('nEl').value = S.n;
  $('ap').value = S.d;
  $('fx').value = S.fx;
  $('fy').value = S.fy;
  $('res').value = String(S.grid);
  sincronizarEtiquetas();
  if (cfg.modo === 'manual') {
    // una config manual sin fuentes explícitas arranca del arreglo: es lo que espera el que
    // guardó la vista en la pestaña "fuentes" sin haber tocado nada
    if (!S.fuentes) generarDesdeArreglo();
    cambiarModo('manual');
  } else {
    cambiarModo('array');
  }
}

$('btnExportarVista').addEventListener('click', () => window.exportarVista());
$('btnGuardar').addEventListener('click', () => {
  const texto = JSON.stringify(instantanea(), null, 2);
  const url = URL.createObjectURL(new Blob([texto], { type: 'application/json' }));
  const a = document.createElement('a');
  a.href = url;
  a.download = `sim_${S.freq}kHz_${S.n}el.json`;
  a.click();
  URL.revokeObjectURL(url);
  $('estado').textContent = 'configuración descargada';
});

$('btnCargar').addEventListener('click', () => $('archivo').click());
$('archivo').addEventListener('change', (ev) => {
  const f = ev.target.files && ev.target.files[0];
  if (!f) return;
  const lector = new FileReader();
  lector.onload = () => {
    try {
      aplicar(JSON.parse(lector.result));
      $('estado').textContent = `configuración cargada: ${f.name}`;
    } catch (e) {
      $('estado').className = 'pie error';
      $('estado').textContent = e.message;
    }
  };
  lector.readAsText(f);
  ev.target.value = '';                                  // permite recargar el mismo archivo
});

/* se recuerda la última configuración: recargar la página no debería perder el trabajo */
const CLAVE_LOCAL = 'ultrasound-sim/config';
function recordar() {
  try {
    localStorage.setItem(CLAVE_LOCAL, JSON.stringify(instantanea()));
  } catch (e) { /* modo privado o cuota llena: no es un error del simulador */ }
}
function olvidar() {
  try { return JSON.parse(localStorage.getItem(CLAVE_LOCAL) || 'null'); } catch (e) { return null; }
}

window.addEventListener('resize', () => { clearTimeout(window.__tAsas); window.__tAsas = setTimeout(pintarAsas, 150); });

/* ---------------------------------------------------------------- cóclea (corte transversal)
 * Estos tres números son los que faltan medir (D20/D30): se editan acá y el corte se resuelve aparte
 * del mapa, porque es otro modelo (cápsula elástica + fluido + membrana), no una apertura radiando.
 * La resolución cuesta ~2.5 s, así que se dispara con el botón y no en cada movimiento del slider.
 */
const COC = { espesor: 2.0, c_s: 1900.0, alpha: 0.0 };

function sincronizarCoclea() {
  $('cocEspVal').textContent = `${COC.espesor.toFixed(1)} mm`;
  $('cocCsVal').textContent = `${COC.c_s.toFixed(0)} m/s`;
  $('cocAlphaVal').textContent = `${COC.alpha.toFixed(1)} Np/m`;
}

async function resolverCoclea() {
  const est = $('cocEstado');
  est.className = 'pie';
  est.textContent = 'resolviendo el corte…';
  try {
    const r = await fetch('/api/cochlea', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ thickness_mm: COC.espesor, c_s: COC.c_s, alpha_np_m: COC.alpha }),
    });
    const d = await r.json();
    if (!r.ok) throw new Error(d.error || `HTTP ${r.status}`);
    $('cocT').textContent = d.t_interior.toFixed(4);
    $('cocN0').textContent = d.t_interior_n0.toFixed(4);
    $('cocAbs').textContent = `${(100.0 * d.absorbido).toFixed(1)} %`;
    $('cocTiempo').textContent = `${d.meta.seconds.toFixed(1)} s`;
    est.textContent = `${d.gdl} nodos · ${d.meta.seconds} s · t=2 mm sin membrana daría 0.9722`;
    est.scrollIntoView({ block: 'nearest' });
  } catch (e) {
    est.className = 'pie error';
    est.textContent = e.message;
  }
}

const enlazarCoclea = (id, clave) => {
  $(id).addEventListener('input', (e) => {
    COC[clave] = parseFloat(e.target.value);
    sincronizarCoclea();
  });
};
enlazarCoclea('cocEsp', 'espesor');
enlazarCoclea('cocCs', 'c_s');
enlazarCoclea('cocAlpha', 'alpha');
$('btnCoclea').addEventListener('click', resolverCoclea);
sincronizarCoclea();

/* ---------------------------------------------------------------- arranque */
const previa = olvidar();
if (previa) {
  try { aplicar(previa); } catch (e) { /* configuración vieja o corrupta: se ignora */ }
} else {
  sincronizarEtiquetas();
  enviar();
}
setInterval(recordar, 3000);
