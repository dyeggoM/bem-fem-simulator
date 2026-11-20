/* Dashboard BEM/FEM — cliente. Sin dependencias salvo Plotly, que se sirve local.
 *
 * Todo el estado vive en S; cada control llama a `pedir()`, que debouncea y aborta el pedido
 * anterior. El servidor devuelve el campo ya en dB y las métricas: acá no se calcula física.
 */

const S = {
  material: 'agua',
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

function pedir() {
  clearTimeout(temporizador);
  temporizador = setTimeout(enviar, 150);
}

async function enviar() {
  if (control) control.abort();
  control = new AbortController();
  const cuerpo = {
    material: S.material,
    frequency_khz: S.freq,
    grid: { x0: 0.004, x1: Math.max(0.06, S.fx / 1000 * 1.6), y0: -Math.max(0.02, S.d / 1000 * 1.1),
            y1: Math.max(0.02, S.d / 1000 * 1.1), nx: S.grid, ny: S.grid },
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
      `${datos.meta.n_sources} fuentes · ${datos.meta.nx}×${datos.meta.ny} · ` +
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
    { k: 'amp', v: String(src.amp), paso: '0.05' },
    { k: 'phase', v: String(Number(src.phase.toFixed(3))), paso: '0.1' },
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

sincronizarEtiquetas();
enviar();
