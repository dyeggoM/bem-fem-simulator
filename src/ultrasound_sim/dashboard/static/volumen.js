/* Escena 3D por rebanadas (Fase U4) ---------------------------------------------------------
 *
 * Dos modos:
 *   - puntual: monopolo en una caja de agua, ABC en las seis caras. El mapa muestra una onda
 *     esférica y la métrica es el error contra el Green 3D (comparado en un anillo).
 *   - plano: capas apiladas con onda plana normal, ABC sólo en las caras de x. La métrica es R y T
 *     contra la matriz de transferencia.
 *
 * El caso plano es físicamente 1D (uniforme en y/z): eso no es un defecto del render, es la física
 * de la incidencia normal. El modo puntual es el que muestra algo realmente 3D.
 */

const VOL = {
  fuente: 'puntual',
  freq: 100,     // kHz
  ppw: 10,
  alto: 30,      // mm ≈ 2 λ a 100 kHz: con 1.5 λ no entra ni un anillo completo
  capas: [
    { material: 'agua', espesor_mm: 22 },
    { material: 'hueso', espesor_mm: 8 },
    { material: 'agua', espesor_mm: 22 },
  ],
  ka: 1.0,
  matEsfera: 'hueso',
  resuelto: null,
  ocupado: false,
};

/* Los criterios (Green al 10 %, R/T al 2 %) los evalúa el servidor: acá no hay umbrales (U6). */

const vol$ = (id) => document.getElementById(id);
const volFmt = (v, d = 2) => (v === null || v === undefined || Number.isNaN(v) ? '—' : Number(v).toFixed(d));

function volPayload() {
  const homogenea = VOL.fuente === 'puntual' || VOL.fuente === 'esfera';
  const capas = homogenea
    ? [{ material: 'agua', espesor_m: VOL.alto * 1e-3 }]
    : VOL.capas.map((c) => ({ material: c.material, espesor_m: c.espesor_mm * 1e-3 }));
  const cuerpo = {
    frecuencia_hz: VOL.freq * 1e3,
    ppw: VOL.ppw,
    alto_m: VOL.alto * 1e-3,
    fuente: VOL.fuente,
    capas,
  };
  if (VOL.fuente === 'esfera') {
    cuerpo.ka = VOL.ka;
    cuerpo.material_esfera = VOL.matEsfera;
  }
  return cuerpo;
}

async function volPost(ruta, cuerpo) {
  const r = await fetch(ruta, {
    method: 'POST',
    headers: { 'content-type': 'application/json' },
    body: JSON.stringify(cuerpo),
  });
  const datos = await r.json();
  if (!datos.ok) throw new Error(datos.error || 'el servidor no devolvió ok');
  return datos;
}

async function volEstimar() {
  vol$('volEstado').textContent = 'estimando…';
  try {
    const r = await volPost('/api/volume3d/cost', volPayload());
    vol$('volNodos').textContent = r.malla.n_nodos;
    vol$('volElems').textContent = r.malla.n_elementos;
    vol$('volEstimado').textContent = `${volFmt(r.estimado_solve_s, 1)} s`;
    vol$('volEstado').textContent =
      `${r.malla.nx}×${r.malla.ny}×${r.malla.nz} nodos · ${volFmt(r.malla.longitud_mm, 1)} mm de largo · ` +
      `techo conocido: ${r.techo_conocido}`;
    return r;
  } catch (e) {
    vol$('volEstado').textContent = `no se pudo estimar: ${e.message}`;
    return null;
  }
}

async function volResolver() {
  if (VOL.ocupado) return;
  VOL.ocupado = true;
  const boton = vol$('btnVolumen');
  boton.disabled = true;
  vol$('volEstado').textContent = 'resolviendo 3D… (puede tardar decenas de segundos)';
  try {
    const r = await volPost('/api/volume3d/solve', volPayload());
    VOL.resuelto = r;
    volPintar(r);
    volDibujar(r);
    vol$('volEstado').textContent =
      `resuelto: ${r.malla.n_nodos} nodos, ${r.malla.n_elementos} tets · ` +
      `malla ${volFmt(r.tiempos.malla_s, 2)} s + solve ${volFmt(r.tiempos.solve_s, 2)} s`;
  } catch (e) {
    vol$('volEstado').textContent = `no se pudo resolver: ${e.message}`;
  } finally {
    VOL.ocupado = false;
    boton.disabled = false;
  }
}

function volPintar(r) {
  const m = r.metricas;
  vol$('volTiempo').textContent = `${volFmt(r.tiempos.solve_s, 2)} s`;
  const esPuntual = r.escena.fuente === 'puntual';

  const esEsfera = r.escena.fuente === 'esfera';
  vol$('volGreen').textContent = esPuntual ? `${volFmt(m.error_green_pct, 1)} %` : '—';
  vol$('volAmp').textContent = esPuntual ? volFmt(m.amplitud_green, 3) : '—';
  vol$('volSerie').textContent = esEsfera ? `${volFmt(m.error_serie_pct, 2)} %` : '—';
  vol$('volRadio').textContent = esEsfera ? `${volFmt(m.radio_mm, 2)} mm` : '—';
  vol$('volR').textContent = esPuntual ? '—' : volFmt(m.R, 4);
  vol$('volRan').textContent = esPuntual ? '—' : volFmt(m.R_analitico, 4);
  vol$('volErrR').textContent = esPuntual ? '—' : `abs ${volFmt(m.error_R_abs, 4)}`;
  vol$('volErrT').textContent = esPuntual ? '—' : `abs ${volFmt(m.error_T_abs, 4)}`;

  mostrarSello('volSello', r.sello, 'esta corrida no reportó criterio.');
}

/* ------------------------------------------------------------------ dibujo */

const VOL_CONFIG = { displayModeBar: false, responsive: true, displaylogo: false };

function volDibujar(r) {
  const s = r.rebanadas;

  // En la esfera la onda viaja en +z: el plano que muestra los lóbulos es el x-z, y el x-y (a z=0,
  // perpendicular a la propagación) sale casi simétrico — una mancha oscura central. Por eso en ese
  // modo se invierte cuál es el mapa principal.
  const esEsfera = r.escena.fuente === 'esfera';
  const principal = esEsfera
    ? { x: s.xz.x_mm, y: s.xz.z_mm, z: s.xz.db, titulo: `campo dispersado, plano de propagación (${s.xz.donde})`,
        ejeY: 'z (mm)', colorbar: true }
    : { x: s.xy.x_mm, y: s.xy.y_mm, z: s.xy.db, titulo: `rebanada x-y (${s.xy.donde})`,
        ejeY: 'y (mm)', colorbar: true };
  const secundaria = esEsfera
    ? { x: s.xy.x_mm, y: s.xy.y_mm, z: s.xy.db, titulo: `rebanada x-y — perpendicular a la propagación (${s.xy.donde})`,
        ejeY: 'y (mm)' }
    : { x: s.xz.x_mm, y: s.xz.z_mm, z: s.xz.db, titulo: `rebanada x-z (${s.xz.donde})`,
        ejeY: 'z (mm)' };

  const interfases = s.interfaces_mm.map((x) => ({
    type: 'line', x0: x, x1: x, yref: 'paper', y0: 0, y1: 1,
    line: { color: '#58a6ff', width: 1, dash: 'dot' },
  }));

  Plotly.react('mapaVolumen', [{
    type: 'heatmap', x: principal.x, y: principal.y, z: principal.z,
    colorscale: 'Inferno', zmin: -30, zmax: 0,
    colorbar: { title: { text: 'dB rel. pico', side: 'top' }, thickness: 11, len: 0.8,
                tickfont: { size: 9, color: '#8b95a1' }, outlinewidth: 0 },
    hovertemplate: 'x %{x:.1f} mm<br>y %{y:.1f} mm<br>%{z:.1f} dB<extra></extra>',
  }], {
    margin: { l: 52, r: 8, t: 26, b: 40 },
    title: { text: principal.titulo, font: { size: 10, color: '#8b95a1' } },
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
    font: { color: '#8b95a1', size: 10 },
    // aspecto 1:1: sin esto el div ancho estira 30 mm sobre 740 px y una onda esférica se ve
    // como un lóbulo horizontal. El campo ES simétrico (el Green cierra al 7 %); era el render.
    xaxis: { title: { text: 'x (mm)' }, gridcolor: '#21262d', zeroline: false,
             scaleanchor: 'y', scaleratio: 1, constrain: 'domain' },
    yaxis: { title: { text: principal.ejeY }, gridcolor: '#21262d', zeroline: false },
    shapes: interfases,
  }, VOL_CONFIG);

  Plotly.react('corteVolumen', [{
    type: 'heatmap', x: secundaria.x, y: secundaria.y, z: secundaria.z,
    colorscale: 'Inferno', zmin: -30, zmax: 0, showscale: false,
    hovertemplate: 'x %{x:.1f} mm<br>y %{y:.1f} mm<br>%{z:.1f} dB<extra></extra>',
  }], {
    ...VOL_CONFIG,
    margin: { l: 52, r: 8, t: 24, b: 34 },
    title: { text: secundaria.titulo, font: { size: 10, color: '#8b95a1' } },
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
    font: { color: '#8b95a1', size: 9 },
    xaxis: { title: { text: 'x (mm)' }, gridcolor: '#21262d',
             scaleanchor: 'y', scaleratio: 1, constrain: 'domain' },
    yaxis: { title: { text: secundaria.ejeY }, gridcolor: '#21262d' },
    shapes: interfases.map((f) => ({ ...f, yref: 'paper' })),
  });
}

/* --------------------------------------------------------------- capas (modo plano) */

function volPintarCapas() {
  const cont = vol$('volCapas');
  cont.innerHTML = '';
  VOL.capas.forEach((capa) => {
    const fila = document.createElement('div');
    fila.className = 'fila-capa';
    const sel = document.createElement('select');
    for (const m of ['agua', 'tejido', 'hueso']) {
      const opt = document.createElement('option');
      opt.value = m;
      opt.textContent = m;
      if (m === capa.material) opt.selected = true;
      sel.appendChild(opt);
    }
    sel.addEventListener('change', () => { capa.material = sel.value; });
    const esp = document.createElement('input');
    esp.type = 'number';
    esp.min = '0.5';
    esp.step = '0.5';
    esp.value = String(capa.espesor_mm);
    esp.addEventListener('change', () => {
      const v = parseFloat(esp.value);
      if (Number.isFinite(v) && v > 0) capa.espesor_mm = v;
    });
    const hueco = document.createElement('span');   // sin botón de borrar: son 3 fijas
    fila.append(sel, esp, hueco);
    cont.appendChild(fila);
  });
}

/* --------------------------------------------------------------- enlaces */

function volModo() {
  vol$('volCapasBox').hidden = VOL.fuente !== 'plano';
  vol$('volEsferaBox').hidden = VOL.fuente !== 'esfera';
}

function volIniciar() {
  volPintarCapas();
  volModo();
  if (!VOL.resuelto) volEstimar();
  setTimeout(() => window.dispatchEvent(new Event('resize')), 30);
}
window.volIniciar = volIniciar;

vol$('volFuente').addEventListener('change', (e) => {
  VOL.fuente = e.target.value;
  volModo();
  volEstimar();
});
vol$('volFreq').addEventListener('input', (e) => {
  VOL.freq = parseInt(e.target.value, 10);
  vol$('volFreqVal').textContent = `${VOL.freq} kHz`;
});
vol$('volPpw').addEventListener('input', (e) => {
  VOL.ppw = parseInt(e.target.value, 10);
  vol$('volPpwVal').textContent = String(VOL.ppw);
});
vol$('volAlto').addEventListener('input', (e) => {
  VOL.alto = parseInt(e.target.value, 10);
  vol$('volAltoVal').textContent = `${VOL.alto} mm`;
});
vol$('volKa').addEventListener('input', (e) => {
  VOL.ka = parseFloat(e.target.value);
  vol$('volKaVal').textContent = volFmt(VOL.ka, 1);
});
vol$('volMatEsfera').addEventListener('change', (e) => { VOL.matEsfera = e.target.value; });
vol$('btnVolEstimar').addEventListener('click', volEstimar);
vol$('btnVolumen').addEventListener('click', volResolver);
