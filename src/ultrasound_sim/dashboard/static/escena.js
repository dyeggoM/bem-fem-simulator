/* Escena FEM 2D (Fase U1) ------------------------------------------------------
 *
 * Vista independiente del explorador de campo: acá no hay apertura ni Rayleigh, hay un apilado de
 * capas con material y un solve FEM de verdad (`helmholtz2d.solve`), validado contra la matriz de
 * transferencia. Todo número que se muestra viene del backend; acá no se calcula física.
 *
 * El sello de la derecha no es decorativo (D33): se arma con el error que devolvió la corrida y el
 * criterio de la Fase 1 (2 % en amplitud). Si el error no está, dice que no está.
 */

const ESC = {
  freq: 100,          // kHz
  ppw: 12,
  alto: 20,           // mm
  // caso por defecto con contraste visible (agua/hueso/agua, R = 0.90): con agua/tejido R vale
  // 1.2 % y el mapa sale casi uniforme, que es correcto pero no muestra nada
  capas: [
    { material: 'agua', espesor_mm: 20 },
    { material: 'hueso', espesor_mm: 10 },
    { material: 'agua', espesor_mm: 20 },
  ],
  rango: 40,          // dB que cubre la escala de color, o 'lineal'
  resuelto: null,
  materiales: [],
  ocupado: false,
};

/* El criterio NO vive acá: lo evalúa el servidor (dashboard/validacion.py) y llega en `r.sello`.
   Si hiciera falta un umbral en el cliente sería una señal de que la validación quedó partida. */
function mostrarSello(id, sello, alternativa) {
  const el = esc$(id);
  if (!sello || !sello.texto) {
    el.textContent = alternativa || 'sin sello: esta corrida no reportó criterio.';
    el.className = 'hint sello-sin_referencia';
    return;
  }
  el.textContent = sello.texto;
  el.className = `hint sello-${sello.estado}`;
}

const esc$ = (id) => document.getElementById(id);
const escFmt = (v, d = 4) => (v === null || v === undefined || Number.isNaN(v) ? '—' : Number(v).toFixed(d));

function escPayload() {
  return {
    frecuencia_hz: ESC.freq * 1e3,
    ppw: ESC.ppw,
    alto_m: ESC.alto * 1e-3,
    capas: ESC.capas.map((c) => ({ material: c.material, espesor_m: c.espesor_mm * 1e-3 })),
  };
}

async function escPost(ruta, cuerpo) {
  const r = await fetch(ruta, {
    method: 'POST',
    headers: { 'content-type': 'application/json' },
    body: JSON.stringify(cuerpo),
  });
  const datos = await r.json();
  if (!datos.ok) throw new Error(datos.error || 'el servidor no devolvió ok');
  return datos;
}

/* ------------------------------------------------------------------ capas */

function escPintarCapas() {
  const caja = esc$('capas');
  caja.innerHTML = '';
  ESC.capas.forEach((capa, i) => {
    const fila = document.createElement('div');
    fila.className = 'fila-capa';

    const sel = document.createElement('select');
    for (const m of ESC.materiales) {
      const opt = document.createElement('option');
      opt.value = m.id;
      opt.textContent = `${m.id} · ${Math.round(m.c_m_s)} m/s`;
      if (m.id === capa.material) opt.selected = true;
      sel.appendChild(opt);
    }
    sel.addEventListener('change', () => { capa.material = sel.value; esc$('capaAviso').textContent = ''; });

    const esp = document.createElement('input');
    esp.type = 'number';
    esp.min = '0.2';
    esp.step = '0.5';
    esp.value = String(capa.espesor_mm);
    esp.addEventListener('change', () => {
      const v = parseFloat(esp.value);
      if (!Number.isFinite(v) || v <= 0) { esc$('capaAviso').textContent = 'el espesor tiene que ser > 0'; return; }
      capa.espesor_mm = v;
      esc$('capaAviso').textContent = '';
    });

    const del = document.createElement('button');
    del.className = 'btn mini';
    del.textContent = '×';
    del.title = 'quitar capa';
    del.disabled = ESC.capas.length <= 2;
    del.addEventListener('click', () => {
      if (ESC.capas.length <= 2) return;
      ESC.capas.splice(i, 1);
      escPintarCapas();
    });

    fila.append(sel, esp, del);
    caja.appendChild(fila);
  });
  esc$('btnCapaAdd').disabled = ESC.capas.length >= 12;
}

/* ------------------------------------------------------- estimar y resolver */

async function escEstimar() {
  esc$('escEstado').textContent = 'estimando…';
  try {
    const r = await escPost('/api/scene/cost', escPayload());
    esc$('escNodos').textContent = r.grilla.n_nodos;
    esc$('escElems').textContent = r.grilla.n_elementos;
    esc$('escPaso').textContent = `${escFmt(r.grilla.paso_y_mm, 2)} mm`;
    esc$('escEstimado').textContent = `${escFmt(r.estimado_solve_s, 2)} s`;
    esc$('escEstado').textContent =
      `malla ${r.grilla.nx}×${r.grilla.ny} · ${escFmt(r.grilla.longitud_mm, 1)} mm de largo · interfaces en ` +
      `${r.grilla.interfaces_mm.map((v) => escFmt(v, 1)).join(', ')} mm`;
    return r;
  } catch (e) {
    esc$('escEstado').textContent = `no se pudo estimar: ${e.message}`;
    return null;
  }
}

async function escResolver() {
  if (ESC.ocupado) return;
  ESC.ocupado = true;
  const boton = esc$('btnEscena');
  boton.disabled = true;
  esc$('escEstado').textContent = 'resolviendo…';
  try {
    const r = await escPost('/api/scene/solve', escPayload());
    ESC.resuelto = r;
    escPintarResultado(r);
    escDibujar(r);
    window.vistaGuardar('escena', r);   // volver a la pestaña no debe costar otro solve
    esc$('escEstado').textContent =
      `resuelto: ${r.grilla.n_elementos} elementos, ${r.grilla.n_nodos} nodos · ` +
      `malla ${escFmt(r.tiempos.malla_s, 3)} s + solve ${escFmt(r.tiempos.solve_s, 3)} s`;
  } catch (e) {
    esc$('escEstado').textContent = `no se pudo resolver: ${e.message}`;
  } finally {
    ESC.ocupado = false;
    boton.disabled = false;
  }
}

function escPintarResultado(r) {
  const m = r.metricas;
  const hay = m.R !== undefined;
  esc$('escR').textContent = hay ? escFmt(m.R) : '—';
  esc$('escRan').textContent = hay ? escFmt(m.R_analitico) : '—';
  esc$('escT').textContent = hay ? escFmt(m.T) : '—';
  esc$('escTan').textContent = hay ? escFmt(m.T_analitico) : '—';
  esc$('escEnergia').textContent = hay ? escFmt(m.energia) : '—';
  esc$('escTiempo').textContent = `${escFmt(r.tiempos.solve_s, 3)} s`;

  esc$('escErrR').textContent = hay ? `abs ${escFmt(m.error_R_abs)}` : '—';
  esc$('escErrT').textContent = hay ? `abs ${escFmt(m.error_T_abs)}` : '—';
  // El sello lo evalúa el servidor con su criterio y su fuente (U6): acá no hay umbrales.
  mostrarSello('escSello', r.sello, m.aviso || 'la corrida no reportó R ni T.');
}

/* ------------------------------------------------------------------ dibujo */

const ESC_CONFIG = { displayModeBar: false, responsive: true, displaylogo: false };

function escDibujar(r) {
  const interfases = r.grilla.interfaces_mm.map((x) => ({
    type: 'line', x0: x, x1: x, yref: 'paper', y0: 0, y1: 1,
    line: { color: '#58a6ff', width: 1, dash: 'dot' },
  }));

  Plotly.react('mapaEscena', [{
    type: 'heatmap',
    x: r.grilla.x_mm,
    y: r.grilla.y_mm,
    z: ESC.rango === 'lineal' ? r.campo.amplitud : r.campo.db,
    colorscale: 'Inferno',
    zmin: ESC.rango === 'lineal' ? 0 : -ESC.rango,
    zmax: ESC.rango === 'lineal' ? Math.max(...r.campo.amplitud.flat()) : 0,
    colorbar: { title: { text: ESC.rango === 'lineal' ? '|p|' : 'dB rel. pico', side: 'top' },
                thickness: 11, len: 0.8,
                tickfont: { size: 9, color: '#8b95a1' }, outlinewidth: 0 },
    hovertemplate: 'x %{x:.1f} mm<br>y %{y:.1f} mm<br>%{z:.1f} dB<extra></extra>',
  }], {
    margin: { l: 52, r: 8, t: 10, b: 40 },
    paper_bgcolor: 'rgba(0,0,0,0)',
    plot_bgcolor: 'rgba(0,0,0,0)',
    font: { color: '#8b95a1', size: 10 },
    xaxis: { title: { text: 'x (mm)', font: { size: 10 } }, gridcolor: '#21262d', zeroline: false },
    yaxis: { title: { text: 'y (mm)', font: { size: 10 } }, gridcolor: '#21262d', zeroline: false },
    shapes: interfases,
    annotations: [{
      x: 0, y: 1.02, xref: 'x', yref: 'paper', showarrow: false,
      text: 'onda plana →', font: { size: 9, color: '#8b95a1' }, xanchor: 'left',
    }],
  }, ESC_CONFIG);

  const c = r.cortes;
  Plotly.react('axialEscena', [{
    type: 'scatter', mode: 'lines', x: c.eje_x_mm, y: c.axial,
    line: { color: '#58a6ff', width: 1.4 },
    hovertemplate: 'x %{x:.1f} mm<br>|p| %{y:.3f}<extra></extra>',
  }], {
    ...ESC_CONFIG,
    margin: { l: 46, r: 8, t: 24, b: 34 },
    title: { text: 'corte axial (|p| en y medio)', font: { size: 10, color: '#8b95a1' } },
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
    font: { color: '#8b95a1', size: 9 },
    xaxis: { title: { text: 'x (mm)' }, gridcolor: '#21262d' },
    yaxis: { gridcolor: '#21262d' },
    shapes: interfases.map((s) => ({ ...s, yref: 'paper' })),
  });

  Plotly.react('lateralEscena', [{
    type: 'scatter', mode: 'lines', x: c.lateral_en_x_pico, y: c.eje_y_mm,
    line: { color: '#f0883e', width: 1.4 },
    hovertemplate: 'y %{y:.1f} mm<br>|p| %{x:.3f}<extra></extra>',
  }], {
    ...ESC_CONFIG,
    margin: { l: 46, r: 8, t: 24, b: 34 },
    title: { text: `corte lateral en x = ${escFmt(c.x_pico_mm, 1)} mm`, font: { size: 10, color: '#8b95a1' } },
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
    font: { color: '#8b95a1', size: 9 },
    xaxis: { title: { text: '|p|' }, gridcolor: '#21262d' },
    yaxis: { title: { text: 'y (mm)' }, gridcolor: '#21262d' },
  });
}

/* ------------------------------------------------------------- vista activa */

/* Tres vistas, un solo conmutador: panel + main por nombre. Agregar una vista es agregarla acá. */
const VISTAS = {
  campo: { panel: 'vistaCampo', main: 'mainCampo' },
  escena: { panel: 'vistaEscena', main: 'vistaEscenaMain' },
  dispersor: { panel: 'vistaDispersor', main: 'vistaDispersorMain' },
  cochlea: { panel: 'vistaCochlea', main: 'vistaCochleaMain' },
  barridos: { panel: 'vistaBarridos', main: 'vistaBarridosMain' },
  volumen: { panel: 'vistaVolumen', main: 'vistaVolumenMain' },
  esfera: { panel: 'vistaEsfera', main: 'vistaEsferaMain' },
  criterios: { panel: null, main: 'vistaCriteriosMain' },
};

function escCambiarVista(vista) {
  document.querySelectorAll('.vista-btn').forEach((b) => b.classList.toggle('activo', b.dataset.vista === vista));
  for (const [nombre, ref] of Object.entries(VISTAS)) {
    const panel = esc$(ref.panel);
    const main = esc$(ref.main);
    if (panel) panel.hidden = nombre !== vista;
    if (main) main.hidden = nombre !== vista;
  }
  if (vista === 'dispersor' && typeof window.disIniciar === 'function') window.disIniciar();
  if (vista === 'cochlea' && typeof window.cocfIniciar === 'function') window.cocfIniciar();
  if (vista === 'barridos' && typeof window.barIniciar === 'function') window.barIniciar();
  if (vista === 'volumen' && typeof window.volIniciar === 'function') window.volIniciar();
  if (vista === 'esfera' && typeof window.esfIniciar === 'function') window.esfIniciar();
  if (vista === 'criterios' && typeof window.critIniciar === 'function') window.critIniciar();
  // el aviso naranja y la línea de estado son del explorador de campo ("32 fuentes · baffle…",
  // habla del foco de la apertura): en la vista de escena confunden. Se ocultan por CSS
  // (`body.vista-escena`), no por `hidden`: `app.js` los re-muestra cuando llega su respuesta.
  document.body.classList.toggle('vista-escena', vista === 'escena');
  document.body.classList.toggle('vista-dispersor', vista === 'dispersor');
  document.body.classList.toggle('vista-cochlea', vista === 'cochlea');
  document.body.classList.toggle('vista-barridos', vista === 'barridos');
  document.body.classList.toggle('vista-volumen', vista === 'volumen');
  document.body.classList.toggle('vista-esfera', vista === 'esfera');
  if (vista === 'escena') {
    // resultado de la sesión anterior: se repinta antes de estimar, así volver a la pestaña no
    // obliga a resolver de nuevo
    const previo = ESC.resuelto ? null : window.vistaCargar('escena');
    if (previo) {
      try {
        ESC.resuelto = previo.resultado;
        escPintarResultado(previo.resultado);
        escDibujar(previo.resultado);
        esc$('escEstado').textContent =
          'restaurado de ' + new Date(previo.fecha).toLocaleString('es');
      } catch (e) {
        console.warn('no se pudo restaurar la escena:', e && e.message);
        window.vistaBorrar('escena');
        ESC.resuelto = null;
      }
    }
    if (!ESC.materiales.length) escCargarMateriales().then(escEstimar);
    else if (!ESC.resuelto) escEstimar();
    setTimeout(() => window.dispatchEvent(new Event('resize')), 30);
  } else if (mainCampo) {
    setTimeout(() => window.dispatchEvent(new Event('resize')), 30);
  }
}

async function escCargarMateriales() {
  try {
    const r = await fetch('/api/materials');
    const datos = await r.json();
    ESC.materiales = datos.scene_materials || [];
  } catch (e) {
    ESC.materiales = [
      { id: 'agua', c_m_s: 1500 }, { id: 'tejido', c_m_s: 1540 }, { id: 'hueso', c_m_s: 3500 },
    ];
  }
  escPintarCapas();
}

/* --------------------------------------------------------------- enlaces */

document.querySelectorAll('.vista-btn').forEach((b) => {
  b.addEventListener('click', () => escCambiarVista(b.dataset.vista));
});

esc$('escFreq').addEventListener('input', (e) => {
  ESC.freq = parseInt(e.target.value, 10);
  esc$('escFreqVal').textContent = `${ESC.freq} kHz`;
});
esc$('escPpw').addEventListener('input', (e) => {
  ESC.ppw = parseInt(e.target.value, 10);
  esc$('escPpwVal').textContent = String(ESC.ppw);
});
esc$('escAlto').addEventListener('input', (e) => {
  ESC.alto = parseInt(e.target.value, 10);
  esc$('escAltoVal').textContent = `${ESC.alto} mm`;
});
esc$('btnCapaAdd').addEventListener('click', () => {
  if (ESC.capas.length >= 12) return;
  // la capa nueva hereda el material de la última: casi siempre se agrega para alargar el medio
  ESC.capas.push({ material: ESC.capas[ESC.capas.length - 1].material, espesor_mm: 10 });
  escPintarCapas();
});
esc$('escRango').addEventListener('change', (e) => {
  ESC.rango = e.target.value === 'lineal' ? 'lineal' : parseInt(e.target.value, 10);
  if (ESC.resuelto) escDibujar(ESC.resuelto);
});
esc$('btnEstimar').addEventListener('click', escEstimar);
esc$('btnEscena').addEventListener('click', escResolver);

/* --------------------------------------------------------------- arranque */
escCargarMateriales();
