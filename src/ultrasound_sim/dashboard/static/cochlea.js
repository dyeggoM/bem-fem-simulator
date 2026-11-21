/* Cóclea de trabajo (Fase U3) ----------------------------------------------------------------
 *
 * El corte transversal de la Fase 5 completo: cápsula elástica, perilinfa, membrana basilar y
 * ventana redonda. El mapa es el campo interior interpolado (no un esquema): sale de `sol.p` del
 * solver acoplado, vía `dashboard/cochlea_field.py`.
 *
 * Q no se muestra: medirlo pide un barrido de frecuencia (≈ 9 solves de 5 s) y eso es U5. Acá se
 * dice en vez de inventar un número.
 */

const COCF = {
  esp: 2.0,
  cs: 1900,
  alpha: 0,
  membrana: true,
  ventana: true,
  resuelto: null,
  ocupado: false,
};

/* Las referencias de validación viven en el servidor (`dashboard/validacion.py`): acá no hay
   ningún número de D30 escrito. `mostrarSello` viene de escena.js. */

const cocf$ = (id) => document.getElementById(id);
const cocfFmt = (v, d = 4) => (v === null || v === undefined || Number.isNaN(v) ? '—' : Number(v).toFixed(d));

function cocfPayload() {
  return {
    thickness_mm: COCF.esp,
    c_s: COCF.cs,
    alpha_np_m: COCF.alpha,
    membrana: COCF.membrana,
    ventana: COCF.ventana,
  };
}

async function cocfPost(ruta, cuerpo) {
  const r = await fetch(ruta, {
    method: 'POST',
    headers: { 'content-type': 'application/json' },
    body: JSON.stringify(cuerpo),
  });
  const datos = await r.json();
  if (!datos.ok) throw new Error(datos.error || 'el servidor no devolvió ok');
  return datos;
}

async function cocfEstimar() {
  cocf$('cocFieldEstado').textContent = 'estimando…';
  try {
    const r = await cocfPost('/api/cochlea/field/cost', cocfPayload());
    cocf$('cocFieldNodos').textContent = r.malla.n_nodos;
    cocf$('cocFieldEstimado').textContent = `${cocfFmt(r.estimado_solve_s, 1)} s`;
    cocf$('cocFieldEstado').textContent =
      `${r.malla.n_elementos} elementos · ${r.malla.n_nodos} nodos · a = ${cocfFmt(r.malla.radio_mm, 1)} mm`;
    return r;
  } catch (e) {
    cocf$('cocFieldEstado').textContent = `no se pudo estimar: ${e.message}`;
    return null;
  }
}

async function cocfResolver() {
  if (COCF.ocupado) return;
  COCF.ocupado = true;
  const boton = cocf$('btnCocField');
  boton.disabled = true;
  cocf$('cocFieldEstado').textContent = 'resolviendo el corte (unos segundos)…';
  try {
    const r = await cocfPost('/api/cochlea/field', cocfPayload());
    COCF.resuelto = r;
    cocfPintar(r);
    cocfDibujar(r);
    cocf$('cocFieldEstado').textContent =
      `resuelto: ${r.malla.n_elementos} elementos · malla ${cocfFmt(r.tiempos.malla_s, 2)} s + ` +
      `solve ${cocfFmt(r.tiempos.solve_s, 2)} s`;
  } catch (e) {
    cocf$('cocFieldEstado').textContent = `no se pudo resolver: ${e.message}`;
  } finally {
    COCF.ocupado = false;
    boton.disabled = false;
  }
}

function cocfPintar(r) {
  const m = r.metricas;
  cocf$('cocFieldT').textContent = cocfFmt(m.t_interior, 5);
  cocf$('cocFieldN0').textContent = cocfFmt(m.t_interior_n0, 5);
  cocf$('cocFieldAbs').textContent = `${cocfFmt(100 * m.absorcion_membrana, 1)} %`;
  cocf$('cocFieldVentana').textContent = `${cocfFmt(100 * m.ventana_real, 2)} %`;
  cocf$('cocFieldBalance').textContent = `${cocfFmt(100 * m.balance_energia, 2)} %`;
  cocf$('cocFieldTiempo').textContent = `${cocfFmt(r.tiempos.solve_s, 2)} s`;

  // El sello viene del servidor: incluye la comparación contra D30 (sólo cuando los parámetros son
  // los del caso base) y el cierre del balance. Acá no hay umbrales ni referencias (U6).
  mostrarSello('cocFieldSello', r.sello,
    'esta corrida no reportó criterio. Q (ancho de banda) sigue pendiente: necesita una definición.');
}

/* ------------------------------------------------------------------ dibujo */

const COCF_CONFIG = { displayModeBar: false, responsive: true, displaylogo: false };

function cocfDibujar(r) {
  const a = r.mapa.radio_mm;
  const t = r.mapa.cascara_mm;
  const formas = [
    { type: 'circle', x0: -a, x1: a, y0: -a, y1: a,
      line: { color: '#58a6ff', width: 1.2 }, fillcolor: 'rgba(0,0,0,0)' },
    { type: 'circle', x0: -(a + t), x1: a + t, y0: -(a + t), y1: a + t,
      line: { color: '#f0883e', width: 1.2, dash: 'dot' }, fillcolor: 'rgba(0,0,0,0)' },
  ];
  if (r.mapa.membrana_angulo_deg !== null) {
    const rad = (r.mapa.membrana_angulo_deg * Math.PI) / 180;
    formas.push({
      type: 'line',
      x0: -a * Math.cos(rad), y0: -a * Math.sin(rad),
      x1: a * Math.cos(rad), y1: a * Math.sin(rad),
      line: { color: '#3fb950', width: 2 },
    });
  }

  Plotly.react('mapaCochlea', [{
    type: 'heatmap',
    x: r.mapa.x_mm, y: r.mapa.y_mm, z: r.mapa.db,
    colorscale: 'Inferno', zmin: -30, zmax: 0,
    colorbar: { title: { text: 'dB rel. pico', side: 'top' }, thickness: 11, len: 0.8,
                tickfont: { size: 9, color: '#8b95a1' }, outlinewidth: 0 },
    hovertemplate: 'x %{x:.1f} mm<br>y %{y:.1f} mm<br>%{z:.1f} dB<extra></extra>',
  }], {
    margin: { l: 52, r: 8, t: 10, b: 40 },
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
    font: { color: '#8b95a1', size: 10 },
    xaxis: { title: { text: 'x (mm)' }, gridcolor: '#21262d', zeroline: false },
    yaxis: { title: { text: 'y (mm)' }, gridcolor: '#21262d', zeroline: false,
             scaleanchor: 'x', scaleratio: 1 },
    shapes: formas,
    annotations: [{
      x: 0.01, y: 0.99, xref: 'paper', yref: 'paper', showarrow: false, xanchor: 'left',
      text: 'azul: pared interior · naranja: pared exterior · verde: membrana',
      font: { size: 9, color: '#8b95a1' },
    }],
  }, COCF_CONFIG);

  // corte radial: |p| sobre el diámetro horizontal, que muestra la onda estacionaria de la cavidad
  const fila = Math.floor(r.mapa.y_mm.length / 2);
  const amplitud = r.mapa.db[fila].map((v) => (v === null ? null : Math.pow(10, v / 20)));
  Plotly.react('corteCochlea', [{
    type: 'scatter', mode: 'lines', x: r.mapa.x_mm, y: amplitud,
    line: { color: '#58a6ff', width: 1.5 }, connectgaps: false,
    hovertemplate: 'x %{x:.1f} mm<br>|p|/pico %{y:.3f}<extra></extra>',
  }], {
    ...COCF_CONFIG,
    margin: { l: 52, r: 8, t: 24, b: 34 },
    title: { text: 'corte por el centro: |p|/pico', font: { size: 10, color: '#8b95a1' } },
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
    font: { color: '#8b95a1', size: 9 },
    xaxis: { title: { text: 'x (mm)' }, gridcolor: '#21262d' },
    yaxis: { title: { text: '|p|/pico' }, gridcolor: '#21262d' },
  });
}

/* --------------------------------------------------------------- enlaces */

function cocfIniciar() {
  if (!COCF.resuelto) cocfEstimar();
  setTimeout(() => window.dispatchEvent(new Event('resize')), 30);
}
window.cocfIniciar = cocfIniciar;

cocf$('cocFieldEsp').addEventListener('input', (e) => {
  COCF.esp = parseFloat(e.target.value);
  cocf$('cocFieldEspVal').textContent = `${cocfFmt(COCF.esp, 1)} mm`;
});
cocf$('cocFieldCs').addEventListener('input', (e) => {
  COCF.cs = parseInt(e.target.value, 10);
  cocf$('cocFieldCsVal').textContent = `${COCF.cs} m/s`;
});
cocf$('cocFieldAlpha').addEventListener('input', (e) => {
  COCF.alpha = parseInt(e.target.value, 10);
  cocf$('cocFieldAlphaVal').textContent = `${COCF.alpha} Np/m`;
});
cocf$('cocFieldMemb').addEventListener('change', (e) => { COCF.membrana = e.target.checked; });
cocf$('cocFieldVent').addEventListener('change', (e) => { COCF.ventana = e.target.checked; });
cocf$('btnCocFieldEstimar').addEventListener('click', cocfEstimar);
cocf$('btnCocField').addEventListener('click', cocfResolver);
