/* Dispersor cilíndrico rígido: BEM, FEM y serie, lado a lado (Fase U2) ------------------------
 *
 * La pestaña que justifica tener dos métodos: el mismo cilindro resuelto por BEM (exterior), FEM
 * (corona en campo dispersado) y la serie de Bessel. La serie no es un adorno: es lo que dice cuál
 * de los dos se equivoca cuando no coinciden.
 *
 * La UI no calcula física: los tres campos y las tres discrepancias vienen de
 * `dashboard/scatterer.py`, que replica la receta ya validada de `scripts/cross_validate.py`.
 */

const DIS = {
  ka: 1.0,
  ppw: 20,
  rOut: 3.0,
  abc: '2',
  ecuacion: 'standard',
  resuelto: null,
  ocupado: false,
};

const CRITERIO_DISP = 2.0;   // Fase 4: discrepancia < 2 % en campo

const dis$ = (id) => document.getElementById(id);
const disFmt = (v, d = 2) => (v === null || v === undefined || Number.isNaN(v) ? '—' : Number(v).toFixed(d));

function disPayload() {
  return {
    ka: DIS.ka,
    ppw: DIS.ppw,
    r_out: DIS.rOut,
    abc_order: parseInt(DIS.abc, 10),
    equation: DIS.ecuacion,
  };
}

async function disPost(ruta, cuerpo) {
  const r = await fetch(ruta, {
    method: 'POST',
    headers: { 'content-type': 'application/json' },
    body: JSON.stringify(cuerpo),
  });
  const datos = await r.json();
  if (!datos.ok) throw new Error(datos.error || 'el servidor no devolvió ok');
  return datos;
}

async function disEstimar() {
  dis$('disEstado').textContent = 'estimando…';
  try {
    const r = await disPost('/api/scatterer/cost', disPayload());
    dis$('disRadio').textContent = `${disFmt(r.geo.radio_mm)} mm`;
    dis$('disElems').textContent = r.malla.n_elementos;
    dis$('disPaneles').textContent = r.paneles_bem;
    dis$('disEstimado').textContent = `${disFmt(r.estimado_solve_s, 1)} s`;
    dis$('disEstado').textContent =
      `a = ${disFmt(r.geo.radio_mm)} mm · R = ${disFmt(r.geo.R_mm)} mm · λ = ${disFmt(r.geo.longitud_onda_mm)} mm`;
    return r;
  } catch (e) {
    dis$('disEstado').textContent = `no se pudo estimar: ${e.message}`;
    return null;
  }
}

async function disResolver() {
  if (DIS.ocupado) return;
  DIS.ocupado = true;
  const boton = dis$('btnDisolver');
  boton.disabled = true;
  dis$('disEstado').textContent = 'resolviendo los tres…';
  try {
    const r = await disPost('/api/scatterer/solve', disPayload());
    DIS.resuelto = r;
    disPintarMetricas(r);
    disDibujar(r);
    dis$('disEstado').textContent =
      `ka = ${disFmt(r.geo.ka, 1)} · ${r.malla.n_elementos} elementos FEM · ` +
      `${r.malla.n_paneles_bem} paneles BEM · total ${disFmt(r.tiempos.total_s, 2)} s`;
  } catch (e) {
    dis$('disEstado').textContent = `no se pudo resolver: ${e.message}`;
  } finally {
    DIS.ocupado = false;
    boton.disabled = false;
  }
}

function disPintarMetricas(r) {
  const m = r.metricas;
  dis$('disFemBem').textContent = `${disFmt(m.fem_vs_bem_pct)} %`;
  dis$('disFemSerie').textContent = `${disFmt(m.fem_vs_serie_pct)} %`;
  dis$('disBemSerie').textContent = `${disFmt(m.bem_vs_serie_pct)} %`;
  dis$('disCond').textContent = disFmt(m.condicion, 2);
  dis$('disTFem').textContent = `${disFmt(r.tiempos.fem_s * 1e3, 0)} ms`;
  dis$('disTBem').textContent = `${disFmt(r.tiempos.bem_s * 1e3, 0)} ms`;

  // sello: sale de los números de esta corrida, no de una constante
  const ok = m.fem_vs_serie_pct < CRITERIO_DISP && m.bem_vs_serie_pct < CRITERIO_DISP;
  dis$('disSello').textContent = ok
    ? `los dos métodos quedan dentro del ${CRITERIO_DISP} % de la serie de Bessel (criterio Fase 4). ` +
      `El que más se aleja es ${m.fem_vs_serie_pct > m.bem_vs_serie_pct ? 'FEM' : 'BEM'}.`
    : `fuera del criterio de la Fase 4 (${CRITERIO_DISP} %): subí ppw o la cantidad de paneles antes de usar estos números.`;
}

/* ------------------------------------------------------------------ dibujo */

const DIS_CONFIG = { displayModeBar: false, responsive: true, displaylogo: false };

function disCirculo(radioMm) {
  const th = [];
  const r = [];
  for (let i = 0; i <= 120; i += 1) {
    const t = (2 * Math.PI * i) / 120;
    th.push(radioMm * Math.cos(t));
    r.push(radioMm * Math.sin(t));
  }
  return {
    type: 'scatter', mode: 'lines', x: th, y: r, hoverinfo: 'skip',
    line: { color: '#f0883e', width: 1.2 }, name: 'cilindro',
  };
}

function disDibujar(r) {
  Plotly.react('mapaDisp', [{
    type: 'heatmap',
    x: r.mapa.x_mm, y: r.mapa.y_mm, z: r.mapa.db,
    colorscale: 'Inferno', zmin: -30, zmax: 0,
    colorbar: { title: { text: 'dB rel. pico', side: 'top' }, thickness: 11, len: 0.8,
                tickfont: { size: 9, color: '#8b95a1' }, outlinewidth: 0 },
    hovertemplate: 'x %{x:.1f} mm<br>y %{y:.1f} mm<br>%{z:.1f} dB<extra></extra>',
  }, disCirculo(r.mapa.radio_mm)], {
    margin: { l: 52, r: 8, t: 10, b: 40 },
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
    font: { color: '#8b95a1', size: 10 },
    // sin scaleanchor/scaleratio: cuando el div queda angosto Plotly falla con
    // "Something went wrong with axis scaling" (D27). El círculo se ve ovalado si el div no es
    // cuadrado, pero el mapa se dibuja siempre.
    xaxis: { title: { text: 'x (mm)' }, gridcolor: '#21262d', zeroline: false },
    yaxis: { title: { text: 'y (mm)' }, gridcolor: '#21262d', zeroline: false },
    showlegend: false,
  }, DIS_CONFIG);

  const c = r.comparacion;
  const traza = (nombre, y, color, dash) => ({
    type: 'scatter', mode: 'lines', x: c.angulo_deg, y,
    line: { color, width: 1.5, dash }, name: nombre,
    hovertemplate: `${nombre}<br>θ %{x:.0f}°<br>|p| %{y:.3f}<extra></extra>`,
  });
  Plotly.react('comparacionDisp', [
    traza('serie (referencia)', c.serie, '#8b95a1', 'dot'),
    traza('BEM', c.bem, '#58a6ff', 'solid'),
    traza('FEM', c.fem, '#f0883e', 'solid'),
  ], {
    ...DIS_CONFIG,
    margin: { l: 52, r: 8, t: 26, b: 40 },
    title: { text: `|p| en r = ${disFmt(c.r_eval_mm, 2)} mm: los tres métodos`, font: { size: 10, color: '#8b95a1' } },
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
    font: { color: '#8b95a1', size: 9 },
    legend: { orientation: 'h', y: 1.14, font: { size: 9 } },
    xaxis: { title: { text: 'θ (grados)' }, gridcolor: '#21262d' },
    yaxis: { title: { text: '|p|' }, gridcolor: '#21262d' },
  });
}

/* --------------------------------------------------------------- enlaces */

function disIniciar() {
  if (!DIS.resuelto) disEstimar();
  setTimeout(() => window.dispatchEvent(new Event('resize')), 30);
}
window.disIniciar = disIniciar;

dis$('disKa').addEventListener('input', (e) => {
  DIS.ka = parseFloat(e.target.value);
  dis$('disKaVal').textContent = disFmt(DIS.ka, 1);
});
dis$('disPpw').addEventListener('input', (e) => {
  DIS.ppw = parseInt(e.target.value, 10);
  dis$('disPpwVal').textContent = String(DIS.ppw);
});
dis$('disRout').addEventListener('input', (e) => {
  DIS.rOut = parseFloat(e.target.value);
  dis$('disRoutVal').textContent = disFmt(DIS.rOut, 1);
});
dis$('disAbc').addEventListener('change', (e) => { DIS.abc = e.target.value; });
dis$('disEcuacion').addEventListener('change', (e) => { DIS.ecuacion = e.target.value; });
dis$('btnDisEstimar').addEventListener('click', disEstimar);
dis$('btnDisolver').addEventListener('click', disResolver);
