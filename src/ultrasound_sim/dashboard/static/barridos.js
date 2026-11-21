/* Barridos con cola, progreso y cancelación (Fase U5) ----------------------------------------
 *
 * El barrido de sensibilidad son 31 solves de ~2.5 s: 88 s en total. Un `POST` síncrono dejaría la
 * pestaña colgada sin señal, así que acá el servidor encola un *job* y esta vista hace polling.
 *
 * El sello compara el ranking medido contra los pesos de D30 que ya están publicados: si la UI no
 * los reproduce, el que está mal es el barrido de la UI, no el solver.
 */

const BAR = {
  jobId: null,
  temporizador: null,
  resultado: null,
};

const bar$ = (id) => document.getElementById(id);
const barFmt = (v, d = 2) => (v === null || v === undefined || Number.isNaN(v) ? '—' : Number(v).toFixed(d));

const BAR_CONFIG = { displayModeBar: false, responsive: true, displaylogo: false };

async function barGet(ruta) {
  const r = await fetch(ruta);
  const datos = await r.json();
  if (datos.ok === false) throw new Error(datos.error || 'error del servidor');
  return datos;
}

async function barInfo() {
  try {
    const r = await barGet('/api/sweeps');
    const b = r.barridos[0];
    bar$('barInfo').textContent =
      `${b.puntos} puntos · estimado ${barFmt(b.estimado_s / 60, 1)} min · referencia ${b.referencia}`;
  } catch (e) {
    bar$('barInfo').textContent = `no se pudo leer la lista de barridos: ${e.message}`;
  }
}

async function barLanzar() {
  bar$('btnBarrido').disabled = true;
  bar$('btnCancelar').disabled = false;
  bar$('barEstado').textContent = 'encolando…';
  bar$('barProgreso').style.width = '0%';
  bar$('barRanking').innerHTML = '';
  bar$('barSello').textContent = '';
  try {
    const r = await fetch('/api/jobs/sensibilidad', {
      method: 'POST',
      headers: { 'content-type': 'application/json' },
      body: '{}',
    }).then((x) => x.json());
    if (!r.ok) throw new Error(r.error || 'no se pudo encolar');
    BAR.jobId = r.job.id;
    barSeguir();
  } catch (e) {
    bar$('barEstado').textContent = `no se pudo lanzar: ${e.message}`;
    bar$('btnBarrido').disabled = false;
    bar$('btnCancelar').disabled = true;
  }
}

async function barCancelar() {
  if (!BAR.jobId) return;
  bar$('btnCancelar').disabled = true;
  try {
    await fetch(`/api/jobs/${BAR.jobId}/cancel`, { method: 'POST' });
    bar$('barEstado').textContent = 'cancelación pedida…';
  } catch (e) {
    bar$('barEstado').textContent = `no se pudo cancelar: ${e.message}`;
  }
}

function barSeguir() {
  clearInterval(BAR.temporizador);
  BAR.temporizador = setInterval(async () => {
    try {
      const j = await barGet(`/api/jobs/${BAR.jobId}`);
      bar$('barProgreso').style.width = `${Math.round(100 * j.progreso)}%`;
      bar$('barEstado').textContent = `${j.estado} · ${Math.round(100 * j.progreso)} % · ${j.mensaje}`;
      if (j.estado === 'listo') {
        clearInterval(BAR.temporizador);
        BAR.resultado = j.resultado;
        barPintar(j.resultado, j.segundos);
        bar$('btnBarrido').disabled = false;
        bar$('btnCancelar').disabled = true;
      } else if (j.estado === 'error' || j.estado === 'cancelado') {
        clearInterval(BAR.temporizador);
        bar$('barEstado').textContent = j.estado === 'error'
          ? `error: ${j.error}` : 'cancelado';
        bar$('btnBarrido').disabled = false;
        bar$('btnCancelar').disabled = true;
      }
    } catch (e) {
      clearInterval(BAR.temporizador);
      bar$('barEstado').textContent = `se cortó el seguimiento: ${e.message}`;
      bar$('btnBarrido').disabled = false;
      bar$('btnCancelar').disabled = true;
    }
  }, 700);
}

function barPintar(r, segundos) {
  const comp = r.comparacion || [];
  bar$('barRanking').innerHTML = comp.map((f) => {
    const ok = f.d30_pct !== null && Math.abs(f.desvio_pct) < 5;
    return `<div class="rank-fila${ok ? '' : ' fuera'}">
      <span class="rank-nombre">${f.parametro}</span>
      <span class="rank-valor">${barFmt(f.rango_pct)} %</span>
      <span class="rank-ref">D30 ${f.d30_pct === null ? '—' : barFmt(f.d30_pct)} %</span>
    </div>`;
  }).join('');

  const desvios = comp.filter((f) => f.d30_pct !== null && f.desvio_pct !== null)
    .map((f) => Math.abs(f.desvio_pct));
  // sin comparación no se puede declarar nada: `Math.max(...[], 0)` da 0 y diría "desvío 0 %"
  // — un falso positivo que ya me comí una vez
  if (!desvios.length) {
    bar$('barSello').textContent =
      'el barrido terminó pero no trajo comparación contra D30: revisá el resultado antes de usarlo.';
    return;
  }
  const peor = Math.max(...desvios);
  bar$('barSello').textContent = peor < 5
    ? `el ranking reproduce el de D30: el peor desvío es ${barFmt(peor, 1)} % ` +
      `(${r.puntos} puntos, ${barFmt(segundos, 0)} s, t_interior base ${barFmt(r.base.t_interior, 5)}).`
    : `el ranking se aparta de D30 hasta un ${barFmt(peor, 1)} %: no usar estos pesos sin revisar el barrido.`;

  const colores = comp.map((f) => (f.d30_pct !== null && Math.abs(f.desvio_pct) < 5 ? '#3fb950' : '#f85149'));
  Plotly.react('barrasBarrido', [{
    type: 'bar', orientation: 'h',
    x: comp.map((f) => f.rango_pct),
    y: comp.map((f) => f.parametro),
    marker: { color: colores },
    text: comp.map((f) => `${barFmt(f.rango_pct)} %`),
    textposition: 'outside',
    hovertemplate: '%{y}: %{x:.2f} % del t_interior base<extra></extra>',
  }], {
    ...BAR_CONFIG,
    margin: { l: 92, r: 40, t: 26, b: 34 },
    title: { text: 'rango relativo de t_interior por parámetro (verde: coincide con D30)',
             font: { size: 10, color: '#8b95a1' } },
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
    font: { color: '#8b95a1', size: 9 },
    xaxis: { title: { text: '% del t_interior base' }, gridcolor: '#21262d', range: [0, Math.max(...comp.map((f) => f.rango_pct)) * 1.25] },
    yaxis: { autorange: 'reversed' },
    showlegend: false,
  });

  // curva de cada parámetro: t_interior normalizado contra el valor del parámetro
  const grupos = {};
  for (const fila of r.filas) {
    (grupos[fila.parametro] = grupos[fila.parametro] || []).push(fila);
  }
  const trazas = Object.entries(grupos).map(([nombre, filas]) => ({
    type: 'scatter', mode: 'lines+markers', name: nombre,
    x: filas.map((_, i) => i + 1), y: filas.map((f) => f.t_interior),
    text: filas.map((f) => f.valor),
    hovertemplate: `${nombre}: %{text}<br>t_interior %{y:.4f}<extra></extra>`,
  }));
  Plotly.react('graficoBarrido', trazas, {
    ...BAR_CONFIG,
    margin: { l: 52, r: 8, t: 26, b: 34 },
    title: { text: 't_interior en cada punto del barrido', font: { size: 10, color: '#8b95a1' } },
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
    font: { color: '#8b95a1', size: 9 },
    legend: { orientation: 'h', y: 1.18, font: { size: 9 } },
    xaxis: { title: { text: 'punto dentro de cada parámetro' }, gridcolor: '#21262d' },
    yaxis: { title: { text: 't_interior' }, gridcolor: '#21262d' },
  });
}

function barIniciar() {
  if (!BAR.jobId) barInfo();
  setTimeout(() => window.dispatchEvent(new Event('resize')), 30);
}
window.barIniciar = barIniciar;

bar$('btnBarrido').addEventListener('click', barLanzar);
bar$('btnCancelar').addEventListener('click', barCancelar);
