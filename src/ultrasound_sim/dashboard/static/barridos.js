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
      window.cargando('graficoBarrido', `barriendo · ${Math.round(100 * j.progreso)} %`, j.mensaje);
      if (j.estado === 'listo') {
        clearInterval(BAR.temporizador);
        BAR.resultado = j.resultado;
        window.terminarCargando('graficoBarrido');
        barPintar(j.resultado, j.segundos);
        // 30 puntos de escena son ~30 solves: guardarlo para no repetirlos al volver a la pestaña
        window.vistaGuardar('barridos', { resultado: j.resultado, segundos: j.segundos });
        bar$('btnBarrido').disabled = false;
        bar$('btnCancelar').disabled = true;
      } else if (j.estado === 'error' || j.estado === 'cancelado') {
        clearInterval(BAR.temporizador);
        window.terminarCargando('graficoBarrido');
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

  // sello del servidor (U6). Ojo: si no viene, NO se declara éxito — el `Math.max(...[], 0)` de la
  // versión anterior devolvía 0 y anunciaba "desvío 0 %" con la tabla vacía.
  mostrarSello('barSello', r.sello,
    `el barrido terminó (${r.puntos} puntos, ${barFmt(segundos, 0)} s) pero no trajo criterio: ` +
    'revisá el resultado antes de usarlo.');

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

  // La curva de sensibilidad va al gráfico de barras sólo si NO se corrió el barrido de frecuencia
  // (ese usa el gráfico principal). Sin esta guarda, el último que corre pisa al otro.
  if (BARESC.jobId) return;

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

/* ------------------------------------------------ barrido de apertura (FWHM contra pitch) */

async function barApertura() {
  // `.value` y no el elemento: `parseInt(<input>, 10)` devuelve NaN, el JSON lo manda como null y el
  // servidor muere con "float() argument must be ... not 'NoneType'".
  const valor = (id) => esfBarId(id).value;
  const cuerpo = {
    desde: parseInt(valor('barApDesde'), 10),
    hasta: parseInt(valor('barApHasta'), 10),
    puntos: parseInt(valor('barApPuntos'), 10),
    apertura_mm: parseFloat(valor('barApApertura')),
    foco_mm: parseFloat(valor('barApFoco')),
    frecuencia_khz: parseFloat(valor('barApFreq')),
    material: valor('barApMaterial'), modelo: 'baffle',
  };
  bar$('btnBarApertura').disabled = true;
  bar$('btnBarAperturaCancelar').disabled = false;
  // el spinner va antes del fetch: encolar tarda y sin esto parece que el click no hizo nada
  window.cargando('graficoApertura', 'encolando el barrido de apertura…', '');
  bar$('barApProgreso').style.width = '0%';
  bar$('barApEstado').textContent = 'encolando…';
  bar$('barApResumen').textContent = '';
  try {
    const r = await barPostJson('/api/jobs/barrido_apertura', cuerpo);
    BAR.apJob = r.job.id;
    barAperturaSeguir();
  } catch (e) {
    bar$('barApEstado').textContent = `no se pudo lanzar: ${e.message}`;
    bar$('btnBarApertura').disabled = false;
    bar$('btnBarAperturaCancelar').disabled = true;
  }
}

function esfBarId(id) { return document.getElementById(id); }

async function barPostJson(ruta, cuerpo) {
  const r = await fetch(ruta, {
    method: 'POST', headers: { 'content-type': 'application/json' }, body: JSON.stringify(cuerpo),
  }).then((x) => x.json());
  if (!r.ok) throw new Error(r.error || 'no se pudo encolar');
  return r;
}

async function barAperturaCancelar() {
  if (!BAR.apJob) return;
  bar$('btnBarAperturaCancelar').disabled = true;
  await fetch(`/api/jobs/${BAR.apJob}/cancel`, { method: 'POST' });
}

function barAperturaSeguir() {
  clearInterval(BAR.apTimer);
  BAR.apTimer = setInterval(async () => {
    try {
      const j = await barGet(`/api/jobs/${BAR.apJob}`);
      bar$('barApProgreso').style.width = `${Math.round(100 * j.progreso)}%`;
      bar$('barApEstado').textContent = `${j.estado} · ${Math.round(100 * j.progreso)} % · ${j.mensaje}`;
      window.cargando('graficoApertura', `barriendo apertura · ${Math.round(100 * j.progreso)} %`,
        j.mensaje);
      if (j.estado === 'listo') {
        clearInterval(BAR.apTimer);
        window.terminarCargando('graficoApertura');
        barPintarApertura(j.resultado);
        bar$('btnBarApertura').disabled = false;
        bar$('btnBarAperturaCancelar').disabled = true;
      } else if (j.estado === 'error' || j.estado === 'cancelado') {
        clearInterval(BAR.apTimer);
        window.terminarCargando('graficoApertura');
        bar$('barApEstado').textContent = j.estado === 'error' ? `error: ${j.error}` : 'cancelado';
        bar$('btnBarApertura').disabled = false;
        bar$('btnBarAperturaCancelar').disabled = true;
      }
    } catch (e) {
      clearInterval(BAR.apTimer);
      bar$('barApEstado').textContent = `se cortó el seguimiento: ${e.message}`;
      bar$('btnBarApertura').disabled = false;
      bar$('btnBarAperturaCancelar').disabled = true;
    }
  }, 600);
}

function barPintarApertura(r) {
  const b = r.barrido;
  const ref = r.referencia_continua;
  const cuatro = b.pitch_over_lambda.map((x) => x * 4);
  Plotly.react('graficoApertura', [
    {
      type: 'scatter', mode: 'lines+markers', name: 'FWHM lateral (mm)',
      x: b.pitch_over_lambda, y: b.fwhm_lateral_mm,
      line: { color: '#58a6ff', width: 1.8 }, marker: { size: 7 },
      hovertemplate: 'pitch %{x:.3f} λ<br>FWHM %{y:.2f} mm<extra></extra>',
    },
    {
      type: 'scatter', mode: 'lines', name: 'apertura continua (referencia)',
      x: b.pitch_over_lambda, y: b.pitch_over_lambda.map(() => ref.fwhm_lateral_mm),
      line: { color: '#3fb950', width: 1.2, dash: 'dash' }, hoverinfo: 'skip',
    },
  ], {
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
    font: { color: '#8b95a1', size: 9 }, margin: { l: 52, r: 12, t: 28, b: 40 },
    title: {
      text: `${b.apertura_mm} mm de apertura · foco ${b.foco_mm} mm · ${b.frecuencia_khz} kHz`
        + ' (verde punteado: apertura continua)',
      font: { size: 10, color: '#8b95a1' },
    },
    xaxis: { title: { text: 'pitch / λ' } },
    yaxis: { title: { text: 'FWHM lateral (mm)' } },
    legend: { orientation: 'h', y: 1.2, font: { size: 9 } },
    shapes: [{
      type: 'line', xref: 'paper', x0: 0, x1: 1, y0: 0, y1: 0, opacity: 0,
    }],
  });

  const s = r.resumen;
  bar$('barApEstado').textContent = 'listo';
  bar$('barApResumen').textContent =
    `referencia continua: FWHM ${ref.fwhm_lateral_mm.toFixed(2)} mm · `
    + (s.pitch_con_2pct === null
      ? 'ningún pitch de este rango empata dentro del 2 % (usá más elementos)'
      : `el pitch más grueso dentro del 2 % es ${s.pitch_con_2pct.toFixed(3)} λ `
        + `(${(s.pitch_con_2pct / 0.25).toFixed(2)}× λ/4, que es el criterio de D25)`)
    + ` · peor desvío ${(100 * s.peor_desvio).toFixed(1)} % con n = ${s.peor_n}`
    + ` · ${s.segundos_totales.toFixed(0)} s`;
  window.vistaGuardar('barrido_apertura', r);
}

function barIniciar() {
  // listeners del barrido de apertura (idempotente: `barIniciar` se llama en cada cambio de vista)
  const ba = bar$('btnBarApertura');
  if (ba && !ba.dataset.listo) {
    ba.dataset.listo = '1';
    ba.addEventListener('click', barApertura);
    bar$('btnBarAperturaCancelar').addEventListener('click', barAperturaCancelar);
    const previoAp = window.vistaCargar('barrido_apertura');
    if (previoAp) {
      try {
        barPintarApertura(previoAp.resultado);
        bar$('barApEstado').textContent =
          'restaurado de ' + new Date(previoAp.fecha).toLocaleString('es');
      } catch (e) { console.warn('no se pudo restaurar el barrido de apertura:', e && e.message); }
    }
  }
  const previo = window.vistaCargar('barridos');
  if (previo && !BAR.resultado && previo.resultado && previo.resultado.resultado) {
    try {
      BAR.resultado = previo.resultado.resultado;
      barPintar(previo.resultado.resultado, previo.resultado.segundos);
      bar$('barEstado').textContent =
        'restaurado de ' + new Date(previo.fecha).toLocaleString('es');
    } catch (e) {
      console.warn('no se pudo restaurar el barrido:', e && e.message);
      window.vistaBorrar('barridos');
    }
  }
  if (!BAR.jobId) barInfo();
  setTimeout(() => window.dispatchEvent(new Event('resize')), 30);
}
window.barIniciar = barIniciar;

bar$('btnBarrido').addEventListener('click', barLanzar);
bar$('btnCancelar').addEventListener('click', barCancelar);

/* --------------------------------------------------- barrido de frecuencia de la escena 2D */

const BARESC = { desde: 100, hasta: 400, puntos: 30, jobId: null, temporizador: null };

function barescEscena() {
  // usa la escena de la pestaña "escena FEM" si el usuario la tocó; si no, una losa por defecto
  if (typeof ESC !== 'undefined' && ESC && ESC.capas) {
    return {
      frecuencia_hz: ESC.freq * 1e3,
      ppw: ESC.ppw,
      alto_m: ESC.alto * 1e-3,
      capas: ESC.capas.map((c) => ({ material: c.material, espesor_m: c.espesor_mm * 1e-3 })),
    };
  }
  return {
    frecuencia_hz: 100e3, ppw: 12, alto_m: 20e-3,
    capas: [{ material: 'agua', espesor_m: 20e-3 },
            { material: 'hueso', espesor_m: 10e-3 },
            { material: 'agua', espesor_m: 20e-3 }],
  };
}

async function barescLanzar() {
  if (BARESC.desde >= BARESC.hasta) {
    bar$('barEscEstado').textContent = 'la banda necesita desde < hasta';
    return;
  }
  bar$('btnBarridoEscena').disabled = true;
  bar$('btnCancelarEscena').disabled = false;
  bar$('barEscProgreso').style.width = '0%';
  bar$('barEscEstado').textContent = 'encolando…';
  const escena = barescEscena();
  bar$('barEscInfo').textContent =
    `${escena.capas.length} capas · ppw ${escena.ppw} · ` +
    escena.capas.map((c) => `${c.material} ${(c.espesor_m * 1e3).toFixed(0)} mm`).join(' / ');
  try {
    const r = await fetch('/api/jobs/barrido_escena', {
      method: 'POST',
      headers: { 'content-type': 'application/json' },
      body: JSON.stringify({ escena, desde_hz: BARESC.desde * 1e3,
                             hasta_hz: BARESC.hasta * 1e3, puntos: BARESC.puntos }),
    }).then((x) => x.json());
    if (!r.ok) throw new Error(r.error || 'no se pudo encolar');
    BARESC.jobId = r.job.id;
    barescSeguir();
  } catch (e) {
    bar$('barEscEstado').textContent = `no se pudo lanzar: ${e.message}`;
    bar$('btnBarridoEscena').disabled = false;
    bar$('btnCancelarEscena').disabled = true;
  }
}

function barescSeguir() {
  clearInterval(BARESC.temporizador);
  BARESC.temporizador = setInterval(async () => {
    try {
      const j = await (await fetch(`/api/jobs/${BARESC.jobId}`)).json();
      bar$('barEscProgreso').style.width = `${Math.round(100 * j.progreso)}%`;
      bar$('barEscEstado').textContent = `${j.estado} · ${Math.round(100 * j.progreso)} % · ${j.mensaje}`;
      if (j.estado === 'listo') {
        clearInterval(BARESC.temporizador);
        barescPintar(j.resultado);
        bar$('btnBarridoEscena').disabled = false;
        bar$('btnCancelarEscena').disabled = true;
      } else if (j.estado === 'error' || j.estado === 'cancelado') {
        clearInterval(BARESC.temporizador);
        bar$('barEscEstado').textContent = j.estado === 'error' ? `error: ${j.error}` : 'cancelado';
        bar$('btnBarridoEscena').disabled = false;
        bar$('btnCancelarEscena').disabled = true;
      }
    } catch (e) {
      clearInterval(BARESC.temporizador);
      bar$('barEscEstado').textContent = `se cortó el seguimiento: ${e.message}`;
      bar$('btnBarridoEscena').disabled = false;
      bar$('btnCancelarEscena').disabled = true;
    }
  }, 600);
}

function barescPintar(r) {
  const validas = r.filas.filter((f) => f.T !== null);
  const f_khz = validas.map((f) => f.f_hz / 1e3);
  const trazas = [
    { type: 'scatter', mode: 'lines+markers', name: 'T (FEM)', x: f_khz,
      y: validas.map((f) => f.T), line: { color: '#58a6ff', width: 1.6 },
      hovertemplate: '%{x:.0f} kHz<br>T = %{y:.3f}<extra></extra>' },
    { type: 'scatter', mode: 'lines', name: 'T (matriz de transferencia)', x: f_khz,
      y: validas.map((f) => f.T_analitico), line: { color: '#8b95a1', width: 1, dash: 'dot' } },
    { type: 'scatter', mode: 'lines', name: 'R (FEM)', x: f_khz,
      y: validas.map((f) => f.R), line: { color: '#f0883e', width: 1.4 } },
  ];
  const formas = [];
  if (r.prediccion) {
    for (const f of r.prediccion.resonancias_hz) {
      formas.push({ type: 'line', x0: f / 1e3, x1: f / 1e3, yref: 'paper', y0: 0, y1: 1,
                    line: { color: '#3fb950', width: 1, dash: 'dash' } });
    }
  }
  Plotly.react('graficoBarrido', trazas, {
    ...BAR_CONFIG,
    margin: { l: 52, r: 8, t: 26, b: 40 },
    title: { text: 'R y T vs frecuencia (verde: resonancias de media onda predichas)',
             font: { size: 10, color: '#8b95a1' } },
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
    font: { color: '#8b95a1', size: 9 },
    legend: { orientation: 'h', y: 1.16, font: { size: 9 } },
    xaxis: { title: { text: 'f (kHz)' }, gridcolor: '#21262d' },
    yaxis: { title: { text: 'R, T' }, gridcolor: '#21262d', range: [0, 1.05] },
    shapes: formas,
  });

  const picos = r.picos.map((p) => `${barFmt(p.f_hz / 1e3, 1)} kHz`).join(', ') || 'ninguno';
  let texto = `picos de transmisión: ${picos}.`;
  if (r.prediccion) {
    const pred = r.prediccion.resonancias_hz.map((f) => barFmt(f / 1e3, 1)).join(', ');
    texto += ` predichas para ${r.prediccion.espesor_mm} mm de ${r.prediccion.material} ` +
             `(c = ${r.prediccion.c_m_s} m/s): ${pred} kHz.`;
  }
  if (r.max_error_T !== null) {
    texto += ` peor error de T del barrido: ${barFmt(100 * r.max_error_T, 2)} %.`;
  }
  bar$('barEscSello').textContent = texto;
  if (r.sello) {
    const el = bar$('barEscSello');
    el.textContent += ' ' + r.sello.texto;
    el.className = `hint sello-${r.sello.estado}`;
  }
}

bar$('barEscDesde').addEventListener('input', (e) => {
  BARESC.desde = parseInt(e.target.value, 10);
  bar$('barEscDesdeVal').textContent = `${BARESC.desde} kHz`;
});
bar$('barEscHasta').addEventListener('input', (e) => {
  BARESC.hasta = parseInt(e.target.value, 10);
  bar$('barEscHastaVal').textContent = `${BARESC.hasta} kHz`;
});
bar$('barEscPuntos').addEventListener('input', (e) => {
  BARESC.puntos = parseInt(e.target.value, 10);
  bar$('barEscPuntosVal').textContent = String(BARESC.puntos);
});
bar$('btnBarridoEscena').addEventListener('click', barescLanzar);
bar$('btnCancelarEscena').addEventListener('click', async () => {
  if (!BARESC.jobId) return;
  bar$('btnCancelarEscena').disabled = true;
  await fetch(`/api/jobs/${BARESC.jobId}/cancel`, { method: 'POST' });
});
