/* Vista "esfera elástica" (U7).
 *
 * La física la resuelve `dashboard/elastic_sphere.py` y el sello lo evalúa `validacion.py`: acá sólo
 * se pide, se muestra el progreso y se pintan los resultados. El único número que se calcula en el
 * cliente es el costo estimado, que también viene del servidor (`/api/esfera/cost`).
 */

const ESF = {
  ka: 1.0, freq_khz: 100, modo: 'precisa', jobId: null, timer: null, listo: false,
  ultimo: null,        // el resultado completo de la última corrida, para exportar
  historial: [],       // corridas guardadas en localStorage
  superpuestas: new Set(),  // índices del historial que se dibujan encima
};
const ESF_CLAVE = 'esf.historial.v1';
const ESF_MAX_HISTORIAL = 12;

function esf$(id) { return document.getElementById(id); }

const ESF_LAYOUT = {
  paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
  font: { color: '#8b95a1', size: 9 },
  margin: { l: 52, r: 10, t: 26, b: 40 },
};

window.esfIniciar = function () {
  if (ESF.listo) return;
  ESF.listo = true;
  if (!ESF.inicializado) {
    esf$('esfKa').addEventListener('input', (e) => {
      ESF.ka = parseFloat(e.target.value);
      esf$('esfKaVal').textContent = ESF.ka.toFixed(2);
      esfEstimar();
    });
    esf$('esfFreq').addEventListener('input', (e) => {
      ESF.freq_khz = parseInt(e.target.value, 10);
      esf$('esfFreqVal').textContent = ESF.freq_khz + ' kHz';
      esfEstimar();
    });
    esf$('esfModo').addEventListener('change', (e) => { ESF.modo = e.target.value; esfEstimar(); });
    // En la dispersada no hay modos rígidos que restringir: pedirlos sería contradictorio y el
    // servidor lo rechaza. Se deshabilita el selector en vez de dejar que el usuario elija en vano.
    esf$('esfFormulacion').addEventListener('change', (e) => {
      const dispersada = e.target.value === 'dispersada';
      esf$('esfRigidos').disabled = dispersada;
      esf$('esfRigidos').title = dispersada
        ? 'No aplica: en la formulación dispersada el desplazamiento medio de u_s ya es el de la física.'
        : '';
    });
    esf$('btnEsfera').addEventListener('click', esfLanzar);
    esf$('btnEsferaCancelar').addEventListener('click', esfCancelar);
    esf$('btnEsfJson').addEventListener('click', () => esfDescargarJson());
    esf$('btnEsfCsv').addEventListener('click', () => esfDescargarCsv());
    esf$('btnEsfGuardar').addEventListener('click', esfGuardar);
    esf$('btnEsfBarrido').addEventListener('click', esfBarrido);
    esf$('btnEsfBarridoCancelar').addEventListener('click', esfBarridoCancelar);
    esfCargarHistorial();
    esfPintarHistorial();
    ESF.inicializado = true;
  }
  esfEstimar();
  // resultado de la sesión anterior: repintarlo en vez de dejar la vista en blanco y obligar a
  // volver a esperar los ~30 s del modo preciso
  const previo = window.vistaCargar('esfera');
  if (previo && !ESF.ultimo) {
    try {
      // `persistir = false`: restaurar no es resolver. Si se volviera a guardar, la fecha se
      // renovaría en cada recarga y el TTL de 12 h no vencería nunca, que es justo lo que evita que
      // un resultado viejo se muestre como reciente.
      esfPintar(previo.resultado, false);
      esf$('esfEstado').textContent =
        'restaurado de ' + new Date(previo.fecha).toLocaleString('es');
    } catch (e) {
      console.warn('no se pudo restaurar la esfera:', e && e.message);
      window.vistaBorrar('esfera');
    }
  }
};

// ---------------------------------------------------------------------------------------------
// exportación
// ---------------------------------------------------------------------------------------------
function esfDescargar(nombre, texto, tipo) {
  const blob = new Blob([texto], { type: tipo });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = nombre;
  document.body.appendChild(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

function esfNombreBase() {
  const m = ESF.ultimo.metrica;
  const sello = ESF.ultimo.sello.estado;
  return `esfera_${m.material}_ka${m.ka.toFixed(2)}_${m.formulacion}_${m.modo}_${sello}`.replace(
    /[^a-zA-Z0-9._-]+/g, '_');
}

function esfDescargarJson() {
  if (!ESF.ultimo) { esf$('esfDescarga').textContent = 'todavía no hay resultados: lanzá una corrida.'; return; }
  // se agrega la fecha y la versión del formato: un JSON sin fecha no se puede ordenar después
  const salida = {
    formato: 'ultrasound-sim/dashboard/esfera/1',
    fecha: new Date().toISOString(),
    ...ESF.ultimo,
  };
  esfDescargar(`${esfNombreBase()}.json`, JSON.stringify(salida, null, 2), 'application/json');
  esf$('esfDescarga').textContent = `descargado ${esfNombreBase()}.json`;
}

function esfDescargarCsv() {
  if (!ESF.ultimo) { esf$('esfDescarga').textContent = 'todavía no hay resultados: lanzá una corrida.'; return; }
  const a = ESF.ultimo.angulos;
  const filas = ['theta_deg,abs_p_s_fem,abs_p_s_faran,delta_fase_deg'];
  for (let i = 0; i < a.theta_deg.length; i += 1) {
    const d = ((a.fem_fase_deg[i] - a.faran_fase_deg[i] + 540) % 360) - 180;
    filas.push(`${a.theta_deg[i].toFixed(4)},${a.fem_modulo[i].toExponential(8)},`
      + `${a.faran_modulo[i].toExponential(8)},${d.toFixed(4)}`);
  }
  esfDescargar(`${esfNombreBase()}.csv`, filas.join('\n'), 'text/csv');
  esf$('esfDescarga').textContent = `descargado ${esfNombreBase()}.csv (${a.theta_deg.length} puntos)`;
}

// ---------------------------------------------------------------------------------------------
// barrido en ka: dónde vale cada formulación
// ---------------------------------------------------------------------------------------------
async function esfBarrido() {
  const cuerpo = {
    desde: parseFloat(esf$('esfBarDesde').value),
    hasta: parseFloat(esf$('esfBarHasta').value),
    puntos: parseInt(esf$('esfBarPuntos').value, 10),
    material: esf$('esfMaterial').value,
    modo: ESF.modo,
    formulacion: esf$('esfFormulacion').value,
    frecuencia_hz: ESF.freq_khz * 1e3,
    puntos_angulo: 61,
  };
  esf$('btnEsfBarrido').disabled = true;
  esf$('btnEsfBarridoCancelar').disabled = false;
  esf$('esfBarResumen').textContent = '';
  esf$('esfBarProgreso').style.width = '0%';
  esf$('esfBarEstado').textContent = 'encolando…';
  try {
    const r = await fetch('/api/jobs/esfera_barrido', {
      method: 'POST', headers: { 'content-type': 'application/json' }, body: JSON.stringify(cuerpo),
    }).then((x) => x.json());
    if (!r.ok) throw new Error(r.error || 'no se pudo encolar');
    ESF.barridoJob = r.job.id;
    esfBarridoSeguir();
  } catch (e) {
    esf$('esfBarEstado').textContent = `no se pudo lanzar: ${e.message}`;
    esf$('btnEsfBarrido').disabled = false;
    esf$('btnEsfBarridoCancelar').disabled = true;
  }
}

async function esfBarridoCancelar() {
  if (!ESF.barridoJob) return;
  esf$('btnEsfBarridoCancelar').disabled = true;
  await fetch(`/api/jobs/${ESF.barridoJob}/cancel`, { method: 'POST' });
}

function esfBarridoSeguir() {
  clearInterval(ESF.barridoTimer);
  ESF.barridoTimer = setInterval(async () => {
    try {
      const j = await (await fetch(`/api/jobs/${ESF.barridoJob}`)).json();
      esf$('esfBarProgreso').style.width = `${Math.round(100 * j.progreso)}%`;
      esf$('esfBarEstado').textContent = `${j.estado} · ${Math.round(100 * j.progreso)} % · ${j.mensaje}`;
      if (j.estado === 'listo') {
        clearInterval(ESF.barridoTimer);
        esfPintarBarrido(j.resultado);
        esf$('btnEsfBarrido').disabled = false;
        esf$('btnEsfBarridoCancelar').disabled = true;
      } else if (j.estado === 'error' || j.estado === 'cancelado') {
        clearInterval(ESF.barridoTimer);
        esf$('esfBarEstado').textContent = j.estado === 'error' ? `error: ${j.error}` : 'cancelado';
        esf$('btnEsfBarrido').disabled = false;
        esf$('btnEsfBarridoCancelar').disabled = true;
      }
    } catch (e) {
      clearInterval(ESF.barridoTimer);
      esf$('esfBarEstado').textContent = `se cortó el seguimiento: ${e.message}`;
      esf$('btnEsfBarrido').disabled = false;
      esf$('btnEsfBarridoCancelar').disabled = true;
    }
  }, 800);
}

function esfPintarBarrido(r) {
  const b = r.barrido;
  const tol = r.resumen.tolerancia;
  const puntos = b.ka.map((ka, i) => ({
    ka, err: 100 * b.error_relativo[i], l2: 100 * b.error_l2_relativo[i], sello: b.sello[i],
  }));
  const color = (s) => (s === 'validado' ? '#3fb950' : (s === 'fuera' ? '#f85149' : '#d29922'));
  Plotly.react('barridoEsf', [
    {
      type: 'scatter', mode: 'lines+markers', name: 'error máx. (%)', x: b.ka,
      y: b.error_relativo.map((v) => 100 * v),
      line: { color: '#58a6ff', width: 1.6 },
      marker: { size: 7, color: puntos.map((p) => color(p.sello)) },
      hovertemplate: 'ka %{x:.2f}<br>error máx. %{y:.2f} %<extra></extra>',
    },
    {
      type: 'scatter', mode: 'lines+markers', name: 'error L2 (%)', x: b.ka,
      y: b.error_l2_relativo.map((v) => 100 * v),
      line: { color: '#d29922', width: 1.3, dash: 'dot' }, marker: { size: 4 },
      hovertemplate: 'ka %{x:.2f}<br>L2 %{y:.2f} %<extra></extra>',
    },
  ], {
    ...ESF_LAYOUT,
    title: {
      text: `barrido en ka · ${b.material} · ${b.formulacion} · ${b.modo}`
        + ` (verde = dentro del criterio, rojo = fuera)`,
      font: { size: 10, color: '#8b95a1' },
    },
    xaxis: { title: { text: 'ka = k₀·a' } },
    yaxis: { title: { text: 'error contra Faran (%)' } },
    legend: { orientation: 'h', y: 1.18, font: { size: 9 } },
    shapes: [{
      type: 'line', xref: 'paper', x0: 0, x1: 1, y0: 100 * tol, y1: 100 * tol,
      line: { color: '#f85149', width: 1, dash: 'dash' },
    }],
  });

  const s = r.resumen;
  esf$('esfBarEstado').textContent = 'listo';
  esf$('esfBarResumen').textContent =
    `${s.n_dentro} de ${b.n_puntos} puntos dentro del criterio (${(100 * s.tolerancia).toFixed(0)} %)`
    + (s.n_dentro ? `: ka ${s.dentro_de.map((v) => v.toFixed(2)).join(', ')}` : '')
    + ` · peor ${(100 * s.peor_error).toFixed(2)} % en ka ${s.peor_ka.toFixed(2)}`
    + ` · mejor ${(100 * s.mejor_error).toFixed(2)} % · ${s.segundos_totales.toFixed(0)} s en total`;
  window.vistaGuardar('esfera_barrido', r);
}

// ---------------------------------------------------------------------------------------------
// historial (localStorage: no hay backend de resultados y no hace falta para comparar en sesión)
// ---------------------------------------------------------------------------------------------
function esfCargarHistorial() {
  try {
    ESF.historial = JSON.parse(localStorage.getItem(ESF_CLAVE) || '[]');
    if (!Array.isArray(ESF.historial)) ESF.historial = [];
  } catch { ESF.historial = []; }
}

function esfGuardarHistorial() {
  try {
    localStorage.setItem(ESF_CLAVE, JSON.stringify(ESF.historial));
  } catch {
    // el navegador puede rechazar por cuota: se avisa y se recorta a la mitad
    ESF.historial = ESF.historial.slice(0, Math.ceil(ESF_MAX_HISTORIAL / 2));
    try { localStorage.setItem(ESF_CLAVE, JSON.stringify(ESF.historial)); } catch { /* nada */ }
  }
}

function esfGuardar() {
  if (!ESF.ultimo) { esf$('esfDescarga').textContent = 'todavía no hay resultados: lanzá una corrida.'; return; }
  const m = ESF.ultimo.metrica;
  ESF.historial.unshift({
    fecha: new Date().toISOString(),
    material: m.material, formulacion: m.formulacion, modo: m.modo, ka: m.ka,
    error_rel: m.error_relativo, error_abs: m.error_absoluto,
    razon: m.razon_amplitudes, contraste: m.contraste_impedancia,
    sello: ESF.ultimo.sello.estado,
    theta_deg: ESF.ultimo.angulos.theta_deg,
    fem_modulo: ESF.ultimo.angulos.fem_modulo,
  });
  ESF.historial = ESF.historial.slice(0, ESF_MAX_HISTORIAL);
  esfGuardarHistorial();
  ESF.superpuestas.clear();
  esfPintarHistorial();
}

function esfPintarHistorial() {
  const caja = esf$('esfHistorial');
  if (!caja) return;
  if (!ESF.historial.length) {
    caja.innerHTML = '<p class="hint">todavía no guardaste ninguna corrida</p>';
    return;
  }
  caja.innerHTML = '';
  ESF.historial.forEach((h, i) => {
    const fila = document.createElement('label');
    fila.className = 'hist-fila';
    const cb = document.createElement('input');
    cb.type = 'checkbox';
    cb.checked = ESF.superpuestas.has(i);
    cb.addEventListener('change', () => {
      if (cb.checked) ESF.superpuestas.add(i); else ESF.superpuestas.delete(i);
      if (ESF.ultimo) esfPintar(ESF.ultimo);
    });
    const t = document.createElement('span');
    const cuando = h.fecha.slice(11, 19);
    t.textContent = `${cuando} · ${h.material} · ${h.formulacion} · ka ${h.ka.toFixed(2)}`
      + ` · ${(100 * h.error_rel).toFixed(2)} % · ${h.sello}`;
    fila.appendChild(cb);
    fila.appendChild(t);
    caja.appendChild(fila);
  });
  const limpiar = document.createElement('button');
  limpiar.className = 'btn mini';
  limpiar.textContent = 'borrar historial';
  limpiar.addEventListener('click', esfBorrarHistorial);
  caja.appendChild(limpiar);
}

function esfBorrarHistorial() {
  ESF.historial = [];
  ESF.superpuestas.clear();
  esfGuardarHistorial();
  esfPintarHistorial();
}

async function esfEstimar() {
  try {
    const r = await fetch(`/api/esfera/cost?modo=${encodeURIComponent(ESF.modo)}`
      + `&ka=${ESF.ka}&frecuencia_hz=${ESF.freq_khz * 1e3}`).then((x) => x.json());
    if (!r.ok) throw new Error(r.error || 'sin dato');
    const n = r.nodos_estimados.toLocaleString('es');
    esf$('esfCosto').textContent =
      `≈ ${n} nodos · radio ${(r.radio_m * 1e3).toFixed(2)} mm · elemento `
      + `${(r.elemento_m * 1e3).toFixed(3)} mm · ≈ ${r.segundos_estimados} s`
      + (r.aviso ? ` — ${r.aviso}` : '');
  } catch (e) {
    esf$('esfCosto').textContent = `no se pudo estimar: ${e.message}`;
  }
}

async function esfLanzar() {
  esf$('btnEsfera').disabled = true;
  esf$('btnEsferaCancelar').disabled = false;
  esf$('esfProgreso').style.width = '0%';
  esf$('esfEstado').textContent = 'encolando…';
  esf$('esfSello').textContent = '';
  esf$('esfMetricas').textContent = '';
  const cuerpo = {
    ka: ESF.ka, frecuencia_hz: ESF.freq_khz * 1e3, modo: ESF.modo,
    material: esf$('esfMaterial').value,
    fijar_rigidos: esf$('esfRigidos').value === 'si',
    formulacion: esf$('esfFormulacion').value,
    puntos_angulo: 121,
  };
  try {
    const r = await fetch('/api/jobs/esfera', {
      method: 'POST', headers: { 'content-type': 'application/json' },
      body: JSON.stringify(cuerpo),
    }).then((x) => x.json());
    if (!r.ok) throw new Error(r.error || 'no se pudo encolar');
    ESF.jobId = r.job.id;
    esfSeguir();
  } catch (e) {
    esf$('esfEstado').textContent = `no se pudo lanzar: ${e.message}`;
    esf$('btnEsfera').disabled = false;
    esf$('btnEsferaCancelar').disabled = true;
  }
}

function esfSeguir() {
  clearInterval(ESF.timer);
  ESF.timer = setInterval(async () => {
    try {
      const j = await (await fetch(`/api/jobs/${ESF.jobId}`)).json();
      esf$('esfProgreso').style.width = `${Math.round(100 * j.progreso)}%`;
      esf$('esfEstado').textContent = `${j.estado} · ${Math.round(100 * j.progreso)} % · ${j.mensaje}`;
      if (j.estado === 'listo') {
        clearInterval(ESF.timer);
        esfPintar(j.resultado);
        esf$('btnEsfera').disabled = false;
        esf$('btnEsferaCancelar').disabled = true;
      } else if (j.estado === 'error' || j.estado === 'cancelado') {
        clearInterval(ESF.timer);
        esf$('esfEstado').textContent = j.estado === 'error' ? `error: ${j.error}` : 'cancelado';
        esf$('btnEsfera').disabled = false;
        esf$('btnEsferaCancelar').disabled = true;
      }
    } catch (e) {
      clearInterval(ESF.timer);
      esf$('esfEstado').textContent = `se cortó el seguimiento: ${e.message}`;
      esf$('btnEsfera').disabled = false;
      esf$('btnEsferaCancelar').disabled = true;
    }
  }, 700);
}

async function esfCancelar() {
  if (!ESF.jobId) return;
  esf$('btnEsferaCancelar').disabled = true;
  await fetch(`/api/jobs/${ESF.jobId}/cancel`, { method: 'POST' });
}

function esfPintar(r, persistir = true) {
  const m = r.metrica;
  ESF.ultimo = r;   // queda disponible para exportar y para volver a pintar al cambiar superposiciones
  if (persistir) window.vistaGuardar('esfera', r);

  // ---- mapa |p_s| en el plano x-z (el interior de la esfera va enmascarado, no rellenado)
  const mapa = r.mapa;
  Plotly.react('mapaEsf', [{
    type: 'heatmap', x: mapa.x_mm, y: mapa.z_mm, z: mapa.modulo,
    colorscale: 'Viridis', hoverongaps: false,
    colorbar: { title: { text: '|p_s|', side: 'right' }, thickness: 10, len: 0.8 },
    hovertemplate: 'x %{x:.2f} mm · z %{y:.2f} mm<br>|p_s| = %{z:.3f}<extra></extra>',
  }], {
    ...ESF_LAYOUT, margin: { l: 52, r: 10, t: 26, b: 44 },
    title: { text: 'campo dispersado, plano x-z (la esfera no tiene fluido adentro)',
             font: { size: 10, color: '#8b95a1' } },
    xaxis: { title: { text: 'x (mm)' }, constrain: 'domain' },
    yaxis: { title: { text: 'z (mm)' }, constrain: 'domain', scaleanchor: 'x' },
    shapes: [{
      type: 'circle', xref: 'x', yref: 'y', x0: -mapa.radio_mm, x1: mapa.radio_mm,
      y0: -mapa.radio_mm, y1: mapa.radio_mm, line: { color: '#f0883e', width: 1.2, dash: 'dot' },
    }],
  });

  // ---- módulo contra ángulo: FEM contra Faran
  const a = r.angulos;
  // corridas guardadas que el usuario marcó: van punteadas y finas para no tapar la actual
  const trazasHist = [...ESF.superpuestas]
    .map((i) => ESF.historial[i])
    .filter(Boolean)
    .map((h, k) => ({
      type: 'scatter', mode: 'lines', name: `guardada ${h.fecha.slice(11, 19)} ${h.material} ${h.formulacion}`,
      x: h.theta_deg, y: h.fem_modulo,
      line: { color: `hsl(${(k * 67 + 280) % 360} 70% 62%)`, width: 1.1, dash: 'dot' },
    }));
  Plotly.react('angEsf', [
    ...trazasHist,
    { type: 'scatter', mode: 'lines+markers', name: 'FEM acoplado', x: a.theta_deg, y: a.fem_modulo,
      line: { color: '#58a6ff', width: 1.6 }, marker: { size: 3 } },
    { type: 'scatter', mode: 'lines', name: 'Faran (exacto)', x: a.theta_deg, y: a.faran_modulo,
      line: { color: '#3fb950', width: 1.4, dash: 'dash' } },
  ], {
    ...ESF_LAYOUT,
    title: { text: '|p_s| contra θ en r = 1.5a', font: { size: 10, color: '#8b95a1' } },
    xaxis: { title: { text: 'θ (°)' } }, yaxis: { title: { text: '|p_s|' } },
    legend: { orientation: 'h', y: 1.18, font: { size: 9 } },
  });

  // ---- diferencia de fase, no las dos fases: las curvas de fase tienen saltos de ±180° que se ven
  // como discontinuidades del modelo y no lo son. La diferencia (envuelta a ±180°) es lo que se
  // quiere leer: si el modelo tuviera la fase corrida, acá se ve como una banda horizontal.
  const difFase = a.fem_fase_deg.map((v, i) => ((v - a.faran_fase_deg[i] + 540) % 360) - 180);
  Plotly.react('faseEsf', [
    { type: 'scatter', mode: 'lines+markers', name: 'FEM − Faran', x: a.theta_deg, y: difFase,
      line: { color: '#d29922', width: 1.4 }, marker: { size: 3 } },
  ], {
    ...ESF_LAYOUT,
    title: { text: 'diferencia de fase FEM − Faran (0° = mismo desfase)',
             font: { size: 10, color: '#8b95a1' } },
    xaxis: { title: { text: 'θ (°)' } },
    yaxis: { title: { text: 'Δ fase (°)' }, range: [-30, 30] },
    showlegend: false,
    shapes: [{ type: 'line', xref: 'paper', x0: 0, x1: 1, y0: 0, y1: 0,
               line: { color: '#3fb950', width: 1, dash: 'dot' } }],
  });

  const resid = m.residuo_rigido === null ? '—' : m.residuo_rigido.toExponential(1);
  const formul = m.formulacion === 'dispersada'
    ? 'dispersada (u_s = u − u_i, operador de contraste)'
    : 'total + traslación media nula';
  esf$('esfMetricas').textContent =
    `${m.material} · ${formul} · contraste ${m.contraste_impedancia.toFixed(2)}× · ka ${m.ka.toFixed(2)} · `
    + `radio ${m.radio_mm.toFixed(2)} mm · ${m.nodos.toLocaleString('es')} nodos `
    + `(${m.elementos_solidos} tets sólidos) · ${m.segundos.toFixed(1)} s · `
    + `error máx. ${(100 * m.error_relativo).toFixed(2)} % · L2 `
    + `${(100 * (m.error_l2_relativo ?? m.error_relativo)).toFixed(2)} % · `
    + `razón de amplitudes ${m.razon_amplitudes.toFixed(3)} · error absoluto `
    + `${m.error_absoluto.toExponential(2)} Pa (incidente = 1 Pa) · `
    + `fuerza neta no cancelada ${resid} N`
    + (m.restriccion_automatica ? ' · restricción de rígidos DESACTIVADA automáticamente' : '');

  const sello = r.sello;
  esf$('esfSello').textContent = sello.texto;
  esf$('esfSello').className = `hint sello-${sello.estado}`;
  // Los materiales de bajo contraste tienen un error de MODELO conocido (D47), no de malla: el
  // cuerpo se mueve parcialmente con el fluido y la condición de traslación media nula deja de ser
  // exacta. Se dice acá en vez de subir el umbral del criterio hasta que entre.
  if (m.elementos_por_lambda_corte < 4) {
    esf$('esfSello').textContent +=
      ` ATENCIÓN: la malla resuelve la onda de corte con sólo `
      + `${m.elementos_por_lambda_corte.toFixed(1)} elementos por longitud de onda (λ_s = `
      + `${m.lambda_corte_mm.toFixed(2)} mm, elemento ${m.elemento_mm.toFixed(3)} mm). El tamaño se`
      + ` elige con la longitud de onda del FLUIDO; un sólido blando necesita malla propia más fina.`;
  }
  if (m.contraste_impedancia < 2.5) {
    esf$('esfSello').textContent +=
      ` Contraste ${m.contraste_impedancia.toFixed(2)}× el agua: acá el cuerpo se mueve en parte CON`
      + ' el fluido, la condición de traslación media nula deja de valer y se desactiva sola. El'
      + ' error RELATIVO no es un indicador útil porque el campo dispersado es muy débil: mirá el'
      + ' error absoluto, que es contra la incidente de 1 Pa.';
  } else if (sello.estado === 'fuera' && m.formulacion === 'dispersada') {
    // D49: la formulación correcta es peor contra Faran en esta malla, y el patrón medido es que el
    // error crece al bajar el contraste (tungsteno 10.70 %, acero 19.49 %, aluminio 44.14 %, hueso
    // 54.53 %): el término de contraste es una diferencia de operadores grandes y el error de
    // interpolación de la incidente no se achica al bajar el contraste.
    esf$('esfSello').textContent +=
      ' Esto es la formulación dispersada: es la correcta en el continuo (traslación 3/(2g+1) al'
      + ' 1-2 %, caso transparente exacto) pero su término de contraste es una resta de operadores'
      + ' grandes y el error de interpolación de la incidente no se achica al bajar el contraste.'
      + ' Medido: tungsteno 10.70 %, acero 19.49 %, aluminio 44.14 %, hueso 54.53 % (D49). El número'
      + ' que se reporta hoy es el de la formulación total.';
  } else if (sello.estado === 'fuera' && m.fijar_rigidos) {
    // D48: la restricción suprime una respuesta dipolar que Faran SÍ incluye (Feuillade y Clay 1999,
    // "Anderson (1950) revisited", §IV). El error crece con 3/(2·g+1), g = ρ_sólido/ρ_fluido: chico
    // en metales, grande en tejido. Es aproximación de la implementación, no del método.
    esf$('esfSello').textContent +=
      ' OJO: la restricción de traslación media nula suprime una respuesta dipolar que Faran sí'
      + ' incluye (Feuillade y Clay 1999, §IV): el error crece con 3/(2g+1), g = ρ_sólido/ρ_fluido,'
      + ' así que se nota más cuanto más liviano es el material. Es una limitación de esta'
      + ' implementación, no del método, y refinar la malla no la baja (D48).';
  }
  if (m.modo === 'rapida' && sello.estado === 'fuera') {
    esf$('esfSello').textContent +=
      ' El modo rápido usa un dominio de 2.5a: el error es del ABC, no del solver (D46). '
      + 'Para un número que se reporte, usá la malla precisa.';
  }
  if (persistir) esf$('esfEstado').textContent = 'listo';
  esfPintarHistorial();
}
