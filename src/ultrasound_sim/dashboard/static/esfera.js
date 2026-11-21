/* Vista "esfera elástica" (U7).
 *
 * La física la resuelve `dashboard/elastic_sphere.py` y el sello lo evalúa `validacion.py`: acá sólo
 * se pide, se muestra el progreso y se pintan los resultados. El único número que se calcula en el
 * cliente es el costo estimado, que también viene del servidor (`/api/esfera/cost`).
 */

const ESF = { ka: 1.0, freq_khz: 100, modo: 'precisa', jobId: null, timer: null, listo: false };

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
    esf$('btnEsfera').addEventListener('click', esfLanzar);
    esf$('btnEsferaCancelar').addEventListener('click', esfCancelar);
    ESF.inicializado = true;
  }
  esfEstimar();
};

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

function esfPintar(r) {
  const m = r.metrica;

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
  Plotly.react('angEsf', [
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
  esf$('esfMetricas').textContent =
    `${m.material} · contraste ${m.contraste_impedancia.toFixed(2)}× · ka ${m.ka.toFixed(2)} · `
    + `radio ${m.radio_mm.toFixed(2)} mm · ${m.nodos.toLocaleString('es')} nodos `
    + `(${m.elementos_solidos} tets sólidos) · ${m.segundos.toFixed(1)} s · `
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
  } else if (m.material === 'hueso' && sello.estado === 'fuera') {
    esf$('esfSello').textContent +=
      ' OJO: el hueso (contraste 4.6×) tiene un error de MODELO, no de malla: con contraste'
      + ' moderado el cuerpo se mueve parcialmente con el fluido y la restricción deja de ser'
      + ' exacta (D47). Refinar no lo baja.';
  }
  if (m.modo === 'rapida' && sello.estado === 'fuera') {
    esf$('esfSello').textContent +=
      ' El modo rápido usa un dominio de 2.5a: el error es del ABC, no del solver (D46). '
      + 'Para un número que se reporte, usá la malla precisa.';
  }
  esf$('esfEstado').textContent = 'listo';
}
