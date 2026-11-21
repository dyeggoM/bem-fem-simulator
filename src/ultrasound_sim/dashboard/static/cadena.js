/* Vista «cadena piel→cóclea»: capas de tejido real acopladas en una malla (U9).
 *
 * Es la única vista donde el modelo mezcla SÓLIDOS (hueso, con cizalla) y FLUIDOS (piel, líquidos) en el
 * mismo ensamblado. Por eso no alcanza con mostrar el mapa: el número que importa es la transmisión
 * medida contra la exacta, y cuando no hay exacta (porque hay cizalla) hay que decirlo en vez de inventar
 * una referencia.
 *
 * El estado vive en CAD y se guarda con la persistencia por vista, igual que las demás pestañas.
 */
(function () {
  'use strict';

  const $ = (id) => document.getElementById(id);
  // `vistaGuardar` ya antepone 'usim.vista.' y envuelve todo en {fecha, resultado}: acá va el nombre
  // pelado y el paquete entero como «resultado», no un objeto con su propio campo `resultado`
  const CLAVE = 'cadena';

  const CAD = {
    medios: {},          // nombre -> {etiqueta, tipo, fuente, verificado, nota}
    capas: [],           // [{medio, espesor_mm}]
    freq: 500,
    ppw: 14,
    // estructura: `alto` 0 = capas infinitas (interfaz plana). Con alto > 0 es una pieza finita que se
    // arrastra sobre el mapa; ahí la onda rodea la pieza y difracta en los bordes.
    estructura: { alto: 0, x: 12, y: 0 },
    excitacion: 'plana',      // 'plana' | 'fuentes' (el arreglo de la vista BEM · fuentes)
    resultado: null,
    resolviendo: false,
  };
  window.CAD = CAD;

  /* ---------------------------------------------------------------- capas */
  function filaCapa(capa, i) {
    const fila = document.createElement('div');
    fila.className = 'fila-capa';

    const sel = document.createElement('select');
    for (const [nombre, medio] of Object.entries(CAD.medios)) {
      const op = document.createElement('option');
      op.value = nombre;
      op.textContent = medio.etiqueta + (medio.tipo === 'solido' ? ' (sólido)' : '');
      if (nombre === capa.medio) op.selected = true;
      sel.appendChild(op);
    }
    sel.addEventListener('change', () => { capa.medio = sel.value; actualizarAviso(); });

    const esp = document.createElement('input');
    esp.type = 'number';
    esp.min = '0.1';
    esp.step = '0.1';
    esp.value = capa.espesor_mm;
    esp.addEventListener('input', () => {
      capa.espesor_mm = Math.max(0.1, parseFloat(esp.value) || 0.1);
    });

    const quitar = document.createElement('button');
    quitar.className = 'btn mini';
    quitar.textContent = '×';
    quitar.title = 'quitar esta capa';
    quitar.addEventListener('click', () => {
      CAD.capas.splice(i, 1);
      pintarCapas();
    });

    fila.append(sel, esp, quitar);
    return fila;
  }

  function pintarCapas() {
    const cont = $('cadenaCapas');
    cont.innerHTML = '';
    CAD.capas.forEach((capa, i) => cont.appendChild(filaCapa(capa, i)));
    actualizarAviso();
  }

  function actualizarAviso() {
    const aviso = $('cadenaAviso');
    if (!aviso) return;
    const solidos = CAD.capas.filter((c) => CAD.medios[c.medio] && CAD.medios[c.medio].tipo === 'solido');
    if (solidos.length) {
      aviso.hidden = false;
      aviso.textContent =
        `${solidos.map((c) => CAD.medios[c.medio].etiqueta).join(', ')} entra como sólido elástico `
        + '(tiene cizalla): no hay coeficiente exacto de capas fluidas contra el cual comparar, así que '
        + 'se reporta el transmitido medido sin error.';
    } else {
      aviso.hidden = true;
    }
  }

  /* -------------------------------------------------------------- resolver */
  function estado(texto, esError) {
    const el = $('cadenaEstado');
    el.textContent = texto;
    el.style.color = esError ? '#f85149' : '';
  }

  async function resolver() {
    if (CAD.resolviendo) return;
    if (!CAD.capas.length) { estado('agregá al menos una capa.', true); return; }
    CAD.resolviendo = true;
    $('cadenaResolver').disabled = true;
    estado('resolviendo… (malla + acople de N regiones)');

    const cuerpo = {
      capas: CAD.capas.map((c) => ({ medio: c.medio, espesor_mm: c.espesor_mm })),
      frecuencia_hz: CAD.freq * 1e3,
      ppw: CAD.ppw,
      nx: 160,
      ny: 100,
      // con "fuentes" se manda el arreglo de la vista de campo tal cual está en pantalla
      // `campoFuentesActuales` devuelve el arreglo del modo activo (manual o arreglo): leer sólo
      // `S.fuentes` daba cero en modo arreglo y la corrida caía en onda plana sin avisar
      fuentes: CAD.excitacion === 'fuentes' && typeof window.campoFuentesActuales === 'function'
        ? window.campoFuentesActuales() : null,
      estructura_alto_mm: CAD.estructura.alto,
      estructura_x_mm: CAD.estructura.x,
      estructura_y_mm: CAD.estructura.y,
    };
    try {
      const alta = await fetch('/api/jobs/cadena', {
        method: 'POST', headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(cuerpo),
      }).then((r) => r.json());
      if (!alta.ok) throw new Error(alta.error || 'no se pudo encolar');
      const id = alta.job.id;
      let job = null;
      for (let i = 0; i < 240; i += 1) {
        await new Promise((r) => setTimeout(r, 600));
        job = await fetch(`/api/jobs/${id}`).then((r) => r.json());
        if (['listo', 'listo_con_aviso', 'error', 'cancelado'].includes(job.estado)) break;
      }
      if (!job || job.estado === 'error') throw new Error((job && job.error) || 'la corrida falló');
      if (job.estado === 'cancelado') throw new Error('cancelado');
      CAD.resultado = job.resultado;
      pintar(job.resultado);
      estado(`listo · ${job.resultado.metricas.nodos} nodos · `
        + `${job.resultado.tiempos.solve_s} s de solve`);
      if (window.vistaGuardar) {
        window.vistaGuardar(CLAVE, {
          corrida: CAD.resultado, capas: CAD.capas, freq: CAD.freq, ppw: CAD.ppw,
          estructura: CAD.estructura,
        });
      }
    } catch (err) {
      estado(String(err.message || err), true);
    } finally {
      CAD.resolviendo = false;
      $('cadenaResolver').disabled = false;
    }
  }

  async function resolverSTL(archivo) {
    const estadoSTL = $('stlEstado');
    estadoSTL.style.color = '';
    estadoSTL.textContent = `leyendo ${archivo.name}…`;
    const buffer = await archivo.arrayBuffer();
    let bin = '';
    const bytes = new Uint8Array(buffer);
    for (let i = 0; i < bytes.length; i += 1) bin += String.fromCharCode(bytes[i]);
    const datos_b64 = btoa(bin);

    const cuerpo = {
      datos_b64,
      nombre: archivo.name,
      material: $('stlMaterial').value,
      frecuencia_hz: CAD.freq * 1e3,
      ppw: CAD.ppw,
      plano: { eje: 'z', valor_mm: parseFloat($('stlPlano').value) },
      nx: 160, ny: 110,
    };
    estadoSTL.textContent = `resolviendo la sección de ${archivo.name}…`;
    try {
      const alta = await fetch('/api/jobs/stl', {
        method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(cuerpo),
      }).then((r) => r.json());
      if (!alta.ok) throw new Error(alta.error || 'no se pudo encolar');
      let job = null;
      for (let i = 0; i < 240; i += 1) {
        await new Promise((r) => setTimeout(r, 600));
        job = await fetch(`/api/jobs/${alta.job.id}`).then((r) => r.json());
        if (['listo', 'listo_con_aviso', 'error', 'cancelado'].includes(job.estado)) break;
      }
      if (!job || job.estado === 'error') throw new Error((job && job.error) || 'la corrida falló');
      pintarSTL(job.resultado);
      estadoSTL.textContent = `${archivo.name} · sección ${job.resultado.stl.seccion.ancho_mm}×`
        + `${job.resultado.stl.seccion.alto_mm} mm · ${job.resultado.metricas.nodos} nodos · `
        + `${job.resultado.tiempos.solve_s} s de solve`;
    } catch (err) {
      estadoSTL.style.color = '#f85149';
      estadoSTL.textContent = `${archivo.name}: ${err.message || err}`;
    }
  }

  function pintarSTL(r) {
    const m = r.metricas;
    const cont = $('cadenaMetricas');
    cont.hidden = false;
    cont.innerHTML = `
      <div class="chip"><span>transmitido detrás</span><b>${m.transmitido === null ? '—' : m.transmitido.toFixed(4)}</b></div>
      <div class="chip"><span>material</span><b>${r.material.etiqueta}</b></div>
      <div class="chip"><span>área de la sección</span><b>${r.stl.seccion.area_mm2} mm²</b></div>
      <div class="chip"><span>malla</span><b>${m.tamano_malla_mm} mm</b></div>
      <div class="chip"><span>dominio</span><b>${m.dominio_mm[0]}×${m.dominio_mm[1]} mm</b></div>`;

    const g = r.grilla;
    const shapes = (g.contorno || []).map((c) => ({
      type: 'path', path: c.map((p, i) => `${i ? 'L' : 'M'}${p[0]},${p[1]}`).join(' ') + ' Z',
      line: { color: '#58a6ff', width: 2 }, fillcolor: 'rgba(88,166,255,.12)', layer: 'above',
    }));
    Plotly.react('mapaCadena', [{
      type: 'heatmap', x: g.x_mm, y: g.y_mm, z: g.db, zmin: -30, zmax: 0,
      colorscale: 'Viridis', colorbar: { title: 'dB', thickness: 12 },
      hovertemplate: 'x %{x:.1f} mm<br>y %{y:.1f} mm<br>%{z:.1f} dB<extra></extra>',
    }], {
      margin: { l: 50, r: 10, t: 30, b: 40 }, height: 340,
      title: { text: `|p| (dB) · sección de ${r.stl.nombre}`, font: { size: 12 } },
      xaxis: { title: 'x [mm]', range: [Math.min(...g.x_mm), Math.max(...g.x_mm)] },
      yaxis: { title: 'y [mm]', range: [Math.min(...g.y_mm), Math.max(...g.y_mm)] },
      shapes,
    }, { displayModeBar: false, responsive: true });
  }

  /* ---------------------------------------------------------------- pintar */
  function pintar(r) {
    const m = r.metricas;
    const cont = $('cadenaMetricas');
    cont.hidden = false;
    const tMed = m.t_medido === null ? '—' : m.t_medido.toFixed(4);
    const tTeo = m.t_teorico === null ? 'sin referencia' : m.t_teorico.toFixed(4);
    const err = m.error_pct === null ? '—' : m.error_pct.toFixed(2) + ' %';
    const pico = m.pico_detras === null || m.pico_detras === undefined
      ? '—' : m.pico_detras.toExponential(2);
    const chips = m.hay_fuentes
      ? `<div class="chip"><span>excitación</span><b>${m.n_fuentes} fuentes</b></div>
         <div class="chip"><span>pico detrás</span><b>${pico}</b></div>
         <div class="chip"><span>espesor total</span><b>${m.espesor_total_mm} mm</b></div>
         <div class="chip"><span>λ en agua</span><b>${m.lambda_agua_mm} mm</b></div>
         <div class="chip"><span>malla</span><b>${m.ppw_efectivo} elem./λ</b></div>`
      : `<div class="chip"><span>transmitido medido</span><b>${tMed}</b></div>
         <div class="chip"><span>exacto (capas fluidas)</span><b>${tTeo}</b></div>
         <div class="chip"><span>error</span><b>${err}</b></div>
         <div class="chip"><span>espesor total</span><b>${m.espesor_total_mm} mm</b></div>
         <div class="chip"><span>λ en agua</span><b>${m.lambda_agua_mm} mm</b></div>`;
    cont.innerHTML = chips;
    // el aviso de malla recortada no puede pasar inadvertido: cambia cuánto resuelve el mapa
    if (m.aviso_malla) {
      estado(`${m.aviso_malla}.`, false);
    }

    const g = r.grilla;
    const e = r.estructura || { x_mm: 0, y_mm: 0, ancho_mm: 0, alto_mm: null, finita: false };
    // índice del shape de la estructura en el layout, para reconocerlo en `plotly_relayout`
    CAD.iEstructura = e.finita ? (r.capas || []).length * 2 : -1;

    // Con la estructura finita, las líneas de interfaz se recortan a su alto: si atraviesan todo el
    // dominio parece una capa infinita y no se entiende que la pieza se puede rodear.
    const yTop = e.finita ? e.y_mm + (e.alto_mm || 0) / 2 : g.y_mm[g.y_mm.length - 1];
    const yBot = e.finita ? e.y_mm - (e.alto_mm || 0) / 2 : g.y_mm[0];

    const shapes = [
      ...(r.capas || []).flatMap((c) => [c.desde_mm, c.hasta_mm].map((x) => ({
        type: 'line', x0: x, x1: x, y0: yBot, y1: yTop,
        line: { color: '#ffffff', width: 2, dash: 'solid' },
      }))),
    ];
    if (e.finita) {
      // El rectángulo se dibuja con relleno oscuro y borde azul: con relleno translúcido se confundía
      // con el mapa de presión (el sólido no tiene gdl de presión y salía negro, igual que un mínimo
      // del colormap) y no se leía como objeto.
      shapes.push({
        type: 'rect', name: 'estructura',
        x0: e.x_mm, x1: e.x_mm + e.ancho_mm,
        y0: e.y_mm - (e.alto_mm || 0) / 2, y1: e.y_mm + (e.alto_mm || 0) / 2,
        line: { color: '#58a6ff', width: 2.5, dash: 'dash' },
        fillcolor: 'rgba(13,17,23,.62)', layer: 'above',
      });
    }

    Plotly.react('mapaCadena', [{
      type: 'heatmap', x: g.x_mm, y: g.y_mm, z: g.db, zmin: -40, zmax: 0,
      colorscale: 'Viridis', colorbar: { title: 'dB', thickness: 12 },
      hovertemplate: 'x %{x:.1f} mm<br>y %{y:.1f} mm<br>%{z:.1f} dB<extra></extra>',
    }], {
      margin: { l: 50, r: 10, t: 30, b: 40 }, height: 320,
      title: { text: `|p| relativo (dB) · ${(m.frecuencia_hz / 1e3).toFixed(0)} kHz`, font: { size: 12 } },
      // Con `scaleanchor` y un div ancho, Plotly estiraba el eje x a −40…80 y recortaba el y: las
      // líneas de interfaz y las etiquetas de capa quedaban FUERA del rango visible. Con rangos fijos
      // y sin anclaje se ve todo el dominio; la proporción física la da el propio dominio (36×48 mm).
      xaxis: { title: 'x [mm]', range: [0, Math.max(...g.x_mm)] },
      yaxis: { title: 'y [mm]', range: [Math.min(...g.y_mm), Math.max(...g.y_mm)] },
      // las interfaces de capa tienen que verse: si no, el patrón de interferencia se confunde con
      // las capas y no se sabe dónde empieza qué
      shapes,
      // las etiquetas van un poco adentro del borde: en y = ymax quedaban en el margen y no se leían
      annotations: (r.capas || []).map((c) => ({
        x: (c.desde_mm + c.hasta_mm) / 2, y: Math.max(...g.y_mm) * 0.92, text: c.etiqueta,
        showarrow: false, yanchor: 'top', font: { size: 10, color: '#ffffff' },
        bgcolor: 'rgba(13,17,23,.55)', borderpad: 1,
      })),
    }, {
      displayModeBar: false, responsive: true,
      // con esto el rectángulo de la estructura se puede arrastrar directamente en el mapa
      edits: { shapePosition: e.finita },
    });

    // Arrastrar el rectángulo emite `plotly_relayout` con claves `shapes[N].x0/x1/y0/y1` (N cambia con
    // la cantidad de capas: es el índice del rect, después de las 2 líneas por interfaz). Se registra
    // acá, después del react, y se quita el anterior para no acumular handlers.
    const gdMapa = $('mapaCadena');
    // sin esto no hay ninguna señal de que el rectángulo se puede arrastrar: se lee como un obstáculo
    if (gdMapa) gdMapa.style.cursor = e.finita ? 'grab' : 'default';
    if (gdMapa && typeof gdMapa.on === 'function') {
      if (gdMapa.__cadHandler) gdMapa.removeListener('plotly_relayout', gdMapa.__cadHandler);
      gdMapa.__cadHandler = (ev) => {
        const tocaShapes = Object.keys(ev).some((k) => k === 'shapes' || k.startsWith('shapes['));
        if (!tocaShapes || CAD.iEstructura < 0) return;
        const sh = (gdMapa.layout.shapes || [])[CAD.iEstructura];
        if (!sh || sh.name !== 'estructura') return;
        // `estructura_x_mm` es el borde izquierdo (el backend apila las capas desde ahí) y el alto no
        // cambia: arrastrar mueve, no redimensiona
        const nx = Math.round(sh.x0 * 10) / 10;
        const ny = Math.round(((sh.y0 + sh.y1) / 2) * 10) / 10;
        if (nx === CAD.estructura.x && ny === CAD.estructura.y) return;
        CAD.estructura.x = nx;
        CAD.estructura.y = ny;
        if (window.cadSincronizarEstructura) window.cadSincronizarEstructura();
        estado(`estructura en x = ${nx.toFixed(1)} mm, y = ${ny.toFixed(1)} mm · resolviendo…`);
        resolver();
      };
      gdMapa.on('plotly_relayout', gdMapa.__cadHandler);
    }

    Plotly.react('axialCadena', [{
      x: r.corte_axial.x_mm, y: r.corte_axial.abs_p, type: 'scatter', mode: 'lines',
      line: { color: '#58a6ff' },
    }], {
      margin: { l: 48, r: 10, t: 28, b: 38 }, height: 210,
      title: { text: 'corte axial · |p| en el eje', font: { size: 11 } },
      xaxis: { title: 'x [mm]' }, yaxis: { title: '|p|' },
    }, { displayModeBar: false, responsive: true });

    // Las capas se dibujan como rectángulos en su posición real en x: una barra por capa con su
    // etiqueta en y no decía a qué distancia estaba cada una.
    const haySolido = (r.capas || []).some((c) => CAD.medios[c.medio] && CAD.medios[c.medio].tipo === 'solido');
    Plotly.react('capasCadena', [{
      x: [0], y: [0], type: 'scatter', mode: 'markers', marker: { size: 0.01, opacity: 0 },
      showlegend: false, hoverinfo: 'skip',
    }], {
      margin: { l: 46, r: 10, t: 28, b: 38 }, height: 210,
      title: {
        text: (haySolido ? 'capas · ámbar = sólido elástico (con cizalla)' : 'capas · todas fluidas')
          + (r.estructura && r.estructura.finita ? ` · pieza finita de ${r.estructura.alto_mm} mm` : ''),
        font: { size: 11 },
      },
      xaxis: { title: 'x [mm]', range: [0, Math.max(...g.x_mm)] },
      yaxis: { visible: false, range: [0, 1] },
      shapes: [
        // el agua de entrada y de salida también son parte de lo que se simula: sin ellas las capas
        // parecían empezar en el medio de la nada
        { type: 'rect', x0: 0, x1: (r.capas[0] || {}).desde_mm || 0, y0: 0.25, y1: 0.75,
          line: { width: 0 }, fillcolor: '#30363d' },
        { type: 'rect', x0: (r.capas[r.capas.length - 1] || {}).hasta_mm || 0,
          x1: Math.max(...g.x_mm), y0: 0.25, y1: 0.75, line: { width: 0 }, fillcolor: '#30363d' },
        ...(r.capas || []).map((c) => ({
          type: 'rect', x0: c.desde_mm, x1: c.hasta_mm, y0: 0.25, y1: 0.75,
          line: { color: '#0d1117', width: 1 },
          fillcolor: CAD.medios[c.medio] && CAD.medios[c.medio].tipo === 'solido' ? '#d29922' : '#3fb950',
        })),
      ],
      annotations: [
        ...(r.capas || []).map((c) => ({
          x: (c.desde_mm + c.hasta_mm) / 2, y: 0.5, text: `${c.etiqueta}\n${c.espesor_mm} mm`,
          showarrow: false, font: { size: 9, color: '#0d1117' },
        })),
        { x: ((r.capas[0] || {}).desde_mm || 0) / 2, y: 0.5, text: 'agua',
          showarrow: false, font: { size: 9, color: '#8f9aa6' } },
        { x: ((r.capas[r.capas.length - 1] || {}).hasta_mm + Math.max(...g.x_mm)) / 2, y: 0.5,
          text: 'agua', showarrow: false, font: { size: 9, color: '#8f9aa6' } },
      ],
    }, { displayModeBar: false, responsive: true });
  }

  /* ----------------------------------------------------------------- init */
  async function init() {
    const op = await fetch('/api/cadena/opciones').then((r) => r.json());
    CAD.medios = Object.fromEntries(op.medios.map((m) => [m.nombre, m]));
    const avisoTop = document.querySelector('#vistaCadena .que-ves');
    if (avisoTop && op.aviso) avisoTop.title = op.aviso;

    const guardado = window.vistaCargar ? window.vistaCargar(CLAVE) : null;
    const previo = guardado && guardado.resultado && guardado.resultado.corrida;
    if (previo && (guardado.resultado.capas || []).length) {
      CAD.capas = guardado.resultado.capas;
      CAD.freq = guardado.resultado.freq || 500;
      CAD.ppw = guardado.resultado.ppw || 14;
      if (guardado.resultado.estructura) CAD.estructura = guardado.resultado.estructura;
      pintarCapas();
      if (window.cadSincronizarEstructura) window.cadSincronizarEstructura();
      $('cadFreq').value = CAD.freq; $('cadFreqVal').textContent = CAD.freq + ' kHz';
      $('cadPpw').value = CAD.ppw; $('cadPpwVal').textContent = CAD.ppw;
      CAD.resultado = previo;
      pintar(previo);
      estado('restaurado de la sesión anterior (no se volvió a resolver).');
    } else {
      CAD.capas = [
        { medio: 'piel', espesor_mm: 2 },
        { medio: 'hueso_cortical', espesor_mm: 5 },
        { medio: 'perilinfa', espesor_mm: 5 },
      ];
      pintarCapas();
    }

    $('cadenaAdd').addEventListener('click', () => {
      const ultimo = CAD.capas[CAD.capas.length - 1];
      CAD.capas.push({ medio: ultimo ? ultimo.medio : 'agua', espesor_mm: 1 });
      pintarCapas();
    });
    $('cadenaEjemplo').addEventListener('click', () => {
      CAD.capas = [
        { medio: 'piel', espesor_mm: 2 },
        { medio: 'grasa', espesor_mm: 3 },
        { medio: 'hueso_cortical', espesor_mm: 5 },
        { medio: 'perilinfa', espesor_mm: 5 },
      ];
      pintarCapas();
    });
    $('cadFreq').addEventListener('input', (e) => {
      CAD.freq = parseInt(e.target.value, 10);
      $('cadFreqVal').textContent = CAD.freq + ' kHz';
      if (window.syncFrecuencia) window.syncFrecuencia(CAD.freq, 'cadena');
    });
    $('cadPpw').addEventListener('input', (e) => {
      CAD.ppw = parseInt(e.target.value, 10);
      $('cadPpwVal').textContent = CAD.ppw;
    });
    $('cadenaResolver').addEventListener('click', resolver);

    $('cadExcitacion').addEventListener('change', (e) => {
      CAD.excitacion = e.target.value;
      const conFuentes = CAD.excitacion === 'fuentes';
      $('cadExcitacionInfo').textContent = conFuentes
        ? 'Se usan las fuentes de la vista BEM · fuentes (posición, ancho, amplitud y fase de cada '
          + 'elemento). El campo enfoca y ya no hay coeficiente exacto contra el cual comparar.'
        : 'Con onda plana la cadena se compara contra el coeficiente exacto de capas.';
      // con fuentes, la referencia analítica no aplica: se avisa antes de correr, no después
      if (conFuentes) {
        const n = (window.campoFuentesActuales ? window.campoFuentesActuales().length : 0);
        estado(n ? `${n} fuentes del arreglo listas. Resolvé para ver el campo.`
                 : 'la vista BEM · fuentes no tiene fuentes cargadas.');
      }
    });

    // ---- STL: se simula su sección (el modelo es 2D) -----------------------
    (async () => {
      try {
        const op = await fetch('/api/stl/opciones').then((r) => r.json());
        const sel = $('stlMaterial');
        sel.innerHTML = '';
        for (const s of op.solidos) {
          const o = document.createElement('option');
          o.value = s.id;
          o.textContent = s.etiqueta + (s.verificado ? '' : ' (VERIFY)');
          if (s.id === 'hueso_cortical') o.selected = true;
          sel.appendChild(o);
        }
      } catch (e) {
        console.warn('no se pudieron cargar los materiales de STL:', e && e.message);
      }
    })();
    $('stlPlano').addEventListener('input', (e) => {
      $('stlPlanoVal').textContent = parseFloat(e.target.value).toFixed(1);
    });
    $('stlCargarBtn').addEventListener('click', () => $('stlArchivo').click());
    $('stlEjemplo').addEventListener('click', async () => {
      try {
        const blob = await fetch('/api/stl/ejemplo').then((r) => r.blob());
        resolverSTL(new File([blob], 'cilindro_ejemplo.stl', { type: 'model/stl' }));
      } catch (e) {
        const st = $('stlEstado');
        st.style.color = '#f85149';
        st.textContent = 'no se pudo traer el ejemplo: ' + (e.message || e);
      }
    });
    $('stlArchivo').addEventListener('change', (e) => {
      const f = e.target.files && e.target.files[0];
      e.target.value = '';
      if (f) resolverSTL(f);
    });

    // ---- estructura: alto, x, y y arrastre en el mapa ----------------------
    function sincronizarControlesEstructura() {
      $('cadAlto').value = CAD.estructura.alto;
      $('cadAltoVal').textContent = CAD.estructura.alto ? CAD.estructura.alto + ' mm' : 'todo el dominio';
      $('cadX').value = CAD.estructura.x;
      $('cadXVal').textContent = CAD.estructura.x.toFixed(1);
      $('cadY').value = CAD.estructura.y;
      $('cadYVal').textContent = CAD.estructura.y.toFixed(1);
    }
    window.cadSincronizarEstructura = sincronizarControlesEstructura;

    $('cadAlto').addEventListener('input', (e) => {
      CAD.estructura.alto = parseInt(e.target.value, 10);
      sincronizarControlesEstructura();
    });
    $('cadX').addEventListener('input', (e) => {
      CAD.estructura.x = parseFloat(e.target.value);
      $('cadXVal').textContent = CAD.estructura.x.toFixed(1);
    });
    $('cadY').addEventListener('input', (e) => {
      CAD.estructura.y = parseFloat(e.target.value);
      $('cadYVal').textContent = CAD.estructura.y.toFixed(1);
    });
    sincronizarControlesEstructura();

    // El listener del arrastre se registra en `pintar()`: acá el div todavía no es un gráfico de Plotly
    // y no tiene `.on`, así que quedaba sin registrar y el arrastre no recalculaba nada.

    window.cadenaPintar = () => { if (CAD.resultado) pintar(CAD.resultado); };
  }

  // Un fallo dentro de `init` es un rechazo de promesa: sin este catch la pestaña quedaba vacía y sin
  // decir por qué (que es exactamente lo que pasó con la clave de la persistencia).
  function arrancar() {
    init().catch((e) => {
      estado('no se pudo preparar la vista: ' + (e && e.message ? e.message : e), true);
      console.error('[cadena] init falló:', e);
    });
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', arrancar);
  } else {
    arrancar();
  }
})();
