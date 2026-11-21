/* Vista "criterios": el registro de validación del servidor, con lo medido en esta sesión al lado.
 *
 * Por qué es una vista y no un tooltip
 * ------------------------------------
 * El sello de cada corrida es la pieza que decide si un número se reporta o no, y hasta ahora sólo
 * se veía de a uno, dentro de la vista que lo produjo. Acá se ve el registro completo: qué se mide,
 * con qué tolerancia y con qué fuente, y cuáles de esos criterios tienen un número medido en esta
 * sesión. Un criterio sin medición no es un aprobado: dice "sin medir".
 *
 * Los números medidos se leen del localStorage (los mismos resultados que guardan las vistas). La
 * UI no mide nada acá: sólo cruza el registro con lo que ya se resolvió.
 */

const CRIT = { criterios: [], listo: false };

function crit$(id) { return document.getElementById(id); }

window.critIniciar = function critIniciar() {
  critRefrescar();
  if (!CRIT.listo) {
    CRIT.listo = true;
    // si alguna vista resuelve mientras esta pestaña está abierta, al volver se ve el número nuevo
    window.addEventListener('focus', critRefrescar);
  }
};

/* Junta lo medido en esta sesión: recorre los resultados guardados por las vistas y se queda con el
 * último sello de cada criterio. `usim.vista.<nombre>` lo escribe app.js en `vistaGuardar`. */
function critMediciones() {
  const porCriterio = {};
  for (let i = 0; i < localStorage.length; i += 1) {
    const clave = localStorage.key(i);
    if (!clave || !clave.startsWith('usim.vista.')) continue;
    const vista = clave.slice('usim.vista.'.length);
    let datos;
    try { datos = JSON.parse(localStorage.getItem(clave)); } catch { continue; }
    const sello = datos && datos.resultado && datos.resultado.sello;
    if (!sello || !Array.isArray(sello.criterios)) continue;
    for (const c of sello.criterios) {
      if (!c || !c.id) continue;
      const previo = porCriterio[c.id];
      // se queda el más reciente
      if (!previo || (datos.fecha && datos.fecha > previo.fecha)) {
        porCriterio[c.id] = {
          medido: c.medido, estado: c.estado, margen: c.margen,
          unidad: c.unidad, vista, fecha: datos.fecha,
        };
      }
    }
  }
  return porCriterio;
}

function critEstadoClase(estado) {
  if (estado === 'validado') return 'ok';
  if (estado === 'fuera') return 'mal';
  if (estado === 'sin_referencia') return 'ojo';
  return 'gris';
}

/* Sólo las unidades que SON un porcentaje se multiplican por 100. `amplitud` es un valor absoluto
 * adimensional (el sello lo dice así: "0.001945 amplitud"), y mostrarlo como 0.19 % inventaría una
 * unidad que la referencia no usa. */
const UNIDADES_PORCENTUALES = new Set(['máx. relativo', 'L2 relativo', '% relativo',
  'de la potencia interceptada']);

function critNumero(v, unidad) {
  if (v === null || v === undefined) return '—';
  if (unidad === 'amplitud' || unidad === 'adimensional') {
    return Math.abs(v) !== 0 && Math.abs(v) < 1e-3 ? v.toExponential(2) : v.toPrecision(4);
  }
  if (UNIDADES_PORCENTUALES.has(unidad)) return `${(100 * v).toFixed(2)} %`;
  const a = Math.abs(v);
  return a !== 0 && (a < 1e-3 || a >= 1e5) ? v.toExponential(2) : v.toPrecision(4);
}

function critFecha(iso) {
  if (!iso) return '';
  try { return new Date(iso).toLocaleString('es', { day: '2-digit', month: '2-digit', hour: '2-digit', minute: '2-digit' }); }
  catch { return iso; }
}

async function critRefrescar() {
  const caja = crit$('tablaCriterios');
  if (!caja) return;
  try {
    if (!CRIT.criterios.length) {
      // `/api/validacion` devuelve `{criterios: [...]}` sin campo `ok` (las otras rutas sí lo traen).
      // Comprobar `r.ok` acá lanzaba siempre y la tabla quedaba vacía sin decir por qué.
      const r = await fetch('/api/validacion').then((x) => x.json());
      if (!r || !Array.isArray(r.criterios)) throw new Error(r && r.error ? r.error : 'sin datos');
      CRIT.criterios = r.criterios;
    }
  } catch (e) {
    caja.innerHTML = `<p class="hint">no se pudo leer el registro: ${e.message}</p>`;
    return;
  }

  const medidos = critMediciones();
  let validados = 0; let fuera = 0; let sinMedir = 0; let sinRef = 0;

  const filas = CRIT.criterios.map((c) => {
    const m = medidos[c.id];
    const estado = m ? m.estado : 'sin medir';
    if (estado === 'validado') validados += 1;
    else if (estado === 'fuera') fuera += 1;
    else if (estado === 'sin_referencia') sinRef += 1;
    else sinMedir += 1;
    const tol = c.modo === 'rel' || UNIDADES_PORCENTUALES.has(c.unidad)
      ? `≤ ${(100 * c.tolerancia).toFixed(1)} %`
      : `≤ ${c.tolerancia} ${c.unidad}`;
    return `
      <tr class="crit-fila crit-${critEstadoClase(estado)}">
        <td class="crit-id">${c.id}</td>
        <td class="crit-estado">${estado === 'sin medir' ? 'sin medir' : estado}</td>
        <td class="crit-medido">${m ? critNumero(m.medido, c.unidad) : '—'}</td>
        <td class="crit-tol">${tol}</td>
        <td class="crit-desc">${c.descripcion}</td>
        <td class="crit-vista">${m ? `${m.vista}<br><span class="crit-fecha">${critFecha(m.fecha)}</span>` : '—'}</td>
        <td class="crit-fuente" title="${(c.fuente || '').replace(/"/g, '&quot;')}">${(c.fuente || '').slice(0, 60)}${(c.fuente || '').length > 60 ? '…' : ''}</td>
      </tr>`;
  }).join('');

  caja.innerHTML = `
    <table>
      <thead><tr>
        <th>criterio</th><th>estado</th><th>medido</th><th>tolerancia</th>
        <th>qué se mide</th><th>vista</th><th>fuente</th>
      </tr></thead>
      <tbody>${filas}</tbody>
    </table>`;

  crit$('criteriosResumen').textContent =
    `${CRIT.criterios.length} criterios en el registro · ${validados} validados con lo medido en esta`
    + ` sesión · ${fuera} fuera · ${sinRef} sin referencia · ${sinMedir} sin medir todavía`;
}
