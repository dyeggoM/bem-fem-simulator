/* Cargar el arreglo de fuentes planeado, y el catálogo de transductores (U10).
 *
 * Dos cosas distintas que se piden juntas:
 *   1. un **transductor del catálogo** (equipo comercial o de referencia) del que salen frecuencia,
 *      cantidad de elementos, apertura y foco;
 *   2. un **arreglo propio** en CSV o JSON, con la posición, el ancho, la amplitud y la fase de cada
 *      elemento.
 *
 * Todo lo que se carga pasa por el servidor antes de entrar al solver: un archivo con dos elementos
 * pisados o con las distancias en metros daría un campo que parece válido y no lo es. Si el servidor
 * rechaza el archivo, se dice por qué y no se toca lo que está en pantalla.
 */
(function () {
  'use strict';

  const $ = (id) => document.getElementById(id);
  let catalogo = [];

  function aviso(texto, esError) {
    const el = $('arregloAviso');
    if (!el) return;
    el.textContent = texto || '';
    el.style.color = esError ? '#f85149' : '';
  }

  /* ------------------------------------------------------------ catálogo */
  async function cargarCatalogo() {
    // el selector de material de la vista de campo es `#mat`: los lóbulos de rejilla y el paso en λ se
    // calculan con la longitud de onda de ESE medio, no con la del agua
    const material = ($('mat') && $('mat').value) || 'agua';
    try {
      const r = await fetch(`/api/transductores?material=${encodeURIComponent(material)}`).then((x) => x.json());
      if (!r.ok) throw new Error(r.error || 'no se pudo leer el catálogo');
      catalogo = r.transductores || [];
      const sel = $('selTransductor');
      sel.innerHTML = '<option value="">— elegir —</option>';
      for (const t of catalogo) {
        const op = document.createElement('option');
        op.value = t.key;
        op.textContent = `${t.name} · ${Math.round(t.frequency_hz / 1e3)} kHz · ${t.n_elements} elem.`
          + (t.verified ? '' : ' (VERIFY)');
        sel.appendChild(op);
      }
      if (!catalogo.length) aviso('el catálogo está vacío: copiá la plantilla de data/transducers.toml', true);
    } catch (e) {
      aviso('catálogo: ' + (e.message || e), true);
    }
  }

  function describir(t) {
    const partes = [
      `${Math.round(t.frequency_hz / 1e3)} kHz`,
      `${t.n_elements} elemento${t.n_elements === 1 ? '' : 's'}`,
      `apertura ${t.width_mm} mm`,
      t.pitch_over_lambda ? `paso ${t.pitch_over_lambda.toFixed(3)} λ` : null,
      t.f_number ? `f/${t.f_number.toFixed(2)}` : null,
      t.grating_lobe_deg !== null && t.grating_lobe_deg !== undefined
        ? `lóbulo de rejilla ${t.grating_lobe_deg.toFixed(1)}°` : 'sin lóbulo de rejilla',
    ].filter(Boolean);
    return partes.join(' · ') + (t.verified ? '' : ' · VERIFY: sin datasheet verificado');
  }

  function aplicarTransductor() {
    const key = $('selTransductor').value;
    if (!key) { $('transductorInfo').textContent = ''; return; }
    const t = catalogo.find((x) => x.key === key);
    if (!t) return;
    $('transductorInfo').textContent = describir(t);
    if (typeof window.campoAplicarTransductor === 'function') {
      window.campoAplicarTransductor(t);
      aviso(`transductor aplicado: ${t.n_elements} elementos, apertura ${t.width_mm} mm.`);
    }
  }

  /* ------------------------------------------------------- carga/exportación */
  async function cargarArchivo(archivo) {
    if (!archivo) return;
    const texto = await archivo.text();
    try {
      const r = await fetch('/api/arreglo/cargar', {
        method: 'POST', headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ texto, formato: archivo.name.toLowerCase().endsWith('.json') ? 'json' : 'csv' }),
      }).then((x) => x.json());
      if (!r.ok) throw new Error(r.error || 'el servidor rechazó el archivo');
      const fuentes = r.fuentes.map((f) => ({
        y: f.y_mm * 1e-3, width: f.ancho_mm * 1e-3, amp: f.amp, phase: f.fase_rad,
      }));
      window.campoAplicarFuentes(fuentes);
      aviso(`${r.n} fuentes cargadas de ${archivo.name} · apertura ${r.apertura_mm.toFixed(1)} mm.`);
    } catch (e) {
      // no se toca lo que está en pantalla: el archivo no pasó la validación
      aviso(`${archivo.name}: ${e.message || e}`, true);
    }
  }

  async function exportar() {
    const fuentes = window.campoLeerFuentes();
    if (!fuentes.length) { aviso('no hay fuentes para exportar', true); return; }
    try {
      const r = await fetch('/api/arreglo/exportar', {
        method: 'POST', headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ fuentes }),
      }).then((x) => x.json());
      if (!r.ok) throw new Error(r.error || 'no se pudo exportar');
      const blob = new Blob([r.texto], { type: 'text/csv' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = `arreglo-${r.n}fuentes.csv`;
      a.click();
      URL.revokeObjectURL(url);
      aviso(`exportadas ${r.n} fuentes a CSV.`);
    } catch (e) {
      aviso('exportar: ' + (e.message || e), true);
    }
  }

  function init() {
    if (!$('selTransductor')) return;
    $('selTransductor').addEventListener('change', aplicarTransductor);
    $('btnCargarArreglo').addEventListener('click', () => $('archivoArreglo').click());
    $('archivoArreglo').addEventListener('change', (e) => {
      cargarArchivo(e.target.files && e.target.files[0]);
      e.target.value = '';
    });
    $('btnExportarArreglo').addEventListener('click', exportar);
    // el catálogo depende del material: los lóbulos de rejilla se calculan con su longitud de onda
    const mat = $('mat');
    if (mat) mat.addEventListener('change', cargarCatalogo);
    cargarCatalogo();
  }

  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init);
  else init();
})();
