"""Genera ``out/report.html``: una página de estado autocontenida del simulador.

No necesita servidor ni dependencias fuera de la biblioteca estándar: lee los CSV de ``out/`` y los
ADR de ``docs/decisions.md``, embebe las figuras como data-URI y escribe un único HTML que se abre
con doble clic y se puede publicar tal cual (GitHub Pages, por ejemplo).

Uso::

    .venv/bin/python scripts/report.py

Qué NO hace: no mide nada. Si un número no está en un CSV o en la tabla ``FASES`` de este archivo,
no aparece. Las cifras de ``FASES`` son las que verifican los tests, y cada fila cita el ADR que las
respalda: si cambia un número, cambia también su ADR.
"""

from __future__ import annotations

import base64
import csv
import datetime as dt
import html
import re
import subprocess
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
OUT = RAIZ / "out"

#: Estado de las fases. ``measured`` es lo que imprimen los tests; ``adr`` es dónde está documentado.
FASES = [
    ("0", "Helmholtz 1D por mano", "error L2 < 1 %", "0.399 %", True, "tests/test_helmholtz1d.py"),
    ("1", "FEM acústico 2D + ABC", "onda plana < 1 %", "0.511 %", True, "tests/test_fem_layered.py"),
    ("1", "ABC: reflexión espuria", "< −40 dB", "−73 dB", True, "tests/test_fem_scattered.py"),
    ("2", "Placa fluido–sólido–fluido", "|T| < 2 %", "0.079 %", True, "tests/test_coupled_fluid_solid.py"),
    ("2", "Snell con conversión de modo", "< 1°", "0.42°", True, "tests/test_elastic2d.py"),
    ("3", "BEM: cilindro vs serie", "< 2 %", "0.13 %", True, "tests/test_bem_cylinder.py"),
    ("4", "Contraste FEM vs BEM", "< 2 %", "0.97 % / 0.59 %", True, "tests/test_cross_validation.py"),
    ("4B", "Arreglo vs apertura continua", "FWHM ≤ +2 %", "λ/4 → +1.4 %", True, "D25, D26"),
    ("5a", "Cóclea: cápsula + fluido", "L2 < 5 %", "12 %", False, "D18, D21"),
    ("5b-2", "Membrana basilar", "romper el Q de la cavidad", "no empezada", None, "D24"),
]

CONCLUSIONES = [
    ("Un arreglo empata a un monolítico",
     "Los cuatro arreglos documentados caen dentro del 2 % de la apertura continua: aire 40 kHz +1.3 %, "
     "transcraneal 250 kHz +0.4 %, 500 kHz +1.6 %, diseño 1 MHz +1.8 %."),
    ("Manda la cantidad de elementos, no el pitch",
     "Con 256 elementos y pitch λ/1.9 el error es +1.8 %; con 16 elementos hay que bajar a λ/4 para el 2 %. "
     "El error es de cuantización del perfil de fase."),
    ("A 40 kHz no hay foco que medir",
     "FWHM 26.6 mm contra una cóclea de ~9 mm: en esa banda se mide iluminación, no focalización (D8)."),
]


def _fecha_branch() -> tuple[str, str]:
    fecha = dt.datetime.now().strftime("%Y-%m-%d %H:%M")
    try:
        out = subprocess.run(
            ["git", "log", "-1", "--format=%h %s"], cwd=RAIZ,
            capture_output=True, text=True, timeout=10, check=False,
        )
        commit = out.stdout.strip() or "(sin git)"
    except (OSError, subprocess.SubprocessError):
        commit = "(sin git)"
    return fecha, commit


def _leer_csv(nombre: str) -> list[dict[str, str]]:
    ruta = OUT / nombre
    if not ruta.exists():
        return []
    with ruta.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


def _img_data_uri(ruta: Path) -> str | None:
    if not ruta.exists():
        return None
    datos = base64.b64encode(ruta.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{datos}"


def _adrs() -> list[tuple[str, str]]:
    ruta = RAIZ / "docs" / "decisions.md"
    if not ruta.exists():
        return []
    patron = re.compile(r"^##\s+(D\d+)\s*[—–-]\s*(.+)$", re.MULTILINE)
    return [(m.group(1), m.group(2).strip()) for m in patron.finditer(ruta.read_text(encoding="utf-8"))]


def _tabla_fases() -> str:
    filas = []
    for fase, nombre, criterio, medido, ok, fuente in FASES:
        if ok is True:
            marca, clase = "cumple", "ok"
        elif ok is False:
            marca, clase = "no cumple", "no"
        else:
            marca, clase = "pendiente", "pend"
        filas.append(
            f"<tr><td class='mono'>{html.escape(fase)}</td><td>{html.escape(nombre)}</td>"
            f"<td class='mono'>{html.escape(criterio)}</td><td class='mono'>{html.escape(medido)}</td>"
            f"<td class='{clase}'>{marca}</td><td class='mono dim'>{html.escape(fuente)}</td></tr>"
        )
    return (
        "<table><thead><tr><th>Fase</th><th>Qué se verifica</th><th>Criterio</th>"
        "<th>Medido</th><th>Estado</th><th>Dónde</th></tr></thead><tbody>"
        + "".join(filas) + "</tbody></table>"
    )


def _tabla_focalizacion(filas: list[dict[str, str]]) -> str:
    if not filas:
        return "<p class='dim'>Falta <code>out/focusing_study.csv</code>: correr <code>scripts/focusing_study.py</code>.</p>"
    cuerpo = []
    for f in filas:
        ratio = float(f["fwhm_ratio"])
        desvio = (ratio - 1.0) * 100.0
        clase = "ok" if desvio <= 2.0 else ("warn" if desvio <= 7.0 else "no")
        cuerpo.append(
            "<tr>"
            f"<td class='mono'>{html.escape(f['nombre'])}</td>"
            f"<td class='mono num'>{float(f['frecuencia_khz']):.0f}</td>"
            f"<td class='mono num'>{f['n_elementos']}</td>"
            f"<td class='mono num'>{float(f['pitch_over_lambda']):.3f}</td>"
            f"<td class='mono num'>{float(f['fwhm_mm']):.2f}</td>"
            f"<td class='mono num {clase}'>{desvio:+.1f} %</td>"
            f"<td class='mono num'>{float(f['correlacion']):.4f}</td>"
            f"<td class='mono num'>{float(f['ganancia_db']):.2f}</td>"
            "</tr>"
        )
    return (
        "<table><thead><tr><th>Configuración</th><th>f [kHz]</th><th>N</th><th>pitch/λ</th>"
        "<th>FWHM [mm]</th><th>vs continua</th><th>corr</th><th>ganancia [dB]</th></tr></thead><tbody>"
        + "".join(cuerpo) + "</tbody></table>"
    )


def _bloque_sensibilidad(filas: list[dict[str, str]]) -> str:
    ranking = next((f for f in filas if f.get("parametro") == "ranking"), None)
    if ranking is None:
        return "<p class='dim'>Falta <code>out/sensitivity_cochlea.csv</code>.</p>"
    base = ranking.get("t_interior", "")
    items = []
    for trozo in ranking["valor"].split(";"):
        if "=" in trozo:
            k, v = trozo.split("=", 1)
            items.append((k.strip(), v.strip()))
    filas_html = "".join(
        f"<tr><td>{html.escape(k)}</td><td class='mono num'>{html.escape(v)}</td></tr>" for k, v in items
    )
    return (
        f"<p>Caso base <code>t_interior = {html.escape(base)}</code>. "
        f"{len(filas) - 1} corridas, rango relativo <code>(máx−mín)/base</code>:</p>"
        "<table><thead><tr><th>Parámetro</th><th>Influye</th></tr></thead><tbody>"
        + filas_html + "</tbody></table>"
        "<p class='nota'>El ranking no es definitivo: el 72 % del espesor es en buena parte resonancia "
        "de cavidad, y la membrana basilar (5b-2) es justo lo que rompe ese Q. Hay que re-correrlo "
        "después (D24).</p>"
    )


CSS = """
:root { --bg:#0e1116; --panel:#161b22; --linea:#242c37; --texto:#d7dde5; --dim:#7d8899;
        --ok:#3fb950; --warn:#d29922; --no:#f85149; --acento:#58a6ff; }
* { box-sizing:border-box; }
body { margin:0; padding:0 20px 60px; background:var(--bg); color:var(--texto);
       font:15px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif; }
.wrap { max-width:1080px; margin:0 auto; }
header { padding:38px 0 22px; border-bottom:1px solid var(--linea); margin-bottom:34px; }
h1 { margin:0 0 6px; font-size:27px; letter-spacing:-.4px; }
h2 { margin:44px 0 14px; font-size:19px; letter-spacing:-.2px; }
h2::before { content:""; display:inline-block; width:3px; height:15px; background:var(--acento);
             margin-right:9px; vertical-align:-2px; border-radius:2px; }
.sub { color:var(--dim); font-size:13.5px; }
.sub code { color:var(--texto); }
table { width:100%; border-collapse:collapse; background:var(--panel); border-radius:8px;
        overflow:hidden; font-size:13.5px; margin:10px 0 6px; }
th { text-align:left; font-weight:600; color:var(--dim); font-size:11.5px; text-transform:uppercase;
     letter-spacing:.6px; padding:9px 12px; background:#11161d; border-bottom:1px solid var(--linea); }
td { padding:8px 12px; border-bottom:1px solid #1b212a; }
tr:last-child td { border-bottom:none; }
tr:hover td { background:#1a212a; }
.mono { font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace; font-size:12.5px; }
.num { text-align:right; }
.dim { color:var(--dim); }
.ok { color:var(--ok); } .warn { color:var(--warn); } .no { color:var(--no); } .pend { color:var(--dim); }
.cards { display:grid; grid-template-columns:repeat(auto-fit,minmax(230px,1fr)); gap:12px; margin:14px 0; }
.card { background:var(--panel); border:1px solid var(--linea); border-radius:8px; padding:14px 16px; }
.card h3 { margin:0 0 6px; font-size:14px; color:var(--acento); }
.card p { margin:0; font-size:13px; color:#b9c2ce; }
figure { margin:18px 0; }
figure img { width:100%; border-radius:8px; border:1px solid var(--linea); display:block; }
figcaption { color:var(--dim); font-size:12.5px; margin-top:7px; }
.nota { color:var(--dim); font-size:12.5px; border-left:2px solid var(--linea); padding-left:11px; }
code { background:#1d242e; padding:1px 5px; border-radius:4px; font-size:12.5px; }
ol.adr { columns:2; column-gap:28px; padding-left:22px; font-size:13.5px; }
ol.adr li { margin-bottom:7px; break-inside:avoid; }
footer { margin-top:52px; padding-top:18px; border-top:1px solid var(--linea); color:var(--dim); font-size:12.5px; }
@media (max-width:720px) { ol.adr { columns:1; } table { font-size:12px; } th,td { padding:7px 8px; } }
"""


def main() -> int:
    OUT.mkdir(exist_ok=True)
    fecha, commit = _fecha_branch()
    focal = _leer_csv("focusing_study.csv")
    sens = _leer_csv("sensitivity_cochlea.csv")
    adrs = _adrs()

    figs = []
    for nombre, titulo in [
        ("focusing_summary.png", "FWHM y equivalencia con la apertura continua, por configuración"),
        ("focus_transcraneal_250k.png", "Mapa de intensidad: arreglo transcraneal de 250 kHz, 256 elementos"),
        ("focus_diseno_1m.png", "Mapa de intensidad: diseño de 1 MHz, 256 elementos, pitch 0.521 λ"),
    ]:
        uri = _img_data_uri(OUT / nombre)
        if uri:
            figs.append(f"<figure><img alt='{html.escape(titulo)}' src='{uri}'>"
                        f"<figcaption>{html.escape(titulo)}</figcaption></figure>")

    tarjetas = "".join(
        f"<div class='card'><h3>{html.escape(t)}</h3><p>{html.escape(d)}</p></div>" for t, d in CONCLUSIONES
    )
    lista_adr = "".join(f"<li><b class='mono'>{n}</b> — {html.escape(t)}</li>" for n, t in adrs)

    partes = [
        "<!doctype html><html lang='es'><head><meta charset='utf-8'>",
        "<meta name='viewport' content='width=device-width,initial-scale=1'>",
        "<title>bem-fem-simulator · estado</title>",
        f"<style>{CSS}</style></head><body><div class='wrap'>",
        "<header><h1>bem-fem-simulator — estado</h1>",
        f"<p class='sub'>Generado {html.escape(fecha)} · commit <code>{html.escape(commit)}</code> · "
        f"{len(adrs)} ADR · {len(focal)} configuraciones de focalización · "
        f"{len(sens) - 1 if sens else 0} corridas de sensibilidad</p></header>",
        "<h2>Estado por fase</h2>",
        "<p class='dim'>Cada fila es un criterio que un test hace cumplir. La columna «Dónde» apunta al "
        "test o al ADR que lo respalda; un número sin test no entra en esta tabla.</p>",
        _tabla_fases(),
        "<h2>Focalización: arreglo contra apertura continua</h2>",
        f"<div class='cards'>{tarjetas}</div>",
        _tabla_focalizacion(focal),
        "<p class='nota'>Modelo 2D: un arreglo real se representa por su corte lineal. FWHM medido en el "
        "plano del pico axial real de cada configuración. «vs continua» es el ensanchamiento del foco "
        "respecto de una apertura continua del mismo diámetro: verde ≤ 2 %, ámbar ≤ 7 %, rojo por encima.</p>",
        "".join(figs),
        "<h2>Sensibilidad del modelo de cóclea</h2>",
        _bloque_sensibilidad(sens),
        "<h2>Decisiones (ADR)</h2>",
        f"<ol class='adr'>{lista_adr}</ol>",
        "<h2>Limitaciones conocidas</h2>",
        "<table><thead><tr><th>Limitación</th><th>Consecuencia práctica</th></tr></thead><tbody>"
        "<tr><td>El material manda sobre la malla</td><td>Bajar el error de discretización por debajo "
        "del 8–9 % no cambia conclusiones hasta que mejoren los datos de material (D21, D24)</td></tr>"
        "<tr><td>Sin atenuación de hueso bajo 270 kHz</td><td>Justo la banda de diseño; dos estudios "
        "serios difieren 2.4×–3.3× (D20)</td></tr>"
        "<tr><td>Cápsula ótica sin propiedades medidas</td><td>Sustituto isotrópico declarado, "
        "conservador en transmisión (D17)</td></tr>"
        "<tr><td>Convergencia de malla no monótona</td><td>Refinar la pared puede empeorar el "
        "resultado (D21)</td></tr>"
        "</tbody></table>",
        "<footer>Regenerar: <code>.venv/bin/python scripts/focusing_study.py</code> (139 s) y "
        "<code>.venv/bin/python scripts/report.py</code>.<br>"
        "Los números de la tabla de fases los imprimen los tests; los de focalización y sensibilidad "
        "salen de los CSV de <code>out/</code>.</footer>",
        "</div></body></html>",
    ]
    destino = OUT / "report.html"
    destino.write_text("\n".join(partes), encoding="utf-8")
    print(f"escrito {destino.relative_to(RAIZ)}  ({destino.stat().st_size / 1024:.0f} kB, "
          f"{len(figs)} figuras embebidas, {len(adrs)} ADR)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
