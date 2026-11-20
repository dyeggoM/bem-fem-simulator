"""Dashboard local interactivo del simulador.

Uso::

    .venv/bin/python -m ultrasound_sim.dashboard          # http://127.0.0.1:8765/
    .venv/bin/python -m ultrasound_sim.dashboard --port 9000
"""

from .field import DEFAULT_QUAD, Source, array_sources, evaluate, metrics
from .server import build_server, serve

__all__ = ["DEFAULT_QUAD", "Source", "array_sources", "build_server", "evaluate", "metrics", "serve"]
