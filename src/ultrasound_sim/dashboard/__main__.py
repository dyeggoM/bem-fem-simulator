"""``python -m ultrasound_sim.dashboard``"""

from __future__ import annotations

import argparse

from .server import serve


def main() -> int:
    p = argparse.ArgumentParser(description="Dashboard local del simulador de ultrasonido")
    p.add_argument("--host", default="127.0.0.1", help="interfaz (por defecto sólo local)")
    p.add_argument("--port", type=int, default=8765, help="puerto (0 = uno libre)")
    args = p.parse_args()
    serve(args.host, args.port)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
