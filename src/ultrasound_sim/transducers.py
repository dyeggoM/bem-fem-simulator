"""Importación de transductores desde un catálogo editable (``data/transducers.toml``).

El catálogo es la única fuente de geometría de transductores reales. El importador no inventa nada:
si un campo falta, falla nombrando el campo y el equipo; si la geometría no es la que el modelo
soporta (elementos con kerf), falla en vez de aproximarla.

Del datasheet al campo: un ``TransducerSpec`` se convierte en un :class:`arrays.LinearArray`, y de
ahí en aperturas con la fase de foco o de apuntamiento. Todo lo demás (presión, métricas) es común
con el resto del paquete.

Ejemplo
-------
>>> cat = load_catalog()
>>> spec = cat["referencia_lambda8_100k"]
>>> cat_aperturas = spec.apertures(materials.WATER, focus_point=(0.04, 0.0))

Regla (docs/decisions.md D5): un equipo con ``verified = False`` se puede simular, pero no entra en
una validación como si sus números fueran medidos.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

from .aperture import Aperture
from .arrays import LinearArray
from .materials import Material

__all__ = [
    "TransducerSpec",
    "CatalogError",
    "default_catalog_path",
    "load_catalog",
    "parse_catalog",
    "catalog_summary",
    "grating_lobe_angle_deg",
]

_REQUIRED = ("name", "frequency_hz", "n_elements", "pitch_m", "source")
_KNOWN = set(_REQUIRED) | {
    "element_width_m",
    "focus_distance_m",
    "amplitude",
    "verified",
    "notes",
}


class CatalogError(ValueError):
    """Error de contenido del catálogo (campo faltante, tipo, geometría no soportada)."""


def default_catalog_path() -> Path:
    """``<repo>/data/transducers.toml`` (resuelto desde este archivo, no desde el cwd)."""
    return Path(__file__).resolve().parents[2] / "data" / "transducers.toml"


def _as_float(data: dict, key: str, entry: str, required: bool) -> float | None:
    if key not in data or data[key] is None:
        if required:
            raise CatalogError(f"transductor '{entry}': falta el campo obligatorio '{key}'")
        return None
    value = data[key]
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise CatalogError(
            f"transductor '{entry}': '{key}' debe ser numérico; recibido {value!r}"
        )
    if not math.isfinite(float(value)):
        raise CatalogError(f"transductor '{entry}': '{key}' debe ser finito")
    return float(value)


def _opt_float(data: dict, key: str, entry: str, default: float) -> float:
    """``_as_float`` para campos opcionales: ``None``/ausente ⇒ ``default`` (no ``0.0``)."""
    value = _as_float(data, key, entry, False)
    return default if value is None else value


def _req_float(data: dict, key: str, entry: str) -> float:
    """``_as_float`` para campos obligatorios: devuelve ``float``, nunca ``None``."""
    value = _as_float(data, key, entry, True)
    assert value is not None  # _as_float con required=True o lanza o devuelve float
    return value


@dataclass(frozen=True)
class TransducerSpec:
    """Un transductor del catálogo.

    ``n_elements = 1`` es un monolítico; ``> 1`` un arreglo lineal de elementos contiguos.
    ``focus_distance_m`` es la distancia focal *geométrica* del datasheet (la fase se calcula con
    ella); recuérdese D13: en f/# bajo el máximo real no cae en la focal geométrica.
    """

    key: str
    name: str
    frequency_hz: float
    n_elements: int
    pitch_m: float
    source: str
    element_width_m: float | None = None
    focus_distance_m: float | None = None
    amplitude: float = 1.0
    verified: bool = False
    notes: str = ""
    raw: dict = field(default_factory=dict, repr=False, compare=False)

    # ---- validación --------------------------------------------------------
    def __post_init__(self) -> None:
        if not (self.frequency_hz > 0.0):
            raise CatalogError(
                f"transductor '{self.key}': frequency_hz debe ser > 0; recibido {self.frequency_hz!r}"
            )
        if int(self.n_elements) < 1:
            raise CatalogError(
                f"transductor '{self.key}': n_elements debe ser >= 1; recibido {self.n_elements!r}"
            )
        if not (self.pitch_m > 0.0):
            raise CatalogError(
                f"transductor '{self.key}': pitch_m debe ser > 0; recibido {self.pitch_m!r}"
            )
        if not str(self.source).strip():
            raise CatalogError(
                f"transductor '{self.key}': 'source' no puede estar vacío (procedencia obligatoria)"
            )
        if self.focus_distance_m is not None and not (self.focus_distance_m > 0.0):
            raise CatalogError(
                f"transductor '{self.key}': focus_distance_m debe ser > 0 si se especifica"
            )
        if not (self.amplitude > 0.0):
            raise CatalogError(f"transductor '{self.key}': amplitude debe ser > 0")
        # El modelo de apertura asume elementos contiguos (sin kerf). Mejor fallar que aproximar.
        if self.element_width_m is not None and not math.isclose(
            self.element_width_m, self.pitch_m, rel_tol=1e-9, abs_tol=1e-15
        ):
            raise CatalogError(
                f"transductor '{self.key}': element_width_m ({self.element_width_m}) != pitch_m "
                f"({self.pitch_m}); el modelo de apertura exige elementos contiguos sin kerf"
            )

    # ---- geometría derivada ------------------------------------------------
    @property
    def aperture_width_m(self) -> float:
        """Ancho total de la cara activa, ``N · pitch`` [m]."""
        return self.n_elements * self.pitch_m

    @property
    def f_number(self) -> float | None:
        """``F/D`` si el transductor trae focal; ``None`` si es plano."""
        if self.focus_distance_m is None:
            return None
        return self.focus_distance_m / self.aperture_width_m

    def wavelength(self, material: Material) -> float:
        return material.wavelength(self.frequency_hz)

    def pitch_over_wavelength(self, material: Material) -> float:
        return self.pitch_m / self.wavelength(material)

    def aperture_over_wavelength(self, material: Material) -> float:
        return self.aperture_width_m / self.wavelength(material)

    # ---- al modelo acústico ------------------------------------------------
    def as_array(self, center_y: float = 0.0) -> LinearArray:
        return LinearArray(center_y=center_y, n_elements=self.n_elements, pitch=self.pitch_m)

    def apertures(
        self,
        material: Material,
        center_y: float = 0.0,
        focus_point: tuple[float, float] | None = None,
        steer_angle: float | None = None,
        frequency: float | None = None,
    ) -> list[Aperture]:
        """Aperturas listas para :func:`aperture.pressure`.

        Si no se pasa ``focus_point`` y el transductor declara ``focus_distance_m``, se enfoca en
        ``(focus_distance_m, 0)``. ``steer_angle`` [rad] apunta el haz (no se combina con foco: el
        foco manda, porque una fase sola no puede hacer las dos cosas sin un frente de onda).
        """
        f = self.frequency_hz if frequency is None else float(frequency)
        array = self.as_array(center_y)
        if focus_point is None and self.focus_distance_m is not None:
            focus_point = (float(self.focus_distance_m), 0.0)
        if focus_point is not None:
            aps = array.focused(focus_point, f, material)
        elif steer_angle is not None:
            aps = array.steered(float(steer_angle), f, material)
        else:
            aps = array.apertures()
        return [Aperture(a.center_y, a.width, a.amplitude * self.amplitude) for a in aps]


def _spec_from_mapping(key: str, data: dict) -> TransducerSpec:
    if not isinstance(data, dict):
        raise CatalogError(f"transductor '{key}': debe ser una tabla")
    unknown = set(data) - _KNOWN
    if unknown:
        raise CatalogError(
            f"transductor '{key}': campos desconocidos {sorted(unknown)} "
            f"(válidos: {sorted(_KNOWN)})"
        )
    for missing in _REQUIRED:
        if missing not in data:
            raise CatalogError(f"transductor '{key}': falta el campo obligatorio '{missing}'")

    n_elements = data["n_elements"]
    if isinstance(n_elements, bool) or not isinstance(n_elements, int):
        raise CatalogError(
            f"transductor '{key}': 'n_elements' debe ser entero; recibido {n_elements!r}"
        )

    return TransducerSpec(
        key=key,
        name=str(data["name"]),
        frequency_hz=_req_float(data, "frequency_hz", key),
        n_elements=int(n_elements),
        pitch_m=_req_float(data, "pitch_m", key),
        source=str(data["source"]),
        element_width_m=_as_float(data, "element_width_m", key, False),
        focus_distance_m=_as_float(data, "focus_distance_m", key, False),
        amplitude=_opt_float(data, "amplitude", key, 1.0),
        verified=bool(data.get("verified", False)),
        notes=str(data.get("notes", "")),
        raw=dict(data),
    )


def parse_catalog(text: str) -> dict[str, TransducerSpec]:
    """Parsea el texto de un catálogo TOML. Devuelve ``{id: TransducerSpec}``."""
    try:
        import tomllib
    except ModuleNotFoundError:  # Python 3.10
        try:
            import tomli as tomllib  # type: ignore[no-redef]
        except ModuleNotFoundError as exc:  # pragma: no cover
            raise CatalogError(
                "hace falta tomllib (Python >= 3.11) o el paquete 'tomli' para leer el catálogo"
            ) from exc

    try:
        data = tomllib.loads(text)
    except Exception as exc:  # tomllib.TOMLDecodeError no existe en 3.10
        raise CatalogError(f"catálogo TOML inválido: {exc}") from exc

    entries = data.get("transducers", {})
    if not isinstance(entries, dict):
        raise CatalogError("la sección [transducers] debe contener tablas, una por equipo")
    return {key: _spec_from_mapping(key, entry) for key, entry in entries.items()}


def load_catalog(path: str | Path | None = None) -> dict[str, TransducerSpec]:
    """Carga el catálogo. Sin argumento usa ``data/transducers.toml`` del repo."""
    p = Path(path) if path is not None else default_catalog_path()
    if not p.is_file():
        raise CatalogError(
            f"no existe el catálogo '{p}'. Copiá la plantilla de data/transducers.toml y completala"
        )
    return parse_catalog(p.read_text(encoding="utf-8"))


def grating_lobe_angle_deg(spec: TransducerSpec, material: Material,
                           steer_angle: float = 0.0) -> float | None:
    """Ángulo del primer lóbulo de grataje [°], o ``None`` si ``pitch < λ/2``.

    ``sinθ_g = sinθ_0 + m·λ/d`` con ``m = ±1`` (docs/plans, Fase 4B). Si el seno se sale de rango,
    no hay lóbulo de grataje real.
    """
    ratio = spec.pitch_over_wavelength(material)
    if ratio < 0.5 - 1e-12:
        return None
    s = math.sin(steer_angle) + material.wavelength(spec.frequency_hz) / spec.pitch_m
    if abs(s) > 1.0:
        return None
    return math.degrees(math.asin(s))


def catalog_summary(catalog: dict[str, TransducerSpec], material: Material) -> list[dict]:
    """Tabla lista para imprimir: geometría en longitudes de onda y lóbulo de grataje."""
    rows = []
    for key, spec in catalog.items():
        rows.append(
            dict(
                key=key,
                name=spec.name,
                frequency_hz=spec.frequency_hz,
                n_elements=spec.n_elements,
                width_mm=spec.aperture_width_m * 1e3,
                pitch_over_lambda=spec.pitch_over_wavelength(material),
                aperture_over_lambda=spec.aperture_over_wavelength(material),
                f_number=spec.f_number,
                grating_lobe_deg=grating_lobe_angle_deg(spec, material),
                verified=spec.verified,
            )
        )
    return rows
