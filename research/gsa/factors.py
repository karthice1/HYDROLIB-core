"""Uncertain input factors for the GSA and their sampling bounds."""

from dataclasses import dataclass
from typing import Dict, List


@dataclass(frozen=True)
class Factor:
    """One uncertain input of the model.

    Attributes:
        name: Unique factor name, used as column name in the sample table.
        lower: Lower bound of the uniform sampling range.
        upper: Upper bound of the uniform sampling range.
        kind: How the factor is applied to the model. One of
            "friction_multiplier", "discharge_multiplier", "waterlevel_offset"
            or "crest_offset".
        targets: Names the factor applies to (friction ids, forcing names or
            structure ids). An empty list means "all of that kind".
        description: Human-readable description for reports.
    """

    name: str
    lower: float
    upper: float
    kind: str
    targets: tuple = ()
    description: str = ""


KINDS = {
    "friction_multiplier",
    "discharge_multiplier",
    "waterlevel_offset",
    "crest_offset",
}


def validate(factors: List[Factor]) -> None:
    """Raise a ValueError when the factor list is inconsistent."""
    names = [f.name for f in factors]
    if len(set(names)) != len(names):
        raise ValueError(f"Factor names must be unique, got {names}.")
    for f in factors:
        if f.kind not in KINDS:
            raise ValueError(f"Unknown factor kind '{f.kind}' for '{f.name}'.")
        if not f.lower < f.upper:
            raise ValueError(f"Factor '{f.name}' needs lower < upper.")


def to_salib_problem(factors: List[Factor]) -> Dict:
    """Convert factors to the problem dictionary used by SALib."""
    validate(factors)
    return {
        "num_vars": len(factors),
        "names": [f.name for f in factors],
        "bounds": [[f.lower, f.upper] for f in factors],
    }


# Provisional factor set for the Periyar August 2018 flood. Ranges are first
# estimates and must be justified from literature before the production runs.
PERIYAR_2018 = [
    Factor("n_channel", 0.7, 1.3, "friction_multiplier", ("Channel",),
           "Main channel roughness, multiplier on calibrated value"),
    Factor("n_floodplain", 0.6, 1.5, "friction_multiplier", ("Floodplain",),
           "Floodplain roughness, multiplier on calibrated value"),
    Factor("q_idukki", 0.8, 1.2, "discharge_multiplier", ("Idukki_release",),
           "Idukki dam release, multiplier on reported hydrograph"),
    Factor("q_idamalayar", 0.8, 1.2, "discharge_multiplier", ("Idamalayar_release",),
           "Idamalayar dam release, multiplier on reported hydrograph"),
    Factor("q_lateral", 0.7, 1.3, "discharge_multiplier", ("Lateral_inflow",),
           "Ungauged lateral inflow from the intermediate catchment"),
    Factor("tide_offset", -0.2, 0.3, "waterlevel_offset", ("Kochi_tide",),
           "Downstream tidal water level offset in metres"),
    Factor("crest_offset", -0.3, 0.3, "crest_offset", (),
           "Crest level error of weirs and regulators in metres"),
]
