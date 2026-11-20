"""Apply sampled factor values to D-Flow FM input files with HYDROLIB-core."""

from typing import Iterable

from hydrolib.core.dflowfm.bc.models import ForcingModel
from hydrolib.core.dflowfm.friction.models import FrictionModel, FrictionType
from hydrolib.core.dflowfm.structure.models import StructureModel


def _selected(name: str, targets: Iterable[str]) -> bool:
    targets = tuple(targets)
    return not targets or name in targets


# Friction types whose coefficient is a conveyance (larger value = smoother bed).
# For these the roughness multiplier is applied as a division, so a multiplier
# above 1 always means a rougher bed regardless of the friction formulation.
_CONVEYANCE_TYPES = {FrictionType.chezy, FrictionType.strickler, FrictionType.debosbijkerk}


def _roughen(value: float, friction_type: FrictionType, multiplier: float) -> float:
    if friction_type in _CONVEYANCE_TYPES:
        return value / multiplier
    return value * multiplier


def scale_friction(model: FrictionModel, multiplier: float,
                   targets: Iterable[str] = ()) -> int:
    """Make the bed rougher (multiplier > 1) or smoother (< 1). Returns the number of edits.

    `targets` are friction ids of the [Global] sections. [Branch] sections of a
    file are scaled when one of its [Global] sections is selected, or when no
    targets are given.
    """
    edits = 0
    file_selected = False
    for section in model.global_:
        if _selected(section.frictionid, targets):
            section.frictionvalue = _roughen(section.frictionvalue, section.frictiontype, multiplier)
            edits += 1
            file_selected = True
    if file_selected or not targets:
        for branch in model.branch:
            if branch.frictionvalues:
                branch.frictionvalues = [
                    _roughen(v, branch.frictiontype, multiplier) for v in branch.frictionvalues
                ]
                edits += 1
    return edits


def _scale_column(model: ForcingModel, quantity_prefix: str, transform,
                  targets: Iterable[str]) -> int:
    edits = 0
    for forcing in model.forcing:
        if not _selected(forcing.name, targets) or not forcing.datablock:
            continue
        quantities = [q.quantity for q in forcing.quantityunitpair]
        columns = [i for i, q in enumerate(quantities) if q.startswith(quantity_prefix)]
        if not columns:
            continue
        for row in forcing.datablock:
            for i in columns:
                row[i] = transform(row[i])
        edits += 1
    return edits


def scale_discharge(model: ForcingModel, multiplier: float,
                    targets: Iterable[str] = ()) -> int:
    """Multiply discharge time series (dischargebnd, lateral_discharge, ...)."""
    return _scale_column(model, "discharge", lambda v: v * multiplier, targets) + \
        _scale_column(model, "lateral_discharge", lambda v: v * multiplier, targets)


def offset_waterlevel(model: ForcingModel, offset: float,
                      targets: Iterable[str] = ()) -> int:
    """Add a constant offset to water level boundary time series."""
    return _scale_column(model, "waterlevel", lambda v: v + offset, targets)


def offset_crest(model: StructureModel, offset: float,
                 targets: Iterable[str] = ()) -> int:
    """Add an offset to every numeric crest level of the selected structures.

    Crest levels given as a time series file are left untouched; those must be
    perturbed through their forcing file instead.
    """
    edits = 0
    for structure in model.structure:
        if not _selected(structure.id, targets):
            continue
        level = getattr(structure, "crestlevel", None)
        if isinstance(level, (int, float)):
            structure.crestlevel = level + offset
            edits += 1
    return edits
