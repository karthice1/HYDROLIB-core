"""Write one perturbed copy of a D-Flow FM model per sample."""

import shutil
from pathlib import Path
from typing import Dict, List

import pandas as pd

from hydrolib.core.dflowfm.bc.models import ForcingModel
from hydrolib.core.dflowfm.friction.models import FrictionModel
from hydrolib.core.dflowfm.structure.models import StructureModel

from . import perturb
from .factors import Factor


def build(base_dir: Path, out_dir: Path, samples: pd.DataFrame,
          factors: List[Factor], files: Dict[str, List[str]]) -> pd.DataFrame:
    """Create run directories `out_dir/run_00000`, ... with perturbed inputs.

    Args:
        base_dir: Directory of the calibrated base model.
        out_dir: Directory in which run directories are created.
        samples: Sample table, one column per factor.
        factors: Factor definitions matching the sample columns.
        files: Relative paths of the files to perturb, keyed by
            "friction", "forcing" and "structure".

    Returns:
        A manifest with the run directory and number of edits per factor.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    by_name = {f.name: f for f in factors}
    manifest = []
    for run, row in samples.iterrows():
        run_dir = out_dir / f"run_{run:05d}"
        if run_dir.exists():
            shutil.rmtree(run_dir)
        shutil.copytree(base_dir, run_dir)
        edits = _apply(run_dir, row, by_name, files)
        manifest.append({"run": run, "path": str(run_dir), **row.to_dict(),
                         **{f"edits_{k}": v for k, v in edits.items()}})
    result = pd.DataFrame(manifest).set_index("run")
    result.to_csv(out_dir / "manifest.csv")
    return result


def _apply(run_dir: Path, row: pd.Series, factors: Dict[str, Factor],
           files: Dict[str, List[str]]) -> Dict[str, int]:
    edits = {name: 0 for name in row.index}
    for rel in files.get("friction", []):
        model = FrictionModel(run_dir / rel)
        for name, value in row.items():
            f = factors[name]
            if f.kind == "friction_multiplier":
                edits[name] += perturb.scale_friction(model, value, f.targets)
        model.save()
    for rel in files.get("forcing", []):
        model = ForcingModel(run_dir / rel)
        for name, value in row.items():
            f = factors[name]
            if f.kind == "discharge_multiplier":
                edits[name] += perturb.scale_discharge(model, value, f.targets)
            elif f.kind == "waterlevel_offset":
                edits[name] += perturb.offset_waterlevel(model, value, f.targets)
        model.save()
    for rel in files.get("structure", []):
        model = StructureModel(run_dir / rel)
        for name, value in row.items():
            f = factors[name]
            if f.kind == "crest_offset":
                edits[name] += perturb.offset_crest(model, value, f.targets)
        model.save()
    return edits
