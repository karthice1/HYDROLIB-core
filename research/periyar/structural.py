"""Structural tests of the Periyar model at nominal parameter values.

Each variant rebuilds the base model with one structural choice changed (terrain,
channel mask, grid size or tide), runs it once, and stores its flood metrics, the
full confusion matrix against the Sentinel-1 map, the flooded land per output time
and the wet/dry state of every cell at the radar pass. Map files are deleted after
scoring.

Usage:
    python -m research.periyar.structural DATA_DIR SAR_TIF OUT_DIR [--variants nominal min_terrain ...]
"""

import argparse
import json
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd

from . import validate
from .build_base_model import build
from .sensitivity import metrics
from .time_series import flooded_land

# Keyword arguments of build_base_model.build that differ from the v5 model.
VARIANTS = {
    "nominal": {},
    "min_terrain": {"terrain": "min"},
    "full_canopy": {"canopy_factor": 1.0},
    "wide_channel": {"water_threshold": 25.0},
    "grid_100m": {"cell_size": 100.0},
    "grid_400m": {"cell_size": 400.0},
    "tide_x1.5": {"tide_scale": 1.5},
    "tide_shift_6h": {"tide_shift_hours": 6.21},
}
BASE = dict(cell_size=200.0, channel_depth=4.0, max_elevation=20.0)


def radar_state(model_dir: Path, map_file: Path, sar_file: Path) -> dict:
    """Per-cell model and radar flooding at the radar pass, for maps."""
    with netCDF4.Dataset(map_file) as nc:
        x, y = nc["mesh2d_face_x"][:], nc["mesh2d_face_y"][:]
        hours = nc["time"][:] / 3600.0
        depth = nc["mesh2d_waterdepth"][:]
    sar_hour = (pd.Timestamp(validate.SAR_TIME) - validate.REFERENCE).total_seconds() / 3600
    step = int(np.argmin(np.abs(hours - sar_hour)))
    fraction = validate._sample(sar_file, x, y, -1.0)
    landcover = validate._sample(model_dir / "landcover.tif", x, y, 0)
    return dict(x=np.asarray(x), y=np.asarray(y), land=np.asarray(depth[0] <= validate.WET_DEPTH),
                model_wet=np.asarray(depth[step] > validate.WET_DEPTH),
                radar_fraction=np.asarray(fraction), landcover=np.asarray(landcover))


def run_variant(name: str, data_dir: Path, sar_file: Path, out_dir: Path, dflowfm: str,
                threads: int) -> dict:
    model_dir = out_dir / name
    summary = model_dir / "summary.json"
    if summary.exists():
        return json.loads(summary.read_text())
    build(data_dir, model_dir, **{**BASE, **VARIANTS[name]})
    with (model_dir / "run.log").open("w") as out:
        subprocess.run([dflowfm, "--autostartstop", "periyar.mdu"], cwd=model_dir, check=True,
                       stdout=out, stderr=subprocess.STDOUT,
                       env={**os.environ, "OMP_NUM_THREADS": str(threads)})
    map_file = model_dir / "DFM_OUTPUT_periyar" / "periyar_map.nc"
    values = {"variant": name, **VARIANTS[name], **metrics(map_file)}
    values.update(validate.scores(model_dir, map_file, sar_file))
    flooded_land(map_file).to_csv(model_dir / "flooded_land.csv")
    np.savez_compressed(model_dir / "radar_state.npz", **radar_state(model_dir, map_file, sar_file))
    map_file.unlink()
    summary.write_text(json.dumps(values, indent=2))
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("data_dir", type=Path)
    parser.add_argument("sar_file", type=Path)
    parser.add_argument("out_dir", type=Path)
    parser.add_argument("--variants", nargs="+", default=list(VARIANTS))
    parser.add_argument("--parallel", type=int, default=4)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--dflowfm", default="dflowfm")
    args = parser.parse_args()

    def task(name):
        return run_variant(name, args.data_dir, args.sar_file, args.out_dir, args.dflowfm,
                           args.threads)

    with ThreadPoolExecutor(max_workers=args.parallel) as pool:
        rows = list(pool.map(task, args.variants))
    table = pd.DataFrame(rows).set_index("variant")
    table.to_csv(args.out_dir / "structural_tests.csv")
    print(table[["max_flooded_area_km2", "model_flooded_km2", "hfl_error_m",
                 "open_critical_success_index", "all_critical_success_index"]].round(3).to_string())


if __name__ == "__main__":
    main()
