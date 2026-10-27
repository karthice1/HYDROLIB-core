"""Flooded land through time for selected runs of a sensitivity ensemble.

The ensemble keeps only scalar metrics per run, so the selected runs (by default the
behavioural ones from analyse_runs.py) are run again and their flooded land on cells
that were dry at the start is written for every map output time.

Usage:
    python -m research.periyar.time_series BASE_DIR RUNS_DIR OUT_CSV [--runs 3 17 ...]
"""

import argparse
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import netCDF4
import pandas as pd

from .build_base_model import REFERENCE
from .sensitivity import WET_DEPTH, prepare


def flooded_land(map_file: Path) -> pd.Series:
    with netCDF4.Dataset(map_file) as nc:
        hours = nc["time"][:] / 3600.0
        depth = nc["mesh2d_waterdepth"][:]
        area = nc["mesh2d_flowelem_ba"][:]
    land = depth[0] <= WET_DEPTH
    km2 = [(area[land & (d > WET_DEPTH)]).sum() / 1e6 for d in depth]
    times = REFERENCE + pd.to_timedelta(hours, unit="h")
    return pd.Series(km2, index=pd.Index(times, name="time_ist"), name="flooded_land_km2")


def run(dflowfm: str, base_dir: Path, run_dir: Path, row: pd.Series) -> pd.Series:
    done = run_dir / "flooded_land.csv"
    if done.exists():  # finished before a restart
        return pd.read_csv(done, index_col=0, parse_dates=True)["flooded_land_km2"]
    prepare(base_dir, run_dir, row)
    with (run_dir / "run.log").open("w") as out:
        subprocess.run([dflowfm, "--autostartstop", "periyar.mdu"], cwd=run_dir, check=True,
                       stdout=out, stderr=subprocess.STDOUT,
                       env={**os.environ, "OMP_NUM_THREADS": "1"})
    map_file = run_dir / "DFM_OUTPUT_periyar" / "periyar_map.nc"
    series = flooded_land(map_file)
    series.to_csv(done)
    map_file.unlink()
    return series


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("base_dir", type=Path)
    parser.add_argument("runs_dir", type=Path)
    parser.add_argument("out_csv", type=Path)
    parser.add_argument("--runs", type=int, nargs="+", help="Default: the behavioural runs")
    parser.add_argument("--parallel", type=int, default=4)
    parser.add_argument("--dflowfm", default="dflowfm")
    args = parser.parse_args()

    samples = pd.read_csv(args.runs_dir / "samples.csv", index_col=0)
    runs = args.runs
    if runs is None:
        judged = pd.read_csv(args.runs_dir / "runs_with_behaviour.csv", index_col=0)
        runs = judged.index[judged["behavioural"]].tolist()
    work = args.runs_dir / "time_series"

    def task(run_id):
        return run_id, run(args.dflowfm, args.base_dir, work / f"run_{run_id:04d}",
                           samples.loc[run_id])

    with ThreadPoolExecutor(max_workers=args.parallel) as pool:
        series = dict(pool.map(task, runs))
    table = pd.DataFrame(series)
    table.columns.name = "run"
    table.to_csv(args.out_csv)
    print(table.max().round(1).to_string())


if __name__ == "__main__":
    main()
