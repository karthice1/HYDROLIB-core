"""Global sensitivity analysis of the Periyar 2018 flood model.

Each run is a copy of the base model (see build_base_model.py) in which the bed
level, roughness and boundary forcing are rewritten from one sample of the
uncertain factors below. Runs are executed with D-Flow FM and every map output is
reduced to a few scalar flood metrics before it is deleted, so hundreds of runs fit
on disk.

Usage:
    python -m research.periyar.sensitivity BASE_DIR RUNS_DIR --method morris --trajectories 4
"""

import argparse
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd
import rasterio

from hydrolib.core.dflowfm.bc.models import ForcingModel
from research.gsa import analysis, perturb, sampling
from research.gsa.factors import Factor

from . import validate
from .build_base_model import MANNING, REFERENCE

DATA = Path(__file__).parent / "data"

# Ranges are first estimates and are discussed in the manuscript. Travel times from
# the dams to Neeleswaram are not published; the lags below are assumptions.
FACTORS = [
    Factor("n_channel", 0.7, 1.3, "custom", (), "Channel Manning n, multiplier"),
    Factor("n_floodplain", 0.6, 1.5, "custom", (), "Floodplain Manning n, multiplier"),
    Factor("channel_depth", 2.0, 8.0, "custom", (), "Channel depth below DEM water surface (m)"),
    Factor("canopy_factor", 0.3, 0.9, "custom", (), "Fraction of canopy height in the DSM"),
    Factor("q_rating", 0.85, 1.15, "custom", (), "Neeleswaram rating-curve error, multiplier"),
    Factor("q_idukki", 0.8, 1.2, "custom", (), "Idukki spill, multiplier"),
    Factor("q_idamalayar", 0.8, 1.2, "custom", (), "Idamalayar release, multiplier"),
    Factor("tide_offset", -0.2, 0.5, "custom", (), "Sea level offset incl. surge (m)"),
]
# Values of the factors that a Latin hypercube run does not vary.
NOMINAL = {"n_channel": 1.0, "n_floodplain": 1.0, "channel_depth": 4.0, "canopy_factor": 0.6,
           "q_rating": 1.0, "q_idukki": 1.0, "q_idamalayar": 1.0, "tide_offset": 0.0}
IDUKKI_LAG_HOURS = 8.0
IDAMALAYAR_LAG_HOURS = 4.0
WET_DEPTH = 0.1  # m, depth above which a cell counts as flooded


def _release_hours(lag_hours):
    """Hourly Idukki spill and Idamalayar release on the model time axis, shifted by lag."""
    idk = pd.read_csv(DATA / "idukki_cwc_fig4.csv", parse_dates=["datetime_ist"])
    ida = pd.read_csv(DATA / "idamalayar_cwc_fig5.csv", parse_dates=["date"])
    t_idk = (idk["datetime_ist"] - REFERENCE).dt.total_seconds() / 3600 + lag_hours[0]
    # Daily values are plotted at midnight; treat them as daily means centred at noon.
    t_ida = (ida["date"] - REFERENCE).dt.total_seconds() / 3600 + 12 + lag_hours[1]
    return (t_idk.to_numpy(), idk["outflow_m3s"].to_numpy(float)), \
        (t_ida.to_numpy(), ida["outflow_m3s"].to_numpy(float))


def _rewrite_rasters(run_dir: Path, base_dir: Path, row: pd.Series) -> None:
    def read(name):
        with rasterio.open(base_dir / name) as src:
            return src.read(1), src.profile
    dsm, profile = read("dem_raw.tif")
    dsm_min, _ = read("dem_min.tif")
    canopy, _ = read("chm.tif")
    channel = read("channel.tif")[0] > 0.5
    landcover, _ = read("landcover.tif")
    valid = dsm > -9000
    canopy = np.where((canopy >= 0) & (canopy < 100), canopy, 0.0)
    bed = np.where(channel, dsm_min - row["channel_depth"], dsm - row["canopy_factor"] * canopy)
    manning = np.vectorize(lambda c: MANNING.get(int(c), 0.04), otypes=["float32"])(landcover)
    manning = np.where(channel, MANNING[80] * row["n_channel"], manning * row["n_floodplain"])
    for name, array in (("bedlevel.tif", bed), ("manning.tif", manning)):
        with rasterio.open(run_dir / name, "w", **profile) as dst:
            dst.write(np.where(valid, array, profile["nodata"]).astype("float32"), 1)


def _rewrite_forcing(run_dir: Path, row: pd.Series, release_model: str = "additive") -> None:
    """Perturb the Neeleswaram inflow and the sea level.

    The Neeleswaram hydrograph Q(t) already contains the dam releases, so the release
    multipliers m only add or remove the share (m - 1) of each lagged release:
    "additive":       q_rating * Q(t) + sum (m - 1) R(t - lag)
    "multiplicative": q_rating * [Q(t) + sum (m - 1) R(t - lag)]
    At m = 1 no water is added in either form.
    """
    model = ForcingModel(run_dir / "boundaries.bc")
    (t_idk, q_idk), (t_ida, q_ida) = _release_hours((IDUKKI_LAG_HOURS, IDAMALAYAR_LAG_HOURS))
    for forcing in model.forcing:
        if forcing.name.startswith("upstream"):
            for record in forcing.datablock:
                t = record[0]
                extra = (row["q_idukki"] - 1) * np.interp(t, t_idk, q_idk) \
                    + (row["q_idamalayar"] - 1) * np.interp(t, t_ida, q_ida)
                if release_model == "multiplicative":
                    record[1] = max((record[1] + extra) * row["q_rating"], 0.0)
                else:
                    record[1] = max(record[1] * row["q_rating"] + extra, 0.0)
    perturb.offset_waterlevel(model, row["tide_offset"])
    model.save()


def prepare(base_dir: Path, run_dir: Path, row: pd.Series, release_model: str = "additive") -> None:
    if run_dir.exists():
        shutil.rmtree(run_dir)
    shutil.copytree(base_dir, run_dir, ignore=shutil.ignore_patterns("DFM_OUTPUT_*", "*.log"))
    _rewrite_rasters(run_dir, base_dir, row)
    _rewrite_forcing(run_dir, row, release_model)


def metrics(map_file: Path, model_dir: Path = None, sar_file: Path = None) -> dict:
    """Reduce a map file to scalar flood metrics.

    Cells that are already wet at the start (sea, backwaters, river channel) are
    excluded, so the metrics describe land that the flood inundated.
    """
    with netCDF4.Dataset(map_file) as nc:
        depth = nc["mesh2d_waterdepth"][:]
        area = nc["mesh2d_flowelem_ba"][:]
    land = depth[0] <= WET_DEPTH
    wet = (depth > WET_DEPTH) & land
    flooded = (wet * area).sum(axis=1) / 1e6
    max_depth = np.where(land, depth.max(axis=0), 0.0)
    ever_wet = max_depth > WET_DEPTH
    values = {
        "max_flooded_area_km2": float(flooded.max()),
        "hours_to_peak_area": int(flooded.argmax()),
        "max_inundation_extent_km2": float((ever_wet * area).sum() / 1e6),
        "mean_max_depth_m": float(max_depth[ever_wet].mean()) if ever_wet.any() else 0.0,
        "max_flood_volume_mcm": float((np.where(land, depth, 0) * area).sum(axis=1).max() / 1e6),
    }
    if sar_file is not None:
        v = validate.scores(model_dir, map_file, sar_file)
        values.update({key: v[key] for key in (
            "neeleswaram_peak_level_m", "hfl_error_m", "model_flooded_km2",
            "all_critical_success_index", "open_critical_success_index", "open_hit_rate",
            "open_false_alarm_ratio")})
    return values


def run_one(dflowfm: str, base_dir: Path, run_dir: Path, row: pd.Series, threads: int,
            sar_file: Path = None, release_model: str = "additive", keep_map: bool = False) -> dict:
    prepare(base_dir, run_dir, row, release_model)
    log = run_dir / "run.log"
    env = {"OMP_NUM_THREADS": str(threads)}
    with log.open("w") as out:
        result = subprocess.run([dflowfm, "--autostartstop", "periyar.mdu"], cwd=run_dir,
                                stdout=out, stderr=subprocess.STDOUT,
                                env={**os.environ, **env})
    map_file = run_dir / "DFM_OUTPUT_periyar" / "periyar_map.nc"
    if result.returncode != 0 or not map_file.exists():
        return {"status": "failed"}
    values = metrics(map_file, run_dir, sar_file)
    if not keep_map:
        map_file.unlink()  # keep the disk free; scalar metrics are all the GSA needs
    return {"status": "ok", **values}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("base_dir", type=Path)
    parser.add_argument("runs_dir", type=Path)
    parser.add_argument("--method", choices=["morris", "sobol", "lhs", "file"], default="morris")
    parser.add_argument("--samples-file", type=Path, help="Sample table for --method file; "
                        "missing factors are held at NOMINAL")
    parser.add_argument("--release-model", choices=["additive", "multiplicative"],
                        default="additive", help="How the release multipliers enter the inflow")
    parser.add_argument("--keep-maps", action="store_true", help="Keep the map files")
    parser.add_argument("--trajectories", type=int, default=4, help="Morris trajectories")
    parser.add_argument("--base-samples", type=int, default=64, help="Sobol base samples")
    parser.add_argument("--samples", type=int, default=100, help="Latin hypercube runs")
    parser.add_argument("--seed", type=int, default=42, help="Latin hypercube seed")
    parser.add_argument("--factors", nargs="+", help="Vary only these factors (lhs, morris); "
                        "the others are held at NOMINAL")
    parser.add_argument("--parallel", type=int, default=2, help="Simultaneous runs")
    parser.add_argument("--threads", type=int, default=2, help="OpenMP threads per run")
    parser.add_argument("--dflowfm", default="dflowfm", help="D-Flow FM executable or wrapper")
    parser.add_argument("--sar", type=Path, help="Sentinel-1 flood map; adds validation scores")
    args = parser.parse_args()

    varied = [f for f in FACTORS if not args.factors or f.name in args.factors]
    if args.method == "morris":
        samples = sampling.morris(varied, trajectories=args.trajectories)
        for name, value in NOMINAL.items():
            if name not in samples:
                samples[name] = value
        samples = samples[[f.name for f in FACTORS]]
    elif args.method == "sobol":
        samples = sampling.sobol(FACTORS, base_samples=args.base_samples)
    elif args.method == "file":
        samples = pd.read_csv(args.samples_file, index_col=0)
        for name, value in NOMINAL.items():
            if name not in samples:
                samples[name] = value
        samples = samples[[f.name for f in FACTORS]]
    else:
        samples = sampling.latin_hypercube(varied, samples=args.samples, seed=args.seed)
        for name, value in NOMINAL.items():
            if name not in samples:
                samples[name] = value
        samples = samples[[f.name for f in FACTORS]]
    args.runs_dir.mkdir(parents=True, exist_ok=True)
    samples.to_csv(args.runs_dir / "samples.csv")

    results_file = args.runs_dir / "results.csv"
    done = pd.read_csv(results_file, index_col="run") if results_file.exists() else pd.DataFrame()

    def task(run):
        if run in done.index and done.loc[run, "status"] == "ok":
            return run, done.loc[run].to_dict()
        out = run_one(args.dflowfm, args.base_dir, args.runs_dir / f"run_{run:04d}",
                      samples.loc[run], args.threads, args.sar, args.release_model,
                      args.keep_maps)
        print(f"run {run}: {out}", flush=True)
        return run, out

    rows = {}
    with ThreadPoolExecutor(max_workers=args.parallel) as pool:
        for run, out in pool.map(task, samples.index):
            rows[run] = out
            pd.DataFrame.from_dict(rows, orient="index").rename_axis("run").sort_index() \
                .to_csv(results_file)

    results = pd.DataFrame.from_dict(rows, orient="index").sort_index()
    ok = results["status"] == "ok"
    if not ok.all():
        print(f"{(~ok).sum()} runs failed; indices are computed only when all runs succeed.")
        return
    if args.method in ("lhs", "file"):
        return  # conditioned on the observations by analyse_runs.py, no indices
    report = []
    for metric in [c for c in results.columns if c != "status"]:
        y = results[metric].astype(float)
        if args.method == "morris":
            table = analysis.morris(varied, samples[[f.name for f in varied]], y)
        else:
            table = analysis.sobol(FACTORS, y)
        table.insert(0, "metric", metric)
        report.append(table)
    pd.concat(report).rename_axis("factor").to_csv(args.runs_dir / f"indices_{args.method}.csv")
    print(pd.concat(report).to_string())


if __name__ == "__main__":
    main()
