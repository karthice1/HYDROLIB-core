"""Morris indices and a behavioural check from a (possibly unfinished) sensitivity run.

Only whole Morris trajectories whose runs all succeeded are used, so indices can be
computed while the ensemble is still running. A run is called behavioural when its
peak level near Neeleswaram is within ``--hfl-tolerance`` of the 2018 HFL and its
flooded land at the Sentinel-1 pass is at least the radar's open-water area (radar
misses water under tree cover, so it is a lower bound) and at most ``--sar-factor``
times it.

Usage:
    python -m research.periyar.analyse_runs RUNS_DIR
"""

import argparse
from pathlib import Path

import pandas as pd

from research.gsa import analysis

from .sensitivity import FACTORS

OUTPUTS = ["max_flooded_area_km2", "max_flood_volume_mcm", "mean_max_depth_m",
           "hours_to_peak_area", "neeleswaram_peak_level_m", "model_flooded_km2",
           "open_critical_success_index"]


def complete_trajectories(samples: pd.DataFrame, results: pd.DataFrame) -> int:
    per = len(FACTORS) + 1
    ok = results["status"].reindex(samples.index) == "ok"
    n = 0
    while (n + 1) * per <= len(samples) and ok.iloc[n * per:(n + 1) * per].all():
        n += 1
    return n


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("runs_dir", type=Path)
    parser.add_argument("--sar-km2", type=float, default=22.0,
                        help="Radar open-water flood area on land at the pass (km2)")
    parser.add_argument("--sar-factor", type=float, default=3.0)
    parser.add_argument("--hfl-tolerance", type=float, default=0.5, help="m")
    args = parser.parse_args()

    samples = pd.read_csv(args.runs_dir / "samples.csv", index_col=0)
    results = pd.read_csv(args.runs_dir / "results.csv", index_col="run")
    n = complete_trajectories(samples, results)
    used = samples.iloc[:n * (len(FACTORS) + 1)]
    print(f"{n} complete trajectories, {len(used)} runs")

    tables = []
    for metric in OUTPUTS:
        table = analysis.morris(FACTORS, used, results.loc[used.index, metric].astype(float))
        table.insert(0, "metric", metric)
        tables.append(table)
    indices = pd.concat(tables).rename_axis("factor")
    indices.to_csv(args.runs_dir / f"indices_morris_r{n}.csv")
    print(indices.round(3).to_string())

    runs = results.loc[results["status"] == "ok"].join(samples)
    hfl_ok = runs["hfl_error_m"].abs() <= args.hfl_tolerance
    extent_ok = runs["model_flooded_km2"].between(args.sar_km2, args.sar_km2 * args.sar_factor)
    runs["behavioural"] = hfl_ok & extent_ok
    runs.to_csv(args.runs_dir / "runs_with_behaviour.csv")
    print(f"\n{len(runs)} runs: {hfl_ok.sum()} within {args.hfl_tolerance} m of the HFL, "
          f"{extent_ok.sum()} within the radar extent bounds, {runs['behavioural'].sum()} both")
    names = [f.name for f in FACTORS]
    summary = pd.concat({"all": runs[names + OUTPUTS].describe(percentiles=[0.05, 0.5, 0.95]),
                         "behavioural": runs.loc[runs["behavioural"], names + OUTPUTS]
                         .describe(percentiles=[0.05, 0.5, 0.95])})
    summary.to_csv(args.runs_dir / "behavioural_summary.csv")
    print(summary.loc[(slice(None), ["min", "5%", "50%", "95%", "max"]), :].round(2).T.to_string())


if __name__ == "__main__":
    main()
