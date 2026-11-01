"""Forcing variants for 14 stratified members of the 50 runs accepted at +-1.0 m but not at
+-0.5 m (round 5). Compares each rerun with its original in the pooled 600-run ensemble.

The 14 runs were chosen from those 50, sorted by peak flooded land: 12 evenly spaced in that
order (including the smallest and largest) plus the runs with the largest and smallest HFL
error. Their parameters are in samples_central.csv and their originals in central_parents.csv.

Usage: python compare_central.py   (run from this folder)
"""
from pathlib import Path

import pandas as pd

HERE = Path(__file__).parent
CENTRAL = HERE.parent / "central"


def accepted(r, tol):
    return (r.hfl_error_m.abs() <= tol) & (r.model_flooded_km2 >= 22) & (r.open_hit_rate >= 0.75)


def compare(name, original):
    rerun = pd.read_csv(CENTRAL / name / "results.csv").sort_values("run").reset_index(drop=True)
    d_area = rerun.max_flooded_area_km2 - original.max_flooded_area_km2
    d_hfl = rerun.hfl_error_m - original.hfl_error_m
    keep = accepted(rerun, 1.0)
    return {
        "variant": name, "runs": len(rerun),
        "abs_dA_lt_1": int((d_area.abs() < 1).sum()), "max_abs_dA_km2": d_area.abs().max(),
        "abs_de_lt_0.01": int((d_hfl.abs() < 0.01).sum()), "max_abs_de_m": d_hfl.abs().max(),
        "de_min_m": d_hfl.min(), "de_max_m": d_hfl.max(),
        "dropouts_1.0": int((~keep).sum()),
        "fail_hfl": int((~keep & (rerun.hfl_error_m.abs() > 1.0)).sum()),
        "fail_radar": int((~keep & ~((rerun.model_flooded_km2 >= 22) & (rerun.open_hit_rate >= 0.75))).sum()),
        "now_within_0.5": int(accepted(rerun, 0.5).sum()),
        "accepted_area_min": rerun.max_flooded_area_km2[keep].min(),
        "accepted_area_max": rerun.max_flooded_area_km2[keep].max(),
    }


if __name__ == "__main__":
    original = pd.read_csv(HERE / "central_parents.csv").sort_values("new_run").reset_index(drop=True)
    print(f"originals: area {original.max_flooded_area_km2.min():.1f}-"
          f"{original.max_flooded_area_km2.max():.1f} km2, HFL error "
          f"{original.hfl_error_m.min():+.2f} to {original.hfl_error_m.max():+.2f} m")
    table = pd.DataFrame([compare(n, original) for n in ["multiplicative", "lag_short", "lag_long"]])
    print(table.round(2).to_string(index=False))
    table.to_csv(CENTRAL / "summary.csv", index=False)
