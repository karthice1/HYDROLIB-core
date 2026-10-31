"""Compare the reruns of the 23 strictly accepted runs (release form, travel times, repeats)
with their originals in the main 600-run ensemble.

Usage: python compare_reruns.py   (run from this folder)
"""
import numpy as np
import pandas as pd

FACTORS = ["n_floodplain", "channel_depth", "canopy_factor", "q_rating", "q_idukki"]
OPERATORS = ["hfl_err_nearest", "hfl_err_r300_median", "hfl_err_r750_median",
             "hfl_err_r1000_median", "hfl_err_r500_max", "hfl_err_r500_min"]


def accepted(r, tol):
    return (r.hfl_error_m.abs() <= tol) & (r.model_flooded_km2 >= 22) & (r.open_hit_rate >= 0.75)


main = pd.read_csv("../uncertainty_600_tol0.5/runs.csv")
parents = pd.read_csv("samples_accepted.csv")
index = [np.abs(main[FACTORS].values - p[FACTORS].values).sum(1).argmin() for _, p in parents.iterrows()]
original = main.iloc[index].reset_index(drop=True)

for name in ["multiplicative", "lag_short", "lag_long"]:
    rerun = pd.read_csv(f"{name}/results.csv").sort_values("run").reset_index(drop=True)
    d_area = rerun.max_flooded_area_km2 - original.max_flooded_area_km2
    d_hfl = rerun.hfl_error_m - original.hfl_error_m
    print(f"{name}: |d area| < 1 km2 in {(d_area.abs() < 1).sum()}/23, "
          f"|d HFL| < 0.01 m in {(d_hfl.abs() < 0.01).sum()}/23, d HFL {d_hfl.min():.2f} to {d_hfl.max():.2f} m")
    for tol in (0.5, 1.0, 1.5):
        a = accepted(rerun, tol)
        print(f"  +-{tol} m: {a.sum()} accepted, {rerun.max_flooded_area_km2[a].min():.0f}-"
              f"{rerun.max_flooded_area_km2[a].max():.0f} km2")

add = pd.read_csv("repeat_additive/results.csv").set_index("run").sort_index()
mul = pd.read_csv("repeat_multiplicative/results.csv").set_index("run").sort_index()
print("repeat additive reproduces original:",
      np.allclose(add.hfl_error_m.values, original.loc[add.index, "hfl_error_m"].values))
print("HFL change, multiplicative minus additive, by operator:")
print((mul[["hfl_error_m"] + OPERATORS] - add[["hfl_error_m"] + OPERATORS]).round(2).to_string())
