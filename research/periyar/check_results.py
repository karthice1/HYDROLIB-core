"""Recompute the manuscript's principal numbers from the stored per-run results and compare.

No model runs are needed: every check reads CSV or JSON files under ``results/``. Each line
prints the manuscript value, the recomputed value and PASS or FAIL. The exit code is the
number of failed checks, so the script can run in CI.

Usage: python -m research.periyar.check_results
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

RESULTS = Path(__file__).parent / "results"
REV2 = RESULTS / "revision2"
SUPP = REV2 / "supplementary"
FIVE_FACTORS = ["n_floodplain", "channel_depth", "canopy_factor", "q_rating", "q_idukki"]

failures = 0


def check(label, expected, value, tol):
    global failures
    ok = abs(value - expected) <= tol
    failures += not ok
    print(f"{'PASS' if ok else 'FAIL'}  {label}: manuscript {expected}, recomputed {value:.3f}")


def accepted(runs, tol, hfl=True, radar=True):
    keep = np.ones(len(runs), bool)
    if hfl:
        keep &= (runs.hfl_error_m.abs() <= tol).values
    if radar:
        keep &= ((runs.model_flooded_km2 >= 22) & (runs.open_hit_rate >= 0.75)).values
    return keep


def contraction(runs, keep):
    """Conditional ensemble-range contraction C (Eq. 3) of peak flooded land."""
    area = runs.max_flooded_area_km2.values
    n = keep.sum()
    coverage = (n - 1) / (n + 1)
    low, high = np.quantile(area, [(1 - coverage) / 2, (1 + coverage) / 2])
    return 1 - (area[keep].max() - area[keep].min()) / (high - low)


def load_design(name):
    folder = RESULTS / name
    return pd.read_csv(folder / "results.csv").merge(pd.read_csv(folder / "samples.csv"), on="run")


def parents_in_pool(pool, parents):
    """Rows of the pooled ensemble with the parameter sets of a rerun design, in its order."""
    index = [np.abs(pool[FIVE_FACTORS].values - p[FIVE_FACTORS].values).sum(1).argmin()
             for _, p in parents.iterrows()]
    return pool.iloc[index].reset_index(drop=True)


def forcing_row(name, rerun, original, expected):
    """One row of a forcing-response table: small changes, largest changes and drop-outs."""
    small_area, max_area, small_hfl, max_hfl, out_05, out_10 = expected
    d_area = (rerun.max_flooded_area_km2 - original.max_flooded_area_km2).abs()
    d_hfl = (rerun.hfl_error_m - original.hfl_error_m).abs()
    base_05, base_10 = accepted(original, 0.5).sum(), accepted(original, 1.0).sum()
    keep_05 = (accepted(original, 0.5) & accepted(rerun, 0.5)).sum()
    keep_10 = (accepted(original, 1.0) & accepted(rerun, 1.0)).sum()
    check(f"{name} runs with |dA| < 1 km2", small_area, (d_area < 1).sum(), 0)
    check(f"{name} max |dA| (km2)", max_area, d_area.max(), 0.5)
    check(f"{name} runs with |de| < 0.01 m", small_hfl, (d_hfl < 0.01).sum(), 0)
    check(f"{name} max |de| (m)", max_hfl, d_hfl.max(), 0.005)
    check(f"{name} drop-outs at +-0.5 m", out_05, base_05 - keep_05, 0)
    check(f"{name} drop-outs at +-1.0 m", out_10, base_10 - keep_10, 0)


def main():
    print("Table 4: Morris mu* for peak flooded land (km2), v5")
    morris = pd.read_csv(RESULTS / "v5_morris_r20" / "indices_morris_r20.csv")
    area = morris[morris.metric == "max_flooded_area_km2"].set_index("factor").mu_star
    for factor, expected in [("channel_depth", 149), ("n_floodplain", 101), ("q_rating", 65)]:
        check(f"mu* {factor}", expected, area[factor], 0.5)

    print("\nTable 5: pooled five-factor ensemble (seed 42 + seed 43)")
    pooled = pd.read_csv(REV2 / "uncertainty_600_tol0.5" / "runs.csv")
    check("runs", 600, len(pooled), 0)
    for tol, n_exp, c_exp in [(1.0, 73, 0.62), (0.5, 23, 0.77), (1.5, 139, 0.53), (2.0, 217, 0.34)]:
        keep = accepted(pooled, tol)
        check(f"accepted at +-{tol} m", n_exp, keep.sum(), 0)
        check(f"C at +-{tol} m", c_exp, contraction(pooled, keep), 0.005)
    check("C, HFL only, +-1.0 m", 0.49, contraction(pooled, accepted(pooled, 1.0, radar=False)), 0.005)
    check("C, radar only", 0.20, contraction(pooled, accepted(pooled, 0, hfl=False)), 0.005)

    print("\nThe two Latin hypercube designs separately (+-1.0 m)")
    for name, n_exp, c_exp in [("v5_lhs200", 23, 0.63), ("v5_lhs400", 50, 0.61)]:
        runs = load_design(name)
        keep = accepted(runs, 1.0)
        check(f"{name} accepted", n_exp, keep.sum(), 0)
        check(f"{name} C", c_exp, contraction(runs, keep), 0.005)

    print("\nTable 11: eight-factor ensemble (seed 44)")
    lhs8 = pd.read_csv(REV2 / "lhs8_tol1.0" / "runs.csv")
    for tol, n_exp, c_exp in [(1.0, 22, 0.73), (1.5, 47, 0.54), (2.0, 73, 0.33)]:
        keep = accepted(lhs8, tol)
        check(f"accepted at +-{tol} m", n_exp, keep.sum(), 0)
        check(f"C at +-{tol} m", c_exp, contraction(lhs8, keep), 0.005)

    print("\nTable 11: second eight-factor ensemble (seed 46) and both designs pooled")
    lhs8b = pd.read_csv(REV2 / "lhs8b_tol1.0" / "runs.csv")
    for tol, n_exp, c_exp in [(1.0, 23, 0.69), (0.5, 7, 0.75)]:
        keep = accepted(lhs8b, tol)
        check(f"seed 46 accepted at +-{tol} m", n_exp, keep.sum(), 0)
        check(f"seed 46 C at +-{tol} m", c_exp, contraction(lhs8b, keep), 0.005)
    both = pd.concat([lhs8, lhs8b], ignore_index=True)
    check("pooled eight-factor accepted at +-1.0 m", 45, accepted(both, 1.0).sum(), 0)
    check("pooled eight-factor C at +-1.0 m", 0.72, contraction(both, accepted(both, 1.0)), 0.005)

    print("\nTable 11: fixed factors around the 23 strict runs")
    fixed = pd.read_csv(SUPP / "fixed" / "results.csv")
    check("fixed factors accepted at +-1.0 m", 57, accepted(fixed, 1.0).sum(), 0)
    check("fixed factors accepted at +-0.5 m", 30, accepted(fixed, 0.5).sum(), 0)

    print("\nTable 10: forcing reruns of the 23 strict runs, every cell")
    original = parents_in_pool(pooled, pd.read_csv(SUPP / "samples_accepted.csv"))
    table10 = {"multiplicative": (22, 37, 13, 0.76, 2, 1),
               "lag_short": (21, 37, 17, 0.29, 2, 0),
               "lag_long": (22, 62, 12, 0.64, 1, 1)}
    for name, expected in table10.items():
        rerun = pd.read_csv(SUPP / name / "results.csv").sort_values("run").reset_index(drop=True)
        forcing_row(name, rerun, original, expected)

    print("\nTable 10, lower block: forcing reruns of 14 stratified central runs")
    central = pd.read_csv(SUPP / "central_parents.csv").sort_values("new_run").reset_index(drop=True)
    check("central originals accepted at +-1.0 m", 14, accepted(central, 1.0).sum(), 0)
    check("central originals accepted at +-0.5 m", 0, accepted(central, 0.5).sum(), 0)
    table10_central = {"multiplicative": (12, 53, 8, 0.51, 0, 1),
                       "lag_short": (13, 41, 11, 0.51, 0, 1),
                       "lag_long": (13, 41, 8, 0.51, 0, 2)}
    for name, expected in table10_central.items():
        rerun = pd.read_csv(REV2 / "central" / name / "results.csv").sort_values("run").reset_index(drop=True)
        forcing_row(f"central {name}", rerun, central, expected)

    print("\nTable 7: contingency table, every cell")
    scores = pd.read_csv(REV2 / "confusion_median_accepted.csv").set_index(["run", "domain"])
    table7 = {("nominal", "all_land"): (278, 272, 4216, 17122, 0.51, 0.94, 0.058, 8.2),
              ("nominal", "open_land"): (85, 10, 993, 872, 0.89, 0.92, 0.078, 11.3),
              ("median_accepted", "all_land"): (241, 310, 2958, 18512, 0.44, 0.92, 0.069, 5.8),
              ("median_accepted", "open_land"): (74, 21, 695, 1175, 0.78, 0.90, 0.094, 8.1),
              ("bare_earth", "all_land"): (340, 200, 5789, 21407, 0.63, 0.94, 0.054, 11.3)}
    columns = ["hits", "misses", "false_alarms", "correct_negatives", "hit_rate",
               "false_alarm_ratio", "csi", "bias"]
    tolerances = [0, 0, 0, 0, 0.005, 0.005, 0.0005, 0.05]
    for key, expected in table7.items():
        for column, value, tol in zip(columns, expected, tolerances):
            check(f"{key[0]} {key[1]} {column}", value, scores.loc[key, column], tol)
    central = pooled[accepted(pooled, 1.0)]
    for label, column, low, high, tol in [
            ("73 accepted, all-land CSI", "all_critical_success_index", 0.059, 0.076, 0.0005),
            ("73 accepted, open-land hit rate", "open_hit_rate", 0.76, 0.88, 0.005),
            ("73 accepted, open-land FAR", "open_false_alarm_ratio", 0.90, 0.92, 0.005),
            ("73 accepted, open-land CSI", "open_critical_success_index", 0.080, 0.098, 0.0005)]:
        check(f"{label} min", low, central[column].min(), tol)
        check(f"{label} max", high, central[column].max(), tol)

    print("\nSection 4.6: drainage diagnosis, nominal run")
    drainage = json.loads((REV2 / "drainage_nominal" / "drainage_summary.json").read_text())
    check("flooded land at the pass (km2)", 179.8, drainage["flooded_land_at_pass_km2"], 0.05)
    check("water stored on land (Mcm)", 355, drainage["land_storage_at_pass_mcm"], 0.5)

    print(f"\n{failures} check(s) failed")
    return failures


if __name__ == "__main__":
    sys.exit(main())
