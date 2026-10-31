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

    print("\nTable 10: eight-factor ensemble (seed 44)")
    lhs8 = pd.read_csv(REV2 / "lhs8_tol1.0" / "runs.csv")
    for tol, n_exp, c_exp in [(1.0, 22, 0.73), (1.5, 47, 0.54), (2.0, 73, 0.33)]:
        keep = accepted(lhs8, tol)
        check(f"accepted at +-{tol} m", n_exp, keep.sum(), 0)
        check(f"C at +-{tol} m", c_exp, contraction(lhs8, keep), 0.005)

    print("\nTable 10: fixed factors and forcing reruns of the 23 strict runs")
    fixed = pd.read_csv(SUPP / "fixed" / "results.csv")
    check("fixed factors accepted at +-1.0 m", 57, accepted(fixed, 1.0).sum(), 0)
    check("fixed factors accepted at +-0.5 m", 30, accepted(fixed, 0.5).sum(), 0)
    for name, n05, n10 in [("multiplicative", 21, 22), ("lag_short", 21, 23), ("lag_long", 22, 22)]:
        runs = pd.read_csv(SUPP / name / "results.csv")
        check(f"{name} accepted at +-0.5 m", n05, accepted(runs, 0.5).sum(), 0)
        check(f"{name} accepted at +-1.0 m", n10, accepted(runs, 1.0).sum(), 0)

    print("\nTable 7: contingency scores of the median accepted run")
    scores = pd.read_csv(REV2 / "confusion_median_accepted.csv").set_index(["run", "domain"])
    check("CSI, all land", 0.069, scores.loc[("median_accepted", "all_land"), "csi"], 0.0005)
    check("CSI, open land", 0.094, scores.loc[("median_accepted", "open_land"), "csi"], 0.0005)

    print("\nSection 4.6: drainage diagnosis, nominal run")
    drainage = json.loads((REV2 / "drainage_nominal" / "drainage_summary.json").read_text())
    check("flooded land at the pass (km2)", 179.8, drainage["flooded_land_at_pass_km2"], 0.05)
    check("water stored on land (Mcm)", 355, drainage["land_storage_at_pass_mcm"], 0.5)

    print(f"\n{failures} check(s) failed")
    return failures


if __name__ == "__main__":
    sys.exit(main())
