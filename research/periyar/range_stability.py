"""Are the conditioned results stable when a factor's sampled range is narrowed?

For each Latin hypercube factor, drops the runs in the lowest or the highest fifth
of its range (as if that range had been specified narrower), then recomputes the
prior 5-95% interval of peak flooded land, the accepted runs and the like-for-like
narrowing of uncertainty.py.

Usage:
    python -m research.periyar.range_stability OUT_CSV RUNS_CSV
"""

import argparse
from pathlib import Path

import pandas as pd

from .sensitivity import FACTORS
from .uncertainty import OUTPUT, accepted, narrowing

LHS_FACTORS = ["n_floodplain", "channel_depth", "canopy_factor", "q_rating", "q_idukki"]
TRIM = 0.2


def summary(label, runs, tol):
    prior = runs[OUTPUT].to_numpy()
    post = runs.loc[accepted(runs, hfl_tol=tol), OUTPUT].to_numpy()
    lo, hi = pd.Series(prior).quantile([0.05, 0.95])
    return {"case": label, "n_runs": len(runs), "prior_q05": lo, "prior_q95": hi,
            **{f"accepted_{k}": v for k, v in narrowing(prior, post).items()}}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("out_csv", type=Path)
    parser.add_argument("runs_csv", type=Path)
    parser.add_argument("--hfl-tolerance", type=float, default=0.5)
    args = parser.parse_args()
    tol = args.hfl_tolerance
    runs = pd.read_csv(args.runs_csv)
    runs = runs[runs["status"] == "ok"]
    bounds = {f.name: (f.lower, f.upper) for f in FACTORS}
    rows = [summary("full ranges", runs, tol)]
    for name in LHS_FACTORS:
        lo, hi = bounds[name]
        cut = TRIM * (hi - lo)
        rows.append(summary(f"{name} without lowest fifth", runs[runs[name] >= lo + cut], tol))
        rows.append(summary(f"{name} without highest fifth", runs[runs[name] <= hi - cut], tol))
    table = pd.DataFrame(rows)
    table.to_csv(args.out_csv, index=False)
    print(table[["case", "n_runs", "prior_q05", "prior_q95", "accepted_n", "accepted_low",
                 "accepted_high", "accepted_narrowing"]].round(2).to_string(index=False))


if __name__ == "__main__":
    main()
