"""How far do the HFL mark and the radar map narrow the predicted flood?

Re-analyses Latin hypercube ensembles (runs_with_behaviour.csv from analyse_runs.py):
empirical quantiles of the prior and conditioned ensembles, bootstrap intervals for
the range ends and for the narrowing, the narrowing under a grid of acceptance
thresholds, the HFL test and the radar test on their own, and an informal GLUE
likelihood weighting as an alternative to the yes/no acceptance rule.

The ensemble range of n accepted runs covers on average (n - 1)/(n + 1) of their
distribution, so the narrowing is reported for like-for-like intervals: the central
5-95% of both ensembles where enough runs are accepted, and the full range of the
accepted runs against the central interval of the prior with the same expected
coverage.

Usage:
    python -m research.periyar.uncertainty OUT_DIR RUNS_DIR [RUNS_DIR ...]
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

SAR_KM2 = 22.0
OUTPUT = "max_flooded_area_km2"
HFL_TOLERANCES = [0.25, 0.5, 0.75, 1.0, 1.5]
HIT_RATES = [0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9]


def accepted(runs, hfl_tol=0.5, hit_rate=0.75, use_hfl=True, use_radar=True, sar_km2=SAR_KM2):
    ok = runs["status"] == "ok"
    if use_hfl:
        ok &= runs["hfl_error_m"].abs() <= hfl_tol
    if use_radar:
        ok &= (runs["model_flooded_km2"] >= sar_km2) & (runs["open_hit_rate"] >= hit_rate)
    return ok


def intervals(prior, post):
    """The interval constructions compared in the paper, side by side."""
    n = len(post)
    coverage = (n - 1) / (n + 1)
    return {
        "prior_range": (prior.min(), prior.max()),
        "prior_5_95": tuple(np.quantile(prior, [0.05, 0.95])),
        "prior_like_for_like": tuple(np.quantile(prior, [(1 - coverage) / 2, (1 + coverage) / 2])),
        "accepted_range": (post.min(), post.max()),
        "accepted_5_95": tuple(np.quantile(post, [0.05, 0.95])),
    }


def narrowing(prior, post):
    """Narrowing of the accepted range against the prior interval of equal expected coverage."""
    n = len(post)
    if n < 2:
        return dict(n=n, low=np.nan, high=np.nan, coverage=np.nan, prior_low=np.nan,
                    prior_high=np.nan, narrowing=np.nan)
    coverage = (n - 1) / (n + 1)
    lo, hi = np.quantile(prior, [(1 - coverage) / 2, (1 + coverage) / 2])
    return dict(n=n, low=post.min(), high=post.max(), coverage=coverage, prior_low=lo,
                prior_high=hi, narrowing=1 - (post.max() - post.min()) / (hi - lo))


def weighted_quantiles(values, weights, q):
    order = np.argsort(values)
    v, w = np.asarray(values)[order], np.asarray(weights)[order]
    cdf = (np.cumsum(w) - 0.5 * w) / w.sum()
    return np.interp(q, cdf, v)


def glue_weights(runs, sigma_hfl, power=4):
    """Gaussian likelihood on the HFL residual times the radar open-land hit rate.

    Runs that flood less land than the radar's open water get zero weight.
    """
    w = np.exp(-0.5 * (runs["hfl_error_m"] / sigma_hfl) ** 2) * runs["open_hit_rate"] ** power
    return np.where(runs["model_flooded_km2"] >= SAR_KM2, w, 0.0)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("out_dir", type=Path)
    parser.add_argument("runs_dirs", type=Path, nargs="+")
    parser.add_argument("--bootstrap", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--hfl-tolerance", type=float, default=0.5, help="m, for the main rule")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    runs = pd.concat([pd.read_csv(d / "runs_with_behaviour.csv") for d in args.runs_dirs],
                     ignore_index=True)
    runs = runs[runs["status"] == "ok"].reset_index(drop=True)
    prior = runs[OUTPUT].to_numpy()
    main_tol = args.hfl_tolerance
    out = {"n_runs": len(runs), "hfl_tolerance_m": main_tol}

    # Empirical distributions.
    q = [0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1.0]
    rows = {"prior": np.quantile(prior, q)}
    for label, kw in [("hfl_only", dict(use_radar=False)), ("radar_only", dict(use_hfl=False)),
                      ("both", {})]:
        post = runs.loc[accepted(runs, hfl_tol=main_tol, **kw), OUTPUT].to_numpy()
        rows[label] = np.quantile(post, q) if len(post) else np.full(len(q), np.nan)
        out[f"n_{label}"] = int(len(post))
    quantiles = pd.DataFrame(rows, index=[f"q{int(x * 100):02d}" for x in q]).T
    quantiles.to_csv(args.out_dir / "quantiles_peak_area.csv")
    post = runs.loc[accepted(runs, hfl_tol=main_tol), OUTPUT].to_numpy()
    out["headline_5_95_narrowing"] = 1 - np.ptp(post) / np.diff(np.quantile(prior, [0.05, 0.95]))[0]
    out["like_for_like"] = narrowing(prior, post)

    # Bootstrap over the whole ensemble: resample runs, re-apply the rule.
    rng = np.random.default_rng(args.seed)
    boot = []
    for _ in range(args.bootstrap):
        sample = runs.iloc[rng.integers(0, len(runs), len(runs))]
        p = sample[OUTPUT].to_numpy()
        a = sample.loc[accepted(sample, hfl_tol=main_tol), OUTPUT].to_numpy()
        r = narrowing(p, a)
        boot.append([r["n"], r["low"], r["high"], r["narrowing"]])
    boot = pd.DataFrame(boot, columns=["n", "low", "high", "narrowing"])
    boot.to_csv(args.out_dir / "bootstrap.csv", index=False)
    out["bootstrap"] = {c: dict(zip(["q05", "q50", "q95"], boot[c].quantile([0.05, 0.5, 0.95])))
                        for c in boot}
    out["bootstrap"]["share_fewer_than_2"] = float((boot["n"] < 2).mean())

    # Threshold grid.
    grid = []
    for tol in HFL_TOLERANCES:
        for hit in HIT_RATES:
            a = runs.loc[accepted(runs, tol, hit), OUTPUT].to_numpy()
            grid.append({"hfl_tolerance_m": tol, "hit_rate": hit, **narrowing(prior, a)})
    pd.DataFrame(grid).to_csv(args.out_dir / "threshold_grid.csv", index=False)

    # Informal likelihood weighting instead of yes/no acceptance.
    glue = {}
    for sigma in (0.25, 0.5, 1.0):
        w = glue_weights(runs, sigma)
        lo, med, hi = weighted_quantiles(prior, w, [0.05, 0.5, 0.95])
        glue[f"sigma_{sigma}"] = dict(q05=lo, q50=med, q95=hi,
                                      effective_runs=float(w.sum() ** 2 / (w ** 2).sum()),
                                      narrowing_5_95=1 - (hi - lo) / np.diff(np.quantile(prior, [0.05, 0.95]))[0])
    out["glue"] = glue

    # Likelihood weighting with other powers of the hit rate, and with the HFL alone.
    width = np.diff(np.quantile(prior, [0.05, 0.95]))[0]
    weights = []
    for power in (0, 1, 4, 8):
        for sigma in (0.25, 0.5, 1.0):
            w = glue_weights(runs, sigma, power)
            lo, hi = weighted_quantiles(prior, w, [0.05, 0.95])
            weights.append({"hit_rate_power": power, "sigma_m": sigma, "q05": lo, "q95": hi,
                            "effective_runs": w.sum() ** 2 / (w ** 2).sum(),
                            "narrowing_5_95": 1 - (hi - lo) / width})
    pd.DataFrame(weights).to_csv(args.out_dir / "likelihood_weights.csv", index=False)

    # What each observation contributes, at two HFL tolerances, with every interval.
    rows = []
    for tol in (0.5, 1.0):
        for label, kw in [("hfl_only", dict(use_radar=False)), ("radar_only", dict(use_hfl=False)),
                          ("both", {})]:
            post = runs.loc[accepted(runs, hfl_tol=tol, **kw), OUTPUT].to_numpy()
            iv = intervals(prior, post)
            rows.append({"hfl_tolerance_m": tol, "test": label, "n": len(post),
                         **{f"{k}_{e}": v[i] for k, v in iv.items() for i, e in enumerate(("low", "high"))},
                         "narrowing_like_for_like": narrowing(prior, post)["narrowing"],
                         "narrowing_vs_prior_5_95": 1 - np.ptp(post) / width,
                         "narrowing_5_95_vs_5_95": 1 - np.diff(iv["accepted_5_95"])[0] / width})
    pd.DataFrame(rows).to_csv(args.out_dir / "contributions.csv", index=False)

    # The radar lower bound at the ends of its classification range.
    out["lower_bound_test"] = {str(km2): int(accepted(runs, hfl_tol=main_tol, sar_km2=km2).sum()) for km2 in (18.0, 22.0, 26.0)}
    out["min_model_flooded_at_pass_km2"] = float(runs["model_flooded_km2"].min())
    out["correlation_with_peak_area"] = {c: float(np.corrcoef(runs[c], runs[OUTPUT])[0, 1]) for c in
                                         ("open_hit_rate", "open_critical_success_index",
                                          "all_critical_success_index")}

    # HFL residuals.
    out["hfl_residuals"] = {
        "prior": dict(zip(["q05", "q50", "q95"], runs["hfl_error_m"].quantile([0.05, 0.5, 0.95]))),
        "accepted": runs.loc[accepted(runs, hfl_tol=main_tol), "hfl_error_m"].describe().to_dict(),
    }
    runs.assign(accepted=accepted(runs, hfl_tol=main_tol)).to_csv(args.out_dir / "runs.csv", index=False)
    text = json.dumps(out, indent=2, default=float)
    (args.out_dir / "summary.json").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
