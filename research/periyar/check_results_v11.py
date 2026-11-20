"""Check the corrected-map values of the manuscript (v1.1-periyar) against the stored results.

Every radar score is recomputed here from the stored per-cell model states (results_v11/model_states) and the
corrected Sentinel-1 map (results_v11/radar). The HFL errors, peak areas and factor values come from the stored
per-run results of v1.0 (results/), which are unchanged. Acceptance, the conditional contraction C (Eq. 3), the
bootstrap, the threshold grid, the likelihood weighting and all spatial scores are recomputed with the code in this
file, which does not import the analysis modules it checks. No model is run.

Each line prints PASS or FAIL, the manuscript value and the recomputed value. Lines marked INFO are
recomputed values with no single manuscript number to compare. The exit code is the number of failed checks.

Usage: python -m research.periyar.check_results_v11        (from the repository root)
"""

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from scipy.stats import ks_2samp, spearmanr

HERE = Path(__file__).parent
R = HERE / "results"
SUPP = R / "revision2" / "supplementary"
CENTRAL = R / "revision2" / "central"
V11 = HERE / "results_v11"
STATES = V11 / "model_states"
RADAR = V11 / "radar"
MAP = RADAR / "sar_flood_20180821_corrected.tif"
MAP_SHA256 = "cda3ee7a00dbe42a5b30ba9a74c3ce121e75bb113f0503580e16cb4542691d8b"

OPEN = (30, 40, 60, 90)          # WorldCover grassland, cropland, bare, herbaceous wetland
TREE, BUILT = 10, 50
SAR_KM2 = 22.0                   # lower bound in the acceptance code (never binding, see check below)
AREA = "max_flooded_area_km2"
FIVE = ["n_floodplain", "channel_depth", "canopy_factor", "q_rating", "q_idukki"]
PRIOR = {"n_floodplain": (0.6, 1.5), "channel_depth": (2, 8), "canopy_factor": (0.3, 0.9),
         "q_rating": (0.85, 1.15), "q_idukki": (0.8, 1.2)}

failures = 0
n_checks = 0


def check(label, expected, value, tol):
    global failures, n_checks
    n_checks += 1
    ok = value is not None and np.isfinite(value) and abs(value - expected) <= tol + 1e-12
    failures += not ok
    shown = "n/a" if value is None or not np.isfinite(value) else f"{value:.4g}"
    print(f"{'PASS' if ok else 'FAIL'}  {label}: manuscript {expected}, recomputed {shown}")


def check_true(label, condition, detail=""):
    global failures, n_checks
    n_checks += 1
    failures += not condition
    print(f"{'PASS' if condition else 'FAIL'}  {label}{': ' + detail if detail else ''}")


def info(label, value):
    print(f"INFO  {label}: {value}")


# ----------------------------------------------------------------------------------------------- radar scoring
with rasterio.open(MAP) as _src:
    _MAP_DATA, _MAP_TF = _src.read(1), _src.transform


def sample(raster_data, transform, x, y):
    rows, cols = rasterio.transform.rowcol(transform, x, y)
    rows, cols = np.asarray(rows), np.asarray(cols)
    inside = (rows >= 0) & (rows < raster_data.shape[0]) & (cols >= 0) & (cols < raster_data.shape[1])
    out = np.full(len(rows), -1.0, dtype=raster_data.dtype)
    out[inside] = raster_data[rows[inside], cols[inside]]
    return out


def table(model, observed):
    h = int((model & observed).sum()); m = int((~model & observed).sum())
    f = int((model & ~observed).sum()); c = int((~model & ~observed).sum())
    return dict(hits=h, misses=m, false_alarms=f, correct_negatives=c, hit_rate=h / max(h + m, 1),
                far=f / max(h + f, 1), csi=h / max(h + m + f, 1), bias=(h + f) / max(h + m, 1))


def score(x, y, landcover, land, wet, raster=None, cut=0.5):
    data, tf = (_MAP_DATA, _MAP_TF) if raster is None else raster
    fr = sample(data, tf, x, y)
    land = land & (fr >= 0)
    observed = fr > cut
    open_land = land & np.isin(landcover, OPEN)
    out = {"all": table(wet[land], observed[land]), "open": table(wet[open_land], observed[open_land])}
    fa = wet & ~observed & land
    out["fa_nonopen_share"] = float((fa & ~open_land).sum() / max(fa.sum(), 1))
    out["fa_tree_share"] = float((fa & (landcover == TREE)).sum() / max(fa.sum(), 1))
    out["radar_cells_on_land"] = int((observed & land).sum())
    return out


GRID = np.load(STATES / "grid.npz")


def group_scores(group, raster=None, cut=0.5):
    """Corrected radar scores of every stored state of a group, by run number."""
    s = np.load(STATES / f"states_{group}.npz")
    n = int(s["n_cells"])
    rows = []
    for name, land, wet in zip(s["run"], s["land"], s["model_wet"]):
        land = np.unpackbits(land)[:n].astype(bool); wet = np.unpackbits(wet)[:n].astype(bool)
        sc = score(GRID["x"], GRID["y"], GRID["landcover"], land, wet, raster, cut)
        rows.append(dict(run=int(str(name)[4:]), open_hit_rate=sc["open"]["hit_rate"],
                         open_false_alarm_ratio=sc["open"]["far"], open_critical_success_index=sc["open"]["csi"],
                         open_bias=sc["open"]["bias"], all_critical_success_index=sc["all"]["csi"],
                         all_far=sc["all"]["far"], all_bias=sc["all"]["bias"], all_hit_rate=sc["all"]["hit_rate"],
                         fa_nonopen_share=sc["fa_nonopen_share"], fa_tree_share=sc["fa_tree_share"],
                         _all=sc["all"], _open=sc["open"]))
    return pd.DataFrame(rows).set_index("run")


RADAR_COLS = ["open_hit_rate", "open_false_alarm_ratio", "open_critical_success_index", "all_critical_success_index"]


def with_corrected_scores(stored, scores):
    """Stored v1.0 results with the four radar columns replaced by corrected scores (NaN where no state exists)."""
    out = stored.drop(columns=[c for c in RADAR_COLS if c in stored]).set_index("run")
    out = out.join(scores.drop(columns=["_all", "_open"]), how="left")
    return out.reset_index()


# ----------------------------------------------------------------------------------------------- conditioning
def accepted(runs, tol, hfl=True, radar=True, hit=0.75, bound=SAR_KM2, err="hfl_error_m"):
    keep = (runs["status"] == "ok").values.copy()
    if hfl:
        keep &= (runs[err].abs() <= tol).values
    if radar:
        keep &= ((runs.model_flooded_km2 >= bound) & (runs.open_hit_rate >= hit)).fillna(False).values
    return keep


def contraction(prior, post):
    n = len(post)
    if n < 2:
        return dict(n=n, C=np.nan)
    cov = (n - 1) / (n + 1)
    lo, hi = np.quantile(prior, [(1 - cov) / 2, (1 + cov) / 2])
    return dict(n=n, low=post.min(), high=post.max(), prior_low=lo, prior_high=hi,
                C=1 - (post.max() - post.min()) / (hi - lo))


def C_of(runs, keep):
    return contraction(runs[AREA].to_numpy(), runs.loc[keep, AREA].to_numpy())


def bootstrap(runs, tol, n_boot=5000, seed=1):
    rng = np.random.default_rng(seed)
    cs, ns = [], []
    for _ in range(n_boot):
        s = runs.iloc[rng.integers(0, len(runs), len(runs))].reset_index(drop=True)
        c = C_of(s, accepted(s, tol))
        if c["n"] >= 2:
            cs.append(c["C"]); ns.append(c["n"])
    return np.quantile(cs, [0.05, 0.95]), np.quantile(ns, [0.05, 0.95])


def weighted_quantiles(values, weights, q):
    order = np.argsort(values)
    v, w = np.asarray(values)[order], np.asarray(weights)[order]
    cdf = (np.cumsum(w) - 0.5 * w) / w.sum()
    return np.interp(q, cdf, v)


def likelihood(runs, s, p):
    w = np.exp(-0.5 * (runs.hfl_error_m / s) ** 2) * runs.open_hit_rate ** p
    w = np.where(runs.model_flooded_km2 >= SAR_KM2, w, 0.0)
    prior = runs[AREA].to_numpy()
    lo, hi = weighted_quantiles(prior, w, [0.05, 0.95])
    width = np.diff(np.quantile(prior, [0.05, 0.95]))[0]
    return 1 - (hi - lo) / width, w.sum() ** 2 / (w ** 2).sum(), lo, hi


def load(folder):
    return pd.read_csv(folder / "results.csv").merge(pd.read_csv(folder / "samples.csv"), on="run", suffixes=("", "_s"))


def main():
    print("Corrected radar map")
    sha = hashlib.sha256(MAP.read_bytes()).hexdigest()
    check_true("corrected map SHA-256", sha == MAP_SHA256, sha)
    qc = json.loads((RADAR / "georef_qc.json").read_text())
    d = qc["20180821"]
    check("descending shift, metres west", 820, -d["image_displacement_m"]["east"], 0)
    check("descending shift, metres north", 180, d["image_displacement_m"]["north"], 0)
    check("descending displacement (m)", 840, d["image_displacement_m"]["magnitude"], 0.5)
    check("look direction (deg)", 282, d["look_direction_deg"], 0.5)
    check("residual against WorldCover built-up (m, largest component)", 20,
          max(abs(v) for v in d["verification"]["corrected"]["builtup_best_shift_m"]), 0)
    check("ascending displacement (m)", 150, qc["20180821_asc"]["image_displacement_m"]["magnitude"], 10)
    check("morning Otsu threshold (dB)", -13.75, qc["otsu_threshold_db"]["morning_corrected"], 0.005)
    check("evening Otsu threshold (dB)", -13.31, qc["otsu_threshold_db"]["evening_corrected"], 0.005)
    info("gcp heights descending (m)", d["annotation"]["gcp_height_m"])
    info("gcp heights ascending (m)", qc["20180821_asc"]["annotation"]["gcp_height_m"])
    info("map hash of the distributed file recorded in georef_qc.json",
         qc["maps"]["sar_flood_20180821_corrected"]["sha256"])
    info("map hash as written by sar_georef_correction.py (same pixels, different DEFLATE byte stream)",
         qc["maps"]["sar_flood_20180821_corrected"]["sha256_as_written_by_script"])

    print("\nStored structural states (Table 9, Table 10)")
    st = {}
    for name in ("nominal", "min_terrain", "full_canopy", "wide_channel", "grid_100m", "tide_x1.5", "tide_shift_6h"):
        z = np.load(R / "revision" / "structural" / f"{name}_radar_state.npz")
        st[name] = score(z["x"], z["y"], z["landcover"], z["land"], z["model_wet"])
    nom = st["nominal"]
    check("radar-flooded land cells (fraction > 0.5)", 571, nom["radar_cells_on_land"], 0)
    check("radar-flooded land (km2)", 23, nom["radar_cells_on_land"] * 0.04, 0.5)
    table9 = {("nominal", "all"): (508, 63, 3986, 17182, 0.89, 0.89, 0.111, 7.9),
              ("nominal", "open"): (455, 55, 623, 827, 0.89, 0.58, 0.402, 2.1),
              ("min_terrain", "all"): (544, 20, 5585, 21440, 0.96, 0.91, 0.088, 10.9)}
    keys = ["hits", "misses", "false_alarms", "correct_negatives", "hit_rate", "far", "csi", "bias"]
    tols = [0, 0, 0, 0, 0.005, 0.005, 0.0005, 0.05]
    for (run, dom), exp in table9.items():
        for k, e, t in zip(keys, exp, tols):
            check(f"Table 9 {run} {dom} {k}", e, st[run][dom][k], t)
    table10 = {"nominal": (0.111, 0.402), "min_terrain": (0.088, 0.424), "full_canopy": (0.100, 0.388),
               "wide_channel": (0.127, 0.395), "grid_100m": (0.116, 0.304), "tide_x1.5": (0.111, 0.401),
               "tide_shift_6h": (0.110, 0.401)}
    for run, (a, o) in table10.items():
        check(f"Table 10 {run} all-land CSI", a, st[run]["all"]["csi"], 0.0005)
        check(f"Table 10 {run} open-land CSI", o, st[run]["open"]["csi"], 0.0005)
    check_true("no structural variant above all-land CSI 0.13", max(v["all"]["csi"] for v in st.values()) <= 0.13)

    print("\nMain ensemble: corrected radar scores recomputed from 600 per-cell states")
    designs = {}
    for name in ("lhs200", "lhs400"):
        stored = load(R / f"v5_{name}")
        sc = group_scores(name)
        check(f"{name}: states available", len(stored), sc.index.isin(stored.run).sum(), 0)
        dd = with_corrected_scores(stored, sc)
        dd["design"] = name
        designs[name] = dd
    m = pd.concat([designs["lhs200"], designs["lhs400"]], ignore_index=True)
    check("pooled runs", 600, len(m), 0)
    check_true("every run has a corrected score", m.open_hit_rate.notna().all())
    ref = pd.read_csv(V11 / "rescoring" / "runs_corrected.csv")
    diff = max(np.nanmax(np.abs(m[c].to_numpy() - ref[c].to_numpy())) for c in RADAR_COLS)
    check("agreement with runs_corrected.csv (max abs difference)", 0, diff, 1e-9)
    A = m[AREA].to_numpy()

    print("\nTable 6: intervals and C (Eq. 3)")
    check("prior minimum (km2)", 88.8, A.min(), 0.05)
    check("prior maximum (km2)", 363.6, A.max(), 0.05)
    q = np.quantile(A, [0.05, 0.5, 0.95])
    check("prior 5% (km2)", 107.3, q[0], 0.05); check("prior median (km2)", 220, q[1], 0.5)
    check("prior 95% (km2)", 329.2, q[2], 0.05)
    width = q[2] - q[0]
    table6 = {(1.0, "hfl"): (131, 91.7, 353.7, 110.5, 243.4, 124.2, 212.2, 0.49, 0.40),
              (1.0, "radar"): (488, 89.1, 360.6, 143.4, 363.6, 162.2, 334.2, 0.19, 0.01),
              (1.0, "both"): (69, 96.6, 350.3, 143.4, 243.4, 152.2, 237.8, 0.61, 0.55),
              (0.5, "hfl"): (55, 97.7, 349.0, 113.9, 198.2, 131.9, 179.4, 0.66, 0.62),
              (0.5, "radar"): (488, 89.1, 360.6, 143.4, 363.6, 162.2, 334.2, 0.19, 0.01),
              (0.5, "both"): (20, 106.8, 331.4, 145.8, 198.2, 149.2, 182.6, 0.77, 0.76),
              (1.5, "hfl"): (217, 90.2, 357.3, 102.3, 268.9, 118.9, 233.2, 0.38, 0.25),
              (1.5, "both"): (135, 91.5, 353.7, 143.4, 268.9, 154.2, 240.6, 0.52, 0.43),
              (2.0, "hfl"): (306, 89.3, 360.1, 101.2, 323.2, 110.7, 270.0, 0.18, 0.00),
              (2.0, "both"): (213, 90.2, 357.1, 143.4, 323.2, 157.3, 278.3, 0.33, 0.19)}
    labels = ["n", "like-for-like low", "like-for-like high", "accepted min", "accepted max", "accepted 5%",
              "accepted 95%", "C", "C vs prior 5-95%"]
    for (tol, test), exp in table6.items():
        keep = accepted(m, tol, hfl=test != "radar", radar=test != "hfl")
        c = C_of(m, keep); post = A[keep]
        got = [c["n"], c["prior_low"], c["prior_high"], c["low"], c["high"], *np.quantile(post, [0.05, 0.95]),
               c["C"], 1 - np.ptp(post) / width]
        for lab, e, g, t in zip(labels, exp, got, [0, .05, .05, .05, .05, .05, .05, .005, .005]):
            check(f"+-{tol} m {test}: {lab}", e, g, t)
    keep10, keep05 = accepted(m, 1.0), accepted(m, 0.5)
    hfl10 = (m.hfl_error_m.abs() <= 1.0).values
    radar_ok = accepted(m, 0, hfl=False)
    info("C at +-1.0 m, +-0.5 m, +-1.5 m, +-2.0 m (unrounded)",
         [round(C_of(m, accepted(m, t))["C"], 4) for t in (1.0, 0.5, 1.5, 2.0)])
    check("coverage of like-for-like interval at n = 69 (%)", 97, 100 * 68 / 70, 0.5)
    check("accepted HFL error min (m)", -0.95, m.hfl_error_m[keep10].min(), 0.005)
    check("accepted HFL error max (m)", 1.00, m.hfl_error_m[keep10].max(), 0.005)
    hq = np.quantile(m.hfl_error_m, [0.05, 0.5, 0.95])
    check("prior HFL residual median (m)", 1.9, hq[1], 0.05)
    check("prior HFL residual 5% (m)", -1.45, hq[0], 0.005); check("prior HFL residual 95% (m)", 3.58, hq[2], 0.005)
    vol = m.max_flood_volume_mcm
    check("accepted volume min +-1.0 (Mm3)", 396, vol[keep10].min(), 0.5)
    check("accepted volume max +-1.0 (Mm3)", 700, vol[keep10].max(), 0.5)
    check("accepted volume min +-0.5 (Mm3)", 411, vol[keep05].min(), 0.5)
    check("accepted volume max +-0.5 (Mm3)", 558, vol[keep05].max(), 0.5)
    vq = np.quantile(vol, [0.05, 0.95])
    check("prior volume 5% (Mm3)", 343, vq[0], 0.5); check("prior volume 95% (Mm3)", 1121, vq[1], 0.5)
    check("+-2.0 m accepted min (km2)", 143, A[accepted(m, 2.0)].min(), 0.5)
    check("+-2.0 m accepted max (km2)", 323, A[accepted(m, 2.0)].max(), 0.5)

    print("\nBootstrap (5,000 resamples, seed 1)")
    (c05, c95), (n05, n95) = bootstrap(m, 1.0)
    check("+-1.0 m bootstrap C 5%", 0.60, c05, 0.005); check("+-1.0 m bootstrap C 95%", 0.64, c95, 0.005)
    check("+-1.0 m bootstrap n 5%", 56, n05, 0.5); check("+-1.0 m bootstrap n 95%", 83, n95, 0.5)
    (c05, c95), (n05, n95) = bootstrap(m, 0.5)
    check("+-0.5 m bootstrap C 5%", 0.75, c05, 0.005); check("+-0.5 m bootstrap C 95%", 0.86, c95, 0.005)
    check("+-0.5 m bootstrap n 5%", 13, n05, 0.5); check("+-0.5 m bootstrap n 95%", 27, n95, 0.5)

    print("\nThe two designs separately")
    for name, tol, n_e, lo_e, hi_e, c_e in [("lhs200", 1.0, 22, 157, 240, 0.63), ("lhs400", 1.0, 47, 143, 243, 0.60),
                                            ("lhs200", 0.5, 7, 158, 182, 0.84), ("lhs400", 0.5, 13, None, None, 0.74)]:
        dd = designs[name]; k = accepted(dd, tol); c = C_of(dd, k)
        check(f"{name} +-{tol} accepted", n_e, c["n"], 0)
        check(f"{name} +-{tol} C", c_e, c["C"], 0.005)
        if lo_e is not None:
            check(f"{name} +-{tol} accepted min (km2)", lo_e, c["low"], 0.5)
            check(f"{name} +-{tol} accepted max (km2)", hi_e, c["high"], 0.5)
    for name, rng_d, rng_n in [("lhs200", (2.4, 7.8), (0.61, 1.31)), ("lhs400", (2.1, 7.8), (0.60, 1.42))]:
        dd = designs[name]; k = accepted(dd, 1.0)
        check(f"{name} accepted depth min", rng_d[0], dd.channel_depth[k].min(), 0.05)
        check(f"{name} accepted depth max", rng_d[1], dd.channel_depth[k].max(), 0.05)
        check(f"{name} accepted roughness min", rng_n[0], dd.n_floodplain[k].min(), 0.005)
        check(f"{name} accepted roughness max", rng_n[1], dd.n_floodplain[k].max(), 0.005)

    print("\nWhat each observation contributes; radar rejects")
    rej = ~radar_ok
    check("radar rejects", 112, rej.sum(), 0)
    check("radar rejects peak min (km2)", 89, A[rej].min(), 0.5); check("radar rejects peak max (km2)", 175, A[rej].max(), 0.5)
    h05 = (m.hfl_error_m.abs() <= 0.5).values
    check("runs passing HFL at +-0.5 m", 55, h05.sum(), 0)
    check("of which rejected by radar", 35, (rej & h05).sum(), 0)
    check_true("radar test rejects no run with peak land above 175 km2", bool(radar_ok[A > 175].all()))
    info("smallest radar-accepted peak / largest radar-rejected peak (km2)", (A[radar_ok].min(), A[rej].max()))
    check("Pearson r, open hit rate vs peak land", 0.97, np.corrcoef(m.open_hit_rate, A)[0, 1], 0.005)
    check("Pearson r, open CSI vs peak land", -0.76, np.corrcoef(m.open_critical_success_index, A)[0, 1], 0.005)
    check("Pearson r, all CSI vs peak land", -0.96, np.corrcoef(m.all_critical_success_index, A)[0, 1], 0.005)
    check("driest run, land flooded at pass (km2)", 77, m.model_flooded_km2.min(), 0.5)
    for b in (18.0, 23.0, 26.0):
        check_true(f"lower bound {b} km2 accepts the same 69 runs", (accepted(m, 1.0, bound=b) == keep10).all())
    check_true("acceptance code constant (22 km2) equivalent to 23 km2 on these runs",
               (accepted(m, 1.0, bound=22.0) == accepted(m, 1.0, bound=22.84)).all())

    print("\nTable 7 and correlations")
    table7 = {"n_floodplain": (0.42, 0.46, 0.39, 0.77), "channel_depth": (0.32, 0.50, 0.17, 0.81),
              "canopy_factor": (0.27, 0.15, 0.35, 0.80), "q_rating": (0.29, 0.33, 0.19, 0.91),
              "q_idukki": (0.15, 0.10, 0.19, 0.78)}
    for f, (dh, dr, db, share) in table7.items():
        v = m[f]; lo, hi = PRIOR[f]
        check(f"D HFL only {f}", dh, ks_2samp(v[hfl10], v[~hfl10]).statistic, 0.005)
        check(f"D radar only {f}", dr, ks_2samp(v[radar_ok], v[~radar_ok]).statistic, 0.005)
        check(f"D both {f}", db, ks_2samp(v[keep10], v[~keep10]).statistic, 0.005)
        check(f"accepted 5-95% share {f}", share, (np.quantile(v[keep10], .95) - np.quantile(v[keep10], .05)) / (hi - lo), 0.005)
    check("p value, D both, channel depth", 0.06, ks_2samp(m.channel_depth[keep10], m.channel_depth[~keep10]).pvalue, 0.005)
    n1, n2 = keep10.sum(), (~keep10).sum()
    check("5% critical D (both)", 0.17, 1.358 * np.sqrt((n1 + n2) / (n1 * n2)), 0.005)
    sp = lambda a, b, k=None: spearmanr(m[a][k] if k is not None else m[a], m[b][k] if k is not None else m[b])[0]
    check("rho depth-roughness, prior", 0.01, sp("channel_depth", "n_floodplain"), 0.005)
    check("rho depth-roughness, 69 accepted", 0.54, sp("channel_depth", "n_floodplain", keep10), 0.005)
    check("rho depth-roughness, 20 accepted", 0.42, sp("channel_depth", "n_floodplain", keep05), 0.005)
    check("rho depth-canopy, 69 accepted", 0.50, sp("channel_depth", "canopy_factor", keep10), 0.005)
    check("rho depth-canopy, 20 accepted", 0.91, sp("channel_depth", "canopy_factor", keep05), 0.005)
    check("rho depth-canopy, HFL only", 0.22, sp("channel_depth", "canopy_factor", hfl10), 0.005)
    check("rho depth-HFL error", -0.58, sp("channel_depth", "hfl_error_m"), 0.005)
    check("rho roughness-HFL error", 0.53, sp("n_floodplain", "hfl_error_m"), 0.005)
    check("rho depth-hit rate", -0.79, sp("channel_depth", "open_hit_rate"), 0.005)
    check("rho canopy-hit rate", 0.17, sp("canopy_factor", "open_hit_rate"), 0.005)
    check("rho gauge level-peak land", 0.86, sp("neeleswaram_peak_level_m", AREA), 0.005)
    check_true("|rho| >= 0.98 for all-land CSI, open hit rate and open FAR vs peak land",
               min(abs(sp(c, AREA)) for c in ("all_critical_success_index", "open_hit_rate", "open_false_alarm_ratio")) >= 0.975,
               str([round(sp(c, AREA), 3) for c in ("all_critical_success_index", "open_hit_rate", "open_false_alarm_ratio")]))
    for f, lo_e, hi_e, med_e in [("channel_depth", 2.1, 7.8, None), ("n_floodplain", 0.60, 1.42, None),
                                 ("canopy_factor", 0.32, 0.90, 0.75), ("q_rating", 0.85, 1.14, 0.96)]:
        t = 0.05 if f == "channel_depth" else 0.005
        check(f"69 accepted {f} min", lo_e, m[f][keep10].min(), t); check(f"69 accepted {f} max", hi_e, m[f][keep10].max(), t)
        if med_e is not None:
            check(f"69 accepted {f} median", med_e, m[f][keep10].median(), 0.005)

    print("\nTable 8 and threshold grid")
    grid = {(1.0, 0.5): 0.49, (1.0, 0.6): 0.49, (1.0, 0.7): 0.56, (1.0, 0.8): 0.68, (1.0, 0.85): 0.87}
    for (tol, hit), e in grid.items():
        check(f"C at +-{tol} m, hit rate {hit}", e, C_of(m, accepted(m, tol, hit=hit))["C"], 0.005)
    check("runs at +-1.0 m, hit rate 0.85", 6, accepted(m, 1.0, hit=0.85).sum(), 0)
    check("runs at +-1.0 m, hit rate 0.5 (= HFL only)", 131, accepted(m, 1.0, hit=0.5).sum(), 0)
    g05 = [C_of(m, accepted(m, 0.5, hit=h))["C"] for h in (0.5, 0.6, 0.7, 0.75, 0.8)]
    check("C at +-0.5 m, hit rates 0.5-0.8, min", 0.66, min(g05), 0.005)
    check("C at +-0.5 m, hit rates 0.5-0.8, max", 0.84, max(g05), 0.005)
    check("runs at +-0.5 m, hit rate 0.85", 0, accepted(m, 0.5, hit=0.85).sum(), 0)
    check("runs at +-1.0 m, hit rate 0.9", 0, accepted(m, 1.0, hit=0.9).sum(), 0)
    for s, p, nar, neff, lo_e, hi_e in [(0.25, 4, 0.80, 51, 135, 180), (0.5, 4, 0.62, 121, 134, 219),
                                        (1.0, 4, 0.33, 302, 135, 283)]:
        a, b, lo, hi = likelihood(m, s, p)
        check(f"likelihood s={s} p={p} narrowing", nar, a, 0.005)
        check(f"likelihood s={s} p={p} effective runs", neff, b, 0.5)
        check(f"likelihood s={s} p={p} weighted 5%", lo_e, lo, 0.5); check(f"likelihood s={s} p={p} weighted 95%", hi_e, hi, 0.5)
    for p, e in [(0, 0.66), (1, 0.66), (8, 0.54)]:
        check(f"likelihood s=0.5 p={p} narrowing", e, likelihood(m, 0.5, p)[0], 0.005)
    effs = [likelihood(m, 0.5, p)[1] for p in (0, 1, 4, 8)]
    check("likelihood s=0.5, p=0-8, effective runs min", 115, min(effs), 0.5)
    check("likelihood s=0.5, p=0-8, effective runs max", 122, max(effs), 0.5)
    cuts = []
    for f in FIVE:
        lo, hi = PRIOR[f]
        for side, keep in (("low", m[f] >= lo + 0.2 * (hi - lo)), ("high", m[f] <= hi - 0.2 * (hi - lo))):
            d = m[keep.values].reset_index(drop=True); k = accepted(d, 1.0); c = C_of(d, k); c.update(f=f, side=side)
            cuts.append(c)
    check("range cut: accepted min", 37, min(c["n"] for c in cuts), 0)
    check("range cut: accepted max", 66, max(c["n"] for c in cuts), 0)
    others = [c["C"] for c in cuts if not (c["f"] == "canopy_factor" and c["side"] == "high")]
    can = [c for c in cuts if c["f"] == "canopy_factor" and c["side"] == "high"][0]
    check("range cut: C min (all cases)", 0.56, min(others), 0.005)
    check("range cut: C max (other than high canopy)", 0.61, max(others), 0.005)
    check("range cut: high canopy removed, C", 0.80, can["C"], 0.005)
    check("range cut: high canopy removed, upper end (km2)", 190, can["high"], 0.5)
    rng = np.random.default_rng(1); cs, ns = [], []
    for _ in range(5000):
        s = m.iloc[rng.choice(600, 200, replace=False)].reset_index(drop=True); c = C_of(s, accepted(s, 1.0))
        if c["n"] >= 2:
            cs.append(c["C"]); ns.append(c["n"])
    cs = np.array(cs)
    check("200-run subsets C 5%", 0.54, np.quantile(cs, .05), 0.005)
    check("200-run subsets C 95%", 0.72, np.quantile(cs, .95), 0.005)
    check("200-run subsets n 5%", 17, np.quantile(ns, .05), 0.5); check("200-run subsets n 95%", 29, np.quantile(ns, .95), 0.5)

    print("\nSpatial agreement (Section 4.4, Table 9)")
    acc = m[keep10]
    for lab, col, lo_e, hi_e, t in [("hit rate", "all_hit_rate", 0.75, 0.87, 0.005), ("FAR", "all_far", 0.84, 0.88, 0.005),
                                    ("CSI", "all_critical_success_index", 0.116, 0.154, 0.0005), ("bias", "all_bias", 4.6, 7.3, 0.05),
                                    ("open hit rate", "open_hit_rate", 0.75, 0.88, 0.005),
                                    ("open FAR", "open_false_alarm_ratio", 0.44, 0.57, 0.005),
                                    ("open CSI", "open_critical_success_index", 0.395, 0.469, 0.0005), ("open bias", "open_bias", 1.4, 2.0, 0.05)]:
        check(f"69 accepted all-land {lab} min" if not lab.startswith("open") else f"69 accepted {lab} min", lo_e, acc[col].min(), t)
        check(f"69 accepted all-land {lab} max" if not lab.startswith("open") else f"69 accepted {lab} max", hi_e, acc[col].max(), t)
    check("69 accepted: false alarms outside open land, min", 0.83, acc.fa_nonopen_share.min(), 0.005)
    check("69 accepted: false alarms outside open land, max", 0.87, acc.fa_nonopen_share.max(), 0.005)
    check("69 accepted: false alarms under tree cover, min", 0.78, acc.fa_tree_share.min(), 0.005)
    check("69 accepted: false alarms under tree cover, max", 0.83, acc.fa_tree_share.max(), 0.005)
    check("600 runs: all-land FAR min", 0.82, m.all_far.min(), 0.005)
    check("600 runs: all-land FAR max", 0.91, m.all_far.max(), 0.005)
    check("600 runs: all-land CSI max", 0.17, m.all_critical_success_index.max(), 0.005)
    check_true("every run floods at least three times the radar's area (all-land bias >= 3)", m.all_bias.min() >= 3,
               f"minimum bias {m.all_bias.min():.2f}")
    check("600 runs: land flooded at the pass, min (km2)", 77, m.model_flooded_km2.min(), 0.5)
    check("600 runs: land flooded at the pass, max (km2)", 250, m.model_flooded_km2.max(), 0.5)
    check("median peak land of the 69 accepted (km2)", 177, np.median(A[keep10]), 0.5)
    k05 = m[keep05].sort_values(AREA)
    med = k05.iloc[len(k05) // 2 - 1: len(k05) // 2 + 1]
    info("the two central runs of the 20 (design, run, peak)", med[["design", "run", AREA]].values.tolist())
    r70 = m[(m.design == "lhs400") & (m.run == 70)].iloc[0]
    check_true("lhs400 run 70 is accepted at +-0.5 m and +-1.0 m", bool(keep05[r70.name]) and bool(keep10[r70.name]))
    check("run 70 peak land (km2)", 160, r70[AREA], 0.5); check("run 70 HFL error (m)", 0.40, r70.hfl_error_m, 0.005)
    s70 = group_scores("lhs400").loc[70]
    for k, e, t in zip(keys, (430, 139, 2604, 18436, 0.76, 0.86, 0.136, 5.3), tols):
        check(f"Table 9 median accepted all {k}", e, s70["_all"][k], t)
    for k, e, t in zip(keys, (383, 126, 369, 1074, 0.75, 0.49, 0.436, 1.5), tols):
        check(f"Table 9 median accepted open {k}", e, s70["_open"][k], t)
    check_true("every one of the 20 runs accepted at +-0.5 m is accepted at +-1.0 m", bool((keep10 | ~keep05).all()))

    print("\nRadar map variants (Section 4.4)")
    V = RADAR / "sar_variants"
    vals, areas = [], {}
    st200 = {g: np.load(STATES / f"states_{g}.npz") for g in ("lhs200", "lhs400")}
    for name in ("base", "thr_m12.25", "thr_m15.25", "drop2", "drop4"):
        with rasterio.open(V / f"{name}.tif") as src:
            rd = (src.read(1), src.transform)
        areas[name] = int(score(*[np.load(R / "revision" / "structural" / "nominal_radar_state.npz")[k] for k in ("x", "y", "landcover", "land", "model_wet")], rd)["radar_cells_on_land"]) * 0.04
        for cut in ((0.5, 0.3) if name == "base" else (0.5,)):
            for g in ("lhs200", "lhs400"):
                sc = group_scores(g, rd, cut)
                ids = m[(m.design == g) & keep10].run
                vals += sc.loc[ids, "all_critical_success_index"].tolist()
    check("variants: all-land CSI of the 69 accepted, min", 0.09, min(vals), 0.005)
    check("variants: all-land CSI of the 69 accepted, max", 0.22, max(vals), 0.005)
    check("variant open-water area on land, min (km2)", 18, min(areas.values()), 0.5)
    check("variant open-water area on land, max (km2) [proposed correction; DOCX says 28, an all-cell count]", 27, max(areas.values()), 0.5)
    info("variant radar land (km2)", {k: round(v, 1) for k, v in areas.items()})

    print("\nFig. 4 runs (7 of the first design)")
    d200 = designs["lhs200"]; k2 = accepted(d200, 0.5)
    ts = [26, 61, 64, 86, 93, 119, 188]
    still = d200.run[k2].isin(ts).sum(); new = sorted(set(d200.run[k2]) - set(ts))
    check("time-series runs still accepted at +-0.5 m", 6, still, 0)
    check("runs newly accepted (no hourly output)", 1, len(new), 0)
    drop = [r for r in ts if r not in set(d200.run[k2])]
    check("hit rate of the time-series run that drops out", 0.749, float(d200.set_index("run").loc[drop[0], "open_hit_rate"]), 0.0005)

    print("\nEight-factor designs (radar scores exist only for |HFL error| <= 1.0 m)")
    eight = {}
    for name in ("lhs8", "lhs8b"):
        stored = load(SUPP / name)
        sc = group_scores(name)
        d = with_corrected_scores(stored, sc)
        have = d.open_hit_rate.notna()
        check_true(f"{name}: every run within +-1.0 m has a state", bool(have[d.hfl_error_m.abs() <= 1.0].all()))
        info(f"{name}: runs without a corrected radar score", int((~have).sum()))
        eight[name] = d
    check("lhs8 runs without corrected radar score", 153, eight["lhs8"].open_hit_rate.isna().sum(), 0)
    check("lhs8b runs without corrected radar score", 155, eight["lhs8b"].open_hit_rate.isna().sum(), 0)
    exp8 = {("lhs8", 1.0): (20, 150, 210, 0.73, 47, 0.62, 0.07), ("lhs8", 0.5): (3, 150, 172, 0.80, None, None, None),
            ("lhs8b", 1.0): (21, 152, 215, 0.71, 45, 0.53, 0.20), ("lhs8b", 0.5): (9, 152, 196, 0.75, None, None, None)}
    for (name, tol), (n_e, lo_e, hi_e, c_e, nh_e, ch_e, rho_e) in exp8.items():
        d = eight[name]; k = accepted(d, tol); c = C_of(d, k)
        check(f"{name} +-{tol} accepted", n_e, c["n"], 0); check(f"{name} +-{tol} C", c_e, c["C"], 0.005)
        check(f"{name} +-{tol} min (km2)", lo_e, c["low"], 0.5); check(f"{name} +-{tol} max (km2)", hi_e, c["high"], 0.5)
        if nh_e is not None:
            h = (d.hfl_error_m.abs() <= tol).values
            check(f"{name} HFL-only passes", nh_e, h.sum(), 0)
            check(f"{name} HFL-only C", ch_e, C_of(d, h)["C"], 0.005)
            check(f"{name} rho depth-roughness, accepted", rho_e, spearmanr(d.channel_depth[k], d.n_floodplain[k])[0], 0.005)
            (b05, b95), (bn05, bn95) = bootstrap(d.fillna({"open_hit_rate": -1.0}), tol)
            exp_b = {"lhs8": (0.70, 0.79, 13, 27), "lhs8b": (0.69, 0.79, None, None)}[name]
            check(f"{name} bootstrap C 5%", exp_b[0], b05, 0.005); check(f"{name} bootstrap C 95%", exp_b[1], b95, 0.005)
            if exp_b[2] is not None:
                check(f"{name} bootstrap n 5%", exp_b[2], bn05, 0.5); check(f"{name} bootstrap n 95%", exp_b[3], bn95, 0.5)
    d = eight["lhs8"]; k = accepted(d, 1.0)
    for f, lo_e, hi_e in [("n_channel", 0.76, 1.30), ("q_idamalayar", 0.80, 1.19), ("tide_offset", -0.14, 0.46)]:
        check(f"lhs8 accepted {f} min", lo_e, d[f][k].min(), 0.005); check(f"lhs8 accepted {f} max", hi_e, d[f][k].max(), 0.005)
    pool = pd.concat([eight["lhs8"], eight["lhs8b"]], ignore_index=True); kp = accepted(pool, 1.0); cp = C_of(pool, kp)
    check("pooled eight-factor accepted", 41, cp["n"], 0); check("pooled eight-factor C", 0.73, cp["C"], 0.005)
    check("pooled eight-factor min (km2)", 150, cp["low"], 0.5); check("pooled eight-factor max (km2)", 215, cp["high"], 0.5)
    c8 = C_of(eight["lhs8"], accepted(eight["lhs8"], 1.0))
    alt = contraction(A, eight["lhs8"].loc[accepted(eight["lhs8"], 1.0), AREA].to_numpy())
    check("eight-factor accepted range against the five-factor prior, C", 0.73, alt["C"], 0.005)
    p8 = eight["lhs8"][AREA]
    check("eight-factor prior median (km2)", 222, p8.median(), 0.5)
    check("eight-factor prior 5% (km2)", 117, p8.quantile(0.05), 0.5); check("eight-factor prior 95% (km2)", 332, p8.quantile(0.95), 0.5)
    check("share of 200-run subsets with C >= 0.725 (%)", 4, 100 * (cs >= c8["C"]).mean(), 0.5)
    c8b = C_of(eight["lhs8b"], accepted(eight["lhs8b"], 1.0))["C"]
    check("share of 200-run subsets with C >= seed-46 C (%)", 6, 100 * (cs >= c8b).mean(), 0.5)
    check_true("eight-factor radar-only C is not computable (runs outside +-1.0 m have no state)",
               bool(eight["lhs8"].open_hit_rate.isna().any() and eight["lhs8b"].open_hit_rate.isna().any()))
    for tol in (1.5, 2.0):
        miss = sum(int(((e.hfl_error_m.abs() <= tol) & e.open_hit_rate.isna()).sum()) for e in eight.values())
        info(f"eight-factor runs passing the HFL at +-{tol} m without a corrected score (so C not computable)", miss)
    ops = ["hfl_err_nearest", "hfl_err_r300_median", "hfl_err_r750_median", "hfl_err_r1000_median", "hfl_err_r500_max"]
    nn, lo_o, hi_o = [], [], []
    for col in ops:
        k = accepted(d, 1.0, err=col); nn.append(k.sum()); lo_o.append(d[AREA][k].min()); hi_o.append(d[AREA][k].max())
    check("lhs8 operators: accepted min", 19, min(nn), 0); check("lhs8 operators: accepted max", 20, max(nn), 0)
    check("lhs8 operators: range min (km2)", 150, min(lo_o), 0.5); check("lhs8 operators: range max (km2)", 210, max(hi_o), 0.5)
    pas = d.hfl_err_r500_min.abs() <= 1.0
    check("lhs8 lowest-cell operator: passing runs without a state", 6, int((pas & d.open_hit_rate.isna()).sum()), 0)

    print("\nFixed factors around the 23 strict runs")
    fx = load(SUPP / "fixed").merge(pd.read_csv(SUPP / "fixed_parents.csv"), on="run")
    fx = with_corrected_scores(fx, group_scores("fixed"))
    check_true("fixed: every run within +-1.0 m has a state", bool(fx.open_hit_rate.notna()[fx.hfl_error_m.abs() <= 1.0].all()))
    k10, k05f = accepted(fx, 1.0), accepted(fx, 0.5)
    check("fixed accepted +-1.0 m", 57, k10.sum(), 0); check("fixed accepted +-0.5 m", 30, k05f.sum(), 0)
    check("fixed parents with an accepted variant", 21, fx[k10].parent.nunique(), 0)
    check("fixed accepted min (km2)", 146, fx[AREA][k10].min(), 0.5); check("fixed accepted max (km2)", 214, fx[AREA][k10].max(), 0.5)
    check("fixed rejected", 35, (~k10).sum(), 0)
    check("fixed rejected by radar", 25, ((fx.hfl_error_m.abs() <= 1.0) & ~k10).sum(), 0)
    check("fixed with HFL error above +1.0 m", 10, (fx.hfl_error_m > 1.0).sum(), 0)
    check("fixed mean channel roughness, rejected", 0.89, fx.n_channel[~k10].mean(), 0.005)
    check("fixed mean channel roughness, accepted", 1.07, fx.n_channel[k10].mean(), 0.005)
    check("fixed rho depth-roughness, accepted", 0.67, spearmanr(fx.channel_depth[k10], fx.n_floodplain[k10])[0], 0.005)
    info("fixed HFL-error change and linear fit are not radar dependent (checked by v1.0 check_results)", "-")

    print("\nForcing variants (Table 12, Table 13, Section 4.8)")
    samp = pd.read_csv(SUPP / "samples_accepted.csv"); cen = pd.read_csv(SUPP / "samples_central.csv")

    def parents(s):
        idx = []
        for _, row in s.iterrows():
            j = np.where(np.isclose(m.channel_depth, row.channel_depth) & np.isclose(m.n_floodplain, row.n_floodplain))[0]
            assert len(j) == 1
            idx.append(int(j[0]))
        return m.iloc[idx].reset_index(drop=True)

    exp12 = {("strict", "mult"): (22, 37, 13, 0.76, 2, 1, 16, 17), ("strict", "lagshort"): (21, 37, 17, 0.29, 2, 0, 15, 17),
             ("strict", "laglong"): (22, 62, 12, 0.64, 1, 1, 16, 16), ("central", "mult"): (12, 53, 8, 0.51, None, 0, 0, None),
             ("central", "lagshort"): (13, 41, 11, 0.51, None, 0, 0, None), ("central", "laglong"): (13, 41, 8, 0.51, None, 1, 0, None)}
    folders = {("strict", "mult"): SUPP / "multiplicative", ("strict", "lagshort"): SUPP / "lag_short",
               ("strict", "laglong"): SUPP / "lag_long", ("central", "mult"): CENTRAL / "multiplicative",
               ("central", "lagshort"): CENTRAL / "lag_short", ("central", "laglong"): CENTRAL / "lag_long"}
    base = {"strict": parents(samp), "central": parents(cen)}
    check("strict 23: accepted against the corrected map (both tolerances)", 17, accepted(base["strict"], 0.5).sum(), 0)
    check("strict 23: accepted at +-1.0 m", 17, accepted(base["strict"], 1.0).sum(), 0)
    check("central 14: accepted at +-1.0 m", 13, accepted(base["central"], 1.0).sum(), 0)
    check("central 14: accepted at +-0.5 m", 0, accepted(base["central"], 0.5).sum(), 0)
    tested = set(zip(base["strict"].design, base["strict"].run)) | set(zip(base["central"].design, base["central"].run))
    acc_ids = set(zip(m.design[keep10], m.run[keep10]))
    check("of the 69 accepted, among the 37 tested runs", 30, len(acc_ids & tested), 0)
    check("of the 69 accepted, not rerun under forcing variants", 39, len(acc_ids - tested), 0)
    for (grp, var), (a1, maxA, e1, maxe, out05, out10, n05, n10) in exp12.items():
        v = load(folders[(grp, var)]).sort_values("run").reset_index(drop=True)
        v = with_corrected_scores(v, group_scores(f"{var}_{grp}")).sort_values("run").reset_index(drop=True)
        b = base[grp]
        dA = (v[AREA] - b[AREA]).abs(); de = (v.hfl_error_m - b.hfl_error_m).abs()
        lab = f"{grp} {var}"
        check(f"{lab} |dA| < 1 km2", a1, (dA < 1).sum(), 0); check(f"{lab} max |dA|", maxA, dA.max(), 0.5)
        check(f"{lab} |de| < 0.01 m", e1, (de < 0.01).sum(), 0); check(f"{lab} max |de|", maxe, de.max(), 0.005)
        for tol, out_e, n_e in ((0.5, out05, n05), (1.0, out10, n10)):
            a0, a1_ = accepted(b, tol), accepted(v, tol)
            if out_e is not None:
                check(f"{lab} drop-outs +-{tol}", out_e, (a0 & ~a1_).sum(), 0)
                check_true(f"{lab} drop-outs +-{tol} all on the HFL", bool(((a0 & ~a1_) <= (v.hfl_error_m.abs() > tol).values).all()))
            if n_e is not None:
                check(f"{lab} accepted after change +-{tol}", n_e, a1_.sum(), 0)
            info(f"{lab} accepted range +-{tol} (km2)", (round(v[AREA][a1_].min(), 1), round(v[AREA][a1_].max(), 1)) if a1_.any() else None)
        if (grp, var) == ("strict", "mult"):
            g = accepted(v, 1.0) & ~accepted(b, 1.0)
            check("strict mult: runs gaining acceptance", 1, g.sum(), 0)
            check("strict mult: hit rate of the gaining run before", 0.749, b.open_hit_rate[g].iloc[0], 0.0005)
            check("strict mult: hit rate of the gaining run after", 0.750, v.open_hit_rate[g].iloc[0], 0.0005)
        if grp == "central":
            k1 = accepted(v, 1.0)
            info(f"central {var} accepted range +-1.0", (round(v[AREA][k1].min(), 1), round(v[AREA][k1].max(), 1)))
    bc = base["central"]
    j = bc.open_hit_rate < 0.75
    check("central run failing radar before change: hit rate", 0.747, bc.open_hit_rate[j].iloc[0], 0.0005)
    top = A[keep10]; top = np.sort(top)[::-1]
    check_true("four other accepted runs between 236 and 243 km2", bool(((top[1:] >= 235.5) & (top[1:] <= 243.5)).sum() == 4),
               str(np.round(top[:6], 1)))
    vm = load(CENTRAL / "multiplicative").sort_values("run").reset_index(drop=True)
    i_top = int(np.argmax(bc[AREA]))
    A2 = A.copy(); idx_top = int(np.where(np.isclose(A, bc[AREA].iloc[i_top]) & keep10)[0][0]); A2[idx_top] = vm[AREA].iloc[i_top]
    k2_ = keep10.copy()
    check("C with the multiplicative result substituted for the largest run", 0.61, contraction(A2, A2[k2_])["C"], 0.005)

    print("\nDrainage diagnosis against the corrected map")
    for name, exp in [("drainage_nominal", dict(cd_fa=115.9, fa=159.4, tree=0.83, csi_all=0.111, csi_open=0.402)),
                      ("drainage_filled", dict(csi_all=0.115, csi_open=0.398))]:
        z = np.load(R / "revision2" / name / "drainage_state.npz")
        fr = sample(_MAP_DATA, _MAP_TF, z["x"], z["y"]); land = z["land"] & (fr >= 0); obs = fr > 0.5
        wet = z["wet"] & land
        fa = wet & ~obs
        sc = score(z["x"], z["y"], z["landcover"], z["land"], z["wet"])
        check(f"{name} all-land CSI", exp["csi_all"], sc["all"]["csi"], 0.0005)
        check(f"{name} open-land CSI", exp["csi_open"], sc["open"]["csi"], 0.0005)
        if "fa" in exp:
            cd = fa & (z["classes"] == "closed depression")
            check(f"{name} false alarms (km2)", exp["fa"], fa.sum() * 0.04, 0.05)
            check(f"{name} closed-depression false alarms (km2)", exp["cd_fa"], cd.sum() * 0.04, 0.05)
            check(f"{name} closed-depression false alarms under trees", exp["tree"], (z["landcover"][cd] == TREE).mean(), 0.005)

    print("\nEvening pass (Section 4.5)")
    E = V11 / "evening"
    z = np.load(R / "revision" / "structural" / "nominal_radar_state.npz")

    def em(name):
        with rasterio.open(E / f"radar_{name}_corrected.tif") as src:
            return sample(src.read(1), src.transform, z["x"], z["y"])

    def overlap(a_n, b_n, dom, cut=0.5):
        a, b = em(a_n), em(b_n); base_ = z["land"].copy()
        if dom == "open":
            base_ &= np.isin(z["landcover"], OPEN)
        v = base_ & (a >= 0) & (b >= 0); wa, wb = (a > cut) & v, (b > cut) & v
        return (wa & wb).sum() / max((wa | wb).sum(), 1), (wa & wb).sum() / max(wa.sum(), 1)
    check("threshold-only overlap, all land", 0.77, overlap("morning_threshold_only", "evening_threshold_only", "all")[0], 0.005)
    check("threshold-only overlap, open land", 0.82, overlap("morning_threshold_only", "evening_threshold_only", "open")[0], 0.005)
    check("share of morning full-rule water also in the evening", 0.84,
          overlap("morning_published_rule", "evening_threshold_only", "all")[1], 0.005)
    ov = [overlap(f"morning_threshold_only_{s}dB", f"evening_threshold_only_{s}dB", dm)[0] for s in ("-1.5", "+1.5") for dm in ("all",)]
    ovo = [overlap(f"morning_threshold_only_{s}dB", f"evening_threshold_only_{s}dB", "open")[0] for s in ("-1.5", "+1.5")]
    check("+-1.5 dB overlap min", 0.76, min(ov), 0.005); check("+-1.5 dB overlap max", 0.78, max(ov), 0.005)
    check("+-1.5 dB open-land overlap min", 0.79, min(ovo), 0.005); check("+-1.5 dB open-land overlap max", 0.86, max(ovo), 0.005)
    info("full-rule overlap (both over either), all land", round(overlap("morning_published_rule", "evening_threshold_only", "all")[0], 3))
    old = pd.read_csv(E / "morning_evening_radar_comparison_uncorrected.csv")
    sel = old[(old.fraction_cut == 0.5) & (old.domain == "all_land") & old.scene_a.isin(["morning_published_rule", "morning_threshold_only"])
              & old.scene_b.str.startswith("evening")]
    check("before correction: overlap min (%)", 4, 100 * sel.overlap_both_over_either.min(), 0.5)
    check("before correction: overlap max (%)", 5, 100 * sel.overlap_both_over_either.max(), 0.5)

    print("\nReruns (Section 3.8)")
    rr = pd.read_csv(V11 / "rerun_reproduction.csv")
    check("reruns", 884, len(rr), 0); check("exact", 874, rr.exact.sum(), 0); check("differing", 10, (~rr.exact).sum(), 0)
    check("largest HFL-error difference (m) [proposed correction; DOCX says 0.21]", 0.20, rr.diff_hfl_error_m.abs().max(), 0.005)
    check("largest peak-land difference (km2)", 0.8, rr.diff_max_flooded_area_km2.abs().max(), 0.05)
    by = rr.groupby("group").size()
    check("main-ensemble reruns", 593, by.get("lhs200", 0) + by.get("lhs400", 0), 0)
    check("eight-factor reruns", 92, by.get("lhs8", 0) + by.get("lhs8b", 0), 0)
    check("fixed-factor reruns", 82, by.get("fixed", 0), 0)
    check("forcing reruns", 111, sum(by.get(g, 0) for g in by.index if g.startswith(("mult", "lag"))), 0)
    check("repeat reruns", 6, by.get("repeat_add", 0) + by.get("repeat_mult", 0), 0)
    rr_main = rr[rr.group.isin(["lhs200", "lhs400"]) & ~rr.exact]
    flips = 0
    for _, row in rr_main.iterrows():
        mm = m[(m.design == row.group) & (m.run == row.run)].iloc[0]
        for tol in (0.5, 1.0, 1.5, 2.0):
            a = abs(row.stored_hfl_error_m) <= tol; b = abs(row.rerun_hfl_error_m) <= tol
            flips += int(a != b)
    rr_e = rr[rr.group.isin(["lhs8", "lhs8b"]) & ~rr.exact]
    for _, row in rr_e.iterrows():
        for tol in (0.5, 1.0):
            flips += int((abs(row.stored_hfl_error_m) <= tol) != (abs(row.rerun_hfl_error_m) <= tol))
    check("HFL-test decisions changed by the 10 differing reruns", 0, flips, 0)

    unchanged()
    print(f"\n{n_checks} checks, {failures} failed")
    return failures


def unchanged():
    """Values that do not depend on the radar map, read from the stored v1.0 outputs."""
    print("\nValues not affected by the radar correction (stored v1.0 outputs)")
    mo = pd.read_csv(R / "v5_morris_r20" / "indices_morris_r20.csv").set_index(["metric", "factor"])
    t5 = {"channel_depth": (149, 57, 516, 271, 12.5, 13.3, 2.92, 0.63), "n_floodplain": (101, 96, 414, 379, 42.1, 52.0, 2.86, 2.46),
          "q_rating": (65, 73, 258, 298, 22.3, 38.4, 2.00, 1.31), "canopy_factor": (41, 38, 56, 118, 10.4, 26.6, 1.74, 1.98),
          "q_idukki": (38, 77, 149, 309, 18.9, 41.1, 1.00, 1.78), "n_channel": (31, 16, 119, 74, 10.4, 20.8, 1.00, 0.17),
          "q_idamalayar": (14, 38, 56, 154, 6.9, 11.2, 0.36, 0.60), "tide_offset": (6, 1, 5, 2, 1.9, 4.1, 0.00, 0.00)}
    mets = ["max_flooded_area_km2", "max_flood_volume_mcm", "hours_to_peak_area", "neeleswaram_peak_level_m"]
    for f, v in t5.items():
        for i, met in enumerate(mets):
            t = 0.5 if i < 2 else (0.05 if i == 2 else 0.005)
            check(f"Table 5 {f} {met} mu*", v[2 * i], mo.loc[(met, f), "mu_star"], t)
            check(f"Table 5 {f} {met} sigma", v[2 * i + 1], mo.loc[(met, f), "sigma"], t)
    a = mo.loc["max_flooded_area_km2", "mu_star"]
    check("top three share of summed mu* for area (%)", 71, 100 * a[["channel_depth", "n_floodplain", "q_rating"]].sum() / a.sum(), 0.5)
    check("three fixed factors share of summed mu* for area (%)", 11, 100 * a[["n_channel", "q_idamalayar", "tide_offset"]].sum() / a.sum(), 0.5)
    mr = pd.read_csv(R / "v5_morris_r20" / "results.csv")
    check("Morris runs", 180, len(mr), 0)
    check("Morris peak land min (km2)", 86, mr[AREA].min(), 0.5); check("Morris peak land max (km2)", 388, mr[AREA].max(), 0.5)
    check("Morris volume min (Mm3)", 281, mr.max_flood_volume_mcm.min(), 0.5); check("Morris volume max (Mm3)", 1385, mr.max_flood_volume_mcm.max(), 0.5)
    m4 = pd.read_csv(R / "v4_morris_r20" / "indices_morris_r20.csv").set_index(["metric", "factor"]).mu_star["max_flooded_area_km2"]
    for f, e in [("channel_depth", 150), ("n_floodplain", 97), ("q_rating", 63)]:
        check(f"v4 mu* {f}", e, m4[f], 0.5)
    m6 = pd.read_csv(SUPP / "morris6" / "indices_morris.csv").set_index(["metric", "factor"])
    for f, e in [("channel_depth", 164), ("n_floodplain", 103), ("canopy_factor", 76), ("q_rating", 56), ("n_channel", 38), ("tide_offset", 6)]:
        check(f"six-factor Morris mu* {f}", e, m6.loc[("max_flooded_area_km2", f), "mu_star"], 0.5)
    r6 = pd.read_csv(SUPP / "morris6" / "results.csv")
    check("six-factor Morris runs", 140, len(r6), 0)
    check("six-factor Morris peak land min", 86, r6[AREA].min(), 0.5); check("six-factor Morris peak land max", 399, r6[AREA].max(), 0.5)
    st = pd.read_csv(R / "revision" / "structural" / "structural_tests.csv").set_index("variant")
    t10 = {"nominal": (263, 180, 855, 2.44), "min_terrain": (353, 245, 934, -0.76), "full_canopy": (292, 201, 860, 1.31),
           "wide_channel": (210, 140, 643, 1.84), "grid_100m": (263, 179, 921, 3.32), "tide_x1.5": (265, 180, 856, 2.44),
           "tide_shift_6h": (264, 182, 854, 2.44)}
    for v, (pk, ps, vol, e) in t10.items():
        check(f"Table 10 {v} peak land", pk, st.loc[v, AREA], 0.5); check(f"Table 10 {v} land at pass", ps, st.loc[v, "model_flooded_km2"], 0.5)
        check(f"Table 10 {v} peak volume", vol, st.loc[v, "max_flood_volume_mcm"], 0.5); check(f"Table 10 {v} HFL error", e, st.loc[v, "hfl_error_m"], 0.005)
    dn = json.loads((R / "revision2" / "drainage_nominal" / "drainage_summary.json").read_text())
    check("nominal flooded land at pass (km2)", 179.8, dn["flooded_land_at_pass_km2"], 0.05)
    check("nominal land storage (Mcm)", 355, dn["land_storage_at_pass_mcm"], 0.5)
    check("already ponded 24 h earlier (km2)", 100.9, dn["ponded"]["already_ponded_24h_before_km2"], 0.05)
    check("gauge-cell peak spread (m)", 0.26, dn["gauge_peak_level_spread_m"], 0.005)
    check("gauge-cell peak min (m)", 14.62, min(dn["gauge_cell_peak_levels_m"]), 0.005)
    check("gauge-cell peak max (m)", 14.88, max(dn["gauge_cell_peak_levels_m"]), 0.005)
    check("sample SD of gauge-cell peaks (m)", 0.10, np.std(dn["gauge_cell_peak_levels_m"], ddof=1), 0.005)
    check("sigma_obs (m) [DOCX: approx. 0.50; exact 0.49]", 0.49, np.sqrt(0.46 ** 2 + 0.13 ** 2 + 0.10 ** 2), 0.005)
    check("depression depth median (m)", 2.9, dn["ponded"]["depression_depth_m_quantiles"]["q50"], 0.05)
    check("depression depth 10% (m)", 0.95, dn["ponded"]["depression_depth_m_quantiles"]["q10"], 0.005)
    check("depression depth 90% (m)", 7.3, dn["ponded"]["depression_depth_m_quantiles"]["q90"], 0.05)
    cl = pd.read_csv(R / "revision2" / "drainage_nominal" / "drainage_classes.csv").set_index("class")
    check("closed depression (km2)", 128.1, cl.loc["closed depression", "flooded_km2"], 0.05)
    check("closed depression storage (Mcm)", 250, cl.loc["closed depression", "stored_mcm"], 0.5)
    check("closed depression mean depth (m)", 1.95, cl.loc["closed depression", "mean_depth_m"], 0.005)
    check("closed depression next to channel (km2)", 34.5, cl.loc["closed depression", "next_to_channel_km2"], 0.05)
    check("draining (km2)", 38.0, cl.loc["draining", "flooded_km2"], 0.05)
    check("connected, slow (km2)", 13.6, cl.loc["connected, slow", "flooded_km2"], 0.05)
    th = pd.read_csv(R / "revision2" / "drainage_nominal" / "drainage_thresholds.csv")
    allc = th[th.cells == "all"]; exc = th[th.cells != "all"]
    check("closed class over thresholds, min (km2)", 124.6, allc.closed_km2.min(), 0.05)
    check("closed class over thresholds, max (km2)", 128.1, allc.closed_km2.max(), 0.05)
    check("closed class without channel-side cells, min", 90.1, exc.closed_km2.min(), 0.05)
    check("closed class without channel-side cells, max", 93.6, exc.closed_km2.max(), 0.05)
    d01 = allc[np.isclose(allc.depression_depth_min_m, 0.1)]
    check("draining over level-fall thresholds (depression depth 0.1 m), min", 21.6, d01.draining_km2.min(), 0.05)
    check("draining over level-fall thresholds (depression depth 0.1 m), max", 40.4, d01.draining_km2.max(), 0.05)
    td = pd.read_csv(R / "revision2" / "drainage_nominal" / "terrain_depressions.csv").set_index("terrain")
    for t, (ar, sto) in {"canopy factor 0.3": (352, 1255), "canopy factor 0.6": (363, 1196), "canopy factor 0.9": (375, 1096),
                         "canopy factor 1.0": (379, 1067), "bare-earth proxy": (298, 981)}.items():
        check(f"Table 11b {t} area", ar, td.loc[t, "depression_area_km2_over_0.1m"], 0.5)
        check(f"Table 11b {t} storage", sto, td.loc[t, "depression_storage_mcm"], 0.5)
    fl = pd.read_csv(R / "revision2" / "drainage_filled" / "drainage_classes.csv").set_index("class")
    df_ = json.loads((R / "revision2" / "drainage_filled" / "drainage_summary.json").read_text())
    check("filled: flooded land (km2)", 153.0, df_["flooded_land_at_pass_km2"], 0.05)
    check("filled: storage (Mcm)", 187, df_["land_storage_at_pass_mcm"], 0.5)
    check("filled: closed (km2)", 48.2, fl.loc["closed depression", "flooded_km2"], 0.05)
    check("filled: closed next to channel (km2)", 31.9, fl.loc["closed depression", "next_to_channel_km2"], 0.05)
    check("filled: connected slow (km2)", 57.6, fl.loc["connected, slow", "flooded_km2"], 0.05)
    check("filled: connected slow mean depth (m)", 0.59, fl.loc["connected, slow", "mean_depth_m"], 0.005)
    check("filled: draining (km2)", 47.2, fl.loc["draining", "flooded_km2"], 0.05)
    rs = json.loads((R / "revision2" / "drainage_filled" / "run_summary.json").read_text())
    check("filled: HFL error (m)", 2.66, rs["hfl_error_m"], 0.005)
    check("filled: area reduction (%)", 15, 100 * (1 - df_["flooded_land_at_pass_km2"] / dn["flooded_land_at_pass_km2"]), 0.5)
    md = V11 / "rescoring" / "drainage_median_lhs400_70"
    ms = json.loads((md / "drainage_summary.json").read_text()); mc = pd.read_csv(md / "drainage_classes.csv").set_index("class")
    check("median run flooded at pass (km2)", 121.4, ms["flooded_land_at_pass_km2"], 0.05)
    check("median run closed (km2)", 97.1, mc.loc["closed depression", "flooded_km2"], 0.05)
    check("median run closed next to channel (km2)", 30.2, mc.loc["closed depression", "next_to_channel_km2"], 0.05)
    check("median run closed storage (Mcm)", 192, mc.loc["closed depression", "stored_mcm"], 0.5)
    check("median run land storage (Mcm)", 246, ms["land_storage_at_pass_mcm"], 0.5)


if __name__ == "__main__":
    sys.exit(main())
