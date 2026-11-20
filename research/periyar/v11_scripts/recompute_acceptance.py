"""Recompute every radar-dependent acceptance statistic of the pooled 600-run ensemble with the corrected radar map,
using the UNCHANGED rule and functions of research/periyar/uncertainty.py (accepted, narrowing, intervals).
Inputs: baseline runs.csv (read-only copy), rerun results (rerun_radar/lhs200, lhs400) and rescored stored states.
Runs with |HFL error| > 2.0 m are never accepted by the joint rule at any reported tolerance (<= 2.0 m), so their
open-land hit rate is not needed for the joint statistics; it is left missing and every statistic that needs it
(radar-only acceptance, likelihood weighting, all-600 score ranges) is flagged as not recomputed.
Usage: python -I recompute_acceptance.py OUT_DIR"""
import sys, json, glob
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, '/home/user/hydrolib-core')
from research.periyar import uncertainty as U
OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
BASE = Path('/mnt/project-files/periyar/EMA_FINAL_REVISION/baseline_2026-11-17')
old = pd.read_csv(BASE / 'results_v1.0/revision2/uncertainty_600_tol1.0/runs.csv')
old['design'] = ['lhs200'] * 200 + ['lhs400'] * 400
SCORES = ['model_flooded_km2', 'all_critical_success_index', 'open_critical_success_index', 'open_hit_rate', 'open_false_alarm_ratio']
new = old.copy(); new[SCORES[1:]] = np.nan; new['source'] = 'not needed (|HFL error| > 2 m)'; new['repro_max_abs_diff'] = np.nan
sar_km2 = set()
for d in ('lhs200', 'lhs400'):
    for f in glob.glob(f'/home/user/ext/rerun_radar/{d}/run_*/result.json'):
        v = json.load(open(f)); run = int(Path(f).parent.name[4:]); i = new.index[(new.design == d) & (new.run == run)][0]
        if v.get('status') != 'ok': new.loc[i, 'source'] = 'rerun failed'; continue
        for k in SCORES[1:]: new.loc[i, k] = v[f'corrected_{k}']
        new.loc[i, 'repro_max_abs_diff'] = max(abs(v[f'published_{k}'] - old.loc[i, k]) for k in ['hfl_error_m', 'model_flooded_km2'] + SCORES[1:])
        new.loc[i, 'source'] = 'rerun'; sar_km2.add(round(v['corrected_sar_flooded_km2'], 2))
st = pd.read_csv('/mnt/project-files/periyar/EMA_FINAL_REVISION/radar_corrected/rescored_stored_states.csv')
st = st[(st.radar_map == 'corrected') & st.state.str.startswith('gsa_v5_lhs200')]
for run in st.state.str[-4:].astype(int).unique():
    i = new.index[(new.design == 'lhs200') & (new.run == run)][0]
    if new.loc[i, 'source'] == 'rerun': continue
    s = st[st.state.str.endswith(f'{run:04d}')].set_index('domain')
    new.loc[i, ['all_critical_success_index', 'open_critical_success_index', 'open_hit_rate', 'open_false_alarm_ratio']] = \
        [s.loc['all_land', 'csi'], s.loc['open_land', 'csi'], s.loc['open_land', 'hit_rate'], s.loc['open_land', 'false_alarm_ratio']]
    new.loc[i, 'source'] = 'stored state rescored'
need = new.hfl_error_m.abs() <= 2.0
missing = need & new.open_hit_rate.isna()
SAR_NEW = max(sar_km2) if sar_km2 else np.nan
prior = old[U.OUTPUT].to_numpy()
res = {'n_runs': len(new), 'n_needed_for_joint_rule': int(need.sum()), 'n_needed_missing': int(missing.sum()),
       'sar_flooded_land_km2': {'published_constant': U.SAR_KM2, 'corrected_map': SAR_NEW,
                                'min_model_flooded_km2': float(new.model_flooded_km2.min()),
                                'note': 'area test is non-binding for every run with either value'},
       'reruns_reproduce_stored_results_max_abs_diff': float(new.repro_max_abs_diff.max())}
rows = []
for tol in (0.25, 0.5, 0.75, 1.0, 1.5, 2.0):
    for lab, kw in (('hfl_only', dict(use_radar=False)), ('both', {})):
        a_old = U.accepted(old, tol, **kw); a_new = U.accepted(new, tol, sar_km2=SAR_NEW, **kw) if lab == 'both' else U.accepted(new, tol, **kw)
        n_old, n_new = U.narrowing(prior, old.loc[a_old, U.OUTPUT].to_numpy()), U.narrowing(prior, new.loc[a_new, U.OUTPUT].to_numpy())
        rows.append(dict(hfl_tolerance_m=tol, test=lab, n_published=n_old['n'], C_published=n_old['narrowing'],
                         n_corrected=n_new['n'], C_corrected=n_new['narrowing'], range_corrected=(n_new['low'], n_new['high']),
                         runs_dropped=';'.join(f'{d}:{r}' for d, r in new.loc[a_old & ~a_new, ['design', 'run']].values),
                         runs_added=';'.join(f'{d}:{r}' for d, r in new.loc[~a_old & a_new, ['design', 'run']].values),
                         complete=bool(lab == 'hfl_only' or not (missing & (new.hfl_error_m.abs() <= tol)).any())))
pd.DataFrame(rows).to_csv(OUT / 'contraction_before_after.csv', index=False)
# threshold grid (joint rule, complete for tol <= 2.0)
grid = [dict(hfl_tolerance_m=t, hit_rate=h, **U.narrowing(prior, new.loc[U.accepted(new, t, h, sar_km2=SAR_NEW), U.OUTPUT].to_numpy()))
        for t in U.HFL_TOLERANCES for h in U.HIT_RATES]
pd.DataFrame(grid).to_csv(OUT / 'threshold_grid_corrected.csv', index=False)
# bootstrap of C at +-1.0 m, same procedure and seed as uncertainty.py
rng = np.random.default_rng(1); boot = []
for _ in range(5000):
    s = new.iloc[rng.integers(0, len(new), len(new))]
    boot.append(U.narrowing(s[U.OUTPUT].to_numpy(), s.loc[U.accepted(s, 1.0, sar_km2=SAR_NEW), U.OUTPUT].to_numpy())['narrowing'])
res['bootstrap_C_tol1.0_q05_q50_q95'] = [float(x) for x in np.nanquantile(boot, [0.05, 0.5, 0.95])]
# spatial scores of accepted runs
for tol in (0.5, 1.0):
    a = U.accepted(new, tol, sar_km2=SAR_NEW); s = new.loc[a]
    res[f'accepted_tol{tol}'] = dict(n=int(a.sum()), all_csi=[float(s.all_critical_success_index.min()), float(s.all_critical_success_index.max())],
                                     open_far=[float(s.open_false_alarm_ratio.min()), float(s.open_false_alarm_ratio.max())],
                                     open_hit=[float(s.open_hit_rate.min()), float(s.open_hit_rate.max())])
a05 = new.loc[U.accepted(new, 0.5, sar_km2=SAR_NEW)].sort_values(U.OUTPUT)
res['median_accepted_run_tol0.5'] = a05.iloc[(len(a05) - 1) // 2][['design', 'run', U.OUTPUT, 'all_critical_success_index', 'open_false_alarm_ratio']].to_dict() if len(a05) else None
sub = new[need & ~missing]
res['over_hfl_pass_runs_within_2m'] = dict(n=len(sub), max_all_csi=float(sub.all_critical_success_index.max()),
                                           open_far=[float(sub.open_false_alarm_ratio.min()), float(sub.open_false_alarm_ratio.max())])
res['not_recomputed'] = ['radar-only acceptance and C (needs the 294 runs with |HFL error| > 2 m)',
                         'likelihood weighting (uses every run)', 'score ranges over all 600 runs']
if new.open_hit_rate.notna().all():
    res['not_recomputed'] = []
    width = np.diff(np.quantile(prior, [0.05, 0.95]))[0]
    ro = []
    for tol in (0.5, 1.0):
        for lab, kw in (('hfl_only', dict(use_radar=False)), ('radar_only', dict(use_hfl=False)), ('both', {})):
            a0 = U.accepted(old, tol, **kw); a1 = U.accepted(new, tol, sar_km2=SAR_NEW, **kw)
            ro.append(dict(hfl_tolerance_m=tol, test=lab, n_published=int(a0.sum()), C_published=U.narrowing(prior, old.loc[a0, U.OUTPUT].to_numpy())['narrowing'],
                           n_corrected=int(a1.sum()), C_corrected=U.narrowing(prior, new.loc[a1, U.OUTPUT].to_numpy())['narrowing']))
    pd.DataFrame(ro).to_csv(OUT / 'contributions_corrected.csv', index=False)
    gl = []
    for power in (0, 1, 4, 8):
        for sigma in (0.25, 0.5, 1.0):
            for lab, r_ in (('published', old), ('corrected', new)):
                w = np.exp(-0.5 * (r_.hfl_error_m / sigma) ** 2) * r_.open_hit_rate ** power
                w = np.where(r_.model_flooded_km2 >= (U.SAR_KM2 if lab == 'published' else SAR_NEW), w, 0.0)
                lo, hi = U.weighted_quantiles(prior, w, [0.05, 0.95])
                gl.append(dict(map=lab, hit_rate_power=power, sigma_m=sigma, q05=lo, q95=hi, effective_runs=w.sum() ** 2 / (w ** 2).sum(), narrowing_5_95=1 - (hi - lo) / width))
    pd.DataFrame(gl).to_csv(OUT / 'likelihood_weights_before_after.csv', index=False)
    res['all_600'] = {m: dict(published=[float(old[m].min()), float(old[m].max())], corrected=[float(new[m].min()), float(new[m].max())])
                      for m in ('all_critical_success_index', 'open_false_alarm_ratio', 'open_hit_rate', 'open_critical_success_index')}
    res['all_600_spearman_with_peak_area'] = {m: float(new[[m, U.OUTPUT]].corr(method='spearman').iloc[0, 1]) for m in ('open_hit_rate', 'open_false_alarm_ratio', 'all_critical_success_index')}
    res['repro_not_exact_runs'] = new.loc[new.repro_max_abs_diff > 1e-9, ['design', 'run', 'repro_max_abs_diff']].values.tolist()
    print(pd.DataFrame(ro).round(3).to_string()); print(pd.DataFrame(gl).round(3).to_string())
new.to_csv(OUT / 'runs_corrected.csv', index=False)
json.dump(res, open(OUT / 'summary_corrected.json', 'w'), indent=1, default=str)
print(json.dumps(res, indent=1, default=str)); print(pd.DataFrame(rows)[['hfl_tolerance_m', 'test', 'n_published', 'C_published', 'n_corrected', 'C_corrected', 'complete']].round(3).to_string())
