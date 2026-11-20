"""Forcing tests (release form, travel times, repeats) rescored on the corrected radar map with the unchanged rule.
Stored HFL errors and areas are used; only the open-land hit rate comes from the reruns (each rerun is also
checked against the stored values). Usage: python -I recompute_forcing.py OUT_DIR"""
import sys, json, glob
from pathlib import Path
import numpy as np, pandas as pd
OUT = Path(sys.argv[1]); RV = Path('/mnt/project-files/periyar/EMA_FINAL_REVISION/baseline_2026-11-17/results_v1.0/revision2')
RR = Path('/home/user/ext/rerun_radar')
def acc(r, tol): return (r.hfl_error_m.abs() <= tol) & (r.model_flooded_km2 >= 22.96) & (r.open_hit_rate >= 0.75)
def rescored(stored_csv, rerun_dir):
    s = pd.read_csv(stored_csv).sort_values('run').reset_index(drop=True); diff = []
    s['open_hit_rate_published'] = s.open_hit_rate
    for i, row in s.iterrows():
        v = json.load(open(rerun_dir / f'run_{int(row.run):04d}' / 'result.json'))
        diff.append(max(abs(v['published_hfl_error_m'] - row.hfl_error_m), abs(v['published_open_hit_rate'] - row.open_hit_rate)))
        s.loc[i, 'open_hit_rate'] = v['corrected_open_hit_rate']; s.loc[i, 'all_critical_success_index'] = v['corrected_all_critical_success_index']
    return s, max(diff), int(sum(d > 1e-9 for d in diff))
rows = []
runs600 = pd.read_csv('/mnt/project-files/periyar/EMA_FINAL_REVISION/rescoring/runs_corrected.csv')
for group, sub, names in (('strict (old 23 accepted at +-0.5 m)', 'supplementary', [('multiplicative', 'mult_strict'), ('lag_short', 'lagshort_strict'), ('lag_long', 'laglong_strict')]),
                          ('central (14 of old 50 accepted at +-1.0 m only)', 'central', [('multiplicative', 'mult_central'), ('lag_short', 'lagshort_central'), ('lag_long', 'laglong_central')])):
    for name, d in names:
        r, md, ne = rescored(RV / sub / name / 'results.csv', RR / d)
        row = dict(group=group, variant=name, runs=len(r), repro_max_diff=md, repro_not_exact=ne)
        for tol in (0.5, 1.0, 1.5):
            a = acc(r, tol); a_old = (r.hfl_error_m.abs() <= tol) & (r.model_flooded_km2 >= 22) & (r.open_hit_rate_published >= 0.75)
            row[f'accepted_{tol}_published'] = int(a_old.sum()); row[f'accepted_{tol}_corrected'] = int(a.sum())
            row[f'area_{tol}_corrected'] = f'{r.max_flooded_area_km2[a].min():.1f}-{r.max_flooded_area_km2[a].max():.1f}' if a.any() else ''
        keep = acc(r, 1.0)
        row.update(dropouts_1_0=int((~keep).sum()), fail_hfl=int((~keep & (r.hfl_error_m.abs() > 1.0)).sum()),
                   fail_radar=int((~keep & (r.open_hit_rate < 0.75)).sum()))
        rows.append(row)
rep = {}
for name, d in (('repeat_additive', 'repeat_add'), ('repeat_multiplicative', 'repeat_mult')):
    files = sorted((RR / d).glob('run_*/result.json'))
    if files:
        s = pd.read_csv(RV / 'supplementary' / name / 'results.csv').set_index('run')
        rep[name] = {int(f.parent.name[4:]): dict(stored_hfl=float(s.loc[int(f.parent.name[4:]), 'hfl_error_m']), rerun_hfl=json.load(open(f))['published_hfl_error_m'],
                                                   corrected_open_hit=json.load(open(f))['corrected_open_hit_rate']) for f in files}
df = pd.DataFrame(rows); df.to_csv(OUT / 'forcing_tests_before_after.csv', index=False)
json.dump(rep, open(OUT / 'repeat_runs.json', 'w'), indent=1)
print(df.to_string()); print(json.dumps(rep, indent=1))
