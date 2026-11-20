"""Corrected-map acceptance for the supplementary ensembles (eight-factor designs; fixed-factor variants),
unchanged rule (uncertainty.accepted / narrowing). Only runs with |HFL error| <= 1.0 m were rerun, so results are
complete for tolerances <= 1.0 m. Usage: python -I recompute_supplementary.py OUT_DIR"""
import sys, json, glob
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, '/home/user/hydrolib-core')
from research.periyar import uncertainty as U
OUT = Path(sys.argv[1]); B = Path('/mnt/project-files/periyar/EMA_FINAL_REVISION/baseline_2026-11-17/results_v1.0/revision2')
K = ['all_critical_success_index', 'open_critical_success_index', 'open_hit_rate', 'open_false_alarm_ratio']
def load(old, rerun_dir):
    new = old.copy(); new[K] = np.nan; diff = []
    for f in glob.glob(f'{rerun_dir}/run_*/result.json'):
        v = json.load(open(f)); run = int(Path(f).parent.name[4:]); i = new.index[new.run == run][0]
        for k in K: new.loc[i, k] = v[f'corrected_{k}']
        diff.append(max(abs(v[f'published_{k}'] - old.loc[i, k]) for k in ['hfl_error_m', 'model_flooded_km2'] + K))
    return new, diff
rows, notes = [], {}
for name in ('lhs8', 'lhs8b'):
    old = pd.read_csv(B / f'{name}_tol1.0/runs.csv'); new, diff = load(old, f'/home/user/ext/rerun_radar/{name}')
    notes[name] = dict(rerun=len(diff), needed=int((old.hfl_error_m.abs() <= 1.0).sum()), max_repro_diff=max(diff), n_not_exact=int(sum(d > 1e-9 for d in diff)))
    prior = old[U.OUTPUT].to_numpy()
    for tol in (0.5, 1.0):
        a0, a1 = U.accepted(old, tol), U.accepted(new, tol, sar_km2=22.96)
        n0, n1 = U.narrowing(prior, old.loc[a0, U.OUTPUT].to_numpy()), U.narrowing(prior, new.loc[a1, U.OUTPUT].to_numpy())
        rows.append(dict(ensemble=name, hfl_tolerance_m=tol, n_published=n0['n'], C_published=n0['narrowing'], n_corrected=n1['n'], C_corrected=n1['narrowing'],
                         range_corrected=f"{n1['low']:.1f}-{n1['high']:.1f}", accepted_all_csi=f"{new.loc[a1, K[0]].min():.3f}-{new.loc[a1, K[0]].max():.3f}"))
    new.to_csv(OUT / f'{name}_runs_corrected.csv', index=False)
pd.DataFrame(rows).to_csv(OUT / 'supplementary_contraction_before_after.csv', index=False)
fx = sorted(glob.glob('/home/user/ext/rerun_radar/fixed/run_*/result.json'))
if fx:
    old = pd.read_csv('/home/user/ext/rev2/fixed/results.csv'); new, diff = load(old, '/home/user/ext/rerun_radar/fixed')
    notes['fixed'] = dict(rerun=len(diff), needed=int((old.hfl_error_m.abs() <= 1.0).sum()), max_repro_diff=max(diff), n_not_exact=int(sum(d > 1e-9 for d in diff)))
    for tol in (0.5, 1.0):
        h = old.hfl_error_m.abs() <= tol
        notes['fixed'][f'tol{tol}'] = dict(published_accepted=int(U.accepted(old, tol).sum()), corrected_accepted=int(U.accepted(new, tol, sar_km2=22.96).sum()),
                                          hfl_pass=int(h.sum()), corrected_fail_radar=int((h & (new.open_hit_rate < 0.75)).sum()),
                                          corrected_range=[float(new.loc[U.accepted(new, tol, sar_km2=22.96), U.OUTPUT].min()), float(new.loc[U.accepted(new, tol, sar_km2=22.96), U.OUTPUT].max())])
    new.to_csv(OUT / 'fixed_runs_corrected.csv', index=False)
json.dump(notes, open(OUT / 'supplementary_notes.json', 'w'), indent=1)
print(pd.DataFrame(rows).round(3).to_string()); print(json.dumps(notes, indent=1))
