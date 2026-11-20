"""Compare every rerun (884) with the stored v1.0 result of the same parameter set: peak flooded land, HFL error,
model flooded land at the pass and the uncorrected-map open-land hit rate. Writes rerun_reproduction.csv.
Usage: python -I rerun_reproduction_v11.py OUT_CSV"""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
R = Path('/home/user/hydrolib-core/research/periyar/results'); S = R / 'revision2/supplementary'; C = R / 'revision2/central'
RR = Path('/home/user/ext/rerun_radar')
STORED = {'lhs200': R / 'v5_lhs200', 'lhs400': R / 'v5_lhs400', 'lhs8': S / 'lhs8', 'lhs8b': S / 'lhs8b', 'fixed': S / 'fixed',
          'mult_strict': S / 'multiplicative', 'lagshort_strict': S / 'lag_short', 'laglong_strict': S / 'lag_long',
          'mult_central': C / 'multiplicative', 'lagshort_central': C / 'lag_short', 'laglong_central': C / 'lag_long',
          'repeat_add': S / 'repeat_additive', 'repeat_mult': S / 'repeat_multiplicative'}
PAIRS = [('max_flooded_area_km2', 'max_flooded_area_km2'), ('hfl_error_m', 'published_hfl_error_m'),
         ('model_flooded_km2', 'published_model_flooded_km2'), ('open_hit_rate', 'published_open_hit_rate')]
rows = []
for g, folder in STORED.items():
    st = pd.read_csv(folder / 'results.csv').set_index('run')
    for f in sorted((RR / g).glob('run_*/result.json')):
        run = int(f.parent.name[4:]); v = json.load(open(f)); o = st.loc[run]
        row = dict(group=g, run=run, stored_table=str(folder.relative_to(R)))
        for a, b in PAIRS:
            row[f'stored_{a}'] = o[a]; row[f'rerun_{a}'] = v[b]; row[f'diff_{a}'] = v[b] - o[a]
        rows.append(row)
t = pd.DataFrame(rows)
d = t[[c for c in t if c.startswith('diff_')]].abs()
t['exact'] = (d.max(axis=1) < 1e-9)
t.to_csv(sys.argv[1], index=False)
print(len(t), 'reruns;', int(t.exact.sum()), 'exact;', int((~t.exact).sum()), 'differ')
print(t[~t.exact][['group', 'run', 'diff_max_flooded_area_km2', 'diff_hfl_error_m', 'diff_model_flooded_km2', 'diff_open_hit_rate']].round(3).to_string())
