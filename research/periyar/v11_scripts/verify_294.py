"""Verify the 294 runs of batch 4 (|HFL error| > 2 m) from their actual outputs, not from markers.
For each run: listed id, stored |HFL error|, result.json status, corrected-score keys present, reproduction
of the stored results (published-map scores vs stored CSV), radar_state.npz present and loadable, independent
rescoring of the state against the corrected map equal to the stored corrected score, and use in
runs_corrected.csv. Usage: python -I verify_294.py OUT_CSV"""
import sys, json, hashlib
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, '/home/user/hydrolib-core')
from research.periyar.validate import _sample, OPEN_CLASSES
from research.periyar.radar_check import confusion
EA = Path('/home/user/ext/ema_analysis'); RR = Path('/home/user/ext/rerun_radar')
MAP = Path('/home/user/ext/sar_flood_20180821_corrected.tif')
REPO = Path('/home/user/hydrolib-core/research/periyar/results')
RC = pd.read_csv('/mnt/project-files/periyar/EMA_FINAL_REVISION/rescoring/runs_corrected.csv')
rows = []
for des in ('lhs200', 'lhs400'):
    ids = [int(x) for x in (EA / f'rest_{des}.txt').read_text().split()]
    stored = pd.read_csv(REPO / f'v5_{des}/runs_with_behaviour.csv').set_index('run')
    for r in ids:
        d = RR / des / f'run_{r:04d}'; res = d / 'result.json'; st = d / 'radar_state.npz'
        row = dict(design=des, run=r, stored_abs_hfl_error=abs(stored.loc[r, 'hfl_error_m']), result_json=res.exists(), state=st.exists())
        if res.exists():
            j = json.loads(res.read_text()); row['status'] = j.get('status')
            row['has_corrected_scores'] = all(f'corrected_{k}' in j for k in ('open_hit_rate', 'all_critical_success_index', 'open_false_alarm_ratio', 'model_flooded_km2'))
            row['repro_hfl_diff'] = abs(j['published_hfl_error_m'] - stored.loc[r, 'hfl_error_m'])
            row['repro_area_diff'] = abs(j['max_flooded_area_km2'] - stored.loc[r, 'max_flooded_area_km2'])
            row['repro_pub_hit_diff'] = abs(j['published_open_hit_rate'] - stored.loc[r, 'open_hit_rate'])
            if st.exists():
                s = np.load(st); f = _sample(MAP, s['x'], s['y'], -1.0); L = s['land'] & (f >= 0); O = L & np.isin(s['landcover'], OPEN_CLASSES)
                c_open = confusion(s['model_wet'][O], (f > 0.5)[O]); c_all = confusion(s['model_wet'][L], (f > 0.5)[L])
                row['state_vs_json_hit_diff'] = abs(c_open['hit_rate'] - j['corrected_open_hit_rate'])
                row['state_vs_json_csi_diff'] = abs(c_all['csi'] - j['corrected_all_critical_success_index'])
            rc = RC[(RC.design == des) & (RC.run == r)].iloc[0]
            row['used_in_runs_corrected'] = rc.source == 'rerun' and abs(rc.open_hit_rate - j['corrected_open_hit_rate']) < 1e-12
            row['radar_pass_corrected'] = bool(j['corrected_open_hit_rate'] >= 0.75 and j['corrected_model_flooded_km2'] >= 22)
        rows.append(row)
df = pd.DataFrame(rows); df.to_csv(sys.argv[1], index=False)
print('runs', len(df), 'by design', df.design.value_counts().to_dict())
print('all |HFL|>2:', bool((df.stored_abs_hfl_error > 2).all()), ' min', df.stored_abs_hfl_error.min().round(3))
for c in ('result_json', 'state', 'has_corrected_scores', 'used_in_runs_corrected'): print(c, int(df[c].sum()))
print('status', df.status.value_counts().to_dict())
for c in ('repro_hfl_diff', 'repro_area_diff', 'repro_pub_hit_diff', 'state_vs_json_hit_diff', 'state_vs_json_csi_diff'): print(c, 'max', df[c].max(), 'n>1e-9', int((df[c] > 1e-9).sum()))
print('pass corrected radar test', int(df.radar_pass_corrected.sum()))
print('map sha256', hashlib.sha256(MAP.read_bytes()).hexdigest())
print('duplicate ids vs earlier batches:', len(set(zip(df.design, df.run)) & set((d, int(x)) for d in ('lhs200', 'lhs400') for x in (EA / f'need_{d}.txt').read_text().split())))
