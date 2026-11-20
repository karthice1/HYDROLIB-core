"""Rerun selected runs of the pooled v5 LHS ensemble (unchanged samples, model, forcing and code) only to
recover the per-cell model state at the radar pass, which the original ensemble deleted after scoring.
For each run: prepare and run D-Flow FM exactly as research/periyar/sensitivity.run_one, then
  (a) score against the PUBLISHED map (must reproduce the stored results: determinism check),
  (b) score against the CORRECTED map,
  (c) save radar_state.npz (x, y, land, model_wet, landcover) so that any map can be rescored later,
and delete the map file. Resumable: a run with result.json is skipped.
Usage: python -I rerun_for_radar.py DESIGN_DIR OUT_DIR CORRECTED_MAP RUN_ID [RUN_ID ...] [--parallel 4]"""
import sys, json, os, subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np, pandas as pd
sys.path.insert(0, '/home/user/hydrolib-core')
from research.periyar import sensitivity as S, validate
from research.periyar.structural import radar_state

args = sys.argv[1:]; par = 4; release = 'additive'
if '--parallel' in args: i = args.index('--parallel'); par = int(args[i + 1]); del args[i:i + 2]
if '--release-model' in args: i = args.index('--release-model'); release = args[i + 1]; del args[i:i + 2]
if '--lags' in args:
    i = args.index('--lags'); S.IDUKKI_LAG_HOURS, S.IDAMALAYAR_LAG_HOURS = float(args[i + 1]), float(args[i + 2]); del args[i:i + 3]
design, out, corrected = Path(args[0]), Path(args[1]), Path(args[2]); runs = [int(a) for a in args[3:]]
BASE = Path('/home/user/ext/periyar_v5'); PUBLISHED = Path('/home/user/ext/sar_flood_20180821_v3.tif')
DFM = '/home/user/ext/run_dflowfm.sh'
samples = pd.read_csv(design / 'samples.csv', index_col=0)
KEYS = ('neeleswaram_peak_level_m', 'hfl_error_m', 'model_flooded_km2', 'all_critical_success_index',
        'open_critical_success_index', 'open_hit_rate', 'open_false_alarm_ratio')

def one(run):
    rd = out / f'run_{run:04d}'; res = rd / 'result.json'
    if res.exists(): return run, json.loads(res.read_text())
    S.prepare(BASE, rd, samples.loc[run], release)
    with (rd / 'run.log').open('w') as f:
        rc = subprocess.run([DFM, '--autostartstop', 'periyar.mdu'], cwd=rd, stdout=f, stderr=subprocess.STDOUT,
                            env={**os.environ, 'OMP_NUM_THREADS': '1'}).returncode
    mf = rd / 'DFM_OUTPUT_periyar' / 'periyar_map.nc'
    if rc != 0 or not mf.exists(): return run, {'status': 'failed'}
    v = {'status': 'ok', **S.metrics(mf)}
    pub = validate.scores(rd, mf, PUBLISHED); cor = validate.scores(rd, mf, corrected)
    v.update({f'published_{k}': pub[k] for k in KEYS}); v.update({f'corrected_{k}': cor[k] for k in KEYS})
    v.update({f'corrected_{k}': cor[k] for k in cor if k.startswith(('all_', 'open_', 'sar_')) and isinstance(cor[k], (int, float))})
    v.update({f'published_{k}': pub[k] for k in pub if k.startswith(('all_', 'open_', 'sar_')) and isinstance(pub[k], (int, float))})
    st = radar_state(rd, mf, corrected); st.pop('radar_fraction')
    np.savez_compressed(rd / 'radar_state.npz', **st)
    mf.unlink()
    for p in rd.glob('*.nc'): p.unlink()        # network copy; the state file is what is kept
    res.write_text(json.dumps(v, indent=1, default=float))
    print(run, 'done', flush=True); return run, v

out.mkdir(parents=True, exist_ok=True)
with ThreadPoolExecutor(max_workers=par) as pool:
    rows = dict(pool.map(one, runs))
pd.DataFrame.from_dict(rows, orient='index').rename_axis('run').sort_index().to_csv(out / 'results.csv')
print('all done', len(rows))
