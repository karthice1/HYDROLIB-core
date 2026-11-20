"""Evening Sentinel-1 hold-out check (21 Aug 2018 13:06 UTC, ascending, relative orbit 173) repeated on the
CO-REGISTERED scenes. Identical to EMA_FINAL_REVISION/analysis/evening_holdout.py (frozen sigma0, Otsu,
aggregation, coverage and channel rules; the same pairs and the same +-1.5 dB sensitivity) except that each
scene is first shifted by the fixed co-registration offset in radar_corrected/georef_qc.json.
The model is not rerun; model states at the MORNING pass are scored against each map as a diagnostic only
(the simulation ends 21 Aug 07:00 IST = 01:30 UTC, before the evening pass).
Usage: python -I evening_holdout_corrected.py OUT_DIR"""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd, rasterio
sys.path.insert(0, '/home/user/hydrolib-core')
from research.periyar.sar_flood_map import sigma0_db, otsu, FINE
from research.periyar.validate import OPEN_CLASSES, _sample
from research.periyar.radar_check import confusion
from research.periyar.uncertainty import accepted
OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
EF = Path('/mnt/project-files/periyar/EMA_FINAL_REVISION')
DATA = Path('/home/user/ext/data/s1'); MODEL = Path('/home/user/ext/periyar_v5')
QC = json.load(open(EF / 'radar_corrected/georef_qc.json'))
CORRECTED = Path('/home/user/ext/sar_flood_20180821_corrected.tif')
STATES = EF / 'baseline_2026-11-17/model_states'; RR = Path('/home/user/ext/rerun_radar')
with rasterio.open(MODEL / 'dem_model.tif') as src: profile = src.profile
with rasterio.open(MODEL / 'channel.tif') as src: channel = src.read(1) > 0.5
fac = int(round(profile['transform'].a / FINE)); h, w = profile['height'], profile['width']
def shift(a, dy, dx):
    out = np.full_like(a, np.nan); H, W = a.shape
    ys, xs = slice(max(0, -dy), min(H, H - dy)), slice(max(0, -dx), min(W, W - dx))
    out[ys, xs] = a[ys.start + dy:ys.stop + dy, xs.start + dx:xs.stop + dx]; return out
coarse = lambda a: a[:h * fac, :w * fac].reshape(h, fac, w, fac).mean(axis=(1, 3))
def fraction_map(fl, valid):
    return np.where((coarse(valid.astype('float32')) > 0.5) & ~channel, coarse(fl.astype('float32')), -1.0).astype('float32')
sc = {k: shift(sigma0_db(DATA / k, profile), *QC[k]['shift_px']) for k in ('20180821', '20180728', '20180821_asc')}
morning, pre, evening = sc['20180821'], sc['20180728'], sc['20180821_asc']
t_m, t_e = otsu(morning), otsu(evening)
maps = {'morning_published_rule': fraction_map((morning < t_m) & (pre - morning > 3.0), np.isfinite(morning) & np.isfinite(pre)),
        'morning_threshold_only': fraction_map(morning < t_m, np.isfinite(morning)),
        'evening_threshold_only': fraction_map(evening < t_e, np.isfinite(evening))}
for d in (-1.5, 1.5):
    maps[f'morning_threshold_only_{d:+.1f}dB'] = fraction_map(morning < t_m + d, np.isfinite(morning))
    maps[f'evening_threshold_only_{d:+.1f}dB'] = fraction_map(evening < t_e + d, np.isfinite(evening))
maps['evening_threshold_only_at_morning_threshold'] = fraction_map(evening < t_m, np.isfinite(evening))
for k, m in maps.items():
    with rasterio.open(OUT / f'radar_{k}_corrected.tif', 'w', **{**profile, 'dtype': 'float32', 'nodata': -1.0, 'count': 1}) as d: d.write(m, 1)
with rasterio.open(CORRECTED) as src: corr = src.read(1)
same = bool(np.array_equal(maps['morning_published_rule'], corr))
s = np.load(STATES / 'structural__nominal/radar_state.npz'); x, y, landm, lc = s['x'], s['y'], s['land'], s['landcover']
cell = profile['transform'].a ** 2 / 1e6
smp = lambda n: _sample(OUT / f'radar_{n}_corrected.tif', x, y, -1.0)
def compare(a_n, b_n, cut, dom):
    a, b = smp(a_n), smp(b_n); base = landm.copy()
    if dom == 'open_land': base &= np.isin(lc, OPEN_CLASSES)
    bv = base & (a >= 0) & (b >= 0); wa, wb = (a > cut) & bv, (b > cut) & bv; both, either = wa & wb, wa | wb
    return dict(scene_a=a_n, scene_b=b_n, domain=dom, fraction_cut=cut, cells_compared=int(bv.sum()),
                water_a_km2=wa.sum() * cell, water_b_km2=wb.sum() * cell, water_both_km2=both.sum() * cell,
                water_to_nonwater_km2=(wa & ~wb).sum() * cell, nonwater_to_water_km2=(~wa & wb).sum() * cell,
                overlap_both_over_either=both.sum() / max(either.sum(), 1),
                share_of_a_also_in_b=both.sum() / max(wa.sum(), 1), share_of_b_also_in_a=both.sum() / max(wb.sum(), 1))
pairs = [('morning_published_rule', 'evening_threshold_only'), ('morning_threshold_only', 'evening_threshold_only'),
         ('morning_threshold_only', 'evening_threshold_only_at_morning_threshold'),
         ('morning_threshold_only_-1.5dB', 'evening_threshold_only_-1.5dB'), ('morning_threshold_only_+1.5dB', 'evening_threshold_only_+1.5dB'),
         ('morning_published_rule', 'morning_threshold_only')]
df = pd.DataFrame([compare(a, b, c, d) for a, b in pairs for c in (0.5, 0.3) for d in ('all_land', 'open_land')])
df.to_csv(OUT / 'morning_evening_radar_comparison_corrected.csv', index=False)
# diagnostic: morning-pass model states vs each map
runs = pd.read_csv(EF / 'rescoring/runs_corrected.csv'); k = accepted(runs, hfl_tol=1.0)
st = {'nominal': STATES / 'structural__nominal/radar_state.npz', 'median_accepted_lhs400_70': RR / 'lhs400/run_0070/radar_state.npz'}
for i in np.where(k)[0]:
    d, r = runs.design[i], int(runs.run[i]); p = RR / f'{d}/run_{r:04d}/radar_state.npz'
    st[f'acc_{d}_{r}'] = p if p.exists() else STATES / f'gsa_v5_{d}__time_series__run_{r:04d}/radar_state.npz'
diag = []
for run, p in st.items():
    z = np.load(p)
    for n in ('morning_published_rule', 'morning_threshold_only', 'evening_threshold_only'):
        fr = _sample(OUT / f'radar_{n}_corrected.tif', z['x'], z['y'], -1.0); L = z['land'] & (fr >= 0)
        for dom, msk in (('all_land', L), ('open_land', L & np.isin(z['landcover'], OPEN_CLASSES))):
            diag.append(dict(run=run, radar_map=n, domain=dom, **confusion(z['model_wet'][msk], fr[msk] > 0.5)))
dg = pd.DataFrame(diag); dg.to_csv(OUT / 'evening_model_diagnostic_corrected.csv', index=False)
acc = dg[dg.run.str.startswith('acc_')]
summ = acc.groupby(['radar_map', 'domain'])[['hit_rate', 'false_alarm_ratio', 'csi']].agg(['min', 'max']).round(3)
summ.to_csv(OUT / 'evening_model_diagnostic_accepted69_summary.csv')
log = dict(otsu_db=dict(morning=round(float(t_m), 2), evening=round(float(t_e), 2)), morning_rule_equals_corrected_map=same,
           shifts_px={k: QC[k]['shift_px'] for k in sc}, cells_over_0_5={k: int((m > 0.5).sum()) for k, m in maps.items()},
           accepted_runs_scored=int(k.sum()))
json.dump(log, open(OUT / 'evening_holdout_corrected_log.json', 'w'), indent=1)
print(json.dumps(log, indent=1)); print(df[df.fraction_cut == 0.5].round(3).to_string()); print(summ.to_string())
