"""Rescore every stored per-cell model state (radar_state.npz) against the published and the corrected radar
maps with the unchanged scoring rule of research/periyar/validate.py (land cells, fraction > 0.5 = water,
channel and low-coverage cells excluded, open land = WorldCover 30/40/60/90)."""
import sys, glob
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, '/home/user/hydrolib-core')
from research.periyar.validate import OPEN_CLASSES, _sample
OUT = Path(sys.argv[1])
MAPS = {'published': '/home/user/ext/sar_flood_20180821_v3.tif', 'corrected': '/home/user/ext/sar_flood_20180821_corrected.tif'}
states = {p.split('/ext/')[1].replace('/radar_state.npz', ''): p for p in sorted(glob.glob('/home/user/ext/gsa_v5_lhs200/time_series/run_*/radar_state.npz') + glob.glob('/home/user/ext/rev2/maps/run_*/radar_state.npz') + glob.glob('/home/user/ext/structural/*/radar_state.npz'))}
rows = []
for name, p in states.items():
    st = np.load(p)
    for mname, mp in MAPS.items():
        fr = _sample(Path(mp), st['x'], st['y'], -1.0); landm = st['land'] & (fr >= 0)
        for dom, mask in (('all_land', landm), ('open_land', landm & np.isin(st['landcover'], OPEN_CLASSES))):
            o = fr[mask] > 0.5; m = st['model_wet'][mask]
            H, M, FA = int((m & o).sum()), int((~m & o).sum()), int((m & ~o).sum())
            rows.append(dict(state=name, radar_map=mname, domain=dom, radar_water_cells=int(o.sum()), hits=H, misses=M, false_alarms=FA,
                             hit_rate=H / max(H + M, 1), false_alarm_ratio=FA / max(H + FA, 1), csi=H / max(H + M + FA, 1)))
df = pd.DataFrame(rows); df.to_csv(OUT / 'rescored_stored_states.csv', index=False)
w = df.pivot_table(index=['state', 'domain'], columns='radar_map', values=['hit_rate', 'false_alarm_ratio', 'csi']).round(3)
print(w.to_string())
