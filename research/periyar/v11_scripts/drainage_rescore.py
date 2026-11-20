"""Rescore the stored drainage-diagnosis states (nominal and depressions-filled runs; original drainage.py output,
unchanged classes) against the CORRECTED radar map: false alarms by class and CSI. Usage: python -I drainage_rescore.py"""
import sys, json
from pathlib import Path
import numpy as np
sys.path.insert(0, '/home/user/hydrolib-core')
from research.periyar.validate import OPEN_CLASSES, _sample
from research.periyar.radar_check import confusion
R = Path('/home/user/hydrolib-core/research/periyar/results/revision2'); MAP = Path('/home/user/ext/sar_flood_20180821_corrected.tif')
EF = Path('/mnt/project-files/periyar/EMA_FINAL_REVISION'); out = {}
for name in ('drainage_nominal', 'drainage_filled', 'drainage_accepted'):
    p = R / name / 'drainage_state.npz'
    if not p.exists(): continue
    s = np.load(p); f = _sample(MAP, s['x'], s['y'], -1.0)
    land = s['land'] & (f >= 0); wet = s['wet'] & land; radar_c = f > 0.5; cell = 0.04
    o = {'published_radar_cells_on_land': int((s['radar'] & s['land']).sum()), 'corrected_radar_cells_on_land': int((radar_c & land).sum())}
    for cls in ('closed depression', 'draining', 'connected, slow'):
        mm = wet & (s['classes'] == cls)
        o[cls] = dict(flooded_km2=round(mm.sum() * cell, 1), false_alarm_km2_published=round((mm & ~s['radar']).sum() * cell, 1),
                      false_alarm_km2_corrected=round((mm & ~radar_c).sum() * cell, 1),
                      tree_share_of_corrected_fa=round(float(np.isin(s['landcover'][mm & ~radar_c], [10]).mean()), 3))
    o['false_alarm_km2_published'] = round((wet & ~s['radar']).sum() * cell, 1); o['false_alarm_km2_corrected'] = round((wet & ~radar_c).sum() * cell, 1)
    for dom, msk in (('all', land), ('open', land & np.isin(s['landcover'], OPEN_CLASSES))):
        o[f'{dom}_corrected'] = {k: round(v, 3) if isinstance(v, float) else v for k, v in confusion(s['wet'][msk], radar_c[msk]).items()}
    out[name] = o
json.dump(out, open(EF / 'rescoring/drainage_rescored.json', 'w'), indent=1); print(json.dumps(out, indent=1))
