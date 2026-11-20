"""Before/after radar comparison on the model's land cells (all land and open land), fraction cut-off 0.5:
published vs corrected morning map, and the morning/evening hold-out with both maps corrected."""
import sys
from pathlib import Path
import numpy as np, pandas as pd, rasterio
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, '/home/user/hydrolib-core')
from research.periyar.validate import OPEN_CLASSES, _sample
OUT = Path(sys.argv[1]); RC = Path('/mnt/project-files/periyar/EMA_FINAL_REVISION/radar_corrected'); AN = Path('/mnt/project-files/periyar/EMA_FINAL_REVISION/analysis')
M = {'morning_published': '/home/user/ext/sar_flood_20180821_v3.tif', 'morning_corrected': RC / 'sar_flood_20180821_corrected.tif',
     'evening_threshold_only_uncorrected': AN / 'radar_evening_threshold_only.tif', 'evening_threshold_only_corrected': RC / 'sar_evening_threshold_only_corrected.tif',
     'morning_threshold_only_corrected': RC / 'sar_morning_threshold_only_corrected.tif', 'morning_threshold_only_uncorrected': AN / 'radar_morning_threshold_only.tif'}
s = np.load('/home/user/ext/structural/nominal/radar_state.npz'); km2 = 0.04
F = {k: _sample(Path(v), s['x'], s['y'], -1.0) for k, v in M.items()}
rows = []
for a, b in (('morning_published', 'morning_corrected'), ('morning_published', 'evening_threshold_only_uncorrected'),
             ('morning_corrected', 'evening_threshold_only_corrected'), ('morning_threshold_only_uncorrected', 'evening_threshold_only_uncorrected'),
             ('morning_threshold_only_corrected', 'evening_threshold_only_corrected')):
    for dom in ('all_land', 'open_land'):
        base = s['land'] & (np.isin(s['landcover'], OPEN_CLASSES) if dom == 'open_land' else True)
        v = base & (F[a] >= 0) & (F[b] >= 0); wa, wb = (F[a] > 0.5) & v, (F[b] > 0.5) & v
        rows.append(dict(map_a=a, map_b=b, domain=dom, cells_compared=int(v.sum()), water_a_km2=wa.sum() * km2, water_b_km2=wb.sum() * km2,
                         both_km2=(wa & wb).sum() * km2, a_only_km2=(wa & ~wb).sum() * km2, b_only_km2=(~wa & wb).sum() * km2,
                         overlap_both_over_either=(wa & wb).sum() / max((wa | wb).sum(), 1), share_of_a_in_b=(wa & wb).sum() / max(wa.sum(), 1)))
df = pd.DataFrame(rows); df.to_csv(OUT / 'radar_before_after_comparison.csv', index=False); print(df.round(3).to_string())
fig, ax = plt.subplots(1, 3, figsize=(13, 4.8), constrained_layout=True)
for k, (a, t) in enumerate((('morning_published', '(a) Published morning map'), ('morning_corrected', '(b) Corrected morning map'), ('evening_threshold_only_corrected', '(c) Evening map (threshold only), corrected'))):
    land = s['land']; ax[k].scatter(s['x'][land] / 1e3, s['y'][land] / 1e3, s=1, c='0.88', lw=0)
    w = land & (F[a] > 0.5); ax[k].scatter(s['x'][w] / 1e3, s['y'][w] / 1e3, s=2.5, c='#1f5fa8', lw=0)
    ax[k].set_title(t, fontsize=10); ax[k].set_aspect('equal'); ax[k].set_xlabel('Easting (km, UTM 43N)')
ax[0].set_ylabel('Northing (km)')
fig.savefig(OUT / 'radar_before_after.png', dpi=600)
