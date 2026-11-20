"""Classification-sensitivity variants of the CORRECTED radar map, built exactly as the published variants
(research/periyar/results/revision/sar_variants: Otsu +-1.5 dB, change threshold 2 and 4 dB) but from the
co-registered scenes (shifts fixed in georef_qc.json; no new search). Frozen sigma0 and aggregation code.
Usage: python -I sar_variants_corrected.py OUT_DIR"""
import sys, json
from pathlib import Path
import numpy as np, rasterio
sys.path.insert(0, '/home/user/hydrolib-core')
import research.periyar.sar_flood_map as frozen
OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
S1 = Path('/home/user/ext/data/s1'); MODEL = Path('/home/user/ext/periyar_v5')
QC = json.load(open('/mnt/project-files/periyar/EMA_FINAL_REVISION/radar_corrected/georef_qc.json'))
with rasterio.open(MODEL / 'dem_model.tif') as src: profile = src.profile
with rasterio.open(MODEL / 'channel.tif') as src: channel = src.read(1) > 0.5
h, w = profile['height'], profile['width']; fac = int(round(profile['transform'].a / frozen.FINE))
def shift(a, dy, dx):
    out = np.full_like(a, np.nan); H, W = a.shape
    ys, xs = slice(max(0, -dy), min(H, H - dy)), slice(max(0, -dx), min(W, W - dx))
    out[ys, xs] = a[ys.start + dy:ys.stop + dy, xs.start + dx:xs.stop + dx]; return out
coarse = lambda a: a[:h * fac, :w * fac].reshape(h, fac, w, fac).mean(axis=(1, 3))
def fraction_map(fl, valid):
    return np.where((coarse(valid.astype('float32')) > 0.5) & ~channel, coarse(fl.astype('float32')), -1.0).astype('float32')
post = shift(frozen.sigma0_db(S1 / '20180821', profile), *QC['20180821']['shift_px'])
pre = shift(frozen.sigma0_db(S1 / '20180728', profile), *QC['20180728']['shift_px'])
t = frozen.otsu(post); valid = np.isfinite(post) & np.isfinite(pre)
variants = {'base': (t, 3.0), f'thr_m{-(t + 1.5):.2f}': (t + 1.5, 3.0), f'thr_m{-(t - 1.5):.2f}': (t - 1.5, 3.0),
            'drop2': (t, 2.0), 'drop4': (t, 4.0)}
info = {}
for k, (thr, drop) in variants.items():
    m = fraction_map((post < thr) & (pre - post > drop), valid)
    with rasterio.open(OUT / f'{k}.tif', 'w', **{**profile, 'dtype': 'float32', 'nodata': -1.0, 'count': 1}) as d: d.write(m, 1)
    info[k] = dict(threshold_db=round(float(thr), 2), min_drop_db=drop, cells_over_0_5=int((m > 0.5).sum()),
                   km2_over_0_5=round(float((m > 0.5).sum()) * 0.04, 2))
json.dump(info, open(OUT / 'variants.json', 'w'), indent=1); print(json.dumps(info, indent=1))
