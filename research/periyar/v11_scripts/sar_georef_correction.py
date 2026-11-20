"""Corrected georeferencing of the Sentinel-1 scenes and rebuilt radar flood map (version 1.1, 17 Nov 2026; 1.1 fixes the sign of the reported displacement only).

Problem: research/periyar/sar_flood_map.sigma0_db places each GRD image from the GeoTIFF ground-control
points (GCPs) without terrain correction. A GCP is exact only at its annotated (WGS84 ellipsoidal) height.
For the descending scenes the GCPs around the domain sit at 558-755 m, whereas the lowland is near sea level,
so the image is displaced by about dh/tan(incidence) along the look direction.

Correction (fixed before any rescoring and independent of the model and of the flood):
  1. sigma0 is computed and geocoded exactly as in the frozen code (sigma0_db, unchanged).
  2. Each scene is co-registered to an independent geographic reference, JRC Global Surface Water
     occurrence (Pekel et al. 2016), by the integer 20 m shift (dy, dx) that maximises the correlation between
     sigma0 (dB) and -occurrence over the model domain (coarse search +-1.6 km at 100 m, then 20 m steps).
     The same rule is applied to all three scenes.
  3. Uniformity of the shift is checked in four quadrants of the domain.
  4. Physical check: the shift direction is compared with the radar look direction (platform heading + 90 deg)
     and the implied height error dh = |shift| * tan(incidence) with the annotated GCP heights.
  5. Held-out verification, references not used in step 2: ESA WorldCover 2021 built-up share and the
     model channel mask (ring contrast).
  6. The flood map is rebuilt with the frozen rules and functions (Otsu on the post scene, drop > 3 dB from
     the same-track pre scene, 200 m aggregation, coverage > 0.5, channel cells removed).
No threshold, rule or model output is changed. Usage: python -I sar_georef_correction.py OUT_DIR
"""
import sys, json, hashlib
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np, rasterio
from rasterio.transform import from_origin
from rasterio.warp import reproject, Resampling
from scipy.ndimage import binary_dilation
sys.path.insert(0, '/home/user/hydrolib-core')
import research.periyar.sar_flood_map as frozen   # unchanged frozen code

OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
DATA = Path('/home/user/ext/data'); S1 = DATA / 's1'; MODEL = Path('/home/user/ext/periyar_v5')
with rasterio.open(MODEL / 'dem_model.tif') as src: profile = src.profile
with rasterio.open(MODEL / 'channel.tif') as src: channel = src.read(1) > 0.5
left, bottom, right, top = frozen._bounds(profile, native=True); F = frozen.FINE
shape = (int(round((top - bottom) / F)), int(round((right - left) / F))); dst_t = from_origin(left, top, F, F)
h, w = profile['height'], profile['width']; fac = int(round(profile['transform'].a / F))

def to_grid(path, fn):
    with rasterio.open(path) as src:
        a = fn(src.read(1)).astype('float32'); out = np.full(shape, np.nan, 'float32')
        reproject(a, out, src_transform=src.transform, src_crs=src.crs, dst_transform=dst_t, dst_crs=profile['crs'],
                  dst_nodata=np.nan, resampling=Resampling.average)
    return out
occ = np.fmax(*(to_grid(DATA / f'water/occurrence_70E_{t}v1_4_2021.tif', lambda a: np.where(a > 100, 0, a)) for t in ('10N', '20N')))
built = to_grid(DATA / 'landcover/ESA_WorldCover_10m_2021_v200_N09E075_Map.tif', lambda a: a == 50)
CH = np.kron(channel, np.ones((fac, fac), bool))
RG = np.kron(binary_dilation(channel, iterations=4) & ~binary_dilation(channel, iterations=2), np.ones((fac, fac), bool))

def shift(a, dy, dx):
    """corrected(i, j) = a(i + dy, j + dx); NaN where undefined."""
    out = np.full_like(a, np.nan); H, W = a.shape
    ys, xs = slice(max(0, -dy), min(H, H - dy)), slice(max(0, -dx), min(W, W - dx))
    out[ys, xs] = a[ys.start + dy:ys.stop + dy, xs.start + dx:xs.stop + dx]; return out

def corr(a, ref, sign, mask=None):
    ok = np.isfinite(a) & np.isfinite(ref) & (True if mask is None else mask)
    return sign * np.corrcoef(a[ok], ref[ok])[0, 1]

def search(sar, ref, sign, mask=None, R=80):
    def score(dy, dx):
        return corr(shift(sar, dy, dx), ref, sign, mask)
    c = max((score(dy, dx), dy, dx) for dy in range(-R, R + 1, 5) for dx in range(-R, R + 1, 5))
    return max((score(dy, dx), dy, dx) for dy in range(c[1] - 5, c[1] + 6) for dx in range(c[2] - 5, c[2] + 6))

def ring_contrast(a):
    return float(np.nanmedian(a[RG]) - np.nanmedian(a[CH]))

def annotation(name):
    r = ET.parse(S1 / name / 'annotation' / 'iw-vv.xml').getroot()
    g = np.array([[float(p.find(k).text) for k in ('latitude', 'longitude', 'height', 'incidenceAngle')]
                  for p in r.iter('geolocationGridPoint')])
    sel = (g[:, 1] > 76.1) & (g[:, 1] < 76.6) & (g[:, 0] > 9.8) & (g[:, 0] < 10.4)
    return dict(pass_=r.find('.//pass').text, heading=float(r.find('.//platformHeading').text),
                gcp_height_m=[float(g[sel, 2].min()), float(g[sel, 2].max())], incidence_deg=float(g[sel, 3].mean()))

qc, corrected = {}, {}
Y, X = np.mgrid[0:shape[0], 0:shape[1]]
quads = {'NW': (Y < shape[0] // 2) & (X < shape[1] // 2), 'NE': (Y < shape[0] // 2) & (X >= shape[1] // 2),
         'SW': (Y >= shape[0] // 2) & (X < shape[1] // 2), 'SE': (Y >= shape[0] // 2) & (X >= shape[1] // 2)}
for name in ('20180821', '20180728', '20180821_asc'):
    raw = frozen.sigma0_db(S1 / name, profile)
    r0, dy, dx = search(raw, occ, -1)
    cor = shift(raw, dy, dx); corrected[name] = cor
    ann = annotation(name)
    north, east = -dy * F, dx * F                   # displacement of the raw image relative to the truth (m)
    disp_bearing = float(np.degrees(np.arctan2(east, north)) % 360)
    look = (ann['heading'] + 90) % 360
    mag = float(np.hypot(north, east))
    q = {k: search(raw, occ, -1, m, R=60)[1:] for k, m in quads.items()}
    ver = {}
    for lab, a in (('raw', raw), ('corrected', cor)):
        b = search(a, built, 1, R=40)
        ver[lab] = dict(jrc_r_at_zero=round(corr(a, occ, -1), 3), builtup_r_at_zero=round(corr(a, built, 1), 3),
                        builtup_best_shift_m=[b[1] * F, b[2] * F], channel_ring_contrast_db=round(ring_contrast(a), 2))
    qc[name] = dict(annotation=ann, shift_px=[dy, dx], image_displacement_m=dict(north=north, east=east, magnitude=mag),
                    displacement_bearing_deg=round(disp_bearing, 1), look_direction_deg=round(look, 1),
                    implied_height_error_m=round(mag * np.tan(np.radians(ann['incidence_deg'])), 0),
                    jrc_r_best=round(r0, 3), quadrant_shifts_px=q, verification=ver)
    print(name, json.dumps(qc[name]), flush=True)

# rebuild the map with the frozen rules (identical to sar_flood_map.main / evening_holdout.py)
coarse = lambda a: a[:h * fac, :w * fac].reshape(h, fac, w, fac).mean(axis=(1, 3))
def fraction_map(fl, valid):
    return np.where((coarse(valid.astype('float32')) > 0.5) & ~channel, coarse(fl.astype('float32')), -1.0).astype('float32')
post, pre, eve = corrected['20180821'], corrected['20180728'], corrected['20180821_asc']
t_m, t_e = frozen.otsu(post), frozen.otsu(eve)
maps = {'sar_flood_20180821_corrected': fraction_map((post < t_m) & (pre - post > 3.0), np.isfinite(post) & np.isfinite(pre)),
        'sar_evening_threshold_only_corrected': fraction_map(eve < t_e, np.isfinite(eve)),
        'sar_morning_threshold_only_corrected': fraction_map(post < t_m, np.isfinite(post))}
for k, m in maps.items():
    with rasterio.open(OUT / f'{k}.tif', 'w', **{**profile, 'dtype': 'float32', 'nodata': -1.0, 'count': 1}) as d: d.write(m, 1)
qc['otsu_threshold_db'] = {'morning_corrected': round(float(t_m), 2), 'evening_corrected': round(float(t_e), 2)}
qc['maps'] = {k: dict(valid_cells=int((m >= 0).sum()), cells_over_0_5=int((m > 0.5).sum()),
                      sha256=hashlib.sha256(open(OUT / f'{k}.tif', 'rb').read()).hexdigest()) for k, m in maps.items()}
qc['script_sha256'] = hashlib.sha256(open(__file__, 'rb').read()).hexdigest()
json.dump(qc, open(OUT / 'georef_qc.json', 'w'), indent=1)
print(json.dumps({k: qc[k] for k in ('otsu_threshold_db', 'maps')}, indent=1))
