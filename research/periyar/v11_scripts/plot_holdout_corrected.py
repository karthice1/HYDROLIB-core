"""Before/after figure for the evening hold-out check: (a,b) uncorrected morning rule vs evening threshold map,
(c) their overlap; (d,e,f) the same after co-registration. Usage: python -I plot_holdout_corrected.py OUT_PNG"""
import sys
from pathlib import Path
import numpy as np, rasterio, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
A = Path('/mnt/project-files/periyar/EMA_FINAL_REVISION/analysis'); C = Path('/mnt/project-files/periyar/EMA_FINAL_REVISION/evening_holdout_corrected')
def rd(p):
    with rasterio.open(p) as s: t = s.transform; a = s.read(1)
    return a, [t.c / 1e3, (t.c + t.a * a.shape[1]) / 1e3, (t.f + t.e * a.shape[0]) / 1e3, t.f / 1e3]
plt.rcParams.update({'font.size': 9, 'font.family': 'DejaVu Sans'})
fig, ax = plt.subplots(2, 3, figsize=(7.2, 6.4), sharex=True, sharey=True, constrained_layout=True)
cm1 = ListedColormap(['white', '#e6e6e6', '#1f4e9c']); cm2 = ListedColormap(['white', '#e6e6e6', '#1f4e9c', '#e08a1e', '#7b3fa0'])
def cls(a):
    c = np.zeros(a.shape, int); c[a >= 0] = 1; c[a > 0.5] = 2; return c
for row, (lab, mp, mt, et) in enumerate([
        ('original georeferencing', A / 'radar_morning_published_rule.tif', A / 'radar_morning_threshold_only.tif', A / 'radar_evening_threshold_only.tif'),
        ('corrected georeferencing', C / 'radar_morning_published_rule_corrected.tif', C / 'radar_morning_threshold_only_corrected.tif', C / 'radar_evening_threshold_only_corrected.tif')]):
    p, ext = rd(mp); m, _ = rd(mt); e, _ = rd(et)
    ax[row, 0].imshow(cls(p), cmap=cm1, vmin=0, vmax=2, extent=ext, interpolation='nearest')
    ax[row, 1].imshow(cls(e), cmap=cm1, vmin=0, vmax=2, extent=ext, interpolation='nearest')
    c = np.zeros(m.shape, int); v = (m >= 0) & (e >= 0); c[v] = 1; mm, ee = (m > 0.5) & v, (e > 0.5) & v
    c[mm & ee] = 2; c[mm & ~ee] = 3; c[~mm & ee] = 4
    ax[row, 2].imshow(c, cmap=cm2, vmin=0, vmax=4, extent=ext, interpolation='nearest')
    k = 'abc' if row == 0 else 'def'
    ax[row, 0].set_title(f'({k[0]}) Morning, full rule\n{lab}', fontsize=8)
    ax[row, 1].set_title(f'({k[1]}) Evening, threshold rule\n{lab}', fontsize=8)
    ax[row, 2].set_title(f'({k[2]}) Threshold rule, both scenes\n{lab}', fontsize=8)
for a in ax[1]: a.set_xlabel('Easting (km, UTM 43N)')
for a in ax[:, 0]: a.set_ylabel('Northing (km)')
fig.legend(handles=[Patch(color='#1f4e9c', label='water (in c, f: both scenes)'), Patch(color='#e08a1e', label='morning only'),
                    Patch(color='#7b3fa0', label='evening only'), Patch(color='#e6e6e6', label='not water'),
                    Patch(facecolor='white', edgecolor='k', label='no data / channel')],
           loc='lower center', ncol=3, fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.07))
fig.savefig(sys.argv[1], dpi=600, bbox_inches='tight'); print('saved')
