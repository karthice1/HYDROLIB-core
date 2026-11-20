"""Manuscript figures from the CORRECTED radar map and recomputed results, drawn with the unchanged functions
of research/periyar/figures.py; only the output resolution is raised to 600 dpi (savefig wrapper).
Manuscript Fig. 1 = fig1, Fig. 2 = fig_uncertainty, Fig. 3 = fig2, Fig. 4 = fig3, Fig. 5 = fig_confusion,
Fig. 6 = fig_drainage. Usage: python -I make_figures_corrected.py OUT_DIR"""
import sys
from pathlib import Path
import matplotlib; matplotlib.use('Agg')
from matplotlib.figure import Figure
_save = Figure.savefig
def save600(self, *a, **k): k['dpi'] = 600; return _save(self, *a, **k)
Figure.savefig = save600
sys.path.insert(0, '/home/user/hydrolib-core')
from research.periyar import figures as F
OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
EF = Path('/mnt/project-files/periyar/EMA_FINAL_REVISION'); MAP = Path('/home/user/ext/sar_flood_20180821_corrected.tif')
ST = EF / 'baseline_2026-11-17/model_states'; R2 = Path('/home/user/hydrolib-core/research/periyar/results')
U1, U05 = EF / 'rescoring/uncertainty_corrected_tol1.0', EF / 'rescoring/uncertainty_corrected_tol0.5'
F.fig1(Path('/home/user/ext/periyar_v5'), MAP, OUT / 'fig1.png', ' v5')
F.fig_uncertainty(U1, OUT / 'fig2.png')
F.fig2(U1 / 'runs.csv', U05 / 'runs.csv', OUT / 'fig3.png')
F.fig3(R2 / 'v5_lhs200/behavioural_time_series.csv', OUT / 'fig4.png', sar_km2=22.84)
F.fig_confusion({'(a) Nominal': ST / 'structural__nominal/radar_state.npz',
                 '(b) Median accepted run': Path('/home/user/ext/rerun_radar/lhs400/run_0070/radar_state.npz'),
                 '(c) Bare-earth proxy': ST / 'structural__min_terrain/radar_state.npz'}, OUT / 'fig5.png', MAP)
F.fig_drainage({'(a) Nominal model': R2 / 'revision2/drainage_nominal/drainage_state.npz',
                '(b) Depressions filled': R2 / 'revision2/drainage_filled/drainage_state.npz'}, OUT / 'fig6.png')
print('ok')
