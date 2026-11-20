"""Every radar-dependent number of the manuscript, recomputed on the CORRECTED radar map from stored results
(rescoring/runs_corrected.csv, lhs8/lhs8b/fixed corrected tables, rerun_radar forcing results and
radar_state.npz files). Acceptance rule, tolerances, factor ranges and code are unchanged
(research.periyar.uncertainty.accepted / narrowing). Writes rescoring/manuscript_numbers.json.
Usage: python -I manuscript_numbers.py"""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import ks_2samp, spearmanr
sys.path.insert(0, '/home/user/hydrolib-core')
from research.periyar.uncertainty import accepted, narrowing, OUTPUT
from research.periyar.radar_check import confusion
from research.periyar.validate import OPEN_CLASSES, _sample
EF = Path('/mnt/project-files/periyar/EMA_FINAL_REVISION'); RS = EF / 'rescoring'
RR = Path('/home/user/ext/rerun_radar'); REPO = Path('/home/user/hydrolib-core/research/periyar/results')
SUP = REPO / 'revision2/supplementary'; MAP = Path('/home/user/ext/sar_flood_20180821_corrected.tif')
FAC = ['n_floodplain', 'channel_depth', 'canopy_factor', 'q_rating', 'q_idukki']
PRIOR = {'n_floodplain': (0.6, 1.5), 'channel_depth': (2, 8), 'canopy_factor': (0.3, 0.9), 'q_rating': (0.85, 1.15), 'q_idukki': (0.8, 1.2)}
out = {}
r = lambda v, k=2: None if v is None else round(float(v), k)
m = pd.read_csv(RS / 'runs_corrected.csv')
assert len(m) == 600
A = m[OUTPUT].to_numpy()
def C(df, mask):
    n = int(mask.sum())
    if n < 2: return dict(n=n)
    d = narrowing(df[OUTPUT].to_numpy(), df.loc[mask, OUTPUT].to_numpy())
    return dict(n=n, low=r(d['low'], 1), high=r(d['high'], 1), prior_low=r(d['prior_low'], 1), prior_high=r(d['prior_high'], 1), C=r(d['narrowing'], 3))
radar = lambda df: accepted(df, hfl_tol=np.inf)
# --- per design
for des in ('lhs200', 'lhs400'):
    d = m[m.design == des].reset_index(drop=True)
    for t in (0.5, 1.0):
        k = accepted(d, hfl_tol=t)
        out[f'design_{des}_tol{t}'] = {**C(d, k), **{f'{f}_range': [r(d.loc[k, f].min()), r(d.loc[k, f].max())] for f in ('channel_depth', 'n_floodplain')}}
# --- pooled per tolerance
for t in (0.5, 1.0, 1.5, 2.0):
    k = accepted(m, hfl_tol=t); h = m.hfl_error_m.abs() <= t
    out[f'pooled_tol{t}'] = dict(both=C(m, k), hfl_only=C(m, h.values), radar_only=C(m, radar(m)),
                                 accepted_5_95=[r(np.quantile(A[k], .05), 1), r(np.quantile(A[k], .95), 1)],
                                 hfl_err_range=[r(m.hfl_error_m[k].min()), r(m.hfl_error_m[k].max())],
                                 volume_range=[r(m.max_flood_volume_mcm[k].min(), 0), r(m.max_flood_volume_mcm[k].max(), 0)])
rad = radar(m); rej = ~rad
out['radar_rejects'] = dict(n=int(rej.sum()), peak_range=[r(A[rej].min(), 1), r(A[rej].max(), 1)],
                            of_hfl05=int((rej & (m.hfl_error_m.abs() <= 0.5)).sum()), hfl05=int((m.hfl_error_m.abs() <= 0.5).sum()),
                            min_accepted_radar_peak=r(A[rad].min(), 1))
out['prior_peak'] = dict(range=[r(A.min(), 1), r(A.max(), 1)], q05_50_95=[r(x, 1) for x in np.quantile(A, [.05, .5, .95])])
# --- identifiability at +-1.0
k = accepted(m, hfl_tol=1.0); h = (m.hfl_error_m.abs() <= 1.0).values
ident = {}
for f in FAC:
    lo, hi = PRIOR[f]; v = m[f]
    ident[f] = dict(D_hfl=r(ks_2samp(v[h], v[~h]).statistic), D_radar=r(ks_2samp(v[rad], v[~rad]).statistic),
                    D_both=r(ks_2samp(v[k], v[~k]).statistic), p_both=r(ks_2samp(v[k], v[~k]).pvalue, 3),
                    share=r((np.quantile(v[k], .95) - np.quantile(v[k], .05)) / (hi - lo)),
                    acc_range=[r(v[k].min()), r(v[k].max())], acc_median=r(v[k].median()))
n1, n2 = k.sum(), (~k).sum(); ident['D_crit_5pct_both'] = r(1.358 * np.sqrt((n1 + n2) / (n1 * n2)))
out['identifiability_tol1.0'] = ident
sp = lambda a, b, msk=None: r(spearmanr(m[a][msk] if msk is not None else m[a], m[b][msk] if msk is not None else m[b])[0])
k05 = accepted(m, hfl_tol=0.5); h05 = (m.hfl_error_m.abs() <= 0.5).values
out['correlations'] = dict(depth_rough_prior=sp('channel_depth', 'n_floodplain'), depth_rough_acc10=sp('channel_depth', 'n_floodplain', k),
    depth_rough_acc05=sp('channel_depth', 'n_floodplain', k05), depth_canopy_acc10=sp('channel_depth', 'canopy_factor', k),
    depth_canopy_acc05=sp('channel_depth', 'canopy_factor', k05), depth_canopy_hflonly10=sp('channel_depth', 'canopy_factor', h),
    depth_hit=sp('channel_depth', 'open_hit_rate'), canopy_hit=sp('canopy_factor', 'open_hit_rate'),
    depth_hflerr=sp('channel_depth', 'hfl_error_m'), rough_hflerr=sp('n_floodplain', 'hfl_error_m'),
    level_area=sp('neeleswaram_peak_level_m', OUTPUT),
    pearson_hit_area=r(np.corrcoef(m.open_hit_rate, A)[0, 1]), pearson_opencsi_area=r(np.corrcoef(m.open_critical_success_index, A)[0, 1]),
    pearson_allcsi_area=r(np.corrcoef(m.all_critical_success_index, A)[0, 1]),
    spearman_allcsi_area=sp('all_critical_success_index', OUTPUT), spearman_hit_area=sp('open_hit_rate', OUTPUT),
    spearman_openfar_area=sp('open_false_alarm_ratio', OUTPUT), spearman_opencsi_area=sp('open_critical_success_index', OUTPUT))
out['spatial_600'] = {c: [r(m[c].min(), 3), r(m[c].max(), 3)] for c in ('all_critical_success_index', 'open_false_alarm_ratio', 'open_critical_success_index', 'open_hit_rate')}
for t, kk in ((1.0, k), (0.5, k05)):
    out[f'spatial_acc{t}'] = {c: [r(m.loc[kk, c].min(), 3), r(m.loc[kk, c].max(), 3)] for c in ('all_critical_success_index', 'open_false_alarm_ratio', 'open_critical_success_index', 'open_hit_rate')}
# --- factor range cut by a fifth
cut = []
for f in FAC:
    lo, hi = PRIOR[f]
    for side, keep in (('low', m[f] >= lo + 0.2 * (hi - lo)), ('high', m[f] <= hi - 0.2 * (hi - lo))):
        d = m[keep].reset_index(drop=True); c = C(d, accepted(d, hfl_tol=1.0)); c.update(factor=f, side=side); cut.append(c)
out['range_cut_tol1.0'] = dict(cases=cut, n=[min(c['n'] for c in cut), max(c['n'] for c in cut)], C=[min(c['C'] for c in cut), max(c['C'] for c in cut)])
# --- 200-run subsets of the pool, C at +-1.0
rng = np.random.default_rng(1); cs, ns = [], []
for i in range(5000):
    s = m.iloc[rng.choice(600, 200, replace=False)].reset_index(drop=True); c = C(s, accepted(s, hfl_tol=1.0))
    if 'C' in c: cs.append(c['C']); ns.append(c['n'])
cs = np.array(cs)
out['subsets200_tol1.0'] = dict(C_5_50_95=[r(x) for x in np.quantile(cs, [.05, .5, .95])], n_5_95=[int(x) for x in np.quantile(ns, [.05, .95])])
# --- eight-factor designs (radar scores exist for every run with |HFL error| <= 1.0 m)
l8s = {}
for name in ('lhs8', 'lhs8b'):
    d = pd.read_csv(RS / f'{name}_runs_corrected.csv')
    have = d.open_hit_rate.notna()
    assert have[d.hfl_error_m.abs() <= 1.0].all()
    dd = d.fillna({'open_hit_rate': -1, 'model_flooded_km2': -1})   # only used where |HFL| > 1, i.e. rejected anyway
    res = {}
    for t in (0.5, 1.0):
        kk = accepted(dd, hfl_tol=t); res[f'tol{t}'] = C(dd, kk)
        res[f'tol{t}']['depth_rough'] = r(spearmanr(dd.channel_depth[kk], dd.n_floodplain[kk])[0]) if kk.sum() > 2 else None
        res[f'hfl_only_tol{t}'] = C(dd, (dd.hfl_error_m.abs() <= t).values)
        # bootstrap: rule applied to resamples; valid because runs without radar scores fail the HFL test
        bs = []
        rg = np.random.default_rng(1)
        for b in range(5000):
            s = dd.iloc[rg.integers(0, len(dd), len(dd))].reset_index(drop=True); c = C(s, accepted(s, hfl_tol=t))
            if 'C' in c: bs.append((c['C'], c['n']))
        bs = np.array(bs); res[f'bootstrap_tol{t}'] = dict(C=[r(np.quantile(bs[:, 0], .05)), r(np.quantile(bs[:, 0], .95))], n=[int(np.quantile(bs[:, 1], .05)), int(np.quantile(bs[:, 1], .95))])
        res[f'tol{t}']['factor_ranges'] = {f: [r(dd.loc[kk, f].min()), r(dd.loc[kk, f].max())] for f in ('n_channel', 'q_idamalayar', 'tide_offset')}
    # observation operators at +-1.0 (only where all runs passing the operator have radar scores)
    ops = {}
    for col in ('hfl_err_nearest', 'hfl_err_r300_median', 'hfl_err_r750_median', 'hfl_err_r1000_median', 'hfl_err_r500_max', 'hfl_err_r500_min'):
        pas = d[col].abs() <= 1.0
        complete = bool(have[pas].all())
        if complete:
            dx = dd.copy(); dx['hfl_error_m'] = d[col]; kk = accepted(dx, hfl_tol=1.0)
            ops[col] = dict(n=int(kk.sum()), range=[r(dx.loc[kk, OUTPUT].min(), 1), r(dx.loc[kk, OUTPUT].max(), 1)])
        else:
            ops[col] = dict(computable=False, runs_without_radar=int((pas & ~have).sum()))
    res['operators_tol1.0'] = ops
    l8s[name] = (dd, res); out[f'eight_{name}'] = res
pool = pd.concat([l8s['lhs8'][0], l8s['lhs8b'][0]], ignore_index=True); kk = accepted(pool, hfl_tol=1.0)
out['eight_pooled_tol1.0'] = C(pool, kk)
n = out['eight_lhs8']['tol1.0']['n']; cov = (n - 1) / (n + 1); lo, hi = np.quantile(A, [(1 - cov) / 2, (1 + cov) / 2])
e = out['eight_lhs8']['tol1.0']; out['lhs8_vs_five_factor_prior'] = r(1 - (e['high'] - e['low']) / (hi - lo))
# --- likelihood weighting (from the unchanged script's output)
lw = pd.read_csv(RS / 'uncertainty_corrected_tol1.0/likelihood_weights.csv')
out['likelihood'] = lw.round(3).to_dict(orient='records')
# --- threshold grid at +-1.0 and +-0.5
tg = pd.read_csv(RS / 'uncertainty_corrected_tol1.0/threshold_grid.csv')
out['threshold_grid'] = tg[tg.hfl_tolerance_m.isin([0.5, 1.0])].round(3).to_dict(orient='records')
# --- fixed-factor variants around the old 23 strict runs
fx = pd.read_csv(RS / 'fixed_runs_corrected.csv'); par = pd.read_csv(SUP / 'fixed_parents.csv')
fx = fx.merge(pd.read_csv(SUP / 'samples_fixed.csv'), on='run', suffixes=('', '_s')).merge(par, on='run')
fd = fx.fillna({'open_hit_rate': -1, 'model_flooded_km2': -1})
for t in (0.5, 1.0):
    kk = accepted(fd, hfl_tol=t)
    out[f'fixed_tol{t}'] = dict(n=int(kk.sum()), range=[r(fd.loc[kk, OUTPUT].min(), 1), r(fd.loc[kk, OUTPUT].max(), 1)],
                                fail_radar=int(((fd.hfl_error_m.abs() <= t) & ~kk).sum()), fail_hfl=int((fd.hfl_error_m.abs() > t).sum()),
                                parents_with_accepted=int(fd[kk].parent.nunique()),
                                n_channel_mean_acc=r(fd.n_channel[kk].mean()), n_channel_mean_rej=r(fd.n_channel[~kk].mean()),
                                depth_rough=r(spearmanr(fd.channel_depth[kk], fd.n_floodplain[kk])[0]))
# --- the 23 old strict runs and 14 central runs: acceptance under original forcing and each variant (corrected map)
samp = pd.read_csv(SUP / 'samples_accepted.csv'); cen = pd.read_csv(SUP / 'samples_central.csv')
def pooled_index(s):
    idx = []
    for _, row in s.iterrows():
        j = np.where(np.isclose(m.channel_depth, row.channel_depth) & np.isclose(m.n_floodplain, row.n_floodplain))[0]
        assert len(j) == 1; idx.append(int(j[0]))
    return idx
forcing = {}
for grp, s, pref in (('strict23', samp, 'strict'), ('central14', cen, 'central')):
    idx = pooled_index(s); base = m.iloc[idx].reset_index(drop=True)
    forcing[grp] = dict(pooled=[f"{m.design[i]}:{m.run[i]}" for i in idx])
    for t in (0.5, 1.0):
        forcing[grp][f'base_accepted_tol{t}'] = int(accepted(base, hfl_tol=t).sum())
        forcing[grp][f'base_fail_radar_tol{t}'] = int(((base.hfl_error_m.abs() <= t) & ~accepted(base, hfl_tol=t)).sum())
    for var in ('mult', 'lagshort', 'laglong'):
        v = pd.read_csv(RR / f'{var}_{pref}/results.csv').sort_values('run').reset_index(drop=True)
        v = v.rename(columns={f'corrected_{c}': c for c in ('hfl_error_m', 'model_flooded_km2', 'open_hit_rate', 'all_critical_success_index', 'open_false_alarm_ratio')})
        res = {}
        for t in (0.5, 1.0):
            a0, a1 = accepted(base, hfl_tol=t), accepted(v, hfl_tol=t)
            res[f'tol{t}'] = dict(accepted=int(a1.sum()), dropouts=int((a0 & ~a1).sum()), gains=int((~a0 & a1).sum()),
                                  dropout_reason_hfl=int((a0 & ~a1 & (v.hfl_error_m.abs() > t)).sum()),
                                  range=[r(v.loc[a1, OUTPUT].min(), 1), r(v.loc[a1, OUTPUT].max(), 1)] if a1.any() else None)
        dA = (v[OUTPUT] - base[OUTPUT]).abs(); de = (v.hfl_error_m - base.hfl_error_m).abs()
        res.update(dA_lt1=int((dA < 1).sum()), max_dA=r(dA.max(), 0), de_lt001=int((de < 0.01).sum()), max_de=r(de.max()))
        forcing[grp][var] = res
out['forcing'] = forcing
# --- spatial scores from stored per-cell states
st_root = Path('/mnt/project-files/periyar/EMA_FINAL_REVISION/baseline_2026-11-17/model_states')
def state(d, i):
    p = RR / f'{d}/run_{i:04d}/radar_state.npz'
    if p.exists(): return p
    q = st_root / f'gsa_v5_{d}__time_series__run_{i:04d}/radar_state.npz'   # 7 runs kept in the original study
    assert q.exists(), (d, i); return q
def score_state(npz, tif=MAP, cut=0.5):
    s = np.load(npz); f = _sample(Path(tif), s['x'], s['y'], -1.0); land = s['land'] & (f >= 0)
    op = land & np.isin(s['landcover'], OPEN_CLASSES); ob = f > cut
    return {dom: {k: (r(v, 3) if isinstance(v, float) else v) for k, v in confusion(s['model_wet'][msk], ob[msk]).items()} for dom, msk in (('all', land), ('open', op))}
out['median_run_tol0.5'] = dict(id='lhs400:70', peak=r(m[(m.design == 'lhs400') & (m.run == 70)][OUTPUT].iloc[0], 2),
                                 hfl_err=r(m[(m.design == 'lhs400') & (m.run == 70)].hfl_error_m.iloc[0]),
                                 scores=score_state(RR / 'lhs400/run_0070/radar_state.npz'))
acc_ids = [(m.design[i], int(m.run[i])) for i in np.where(k)[0]]
out['median_tol1.0_peak'] = r(np.median(A[k]), 1)
V = EF / 'radar_corrected/sar_variants'
var_csi = {}
for name in ('base', 'thr_m12.25', 'thr_m15.25', 'drop2', 'drop4'):
    for cutv in (0.5, 0.3):
        if name != 'base' and cutv == 0.3: continue
        vals = [score_state(state(d, i), V / f'{name}.tif', cutv)['all']['csi'] for d, i in acc_ids]
        var_csi[f'{name}_cut{cutv}'] = [min(vals), max(vals)]
out['variant_csi_acc69_all_land'] = var_csi
out['variant_csi_all'] = [min(v[0] for v in var_csi.values()), max(v[1] for v in var_csi.values())]
# bias of accepted runs (all land)
bias = [score_state(state(d, i))['all']['bias'] for d, i in acc_ids]
out['bias_acc69_all_land'] = [min(bias), max(bias)]
# stored structural / nominal states
for p in sorted(st_root.rglob('radar_state.npz')):
    out.setdefault('stored_states', {})[p.parent.name] = score_state(p)
# --- the 7 runs of the first design accepted at +-0.5 m (time-series runs)
d200 = m[m.design == 'lhs200']; k2 = accepted(d200.reset_index(drop=True), hfl_tol=0.5)
out['lhs200_tol0.5_accepted_runs'] = d200.run.to_numpy()[k2].tolist()
out['lhs200_tol1.0_accepted_ts_runs'] = [int(x) for x in [26, 61, 64, 86, 93, 119, 188] if accepted(d200.reset_index(drop=True), hfl_tol=1.0)[d200.run.to_numpy().tolist().index(x)]]
json.dump(out, open(RS / 'manuscript_numbers.json', 'w'), indent=1, default=float)
print(json.dumps(out, indent=1, default=float))
