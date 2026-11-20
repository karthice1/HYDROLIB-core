"""[corrected-map version: python -I identifiability_corrected.py runs_corrected.csv OUT_DIR] Post-hoc analysis of the stored 600-run pooled ensemble (no new model runs).

Reads research/periyar/results/revision2/uncertainty_600_tol1.0/runs.csv and writes
identifiability.csv, equifinality.csv and metric_correlations.csv next to this script.
(a) Regional sensitivity: two-sample Kolmogorov-Smirnov distance between the factor
    values of accepted and rejected runs, for each acceptance test separately.
(b) Equifinality: Spearman correlation between factors within the accepted set
    (the prior design is uncorrelated by construction, so structure here is induced by conditioning).
(c) How the scores co-vary with peak flooded land across all 600 runs.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
RUNS = Path(sys.argv[1]) if len(sys.argv) > 1 else Path('/home/user/hydrolib-core/research/periyar/results/revision2/uncertainty_600_tol1.0/runs.csv')
OUT = Path(sys.argv[2])
r = pd.read_csv(RUNS)
F = ['n_floodplain', 'channel_depth', 'canopy_factor', 'q_rating', 'q_idukki']
radar = (r.model_flooded_km2 >= 22.96) & (r.open_hit_rate >= 0.75)  # corrected-map radar water on land; non-binding
sets = {}
for tol in (0.5, 1.0, 2.0):
    h = r.hfl_error_m.abs() <= tol
    sets[f'hfl_only_{tol}'] = h; sets[f'both_{tol}'] = h & radar
# radar-only set not recomputable: hit rates of runs with |HFL error| > 2 m are not available
rows = []
for name, m in sets.items():
    for f in F:
        a, b = r.loc[m, f], r.loc[~m, f]
        ks = stats.ks_2samp(a, b)
        lo, hi = r[f].min(), r[f].max()
        rows.append(dict(test=name, n_accepted=int(m.sum()), factor=f, ks_D=ks.statistic, ks_p=ks.pvalue,
                         prior_min=lo, prior_max=hi, acc_min=a.min(), acc_max=a.max(),
                         acc_q05=a.quantile(.05), acc_q95=a.quantile(.95),
                         acc_5_95_share_of_prior=(a.quantile(.95) - a.quantile(.05)) / (hi - lo)))
pd.DataFrame(rows).to_csv(OUT / 'identifiability.csv', index=False)
eq = []
for name in ('both_1.0', 'both_0.5', 'hfl_only_1.0'):
    m = sets[name]; c = r.loc[m, F].corr(method='spearman')
    for i, a in enumerate(F):
        for b in F[i + 1:]:
            rho, p = stats.spearmanr(r.loc[m, a], r.loc[m, b])
            prho, _ = stats.spearmanr(r[a], r[b])
            eq.append(dict(set=name, n=int(m.sum()), factor_a=a, factor_b=b, spearman_accepted=rho, p_value=p, spearman_prior=prho))
pd.DataFrame(eq).to_csv(OUT / 'equifinality.csv', index=False)
mc = []
for col in ['hfl_error_m', 'neeleswaram_peak_level_m', 'model_flooded_km2', 'open_hit_rate', 'open_false_alarm_ratio',
            'all_critical_success_index', 'open_critical_success_index']:
    for name, m in (('hfl_within_2m_306', (r.hfl_error_m.abs() <= 2.0).values), ('both_1.0', sets['both_1.0'].values)):
        rho, p = stats.spearmanr(r.loc[m, 'max_flooded_area_km2'], r.loc[m, col])
        mc.append(dict(set=name, metric=col, spearman_with_peak_area=rho, p_value=p,
                       min=r.loc[m, col].min(), max=r.loc[m, col].max()))
pd.DataFrame(mc).to_csv(OUT / 'metric_correlations.csv', index=False)
print('done')
