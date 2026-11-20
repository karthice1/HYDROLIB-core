"""Write corrected copies of runs_with_behaviour.csv (radar scores replaced by corrected-map scores from the reruns;
every other column, including HFL errors and areas, unchanged) so that the unchanged uncertainty.py can be rerun.
Originals are not touched. Usage: python -I build_corrected_inputs.py"""
import json, glob
from pathlib import Path
import pandas as pd, numpy as np
K = ['all_critical_success_index', 'open_critical_success_index', 'open_hit_rate', 'open_false_alarm_ratio']
for src, rr, out in (('gsa_v5_lhs200', 'lhs200', 'lhs200'), ('gsa_v5_lhs400', 'lhs400', 'lhs400'), ('rev2/lhs8', 'lhs8', 'lhs8'), ('rev2/lhs8b', 'lhs8b', 'lhs8b')):
    r = pd.read_csv(f'/home/user/ext/{src}/runs_with_behaviour.csv'); idcol = r.columns[0] if 'run' not in r.columns else 'run'
    r[K] = np.nan; n = 0
    for i, row in r.iterrows():
        f = Path(f'/home/user/ext/rerun_radar/{rr}/run_{int(row[idcol]):04d}/result.json')
        if f.exists():
            v = json.load(open(f)); r.loc[i, K] = [v[f'corrected_{k}'] for k in K]; n += 1
    if rr == 'lhs200':   # 7 runs rescored from stored states
        st = pd.read_csv('/mnt/project-files/periyar/EMA_FINAL_REVISION/radar_corrected/rescored_stored_states.csv')
        st = st[(st.radar_map == 'corrected') & st.state.str.startswith('gsa_v5_lhs200')]
        for name, g in st.groupby('state'):
            run = int(name[-4:]); i = r.index[r[idcol] == run][0]
            if pd.isna(r.loc[i, 'open_hit_rate']):
                g = g.set_index('domain'); r.loc[i, K] = [g.loc['all_land', 'csi'], g.loc['open_land', 'csi'], g.loc['open_land', 'hit_rate'], g.loc['open_land', 'false_alarm_ratio']]; n += 1
    Path(f'/home/user/ext/corrected_inputs/{out}').mkdir(exist_ok=True)
    r.to_csv(f'/home/user/ext/corrected_inputs/{out}/runs_with_behaviour.csv', index=False)
    print(out, len(r), 'with corrected scores:', n)
