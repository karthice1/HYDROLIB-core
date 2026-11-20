"""Rerun lhs400 run 70 (the median accepted run at +-0.5 m on the corrected map) with unchanged samples,
model, forcing and code, KEEPING the map file, so that research/periyar/drainage.py can be applied to it
(as was done for the old median run). Resumable: skips D-Flow FM if the map file exists."""
import sys, os, subprocess
from pathlib import Path
import pandas as pd
sys.path.insert(0, '/home/user/hydrolib-core')
from research.periyar import sensitivity as S
BASE = Path('/home/user/ext/periyar_v5'); DFM = '/home/user/ext/run_dflowfm.sh'
rd = Path('/home/user/ext/median70/run_0070')
samples = pd.read_csv('/home/user/hydrolib-core/research/periyar/results/v5_lhs400/samples.csv', index_col=0)
mf = rd / 'DFM_OUTPUT_periyar' / 'periyar_map.nc'
if not (mf.exists() and (rd / 'DONE').exists()):
    S.prepare(BASE, rd, samples.loc[70], 'additive')
    with (rd / 'run.log').open('w') as f:
        rc = subprocess.run([DFM, '--autostartstop', 'periyar.mdu'], cwd=rd, stdout=f, stderr=subprocess.STDOUT,
                            env={**os.environ, 'OMP_NUM_THREADS': '4'}).returncode
    assert rc == 0 and mf.exists(), rc
    (rd / 'DONE').write_text('ok')
print('dflowfm done', flush=True)
