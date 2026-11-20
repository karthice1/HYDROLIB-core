"""Pack every per-cell model state used for the corrected radar scores into one compact file per group.

Each state holds, for the 32,337 cells of the 200 m grid, whether the cell is land (dry at the start) and whether it
is wet at the radar pass (21 Aug 2018 06:10 IST). x, y and WorldCover class are common to all runs and stored once.
Nothing is recomputed; the bits are copied from the radar_state.npz files written by the reruns (and, for 7 lhs200
runs and the structural variants, from the states kept on 17 Nov 2026). Usage: python -I pack_states_v11.py OUT_DIR"""
import sys, glob, json, hashlib
from pathlib import Path
import numpy as np
OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
RR = Path('/home/user/ext/rerun_radar'); BASE = Path('/mnt/project-files/periyar/EMA_FINAL_REVISION/baseline_2026-11-17/model_states')
ref = np.load(sorted(RR.glob('lhs400/run_*/radar_state.npz'))[0])
np.savez_compressed(OUT / 'grid.npz', x=ref['x'], y=ref['y'], landcover=ref['landcover'])
groups = {g.name: sorted(g.glob('run_*/radar_state.npz')) for g in sorted(RR.iterdir()) if g.is_dir()}
groups['lhs200'] += sorted(BASE.glob('gsa_v5_lhs200__time_series__run_*/radar_state.npz'))
manifest = {}
for g, files in groups.items():
    if not files: continue
    names, land, wet, src = [], [], [], []
    for f in files:
        d = np.load(f)
        for k in ('x', 'y', 'landcover'):
            assert np.array_equal(d[k], ref[k]), (f, k)
        name = f.parent.name.split('__')[-1]
        names.append(name); land.append(np.packbits(d['land'])); wet.append(np.packbits(d['model_wet']))
        src.append('rerun' if str(f).startswith(str(RR)) else 'stored_2026-11-17')
    order = np.argsort(names)
    p = OUT / f'states_{g}.npz'
    np.savez_compressed(p, run=np.array(names)[order], land=np.array(land)[order], model_wet=np.array(wet)[order],
                        source=np.array(src)[order], n_cells=len(ref['x']))
    manifest[g] = dict(runs=len(names), file=p.name, sha256=hashlib.sha256(p.read_bytes()).hexdigest())
    print(g, len(names), flush=True)
# structural variants (own grids) are kept as full states
for f in sorted(BASE.glob('structural__*/radar_state.npz')) + sorted(Path('/home/user/ext/median70').glob('run_*/radar_state.npz')):
    d = np.load(f); name = f.parent.name.replace('structural__', 'structural_') if 'structural' in str(f) else 'median_lhs400_' + f.parent.name
    p = OUT / f'state_{name}.npz'; np.savez_compressed(p, **{k: d[k] for k in d})
    manifest[name] = dict(runs=1, file=p.name, sha256=hashlib.sha256(p.read_bytes()).hexdigest()); print(name, flush=True)
json.dump(manifest, open(OUT / 'states_manifest.json', 'w'), indent=1)
