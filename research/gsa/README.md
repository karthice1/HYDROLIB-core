# Global sensitivity analysis for D-Flow FM

Research workflow for the Periyar (Kerala) August 2018 flood study. It uses
HYDROLIB-core to write a perturbed copy of a calibrated D-Flow FM model for
every sample of a Morris or Sobol design, and SALib to compute sensitivity
indices from the simulated outputs.

| Module | Purpose |
|---|---|
| `factors.py` | Uncertain inputs, their ranges, and the provisional Periyar 2018 factor set |
| `sampling.py` | Morris and Sobol (Saltelli) sample designs |
| `perturb.py` | Edits to friction, boundary forcing and structure files |
| `ensemble.py` | Copies the base model per sample, applies edits, writes `manifest.csv` |
| `analysis.py` | Morris (mu*, sigma) and Sobol (S1, ST) indices for a scalar output |

Roughness factors are multipliers on bed roughness: a value above 1 means a
rougher bed for every friction type (for Chezy, Strickler and De Bos-Bijkerk the
coefficient is divided, for Manning and Nikuradse it is multiplied).

## Usage

```python
from pathlib import Path
from research.gsa import ensemble, sampling
from research.gsa.factors import PERIYAR_2018

samples = sampling.morris(PERIYAR_2018, trajectories=20)
ensemble.build(Path("base_model"), Path("runs"), samples, PERIYAR_2018, {
    "friction": ["roughness/channel.ini", "roughness/floodplain.ini"],
    "forcing": ["boundaries.bc"],
    "structure": ["structures.ini"],
})
```

Then run D-Flow FM in each `runs/run_*` directory, reduce each map output to a
scalar (for example flooded area or peak water level at a gauge), and pass the
resulting series to `analysis.morris` or `analysis.sobol`.

## Tests

```bash
pip install -e . SALib pytest
python -m pytest research/gsa/tests
```

The tests use the Moergestels Broek model from `tests/data` as a stand-in until
the Periyar base model exists.
