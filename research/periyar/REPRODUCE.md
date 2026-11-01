# Reproducing and reusing the Periyar workflow

This guide covers three tasks: installing the workflow, checking the published numbers
without running the model, and rerunning or adapting the model ensembles. The command
behind each table and figure is listed in `README.md`.

## What is HYDROLIB-core and what is research code

| Layer | Location | Scope |
|---|---|---|
| HYDROLIB-core 1.5.0 | `hydrolib/` (upstream Deltares library, unchanged) | Reads and writes D-Flow FM model files (MDU, network, boundary, friction, external forcing) |
| GSA toolkit | `research/gsa/` | D-Flow FM-domain-independent (any D-Flow FM model that HYDROLIB-core reads; not other model formats): factor definitions, Morris, Sobol and LHS designs, file perturbation through HYDROLIB-core, Morris and Sobol indices. Seven unit tests use the separate Moergestels Broek model in `tests/data` and cover factor generation, file perturbations and known rankings on synthetic outputs; they do not run D-Flow FM |
| Periyar study | `research/periyar/` | Builds the open-data Periyar model, runs the ensembles and scores them against the HFL and the Sentinel-1 map. It also computes the conditional contraction, the drainage diagnosis, the structural tests and the figures |

## 1. Install

```bash
git clone https://github.com/karthice1/HYDROLIB-core.git
cd HYDROLIB-core
git checkout claude/project-thread-hlha70         # or the archived release tag
python3.13 -m venv .venv && . .venv/bin/activate
pip install -e .
pip install -r research/periyar/requirements.txt  # pinned versions of the published runs
python -m pytest research/gsa/tests               # 7 tests, about 5 s
```

`bash research/periyar/release_checks.sh [DATA_DIR]` does all of this in a new virtual
environment, runs the checks of step 2 and writes a report to `research/periyar/test_report/`.
The Zenodo archive contains this report for the archived commit itself;
`research/periyar/zenodo/build_archive.sh` builds the archive, COMMIT.txt and the report from one
commit and verifies the result.

Only steps 3 and 4 need D-Flow FM. The published runs used computational kernel 1.2.184,
built from the Deltares Delft3D source at commit
`9c5ef6b8b0bddb3a2fe3b987338c34757f30b38a`. Pass the executable, or a wrapper script that
sets its library path, as `--dflowfm`. The model is deterministic: rerunning a parameter
set with the same build reproduces its results to the last digit.

## 2. Check the published numbers (no model runs, under a minute)

```bash
python -m research.periyar.check_results
```

This recomputes 120 published values from the per-run results stored in `results/`. Each
check prints PASS or FAIL. The exit code is the number of failures, so it can run in CI.

| Checked automatically (no model runs) | Needs D-Flow FM runs |
|---|---|
| Morris mu* of the three leading factors (Table 4) | The per-run results themselves |
| Accepted counts and contractions at ±0.5, ±1.0, ±1.5 and ±2.0 m, pooled and per design (Table 5) | Morris indices of other outputs, Sobol and LHS ensembles from scratch |
| Every cell of the contingency table and the 73-run score ranges (Table 7) | The radar flood map (needs the Sentinel-1 scenes) |
| Every cell of the forcing-response table, 23 strict runs and 14 central runs (Table 10) | All figures |
| Eight-factor, fixed-factor and pooled results (Table 11) | Structural tests (Table 8) and the drainage diagnosis maps |
| Drainage summary (Section 4.6) | |

`python -m research.periyar.verify_data DATA_DIR` checks the downloaded source data against
`data_manifest.csv` (product, version, tile or scene, download location, processing, size and
SHA-256 for every file).

`results/revision2/supplementary/compare_reruns.py` reproduces the rerun comparisons of
Section 4.7.

## 3. Rerun the study

| Step | Command (see `README.md` for arguments) | Runs | Time on 4 cores |
|---|---|---|---|
| Build model v5 and radar map | `build_base_model`, `sar_flood_map` | – | under 1 h |
| Morris screening | `sensitivity --method morris --trajectories 20` | 180 | about 12 h |
| Conditioning ensembles | `sensitivity --method lhs --samples 200 --seed 42` and `--samples 400 --seed 43` | 600 | about 1.5 days |
| Contraction, bootstrap, weighting | `uncertainty OUT RUNS_LHS RUNS_LHS2` | – | minutes |
| Supplementary ensembles | `results/revision2/supplementary/chain_rev[2-6].sh` | 751 | about 2.7 days |
| Drainage diagnosis | `drainage RUN_DIR SAR.tif OUT` | 2 | under 1 h |

Every `sensitivity` call resumes from its `results.csv` when it is restarted. Runs take
4 to 60 min each; low roughness makes runs slower.

## 4. Use the workflow on another D-Flow FM domain

The method is not tied to the Periyar. Using it on another domain needs four
case-specific inputs:

This has not yet been tried end to end on a second domain; the unit tests cover the reusable
layer on a separate D-Flow FM model but do not run a hydrodynamic simulation there.

1. **A base model.** Any D-Flow FM model that HYDROLIB-core can read. `build_base_model.py`
   shows one way to build one from open global data, but an existing model can be used
   directly.
2. **Uncertain factors.** Define them as `research.gsa.factors.Factor` entries: multipliers
   on friction, boundary forcing or structures, or offsets. `research/periyar/sensitivity.py`
   shows how to add model-specific factors such as channel depth or a canopy correction.
3. **Observation operators.** `validate.py` reduces each run to a peak level near a gauge and
   a contingency table against a flood map. Replace the gauge coordinates, the HFL
   (`HFL_NEELESWARAM`) and the radar map with local equivalents. A high-water mark, a gauge
   record or another satellite scene can each be the constraint.
4. **An acceptance rule and tolerances.** `uncertainty.accepted()` holds them. The
   contraction (Eq. 3), the bootstrap, the tolerance grid and the likelihood weighting in
   `uncertainty.py` work on any table of runs with an output column and the score columns.

The drainage diagnosis (`drainage.py`) needs only a D-Flow FM map file, the bed level and
a flood map, so it applies to any 2D model that keeps water on the land longer than
observed.
