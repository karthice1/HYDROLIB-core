# v11_scripts: provenance and reconstruction scripts

These are the scripts that produced the v1.1 corrected radar map, the 884 reruns and the corrected-map
rescoring, kept **exactly as they were run** (byte-identical to the local analysis repository at commit
`ac63eb4`). They are a record of what was executed, not portable tools.

**Every script and batch file in this folder contains absolute paths of the computing environment in which
it ran** (`/home/user/ext/...`, `/home/user/hydrolib-core`, `/mnt/project-files/periyar/...`). They were
not rewritten, so that the record stays exact. To run one elsewhere, edit the path constants at the top of
the script to point at your own copies of the inputs. Several inputs are not in this archive (see the last
column): the D-Flow FM model directory, the D-Flow FM launcher, the raw rerun output folders and the
Sentinel-1 files, which are listed with checksums in `../data_manifest.csv`.

The portable commands of v1.1 are `../check_results_v11.py` and `../release_checks_v11.sh`; they use paths
relative to the repository and are tested in a new environment (the report is in the Zenodo archive, `test_report_v11/`). See
`../REPRODUCE_v11.md`.

| Script | What it did | Inputs outside this archive |
|---|---|---|
| `sar_georef_correction.py` | Estimated the shift of each Sentinel-1 scene against JRC occurrence and rebuilt the flood map with the frozen v1.0 rules. Takes `OUT_DIR` as argument; data and model paths are fixed in the script. | Sentinel-1 scenes and annotations, JRC occurrence, WorldCover, model grid |
| `sar_variants_corrected.py` | The five classification variants of the corrected map | as above |
| `compare_maps.py` | Before/after comparison of the uncorrected and corrected maps | uncorrected map, nominal state |
| `rerun_for_radar.py` | Reran stored parameter sets with D-Flow FM 1.2.184 to keep the map output at the radar pass | model directory `periyar_v5`, D-Flow FM launcher |
| `run_all.sh`, `run_batch2.sh`, `run_batch3.sh`, `run_batch4.sh` | The four rerun batches (884 reruns); run lists in `need_*.txt` and `rest_*.txt` | as above |
| `rerun_median_full.py` | Full-output rerun of lhs400 run 70 for the drainage diagnosis | as above |
| `verify_294.py` | Verified the 294 batch-4 runs from their outputs | raw rerun folders |
| `rerun_reproduction_v11.py` | Compared each rerun with the stored v1.0 result (`../results_v11/rerun_reproduction.csv`) | raw rerun folders |
| `pack_states_v11.py` | Packed the per-cell states into `../results_v11/model_states/` | raw rerun folders |
| `rescore_stored.py` | Rescored the stored v1.0 per-cell states against both maps | stored v1.0 state files |
| `build_corrected_inputs.py` | Wrote corrected copies of the per-run tables | raw rerun folders |
| `recompute_acceptance.py`, `recompute_supplementary.py`, `recompute_forcing.py` | Corrected-map acceptance and contraction for the main, supplementary and forcing ensembles, with the unchanged `uncertainty.py` | corrected tables (in `../results_v11/rescoring/`) at their working paths |
| `identifiability_corrected.py` | Identifiability analysis on the corrected tables | as above |
| `manuscript_numbers.py` | Wrote `../results_v11/rescoring/manuscript_numbers.json` with the project's own functions (a record, not an independent check) | as above |
| `drainage_rescore.py` | Rescored the drainage-diagnosis states | stored drainage states |
| `evening_holdout_corrected.py`, `plot_holdout_corrected.py` | Evening-pass check of the radar map, before and after the correction | evening scene, model grid |
| `make_figures_corrected.py` | Figures of the corrected-map results | model directory, rerun states |

`sar_georef_correction.py` here is the later local version (commit `1ef8945`), which changes only the sign of
the reported displacement. The maps were written by the earlier version (`3d1ce67`, SHA-256 `8b07b5de…`,
recorded in `../results_v11/radar/georef_qc.json`). The shifts and maps are the same.
