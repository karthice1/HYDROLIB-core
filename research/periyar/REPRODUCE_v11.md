# Reproducing v1.1-periyar (proposed)

v1.1 adds two things to v1.0-periyar:

- the radar geolocation correction;
- the corrected-map rescoring of every run.

Nothing in v1.0 is removed or overwritten. The v1.0 per-run results under `results/` are unchanged. `check_results.py` still reproduces the v1.0 numbers, which were computed with the uncorrected radar map.

## What is new

| Path | Content |
|---|---|
| `results_v11/radar/` | Corrected flood-fraction map (`sar_flood_20180821_corrected.tif`), evening and morning threshold-only maps, five classification variants (`sar_variants/`) and the co-registration record (`georef_qc.json`) |
| `results_v11/model_states/` | Per-cell model state at the radar pass of every rescored run. For each run this is whether the cell is land (dry at the start) and whether it is wet, bit-packed per group, with `grid.npz` holding x, y and WorldCover class. 884 rerun states and 7 states kept from the original study; manifest with SHA-256 in `states_manifest.json` |
| `results_v11/rerun_reproduction.csv` | Every rerun compared with the stored v1.0 result: 874 of 884 exact, 10 differ by at most 0.205 m in HFL error and 0.84 km² in peak land |
| `results_v11/rescoring/` | Corrected per-run tables, outputs of the unchanged `uncertainty.py` on corrected inputs, `manuscript_numbers.json`, drainage rescoring, the median accepted run's drainage diagnosis |
| `results_v11/evening/` | Evening-pass check of the radar map, before and after the correction |
| `check_results_v11.py` | Independent check of the corrected-map values (658 checks). It recomputes every radar score from the states and the map and does not import the analysis modules. |
| `release_checks_v11.sh` | Installs the release in a new environment and runs the unit tests, both checkers and the source-data check; writes `test_report_v11/` |
| `data_manifest.csv` | The v1.0 rows unchanged, plus the annotation files of the two descending scenes, the evening scene (measurement, calibration, annotation) and the two derived corrected maps |
| `v11_scripts/` | The scripts that made the correction, the reruns and the rescoring, with the run lists. They are kept unchanged as provenance and contain the absolute paths of the computing environment in which they ran; `v11_scripts/README_PROVENANCE.md` lists each one and the inputs it needs |

## Portable commands: checks without model runs (minutes)

These use paths relative to the repository and were tested in a new Python 3.13 environment built from the v1.0 tree plus the v1.1 files. The report of the archived commit, `test_report_v11/REPORT.md`, is added to the Zenodo archive by `zenodo/build_archive_v11.sh`; it is not stored in the Git repository, as for v1.0. They need the `results/` folder of v1.0 and the `results_v11/` folder of this release, both in the archive. `DATA_DIR` (the downloaded source files) is needed only for `verify_data`.

```
bash research/periyar/release_checks_v11.sh DATA_DIR      # from the repository root; DATA_DIR optional
```

or individually:

```
python -m research.periyar.check_results_v11      # corrected-map values of the manuscript
python -m research.periyar.check_results          # v1.0 values (uncorrected map), historical record
python -m research.periyar.verify_data DATA_DIR   # source data against data_manifest.csv
```

## Environment-specific reconstruction scripts

The scripts below are provenance. They are kept exactly as they were run and contain absolute paths of the original computing environment (`/home/user/ext/...`, `/home/user/hydrolib-core`, `/mnt/project-files/...`). They do not run unchanged on another computer: edit the path constants at the top of each script first. See `v11_scripts/README_PROVENANCE.md`.

### Rebuilding the corrected map (needs the Sentinel-1 files of the manifest)

```
python -I v11_scripts/sar_georef_correction.py OUT_DIR      # after editing DATA, MODEL and the sys.path line
```

The script reads the three scenes and their annotation files. It estimates the shift of each scene against JRC occurrence. It then rebuilds the map with the frozen v1.0 rules: Otsu threshold, a drop of more than 3 dB, coverage above 0.5 and the channel removed.

The canonical map is the file distributed here, `results_v11/radar/sar_flood_20180821_corrected.tif`, SHA-256 `cda3ee7a00dbe42a5b30ba9a74c3ce121e75bb113f0503580e16cb4542691d8b` (22,533 bytes). Every v1.1 score was computed with it, it is listed in `data_manifest.csv`, and the checker tests it. Running the script writes a file with SHA-256 `faf9e9ad28ca7e277fe27fa6d52e20d3c3d51b711a450276cbf6cefc3872af95` (16,520 bytes). The two files have identical pixels, size, transform, CRS, nodata, data type, DEFLATE compression and strip layout; only the compressed byte stream differs. Compare the pixel arrays, not the file hashes. The same holds for the evening and morning threshold-only maps. `georef_qc.json` gives both hashes for each map (`sha256` for the distributed file, `sha256_as_written_by_script`).

### Rerunning the model (days; needs D-Flow FM 1.2.184 and the v5 model directory)

Each rescored run was rerun with D-Flow FM 1.2.184 from its stored v1.0 parameter set, keeping the map output at the radar pass:

```
python -I v11_scripts/rerun_for_radar.py ENSEMBLE_DIR OUT_DIR MAP RUN_IDS...   # after editing BASE, PUBLISHED and DFM
```

The batches are `run_all.sh` and `run_batch2.sh` to `run_batch4.sh`, with run lists `need_*.txt` and `rest_*.txt`. In total there were 884 reruns:

- 593 main-ensemble runs;
- 92 eight-factor runs, those within ±1.0 m of the HFL;
- 82 fixed-factor runs, those within ±1.0 m of the HFL;
- 111 forcing-variant runs;
- 6 repeats.

The 7 remaining main-ensemble states are those kept in v1.0 (`results/revision/behavioural_radar`).

## What cannot be reproduced from this archive alone

The eight-factor radar-only results and the eight-factor ±1.5 and ±2.0 m results are not reproducible, because 153 + 155 eight-factor runs have no state. The manuscript does not report these values.

Everything that needs D-Flow FM runs, as in v1.0, can only be reproduced by rerunning the model: per-run results, figures, structural tests and the drainage diagnosis.
