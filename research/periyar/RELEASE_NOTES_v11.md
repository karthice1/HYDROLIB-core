# v1.1-periyar release notes (proposed; not yet published)

## Why this version exists

In v1.0, the two descending Sentinel-1 scenes (21 Aug and 28 Jul 2018) were placed on the map with ground-control points at hill heights of 558–755 m. As a result, they sat 840 m from their true position (820 m W, 180 m N, along the look direction). The ascending evening scene was off by 150 m.

v1.1 corrects this with a rigid shift of each scene, estimated against JRC water occurrence. It then rebuilds the flood map with the unchanged v1.0 rules and rescores every run against the corrected map. No threshold, acceptance rule, tolerance, factor range or model was changed.

## Effect on the results

| | v1.0 (uncorrected map) | v1.1 (corrected map) |
|---|---|---|
| accepted at ±0.5 / ±1.0 / ±1.5 / ±2.0 m | 23 / 73 / 139 / 217 | 20 / 69 / 135 / 213 |
| C at ±0.5 / ±1.0 / ±1.5 / ±2.0 m | 0.771 / 0.616 / 0.531 / 0.336 | 0.767 / 0.605 / 0.521 / 0.326 |
| radar-only C (five factors) | 0.198 (492 runs) | 0.188 (488 runs) |
| all-land CSI, accepted runs at ±1.0 m | 0.06–0.08 | 0.116–0.154 |
| open-land FAR, accepted runs at ±1.0 m | 0.90–0.92 | 0.44–0.57 |
| eight-factor C, seed 44 / seed 46 / pooled | 0.73 / 0.69 / 0.72 | 0.725 / 0.711 / 0.734 |

The conditional contraction changes by 0.004–0.011. The spatial agreement on open land changes substantially. The all-land excess of model flooding over the radar map is present against both maps.

## Kept unchanged

- **The v1.0 record:** the v1.0 tag, its Zenodo version (10.5281/zenodo.24246263), every v1.0 result file and `check_results.py`.
- **Historical checks still pass:** `check_results.py` still reproduces all 120 v1.0 values.

## Known limits

- **Eight-factor results:** the radar-only result and the ±1.5/±2.0 m results are not computable from this archive, because 308 eight-factor runs were not rerun and have no state. They are not reported in the manuscript.
- **Rerun differences:** 10 of 884 reruns differ slightly from their stored v1.0 results. The differences are at most 0.205 m in HFL error and do not change any HFL-test decision. The stored results are used throughout.
- **The correction is rigid:** terrain-dependent distortion within the scene is not corrected.
- **No ground reference for the correction:** the shift was estimated against JRC water occurrence and checked against WorldCover built-up land, both global products. The residual check is limited to one 20 m reference pixel.
- **Map file hash:** the distributed corrected map has SHA-256 `cda3ee7a…691d8b`; rerunning `v11_scripts/sar_georef_correction.py` writes a file with SHA-256 `faf9e9ad…72af95`. The pixels and georeferencing are identical; only the compressed byte stream differs (see `REPRODUCE_v11.md`).
- **Provenance scripts:** the scripts in `v11_scripts/` contain absolute paths of the original computing environment and do not run unchanged elsewhere (`v11_scripts/README_PROVENANCE.md`). The portable commands are `check_results_v11.py` and `release_checks_v11.sh`.
