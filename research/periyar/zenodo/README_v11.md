# Periyar 2018 flood: open-data D-Flow FM model, sensitivity and uncertainty workflow (v1.1-periyar)

Research code, per-run results and manuscript sources for:

> Karthik V and Roopa D. *How far can one gauge mark and one radar image narrow parameter
> uncertainty in an open-data flood model? Sensitivity analysis of the 2018 Periyar flood.*
> Manuscript prepared for Environmental Modeling & Assessment.

Authors: Dr Karthik V (corresponding, karthik.vasanthakumar@kahedu.edu.in) and Dr Roopa D,
Karpagam Academy of Higher Education, Coimbatore, India.

## Versions

- **v1.1-periyar (this archive)** corrects an 840 m geolocation error of the Sentinel-1 radar
  flood map used in v1.0. It rebuilds the map with the unchanged v1.0 rules and rescores every run
  against the corrected map. No threshold, acceptance rule, tolerance, factor range, seed or model
  was changed. These are the results of the revised manuscript.
- **v1.0-periyar** (Zenodo DOI 10.5281/zenodo.24246263) is the record before the correction. Its
  code and results are still contained here, unchanged; only `README.md`, `data_manifest.csv` (rows
  added) and `.zenodo.json` were edited.

## How to cite

- This version (v1.1-periyar, the corrected analysis): https://doi.org/10.5281/zenodo.24390466
- All versions (concept DOI, resolves to the latest version): https://doi.org/10.5281/zenodo.24246262
- v1.0-periyar, the record before the correction only: https://doi.org/10.5281/zenodo.24246263

See `research/periyar/RELEASE_NOTES_v11.md` for what changed and `research/periyar/REPRODUCE_v11.md`
for how to check it.

## What is in this archive

A snapshot of the GitHub repository <https://github.com/karthice1/HYDROLIB-core> at tag
`v1.1-periyar` (the full commit hash is in `COMMIT.txt`), a fork of Deltares HYDROLIB-core 1.5.0.

| Folder | Contents |
|---|---|
| `hydrolib/` | HYDROLIB-core 1.5.0 by Deltares, unchanged (MIT licence) |
| `research/gsa/` | General sensitivity-analysis toolkit for any D-Flow FM model (Morris, Sobol, LHS), with tests |
| `research/periyar/` | The study (v1.0 code and results, unchanged) |
| `research/periyar/results/` | v1.0 per-run results (computed with the uncorrected radar map) |
| `research/periyar/results_v11/` | v1.1: corrected radar maps, per-cell model states, corrected rescoring, evening-pass check, rerun reproduction |
| `research/periyar/v11_scripts/` | v1.1 provenance scripts, kept as run, with absolute paths of the original environment |
| `research/periyar/test_report_v11/` | Clean-environment test report of the commit in `COMMIT.txt` |
| `research/periyar/manuscript/` | LaTeX source of the v1.0 manuscript (uncorrected numbers) |

The large open input data sets are not included. `research/periyar/data_manifest.csv` gives the
product, tile or scene, download location and checksum of each file.

## Quick start

```bash
pip install -e . && pip install -r research/periyar/requirements.txt
python -m research.periyar.check_results_v11   # 658 corrected-map values (v1.1)
python -m research.periyar.check_results       # 120 v1.0 values, uncorrected map (historical)
python -m research.periyar.verify_data DATA_DIR
```

These checks do not run D-Flow FM. Model runs need D-Flow FM 1.2.184 (Deltares Delft3D source,
commit `9c5ef6b8b0bddb3a2fe3b987338c34757f30b38a`).

## Limits stated in the release notes

The model is not independently validated. The radar correction is a rigid shift estimated against
global reference products, with no ground survey. The eight-factor radar-only and ±1.5/±2.0 m results
cannot be computed from this archive (308 runs have no state). 10 of 884 reruns differ slightly from
their stored v1.0 results (at most 0.205 m).

## Licence

See `NOTICE`. HYDROLIB-core 1.5.0 (everything outside `research/`) is copyright Deltares, MIT
licence (`LICENSE`), included unchanged. The research code and data products in `research/gsa` and
`research/periyar`, including `results_v11/`, are copyright Karthik V and Roopa D, MIT licence
(`research/LICENSE`). The manuscript source in `research/periyar/manuscript` is not covered by the
MIT licence. Third-party source data are not redistributed; their licences are listed in
`research/periyar/data_manifest.csv`.
