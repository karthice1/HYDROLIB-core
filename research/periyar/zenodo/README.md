# Periyar 2018 flood: open-data D-Flow FM model, sensitivity and uncertainty workflow

Research code, per-run results and manuscript sources for:

> Karthik V and Roopa D. *How far can one gauge mark and one radar image narrow parameter
> uncertainty in an open-data flood model? Sensitivity analysis of the 2018 Periyar flood.*
> Submitted to Environmental Modelling & Software.

Authors: Dr Karthik V (corresponding, karthik.vasanthakumar@kahedu.edu.in) and Dr Roopa D,
Karpagam Academy of Higher Education, Coimbatore, India.

## What is in this archive

This archive is a snapshot of the GitHub repository
<https://github.com/karthice1/HYDROLIB-core> (tag `v1.0-periyar`; the full commit hash is in
`COMMIT.txt`). It is a fork of Deltares HYDROLIB-core 1.5.0.

| Folder | Contents |
|---|---|
| `hydrolib/` | HYDROLIB-core 1.5.0 by Deltares, unchanged (MIT licence) |
| `research/gsa/` | General sensitivity-analysis toolkit for any D-Flow FM model (Morris, Sobol, LHS), with tests |
| `research/periyar/` | The study: model build, ensembles, radar flood map, validation, uncertainty, drainage diagnosis, figures |
| `research/periyar/data/` | Discharges and dam releases digitised from the Central Water Commission 2018 reports |
| `research/periyar/results/` | Per-run results (CSV) of every ensemble in the paper |
| `research/periyar/manuscript/` | LaTeX source of the submitted manuscript |

The large open input data sets (Copernicus DEM, canopy height, JRC water, ESA WorldCover,
Sentinel-1 scenes) are not included. `research/periyar/data_manifest.csv` gives the exact
product, tile or scene, download location and checksum of each file.

## Quick start

1. Check the published numbers without running the model (under a minute):

   ```bash
   pip install -e . && pip install -r research/periyar/requirements.txt
   python -m research.periyar.check_results
   ```

   It recomputes 120 published values from the stored results and prints PASS or FAIL for each.

   `python -m research.periyar.verify_data DATA_DIR` checks downloaded source data against
   `research/periyar/data_manifest.csv` (product, tile or scene, URL, processing, SHA-256).
   `research/periyar/test_report/` holds the clean-environment test report of the same commit
   as `COMMIT.txt`.
   These checks do not run D-Flow FM.

2. Rerun or reuse the workflow: see `research/periyar/REPRODUCE.md`. Model runs need
   D-Flow FM 1.2.184 (Deltares Delft3D source, commit
   `9c5ef6b8b0bddb3a2fe3b987338c34757f30b38a`).

## Licence

See `NOTICE`. HYDROLIB-core 1.5.0 (everything outside `research/`) is copyright Deltares,
MIT licence (`LICENSE`), included unchanged. The research code in `research/gsa` and
`research/periyar` is copyright Karthik V and Roopa D, MIT licence (`research/LICENSE`).
The manuscript source in `research/periyar/manuscript` is not covered by the MIT licence.
Third-party source data are not redistributed; their licences are listed in
`research/periyar/data_manifest.csv`.
