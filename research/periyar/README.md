# Periyar 2018 flood model and sensitivity analysis

Scripts, inputs and results behind the manuscript *Global sensitivity analysis of
the 2018 Periyar flood with open data* (`manuscript/manuscript.tex`).

## Inputs

| Input | Source | Use |
|---|---|---|
| Terrain | Copernicus DEM GLO-30 (`s3://copernicus-dem-30m`) | Bed level (a surface model: trees and buildings are included) |
| Canopy height | Meta/WRI global canopy height, 1 m (`s3://dataforgood-fb-data`) | DSM correction |
| Permanent water | JRC Global Surface Water occurrence | Channel cells, lowered by an assumed depth (no public bathymetry) |
| Land cover | ESA WorldCover 2021 | Manning roughness per class |
| Discharges, dam releases, HFL | CWC flood report, digitised (`data/README.md`) | Boundaries and validation |
| Sentinel-1A GRD, 21 Aug and 28 Jul 2018 | Copernicus, AWS open data | Validation flood map |

`DATA_DIR` holds `dem/`, `chm/`, `water/`, `landcover/` and `s1/` (the two SAFE
folders, `20180821` and `20180728`).

## Software

D-Flow FM computational kernel 1.2.184, built from the Deltares Delft3D source at
commit `9c5ef6b8b0bddb3a2fe3b987338c34757f30b38a`. Python 3.13 with HYDROLIB-core
1.5.0, meshkernel 8.3.0, SALib 1.6.0, numpy 2.5.3, pandas 3.0.5, rasterio 1.5.2,
netCDF4 1.7.5 and matplotlib 3.11.2. `DFLOWFM` below is the D-Flow FM executable or a
wrapper script that sets its library path.

## Reproducing the results

```bash
# Model v5 (200 m) and the radar flood map
python -m research.periyar.build_base_model DATA_DIR MODEL_DIR --cell-size 200 --channel-depth 4
python -m research.periyar.sar_flood_map MODEL_DIR DATA_DIR/s1/20180821 DATA_DIR/s1/20180728 SAR.tif

# Morris screening, 20 trajectories (180 runs), and Morris indices
python -m research.periyar.sensitivity MODEL_DIR RUNS_MORRIS --method morris --trajectories 20 \
    --parallel 4 --threads 1 --dflowfm DFLOWFM --sar SAR.tif
python -m research.periyar.analyse_runs RUNS_MORRIS

# Conditioning ensembles over the five leading factors (200 runs, seed 42; 400 runs, seed 43)
# (the 400-run ensemble uses --samples 400 --seed 43)
python -m research.periyar.sensitivity MODEL_DIR RUNS_LHS --method lhs --samples 200 --seed 42 \
    --factors n_floodplain channel_depth canopy_factor q_rating q_idukki \
    --parallel 4 --threads 1 --dflowfm DFLOWFM --sar SAR.tif
python -m research.periyar.analyse_runs RUNS_LHS --no-indices
python -m research.periyar.uncertainty OUT_DIR RUNS_LHS [RUNS_LHS2 ...]
python -m research.periyar.range_stability OUT_DIR/range_stability.csv OUT_DIR/runs.csv

# Behavioural runs again: flooded land per hour, confusion matrices, cell states
python -m research.periyar.time_series MODEL_DIR RUNS_LHS OUT.csv --dflowfm DFLOWFM --sar SAR.tif

# Structural tests at nominal parameters (terrain, canopy, channel mask, grid, tide)
python -m research.periyar.structural DATA_DIR SAR.tif OUT_DIR --dflowfm DFLOWFM

# Scores against radar maps made with other thresholds
python -m research.periyar.radar_check OUT.csv --states NAME=radar_state.npz ... --sar NAME=map.tif ...

# Figures
python -m research.periyar.figures MODEL_DIR SAR.tif FIG_DIR --ensemble RUNS_LHS [RUNS_LHS2 ...] --time-series OUT.csv \
    --uncertainty OUT_DIR
```

The digitised CWC series are in `data/`; per-run results of every ensemble and the
revision tests are in `results/` (`results/revision` for the tests added after review).
Map files are deleted after each run is scored, so the per-run CSV and JSON files and
the `radar_state.npz` cell states are the stored outputs.

**Known limitations.** The sea level is a synthetic M2 + K1 tide, no rainfall falls on
the model domain, and gauge coordinates are inferred from the water mask. See the
manuscript, Section 5.
