# Periyar 2018 base model

`build_base_model.py` builds a first 2D D-Flow FM model of the lower Periyar
River (Neeleswaram to the sea) with HYDROLIB-core, from open global data only:

| Input | Source | Use |
|---|---|---|
| Terrain | Copernicus DEM GLO-30 (`s3://copernicus-dem-30m`) | Bed level (a surface model: trees and buildings are included) |
| Permanent water | JRC Global Surface Water occurrence | Channel cells, lowered by an assumed depth (no public bathymetry) |
| Land cover | ESA WorldCover 2021 | Manning roughness per class |

```bash
python -m research.periyar.build_base_model DATA_DIR OUT_DIR --cell-size 200 --channel-depth 4
```

`DATA_DIR` holds `dem/`, `water/` and `landcover/` GeoTIFF tiles. The model is a
rectilinear grid clipped to the connected area below `--max-elevation` (20 m).

**The boundary forcing is a placeholder.** The upstream discharge is a synthetic
flood wave and the sea level is a synthetic M2 tide. Replace them with the CWC
Neeleswaram record, the Idukki and Idamalayar releases and the Kochi tide before
any result is used. Boundary locations are approximate and must also be checked.

The assumed channel depth, the roughness table and the boundary forcing are the
inputs the sensitivity analysis in `research/gsa` will vary.
