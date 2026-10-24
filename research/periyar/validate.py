"""Compare a Periyar model run with the observed 2018 flood.

Two checks:
* flood extent against a Sentinel-1 flood map from sar_flood_map.py, as a
  contingency table on land cells (cells dry at the start of the run). C-band radar
  sees open water only, and flooded tree cover and built-up land often look bright
  rather than dark, so the scores are also given for the open land-cover classes
  where radar can see floodwater;
* the peak water level in the channel near Neeleswaram against the 2018 HFL of
  EL 12.4 m (CWC).

Usage:
    python -m research.periyar.validate MODEL_DIR SAR_TIF [--sar-time "2018-08-21 06:10"]
"""

import argparse
import json
from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd
import rasterio
from pyproj import Transformer

from .build_base_model import CRS, REFERENCE, UPSTREAM

HFL_NEELESWARAM = 12.4  # m, CWC report "Kerala Floods August 2018"
SAR_TIME = "2018-08-21 06:10"  # IST, Sentinel-1A track 165, 00:40 UTC
WET_DEPTH = 0.1  # m, as in sensitivity.py
OPEN_CLASSES = (30, 40, 60, 90)  # WorldCover grassland, cropland, bare, herbaceous wetland
GAUGE_RADIUS = 500.0  # m


def extent_scores(model_wet: np.ndarray, observed_wet: np.ndarray) -> dict:
    hits = int((model_wet & observed_wet).sum())
    misses = int((~model_wet & observed_wet).sum())
    false_alarms = int((model_wet & ~observed_wet).sum())
    return {
        "hits": hits, "misses": misses, "false_alarms": false_alarms,
        "hit_rate": hits / max(hits + misses, 1),
        "false_alarm_ratio": false_alarms / max(hits + false_alarms, 1),
        "critical_success_index": hits / max(hits + misses + false_alarms, 1),
        "bias": (hits + false_alarms) / max(hits + misses, 1),
    }


def _sample(raster: Path, x, y, fill):
    with rasterio.open(raster) as src:
        rows, cols = rasterio.transform.rowcol(src.transform, x, y)
        rows, cols = np.asarray(rows), np.asarray(cols)
        inside = (rows >= 0) & (rows < src.height) & (cols >= 0) & (cols < src.width)
        data = src.read(1)
    out = np.full(len(rows), fill, dtype=data.dtype)
    out[inside] = data[rows[inside], cols[inside]]
    return out


def scores(model_dir: Path, map_file: Path, sar_file: Path, sar_time: str = SAR_TIME,
           sar_fraction: float = 0.5) -> dict:
    with netCDF4.Dataset(map_file) as nc:
        x, y = nc["mesh2d_face_x"][:], nc["mesh2d_face_y"][:]
        hours = nc["time"][:] / 3600.0
        depth = nc["mesh2d_waterdepth"][:]
        level = nc["mesh2d_s1"][:]
        area = nc["mesh2d_flowelem_ba"][:]

    sar_hour = (pd.Timestamp(sar_time) - REFERENCE).total_seconds() / 3600
    step = int(np.argmin(np.abs(hours - sar_hour)))
    fraction = _sample(sar_file, x, y, -1.0)
    landcover = _sample(model_dir / "landcover.tif", x, y, 0)
    channel = _sample(model_dir / "channel.tif", x, y, 0.0) > 0.5
    land = (depth[0] <= WET_DEPTH) & (fraction >= 0)
    open_land = land & np.isin(landcover, OPEN_CLASSES)
    model_wet = depth[step] > WET_DEPTH
    observed_wet = fraction > sar_fraction
    result = {"model_hour": float(hours[step]), "sar_hour": sar_hour,
              "model_flooded_km2": float(area[land & model_wet].sum() / 1e6),
              "sar_flooded_km2": float(area[land & observed_wet].sum() / 1e6)}
    result.update({f"all_{k}": v for k, v in extent_scores(model_wet[land], observed_wet[land]).items()})
    result.update({f"open_{k}": v for k, v in
                   extent_scores(model_wet[open_land], observed_wet[open_land]).items()})

    to_utm = Transformer.from_crs("EPSG:4326", CRS, always_xy=True)
    ux, uy = to_utm.transform(*np.mean(UPSTREAM, axis=0))
    gauge = channel & (np.hypot(x - ux, y - uy) < GAUGE_RADIUS)
    peak = np.ma.median(level[:, gauge], axis=1)
    result.update({
        "neeleswaram_peak_level_m": float(peak.max()),
        "neeleswaram_peak_hour": float(hours[int(peak.argmax())]),
        "hfl_error_m": float(peak.max() - HFL_NEELESWARAM),
    })
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("model_dir", type=Path)
    parser.add_argument("sar_file", type=Path)
    parser.add_argument("--map-file", type=Path, help="Default: MODEL_DIR/DFM_OUTPUT_periyar/periyar_map.nc")
    parser.add_argument("--sar-time", default=SAR_TIME, help="Acquisition time, IST")
    parser.add_argument("--sar-fraction", type=float, default=0.5,
                        help="Flooded fraction above which a SAR cell counts as wet")
    parser.add_argument("--out", type=Path, help="Write the scores as JSON")
    args = parser.parse_args()
    map_file = args.map_file or args.model_dir / "DFM_OUTPUT_periyar" / "periyar_map.nc"
    result = scores(args.model_dir, map_file, args.sar_file, args.sar_time, args.sar_fraction)
    print(json.dumps(result, indent=2))
    if args.out:
        args.out.write_text(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
