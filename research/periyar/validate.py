"""Compare a Periyar model run with the observed 2018 flood.

Two checks:
* flood extent against a Sentinel-1 flood map from sar_flood_map.py, as a
  contingency table on land cells (cells dry at the start of the run);
* the peak water level near Neeleswaram against the 2018 HFL of EL 12.4 m (CWC).

Usage:
    python -m research.periyar.validate MAP_NC SAR_TIF --sar-time "2018-08-21 06:10"
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
WET_DEPTH = 0.1  # m, as in sensitivity.py


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


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("map_file", type=Path)
    parser.add_argument("sar_file", type=Path)
    parser.add_argument("--sar-time", default="2018-08-21 06:10", help="Acquisition, IST")
    parser.add_argument("--sar-fraction", type=float, default=0.5,
                        help="Flooded fraction above which a SAR cell counts as wet")
    parser.add_argument("--out", type=Path, help="Write the scores as JSON")
    args = parser.parse_args()

    with netCDF4.Dataset(args.map_file) as nc:
        x, y = nc["mesh2d_face_x"][:], nc["mesh2d_face_y"][:]
        hours = nc["time"][:] / 3600.0
        depth = nc["mesh2d_waterdepth"][:]
        level = nc["mesh2d_s1"][:]

    sar_hour = (pd.Timestamp(args.sar_time) - REFERENCE).total_seconds() / 3600
    step = int(np.argmin(np.abs(hours - sar_hour)))
    with rasterio.open(args.sar_file) as src:
        rows, cols = rasterio.transform.rowcol(src.transform, x, y)
        fraction = src.read(1)[np.clip(rows, 0, src.height - 1), np.clip(cols, 0, src.width - 1)]
    land = (depth[0] <= WET_DEPTH) & (fraction >= 0)
    model_wet = depth[step] > WET_DEPTH
    observed_wet = fraction > args.sar_fraction
    scores = extent_scores(model_wet[land], observed_wet[land])
    cell_km2 = float(np.median(np.diff(np.unique(np.round(x)))) ** 2 / 1e6)
    scores.update({
        "model_hour": float(hours[step]), "sar_hour": sar_hour,
        "model_flooded_km2": float(model_wet[land].sum() * cell_km2),
        "sar_flooded_km2": float(observed_wet[land].sum() * cell_km2),
    })

    to_utm = Transformer.from_crs("EPSG:4326", CRS, always_xy=True)
    ux, uy = to_utm.transform(*np.mean(UPSTREAM, axis=0))
    near = np.hypot(x - ux, y - uy) < 500
    peak = level[:, near].max(axis=1)
    scores.update({
        "neeleswaram_peak_level_m": float(peak.max()),
        "neeleswaram_peak_hour": float(hours[int(peak.argmax())]),
        "hfl_observed_m": HFL_NEELESWARAM,
        "hfl_error_m": float(peak.max() - HFL_NEELESWARAM),
    })
    print(json.dumps(scores, indent=2))
    if args.out:
        args.out.write_text(json.dumps(scores, indent=2))


if __name__ == "__main__":
    main()
