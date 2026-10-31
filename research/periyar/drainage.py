"""Why does the model stay flooded after the peak? A diagnosis from one map file.

Classifies every land cell that the model floods at the radar pass:
- "closed depression": the water surface lies below the cell's spill level, the
  lowest level at which water could flow to the sea over the model terrain
  (priority-flood fill of the bed-level raster from the sea edges), so the water
  can leave only by infiltration or evaporation, neither of which is modelled;
- "draining": connected to an outlet and the level fell by more than 0.1 m in the
  24 h before the pass;
- "connected, slow": connected but falling by 0.1 m or less in those 24 h.
It also reports the spread of peak water level over the gauge cells used for
the HFL comparison, and the water still stored on land.

Cells next to a channel cell are counted separately, because there the face bed
level that the model interpolates from its nodes mixes channel and floodplain.

Usage:
    python -m research.periyar.drainage RUN_DIR SAR_TIF OUT_DIR
    python -m research.periyar.drainage RUN_DIR SAR_TIF OUT_DIR --fill-from BASE_DIR --dflowfm EXE
    python -m research.periyar.drainage --terrain-scan BASE_DIR OUT_CSV

With --fill-from, RUN_DIR is first rebuilt from BASE_DIR with every closed
floodplain depression filled to its spill level and run (filled_variant), and its
flood metrics and scores are written to RUN_DIR/summary.json.
"""

import argparse
import heapq
import json
import sys
from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd
import rasterio
from pyproj import Transformer
from scipy.ndimage import binary_dilation

from . import validate
from .build_base_model import CRS, REFERENCE, UPSTREAM

SEA_EDGES = ("west", "south")


def spill_levels(bed: np.ndarray, valid: np.ndarray) -> np.ndarray:
    """Priority-flood: lowest level from which each cell can drain to the sea edges."""
    h, w = bed.shape
    filled = np.full(bed.shape, np.inf)
    heap = []
    seeds = np.zeros(bed.shape, bool)
    seeds[:, 0] = True   # west edge
    seeds[-1, :] = True  # south edge
    # Cells next to nodata on those edges are outlets too (the coast is ragged).
    for i, j in zip(*np.nonzero(seeds & valid)):
        filled[i, j] = bed[i, j]
        heapq.heappush(heap, (bed[i, j], i, j))
    while heap:
        level, i, j = heapq.heappop(heap)
        if level > filled[i, j]:
            continue
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            a, b = i + di, j + dj
            if 0 <= a < h and 0 <= b < w and valid[a, b]:
                new = max(level, bed[a, b])
                if new < filled[a, b]:
                    filled[a, b] = new
                    heapq.heappush(heap, (new, a, b))
    return filled


def main():
    if len(sys.argv) == 4 and sys.argv[1] == "--terrain-scan":
        terrain_scan(Path(sys.argv[2]), Path(sys.argv[3]))
        return
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("sar_file", type=Path)
    parser.add_argument("out_dir", type=Path)
    parser.add_argument("--fill-from", type=Path, help="base model to fill and run first")
    parser.add_argument("--dflowfm", default="dflowfm")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    if args.fill_from:
        from .sensitivity import metrics
        map_file = filled_variant(args.fill_from, args.run_dir, args.dflowfm)
        scores = {**metrics(map_file), **validate.scores(args.run_dir, map_file, args.sar_file)}
        (args.run_dir / "summary.json").write_text(json.dumps(scores, indent=2, default=float))
    map_file = args.run_dir / "DFM_OUTPUT_periyar" / "periyar_map.nc"
    with netCDF4.Dataset(map_file) as nc:
        x, y = np.asarray(nc["mesh2d_face_x"][:]), np.asarray(nc["mesh2d_face_y"][:])
        hours = np.asarray(nc["time"][:]) / 3600.0
        depth = np.asarray(nc["mesh2d_waterdepth"][:])
        level = np.asarray(nc["mesh2d_s1"][:])
        area = np.asarray(nc["mesh2d_flowelem_ba"][:])
        bed_face = np.asarray(nc["mesh2d_flowelem_bl"][:])

    with rasterio.open(args.run_dir / "bedlevel.tif") as src:
        bed = src.read(1).astype(float)
        nodata = src.nodata
        transform = src.transform
    valid = bed != nodata if nodata is not None else np.isfinite(bed)
    filled = spill_levels(np.where(valid, bed, np.inf), valid)
    col = ((x - transform.c) / transform.a).astype(int)
    row = ((y - transform.f) / transform.e).astype(int)
    inside = (row >= 0) & (row < bed.shape[0]) & (col >= 0) & (col < bed.shape[1])
    spill = np.full(x.shape, np.nan)
    spill[inside] = filled[row[inside], col[inside]]

    sar_hour = (pd.Timestamp(validate.SAR_TIME) - REFERENCE).total_seconds() / 3600
    step = int(np.argmin(np.abs(hours - sar_hour)))
    before = int(np.argmin(np.abs(hours - (sar_hour - 24))))
    fraction = validate._sample(args.sar_file, x, y, -1.0)
    landcover = validate._sample(args.run_dir / "landcover.tif", x, y, 0)
    channel = validate._sample(args.run_dir / "channel.tif", x, y, 0.0) > 0.5
    with rasterio.open(args.run_dir / "channel.tif") as src:
        near = binary_dilation(src.read(1) > 0.5, iterations=1)
    next_to_channel = np.zeros(x.shape, bool)
    next_to_channel[inside] = near[row[inside], col[inside]]
    land = (depth[0] <= validate.WET_DEPTH) & (fraction >= 0)
    wet = land & (depth[step] > validate.WET_DEPTH)
    radar = fraction > 0.5
    trapped = wet & (level[step] <= spill + 0.01) & (spill - bed_face > validate.WET_DEPTH)
    falling = (level[before] - level[step]) > 0.1
    classes = np.where(trapped, "closed depression",
                       np.where(falling, "draining", "connected, slow"))
    # Sensitivity of the classes to their two thresholds, with and without the
    # cells next to a channel cell.
    scan = []
    for depth_min in (0.05, 0.1, 0.2, 0.5):
        for fall_min in (0.05, 0.1, 0.2):
            closed = wet & (level[step] <= spill + 0.01) & (spill - bed_face > depth_min)
            drain = wet & ~closed & ((level[before] - level[step]) > fall_min)
            for label, keep in (("all", wet), ("excluding next to channel", wet & ~next_to_channel)):
                scan.append({"depression_depth_min_m": depth_min, "level_fall_min_m": fall_min,
                             "cells": label,
                             "closed_km2": area[closed & keep].sum() / 1e6,
                             "draining_km2": area[drain & keep].sum() / 1e6,
                             "slow_km2": area[keep & ~closed & ~drain].sum() / 1e6,
                             "closed_stored_mcm": (depth[step] * area)[closed & keep].sum() / 1e6})
    pd.DataFrame(scan).to_csv(args.out_dir / "drainage_thresholds.csv", index=False)

    rows = []
    for name in ("closed depression", "draining", "connected, slow"):
        m = wet & (classes == name)
        rows.append({"class": name, "flooded_km2": area[m].sum() / 1e6,
                     "false_alarm_km2": area[m & ~radar].sum() / 1e6,
                     "open_land_false_alarm_km2": area[m & ~radar & np.isin(landcover, validate.OPEN_CLASSES)].sum() / 1e6,
                     "next_to_channel_km2": area[m & next_to_channel].sum() / 1e6,
                     "stored_mcm": (depth[step] * area)[m].sum() / 1e6,
                     "mean_depth_m": float(depth[step][m].mean()) if m.any() else 0.0})
    table = pd.DataFrame(rows)
    table.to_csv(args.out_dir / "drainage_classes.csv", index=False)

    # The ponded cells: how deep their depressions are and what covers the land.
    m = wet & (classes == "closed depression")
    names = {10: "tree cover", 20: "shrubland", 30: "grassland", 40: "cropland", 50: "built-up",
             60: "bare", 80: "water", 90: "herbaceous wetland", 95: "mangroves"}
    cover = pd.Series(area[m], index=[names.get(int(c), "other") for c in landcover[m]])
    ponded = {
        "depression_depth_m_quantiles": dict(zip(["q10", "q50", "q90"], np.round(
            np.quantile((spill - bed_face)[m], [0.1, 0.5, 0.9]), 2).tolist())),
        "water_below_spill_m_quantiles": dict(zip(["q10", "q50", "q90"], np.round(
            np.quantile((spill - level[step])[m], [0.1, 0.5, 0.9]), 2).tolist())),
        "land_cover_share": (cover.groupby(level=0).sum() / cover.sum()).round(3).to_dict(),
        "already_ponded_24h_before_km2": float(area[m & (level[before] <= spill + 0.01)].sum() / 1e6),
    }

    # Spread of the peak level over the gauge cells.
    to_utm = Transformer.from_crs("EPSG:4326", CRS, always_xy=True)
    ux, uy = to_utm.transform(*np.mean(UPSTREAM, axis=0))
    gauge = channel & (np.hypot(x - ux, y - uy) < validate.GAUGE_RADIUS)
    peaks = level[:, gauge].max(axis=0)
    median_series = np.median(level[:, gauge], axis=1)
    summary = {
        "radar_hour": float(hours[step]),
        "flooded_land_at_pass_km2": float(area[wet].sum() / 1e6),
        "false_alarm_km2": float(area[wet & ~radar].sum() / 1e6),
        "land_storage_at_pass_mcm": float((depth[step] * area)[land].sum() / 1e6),
        "land_storage_at_end_mcm": float((depth[-1] * area)[land].sum() / 1e6),
        "gauge_cells": int(gauge.sum()),
        "gauge_cell_peak_levels_m": sorted(np.round(peaks, 2).tolist()),
        "gauge_peak_level_spread_m": float(peaks.max() - peaks.min()),
        "gauge_median_peak_m": float(median_series.max()),
        "gauge_median_peak_hour": float(hours[int(median_series.argmax())]),
        "ponded": ponded,
    }
    (args.out_dir / "drainage_summary.json").write_text(json.dumps(summary, indent=2))
    np.savez_compressed(args.out_dir / "drainage_state.npz", x=x, y=y, land=land, wet=wet,
                        radar=radar, classes=classes.astype("U20"), landcover=landcover,
                        depression_depth=spill - bed_face)
    print(table.round(2).to_string(index=False))
    print(json.dumps(summary, indent=2))


def terrain_scan(base_dir: Path, out_csv: Path) -> pd.DataFrame:
    """Closed-depression storage of the model terrain for several terrain choices.

    Builds the floodplain bed as sensitivity.py does (mean DSM minus a fraction of
    the canopy height; channel cells at the cell-minimum DSM minus 4 m) and, for
    the bare-earth proxy, from the cell-minimum DSM. Depression depth is the
    priority-flood spill level minus the bed, over cells below 20 m.
    """
    def read(name):
        with rasterio.open(base_dir / name) as src:
            return src.read(1).astype(float)
    dsm, dsm_min, canopy = read("dem_raw.tif"), read("dem_min.tif"), read("chm.tif")
    channel = read("channel.tif") > 0.5
    valid = dsm > -9000
    canopy = np.where((canopy >= 0) & (canopy < 100), canopy, 0.0)
    cell_km2 = 0.04
    rows = []
    cases = [(f"canopy factor {c}", np.where(channel, dsm_min - 4.0, dsm - c * canopy))
             for c in (0.3, 0.6, 0.9, 1.0)]
    cases.append(("bare-earth proxy", np.where(channel, dsm_min - 4.0, dsm_min)))
    for label, bed in cases:
        filled = spill_levels(np.where(valid, bed, np.inf), valid)
        depth = np.where(valid & (bed < 20.0) & ~channel, filled - bed, 0.0)
        rows.append({"terrain": label,
                     "depression_area_km2_over_0.1m": float((depth > 0.1).sum() * cell_km2),
                     "depression_area_km2_over_1m": float((depth > 1.0).sum() * cell_km2),
                     "depression_storage_mcm": float(depth.sum() * cell_km2)})
    table = pd.DataFrame(rows)
    table.to_csv(out_csv, index=False)
    print(table.round(1).to_string(index=False))
    return table


def filled_variant(base_dir: Path, run_dir: Path, dflowfm: str) -> Path:
    """Diagnostic run: the nominal model with every closed floodplain depression filled
    to its spill level, so that all land drains by gravity. Not a physical model of
    the floodplain (real depressions drain through canals and sluices), but a test of
    whether the ponding explains the excess flooding at the radar pass."""
    import shutil
    import subprocess
    if run_dir.exists():
        shutil.rmtree(run_dir)
    shutil.copytree(base_dir, run_dir, ignore=shutil.ignore_patterns("DFM_OUTPUT_*", "*.log"))
    with rasterio.open(base_dir / "bedlevel.tif") as src:
        bed, profile = src.read(1).astype(float), src.profile
    with rasterio.open(base_dir / "channel.tif") as src:
        channel = src.read(1) > 0.5
    valid = bed != profile["nodata"]
    filled = spill_levels(np.where(valid, bed, np.inf), valid)
    new = np.where(valid & ~channel, filled, bed)
    with rasterio.open(run_dir / "bedlevel.tif", "w", **profile) as dst:
        dst.write(np.where(valid, new, profile["nodata"]).astype("float32"), 1)
    with (run_dir / "run.log").open("w") as out:
        subprocess.run([dflowfm, "--autostartstop", "periyar.mdu"], cwd=run_dir, check=True,
                       stdout=out, stderr=subprocess.STDOUT, env={**__import__("os").environ,
                                                                  "OMP_NUM_THREADS": "1"})
    return run_dir / "DFM_OUTPUT_periyar" / "periyar_map.nc"


if __name__ == "__main__":
    main()
