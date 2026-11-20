"""Flood extent from a Sentinel-1 GRD scene, on the model's 200 m raster grid.

The scene's VV band is calibrated to sigma0, median filtered against speckle,
georeferenced from its ground control points and averaged to 20 m. A pixel is
flooded when its post-event backscatter is below a threshold (Otsu on the window,
unless given) and it dropped by more than ``--min-drop`` dB from a pre-event scene of
the same track. Permanent water (the model's channel mask) is removed. The result is
the flooded fraction of each 200 m cell of ``dem_model.tif``.

Usage:
    python -m research.periyar.sar_flood_map MODEL_DIR POST_SAFE PRE_SAFE OUT.tif
where POST_SAFE and PRE_SAFE contain measurement/iw-vv.tiff and
annotation/calibration/calibration-iw-vv.xml.
"""

import argparse
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import rasterio
from rasterio.control import GroundControlPoint
from rasterio.transform import from_origin
from rasterio.warp import Resampling, reproject, transform_bounds
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import median_filter

FINE = 20.0  # m, grid on which the SAR scene is resampled before aggregation


def _calibration(xml_file: Path) -> RegularGridInterpolator:
    """sigmaNought look-up table as an interpolator over (line, pixel)."""
    vectors = ET.parse(xml_file).getroot().iter("calibrationVector")
    lines, table, pixels = [], [], None
    for v in vectors:
        lines.append(int(v.find("line").text))
        pixels = np.array(v.find("pixel").text.split(), dtype=float)
        table.append(np.array(v.find("sigmaNought").text.split(), dtype=float))
    return RegularGridInterpolator((np.array(lines, float), pixels), np.array(table),
                                   bounds_error=False, fill_value=None)


def _radar_window(gcps, bounds_lonlat, margin=200):
    """Row/column window of the scene that covers a lon/lat box (affine fit to GCPs)."""
    g = np.array([(p.x, p.y, p.row, p.col) for p in gcps])
    design = np.c_[g[:, 0], g[:, 1], np.ones(len(g))]
    coef_row = np.linalg.lstsq(design, g[:, 2], rcond=None)[0]
    coef_col = np.linalg.lstsq(design, g[:, 3], rcond=None)[0]
    w, s, e, n = bounds_lonlat
    corners = np.array([(w, s, 1), (w, n, 1), (e, s, 1), (e, n, 1)])
    rows, cols = corners @ coef_row, corners @ coef_col
    return (int(rows.min()) - margin, int(rows.max()) + margin,
            int(cols.min()) - margin, int(cols.max()) + margin)


def sigma0_db(safe: Path, grid_profile: dict, filter_size: int = 5) -> np.ndarray:
    """VV sigma0 in dB, resampled (average in linear units) to a FINE grid over the model."""
    with rasterio.open(safe / "measurement" / "iw-vv.tiff") as src:
        gcps, gcp_crs = src.gcps
        r0, r1, c0, c1 = _radar_window(gcps, _bounds(grid_profile))
        r0, c0 = max(r0, 0), max(c0, 0)
        r1, c1 = min(r1, src.height), min(c1, src.width)
        dn = src.read(1, window=((r0, r1), (c0, c1))).astype("float64")
    rows, cols = np.mgrid[r0:r1, c0:c1]
    a = _calibration(safe / "annotation" / "calibration" / "calibration-iw-vv.xml")(
        np.c_[rows.ravel(), cols.ravel()]).reshape(dn.shape)
    sigma0 = np.where(dn > 0, dn ** 2 / a ** 2, np.nan)
    sigma0 = median_filter(np.nan_to_num(sigma0, nan=0.0), size=filter_size)
    local = [GroundControlPoint(row=p.row - r0, col=p.col - c0, x=p.x, y=p.y, z=p.z)
             for p in gcps]
    left, bottom, right, top = _bounds(grid_profile, native=True)
    width, height = int(round((right - left) / FINE)), int(round((top - bottom) / FINE))
    fine = np.full((height, width), np.nan, dtype="float32")
    reproject(sigma0.astype("float32"), fine, gcps=local, src_crs=gcp_crs, src_nodata=0.0,
              dst_transform=from_origin(left, top, FINE, FINE), dst_crs=grid_profile["crs"],
              dst_nodata=np.nan, resampling=Resampling.average)
    with np.errstate(divide="ignore", invalid="ignore"):
        return 10 * np.log10(fine)


def _bounds(profile, native=False):
    t = profile["transform"]
    left, top = t.c, t.f
    right, bottom = left + t.a * profile["width"], top + t.e * profile["height"]
    if native:
        return left, bottom, right, top
    return transform_bounds(profile["crs"], "EPSG:4326", left, bottom, right, top)


def otsu(values: np.ndarray, bins: int = 256) -> float:
    values = values[np.isfinite(values)]
    hist, edges = np.histogram(values, bins=bins)
    centres = (edges[:-1] + edges[1:]) / 2
    w0 = np.cumsum(hist)
    w1 = w0[-1] - w0
    m0 = np.cumsum(hist * centres) / np.maximum(w0, 1)
    m1 = (np.sum(hist * centres) - np.cumsum(hist * centres)) / np.maximum(w1, 1)
    return float(centres[np.argmax(w0 * w1 * (m0 - m1) ** 2)])


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("model_dir", type=Path)
    parser.add_argument("post", type=Path)
    parser.add_argument("pre", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--threshold", type=float, help="dB; Otsu on the post scene if omitted")
    parser.add_argument("--min-drop", type=float, default=3.0, help="dB drop from pre to post")
    args = parser.parse_args()

    with rasterio.open(args.model_dir / "dem_model.tif") as src:
        profile = src.profile
    with rasterio.open(args.model_dir / "channel.tif") as src:
        channel = src.read(1) > 0.5

    post = sigma0_db(args.post, profile)
    pre = sigma0_db(args.pre, profile)
    threshold = args.threshold if args.threshold is not None else otsu(post)
    flooded = (post < threshold) & (pre - post > args.min_drop)
    valid = np.isfinite(post) & np.isfinite(pre)

    factor = int(round(profile["transform"].a / FINE))
    h, w = profile["height"], profile["width"]
    def coarse(a):
        return a[:h * factor, :w * factor].reshape(h, factor, w, factor).mean(axis=(1, 3))
    fraction = coarse(flooded.astype("float32"))
    coverage = coarse(valid.astype("float32"))
    fraction = np.where((coverage > 0.5) & ~channel, fraction, -1.0).astype("float32")

    out_profile = {**profile, "dtype": "float32", "nodata": -1.0, "count": 1}
    with rasterio.open(args.out, "w", **out_profile) as dst:
        dst.write(fraction, 1)
    land = fraction >= 0
    print(f"threshold {threshold:.1f} dB; flooded land (fraction > 0.5): "
          f"{(fraction[land] > 0.5).sum() * profile['transform'].a ** 2 / 1e6:.0f} km2")


if __name__ == "__main__":
    main()
