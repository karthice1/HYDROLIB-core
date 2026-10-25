"""Figures for the Periyar 2018 manuscript.

fig1: study area (terrain, channel, boundaries), the Sentinel-1 flood map, and the
nominal model's flooded land at the radar pass.
fig2: prior and conditioned ranges of peak flooded land from the ensembles.

Usage:
    python -m research.periyar.figures MODEL_DIR SAR_TIF OUT_DIR --ensemble RUNS_DIR [RUNS_DIR ...]
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import pandas as pd
import rasterio

from hydrolib.core.dflowfm.polyfile.models import PolyFile

from .validate import SAR_TIME, WET_DEPTH, _sample
from .build_base_model import REFERENCE


def _extent(profile):
    t = profile["transform"]
    return [t.c / 1e3, (t.c + t.a * profile["width"]) / 1e3,
            (t.f + t.e * profile["height"]) / 1e3, t.f / 1e3]


def _polyline(path):
    pts = PolyFile(path).objects[0].points
    return np.array([p.x for p in pts]) / 1e3, np.array([p.y for p in pts]) / 1e3


def fig1(model_dir: Path, sar_file: Path, out: Path) -> None:
    with rasterio.open(model_dir / "dem_model.tif") as src:
        dem, profile = src.read(1), src.profile
    with rasterio.open(model_dir / "channel.tif") as src:
        channel = src.read(1) > 0.5
    with rasterio.open(sar_file) as src:
        sar = src.read(1)
    map_file = model_dir / "DFM_OUTPUT_periyar" / "periyar_map.nc"
    with netCDF4.Dataset(map_file) as nc:
        x, y = nc["mesh2d_face_x"][:], nc["mesh2d_face_y"][:]
        hours = nc["time"][:] / 3600.0
        step = int(np.argmin(np.abs(hours - (pd.Timestamp(SAR_TIME) - REFERENCE)
                                    .total_seconds() / 3600)))
        depth0 = nc["mesh2d_waterdepth"][0, :]
        depth = nc["mesh2d_waterdepth"][step, :]
    extent = _extent(profile)
    land_wet = (depth > WET_DEPTH) & (depth0 <= WET_DEPTH)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5.2), sharey=True, constrained_layout=True)
    ax = axes[0]
    im = ax.imshow(np.ma.masked_less(dem, -100), extent=extent, cmap="terrain", vmin=-5, vmax=40)
    ax.imshow(np.ma.masked_where(~channel, channel), extent=extent, cmap="Blues", vmin=0, vmax=1.5)
    for name, colour in (("upstream.pli", "red"), ("sea.pli", "black")):
        px, py = _polyline(model_dir / name)
        ax.plot(px, py, color=colour, lw=2.5, label=name.split(".")[0])
    ax.legend(loc="lower right", fontsize=8)
    fig.colorbar(im, ax=ax, shrink=0.7, label="Canopy-corrected terrain (m)")
    ax.set_title("(a) Model terrain, channel and boundaries")

    ax = axes[1]
    ax.imshow(np.ma.masked_less(dem, -100), extent=extent, cmap="Greys", vmin=-5, vmax=60, alpha=0.6)
    im = ax.imshow(np.ma.masked_less_equal(sar, 0.05), extent=extent, cmap="Blues", vmin=0, vmax=1)
    fig.colorbar(im, ax=ax, shrink=0.7, label="Flooded fraction of 200 m cell")
    ax.set_title("(b) Sentinel-1 open-water flood, 21 Aug 06:10 IST")

    ax = axes[2]
    ax.imshow(np.ma.masked_less(dem, -100), extent=extent, cmap="Greys", vmin=-5, vmax=60, alpha=0.6)
    observed = _sample(sar_file, x, y, -1.0) > 0.5
    for mask, colour, label in ((land_wet & ~observed, "tab:orange", "model only"),
                                (land_wet & observed, "tab:green", "both"),
                                (~land_wet & observed & (depth0 <= WET_DEPTH), "tab:blue", "radar only")):
        ax.scatter(x[mask] / 1e3, y[mask] / 1e3, s=1.5, c=colour, marker="s", label=label,
                   linewidths=0)
    ax.legend(loc="lower right", fontsize=8, markerscale=6)
    ax.set_title(f"(c) Nominal model v3, hour {hours[step]:.0f}")
    for ax in axes:
        ax.set_xlabel("UTM 43N easting (km)")
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])
    axes[0].set_ylabel("Northing (km)")
    fig.savefig(out, dpi=200)
    plt.close(fig)


def fig2(ensembles, out: Path) -> None:
    runs = pd.concat([pd.read_csv(d / "runs_with_behaviour.csv") for d in ensembles])
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    ax = axes[0]
    bins = np.linspace(0, runs["max_flooded_area_km2"].max() * 1.05, 30)
    ax.hist(runs["max_flooded_area_km2"], bins=bins, color="0.75", label=f"all runs (n={len(runs)})")
    b = runs[runs["behavioural"]]
    ax.hist(b["max_flooded_area_km2"], bins=bins, color="tab:blue",
            label=f"match HFL and radar (n={len(b)})")
    ax.set_xlabel("Peak flooded land (km²)")
    ax.set_ylabel("Runs")
    ax.legend(fontsize=8)
    ax.set_title("(a) Prior and conditioned peak flooded land")
    ax = axes[1]
    ax.scatter(runs["channel_depth"], runs["n_floodplain"], s=12, c="0.75", label="all runs")
    ax.scatter(b["channel_depth"], b["n_floodplain"], s=24, c="tab:blue", label="behavioural")
    ax.set_xlabel("Channel depth (m)")
    ax.set_ylabel("Floodplain roughness multiplier")
    ax.legend(fontsize=8)
    ax.set_title("(b) Behavioural parameter combinations")
    fig.savefig(out, dpi=200)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("model_dir", type=Path)
    parser.add_argument("sar_file", type=Path)
    parser.add_argument("out_dir", type=Path)
    parser.add_argument("--ensemble", type=Path, nargs="+", default=[])
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    fig1(args.model_dir, args.sar_file, args.out_dir / "fig1_study_area_validation.png")
    if args.ensemble:
        fig2(args.ensemble, args.out_dir / "fig2_conditioning.png")


if __name__ == "__main__":
    main()
