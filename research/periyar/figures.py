"""Figures for the Periyar 2018 manuscript.

fig1: study area (terrain, channel, boundaries), the Sentinel-1 flood map, and the
nominal model's flooded land at the radar pass.
fig2: prior and conditioned ranges of peak flooded land from the ensembles.
fig3: flooded land through time for the behavioural runs (from time_series.py).

Usage:
    python -m research.periyar.figures MODEL_DIR SAR_TIF OUT_DIR --ensemble RUNS_DIR [RUNS_DIR ...]
        [--time-series CSV]
"""

import argparse
from pathlib import Path

import matplotlib
import matplotlib.dates

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import pandas as pd
import rasterio

from hydrolib.core.dflowfm.polyfile.models import PolyFile

from .validate import SAR_TIME, WET_DEPTH, _sample
from .build_base_model import CHALAKUDY_DISCHARGE, DEFAULT_DISCHARGE, REFERENCE


def _extent(profile):
    t = profile["transform"]
    return [t.c / 1e3, (t.c + t.a * profile["width"]) / 1e3,
            (t.f + t.e * profile["height"]) / 1e3, t.f / 1e3]


def _polyline(path):
    pts = PolyFile(path).objects[0].points
    return np.array([p.x for p in pts]) / 1e3, np.array([p.y for p in pts]) / 1e3


def fig1(model_dir: Path, sar_file: Path, out: Path, label: str = "") -> None:
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
    for name, colour in (("upstream.pli", "red"), ("chalakudy.pli", "magenta"), ("sea.pli", "black")):
        if not (model_dir / name).exists():
            continue
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
    for mask, colour, kind in ((land_wet & ~observed, "tab:orange", "model only"),
                                (land_wet & observed, "tab:green", "both"),
                                (~land_wet & observed & (depth0 <= WET_DEPTH), "tab:blue", "radar only")):
        ax.scatter(x[mask] / 1e3, y[mask] / 1e3, s=1.5, c=colour, marker="s", label=kind,
                   linewidths=0)
    ax.legend(loc="lower right", fontsize=8, markerscale=6)
    ax.set_title(f"(c) Nominal model{label}, hour {hours[step]:.0f}")
    for ax in axes:
        ax.set_xlabel("UTM 43N easting (km)")
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])
    axes[0].set_ylabel("Northing (km)")
    fig.savefig(out, dpi=300)
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
    fig.savefig(out, dpi=300)
    plt.close(fig)


def fig3(time_series_csv: Path, out: Path, sar_km2: float = 22.0) -> None:
    table = pd.read_csv(time_series_csv, index_col=0, parse_dates=True)
    fig, ax = plt.subplots(figsize=(9, 4.2), constrained_layout=True)
    for i, run in enumerate(table.columns):
        ax.plot(table.index, table[run], color="tab:blue", lw=1.2, alpha=0.8,
                label=f"behavioural runs (n={table.shape[1]})" if i == 0 else None)
    sar_time = pd.Timestamp(SAR_TIME)
    ax.axvline(sar_time, color="k", ls=":", lw=1)
    ax.annotate("Sentinel-1 pass", (sar_time, ax.get_ylim()[1] * 0.95), ha="right",
                fontsize=8, xytext=(-4, 0), textcoords="offset points")
    ax.plot(sar_time, sar_km2, "kD", ms=6, label=f"radar open-water flooding ({sar_km2:.0f} km², lower bound)")
    ax.set_ylabel("Flooded land (km²)")
    ax.set_xlabel("Time (IST), August 2018")
    q = ax.twinx()
    for csv, name, style in [(DEFAULT_DISCHARGE, "Neeleswaram", "-"), (CHALAKUDY_DISCHARGE, "Arangaly (Chalakudy)", "--")]:
        obs = pd.read_csv(csv, parse_dates=["datetime_ist"])
        obs = obs[obs["datetime_ist"].between(table.index[0], table.index[-1])]
        q.plot(obs["datetime_ist"], obs["discharge_m3s"], color="0.5", ls=style, lw=1, label=f"{name} discharge (CWC)")
    q.set_ylabel("Observed discharge (m³/s)", color="0.4")
    q.set_ylim(bottom=0)
    handles = ax.get_legend_handles_labels()
    more = q.get_legend_handles_labels()
    ax.legend(handles[0] + more[0], handles[1] + more[1], fontsize=8, loc="upper left")
    ax.set_xlim(table.index[0], table.index[-1])
    ax.set_ylim(bottom=0)
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%d %b"))
    fig.savefig(out, dpi=300)
    plt.close(fig)


def fig_confusion(states, out: Path, sar_file: Path) -> None:
    """Model against radar at the radar pass, cell by cell, for several runs."""
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch
    from .validate import _sample
    fig, axes = plt.subplots(1, len(states), figsize=(4.6 * len(states), 5.4),
                             constrained_layout=True, sharex=True, sharey=True)
    colours = ["#d9d9d9", "#e66101", "#5e3c99", "#1b7837"]  # dry in both, model only, radar only, both
    for ax, (title, npz) in zip(np.atleast_1d(axes), states.items()):
        s = np.load(npz)
        fraction = _sample(sar_file, s["x"], s["y"], -1.0)
        land = s["land"] & (fraction >= 0)
        model, radar = s["model_wet"], fraction > 0.5
        code = np.where(model & radar, 3, np.where(radar, 2, np.where(model, 1, 0)))
        ax.scatter(s["x"][land] / 1e3, s["y"][land] / 1e3, c=code[land], s=1.2, marker="s",
                   cmap=ListedColormap(colours), vmin=0, vmax=3, linewidths=0)
        hits, misses = int((code[land] == 3).sum()), int((code[land] == 2).sum())
        false_alarms = int((code[land] == 1).sum())
        ax.set_title(f"{title}\nCSI {hits / (hits + misses + false_alarms):.3f}", fontsize=11)
        ax.set_aspect("equal")
        ax.tick_params(labelsize=9)
        ax.set_xlabel("Easting (km, UTM 43N)", fontsize=10)
    np.atleast_1d(axes)[0].set_ylabel("Northing (km)", fontsize=10)
    ax0 = np.atleast_1d(axes)[0]
    x0, y0 = ax0.get_xlim()[0] + 2, ax0.get_ylim()[0] + 2
    ax0.plot([x0, x0 + 10], [y0, y0], "k-", lw=3)
    ax0.text(x0 + 5, y0 + 1, "10 km", ha="center", fontsize=9)
    fig.legend(handles=[Patch(color=colours[i], label=l) for i, l in
                        [(3, "flooded in both"), (1, "model only"), (2, "radar only"), (0, "dry in both")]],
               loc="lower center", ncol=4, fontsize=10, bbox_to_anchor=(0.5, -0.06))
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)


def fig_uncertainty(uncertainty_dir: Path, out: Path) -> None:
    """Prior and conditioned peak flooded land, threshold grid and HFL residuals."""
    from .uncertainty import accepted
    runs = pd.read_csv(uncertainty_dir / "runs.csv")
    grid = pd.read_csv(uncertainty_dir / "threshold_grid.csv")
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4), constrained_layout=True)
    ax = axes[0]
    sets = [("all runs", np.ones(len(runs), bool), "0.55"),
            ("HFL test only", accepted(runs, use_radar=False), "#eb6834"),
            ("radar test only", accepted(runs, use_hfl=False), "#1baf7a"),
            ("both tests", accepted(runs), "#2a78d6")]
    for label, mask, colour in sets:
        v = np.sort(runs.loc[mask, "max_flooded_area_km2"].to_numpy())
        ax.step(v, np.arange(1, len(v) + 1) / len(v), where="post", color=colour, lw=2,
                label=f"{label} (n={len(v)})")
    ax.set_xlabel("Peak flooded land (km²)")
    ax.set_ylabel("Cumulative share of runs")
    ax.legend(fontsize=9, loc="lower right")
    ax.set_title("(a) Prior and conditioned distributions")
    ax = axes[1]
    table = grid.pivot(index="hfl_tolerance_m", columns="hit_rate", values="narrowing")
    counts = grid.pivot(index="hfl_tolerance_m", columns="hit_rate", values="n")
    im = ax.imshow(table.to_numpy(), cmap="Blues", vmin=0, vmax=1, aspect="auto", origin="lower")
    for i in range(table.shape[0]):
        for j in range(table.shape[1]):
            value, n = table.iat[i, j], int(counts.iat[i, j])
            text = f"{value:.2f}\nn={n}" if np.isfinite(value) else f"n={n}"
            ax.text(j, i, text, ha="center", va="center", fontsize=7.5,
                    color="white" if np.isfinite(value) and value > 0.6 else "black")
    ax.set_xticks(range(table.shape[1]), [f"{c:.2f}" for c in table.columns])
    ax.set_yticks(range(table.shape[0]), [f"{r:.2f}" for r in table.index])
    ax.set_xlabel("Minimum open-land hit rate")
    ax.set_ylabel("HFL tolerance (m)")
    fig.colorbar(im, ax=ax, label="Narrowing of peak flooded land", shrink=0.9)
    ax.set_title("(b) Narrowing under other acceptance rules")
    ax = axes[2]
    bins = np.arange(-3.0, 5.25, 0.25)
    ax.hist(runs["hfl_error_m"], bins=bins, color="0.75", label=f"all runs (n={len(runs)})")
    ax.hist(runs.loc[accepted(runs), "hfl_error_m"], bins=bins, color="#2a78d6",
            label=f"both tests (n={int(accepted(runs).sum())})")
    ax.axvspan(-0.5, 0.5, color="0.9", zorder=0)
    ax.set_xlabel("Modelled minus observed HFL (m)")
    ax.set_ylabel("Runs")
    ax.legend(fontsize=9)
    ax.set_title("(c) HFL residuals")
    fig.savefig(out, dpi=300)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("model_dir", type=Path)
    parser.add_argument("sar_file", type=Path)
    parser.add_argument("out_dir", type=Path)
    parser.add_argument("--ensemble", type=Path, nargs="+", default=[],
                        help="Latin hypercube run directories (Morris samples are not random)")
    parser.add_argument("--label", default="", help="Model version shown in panel (c)")
    parser.add_argument("--time-series", type=Path, help="Output CSV of time_series.py")
    parser.add_argument("--uncertainty", type=Path, help="Output directory of uncertainty.py")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    fig1(args.model_dir, args.sar_file, args.out_dir / "fig1_study_area_validation.png", args.label)
    if args.ensemble:
        fig2(args.ensemble, args.out_dir / "fig2_conditioning.png")
    if args.time_series:
        fig3(args.time_series, args.out_dir / "fig3_flooded_land_time_series.png")
    if args.uncertainty:
        fig_uncertainty(args.uncertainty, args.out_dir / "fig5_uncertainty.png")


if __name__ == "__main__":
    main()
