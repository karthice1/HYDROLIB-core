"""Figures for the Periyar 2018 manuscript.

fig1: study area (terrain, channel, boundaries), the Sentinel-1 flood map, and the
nominal model's flooded land at the radar pass.
fig2: prior and conditioned peak flooded land, and depth against roughness.
fig3: flooded land through time for the behavioural runs (from time_series.py).
fig4: model against radar, cell by cell (fig_confusion).
fig5: prior and conditioned distributions, threshold grid and HFL residuals.
fig6: drainage state of the land flooded at the radar pass (from drainage.py).

Usage:
    python -m research.periyar.figures MODEL_DIR SAR_TIF OUT_DIR
        [--time-series CSV] [--uncertainty DIR] [--strict DIR] [--drainage NPZ NPZ]
"""

import argparse
from pathlib import Path

import matplotlib
import matplotlib.dates

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Large enough to read when a full-width figure is printed at journal-page width.
plt.rcParams.update({"font.size": 15, "axes.titlesize": 15, "axes.labelsize": 15,
                     "legend.fontsize": 13, "xtick.labelsize": 13, "ytick.labelsize": 13})
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

    fig, axes = plt.subplots(1, 3, figsize=(12.5, 5.6), sharey=True, constrained_layout=True)
    ax = axes[0]
    im = ax.imshow(np.ma.masked_less(dem, -100), extent=extent, cmap="terrain", vmin=-5, vmax=40)
    ax.imshow(np.ma.masked_where(~channel, channel), extent=extent, cmap="Blues", vmin=0, vmax=1.5)
    for name, colour in (("upstream.pli", "red"), ("chalakudy.pli", "magenta"), ("sea.pli", "black")):
        if not (model_dir / name).exists():
            continue
        px, py = _polyline(model_dir / name)
        ax.plot(px, py, color=colour, lw=2.5, label=name.split(".")[0])
    ax.legend(loc="lower right", fontsize=12)
    fig.colorbar(im, ax=ax, shrink=0.9, orientation="horizontal", label="Terrain (m)")
    ax.set_title("(a) Terrain and boundaries")

    ax = axes[1]
    ax.imshow(np.ma.masked_less(dem, -100), extent=extent, cmap="Greys", vmin=-5, vmax=60, alpha=0.6)
    im = ax.imshow(np.ma.masked_less_equal(sar, 0.05), extent=extent, cmap="Blues", vmin=0, vmax=1)
    fig.colorbar(im, ax=ax, shrink=0.9, orientation="horizontal", label="Flooded fraction of cell")
    ax.set_title("(b) Sentinel-1, 21 Aug")

    ax = axes[2]
    ax.imshow(np.ma.masked_less(dem, -100), extent=extent, cmap="Greys", vmin=-5, vmax=60, alpha=0.6)
    observed = _sample(sar_file, x, y, -1.0) > 0.5
    for mask, colour, kind in ((land_wet & ~observed, "tab:orange", "model only"),
                                (land_wet & observed, "tab:green", "both"),
                                (~land_wet & observed & (depth0 <= WET_DEPTH), "tab:blue", "radar only")):
        ax.scatter(x[mask] / 1e3, y[mask] / 1e3, s=1.5, c=colour, marker="s", label=kind,
                   linewidths=0)
    ax.legend(loc="lower right", fontsize=12, markerscale=6)
    ax.set_title(f"(c) Nominal{label} at the pass")
    for ax in axes:
        ax.set_xlabel("Easting (km)")
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])
    axes[0].set_ylabel("Northing (km, UTM 43N)")
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)


def fig2(central_runs: Path, strict_runs: Path, out: Path) -> None:
    """Prior and accepted peak flooded land, and depth against roughness, at two tolerances."""
    runs = pd.read_csv(central_runs)
    strict = pd.read_csv(strict_runs)["accepted"].to_numpy()
    central = runs["accepted"].to_numpy()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), constrained_layout=True)
    ax = axes[0]
    bins = np.linspace(0, runs["max_flooded_area_km2"].max() * 1.05, 30)
    ax.hist(runs["max_flooded_area_km2"], bins=bins, color="0.75", label=f"all runs (n={len(runs)})")
    ax.hist(runs.loc[central, "max_flooded_area_km2"], bins=bins, color="#9ec5ec",
            label=f"accepted, ±1.0 m (n={central.sum()})")
    ax.hist(runs.loc[strict, "max_flooded_area_km2"], bins=bins, color="#2a78d6",
            label=f"accepted, ±0.5 m (n={strict.sum()})")
    ax.set_xlabel("Peak flooded land (km²)")
    ax.set_ylabel("Runs")
    ax.legend(fontsize=12)
    ax.set_title("(a) Peak flooded land")
    ax = axes[1]
    ax.scatter(runs["channel_depth"], runs["n_floodplain"], s=12, c="0.8", label="all runs")
    ax.scatter(runs.loc[central, "channel_depth"], runs.loc[central, "n_floodplain"], s=24,
               c="#9ec5ec", label="accepted, ±1.0 m")
    ax.scatter(runs.loc[strict, "channel_depth"], runs.loc[strict, "n_floodplain"], s=28,
               c="#2a78d6", label="accepted, ±0.5 m")
    ax.set_xlabel("Channel depth (m)")
    ax.set_ylabel("Floodplain roughness multiplier")
    ax.legend(fontsize=12, loc="upper right")
    ax.set_title("(b) Accepted parameter combinations")
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
                fontsize=13, xytext=(-4, 0), textcoords="offset points")
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
    ax.legend(handles[0] + more[0], handles[1] + more[1], fontsize=13, loc="upper left")
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
        ax.set_title(f"{title}\nCSI {hits / (hits + misses + false_alarms):.3f}", fontsize=13)
        ax.set_aspect("equal")
        ax.tick_params(labelsize=12)
        ax.set_xlabel("Easting (km, UTM 43N)", fontsize=13)
    np.atleast_1d(axes)[0].set_ylabel("Northing (km)", fontsize=13)
    ax0 = np.atleast_1d(axes)[0]
    x0, y0 = ax0.get_xlim()[0] + 2, ax0.get_ylim()[0] + 2
    ax0.plot([x0, x0 + 10], [y0, y0], "k-", lw=3)
    ax0.text(x0 + 5, y0 + 1, "10 km", ha="center", fontsize=13)
    fig.legend(handles=[Patch(color=colours[i], label=l) for i, l in
                        [(3, "flooded in both"), (1, "model only"), (2, "radar only"), (0, "dry in both")]],
               loc="lower center", ncol=4, fontsize=13, bbox_to_anchor=(0.5, -0.12))
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)


def fig_uncertainty(uncertainty_dir: Path, out: Path, hfl_tol: float = 1.0) -> None:
    """Prior and conditioned peak flooded land, threshold grid and HFL residuals."""
    from .uncertainty import accepted
    runs = pd.read_csv(uncertainty_dir / "runs.csv")
    grid = pd.read_csv(uncertainty_dir / "threshold_grid.csv")
    fig = plt.figure(figsize=(12, 10.5), constrained_layout=True)
    grid_spec = fig.add_gridspec(2, 2)
    axes = [fig.add_subplot(grid_spec[0, 0]), fig.add_subplot(grid_spec[1, :]),
            fig.add_subplot(grid_spec[0, 1])]
    ax = axes[0]
    sets = [("all runs", np.ones(len(runs), bool), "0.55"),
            ("HFL test only", accepted(runs, hfl_tol=hfl_tol, use_radar=False), "#eb6834"),
            ("radar test only", accepted(runs, hfl_tol=hfl_tol, use_hfl=False), "#1baf7a"),
            ("both tests", accepted(runs, hfl_tol=hfl_tol), "#2a78d6")]
    for label, mask, colour in sets:
        v = np.sort(runs.loc[mask, "max_flooded_area_km2"].to_numpy())
        ax.step(v, np.arange(1, len(v) + 1) / len(v), where="post", color=colour, lw=2,
                label=f"{label} (n={len(v)})")
    ax.set_xlabel("Peak flooded land (km²)")
    ax.set_ylabel("Cumulative share of runs")
    ax.legend(fontsize=12, loc="lower right")
    ax.set_title(f"(a) Prior and conditioned, ±{hfl_tol:.1f} m")
    ax = axes[1]
    table = grid.pivot(index="hfl_tolerance_m", columns="hit_rate", values="narrowing")
    counts = grid.pivot(index="hfl_tolerance_m", columns="hit_rate", values="n")
    im = ax.imshow(table.to_numpy(), cmap="Blues", vmin=0, vmax=1, aspect="auto", origin="lower")
    for i in range(table.shape[0]):
        for j in range(table.shape[1]):
            value, n = table.iat[i, j], int(counts.iat[i, j])
            text = f"{value:.2f}\nn={n}" if np.isfinite(value) else f"n={n}"
            ax.text(j, i, text, ha="center", va="center", fontsize=12,
                    color="white" if np.isfinite(value) and value > 0.6 else "black")
    ax.set_xticks(range(table.shape[1]), [f"{c:.2f}" for c in table.columns])
    ax.set_yticks(range(table.shape[0]), [f"{r:.2f}" for r in table.index])
    ax.set_xlabel("Minimum open-land hit rate")
    ax.set_ylabel("HFL tolerance (m)")
    fig.colorbar(im, ax=ax, label="Contraction C", shrink=0.9)
    ax.set_title("(b) Contraction C under other acceptance rules")
    ax = axes[2]
    bins = np.arange(-3.0, 5.25, 0.25)
    ax.hist(runs["hfl_error_m"], bins=bins, color="0.75", label=f"all runs (n={len(runs)})")
    both = accepted(runs, hfl_tol=hfl_tol)
    ax.hist(runs.loc[both, "hfl_error_m"], bins=bins, color="#2a78d6",
            label=f"both tests, ±{hfl_tol:.1f} m (n={int(both.sum())})")
    ax.axvspan(-1.0, 1.0, color="0.93", zorder=0)
    ax.axvspan(-0.5, 0.5, color="0.83", zorder=0)
    ax.set_xlabel("Modelled minus observed HFL (m)")
    ax.set_ylabel("Runs")
    ax.legend(fontsize=12)
    ax.set_title("(c) HFL residuals")
    fig.savefig(out, dpi=300)
    plt.close(fig)


def fig_drainage(states, out: Path) -> None:
    """Drainage class of each land cell flooded at the radar pass (drainage.py states)."""
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch
    names = ["closed depression", "draining", "connected, slow"]
    colours = ["#d9d9d9", "#b2182b", "#2166ac", "#f4a582"]  # dry, then the three classes
    fig, axes = plt.subplots(1, len(states), figsize=(5.6 * len(states), 6.0),
                             constrained_layout=True, sharex=True, sharey=True)
    for ax, (title, npz) in zip(np.atleast_1d(axes), states.items()):
        s = np.load(npz)
        code = np.zeros(s["x"].shape, int)
        for i, name in enumerate(names, start=1):
            code[s["wet"] & (s["classes"] == name)] = i
        land = s["land"]
        ax.scatter(s["x"][land] / 1e3, s["y"][land] / 1e3, c=code[land], s=1.2, marker="s",
                   cmap=ListedColormap(colours), vmin=0, vmax=3, linewidths=0)
        area = {n: (code[land] == i).sum() * 0.04 for i, n in enumerate(names, start=1)}
        ax.set_title(f"{title}\nflooded {sum(area.values()):.0f} km², ponded {area[names[0]]:.0f} km²",
                     fontsize=13)
        ax.set_aspect("equal")
        ax.tick_params(labelsize=12)
        ax.set_xlabel("Easting (km, UTM 43N)", fontsize=13)
    np.atleast_1d(axes)[0].set_ylabel("Northing (km)", fontsize=13)
    fig.legend(handles=[Patch(color=colours[i], label=l) for i, l in
                        [(1, "closed depression"), (2, "draining"), (3, "connected, slow"), (0, "dry land")]],
               loc="lower center", ncol=4, fontsize=13, bbox_to_anchor=(0.5, -0.12))
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("model_dir", type=Path)
    parser.add_argument("sar_file", type=Path)
    parser.add_argument("out_dir", type=Path)
    parser.add_argument("--label", default="", help="Model version shown in panel (c)")
    parser.add_argument("--time-series", type=Path, help="Output CSV of time_series.py")
    parser.add_argument("--uncertainty", type=Path, help="Output directory of uncertainty.py (central tolerance)")
    parser.add_argument("--strict", type=Path, help="Output directory of uncertainty.py (strict tolerance)")
    parser.add_argument("--drainage", type=Path, nargs=2, metavar=("NOMINAL_NPZ", "FILLED_NPZ"),
                        help="drainage_state.npz of the nominal and the filled run")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    fig1(args.model_dir, args.sar_file, args.out_dir / "fig1_study_area_validation.png", args.label)
    if args.uncertainty and args.strict:
        fig2(args.uncertainty / "runs.csv", args.strict / "runs.csv", args.out_dir / "fig2_conditioning.png")
    if args.time_series:
        fig3(args.time_series, args.out_dir / "fig3_flooded_land_time_series.png")
    if args.uncertainty:
        fig_uncertainty(args.uncertainty, args.out_dir / "fig5_uncertainty.png")
    if args.drainage:
        fig_drainage({"(a) Nominal model": args.drainage[0],
                      "(b) Depressions filled": args.drainage[1]}, args.out_dir / "fig6_drainage.png")


if __name__ == "__main__":
    main()
