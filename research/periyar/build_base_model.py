"""Build a first 2D D-Flow FM base model of the lower Periyar River with HYDROLIB-core.

The model covers the floodplain from Neeleswaram (upstream gauge) to the sea and
uses only open global data:

* Copernicus DEM GLO-30 for terrain (a surface model, so vegetation and buildings
  are still included; a known limitation),
* JRC Global Surface Water occurrence to find the permanent channel, whose bed is
  lowered by an assumed depth because no bathymetry survey is public,
* ESA WorldCover 2021 to assign Manning roughness per land-cover class.

The upstream discharge is the observed Neeleswaram record for 14-20 August 2018,
digitised from the CWC flood report (see data/README.md). The sea level is still a
synthetic M2 + K1 tide and must be replaced with an observed or predicted Kochi tide.
"""

import argparse
from pathlib import Path

import meshkernel as mk
import numpy as np
import pandas as pd
import rasterio
from rasterio import features
from rasterio.warp import Resampling, calculate_default_transform, reproject
from shapely.geometry import shape

from hydrolib.core.dflowfm.bc.models import ForcingModel, QuantityUnitPair, TimeSeries
from hydrolib.core.dflowfm.ext.models import Boundary, ExtModel
from hydrolib.core.dflowfm.inifield.models import IniFieldModel, InitialField, ParameterField
from hydrolib.core.dflowfm.mdu.models import FMModel
from hydrolib.core.dflowfm.net.models import Network
from hydrolib.core.dflowfm.polyfile.models import Point, PolyFile, PolyObject, Metadata

REFERENCE = pd.Timestamp("2018-08-13 00:00")  # model time zero, IST
DURATION_HOURS = 192.0  # 13 Aug 00:00 to 21 Aug 00:00
DEFAULT_DISCHARGE = Path(__file__).parent / "data" / "neeleswaram_discharge_cwc_fig3.csv"
CRS = "EPSG:32643"  # WGS 84 / UTM zone 43N
WINDOW = (76.15, 10.0, 76.55, 10.3)  # lon/lat bounds of the lower Periyar

# Manning n per ESA WorldCover class, mid-range values from Chow (1959) and
# Arcement & Schneider (1989). These are GSA factor baselines, not calibrated.
MANNING = {
    10: 0.10,   # tree cover
    20: 0.07,   # shrubland
    30: 0.035,  # grassland
    40: 0.04,   # cropland
    50: 0.08,   # built-up
    60: 0.025,  # bare / sparse vegetation
    80: 0.03,   # permanent water bodies
    90: 0.05,   # herbaceous wetland
    95: 0.10,   # mangroves
}

# Placeholder boundary locations (lon, lat), to be checked against the gauge
# metadata and satellite imagery.
UPSTREAM = [(76.495, 10.215), (76.505, 10.185)]   # across the river near Neeleswaram
SEA = [(76.155, 10.05), (76.155, 10.29)]          # along the coast, west model edge


def warp(src_paths, dst_path, resolution, resampling, dtype="float32", nodata=-9999.0):
    """Mosaic and reproject source rasters to the model grid in UTM."""
    from rasterio.merge import merge

    sources = [rasterio.open(p) for p in src_paths]
    mosaic, transform = merge(sources, bounds=WINDOW)
    src_crs = sources[0].crs
    height, width = mosaic.shape[1:]
    dst_transform, dst_w, dst_h = calculate_default_transform(
        src_crs, CRS, width, height, *WINDOW, resolution=resolution)
    out = np.full((dst_h, dst_w), nodata, dtype=dtype)
    reproject(mosaic[0], out, src_transform=transform, src_crs=src_crs,
              dst_transform=dst_transform, dst_crs=CRS, resampling=resampling,
              dst_nodata=nodata)
    profile = dict(driver="GTiff", height=dst_h, width=dst_w, count=1, dtype=dtype,
                   crs=CRS, transform=dst_transform, nodata=nodata, compress="deflate")
    with rasterio.open(dst_path, "w", **profile) as dst:
        dst.write(out, 1)
    for s in sources:
        s.close()
    return out, dst_transform


def to_utm(lonlat):
    from pyproj import Transformer

    t = Transformer.from_crs("EPSG:4326", CRS, always_xy=True)
    return [t.transform(lon, lat) for lon, lat in lonlat]


def domain_polygon(dem, transform, max_elevation):
    """Largest connected area below `max_elevation`, as a simplified polygon."""
    mask = (dem > -100) & (dem < max_elevation)
    polygons = [shape(geom) for geom, value in
                features.shapes(mask.astype("uint8"), mask=mask, transform=transform)
                if value == 1]
    largest = max(polygons, key=lambda p: p.area)
    return largest.buffer(0).simplify(transform.a * 2).exterior


def write_polyline(path, name, points):
    pf = PolyFile(objects=[PolyObject(
        metadata=Metadata(name=name, n_rows=len(points), n_columns=2),
        points=[Point(x=x, y=y, data=[]) for x, y in points])])
    pf.save(path)
    return pf


# HYDROLIB-core 1.5.0 always writes these MDU keywords, but recent D-Flow FM
# kernels reject them as obsolete and stop.
OBSOLETE_MDU_KEYWORDS = ("qhrelax", "wrishp_enc")


def _drop_obsolete_keywords(mdu_path: Path) -> None:
    lines = Path(mdu_path).read_text().splitlines(keepends=True)
    keep = [line for line in lines
            if line.split("=")[0].strip().lower() not in OBSOLETE_MDU_KEYWORDS]
    Path(mdu_path).write_text("".join(keep))


def build(data_dir: Path, out_dir: Path, cell_size: float, channel_depth: float,
          max_elevation: float, discharge_csv: Path = DEFAULT_DISCHARGE):
    out_dir.mkdir(parents=True, exist_ok=True)
    dem, transform = warp(sorted((data_dir / "dem").glob("*.tif")), out_dir / "dem_raw.tif",
                          cell_size, Resampling.average)
    water, _ = warp(sorted((data_dir / "water").glob("*.tif")), out_dir / "water.tif",
                    cell_size, Resampling.average)
    landcover, _ = warp(sorted((data_dir / "landcover").glob("*.tif")), out_dir / "landcover.tif",
                        cell_size, Resampling.mode, dtype="int16", nodata=0)

    # Burn an assumed channel depth where water is present most of the time.
    channel = (water > 50) & (water <= 100) & (landcover != 0)
    bed = np.where(channel, dem - channel_depth, dem).astype("float32")
    manning = np.vectorize(lambda c: MANNING.get(int(c), 0.04), otypes=["float32"])(landcover)
    manning[channel] = MANNING[80]
    profile = dict(driver="GTiff", height=bed.shape[0], width=bed.shape[1], count=1,
                   dtype="float32", crs=CRS, transform=transform, nodata=-9999.0,
                   compress="deflate")
    for name, array in (("bedlevel.tif", bed), ("manning.tif", manning)):
        with rasterio.open(out_dir / name, "w", **profile) as dst:
            dst.write(array, 1)

    # 2D mesh: rectilinear grid clipped to the low-lying floodplain.
    ring = domain_polygon(dem, transform, max_elevation)
    xs, ys = map(np.array, ring.xy)
    network = Network()
    network.mesh2d_create_rectilinear_within_extent(ring.bounds, cell_size, cell_size)
    network.mesh2d_clip_mesh(mk.GeometryList(xs, ys), inside=False)

    fm = FMModel()
    fm.filepath = out_dir / "periyar.mdu"
    fm.geometry.netfile.network = network
    fm.geometry.netfile.filepath = "periyar_net.nc"

    inifield = IniFieldModel(
        initial=[InitialField(quantity="bedlevel", datafile="bedlevel.tif",
                              datafiletype="GeoTIFF", interpolationmethod="averaging",
                              averagingtype="mean")],
        parameter=[ParameterField(quantity="frictioncoefficient", datafile="manning.tif",
                                  datafiletype="GeoTIFF", interpolationmethod="averaging",
                                  averagingtype="mean")])
    inifield.filepath = Path("initialfields.ini")
    fm.geometry.inifieldfile = inifield
    fm.physics.uniffricttype = 1  # Manning
    fm.physics.uniffrictcoef = 0.04

    write_polyline(out_dir / "upstream.pli", "upstream", to_utm(UPSTREAM))
    write_polyline(out_dir / "sea.pli", "sea", to_utm(SEA))

    # Upstream: observed Neeleswaram discharge (digitised from the CWC report),
    # held constant before the first and after the last record. Sea: a
    # synthetic M2 + K1 tide until an observed Kochi record is available.
    observed = pd.read_csv(discharge_csv, parse_dates=["datetime_ist"])
    hours_obs = (observed["datetime_ist"] - REFERENCE).dt.total_seconds() / 3600.0
    hours = np.arange(0.0, DURATION_HOURS + 1.0, 1.0)
    q = np.interp(hours, hours_obs, observed["discharge_m3s"])
    h = 0.35 * np.sin(2 * np.pi * hours / 12.42) + 0.2 * np.sin(2 * np.pi * hours / 23.93)
    time_unit = f"hours since {REFERENCE:%Y-%m-%d %H:%M:%S}"
    forcing = ForcingModel(forcing=[
        TimeSeries(name="upstream_0001", function="timeseries", timeinterpolation="linear",
                   quantityunitpair=[QuantityUnitPair(quantity="time", unit=time_unit),
                                     QuantityUnitPair(quantity="dischargebnd", unit="m3/s")],
                   datablock=[[t, v] for t, v in zip(hours, q)]),
        *[TimeSeries(name=f"sea_{i:04d}", function="timeseries", timeinterpolation="linear",
                     quantityunitpair=[QuantityUnitPair(quantity="time", unit=time_unit),
                                       QuantityUnitPair(quantity="waterlevelbnd", unit="m")],
                     datablock=[[t, v] for t, v in zip(hours, h)]) for i in (1, 2)],
    ])
    forcing.filepath = Path("boundaries.bc")
    ext = ExtModel(boundary=[
        Boundary(quantity="dischargebnd", locationfile="upstream.pli", forcingfile=forcing),
        Boundary(quantity="waterlevelbnd", locationfile="sea.pli", forcingfile=forcing),
    ])
    ext.filepath = Path("periyar.ext")
    fm.external_forcing.extforcefilenew = ext

    fm.time.refdate = int(f"{REFERENCE:%Y%m%d}")
    fm.time.tunit = "H"
    fm.time.tstart = 0.0
    fm.time.tstop = DURATION_HOURS
    fm.output.mapinterval = [3600.0]
    fm.output.hisinterval = [600.0]
    fm.save(recurse=True)
    _drop_obsolete_keywords(fm.filepath)

    n_faces = network._mesh2d.get_mesh2d().face_x.size
    print(f"Wrote {fm.filepath} with {n_faces} cells of {cell_size:.0f} m.")
    return fm


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("data_dir", type=Path, help="Folder with dem/, water/ and landcover/ rasters")
    parser.add_argument("out_dir", type=Path)
    parser.add_argument("--cell-size", type=float, default=200.0, help="Grid cell size in metres")
    parser.add_argument("--channel-depth", type=float, default=4.0,
                        help="Assumed channel depth below the DEM water surface in metres")
    parser.add_argument("--max-elevation", type=float, default=20.0,
                        help="Upper terrain elevation of the model domain in metres")
    parser.add_argument("--discharge-csv", type=Path, default=DEFAULT_DISCHARGE,
                        help="Upstream discharge series with datetime_ist and discharge_m3s columns")
    args = parser.parse_args()
    build(args.data_dir, args.out_dir, args.cell_size, args.channel_depth, args.max_elevation,
          args.discharge_csv)
