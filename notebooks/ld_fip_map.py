"""Render one regional FIP map in a fresh Python/GMT process.

Invoked by LD-performance.ipynb with one JSON request path. Never import GMT
in the notebook process. Only diff/lon/lat from a regional tile are loaded.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import numpy as np
import xarray as xr


def regional_tile(source, bounds, row_batch=128):
    """Find the regional bounding box using bounded coordinate strips.

    Keep the curvilinear grid intact; the plotter applies the exact geographic
    mask within the tile. Do not stride, interpolate, or change persistence.
    """
    fields = source[["diff", "lon", "lat"]]
    dims = fields["diff"].dims
    if len(dims) != 2 or any(set(fields[v].dims) != set(dims) for v in ("lon", "lat")):
        raise ValueError("Expected diff, lon and lat on the same two-dimensional grid")
    if row_batch < 1:
        raise ValueError("row_batch must be positive")
    west, east, south, north = bounds
    row0, row1, col0, col1 = fields.sizes[dims[0]], 0, fields.sizes[dims[1]], 0
    for start in range(0, fields.sizes[dims[0]], row_batch):
        strip = fields[["lon", "lat"]].isel({dims[0]: slice(start, start + row_batch)}).transpose(*dims).compute()
        lon = (strip.lon.values + 180) % 360 - 180
        lat = strip.lat.values
        inside_lon = ((lon >= west) & (lon <= east)) if west <= east else ((lon >= west) | (lon <= east))
        if east == 180:
            inside_lon |= lon == -180
        rows, cols = np.nonzero(inside_lon & (lat >= south) & (lat <= north))
        if rows.size:
            row0, row1 = min(row0, start + int(rows.min())), max(row1, start + int(rows.max()) + 1)
            col0, col1 = min(col0, int(cols.min())), max(col1, int(cols.max()) + 1)
    if row1 == 0:
        raise ValueError(f"No grid cells in region {bounds}")
    return fields.isel({dims[0]: slice(row0, row1), dims[1]: slice(col0, col1)}).transpose(*dims).compute()


def plot_difference_file(plotter, tile, request, output):
    """Bypass PyGMT's array-backed virtual table for the FIP symbol layer."""
    pygmt = plotter._require_pygmt()
    print(f"PyGMT {pygmt.__version__}; shuga plotter: {sys.modules[plotter.__class__.__module__].__file__}", flush=True)
    bounds = request["bounds"]
    data = plotter.pygmt_da_prep(tile["diff"], lon=tile.lon, lat=tile.lat,
                               mask_zero=False, region=bounds)
    if data.empty:
        raise ValueError("No finite FIP difference points within the plot region")
    # Explicit column order and 17 significant digits preserve float64 values,
    # including zero/negative differences. No resampling or persistence changes.
    with TemporaryDirectory(prefix="ld_fip_points_") as directory:
        points = Path(directory) / "difference.xyz"
        cpt = Path(directory) / "difference.cpt"
        data[["lon", "lat", "z"]].to_csv(points, sep="\t", header=False,
                                        index=False, float_format="%.17g")
        print(f"Wrote {len(data)} lon/lat/difference rows; file-backed plotting", flush=True)
        fig = pygmt.Figure()
        projection = plotter.projection_from_region(bounds, fig_size=plotter.plt_cfg.fip_fig_size)
        print("GMT: makecpt", flush=True)
        pygmt.makecpt(cmap=str(Path(plotter.pth_cfg.fip_diff_con_cmap).expanduser()),
                      series=[-1.0, 1.0, 0.05], continuous=True, output=str(cpt))
        print("GMT: basemap/coast background", flush=True)
        plotter.pygmt_base_layer(fig, bounds, projection, title=request["title"])
        print("GMT: plot file (no x/y/fill arrays)", flush=True)
        fig.plot(data=str(points), style=plotter.plt_cfg.grid_style,
                 cmap=str(cpt), region=bounds, projection=projection)
        print("GMT: coastline overlay", flush=True)
        fig.coast(region=bounds, projection=projection, shorelines=plotter.plt_cfg.shorelines)
        print("GMT: colorbar", flush=True)
        fig.colorbar(cmap=str(cpt), position="JMB+w8c/0.4c+o0.8c/0c",
                      frame=["xaf+lFIP difference (model - AF2020)"])
        print("GMT: savefig", flush=True)
        fig.savefig(output)


def render(request):
    # Configure import paths before loading this standalone shuga checkout.
    sys.path.insert(0, str(Path(request["shuga_dir"]).expanduser().resolve().parent))
    import dask
    from shuga.core.types import RunSpec, ClassificationSpec, MetricsSpec
    from shuga.core.paths import ShugaPaths
    from shuga.plotting.cice import CICEPlotter

    cfg = request["config"]
    method = request["method"]
    run = RunSpec(request["name"], request["start"], request["end"],
                  project=cfg["project"], user=cfg["user"], experiment_group=request["group"])
    cls = ClassificationSpec(ice_type="FI", methods=method, grid_type=cfg["grid_type"],
                             ispd_thresh=cfg["ispd_thresh"], aice_thresh=cfg["aice_thresh"],
                             bin_window=cfg["bin_window"], bin_min_days=cfg["bin_min_days"],
                             roll_window=cfg["roll_window"])
    met = MetricsSpec(methods=method)
    paths = ShugaPaths(run_cfg=run, cls_cfg=cls, met_cfg=met,
                       afim_output_root=Path(request["bank"]), archive_root=Path(request["archive"]),
                       graphics_root=cfg["graphics_root"])
    print("Opening comparison store and selecting regional tile", flush=True)
    with dask.config.set(scheduler="synchronous"):
        with xr.open_zarr(request["store"], consolidated=False) as source:
            tile = regional_tile(source, request["bounds"])
    print(f"Loaded tile {dict(tile.sizes)} ({tile.nbytes / 2**20:.1f} MiB); starting GMT", flush=True)
    output = Path(request["output"])
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.stem + ".partial.png")
    try:
        plotter = CICEPlotter(run_cfg=run, cls_cfg=cls, met_cfg=met, pth_cfg=paths)
        plot_difference_file(plotter, tile, request, temporary)
        if not temporary.is_file() or temporary.stat().st_size == 0:
            raise RuntimeError("GMT returned without a nonempty PNG")
        temporary.replace(output)
    finally:
        temporary.unlink(missing_ok=True)
    print(f"Saved: {output}", flush=True)


if __name__ == "__main__":
    import faulthandler
    faulthandler.enable()
    render(json.loads(Path(sys.argv[1]).read_text()))
