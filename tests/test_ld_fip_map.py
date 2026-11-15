"""Regional cropping checks; no shuga, GMT or production data needed."""
import importlib.util
from pathlib import Path

import dask.array as da
import numpy as np
import pytest
import xarray as xr

spec = importlib.util.spec_from_file_location("ld_fip_map", Path(__file__).parents[1] / "notebooks/ld_fip_map.py")
worker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(worker)


def grid():
    lon, lat = np.meshgrid([0, 10, 20, 170, 180, 190, 350], [-80, -75, -70, -65, -60])
    diff = np.arange(lon.size, dtype=float).reshape(lon.shape) / 100
    diff[2, 2] = np.nan
    return xr.Dataset({v: (("y", "x"), da.from_array(a, chunks=(2, 3)))
                       for v, a in {"lon": lon, "lat": lat, "diff": diff}.items()})


@pytest.mark.parametrize("bounds", [[0, 20, -75, -65], [170, -170, -80, -70], [170, 180, -80, -70]])
def test_tile_preserves_values_and_all_selected_cells(bounds):
    source = grid()
    tile = worker.regional_tile(source, bounds, row_batch=2)
    lon = (source.lon.compute().values + 180) % 360 - 180
    lat = source.lat.compute().values
    west, east, south, north = bounds
    inside = ((lon >= west) & (lon <= east)) if west <= east else ((lon >= west) | (lon <= east))
    if east == 180:
        inside |= lon == -180
    rows, cols = np.nonzero(inside & (lat >= south) & (lat <= north))
    expected = source.isel(y=slice(rows.min(), rows.max()+1), x=slice(cols.min(), cols.max()+1)).compute()
    xr.testing.assert_equal(tile, expected)
    assert not hasattr(tile["diff"].data, "compute")


def test_field_reads_are_limited_to_the_tile():
    from dask import delayed
    def forbidden():
        raise AssertionError("Read diff outside the regional crop")
    good = da.from_array(np.ones((2, 7)), chunks=(2, 7))
    bad = da.from_delayed(delayed(forbidden)(), shape=(3, 7), dtype=float)
    source = grid().assign(diff=(("y", "x"), da.concatenate([good, bad], axis=0)))
    tile = worker.regional_tile(source, [0, 20, -80, -75], row_batch=2)
    assert tile.sizes == {"y": 2, "x": 3}
    assert (tile["diff"] == 1).all()


def test_no_cells_and_incompatible_dimensions():
    with pytest.raises(ValueError, match="No grid cells"):
        worker.regional_tile(grid(), [30, 40, -80, -70])
    with pytest.raises(ValueError, match="same two-dimensional grid"):
        worker.regional_tile(grid().assign(diff=(("time", "y", "x"), np.zeros((2, 5, 7)))), [0, 20, -80, -70])


@pytest.mark.parametrize("crash", [False, True])
def test_notebook_registers_only_successful_worker_outputs(tmp_path, crash):
    import json
    import os
    import subprocess
    import sys
    # An actual child process exits abnormally; the test/notebook parent survives.
    stub = ('import os, signal, resource\nresource.setrlimit(resource.RLIMIT_CORE, (0, 0))\n'
            'os.kill(os.getpid(), signal.SIGSEGV)\n') if crash else (
                'import json, sys\nfrom pathlib import Path\n'
                'r=json.loads(Path(sys.argv[1]).read_text())\nPath(r["output"]).write_bytes(b"PNG")\n')
    (tmp_path / "ld_fip_map.py").write_text(stub)
    notebook = json.loads((Path(__file__).parents[1] / "notebooks/LD-performance.ipynb").read_text())
    cell = next("".join(c["source"]) for c in notebook["cells"] if "FIP_MAP_TIMEOUT_S =" in "".join(c["source"]))
    registrations = []
    ns = dict(Path=Path, os=os, sys=sys, json=json, subprocess=subprocess,
              RUN_FIP=True, FIP_DATA={"case": "store"},
              FIP_PRODUCTS={"case": {"store": "store", "plot_root": tmp_path}},
              D_NB_ROOT=tmp_path, D_shuga=tmp_path, CFG={}, GROUP="group", METHOD="binary-days",
              START="2000-01-01", END="2005-12-31", OBS_START="2000-01-01", OBS_END="2005-12-31",
              BANK=tmp_path, ARCHIVE=tmp_path, REGIONS=["Aus"],
              ANTARCTIC_8_REGIONS={"Aus": {"plot_region": [97, 142, -68, -63]}},
              fip_env=os.environ.copy(), FIP_SHOW_MAPS=False,
              register=lambda *args: registrations.append(args))
    if crash:
        # A stale PNG must never be registered after this worker fails.
        stale = tmp_path / "Aus/2000-01-01_2005-12-31_case_Aus.png"
        stale.parent.mkdir()
        stale.write_bytes(b"OLD PNG")
        with pytest.raises(RuntimeError, match="SIGSEGV"):
            exec(compile(cell, "<plotting cell>", "exec"), ns)
        assert registrations == []
        assert ns["FIP_MAPS"] == {}
    else:
        exec(compile(cell, "<plotting cell>", "exec"), ns)
        assert len(registrations) == 1
        assert ns["FIP_MAPS"][("case", "Aus")].read_bytes() == b"PNG"
