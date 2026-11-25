"""Verify file-backed GMT input without importing the native GMT library."""
import ast
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

import numpy as np
import pandas as pd


source = Path(__file__).parents[1] / "notebooks/ld_fip_map.py"
tree = ast.parse(source.read_text())
function = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "plot_difference_file")
namespace = dict(Path=Path, TemporaryDirectory=TemporaryDirectory, sys=sys)
exec(compile(ast.Module(body=[function], type_ignores=[]), str(source), "exec"), namespace)


class FilePlotTests(unittest.TestCase):
    def test_points_and_palette_pass_as_files(self):
        expected = pd.DataFrame({"lon": [1.2345678901234567, 2., 3.],
                                 "lat": [-71., -70., -69.], "z": [-.7, 0., .9]})
        events = []
        class Figure:
            def plot(self, **kwargs):
                self_check = set(kwargs) & {"x", "y", "fill"}
                assert not self_check, "No arrays may be passed to the GMT symbol layer"
                actual = np.loadtxt(kwargs["data"])
                np.testing.assert_array_equal(actual, expected[["lon", "lat", "z"]].to_numpy())
                assert Path(kwargs["cmap"]).is_file()
                events.append("plot")
            def coast(self, **kwargs): events.append("coast")
            def colorbar(self, **kwargs):
                assert Path(kwargs["cmap"]).is_file()
                events.append("colorbar")
            def savefig(self, path): events.append("savefig")
        def makecpt(**kwargs):
            assert kwargs["series"] == [-1., 1., .05]
            Path(kwargs["output"]).write_text("CPT placeholder")
        class Plotter:
            plt_cfg = SimpleNamespace(fip_fig_size=20, grid_style="s0.05c", shorelines="0.5p")
            pth_cfg = SimpleNamespace(fip_diff_con_cmap="difference.cpt")
            def _require_pygmt(self):
                return SimpleNamespace(__version__="test", Figure=Figure, makecpt=makecpt)
            def pygmt_da_prep(self, *args, **kwargs):
                assert kwargs["mask_zero"] is False
                return expected[["z", "lat", "lon"]]
            def projection_from_region(self, *args, **kwargs): return "S0/-90/20c"
            def pygmt_base_layer(self, *args, **kwargs): events.append("background")
        tile = SimpleNamespace(lon=None, lat=None)
        class Tile:
            lon = lat = None
            def __getitem__(self, key): return None
        namespace["plot_difference_file"](Plotter(), Tile(), {"bounds": [-30, 10, -75, -68], "title": "Test"}, "unused.png")
        self.assertEqual(events, ["background", "plot", "coast", "colorbar", "savefig"])


if __name__ == "__main__":
    unittest.main()
