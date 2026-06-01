import os, sys, warnings
import numpy             as np
import pandas            as pd
import xarray            as xr
import matplotlib.pyplot as plt
from IPython.display     import Video, Image, display
from dataclasses         import dataclass, replace 
from pathlib             import Path
from typing              import Mapping, Sequence
from math                import ceil
xr.set_options(keep_attrs=True)
warnings.filterwarnings(\"default\")
repo_root = Path.home() / "AFIM" / "src" / "mawsons-chest"
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
from shuga import (RunSpec,
                   ClassificationSpec,
                   MetricsSpec,
                   CICEMetrics,
                   ShugaPaths,
                   CICEPlotter,
                   PlottingSpec,
                   CICEPlotter,
                   ObservationSpec,
                   SeaIceObservations,
                   load_cice,
                   load_classified,
                   load_metrics,
                   report_sim_status)
from shuga.grid.lateral_drag import FormFactors
from shuga.core.regions      import ANTARCTIC_8_REGIONS
