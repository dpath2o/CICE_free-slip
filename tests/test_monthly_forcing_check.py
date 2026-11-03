"""Synthetic full workflow checks; no Gadi or production archive needed."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest
import xarray as xr

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('forcing_compare', ROOT / 'scripts/compare_monthly_forcing_check.py')
comparison = importlib.util.module_from_spec(spec)
spec.loader.exec_module(comparison)
spec = importlib.util.spec_from_file_location('forcing_setup', ROOT / 'scripts/setup_monthly_forcing_check.py')
setup = importlib.util.module_from_spec(spec)
spec.loader.exec_module(setup)


def dataset():
    times = pd.date_range('2003-09-02', periods=30)
    return xr.Dataset({'aice': (('time', 'nj', 'ni'), np.ones((30, 2, 3)), {'units': '1'})},
                      coords={'time': times, 'nj': [0, 1], 'ni': [0, 1, 2]})


def run_comparison(tmp_path, reference):
    history = tmp_path / 'history'
    history.mkdir()
    ds = dataset()
    ds.to_netcdf(history / 'iceh.2003-09.nc')
    root = tmp_path / 'reference'
    for month, indices in [('2003-09', slice(0, 29)), ('2003-10', slice(29, 30))]:
        reference.isel(time=indices).to_zarr(root / month, consolidated=False)
    output = tmp_path / 'report'
    r = subprocess.run([sys.executable, str(ROOT / 'scripts/compare_monthly_forcing_check.py'),
                        '--history', str(history), '--reference', str(root),
                        '--variables', 'aice', '--output', str(output)], capture_output=True, text=True)
    return r, output


def test_month_shards_and_cli_report(tmp_path):
    r, output = run_comparison(tmp_path, dataset())
    assert r.returncode == 0, r.stderr
    summary = json.loads((output / 'summary.json').read_text())
    assert summary['passed']
    assert len(summary['reference_shards']) == 2
    assert len(pd.read_csv(output / 'daily_metrics.csv')) == 30


def test_difference_produces_nonzero_exit_and_first_failure(tmp_path):
    ref = dataset()
    ref.aice.values[5, 0, 0] = 0.8
    r, output = run_comparison(tmp_path, ref)
    assert r.returncode == 1, r.stderr
    first = json.loads((output / 'summary.json').read_text())['first_failure']
    assert first['time'] == '2003-09-07T00:00:00'
    assert first['outside_tolerance'] == 1


@pytest.mark.parametrize('problem', ['missing_day', 'duplicate', 'units', 'grid', 'nan_mask', 'variable'])
def test_invalid_comparisons_cannot_pass(problem):
    model, ref = dataset(), dataset()
    if problem == 'missing_day':
        ref = ref.isel(time=slice(1, None))
    elif problem == 'duplicate':
        ref = xr.concat([ref, ref.isel(time=[0])], dim='time')
    elif problem == 'units':
        ref.aice.attrs['units'] = '%'
    elif problem == 'grid':
        ref = ref.assign_coords(ni=[2, 3, 4])
    elif problem == 'nan_mask':
        ref.aice.values[0, 0, 0] = np.nan
    elif problem == 'variable':
        ref = ref.rename(aice='other')
    if problem == 'nan_mask':
        rows = comparison.compare(model, ref, ['aice'], pd.Timestamp('2003-09-01'), 30, 1e-8, 1e-6)
        assert not rows[0]['passed']
        assert rows[0]['nan_mask_mismatch'] == 1
    else:
        with pytest.raises(ValueError):
            comparison.compare(comparison.canonical_time(model), comparison.canonical_time(ref),
                               ['aice'], pd.Timestamp('2003-09-01'), 30, 1e-8, 1e-6)


def test_restart_metadata_checked_before_forcing_access(tmp_path):
    path = tmp_path / 'iced.2003-09-01-00000.nc'
    xr.Dataset(attrs={'myear': 2009, 'mmonth': 9, 'mday': 1, 'msec': 0}).to_netcdf(path)
    with pytest.raises(ValueError, match='2003-09-01'):
        setup.check_inputs(path, tmp_path, tmp_path)


def test_test_case_preserves_cs_high_physics():
    import re
    def nml(path):
        return {m[1].lower(): m[2].strip() for m in re.finditer(
            r'(?m)^\s*(\w+(?:\(\d+\))?)\s*=\s*(.*)$', path.read_text())}
    original = nml(ROOT / 'tests/forcing_monthly/ice_in.Cs-high')
    actual = nml(ROOT / 'forcing-monthly-check/ice_in')
    allowed = {'year_init', 'month_init', 'day_init', 'sec_init', 'istep0', 'npt_unit', 'npt',
               'runtype', 'ice_ic', 'use_restart_time', 'dumpfreq', 'dump_last',
               'atm_data_dir', 'ocn_data_dir'}
    assert original.keys() == actual.keys()
    assert {k for k in original if original[k] != actual[k]} <= allowed


def test_partial_year_resume_and_leap_year():
    from datetime import datetime
    sys.path.insert(0, str(ROOT / 'scripts'))
    from resume_monthly_forcing import year_segments
    assert list(year_segments(datetime(1999, 12, 1), 2000)) == [
        (1999, 31 * 86400), (2000, 366 * 86400)]
    assert list(year_segments(datetime(1996, 9, 1, 0, 30), 1996)) == [
        (1996, 122 * 86400 - 1800)]
