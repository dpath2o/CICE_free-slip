#!/usr/bin/env python3
"""Compare daily CICE NetCDF history with Cs-high monthly Zarr shards."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

DEFAULT_VARS = 'aice,hi,uvel,vvel,uatm,vatm,sst,sss,uocn,vocn,rain,snow'


def canonical_time(ds, shift_hours=0):
    if 'time' not in ds.coords:
        raise ValueError('Dataset has no time coordinate')
    # Compare decoded timestamps explicitly, including time-of-day. Never
    # align by array position or discard unmatched records through inner joins.
    times = pd.DatetimeIndex([str(t) for t in ds.time.values])
    times = times + pd.Timedelta(hours=shift_hours)
    if times.has_duplicates or times.isna().any():
        raise ValueError('Duplicate or invalid daily timestamps')
    return ds.assign_coords(time=times).sortby('time')


def open_reference(root, start, end):
    # Shards are independent stores, not one root array. Include both
    # endpoint months because end-labelled September means can end on Oct 1.
    months = pd.period_range(start.to_period('M'), end.to_period('M'), freq='M')
    paths = [root / str(m) for m in months if (root / str(m)).is_dir()]
    if not paths:
        raise ValueError(f'No monthly Zarr shards in {root}')
    parts = [xr.open_zarr(p, consolidated=False, chunks={}) for p in paths]
    return xr.concat(parts, dim='time', data_vars='minimal', coords='minimal', compat='equals'), paths


def compare(model, reference, variables, start, days, atol, rtol):
    end = start + pd.Timedelta(days=days)
    expected = pd.date_range(start + pd.Timedelta(days=1), end, freq='D')
    def select(ds):
        ds = ds.sel(time=(ds.time > start) & (ds.time <= end))
        if not pd.DatetimeIndex(ds.time.values).equals(expected):
            raise ValueError('Need exactly one end-labelled mean per day, '
                             f'{expected[0]} through {expected[-1]}; '
                             f'got {list(ds.time.values)}. Check time conventions; '
                             'use an explicit reference shift only if justified.')
        return ds
    model, reference = select(model), select(reference)
    rows = []
    for var in variables:
        if var not in model or var not in reference:
            raise ValueError(f'Required comparison variable missing: {var}')
        m, r = model[var], reference[var]
        if set(m.dims) != set(r.dims):
            raise ValueError(f'{var}: dimension mismatch {m.dims} versus {r.dims}')
        if m.attrs.get('units') != r.attrs.get('units'):
            raise ValueError(f'{var}: units differ: {m.attrs.get("units")} / {r.attrs.get("units")}')
        m, r = xr.align(m, r.transpose(*m.dims), join='exact')
        for i, date in enumerate(expected):
            # Materialize one field and day at a time on the full 0.25-degree grid.
            a = np.asarray(m.isel(time=i).values, dtype=np.float64)
            b = np.asarray(r.isel(time=i).values, dtype=np.float64)
            mask_a, mask_b = np.isnan(a), np.isnan(b)
            mask_mismatch = int(np.count_nonzero(mask_a != mask_b))
            if np.isinf(a).any() or np.isinf(b).any():
                raise ValueError(f'{var} {date}: infinite values')
            valid = ~mask_a & ~mask_b
            if not valid.any():
                raise ValueError(f'{var} {date}: no comparable finite cells')
            delta = a[valid] - b[valid]
            bad = np.abs(delta) > atol + rtol * np.abs(b[valid])
            rows.append({'variable': var, 'time': date.isoformat(),
                         'count': int(valid.sum()), 'nan_mask_mismatch': mask_mismatch,
                         'bias': float(delta.mean()), 'rmse': float(np.sqrt(np.mean(delta**2))),
                         'max_abs': float(np.max(np.abs(delta))),
                         'outside_tolerance': int(bad.sum()),
                         'passed': bool(not bad.any() and mask_mismatch == 0)})
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--history', type=Path, required=True)
    p.add_argument('--reference', type=Path,
                   default=Path.home() / 'AFIM_archive/lateral-drag/Cs-high/zarr/iceh_daily.zarr')
    p.add_argument('--start', default='2003-09-01')
    p.add_argument('--days', type=int, default=30)
    p.add_argument('--variables', default=DEFAULT_VARS)
    p.add_argument('--atol', type=float, default=1e-8)
    p.add_argument('--rtol', type=float, default=1e-6)
    p.add_argument('--reference-time-shift-hours', type=float, default=0,
                   help='Explicit correction for documented archive timestamp convention')
    p.add_argument('--output', type=Path, required=True, help='New report directory')
    a = p.parse_args()
    if a.days < 1 or not np.isfinite([a.atol, a.rtol]).all() or a.atol < 0 or a.rtol < 0:
        p.error('days must be positive and tolerances nonnegative')
    if a.output.exists():
        p.error('Output already exists; choose a new report directory')
    files = sorted(a.history.glob('iceh.*.nc'))
    if not files:
        p.error(f'No daily iceh.*.nc files in {a.history}')
    start = pd.Timestamp(a.start)
    variables = [v.strip() for v in a.variables.split(',') if v.strip()]
    def keep_vars(ds):
        missing = set(variables) - set(ds.data_vars)
        if missing:
            raise ValueError(f'History missing variables: {sorted(missing)}')
        return ds[variables]
    model = xr.open_mfdataset(files, combine='by_coords', chunks={'time': 1}, preprocess=keep_vars)
    ref, shards = open_reference(a.reference, start, start + pd.Timedelta(days=a.days))
    try:
        model = canonical_time(model)
        ref = canonical_time(ref, a.reference_time_shift_hours)
        rows = compare(model, ref, variables, start, a.days, a.atol, a.rtol)
        passed = all(row['passed'] for row in rows)
        summary = {'passed': passed, 'atol': a.atol, 'rtol': a.rtol,
                   'start': str(start), 'days': a.days, 'variables': variables,
                   'history_files': [str(f.resolve()) for f in files],
                   'reference_shards': [str(f.resolve()) for f in shards],
                   'reference_time_shift_hours': a.reference_time_shift_hours,
                   'metrics': 'Unweighted grid-cell statistics; no spatial interpolation',
                   'first_failure': next((r for r in rows if not r['passed']), None)}
        a.output.mkdir(parents=True)
        with (a.output / 'daily_metrics.csv').open('w') as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        (a.output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
        print(f'{"PASS" if passed else "DIFFERENCES"}: {a.output / "summary.json"}')
        raise SystemExit(0 if passed else 1)
    finally:
        model.close()
        ref.close()


if __name__ == '__main__':
    main()
