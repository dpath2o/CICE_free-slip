#!/usr/bin/env python3
"""Plan restart-based calendar-year jobs for the three existing free-slip cases."""
import argparse
from datetime import datetime, timedelta
import re
import shlex
import subprocess
from pathlib import Path

from setup_monthly_forcing_check import set_nml


def value(text, key):
    m = re.search(rf'(?mi)^\s*{re.escape(key)}\s*=\s*([^!\n]+)', text)
    if not m:
        raise ValueError(f'Missing {key}')
    return m[1].strip().strip("'\"")


def year_segments(start, end_year):
    if start >= datetime(end_year + 1, 1, 1):
        raise ValueError('Restart is already beyond the requested end year')
    while start.year <= end_year:
        stop = datetime(start.year + 1, 1, 1)
        yield start.year, int((stop - start).total_seconds())
        start = stop


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case', type=Path, required=True)
    p.add_argument('--end', type=int, required=True, help='Last calendar year to process')
    p.add_argument('--submit', action='store_true', help='Write snapshots and submit chained PBS jobs')
    a = p.parse_args()
    case = a.case.resolve()
    if case.name not in ('free-slip', 'free-slip01', 'free-slip02'):
        p.error('This helper is restricted to free-slip, free-slip01 and free-slip02')
    settings = (case / 'cice.settings').read_text()
    run = Path(re.search(r'(?m)^setenv\s+ICE_RUNDIR\s+(\S+)', settings)[1])
    nml = (case / 'ice_in').read_text()
    if value(nml, 'use_leap_years').lower() != '.true.':
        p.error('This helper requires the existing Gregorian/leap-year calendar')
    def in_run(path):
        path = Path(path)
        return path if path.is_absolute() else run / path
    pointer = in_run(value(nml, 'pointer_file'))
    rst_name = pointer.read_text().strip().strip("'\"")
    restart = in_run(rst_name)
    if not restart.is_file():
        restart = in_run(value(nml, 'restart_dir')) / Path(rst_name).name
    import xarray as xr
    with xr.open_dataset(restart, decode_times=False) as ds:
        start = datetime(*(int(ds.attrs[k]) for k in ('myear', 'mmonth', 'mday')))
        start += timedelta(seconds=int(ds.attrs['msec']))
    dt = float(value(nml, 'dt'))
    segments = list(year_segments(start, a.end))
    if any(seconds % dt for _, seconds in segments):
        p.error('Calendar boundary does not fall on a complete model timestep')
    atm = Path(value(nml, 'atm_data_dir'))
    ocn = Path(value(nml, 'ocn_data_dir'))
    first_year = int(value(nml, 'fyear_init'))
    cycle = int(value(nml, 'ycycle'))
    # CICE samples end-of-step forcing, including Jan 1 of the next year.
    for year, _ in segments:
        fyr = year if cycle == 0 else first_year + (year - first_year) % cycle
        first_month = start.month if year == start.year else 1
        paths = [atm / f'era5_for_cice6_{fyr:04d}_{m:02d}.nc' for m in range(first_month, 13)]
        paths += [ocn / f'ORAS_{fyr:04d}.nc']
        next_fyr = year + 1 if cycle == 0 else first_year + (year + 1 - first_year) % cycle
        paths += [atm / f'era5_for_cice6_{next_fyr:04d}_01.nc', ocn / f'ORAS_{next_fyr:04d}.nc']
        for path in paths:
            if not path.is_file():
                raise FileNotFoundError(path)
    if not (run / 'cice').is_file():
        p.error(f'No executable at {run / "cice"}; run ./cice.build first')
    print(f'Restart: {restart}\nResume time: {start}\nEnd: {a.end + 1}-01-01')
    print('Rebuild with the upgraded source before submission; this helper does not build.')
    dest = case / f'resume-monthly-{datetime.now():%Y%m%dT%H%M%S%f}'
    if a.submit:
        dest.mkdir()
    prev = None
    for year, seconds in segments:
        snapshot = nml
        for key, rhs in {'runtype': "'continue'", 'ice_ic': "'none'", 'use_restart_time': '.true.',
                         'npt_unit': "'s'", 'npt': str(seconds), 'dump_last': '.true.'}.items():
            snapshot = set_nml(snapshot, key, rhs)
        path = dest / f'ice_in.{year}'
        cmd = ['qsub']
        if prev:
            cmd += ['-W', f'depend=afterok:{prev}']
        cmd += ['-N', f'{case.name}-{year}', '-v', f'ICE_IN_SRC={path}', 'cice.run.continuous']
        print(shlex.join(cmd))
        if a.submit:
            path.write_text(snapshot)
            result = subprocess.run(cmd, cwd=case, check=True, text=True, capture_output=True)
            prev = result.stdout.strip().split()[0]
            print(f'Submitted {prev}')
        else:
            prev = f'DRYRUN.{year}'
    if not a.submit:
        print('Plan only: no files written and no jobs submitted. Add --submit to execute.')


if __name__ == '__main__':
    main()
