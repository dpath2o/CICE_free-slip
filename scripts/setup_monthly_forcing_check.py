#!/usr/bin/env python3
"""Create an isolated Cs-high forcing check through cice.setup, on Gadi."""
import argparse
import json
import re
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ATM = Path('/g/data/gv90/da1339/afim_input/ERA5/0p25/bilinear/monthly_cice6_fast')
ARCHIVE = Path.home() / 'AFIM_archive/lateral-drag/Cs-high'


def set_nml(text, key, value):
    text, count = re.subn(rf'(?mi)^(\s*{re.escape(key)}\s*=).*$',
                          lambda m: f'{m[1]} {value}', text)
    if count != 1:
        raise ValueError(f'Expected exactly one {key}; found {count}')
    return text


def check_inputs(restart, atm, ocn):
    import xarray as xr
    with xr.open_dataset(restart, decode_times=False) as ds:
        date = tuple(int(ds.attrs[k]) for k in ('myear', 'mmonth', 'mday', 'msec'))
        if date != (2003, 9, 1, 0):
            raise ValueError(f'Restart date {date} is not 2003-09-01 00Z. '
                             'A different initial state cannot reproduce Cs-high September 2003.')
    for path, variables, records in [
        (atm / 'era5_for_cice6_2003_09.nc',
         ['airtmp', 'wndewd', 'wndnwd', 'spchmd', 'glbrad', 'dlwsfc', 'ttlpcp'], 720),
        # CICE samples forcing at the end of a step, including Oct 1 00Z.
        (atm / 'era5_for_cice6_2003_10.nc',
         ['airtmp', 'wndewd', 'wndnwd', 'spchmd', 'glbrad', 'dlwsfc', 'ttlpcp'], 2),
        (ocn / 'ORAS_2003.nc', ['sst', 'sss', 'u', 'v'], 365),
    ]:
        with xr.open_dataset(path, decode_times=False) as ds:
            for var in variables:
                if var not in ds or ds[var].shape[0] < records:
                    raise ValueError(f'{path}: missing/incomplete {var}; need {records} records')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case', default='forcing-monthly-check')
    p.add_argument('--restart', type=Path, default=ARCHIVE / 'restart/iced.2003-09-01-00000.nc')
    p.add_argument('--atm-dir', type=Path, default=ATM)
    p.add_argument('--ocn-dir', type=Path,
                   default=Path('/g/data/gv90/da1339/afim_input/ORAS/daily/sfc/fin'))
    p.add_argument('--run-root', type=Path, default=Path('/g/data/gv90/da1339/cice-dirs/runs'))
    p.add_argument('--check-only', action='store_true', help='Validate Gadi inputs without creating a case')
    p.add_argument('--skip-input-check', action='store_true',
                   help='Generate case only; validate inputs separately before building/running')
    a = p.parse_args()
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*', a.case) or a.case.startswith('free-slip'):
        p.error('Use a separate, simple case name (not free-slip*)')
    case = ROOT / a.case
    rundir = a.run_root.resolve() / a.case
    restart = a.restart.expanduser().resolve()
    if a.check_only:
        check_inputs(restart, a.atm_dir, a.ocn_dir)
        print('Restart date and forcing file checks passed')
        return
    if case.exists() or rundir.exists():
        p.error('Case or run directory already exists; choose a new --case (nothing overwritten)')
    if not a.skip_input_check:
        check_inputs(restart, a.atm_dir, a.ocn_dir)
    cmd = ['./cice.setup', '--case', str(case), '--mach', 'gadi1', '--env', 'intel',
           '--grid', 'gbox12', '--pes', '1232x1x20x20x10', '--acct', 'jk72',
           '--queue', 'normalbw', '--ignore-user-set']
    subprocess.run(cmd, cwd=ROOT, check=True)
    if not (case / 'cice.settings').is_file():
        raise RuntimeError('cice.setup did not produce the expected case')
    # Match the working free-slip build and PBS launch scripts. Keep project jk72.
    for name in ['cice.build', 'Makefile', 'makdep.c', 'Macros.gadi1_intel',
                 'env.gadi1_intel', 'setup_run_dirs.csh', 'cice.run', 'cice.run.continuous']:
        shutil.copy2(ROOT / 'free-slip' / name, case / name)
    settings = (ROOT / 'free-slip/cice.settings').read_text()
    replacements = {'ICE_CASENAME': a.case, 'ICE_SANDBOX': str(ROOT),
                    'ICE_SCRIPTS': str(ROOT / 'configuration/scripts'),
                    'ICE_CASEDIR': str(case), 'ICE_RUNDIR': str(rundir),
                    'ICE_ACCOUNT': 'jk72', 'ICE_CLEANBUILD': 'true'}
    for key, value in replacements.items():
        settings, count = re.subn(rf'(?m)^setenv\s+{key}\s+.*$',
                                  f'setenv {key} {value}', settings)
        if count != 1:
            raise ValueError(f'Missing setting {key}')
    (case / 'cice.settings').write_text(settings)
    for name in ['cice.run', 'cice.run.continuous']:
        t = (case / name).read_text()
        t = t.replace('/home/581/da1339/AFIM/src/CICE_free-slip/free-slip', str(case))
        t = re.sub(r'(?m)^#PBS -N .*$', f'#PBS -N {a.case}', t)
        (case / name).write_text(t)
    nml = (ROOT / 'tests/forcing_monthly/ice_in.Cs-high').read_text()
    for key, value in {'year_init': '2003', 'month_init': '9', 'day_init': '1',
                       'sec_init': '0', 'istep0': '0', 'npt_unit': "'d'", 'npt': '30',
                       'runtype': "'initial'", 'ice_ic': repr(str(restart)),
                       'use_restart_time': '.true.', 'dumpfreq': "'m','x','x','x','x'",
                       'dump_last': '.true.', 'atm_data_dir': repr(str(a.atm_dir.resolve()) + '/'),
                       'ocn_data_dir': repr(str(a.ocn_dir.resolve()) + '/')}.items():
        nml = set_nml(nml, key, value)
    (case / 'ice_in').write_text(nml)
    # cice.continuous in production chains full years; this check is one 30-day job.
    (case / 'cice.continuous').write_text(
        '#!/usr/bin/env bash\nset -euo pipefail\n'
        'cd "$(dirname "$0")"\n'
        'if (( $# )); then echo "This case runs one 30-day job; no yearly options." >&2; exit 2; fi\n'
        'qsub cice.run.continuous\n')
    (case / 'cice.continuous').chmod(0o755)
    (case / 'forcing-check.json').write_text(json.dumps(
        {'setup_command': cmd, 'restart': str(restart), 'start': '2003-09-01',
         'days': 30, 'inputs_checked': not a.skip_input_check,
         'rundir': str(rundir), 'reference': str(ARCHIVE / 'zarr/iceh_daily.zarr')}, indent=2) + '\n')
    print(f'Created {case}\nRun directory: {rundir}\nPBS project: jk72')
    print(f'Next: cd {case} && ./cice.build && qsub cice.run')


if __name__ == '__main__':
    main()
