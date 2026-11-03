# Monthly ERA5 forcing upgrade

The `develop` forcing reader now opens `era5_for_cice6_YYYY_MM.nc`.
The monthly filename, record indexing and two-slot cache are adapted from
`dpath2o/CICE_dyntens` commit `a0de0cf07ac920b11e3c146b7eda44931d303bd1`.

The port is deliberately limited to I/O. Total precipitation (`ttlpcp`) still
enters the existing temperature-based rain/snow partition. The newer branch's
separate ERA5 snowfall/rainfall pathway, precipitation perturbations, tides,
waves and dynamic tensile-strength changes are not included. The ocean reader
is already identical in both branches: it still requires annual `ORAS_YYYY.nc`
files. `ice_init.F90`, `ice_history.F90`, dynamics and Icepack are unchanged.

The source reader holds the final hourly state record at month end, rather
than opening the next month for interpolation. This behaviour is retained.
Consequently the last sub-hour of each month can differ from the annual reader.
Unlike the old annual reader, the monthly reader interpolates the first hour
of January. These are time-sampling differences, not changes to the model's
physical parameterisations. The cache includes year, month and record number,
so month/year rollover cannot reuse a record from a different file.

The base `ice_in` files for `free-slip`, `free-slip01`, `free-slip02` and the
tracked `free-slip/continuous/ice_in/ice_in.*` snapshots point to:

```
/g/data/gv90/da1339/afim_input/ERA5/0p25/bilinear/monthly_cice6_fast/
```

Other experiment directories have not been migrated. Any other case using
`atm_data_type='ERA5'` with this executable must also point to a directory of
monthly files. There is no fallback to annual ERA5 files.

## Thirty-day Cs-high check

The case is `forcing-monthly-check/`, generated with `cice.setup`, then given
the working `free-slip` build/PBS scripts and the supplied Cs-high namelist.
The original pasted namelist is preserved at
`tests/forcing_monthly/ice_in.Cs-high`. Static paths target
`/home/581/da1339/AFIM/src/CICE_free-slip`; the run directory is separate from
all production cases and archives. PBS project **jk72**, queue **normalbw**,
1232 ranks, 700 GB and working compiler/MPI/NetCDF modules are retained.
The build requests a clean rebuild.

The default restart is `~/AFIM_archive/lateral-drag/Cs-high/restart/iced.2003-09-01-00000.nc`.
The request mentioned 2009, but the supplied listing contains 2000–2005 and
includes the matching 2003 restart. Reusing a 2009 state would invalidate a
like-for-like comparison with the September 2003 archive. Input validation
checks the date inside the restart, not just its filename.

The run starts 1 September 2003, lasts 30 days and ends 1 October 00Z.
`runtype='initial'` with `ice_ic` set to the restart path reads that state;
`use_restart_time=.true.` retains its model time. This avoids modifying the
archive or installing an archive file into a writable restart directory.
All physics and history-field settings are taken from the supplied namelist.
Only run dates/length, restart source, forcing paths and final restart output
are adjusted. Final restart/history/log output goes to the new case/run only.

On Gadi, after pulling `develop`:

```bash
cd ~/AFIM/src/CICE_free-slip
module use /g/data/xp65/public/modules
module load conda/analysis3-26.02
python scripts/setup_monthly_forcing_check.py --check-only
cd forcing-monthly-check
./cice.build
qsub cice.run
```

If you pulled the initial upgrade and this case has no `ice_in`, the file
was omitted from that commit because of the repository's `ice_in` ignore rule.
It is now explicitly tracked. From the repository root, retrieve just the
missing file without merging or changing your other local edits:

```bash
git fetch free develop
git restore --source=FETCH_HEAD -- forcing-monthly-check/ice_in
ls -l forcing-monthly-check/ice_in
```

This restores the prepared 30-day Cs-high namelist. If you have since created
or edited that file, save your version before restoring it. **An existing
successful build does not need to be repeated for this namelist-only repair.**
The namelist reads the 2003 restart directly from its archive path; a copied
restart under the run directory is not required.

`qsub cice.run.continuous` is an equivalent single-job launch here.
`./cice.continuous` also submits just this 30-day job; it rejects annual-chain
options. Do not submit both launch scripts for the same run.

The preflight first requires a nonempty case `ice_in`, then checks restart
metadata, September ERA5 records/variables,
October ERA5 availability and annual ORAS variables/records. CICE samples
forcing at the end of each timestep, so the final step needs October 1 forcing.
It does not establish that the monthly reprocessed data equal the old inputs.
Check available disk space before running: 30 daily full-grid files are substantial.

If your checkout is elsewhere, or this run directory already contains output,
create a fresh case through the setup helper. It refuses to overwrite an
existing case or run directory and does not submit jobs:

```bash
cd /path/to/CICE_free-slip
python scripts/setup_monthly_forcing_check.py --case forcing-monthly-check-02
cd forcing-monthly-check-02
./cice.build
qsub cice.run
```

After successful completion:

```bash
cd ~/AFIM/src/CICE_free-slip
module use /g/data/xp65/public/modules
module load conda/analysis3-26.02
python scripts/compare_monthly_forcing_check.py \
  --history /g/data/gv90/da1339/cice-dirs/runs/forcing-monthly-check/history \
  --reference ~/AFIM_archive/lateral-drag/Cs-high/zarr/iceh_daily.zarr \
  --reference-time-shift-hours 24 \
  --output forcing-check-report-2003-09
```

The comparison opens independent monthly Zarr stores (`2003-09`, and
`2003-10` if present), including October 1 when daily means are end-labelled.
Default variables are `aice,hi,uvel,vvel,uatm,vatm,sst,sss,uocn,vocn,rain,snow`.
Use `--variables` to request a different explicit set. All requested variables
must exist in both datasets. No variables are silently omitted.

`daily_metrics.csv` contains per-variable/day bias, RMSE, maximum absolute
difference, finite-cell count, NaN-mask differences and tolerance failures.
Statistics are unweighted grid-cell comparisons, processed one day and field
at a time. `summary.json` records file paths, tolerances and the first failure.
Exit status is 0 for a numeric pass, 1 for differences, and nonzero on invalid
inputs. An existing report directory is never overwritten.

The default numerical criterion is `abs(new-old) <= 1e-8 + 1e-6*abs(old)`;
these are configurable screening thresholds, not validated scientific
acceptance limits. Units, dimensions, indexed coordinates and all 30 daily
end-labelled timestamps must match exactly. Missing days, duplicates,
infinities, mismatched NaN masks or missing fields cannot pass. For archives written by shuga's daily `NC2Zarr` converter, use
`--reference-time-shift-hours 24`, as in the command above. Its
[conversion code](https://github.com/dpath2o/shuga/blob/main/core/data_conversion.py#L538)
subtracts one day from daily NetCDF timestamps before writing Zarr. A NetCDF
mean ending September 2 00Z is therefore labelled September 1 00Z in that
archive. Adding 24 hours to reference coordinates restores end-labelled
alignment; values are unchanged. With zero shift, equal timestamps compare
adjacent daily intervals and can produce large wind, current and ice-state
differences. The script defaults to zero for general archives; verify the
conversion convention when using a different reference. A shifted comparison
must still pass the original tolerances. Never shift dates merely to obtain agreement. No spatial interpolation is performed; confirm
that both products use the same grid and masks.

A bitwise match is not promised: the monthly ERA5 files were reprocessed,
temporal behaviour at month/year edges differs as described above, and the
archived history is single precision. Start by checking the atmospheric
fields and then the ocean fields before interpreting ice-state divergence.

## Resume the existing three runs

Rebuild each case separately after the short test is satisfactory. Retain
its own `ice_in`, restart pointer and physics settings; do not copy the
Cs-high test namelist into production cases.

The existing `cice.continuous` wrapper can reset to initial conditions if it
fails to recognise the restart configuration, and its first job always uses
the base run length. A December or September checkpoint therefore needs a
partial-year first segment. The new helper reads the pointer and NetCDF
restart metadata, requires a real checkpoint/executable and all required
forcing filenames, and plans only the remaining seconds to January 1 before
full calendar-year segments. It never changes the base namelist or pointer.

For example, to finish through December 2005:

```bash
cd ~/AFIM/src/CICE_free-slip/free-slip
./cice.build
cd ..
python scripts/resume_monthly_forcing.py --case free-slip --end 2005
# Inspect the checkpoint date and planned jobs, then submit:
python scripts/resume_monthly_forcing.py --case free-slip --end 2005 --submit
```

Repeat independently with `free-slip01`, and with `free-slip02` after its
September 1996 failure has been resolved. Change `--end` to the intended final
year. The helper uses the existing `cice.run.continuous` (including PBS project
jk72) and chains jobs with `afterok`. New snapshots go in a uniquely named
`resume-monthly-*` directory within that case. It checks the following January's
ERA5 and ORAS files too because the last timestep samples that boundary.
Do not run a second chain concurrently in the same run directory.

Restart from the latest *valid checkpoint*, not the last history date. History
after that checkpoint may be regenerated. Review/archive any partial output
in the production run directory before restarting. No script here writes to
Cs-high archives or forcing data.

## Validation status

Local checks cover Fortran 2008 syntax, preservation of all other source
routines/physics files, test-namelist equivalence, case generation through
`cice.setup`, shell syntax, and synthetic NetCDF-to-monthly-Zarr comparisons
including missing data, time/grid/unit mismatches and numeric failures.
The case generator was exercised with a local module-command stub; this
validates generation, not the Gadi toolchain. No Gadi compile or model run has
been performed here. Those checks and the Cs-high result comparison remain
required before continuing production experiments.
