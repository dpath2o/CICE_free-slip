# Cs-high with monthly ERA5 forcing

Case directory: `free-slip04`. Experiment: `Cs-high-mnthly-frcg`.

The namelist is based on `tests/forcing_monthly/ice_in.Cs-high`, the saved
Cs-high reference. Its only changes are `runtype='initial'`,
`use_restart_time=.false.`, and the monthly ERA5 directory. In particular,
`kstrength=0`, `Cs=1e-3`, `Ktens=0.2`, static lateral drag, the form-factor
file, ORAS forcing and 14-day ocean restoring are retained. The previous
free-slip04 namelist had kstrength=1; it must not be used for this control.

Start: 1995-01-01 00:00, with `ice_ic='none'`. End: 2006-01-01 00:00
(after all of 2005-12-31). This is 11 calendar years: 1995–1999 spin-up
and 2000–2005 assessment, with leap years enabled. Each job runs 12 months;
the annual wrapper, not the base namelist alone, schedules the full period.

The run, compilation, history and restart directories are isolated under:

```text
/g/data/gv90/da1339/cice-dirs/runs/Cs-high-mnthly-frcg
```

The build case name remains free-slip04. Both PBS launchers use jk72.
Existing free-slip04 outputs are not moved or deleted. Existing yearly
namelists in free-slip04/continuous are legacy snapshots; the wrapper
regenerates 1995–2005 from the new base namelist before submission.

## Reconcile Gadi's local commit

From the repository root, with no uncommitted edits:

```bash
git branch backup/monthly-comparisons-301176d 301176d
git fetch free
git rebase free/develop
git push free develop
```

If a conflict occurs, stop and inspect `git status`. Resolve deliberately,
stage the resolved files and use `git rebase --continue`, or use
`git rebase --abort` to return to the starting point. A force push is
unnecessary. If the backup branch already exists, retain it and omit that
first command.

## Preflight and launch

From the repository root, check all 132 monthly input files and annual
ORAS inputs before building. This checks presence, not file contents.

```bash
python - <<'PY'
from pathlib import Path
root = Path("/g/data/gv90/da1339")
era = root / "afim_input/ERA5/0p25/bilinear/monthly_cice6_fast"
oras = root / "afim_input/ORAS/daily/sfc/fin"
required = [
    era / f"era5_for_cice6_{y}_{m:02d}.nc"
    for y in range(1995, 2006) for m in range(1, 13)
]
required += [oras / f"ORAS_{y}.nc" for y in range(1995, 2006)]
missing = [str(p) for p in required if not p.is_file()]
if missing:
    raise SystemExit("Missing inputs:\n" + "\n".join(missing))
run = root / "cice-dirs/runs/Cs-high-mnthly-frcg"
for name in ("history", "restart"):
    directory = run / name
    if directory.exists() and any(directory.iterdir()):
        raise SystemExit(f"Existing outputs in {directory}; inspect before a cold start")
if (run / "ice.restart_file").exists():
    raise SystemExit(f"Existing restart pointer in {run}; inspect before a cold start")
print("Inputs present; no existing outputs detected for this experiment.")
PY
```

Only after preflight succeeds:

```bash
cd free-slip04
./cice.build
bash ./cice.continuous --start 1995 --end 2005 --force-initial --dry-run
```

The dry run writes yearly namelists and prints 11 job submissions without
calling qsub. Check 1995 is initial with no ice and 1996–2005 continue from
the preceding year's restart. It retains 12-month segments throughout.
Launch once:

```bash
bash ./cice.continuous --start 1995 --end 2005 --force-initial
```

Do not also submit cice.run: that would launch another job into the same
directory. Do not repeat the cold-start command after the chain has begun.
For a failed year, inspect its checkpoint before choosing how to resume.

## Scientific comparison

Compare 2000–2005 against the archived Cs-high run using the same fast-ice
classification, temporal filters, observation masks and time conventions.
Inspect regional extent, seasonality, persistence and spatial changes as
well as the hemispheric mean. Keep the 1995–1999 spin-up separate.

This measures the effect of the monthly forcing workflow as a whole:
bilinear preprocessing and monthly-reader sampling, including month-end
handling. It does not isolate interpolation method alone. A small response
in Cs-high would support, but cannot prove, small responses in other
parameterisations. Use its regional differences relative to experiment
effects to decide which reruns are warranted.

The other three production cases are not reset by this configuration.
No jobs are submitted by this repository change. Gadi compilation and
scientific validation remain to be performed.
