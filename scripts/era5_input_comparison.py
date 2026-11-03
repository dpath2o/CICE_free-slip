from pathlib import Path
from datetime import datetime, timedelta
import numpy as np
from netCDF4 import Dataset, num2date

DAY = datetime(1999, 9, 15)
ROOT = Path("/g/data/gv90/da1339/afim_input/ERA5/0p25")

annual = ROOT / f"era5_for_cice6_{DAY.year}.nc"
monthly = ROOT / "bilinear/monthly_cice6_fast" / (
    f"era5_for_cice6_{DAY.year}_{DAY.month:02d}.nc"
)

fields = [
    "airtmp", "wndewd", "wndnwd", "spchmd",
    "glbrad", "dlwsfc", "ttlpcp",
]

unit_aliases = {
    "m s**-1": "m/s",
    "W m**-2": "W/m^2",
    "kg m**-2 s**-1": "kg/m^2/s",
}

# Zero-based hourly offsets.
ia = (DAY - datetime(DAY.year, 1, 1)).days * 24
im = (DAY - datetime(DAY.year, DAY.month, 1)).days * 24


def as_float(values):
    return np.ma.asarray(values, dtype=np.float64).filled(np.nan)


def check_time(ds, offset, label):
    dim = ds.variables[fields[0]].dimensions[0]
    coord = ds.variables.get(dim)
    units = getattr(coord, "units", "") if coord is not None else ""

    if coord is None or "since" not in units.lower():
        raise ValueError(
            f"{label}: no decodable CF time coordinate; "
            "cannot verify hourly alignment."
        )

    dates = num2date(
        coord[offset:offset + 24],
        units,
        calendar=getattr(coord, "calendar", "standard"),
    )
    actual = [d.strftime("%Y-%m-%d %H:%M:%S") for d in dates]
    expected = [
        (DAY + timedelta(hours=h)).strftime("%Y-%m-%d %H:%M:%S")
        for h in range(24)
    ]

    if actual != expected:
        raise ValueError(
            f"{label}: unexpected timestamps:\n{actual}\n"
            "Investigate time indexing before comparing values."
        )

    print(f"{label}: timestamps verified: {actual[0]} to {actual[-1]}")
    return dim


print("ANNUAL :", annual)
print("MONTHLY:", monthly)
print("DATE   :", DAY.date(), "\n")

with Dataset(annual) as old, Dataset(monthly) as new:
    ta = check_time(old, ia, "Annual")
    tm = check_time(new, im, "Monthly")

    for name in ("LAT", "LON"):
        x = as_float(old.variables[name][:])
        y = as_float(new.variables[name][:])

        if x.shape != y.shape:
            raise ValueError(
                f"{name}: grid shapes differ: {x.shape} / {y.shape}"
            )

        if not np.array_equal(np.isfinite(x), np.isfinite(y)):
            raise ValueError(f"{name}: finite-value masks differ")

        valid = np.isfinite(x) & np.isfinite(y)
        if not valid.any():
            raise ValueError(f"{name}: no finite coordinates")

        delta = y[valid] - x[valid]

        if name == "LON":
            # Shortest angular difference: equivalent longitudes give zero.
            delta = (delta + 180.0) % 360.0 - 180.0

        maximum = float(np.abs(delta).max())
        print(f"{name}: maximum coordinate difference = {maximum:.6e} degrees")

        if maximum > 1e-8:
            raise ValueError(f"{name}: annual and monthly grids differ")

    print("Spatial coordinates match, allowing longitude wrapping")

    print("\nDifferences = monthly minus annual")
    print(
        f"{'field':<9} {'bias':>13} {'RMSE':>13} {'max_abs':>13}"
        f" {'different':>12} {'mask_diff':>11}"
    )

    for name in fields:
        a = old.variables[name]
        b = new.variables[name]

        if a.dimensions[0] != ta or b.dimensions[0] != tm:
            raise ValueError(f"{name}: unexpected time dimension")

        if a.dimensions[1:] != b.dimensions[1:]:
            raise ValueError(
                f"{name}: spatial dimension order differs: "
                f"{a.dimensions} / {b.dimensions}"
            )

        if a.shape[1:] != b.shape[1:]:
            raise ValueError(
                f"{name}: spatial shapes differ: {a.shape} / {b.shape}"
            )

        if a.shape[0] < ia + 24 or b.shape[0] < im + 24:
            raise ValueError(f"{name}: insufficient hourly records")

        ua = getattr(a, "units", "").strip()
        ub = getattr(b, "units", "").strip()

        if not ua or not ub:
            raise ValueError(f"{name}: missing units")

        if unit_aliases.get(ua, ua) != unit_aliases.get(ub, ub):
            raise ValueError(
                f"{name}: incompatible units: {ua!r} / {ub!r}"
            )

        total = squared = maximum = 0.0
        count = different = mask_diff = 0

        for hour in range(24):
            x = as_float(a[ia + hour])
            y = as_float(b[im + hour])

            fx = np.isfinite(x)
            fy = np.isfinite(y)

            mask_diff += int(np.count_nonzero(fx != fy))
            valid = fx & fy
            delta = y[valid] - x[valid]

            if delta.size:
                count += delta.size
                different += int(np.count_nonzero(delta))
                total += float(delta.sum())
                squared += float(np.square(delta).sum())
                maximum = max(maximum, float(np.abs(delta).max()))

        if not count:
            raise ValueError(f"{name}: no comparable finite values")

        print(
            f"{name:<9} {total/count:13.6e}"
            f" {np.sqrt(squared/count):13.6e} {maximum:13.6e}"
            f" {different:12d} {mask_diff:11d}"
        )
        print(
            f"          units annual={ua!r}, monthly={ub!r}; "
            f"dtype annual={a.dtype}, monthly={b.dtype}"
        )

print("\nCompleted: all seven fields compared over 24 hourly records.")
