"""Shared, unit-checked reductions for the two lateral-drag publication notebooks.

This module does not submit jobs, change model inputs, or infer unverified units.
"""
from __future__ import annotations
import hashlib
import json
import re
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr


def signature(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, default=str).encode()).hexdigest()


def file_sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def finite_check(da):
    if not bool(np.isfinite(da).all().compute().item()):
        raise ValueError(f'{da.name}: non-finite values in the comparison support')
    return da


def daily_check(obj, start, end):
    """Require exactly one sample per requested Gregorian day; retain sample times."""
    if 'time' not in obj.coords:
        raise ValueError('A time coordinate is required')
    dates = pd.DatetimeIndex(obj.time.values)
    expected = pd.date_range(start, end, freq='D')
    if dates.has_duplicates or dates.normalize().has_duplicates or not dates.is_monotonic_increasing:
        raise ValueError('Duplicate days or unsorted timestamps')
    if not dates.normalize().equals(expected):
        missing = expected.difference(dates.normalize())
        extra = dates.normalize().difference(expected)
        raise ValueError(f'Daily support mismatch: missing={len(missing)}, extra={len(extra)}')
    return obj


def align_daily(*arrays, start, end):
    for arr in arrays:
        daily_check(arr, start, end)
    # Validate complete coverage before deliberately comparing calendar days.
    normalized = [a.assign_coords(time=pd.DatetimeIndex(a.time.values).normalize()) for a in arrays]
    return xr.align(*normalized, join='exact')


def _unit(da):
    return str(da.attrs.get('units', '')).lower().replace(' ', '').replace('²', '2').replace('³', '3')


def area_to_kkm2(da):
    units = _unit(da)
    scales = {'m2': 1e-9, 'm^2': 1e-9, 'km2': 1e-3, 'km^2': 1e-3,
              '10^3km2': 1., '10^3km^2': 1., '10**3km2': 1., '10**3km^2': 1.,
              '1000km2': 1., '1000km^2': 1.}
    if units not in scales:
        raise ValueError(f'{da.name}: unverified area units {da.attrs.get("units")!r}')
    out = da * scales[units]
    out.attrs = dict(da.attrs, units='10^3 km^2')
    return out


def area_m2(area):
    if _unit(area) not in {'m2', 'm^2'}:
        raise ValueError('Actual grid-cell area in m² is required')
    finite_check(area)
    if bool((area <= 0).any().compute().item()):
        raise ValueError('Cell areas must be positive')
    return area


def binary_mask(mask):
    values = mask.where(mask.notnull(), 0)
    if bool((~values.isin([0, 1])).any().compute().item()):
        raise ValueError('Expected binary occurrence; review fractional FIC semantics first')
    return mask.astype(float)


def daily_observation(obs, start, end, max_gap_days=20):
    """Explicit linear FIA interpolation, no extrapolation or bridging large gaps."""
    times = pd.DatetimeIndex(obs.time.values)
    if times.has_duplicates or not times.is_monotonic_increasing:
        raise ValueError('Observation dates must be unique and sorted')
    target = pd.date_range(start, end)
    if len(times) < 2 or target[0] < times[0] or target[-1] > times[-1]:
        raise ValueError('Observations do not bracket the evaluation period')
    if np.any(np.diff(times) > pd.Timedelta(days=max_gap_days)):
        raise ValueError('Observation gap exceeds the reviewed interpolation limit')
    finite_check(obs)
    out = obs.interp(time=target)
    out.attrs = dict(obs.attrs, temporal_support='linear daily interpolation of native composites')
    return out


def fia_table(model, obs, start, end):
    model, obs = align_daily(area_to_kkm2(model), area_to_kkm2(obs), start=start, end=end)
    finite_check(model); finite_check(obs)
    if bool(((model < 0) | (obs < 0)).any().compute().item()):
        raise ValueError('Negative area: check units, coverage and transcription')
    rows = []
    regions = list(model.region.values) if 'region' in model.dims else [None]
    for region in regions:
        m = model.sel(region=region) if region is not None else model
        o = obs.sel(region=region) if region is not None else obs
        if m.dims != ('time',) or o.dims != ('time',):
            raise ValueError('FIA must have only time and optional region dimensions')
        delta = m-o
        row = dict(region=str(region) if region is not None else 'circum_antarctic', n_days=m.sizes['time'],
                   mean=float(m.mean()), min=float(m.min()), max=float(m.max()), p95=float(m.quantile(.95)),
                   obs_mean=float(o.mean()), bias=float(delta.mean()), mae=float(abs(delta).mean()),
                   rmse=float(np.sqrt((delta**2).mean())), units='10^3 km^2')
        for label, data in [('model', m), ('obs', o)]:
            row[label+'_peak_date'] = str(pd.Timestamp(data.time.values[int(data.argmax('time'))]).date())
            row[label+'_minimum_date'] = str(pd.Timestamp(data.time.values[int(data.argmin('time'))]).date())
        rows.append(row)
    return pd.DataFrame(rows)


def pair_table(target, reference, start, end):
    """Target minus reference on exact daily and regional support, not obs skill."""
    table = fia_table(target, reference, start, end)
    return table.rename(columns={'obs_mean': 'reference_mean', 'bias': 'mean_difference',
                                 'obs_peak_date': 'reference_peak_date', 'obs_minimum_date': 'reference_minimum_date'})


def seasonal_slopes(data, months=(4, 5, 6, 7), label='April_July'):
    rows = []
    for year in sorted(set(data.time.dt.year.values.tolist())):
        dates = pd.DatetimeIndex(data.time.values)
        select = (dates.year == year) & dates.month.isin(months)
        part = data.isel(time=select)
        expected = pd.date_range(f'{year}-{min(months):02d}-01',
                                 pd.Timestamp(year, max(months), 1)+pd.offsets.MonthEnd(0))
        if part.sizes.get('time', 0) != len(expected):
            rows.append(dict(year=year, season=label, complete=False, slope=np.nan)); continue
        finite_check(part)
        slope = float(np.polyfit(np.arange(len(expected)), part.values, 1)[0])
        rows.append(dict(year=year, season=label, complete=True, slope=slope, units='10^3 km^2/day'))
    return pd.DataFrame(rows)


def common_fip(model, obs, start, end):
    model, obs = align_daily(binary_mask(model), binary_mask(obs), start=start, end=end)
    valid = model.notnull() & obs.notnull()
    count = valid.sum('time')
    m = model.where(valid).sum('time', min_count=1) / count.where(count > 0)
    o = obs.where(valid).sum('time', min_count=1) / count.where(count > 0)
    return xr.Dataset({'model': m, 'obs': o, 'diff': m-o, 'valid_days': count})


def spatial_scores(pair, area, region_mask, threshold=.5):
    area = area_m2(area)
    pair, area, region_mask = xr.align(pair, area, binary_mask(region_mask), join='exact')
    valid = (region_mask == 1) & pair.model.notnull() & pair.obs.notnull()
    weights = area.where(valid, 0)
    total = float(weights.sum())
    if total <= 0:
        raise ValueError('No common valid area in region')
    diff = pair.model-pair.obs
    m, o = pair.model >= threshold, pair.obs >= threshold
    hit = float(area.where(valid & m & o, 0).sum())
    miss = float(area.where(valid & ~m & o, 0).sum())
    false = float(area.where(valid & m & ~o, 0).sum())
    return dict(bias=float((diff.fillna(0)*weights).sum()/total),
                mae=float((abs(diff).fillna(0)*weights).sum()/total),
                rmse=float(np.sqrt(((diff**2).fillna(0)*weights).sum()/total)),
                hit_area_m2=hit, omission_area_m2=miss, commission_area_m2=false,
                jaccard=hit/(hit+miss+false) if hit+miss+false else np.nan,
                common_area_m2=total, threshold=threshold)


def crop_region(obj, lon, lat, bounds):
    west, east, south, north = bounds
    lon = (lon+180) % 360-180
    inside_lon = ((lon >= west) & (lon <= east)) if west <= east else ((lon >= west) | (lon <= east))
    inside = inside_lon & (lat >= south) & (lat <= north)
    # Longitude 180 is represented as -180; treat an eastern boundary of 180 consistently.
    if east == 180:
        inside = ((lon >= west) | (lon == -180)) & (lat >= south) & (lat <= north)
    if not bool(inside.any().compute().item()):
        raise ValueError('Region has no cells')
    return obj.where(inside, drop=True)


def fi_fields(field, mask):
    field, mask = xr.align(field, binary_mask(mask), join='exact')
    valid = field.notnull() & mask.notnull()
    occurrence = (mask == 1) & valid
    count = valid.sum('time')
    fi_count = occurrence.sum('time')
    total = field.where(occurrence).sum('time', min_count=1)
    out = xr.Dataset({'hi_conditional': total/fi_count.where(fi_count > 0),
                      'hi_weighted': total.fillna(0)/count.where(count > 0),
                      'fip': fi_count/count.where(count > 0)})
    for name in ['hi_conditional', 'hi_weighted']:
        out[name].attrs = dict(field.attrs)
    return out


def strain_q90(history):
    factors = {'s-1': 1., 's^-1': 1., '1/s': 1., '%/day': .01/86400,
               '%day-1': .01/86400, '1/day': 1/86400, 'day-1': 1/86400}
    arrays = []
    for v in ['divu', 'shear']:
        units = _unit(history[v])
        if units not in factors:
            raise ValueError(f'{v}: verify strain units {units!r}')
        arrays.append(history[v]*factors[units])
    eps = np.hypot(*arrays).where(history.aice > .15)
    if eps.chunks is not None:
        eps = eps.chunk({'time': -1})
    out = eps.quantile(.9, dim='time')*1e6
    out.attrs = {'units': '1e-6 s^-1', 'sampling': 'Q90 of daily-history magnitude, ice-covered days'}
    return out


def regional_means(field, area, domain, selection):
    field, area, domain, selection = xr.align(field, area_m2(area), binary_mask(domain), binary_mask(selection), join='exact')
    valid = field.notnull() & (domain == 1) & (selection == 1)
    w = area.where(valid, 0)
    dims = list(area.dims)
    result = (field.fillna(0)*w).sum(dims) / w.sum(dims).where(w.sum(dims) > 0)
    result.attrs = dict(field.attrs)
    return result


def ice_thickness(history, convention):
    if _unit(history.hi) not in {'m', 'meter', 'metre', 'meters', 'metres'}:
        raise ValueError('Verify hi units in metres')
    if convention == 'ice-covered-area':
        return history.hi
    if convention == 'grid-cell-volume-per-area':
        out = history.hi / history.aice.where(history.aice > 0)
        out.attrs = dict(history.hi.attrs, convention='ice-covered-area thickness derived from hi/aice')
        return out
    raise ValueError('Set the reviewed HI_CONVENTION explicitly')


def save_cached(dataset, path, sig):
    dataset = dataset.copy()
    dataset.attrs = dict(dataset.attrs, analysis_signature=sig)
    dataset.to_netcdf(path, engine='scipy')


def read_cached(path, sig):
    ds = xr.open_dataset(path)
    if ds.attrs.get('analysis_signature') != sig:
        ds.close()
        raise ValueError('Cache provenance changed: recompute the reduction')
    return ds


def namelist_values(path):
    """Read scalar fields used by this audit; this is not a general Fortran parser."""
    result = {}
    for line in Path(path).read_text().splitlines():
        m = re.match(r'^\s*([A-Za-z]\w*)\s*=\s*([^!]+)', line)
        if m:
            result[m[1].lower()] = m[2].strip().rstrip(',').strip("'\"")
    return result


def expected_experiments():
    """Settings checked against the archived 30 October namelist bundle."""
    rows = {}
    def add(name, family, **values):
        rows[name] = dict(experiment=name, family=family, boundary_condition='free_slip',
                          lateral_drag=True, cs=1e-3, cq=0., ktens=.2,
                          e_yieldcurve=1.5, e_plasticpot=1.5, kstrength=0, ndte=360,
                          new_run=False)
        rows[name].update(values)
    add('no-slip-LFI-match', 'controls', boundary_condition='no_slip', lateral_drag=False, new_run=True)
    add('free-slip-no-drag', 'controls', lateral_drag=False, new_run=True)
    add('no-slip-LFI', 'reference', boundary_condition='no_slip', cs=0., e_yieldcurve=1.2, e_plasticpot=1.2)
    add('no-slip-def', 'reference', boundary_condition='no_slip', cs=0., lateral_drag=False, ktens=0., e_yieldcurve=2., e_plasticpot=2.)
    for name, cs in [('Cs-low', 1e-4), ('Cs-mid', 5e-4), ('Cs-high', 1e-3)]: add(name, 'static', cs=cs)
    for name, cq in [('Cq-low', 1), ('Cq-mid', 75), ('Cq-high', 350)]: add(name, 'quadratic', cs=0., cq=cq)
    for name in ['cst-drag', 'GIB-drag']: add(name, 'geometry')
    for name, kt in [('Cs-high-ktens-zero', 0.), ('Cs-high-ktens-low', .05), ('Cs-high-ktens-mid', .1)]:
        add(name, 'tensile', ktens=kt, new_run=name.endswith('zero'))
    add('Cs-high-eDef', 'rheology', e_yieldcurve=2., e_plasticpot=2.)
    add('Cs-high-def-rhe', 'rheology', ktens=0., e_yieldcurve=2., e_plasticpot=2., ndte=720)
    add('Cs-high-roth-def', 'rheology', ktens=0., e_yieldcurve=2., e_plasticpot=2., kstrength=1)
    add('Cs-high-roth-rhe-high', 'rheology', kstrength=1)
    return pd.DataFrame(rows.values())


def ld_face_stress(history):
    """Actual face stress magnitudes, never rho*h times an unverified proxy."""
    for name in ['KuxE', 'KuyE', 'KuxN', 'KuyN']:
        if _unit(history[name]) not in {'pa', 'n/m2', 'n/m^2', 'nm-2', 'nm^-2'}:
            raise ValueError(f'{name}: verify stress definition and N/m² units')
    out = xr.Dataset({'ld_stress_E': np.hypot(history.KuxE, history.KuyE),
                      'ld_stress_N': np.hypot(history.KuxN, history.KuyN)})
    for value in out.data_vars.values():
        value.attrs = {'units': 'Pa', 'definition': 'magnitude at original velocity face; not T-grid vector'}
    return out


def anchor_distance_km(lon, lat, anchor):
    """Great-circle distance to an anchor-bearing cell centre, not coastline distance."""
    from scipy.spatial import cKDTree
    lon, lat, anchor = xr.align(lon, lat, binary_mask(anchor), join='exact')
    good = np.isfinite(lon.values) & np.isfinite(lat.values)
    def xyz(lo, la):
        lo, la = np.deg2rad(lo), np.deg2rad(la)
        return np.column_stack([np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)])
    source = good & (anchor.values == 1)
    if not source.any():
        raise ValueError('No valid anchors')
    tree = cKDTree(xyz(lon.values[source], lat.values[source]))
    chord = tree.query(xyz(lon.values[good], lat.values[good]))[0]
    values = np.full(lon.shape, np.nan)
    values[good] = 2*6371.*np.arcsin(np.clip(chord/2, 0, 1))
    return xr.DataArray(values, coords=lon.coords, dims=lon.dims, name='anchor_cell_distance', attrs={'units': 'km'})


def audit_region_totals(regional, total, tolerance=1e-5):
    regional, total = xr.align(regional, total, join='exact')
    delta = regional.sum('region')-total
    return pd.DataFrame({'time': total.time.values, 'regional_minus_total': delta.values,
                         'within_tolerance': np.abs(delta.values) <= tolerance})
