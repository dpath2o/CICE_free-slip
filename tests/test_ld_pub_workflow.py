from pathlib import Path
import sys
import numpy as np
import pandas as pd
import pytest
import xarray as xr
import nbformat
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'notebooks'))
import ld_pub_workflow as w


def series(values,units='m2'):
    return xr.DataArray(values,dims='time',coords={'time':pd.date_range('2000-04-01',periods=len(values))},attrs={'units':units})


def test_daily_support_and_units():
    da=series([1e9,2e9,0])
    assert w.area_to_kkm2(da).values.tolist()==[1,2,0]
    with pytest.raises(ValueError):w.daily_check(da.isel(time=[0,2]),'2000-04-01','2000-04-03')
    with pytest.raises(ValueError):w.area_to_kkm2(da.assign_attrs(units='unknown'))
    with pytest.raises(ValueError):w.daily_check(xr.concat([da,da],dim='time'),'2000-04-01','2000-04-03')


def test_signed_matched_difference():
    out=w.pair_table(series([1e9,2e9]),series([2e9,4e9]),'2000-04-01','2000-04-02')
    assert out.mean_difference.iloc[0]==-1.5
    assert out.rmse.iloc[0]==pytest.approx(np.sqrt(2.5))


def test_common_support_preserves_zero_and_missing():
    m=xr.DataArray([[1,0,np.nan],[0,0,np.nan],[np.nan,1,np.nan]],dims=['time','x'],coords={'time':pd.date_range('2000-04-01',periods=3)})
    o=xr.DataArray([[0,0,np.nan],[1,np.nan,np.nan],[1,1,np.nan]],dims=m.dims,coords=m.coords)
    out=w.common_fip(m,o,'2000-04-01','2000-04-03')
    assert out.valid_days.values.tolist()==[2,2,0]
    np.testing.assert_allclose(out.model.values[:2],[.5,.5])
    assert np.isnan(out.model.values[2])
    with pytest.raises(ValueError):w.binary_mask(m.fillna(.5))


def test_area_weighted_skill_and_thickness_convention():
    pair=xr.Dataset({'model':('x',[1,0]),'obs':('x',[0,1])})
    area=xr.DataArray([1.,3.],dims='x',attrs={'units':'m2'})
    score=w.spatial_scores(pair,area,xr.ones_like(area))
    assert score['bias']==-.5
    assert score['rmse']==1
    assert score['jaccard']==0
    h=xr.Dataset({'hi':xr.DataArray([1.,0],dims='x',attrs={'units':'m'}),'aice':('x',[.5,0])})
    assert w.ice_thickness(h,'grid-cell-volume-per-area').values[0]==2
    assert np.isnan(w.ice_thickness(h,'grid-cell-volume-per-area').values[1])


def test_stress_units_magnitude_and_dateline():
    h=xr.Dataset({v:xr.DataArray([a],dims='x',attrs={'units':'N/m^2'}) for v,a in [('KuxE',3),('KuyE',4),('KuxN',0),('KuyN',0)]})
    assert float(w.ld_face_stress(h).ld_stress_E.item())==5
    h.KuxE.attrs['units']='s-1'
    with pytest.raises(ValueError):w.ld_face_stress(h)
    da=xr.DataArray([1,2,3],dims='x');lon=xr.DataArray([175,-175,0],dims='x');lat=xr.full_like(lon,-70)
    assert w.crop_region(da,lon,lat,(170,-170,-80,-60)).values.tolist()==[1,2]


def test_cache_provenance_and_experiment_scope(tmp_path):
    p=tmp_path/'cache.nc';w.save_cached(xr.Dataset({'x':('i',[1.])}),p,'one')
    with w.read_cached(p,'one') as d:assert d.x.item()==1
    with pytest.raises(ValueError):w.read_cached(p,'two')
    exp=w.expected_experiments()
    assert len(exp)==19
    assert set(exp.loc[exp.new_run,'experiment'])=={'free-slip-no-drag','no-slip-LFI-match','Cs-high-ktens-zero'}


def test_two_clean_notebooks_compile():
    paths=list((Path(__file__).resolve().parents[1]/'notebooks').glob('*.ipynb'))
    assert len(paths)==2
    for p in paths:
        n=nbformat.read(p,4);nbformat.validate(n)
        for i,c in enumerate(n.cells):
            if c.cell_type=='code':
                compile(c.source,f'{p}:{i}','exec')
                assert c.execution_count is None and not c.outputs
