#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Autopep8: https://pypi.org/project/autopep8/
# Check with http://pep8online.com/
##Plotting
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.mpl.ticker as cticker
from cartopy.util import add_cyclic_point
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import *
import matplotlib.colors as mcolors

##Basic libraries
import pandas as pd
from netCDF4 import Dataset
import numpy as np
import xarray as xr
import glob
import math
import numpy.ma as ma


##Statistics
import pymannkendall as mk
from eofs.xarray import Eof
import statsmodels.api as sm
import scipy
from scipy import stats
from scipy.stats import t
from scipy.stats import f
import random

import inspect
import os

def show(function):
    args=inspect.getargspec(function)
    print(args)
    return args

def mon2month(seasonInt):
    monthList=['mon01','mon02','mon03','mon04','mon05','mon06','mon07','mon08','mon09','mon10','mon11','mon12']
    month=monthList[seasonInt-1]
    return month

def latlon(ds):
    lat_str = ''
    lon_str = ''
    other_dims_str = []
    for dim in ds.dims:
        if dim in ['lat', 'latitude']:
            lat_str = dim
        elif dim in ['lon', 'longitude']:
            lon_str = dim
        else:
            other_dims_str.append(dim)
    return lat_str,lon_str

def clim(ds,season='annual',imon=1,iyr=1979,fmon=12,fyr=2005):
    """
        Compute the climatology from monthly data. For seasonal climatology,
        it is possible to shift the start and the end of you
        period in order to select the full winter season (e.g. Dec 1979, Jan 1980, Feb 1980) 
        or select individual winter months
        without removing any data (e.g. for 'DJF' it will take into
        account the first D and last JF).
        The time dimension needs to be named
        'time'.
        Parameters
        ----------
        ds : xarray.core.dataarray.DataArray, xarray.core.dataset.Dataset
            Monthly data to process.
        
        season : int, str, optional
            Season on wchich to compute the climatology. Default is 'annual'.
            Options are:
            - 'annual'
            - single month: int (ex: 1 for January, 2 for February, etc.)
            - any character string (ex: 'DJF', 'JJAS', etc.)
        imon_obs,fmon_obs,iyr_obs,fyr_obs : int,optional
        fyr_obs is set to 2005 in order to avoid erros when data extension is short
        
          Returns
        -------
        clim : xarray.core.dataarray.DataArray, xarray.core.dataset.Dataset
            Weighted climatology.
        -------------------------------------
         Example
        -------
         import xarray as xr
         import sys
         sys.path.insert(1, '/home/maria/Documents/MyPythonLibrary/')
         import climbasis as climb
         da = xr.open_dataarray(...)
         clim = climb.clim(da, season='annual', imon=1,iyr=1979,fmon=12,fyr=2005)
         vals=field0.where(field0['time.month'].isin([12,1,2]).resample(time='AS-Dec')).mean('time') 
    """
    field0=ds.sel(time=slice(str(iyr)+"-"+str(imon), str(fyr)+"-"+str(fmon)))
    month = field0['time.month']
    
    if isinstance(season, int):
        season_sel = (month == season)
        print('Use monthly_selection(ds,mon,iyr,fyr)')
        #print(season_sel)
    elif isinstance(season, str) and len(season) > 1:
        season_str = 'JFMAMJJASONDJFMAMJJASOND'
        #print(season_str.index(season)) ## Index gives the position of the season letter in the season str; then we advance
        #print(len(season))
        month_start = season_str.index(season) + 1
        month_end = month_start + len(season) - 1
        #print(month_start)
        #print(month_end)
        if month_end > 12:
            month_end -= 12  #x-=12 equivalent x=x*12 (multiple)
            season_sel = (month >= month_start) | (month <= month_end)
            #print('Check monthly/seasonal selection')
            #print(season_sel)
        else:
            season_sel = (month >= month_start) & (month <= month_end)
            #print('Check monthly/seasonal selection')
            #print(season_sel)

    else:
        raise ValueError(
                    f"""Invalid season argument: '{season}'. Valid values are:
                - 'annual'"
                - single month: int (ex: 1 for January, 2 for February, etc.)
                - any character string (ex: 'DJF', 'JJAS', etc.)")
                """
                )
    clim=field0.sel(time=season_sel).mean('time')
    #clim.attrs['period'] = str(field0[0]['time.year'].values) + '-' + str(ds[-1]['time.year'].values)
    return clim

def calendar(ds,imon,iyr,fmon,fyr):
    """
    Spatial monthly climatology on a field. The ouput is a 'list' of 12 monthly climatology
    climatology = ds.groupby('time.month').mean('time')
    anomalies = ds.groupby('time.month') - climatology
    """
    field0=ds.sel(time=slice(str(iyr)+"-"+str(imon), str(fyr)+"-"+str(fmon)))
    calendar_clim=field0.groupby('time.month').mean('time')
    calendar_std=field0.groupby('time.month').std('time')
    return calendar_clim,calendar_std

######---------------------------------------------
##
## Monthly and seasonal selection
##
#####-----------------------------------------------

def monthly_selection(ds,season,iyr,fyr):
    #https://stackoverflow.com/questions/60791186/select-xarray-dataset-based-on-month
    # Use .groupby('time.month') to organize the data into months
    # then use .groups to extract the indices for each month
    if isinstance(season, int):
        rmon=season
        #print(season_sel)
    elif isinstance(season, str) and season[0]=='m':
        rmon=int(season.split('mon')[1])
    lat_str,lon_str=latlon(ds)
    print(lat_str,lon_str)
    field_period=ds.sel(time=slice(str(iyr)+"-"+str(1), str(fyr)+"-"+str(12)))
    time_ar=pd.date_range(start=str(iyr)+'-01',end=str(fyr+1)+'-01',freq='Y')
    month_idxs=field_period.groupby('time.month').groups
    # Extract the time indices corresponding to all the Januarys
    selmon_idxs=month_idxs[rmon]
    # Extract the january months by selecting
    valsmon=field_period.isel(time=selmon_idxs)
    anomsmon=valsmon-valsmon.mean('time')
    if len(ds.dims)==3:
        valsmon= xr.DataArray(data=valsmon, dims=["time",ds.dims[1],ds.dims[2]],
            coords=[valsmon.time,valsmon.coords[ds.dims[1]],valsmon.coords[ds.dims[2]]])
        anomsmon= xr.DataArray(data=anomsmon, dims=["time",ds.dims[1],ds.dims[2]],
            coords=[anomsmon.time,anomsmon.coords[ds.dims[1]],anomsmon.coords[ds.dims[2]]])
    elif len(ds.dims)==4:
        valsmon= xr.DataArray(data=valsmon, dims=["time",ds.dims[1],ds.dims[2],ds.dims[3]],
            coords=[valsmon.time,valsmon.coords[ds.dims[1]],valsmon.coords[ds.dims[2]],valsmon.coords[ds.dims[3]]])
        anomsmon= xr.DataArray(data=anomsmon, dims=["time",ds.dims[1],ds.dims[2],ds.dims[3]],
            coords=[anomsmon.time,anomsmon.coords[ds.dims[1]],anomsmon.coords[ds.dims[2]],anomsmon.coords[ds.dims[3]]])
    #time_ar=pd.date_range(start=str(iyr)+'-01',end=str(fyr+1)+'-01',freq='M')
    return valsmon,anomsmon    #time_ar=pd.date_range(start=str(iyr)+'-01',end=str(fyr+1)+'-01',freq='M')


##### select only daily data from June, July, and August
###da_jja_only = ds.sel(time=ds.time.dt.month.isin([6, 7, 8]))
def seasonal_selection(ds,season='MA',iyr=1979,fyr=2005):
    """
        Output of field of seasonal/monthly values 3-D Xarray (nyr,lat,lon). For seasonal ,it is possible to shift the start and the end of you
        period in order to select the full winter season (e.g. Dec 1979, Jan 1980, Feb 1980) or select individual winter months
        without removing any data (e.g. for 'DJF' it will take into
      account the first D and last JF).
        The time dimension needs to be named
        'time'.
        Parameters
        ----------
        ds : xarray.core.dataarray.DataArray, xarray.core.dataset.Dataset
            Monthly data to process.
        
        season : int, str, optional
            Season on wchich to compute the climatology. Default is 'annual'.
            Options are:
            - 'annual'
            - single month: int (ex: 1 for January, 2 for February, etc.)
            - any character string (ex: 'DJF', 'JJAS', etc.)
        imon_obs,fmon_obs,iyr_obs,fyr_obs : int,optional
        fyr_obs is set to 2005 in order to avoid erros when data extension is short
        
          Returns
        -------
     seasonal_values: xarray.core.dataarray.DataArray, xarray.core.dataset.Dataset
            Weighted climatology.
    """
    if season[0]=='D':
        imon=6;fmon=6
    else:
        imon=1;fmon=12
    field0=ds.sel(time=slice(str(iyr)+"-"+str(imon), str(fyr+1)+"-"+str(fmon)))
    #print(field0.time)
    month = field0['time.month']
    if isinstance(season, int):
        season_sel = (month == season)
        #print(season_sel)
        nyr=fyr-iyr+1
        len_season=int(1)
    elif isinstance(season, str) and len(season) > 1:
        season_str = 'JFMAMJJASONDJFMAMJJASOND'
        #print(season_str.index(season)) ## Index gives the position of the season letter in the season str; then we advance
        len_season=len(season) ##I need this for monthly seasonal selection
        month_start = season_str.index(season) + 1
        month_end = month_start + len(season) - 1
            #print(month_start)
            #print(month_end)
        if month_end > 12:
            month_end -= 12  #x-=12 equivalent x=x*12 (multiple)
            season_sel = (month >= month_start) | (month <= month_end)
            #print('Check monthly/seasonal selection: winter')
            #print(season_sel)
            nyr=fyr-iyr+1
        else:
            season_sel = (month >= month_start) & (month <= month_end)
            #print('Check monthly/seasonal selection')
             #print(season_sel)
            nyr=fyr-iyr+1
    else:
        raise ValueError(
                    f"""Invalid season argument: '{season}'. Valid values are:
                - 'annual'"
                - single month: int (ex: 1 for January, 2 for February, etc.)
                - any character string (ex: 'DJF', 'JJAS', etc.)")
                """
                )
    clim=field0.sel(time=season_sel).mean('time')
    #clim.attrs['period'] = str(field0[0]['time.year'].values) + '-' + str(field0[-1]['time.year'].values)
    coords={'time': np.arange(int(iyr),int(fyr)+1,1), field0.dims[1]: field0.coords[field0.dims[1]], field0.dims[2]: field0.coords[field0.dims[2]]}
    zero = xr.DataArray(np.zeros((nyr,field0.shape[1],field0.shape[2])),coords=coords,dims=[field0.dims[0],field0.dims[1], field0.dims[2]])
    seasonal_values=xr.zeros_like(zero)
    seasonal_anomalies=xr.zeros_like(zero)
    if season[0]=='D':
        print('(D-',iyr,' JF-',int(iyr)+1,' to D-',int(fyr),'JF-',int(fyr)+1)
    for i in range(nyr):
        tmp=field0.sel(time=season_sel).values[i*len_season:(i*len_season+len_season),:,:]
        #tmp=field[12*i+(imon-1):12*i+(fmon),:,:]
        seasonal_values[i,:,:]=np.mean(tmp,axis=0)
        seasonal_anomalies[i,:,:]=seasonal_values[i,:,:]-clim[:,:]   
    #seasonal_values.attrs['period'] = str(field0[0]['time.year'].values) + '-' + str(field0[-1]['time.year'].values)
    #seasonal_anomalies.attrs['period'] = str(field0[0]['time.year'].values) + '-' + str(field0[-1]['time.year'].values)
    return seasonal_values,seasonal_anomalies

#############-----------------------------------------------------------
##
## Trends and regression
##
#############------------------------------------------------------------

from scipy import signal
'''
https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.detrend.html
scipy.signal.detrend¶

scipy.signal.detrend(data, axis=- 1, type='linear', bp=0, overwrite_data=False)[source]

    Remove linear trend along axis from data.

    Parameters

        dataarray_like

            The input data.
        axisint, optional

            The axis along which to detrend the data. By default this is the last axis (-1).
        type{‘linear’, ‘constant’}, optional

            The type of detrending. If type == 'linear' (default), the result of a linear least-squares fit to data is subtracted from data. If type == 'constant', only the mean of data is subtracted.
        bparray_like of ints, optional

            A sequence of break points. If given, an individual linear fit is performed for each part of data between two break points. Break points are specified as indices into data. This parameter only has an effect when type == 'linear'.
        overwrite_databool, optional

            If True, perform in place detrending and avoid a copy. Default is False

    Returns

        retndarray

            The detrended input data.


'''
#########----------------
## Detrending
#########----------------
def signal_linear_detrend(seasonal_values,axis=0,type='linear'):
    coords={'time': np.arange(0,seasonal_values.shape[0],1), seasonal_values.dims[1]: seasonal_values.coords[seasonal_values.dims[1]], seasonal_values.dims[2]: seasonal_values.coords[seasonal_values.dims[2]]}
    zero = xr.DataArray(np.zeros((seasonal_values.shape[0],seasonal_values.shape[1],seasonal_values.shape[2])),coords=coords,dims=[seasonal_values.dims[0],seasonal_values.dims[1], seasonal_values.dims[2]])
    detrend=signal.detrend(seasonal_values,axis=0,type='linear')
    detrend_xr=xr.DataArray(data=detrend,coords=coords,dims=[seasonal_values.dims[0],seasonal_values.dims[1], seasonal_values.dims[2]])
    return detrend_xr

def detrend_dim(da, dim, deg=2):
    #https://climate.usu.edu/people/yoshi/pyclm101/monthly.html
    # detrend along a single dimension
    p = da.polyfit(dim=dim, deg=deg)
    fit = xr.polyval(da[dim], p.polyfit_coefficients)
    return da - fit

def detrend_dim_fit(da, dim, deg=2):
    #https://climate.usu.edu/people/yoshi/pyclm101/monthly.html
    # detrend along a single dimension
    p = da.polyfit(dim=dim, deg=deg)
    fit = xr.polyval(da[dim], p.polyfit_coefficients)
    return da - fit,fit

#--------------------------------------------------------------------------------------------------------------------------------
def target_std(ds,rmon,iyr,fyr):
    '''  
    Compute detrended monthly anomalies by removing trend from different period
    -----
    ds: xarray
    rmon: target month
 
    return: 
    field of detrended anomalies
    '''
    f=ds.sel(time=slice(str(iyr)+"-"+str(1), str(fyr)+"-"+str(12)))
    f_month_idxs=f.groupby('time.month').groups
    f_selmon_idxs=f_month_idxs[rmon]
    f_vals=f.isel(time=f_selmon_idxs)
      
    f_anoms=f_vals-f_vals.mean('time') #Anomalies of target period
    f_anoms_std=f_anoms/f_anoms.std('time') #normalised each month as in dyn adj
    
    return f_anoms_std


#####
##Detrend for monthlty anomalies
#####
def target_det(ds,rmon,iyrLong,fyrLong,iyrTarget,fyrTarget,deg):
    '''  
    Compute detrended monthly anomalies by removing trend from different period
    -----
    ds: xarray
    rmon: target month
    iyrLong,fyrLong: long period of ds
    iyrTarget,fyrTarget: shorter period considered to remove clim and trend
    deg:degree for polyfit 
    
    return: 
    field of detrended anomalies
    '''
    fLong=ds.sel(time=slice(str(iyrLong)+"-"+str(1), str(fyrLong)+"-"+str(12)))
    fLong_month_idxs=fLong.groupby('time.month').groups
    fLong_selmon_idxs=fLong_month_idxs[rmon]
    fLong_vals=fLong.isel(time=fLong_selmon_idxs)
    
    fTarget=ds.sel(time=slice(str(iyrTarget)+"-"+str(1), str(fyrTarget)+"-"+str(12)))
    fTarget_month_idxs=fTarget.groupby('time.month').groups
    fTarget_selmon_idxs=fTarget_month_idxs[rmon]
    fTarget_vals=fTarget.isel(time=fTarget_selmon_idxs)
    
    fTarget_anoms=fTarget_vals-fTarget_vals.mean('time') #Anomalies of target period
    fLong_anoms_target=fLong_vals-fTarget_vals.mean('time') #Anomalies Long by removing climatology from target period
    
    p = fTarget_anoms.polyfit(dim='time', deg=deg) #Coeffient of the trend in target period
    fit = xr.polyval(fLong_anoms_target['time'], p.polyfit_coefficients) # fit of long term anoms using p coefficients of target period
    anoms_det=fLong_anoms_target-fit

    return anoms_det

############---------------------
## Spatial trend
############---------------------

def trend_vect(x,y,dim):
    '''
    Compute the spatial trend vectorized, instead of grid by grid.
    ex:
    par=trend_vect(vals.time,vals,'time')

    Source:https://stackoverflow.com/questions/52094320/with-xarray-how-to-parallelize-1d-operations-on-a-multidimensional-dataset
    https://github.com/mickaellalande/MC-Toolkit/blob/master/conda_environment_xarray_xesmf_proplot/xarray/advanced-analysis.ipynb
    
    You can also use a loop on lon/lat but way longer!-> spatial_regression()
    '''
    print('trend-0','intercept-1','rvalue-2','pvalue-3','stderr-4')
    return xr.apply_ufunc(
        stats.linregress, x, y,
        input_core_dims=[[dim], [dim]],
        output_core_dims=[[], [], [], [], []],
        vectorize=True
    )

def trend_theilsen_vect(y,x,ci):
    '''
    Compute the spatial trend vectorized, instead of grid by grid.
    ex:
    par=trend_vect(vals.time,vals,'time')

    Source:https://stackoverflow.com/questions/52094320/with-xarray-how-to-parallelize-1d-operations-on-a-multidimensional-dataset
    https://github.com/mickaellalande/MC-Toolkit/blob/master/conda_environment_xarray_xesmf_proplot/xarray/advanced-analysis.ipynb
    
    scipy.stats.mstats.theilslopes(y, x=None, alpha=0.95)[source]
    Computes the Theil-Sen estimator for a set of points (x, y).

    theilslopes implements a method for robust linear regression. It computes the slope as the median of all slopes between paired values.

     Parameters:
     y : array_like Dependent variable.
     x : {None, array_like}, optional Independent variable. If None, use arange(len(y)) instead.
     alpha : float Confidence degree between 0 and 1. Default is 95% confidence. Note that alpha is symmetric around 0.5, i.e. both 0.1 and 0.9 are interpreted as “find the 90% confidence interval”.

     Returns:
     medslope : float Theil slope. 
     medintercept : float Intercept of the Theil line, as median(y) - medslope*median(x).
     lo_slope : float Lower bound of the confidence interval on medslope.
     up_slope : float Upper bound of the confidence interval on medslope.
    '''
    print('medslope-0','medintercept-1','loslope-2','upslope-3')
    return xr.apply_ufunc(
        stats.theilslopes,y, x, 0.95,
        input_core_dims=[[dim], [dim]],
        output_core_dims=[[], [], [], [], []],
        vectorize=True
    )

###########-----------------------
##Correlation/regression lat,lon (longer time)
###########-----------------------

def spatial_correlation_2fields(field1,field2):
    '''
    Compute the spatial correlation of 2 fields.
    '''
    coords2D={field1.dims[1]: field1.coords[field1.dims[1]], field1.dims[2]: field1.coords[field1.dims[2]]}
    zero2D= xr.DataArray(np.zeros((field1.shape[1],field1.shape[2])),coords=coords2D,dims=[field1.dims[1], field1.dims[2]])
    trend=xr.zeros_like(zero2D)
    intercept=xr.zeros_like(zero2D)
    rvalue=xr.zeros_like(zero2D)
    pvalue=xr.zeros_like(zero2D)
    stderr=xr.zeros_like(zero2D)
    parspatial=[]
    for j in range(field1.shape[1]):
        for i in range(field1.shape[2]):
            par=stats.linregress(field1[:,j,i],field2[:,j,i])
            trend[j,i]=par[0]
            intercept[j,i]=par[1]
            rvalue[j,i]=par[2]
            pvalue[j,i]=par[3]
            stderr[j,i]=par[4]
    parspatial=[trend,intercept,rvalue,pvalue,stderr]
    print('0-trend','1-intercept','2-rvalue','3-pvalue','4-stderr')
    return parspatial


def spatial_regression(seasonal_values,index):
    '''
    Compute the spatial trend from seasonal field. To use after applying seasonal_selection
    '''
    if seasonal_values.shape[0]==len(index):
        vals=seasonal_values.values
        years=index
        coords2D={seasonal_values.dims[1]: seasonal_values.coords[seasonal_values.dims[1]], seasonal_values.dims[2]: seasonal_values.coords[seasonal_values.dims[2]]}
        zero2D= xr.DataArray(np.zeros((seasonal_values.shape[1],seasonal_values.shape[2])),coords=coords2D,dims=[seasonal_values.dims[1], seasonal_values.dims[2]])
        trend=xr.zeros_like(zero2D)
        intercept=xr.zeros_like(zero2D)
        rvalue=xr.zeros_like(zero2D)
        pvalue=xr.zeros_like(zero2D)
        stderr=xr.zeros_like(zero2D)
        parspatial=[]
        for j in range(seasonal_values.shape[1]):
            for i in range(seasonal_values.shape[2]):
                yd=vals[:,j,i]
                par=stats.linregress(years,yd)
                trend[j,i]=par[0]
                intercept[j,i]=par[1]
                rvalue[j,i]=par[2]
                pvalue[j,i]=par[3]
                stderr[j,i]=par[4]
        parspatial=[trend,intercept,rvalue,pvalue,stderr]
        print('0-trend','1-intercept','2-rvalue','3-pvalue','4-stderr')
    else:
        raise ValueError(
                    f"""Check dimensions'"
                """
                )
    return parspatial


def spatial_regression_nans(seasonal_values):
    '''
    Compute the spatial trend from seasonal field. To use after applying seasonal_selection
    '''
    zero2D=np.empty((seasonal_values.shape[1],seasonal_values.shape[2]))
    trend=np.empty((seasonal_values.shape[1],seasonal_values.shape[2]))
    intercept=np.empty((seasonal_values.shape[1],seasonal_values.shape[2]))
    rvalue=np.empty((seasonal_values.shape[1],seasonal_values.shape[2]))
    pvalue=np.empty((seasonal_values.shape[1],seasonal_values.shape[2]))
    stderr=np.empty((seasonal_values.shape[1],seasonal_values.shape[2]))
    parspatial=[]
    for j in range(seasonal_values.shape[1]):
        for i in range(seasonal_values.shape[2]):
            Y=seasonal_values[:,j,i]
            Y=Y[np.logical_not(np.isnan(Y))]
            par=stats.linregress(np.arange(0,Y.shape[0],1),Y)
            trend[j,i]=par[0]
            intercept[j,i]=par[1]
            rvalue[j,i]=par[2]
            pvalue[j,i]=par[3]
            stderr[j,i]=par[4]
    parspatial=[trend,intercept,rvalue,pvalue,stderr]
    print('0 trend','1 intercept','2 rvalue','3 pvalue','4 stderr')
    return parspatial

#############--------------------------------
##
##Spatial average / time series / annual cycle
##
############--------------------------------

def spatial_average(ds):
    # Get the longitude and latitude names + other dimensions to test that the
    # sum of weights is right
    lat_str = ''
    lon_str = ''
    other_dims_str = []
    for dim in ds.dims:
        if dim in ['lat', 'latitude']:
            lat_str = dim
        elif dim in ['lon', 'longitude']:
            lon_str = dim
        else:
            other_dims_str.append(dim)

    # Compute the weights
    coslat = np.cos(np.deg2rad(ds[lat_str])).where(~ds.isnull())
    weights = coslat / coslat.sum(dim=(lat_str, lon_str))

    # Test that the sum of weights equal 1
    np.testing.assert_allclose(
        weights.sum(dim=(lat_str, lon_str)).values,
        np.ones([ds.coords[dim_str].size for dim_str in other_dims_str]),
        rtol=1e-06
    )

    with xr.set_options(keep_attrs=True):
        return (ds * weights).sum(dim=(lat_str, lon_str))
#
#
def annual_cycle(ds_spaavg,imon,iyr,fmon,fyr):
    '''
    Compute the annual cycle of a field.
    First, select starting and ending months and period. This option is set in case to start snow annual cycle from october.
    Second,compute spatial average of field domain weighting latitude of the field domain.
    Finally, group by month and do the climatology.
    '''
   
    field=ds_spaavg.sel(time=slice(str(iyr)+"-"+str(imon), str(fyr)+"-"+str(fmon)))
    annual_cycle=field.groupby('time.month').mean('time')
    return annual_cycle
#
#
def time_series(seasonal_values,field):
    '''
    Compute the time series from seasonal field,i.e., spatial average of seasonal values. 
    Need to call field in order to keep dimensions of field domain (This can be improved and take it from seasonal_values)
    To use after applying seasonal_selection
    '''
    tm=seasonal_values.mean((field.dims[1],field.dims[2]))
    return tm

########################----------------------------------------------------------------
##
##Pattern correlation and RMS difference
##
#https://pcmdi.llnl.gov/staff/taylor/CV/Taylor_diagram_primer.pdf
##https://www.ncl.ucar.edu/Document/Functions/Contributed/pattern_cor.shtml
#https://archive.ipcc.ch/ipccreports/tar/wg1/458.htm
########################----------------------------------------------------------------

def mean_aux(x):
    y = np.nansum(x) / np.size(x);
    return y

def corr2D(a,b):
    a = a - mean_aux(a)
    b = b - mean_aux(b)
    r = np.nansum((a*b))/ math.sqrt(np.nansum((a*a))* np.nansum((b*b)));
    return r

def rms(a,b):
    a=a - mean_aux(a)
    b = b - mean_aux(b)
    rms=math.sqrt(np.nansum((a-b)**2))
    return rms

###################---------------------------
##
## Time lagged cross correlation
##
###################---------------------------

def crosscorr(datax, datay, lag=0):
    """ Lag-N cross correlation. 
    Shifted data filled with NaNs 

    Parameters
    ----------
    lag : int, default 0
    datax, datay : pandas.Series objects of equal length
    Returns
    ----------
    crosscorr : float
    """
    return datax.corr(datay.shift(lag))

####################------------------------------
##
##Running mean
##
####################-------------------------------

def running_mean(x, N):
    """
    ##most efficiennt in time https://stackoverflow.com/questions/13728392/moving-average-or-running-mean
    x: np.array
    N: window
    """
    cumsum = np.cumsum(np.insert(x, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / float(N)


######################---------------------------------------------
##
##EOF analysis
##
#https://ajdawson.github.io/eofs/latest/api/eofs.xarray.html
#https://atmos.washington.edu/~dennis/552_Notes_4.pdf
######################---------------------------------------------

def eofanalysis(anoms,modes=5):
    '''
    Better to use directly:
    coslat = np.cos(np.deg2rad(anomsxO.coords[latxO].values)) ##ydom are Xarray so you have to especify .values
solver = Eof(anomsxO) #, weights=wgts)
    eofs= solver.eofs(neofs=5)
    eofsCov = solver.eofsAsCovariance(neofs=5)
    eofsCor = solver.eofsAsCorrelation(neofs=5)
    pcs = solver.pcs(npcs=5, pcscaling=1)
    fvars= solver.varianceFraction(neigs=5)

    '''
    lat,lon=climb.latlon(anoms)
    coslat = np.cos(np.deg2rad(anoms.coords[lat].values)) ##ydom are Xarray so you have to especify .values
    solver = Eof(anoms) #, weights=wgts)
    eofs= solver.eofs(neofs=modes)
    eofsCov = solver.eofsAsCovariance(neofs=modes)
    eofsCor = solver.eofsAsCorrelation(neofs=modes)
    pcs = solver.pcs(npcs=modes, pcscaling=1)
    fvars= solver.varianceFraction(neigs=modes)
    return solver,eofs,eofsCov,eofsCor,pcs,fvar


##################---------------------------------
##
## Angle between two lines
##
##https://stackoverflow.com/questions/28260962/calculating-angles-between-line-segments-python-with-math-atan2
##slope is for  the sense of the intersection (clockwise positive using a line as reference)
##################---------------------------------

def slope(x1, y1, x2, y2):
    return (y2-y1)/(x2-x1)

def dot(vA, vB):
    return vA[0]*vB[0]+vA[1]*vB[1]

def ang(lineA, lineB):
    # Get nicer vector form
    vA = [(lineA[0][0]-lineA[1][0]), (lineA[0][1]-lineA[1][1])]
    vB = [(lineB[0][0]-lineB[1][0]), (lineB[0][1]-lineB[1][1])]
    # Get dot prod
    dot_prod = dot(vA, vB)
    # Get magnitudes
    magA = dot(vA, vA)**0.5
    magB = dot(vB, vB)**0.5
    # Get cosine value
    cos_ = dot_prod/magA/magB
    # Get angle in radians and then convert to degrees
    angle = math.acos(dot_prod/magB/magA)
    # Basically doing angle <- angle mod 360
    ang_deg = math.degrees(angle)%360

    if ang_deg-180>=0:
        # As in if statement
        return 360 - ang_deg
    else:

        return ang_deg
#################################----------------------------------------
##
## Derivatives
##
#################################-----------------------------------------
def first_derivative(f, axis=None, x=None, delta=None):
    """
    Source: https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.first_derivative.html?highlight=first_derivative
    Calculate the first derivative of a grid of values.

    Works for both regularly-spaced data and grids with varying spacing.

    Either `x` or `delta` must be specified, or `f` must be given as an `xarray.DataArray` with
    attached coordinate and projection information. If `f` is an `xarray.DataArray`, and `x` or
    `delta` are given, `f` will be converted to a `pint.Quantity` and the derivative returned
    as a `pint.Quantity`, otherwise, if neither `x` nor `delta` are given, the attached
    coordinate information belonging to `axis` will be used and the derivative will be returned
    as an `xarray.DataArray`.

    This uses 3 points to calculate the derivative, using forward or backward at the edges of
    the grid as appropriate, and centered elsewhere. The irregular spacing is handled
    explicitly, using the formulation as specified by [Bowen2005]_.

    Parameters
    ----------
    f : array-like
        Array of values of which to calculate the derivative
    axis : int or str, optional
        The array axis along which to take the derivative. If `f` is ndarray-like, must be an
        integer. If `f` is a `DataArray`, can be a string (referring to either the coordinate
        dimension name or the axis type) or integer (referring to axis number), unless using
        implicit conversion to `pint.Quantity`, in which case it must be an integer. Defaults
        to 0. For reference, the current standard axis types are 'time', 'vertical', 'y', and
        'x'.
    x : array-like, optional
        The coordinate values corresponding to the grid points in `f`
    delta : array-like, optional
        Spacing between the grid points in `f`. Should be one item less than the size
        of `f` along `axis`.

    Returns
    -------
    array-like
        The first derivative calculated along the selected axis


    .. versionchanged:: 1.0
       Changed signature from ``(f, **kwargs)``

    See Also
    --------
    second_derivative

    """
    n, axis, delta = _process_deriv_args(f, axis, x, delta)
    take = make_take(n, axis)

    # First handle centered case
    slice0 = take(slice(None, -2))
    slice1 = take(slice(1, -1))
    slice2 = take(slice(2, None))
    delta_slice0 = take(slice(None, -1))
    delta_slice1 = take(slice(1, None))

    combined_delta = delta[delta_slice0] + delta[delta_slice1]
    delta_diff = delta[delta_slice1] - delta[delta_slice0]
    center = (- delta[delta_slice1] / (combined_delta * delta[delta_slice0]) * f[slice0]
              + delta_diff / (delta[delta_slice0] * delta[delta_slice1]) * f[slice1]
              + delta[delta_slice0] / (combined_delta * delta[delta_slice1]) * f[slice2])

    # Fill in "left" edge with forward difference
    slice0 = take(slice(None, 1))
    slice1 = take(slice(1, 2))
    slice2 = take(slice(2, 3))
    delta_slice0 = take(slice(None, 1))
    delta_slice1 = take(slice(1, 2))

    combined_delta = delta[delta_slice0] + delta[delta_slice1]
    big_delta = combined_delta + delta[delta_slice0]
    left = (- big_delta / (combined_delta * delta[delta_slice0]) * f[slice0]
            + combined_delta / (delta[delta_slice0] * delta[delta_slice1]) * f[slice1]
            - delta[delta_slice0] / (combined_delta * delta[delta_slice1]) * f[slice2])

    # Now the "right" edge with backward difference
    slice0 = take(slice(-3, -2))
    slice1 = take(slice(-2, -1))
    slice2 = take(slice(-1, None))
    delta_slice0 = take(slice(-2, -1))
    delta_slice1 = take(slice(-1, None))

    combined_delta = delta[delta_slice0] + delta[delta_slice1]
    big_delta = combined_delta + delta[delta_slice1]
    right = (delta[delta_slice1] / (combined_delta * delta[delta_slice0]) * f[slice0]
             - combined_delta / (delta[delta_slice0] * delta[delta_slice1]) * f[slice1]
             + big_delta / (combined_delta * delta[delta_slice1]) * f[slice2])

    return concatenate((left, center, right), axis=axis)

