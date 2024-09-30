# Get observations and reanalyses

import numpy as np
import xarray as xr

def show(function):
    args=inspect.getargspec(function)
    print(args)
    return args

def get_obs(var,data_shortname):
    '''
    Get info of observations to easily open in script
    '''
    #print('var','filename','units','gridlon','gridlat')
    ###
    ##Mean sea level pressure
    ###
    if var in ['psl','msl','prmsl','slp','MSL']:
        if data_shortname=='hadslp2':
            var='slp'
            filename='slp.mon_hadslp2_newtime_185001-201912.nc'
            units='hPa'
            gridlon='0_360'
            gridlat='90_-90'
        elif data_shortname=='era5':
            var='msl'
            filename='mslp_era5_NH_mon_1979-2020.nc'
            units='Pa'
            gridlon='0_360'
            gridlat='90_0'
        elif data_shortname=='era20c':
            var='psl'
            filename='mslp_mon.era20cr_1900-2010.nc'
            units='Pa'
            gridlon='0_360'
            gridlat='90_-90'
        elif data_shortname=='ncep-ncar':
            var='slp'
            filename='slp.mon.mean.ncep-ncar_194801-202203.nc'
            units='mbar'
            gridlon='0_360'
            gridlat='90_-90'
        elif data_shortname=='eraint':
            var='msl'
            filename='mslp.mon.eraint_197901_201512.nc'
            units='Pa'
            gridlon='0_360'
            gridlat='90_-90'
        elif data_shortname=='noaaV2c':
            var='prmsl'
            filename='prmsl_only.mon.mean.noaaV2c_185101-201412_2.0x2.0.nc'
            units='Pa'
            gridlon='0_360'
            gridlat='90_-90'
        elif data_shortname=='noaaV3':
            var='prmsl'
            filename='prmsl_only.mon.mean.noaaV3_185101-201412_1.0x1.0.nc'
            units='Pa'
            gridlon='0_360'
            gridlat='90_-90'
    ###
    ##Near surface temperature
    ###
    elif var in ['tmp','tas','air','t2m']:
        if data_shortname=='cru':
            var='tmp'
            filename='cru_ts4_mon_1901.2020.tmp.dat.nc'
            units='degrees celsius'
            gridlon='-180_180'
            gridlat='-90_90'
        elif data_shortname=='era5':
            var='t2m'
            filename='t2m_era5_NH_mon_1979-2020.nc'
            units='K'
            gridlon='-180_180'
            gridlat='90_0'
        elif data_shortname=='noaaV2c':
            var='air'
            filename='air.2m.mon.mean.noaaV2c_185101-201412_2.0x2.0.nc'
            units='K'
            gridlon='0_360'
            gridlat='90_-90'
    ###
    ##Sea ice concentration
    ###
    elif var in ['sic']:
        if data_shortname=='HadISST':
            var='sic'
            filename='HadISST_ice.nc'
            units='%'
            gridlon='-180_180'
            gridlat='90_-90'
            print('Need to fill values and ds=ds*100 to get %')
            print('ds = ds0[variable].where(ds0[variable] != -1000.)')
    ####
    ##Sea surface temperature
    ###
    elif var in ['sst']:
        if data_shortname=='HadISST':
            var='sst'
            filename='HadISST_sst_v1.1_187001_202105.nc'
            units='degrees celsius'
            gridlon='-180_180'
            gridlat='90_-90'
            print('Need to fill values')
            print('ds = ds0[variable].where(ds0[variable] != -1000.)')
    ###
    ##Precipitation
    ###
    elif var in ['pr','precip','pre']:
        if data_shortname=='aphro':
            var='precip'
            filename='APHRO_mon_MA_025deg_V1101_EXR1.1951-2015.nc'
            units='mm/day'
            gridlon='60 to 150'
            gridlat='-15 to 55'
        elif data_shortname=='gpcc':
            var='precip'
            filename='precip.monthly_v2022_1891_2020_10.nc'
            #filename='precip.mon.total.gpcc.v2018.nc'
            units='mm'
            gridlon='0_360'
            gridlat='90_-90'
            print('month_length = field0dom.time.dt.days_in_month')
            print('field=field0dom/month_length')
        elif data_shortname=='cru':
            var='pre'
            filename='cru_ts4.05.1901.2020.pre.dat.nc'
            units='mm/month'
            gridlon='-180_180'
            gridlat='-90_90'
            print('month_length = field0dom.time.dt.days_in_month')
            print('field=field0dom/month_length')
    ###
    ##Snow cover fraction
    ###
    elif var in ['snc','snowc','snow_cover_extent']:
        if data_shortname=='noaacdr':
            var='snow_cover_extent'
            filename='nhsce_mon.v01r01_19661004_20210503_1.0x1.0.nc'
            units='%'
            gridlon='-180_180'
            gridlat='0_90'
    ##Sensible heat flux
    elif var in ['shtfl']:
        if data_shortname=='noaaV2c':
            var='shtfl'
            filename='shtfl.mon.mean.noaaV2c_185101-201412.nc'
            units='W/m^2 '
            gridlon='0_360'
            gridlat='90_-90'
        if data_shortname=='ncep-ncar':
            var='shtfl'
            filename='shtfl.sfc.mon.mean.ncep-ncar_194801-202210.nc'
            units='W/m^2 '
            gridlon='0_360'
            gridlat='90_-90' 
    ##Latent heat flux
    elif var in ['lhtfl']:
        if data_shortname=='noaaV2c':
            var='lhtfl'
            filename='lhtfl.mon.mean.noaaV2c_185101-201412.nc'
            units='W/m^2 '
            gridlon='0_360'
            gridlat='90_-90'
        if data_shorname=='ncep-ncar':
            var='lhtfl'
            filename='lhtfl.sfc.mon.mean.ncep-ncar_194801-202210.nc'
            units='W/m^2 '
            gridlon='0_360'
            gridlat='90_-90'            
    elif var in ['thf']:
        if data_shortname=='noaaV2c':
            var='thf'      
            filename=['shtfl.mon.mean.noaaV2c_185101-201412.nc','lhtfl.mon.mean.noaaV2c_185101-201412.nc']
            units='W/m^2 '
            gridlon='0_360'
            gridlat='90_-90'
        if data_shortname=='ncep-ncar':
            var='thf'      
            filename=['shtfl.sfc.mon.mean.ncep-ncar_194801-202210.nc','lhtfl.sfc.mon.mean.ncep-ncar_194801-202210.nc']
            units='W/m^2 '
            gridlon='0_360'
            gridlat='90_-90'
    ####
    ##Geopotential height
    ####
    elif var in ['hgt','z','g300']:
        if data_shortname=='noaaV2c':
            print('ds=dsall.sel(level=level')
            var='hgt'
            filename='hgt.mon.mean.noaaV2c_185101-201412.nc'
            units='m'
            gridlon='0_360'
            gridlat='90_-90'
        if data_shortname=='ncep-ncar':
            print('ds=dsall.sel(level=level')
            var='hgt'
            filename='hgtLevs.mon.mean.ncep-ncar_194801-202210.nc'
            units='m'
            gridlon='0_360'
            gridlat='90_-90'
        if data_shortname=='era20c':
            print('Change units: z= G/g')
            var='g300'
            filename='z300_mon.era20cr_1900-2010.nc'
            units='m**2 s**-2'
            gridlon='0_360'
            gridlat='90_-90'
    elif var in ['uwnd']:
        if data_shortname=='noaaV2c':
            var='uwnd'
            filename='uwnd.10m.mon.mean.noaaV2c_185101-201412.nc'
            units='m/s'
            gridlon='0_360'
            gridlat='90_-90'

    info={'var':var,'filename':filename,'units':units,
            'gridlon':gridlon,'gridlat':gridlat}
    print(info)
    return info

