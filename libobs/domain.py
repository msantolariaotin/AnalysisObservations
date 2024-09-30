import xarray as xr
import inspect
##test local message
def show(function):
    args=inspect.getargspec(function)
    print(args)
    return args

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
##modi
##
#Shifing ds from 0 to 360 -> -180 to 180
#https://stackoverflow.com/questions/53345442/about-changing-longitude-array-from-0-360-to-180-to-180-with-python-xarray
##Maybe faster: cdo sellonlatbox,-180,180,-90,90 ifile ofile
###
def shifting_grid(ds):
    lat_name,lon_name=latlon(ds)
    # Adjust lon values to make sure they are within (-180, 180)
    ds['_longitude_adjusted'] = xr.where(
        ds[lon_name] > 180,
        ds[lon_name] - 360,
        ds[lon_name])
    # reassign the new coords to as the main lon coords
    # and sort DataArray using new coordinate values
    ds = (ds.swap_dims({lon_name: '_longitude_adjusted'})
    .sel(**{'_longitude_adjusted': sorted(ds._longitude_adjusted)})
    .drop(lon_name))
    ds = ds.rename({'_longitude_adjusted': lon_name})
    return ds

def shifting_grid_ALLOW(ds):
    lat_name,lon_name=latlon(ds)
    # Adjust lon values to make sure they are within (-180, 180)
    ds['_longitude_adjusted'] = xr.where(
        ds[lon_name] > 90,
        ds[lon_name] - (360),
        ds[lon_name])
    # reassign the new coords to as the main lon coords
    # and sort DataArray using new coordinate values
    ds = (ds.swap_dims({lon_name: '_longitude_adjusted'})
    .sel(**{'_longitude_adjusted': sorted(ds._longitude_adjusted)})
    .drop(lon_name))
    ds = ds.rename({'_longitude_adjusted': lon_name})
    return ds

def shifting_grid_ALLOWback(ds):
    lat_name,lon_name=latlon(ds)
    # Adjust lon values to make sure they are within (-180, 180)
    ds['_longitude_adjusted'] = xr.where(
        ds[lon_name] < -180,
        ds[lon_name] + (360),
        ds[lon_name])
    # reassign the new coords to as the main lon coords
    # and sort DataArray using new coordinate values
    ds = (ds.swap_dims({lon_name: '_longitude_adjusted'})
    .sel(**{'_longitude_adjusted': sorted(ds._longitude_adjusted)})
    .drop(lon_name))
    ds = ds.rename({'_longitude_adjusted': lon_name})
    return ds

def reversing_lat(ds):
    lat_str = ''
    lon_str = ''
    other_dims_str = []
    for dim in ds.dims:
        if dim in ['lat', 'lon']:
            dsO= ds.reindex(lat=list(reversed(ds.lat)))
        elif dim in ['latitude', 'longitude']:
            dsO= ds.reindex(latitude=list(reversed(ds.latitude)))
        else:
            other_dims_str.append(dim)
    return dsO
##ALLOW lon[130E,110W];lat[20,90]
#latS,latN,lonW,lonE=20,90,-230,-110
#field = ds90.where((latS < ds90.coords[lat]) & ( ds90.coords[lat] < latN) & (lonW < ds90.coords[lon]) & ( ds90.coords[lon] < lonE),drop=True)
# =============================================================================
# Zones
# =============================================================================


def coord_domain(domain):
    """
        Get the latitude and longitude limits for the corresponding zone.
        Parameters
        ----------
        zone : str
            Zone name. Options are:
            - 'GLOB', 'global', 'GLOBAL', 'GLOB-land'
            - 'NH' : North Hemisphere
            - 'HMA' : High Mountain of Asia
            - 'NA' : North America
        Returns
        -------
        latlim, lonlim : slice
            Latitude and longitude limits of the zone.
        -------
    """
    if domain in ['global','GLOB']:
        latS=-90;latN=90; lonW=-180;lonE=180
        latlim=slice(None);lonlim=slice(None)
    # Arctic
    elif domain in ['ARC','Arctic']:
        latS = 50;latN = 90;lonW=-179;lonE=180
        latlim=slice(90,50);lonlim=slice(None)
    # North Hemisphere
    elif domain in ['NH']:
        latS = 0;latN = 90;lonW=-180;lonE=180
        latlim=slice(90,0);lonlim=slice(None)
    elif domain in ['NHpolar']:
        latS = 28;latN = 90;lonW=-180;lonE=180
        latlim=slice(90,28);lonlim=slice(None)
        # North Hemisphere
    elif domain in ['NHeq']:
        latS = 0;latN = 90;lonW=-179;lonE=180
        latlim=slice(90,0);lonlim=slice(None)
    elif domain in ['SH']:
        latS=-90;latN=-20;lonW=-180;lonE=180
        latlim=slice(-90,-20);lonlim=slice(None)
    elif domain in ['SHWest']:
        latS=-90;latN=-20;lonW=-170;lonE=-10
        latlim=slice(-90,-20);lonlim=slice(-170,-10)
    elif domain in ['NINO34']:
        latS=-5;latN=5;lonW=-170;lonE=-120
        latlim=slice(latS,latN);lonlim=slice(None)
    elif domain in ['EQPAC']:
        latS=-15;latN=15;lonW=-170;lonE=-100
        latlim=slice(latS,latN);lonlim=slice(None)
    elif domain in ['SRP']:
        latS = 20;latN = 60;lonW=30;lonE=130
        latlim=slice(latS,latN);lonlim=slice(lonW,lonE)
    elif domain in ['MA']:
        latS = 10;latN = 60;lonW=40;lonE=130
        latlim=slice(latS,latN);lonlim=slice(lonW,lonE)
    elif domain in ['MA20W']:
        latS = 10;latN = 60;lonW=40-20;lonE=130
        latlim=slice(latS,latN);lonlim=slice(lonW,lonE)
    elif domain in ['MA20N']:
        latS = 10;latN = 60+20;lonW=40;lonE=130
        latlim=slice(latS,latN);lonlim=slice(lonW,lonE)
    elif domain in ['MA20E']:
        latS = 10;latN = 60;lonW=40;lonE=130+20
        latlim=slice(latS,latN);lonlim=slice(lonW,lonE)
    elif domain in ['MA20S']:
        latS = 10-20;latN = 60;lonW=40;lonE=130
        latlim=slice(latS,latN);lonlim=slice(lonW,lonE)
    elif domain in ['MA20']:
        latS = 10-20;latN = 60+20;lonW=40-20;lonE=130+20
        latlim=slice(latS,latN);lonlim=slice(lonW,lonE)
    # High Mountain of Asia (HMA)
    elif domain in ['HMA']:
        latS = 20;latN = 45;lonW=60;lonE=110
        latlim=slice(latS,latN);lonlim=slice(lonW,lonE)
 # HK: Hindu-Kush / Karakoram / Western Himalay
# CH: Central and Est Himalaya
# TB: Tibetan Plateau
    elif domain in ['HK']:
        latS = 31;latN = 40;lonW=70;lonE=81
        lonlim = slice(70, 81);latlim = slice(31, 40)
    elif domain in ['CH']:
        latS = 26;latN = 31;lonW=79;lonE=98
        lonlim= slice(79, 98);latlim= slice(26, 31)
    elif domain in ['TP']:
        latS = 31;latN = 39;lonW=81;lonE=104
        lonlim= slice(81, 104);latlim= slice(31, 39)
    # North America
    elif domain in ['NAm']:
        latS = 20;latN = 90;lonW=-175;lonE=-2.5
        latlim=slice(latS,latN);lonlim=slice(lonW,lonE)
    #Europe
    elif domain in ['NHEast']:
        latS = 0;latN = 90;lonW=-20;lonE=180
        latlim=slice(latN,latS);lonlim=slice(lonW,lonE)
    # Eurasia
    elif domain in ['EURA']:
        latS = 50;latN = 90;lonW=0;lonE=180
        latlim=slice(90,50);lonlim=slice(0,180)
        # Eurasia
    elif domain in ['EURAsouth']:
        latS = 40;latN = 90;lonW=0;lonE=180
        latlim=slice(90,40);lonlim=slice(0,180)
        # EURO-ASIA 
    elif domain in ['EAA']:
        latS = 10;latN = 70;lonW=-40;lonE=90
        latlim=slice(70,10);lonlim=slice(-40,90)
    # Scandinavian Pattern
    elif domain in ['SCAND']:
        latS = 60;latN = 90;lonW=0;lonE=150
        latlim=slice(90,60);lonlim=slice(0,150)
    elif domain in ['SCAND20120']:
        latS = 60;latN = 90;lonW=20;lonE=120
        latlim=slice(latN,latS);lonlim=slice(lonW,lonE)
    elif domain in ['SCANDCONT']:
        latS = 50;latN = 70;lonW=40;lonE=120
        latlim=slice(latN,latS);lonlim=slice(lonW,lonE)
       # (30 °W–150 °E, 20° –80 °N)Bueh and Nakamura
    elif domain in ['NAEext']:
        latS=20;latN=90;lonW=-90;lonE=40
        latlim=slice(latN,latS);lonlim=slice(lonW,lonE)
    elif domain in ['NAElim']:
        latS=30;latN=80;lonW=-90;lonE=40
        latlim=slice(latN,latS);lonlim=slice(-90,40)
    elif domain in ['NAEISO']:
        latS=20;latN=90;lonW=-50;lonE=40
        latlim=slice(latN,latS);lonlim=slice(lonW,lonE)
    elif domain in ['NAEz200']:
        latS=5;latN=90;lonW=-50;lonE=40
        latlim=slice(90,5);lonlim=slice(-50,40)
    elif domain in ['NAE']:
        latS=30;latN=90;lonW=-90;lonE=40
        latlim=slice(90,30);lonlim=slice(-90,40)
    elif domain in ['ALLOW']:
        latS=20;latN=80;lonW=-230;lonE=-110
        latlim=slice(latN,latS);lonlim=slice(lonW,lonE)
    elif domain in ['TROPAC']: ##grid de 0 a 360
        latS=-15;latN=15;lonW=260;lonE=100 #100E,100W
        latlim=slice(latN,latS);lonlim=slice(lonW,lonE)
    elif domain in ['PIR']:
        lonW=-3.;lonE=4;latS=41;latN=44
        latlim=slice(latS,latN);lonlim=slice(lonW,lonE)
    elif domain in ['SibHHas']:
        lonW=80;lonE=120;latS=40;latN=65
        latlim=slice(latN,latS);lonlim=slice(lonW,lonE)
    elif domain in ['SibH']:
        lonW=80;lonE=120;latS=30;latN=60
        latlim=slice(latN,latS);lonlim=slice(lonW,lonE)
    elif domain in ['SubPolar']:
        lonW=-62; lonE=-40; latS=50;latN=58
        latlim=slice(50,58)
        lonlim=slice(-62,-40)
    elif domain in ['SubTropic']:
        lonW=-28; lonE=-16; latS=15;latN=32
        latlim=slice(latS,latN)
        lonlim=slice(lonW,lonE)
    elif domain in ['SubCenTropic']:
        latS=16;latN=25;lonW=-55;lonE=-35
        latlim=slice(latS,latN)
        lonlim=slice(lonW,lonE)
    elif domain in ['MED','Med']:
        latS=32;latN=44;lonW=-10;lonE=20
        latlim=slice(latS,latN)
        lonlim=slice(lonW,lonE)
    elif domain in ['SNA']:
        latS=15;latN=25;lonW=20;lonE=30
        latlim=slice(latS,latN)
        lonlim=slice(lonW,lonE)
    elif domain in ['wG']:
        latS=45;latN=80;lonW=-80;lonE=-50
        latlim=slice(latS,latN)
        lonlim=slice(lonW,lonE)
    elif domain in ['BK']:
        latS=50;latN=90;lonW=-30;lonE=120
        latlim=slice(latN,latS)
        lonlim=slice(lonW,lonE)
    elif domain in ['GINBK']:
        latS=-2.5;latN=2.5;lonW=None;lonE=None
        latlim=slice(latS,latN)
        lonlim=slice(lonW,lonE)
    elif domain in ['SWJ']:
        latS=35;latN=45;lonW=40;lonE=120
        latlim=slice(latS,latN)
        lonlim=slice(lonW,lonE)
    elif domain in ['IPOS']:
        latS=-10;latN=10;lonW=-90;lonE=170
        latlim=slice(latS,latN)
        lonlim=slice(lonW,lonE)
    elif domain in ['IPON']:
        latS=25;latN=45;lonW=-150;lonE=150
        latlim=slice(latS,latN)
        lonlim=slice(lonW,lonE)
    else:
        raise ValueError(
            f"""Invalid zone argument: '{zone}'. Valid zones are:
                - 'GLOB', 'global', 'GLOBAL'
                - 'NH' : North Hemisphere
                - 'HMA' : High Mountain of Asia
                - 'NA' : North America
             """
        )
    return latS,latN,lonW,lonE,latlim,lonlim

def field_dom(ds,domain):
    latS,latN,lonW,lonE,latlim,lonlim=coord_domain(domain)
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
    field = ds.where((latS < ds.coords[lat_str]) & ( ds.coords[lat_str] < latN) & (lonW < ds.coords[lon_str]) & ( ds.coords[lon_str] < lonE),drop=True)
    #print('Domain;latS,latN,lonW,lonE:',domain,latS,latN,lonW,lonE)
    return field

def field_dom_noaaCDR(ds,domain):
    latS,latN,lonW,lonE,latlim,lonlim=coord_domain(domain)
    field=ds.sel(lat=slice(latS,latN),lon=slice(lonW,lonE))
    return field

def field_sel_dom(ds,domain):
        # Get the longitude and latitude names + other dimensions to test that the
    latS,latN,lonW,lonE,latlim,lonlim=coord_domain(domain)
    lat_str = ''
    lon_str = ''
    other_dims_str = []
    for dim in ds.dims:
        if dim in ['lat', 'lon']:
            field=ds.sel(lat=latlim,lon=lonlim)
        elif dim in ['latitude', 'longitude']:
            field=ds.sel(latitude=latlim,longitude=lonlim)
        else:
            other_dims_str.append(dim)
    print(domain,latS,latN,lonW,lonE,latlim,lonlim)
    return field

def field_sel_dom_lat_inverse(ds,domain):
    latS,latN,lonW,lonE,latlim,lonlim=coord_domain(domain)
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
    latS,latN,lonW,lonE,latlim,lonlim=coord_domain(domain)
    #print('Domain;Domain;  latS,latN,lonW,lonE:',domain,latS,latN,lonW,lonE)
    field=ds.sel(lat_str=slice(latN,latS),lon_str=slice(lonW,lonE))
    return field

