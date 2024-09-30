import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.mpl.ticker as cticker
from cartopy.util import add_cyclic_point
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import *
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib.ticker as mticker
import matplotlib.lines as mlines
from matplotlib.legend_handler import HandlerTuple
from fpdf import FPDF

from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from shapely import geometry
from collections import namedtuple
from shapely.geometry.polygon import LinearRing


import sys
sys.path.insert(1, '/home/maria/Documents/MyPythonLibrary/')
import climbasis as climb
import domain as dom

import inspect
####
###fieldY_seasons.plot(x="lon", y="lat", col="season", col_wrap=4,cmap='RdBu')
####
def show(function):
    args=inspect.getargspec(function)
    print(args)
    return args

def zerowhite(interval):
    '''
    CS1=axs.contourf(lons,lats,eof[:,:],levs,
            transform=ccrs.PlateCarree(),
            cmap=cmap,extend='both')
    cbar=fig.colorbar(CS1, cax=cbar_ax,ticks=levs_ticks,orientation='horizontal',label='hPa')
    '''
    clevs=np.arange(-2,2.25,0.25).tolist()#.remove(0))
    clevs.remove(0)
    levs_ticks=clevs
    assert len(clevs) % 2 == 0, 'N levels must be even.'
    cmap = mcolors.LinearSegmentedColormap.from_list(name='red_white_blue', 
                                                 colors =[(0, 0, 1), 
                                                          (1, 1., 1), 
                                                          (1, 0, 0)],
                                                 N=len(clevs)-1,
                                                 )
    return clevs,clevs_ticks,cmap
'''
clevs=np.arange(-20,22,2).tolist()#.remove(0))
clevs.remove(0)
#levs_ticks=clevs
assert len(clevs) % 2 == 0, 'N levels must be even.'
cmap = mcolors.LinearSegmentedColormap.from_list(name='red_white_blue', 
                                                 colors =[(0, 0, 1), 
                                                          (1, 1., 1), 
                                                          (1, 0, 0)],
                                                 N=len(clevs)-1,
                                                 )
clevsReg=0.1*np.asarray(clevs)
levs_ticks=clevsReg
'''
#####
###ZERO WHITE CUSTOM
######
'''
When making colormaps with LinearSegmentedColormap.from_list, you can supply a list of tuples of the form (value, color) (as opposed to simply a list of colors) where the values correspond to the relative positions of colors. The values must range from 0 to 1 so you will have to supply an intermediate color. In your case I might try this,

cmap = clr.LinearSegmentedColormap.from_list('custom blue',
                                             [(0,    '#ffff00'),
                                              (0.25, '#002266'),
                                              (1,    '#002266')], N=256)
and tweak color/value until satisfied. Credit goes to https://stackoverflow.com/a/25000108/5285918

##https://www.webucator.com/article/python-color-constants-module/
'''
'''
clevs=[-2.  , -1.75, -1.5 , -1.25 ,-1.,   -0.75 ,-0.5 , -0.25,      0.25 , 0.5,   0.75,
  1.   , 1.25 , 1.5 ,  1.75 , 2.   ]

levs_ticks=[-2.  , -1.5 , -1.   , -0.5  ,-0.25 ,   0.25 , 0.5  ,
  1.  , 1.5   ,  2.  ]
#levs = np.arange(-6,7,0.5)
assert len(clevs) % 2 == 0, 'N levels must be even.'
clevs=np.arange(-7,8,1).tolist()#.remove(0))
clevs.remove(0)
levs_ticks=clevs
assert len(clevs) % 2 == 0, 'N levels must be even.'
cmap= mcolors.LinearSegmentedColormap.from_list(name='blue_lightblue_yellow_red',
                                                         colors =[(0, '#4169E1'),
                                                                 (0.4,'#B0E0E6'),
                                                                 (0.5,'#FFFFFF'),
                                                                 (0.6,'#FFFF00'),
                                                             (1,'#FF0000')],
                                                 N=len(clevs)-1,
                                                 )
'''
######
## Hatched options from 'root'
######
#https://wikidocs.net/142073
#plt.rcParams['hatch.linewidth'] =1.5
###
#fig.add_axes([xmin, ymin, dx, dy])
############------------------------------------------
##
##Custom cmap
##
##https://towardsdatascience.com/beautiful-custom-colormaps-with-matplotlib-5bab3d1f0e72
##https://www.webucator.com/article/python-color-constants-module/
##
##hex_list = ['#0091ad', '#d6f6eb', '#fdf1d2', '#faaaae', '#ff57bb']
##cmap=get_continuous_cmap(hex_list)
############-----------------------------------------
#cmap = mpl.cm.get_cmap("GnBu").copy()
#cmap.set_under('w')
'''
Example for STD:
    clevs=np.arange(2,7.25,0.25)
    cmap = plt.get_cmap('CMRmap_r')
    new_cmap = truncate_colormap(cmap, 0.0, 0.7)
    lons, lats = np.meshgrid(eofCovList[0].coords['lon'],eofCovList[0].coords['lat'])
    new_cmap.set_under('w')
    CS1=axs[i].contourf(lons,lats,plotList[i][:,:],clevs,
            transform=ccrs.PlateCarree(),
            cmap=new_cmap,extend='both')
'''
'''
https://towardsdatascience.com/creating-colormaps-in-matplotlib-4d4de78a04b8
import numpy as np
import matplotlib.pyplot as plt
# create mock data
data = np.random.random([100, 100]) * 10

# define top and bottom colormaps
top = cm.get_cmap('Greens_r', 128) # r means reversed version
bottom = cm.get_cmap('Blues', 128)
# combine it all
newcolors = np.vstack((top(np.linspace(0, 1, 128)),
                       bottom(np.linspace(0, 1, 128))))
# create a new colormaps with a name of OrangeBlue
orange_blue = ListedColormap(newcolors, name='OrangeBlue')

'''
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    '''
    https://stackoverflow.com/questions/18926031/how-to-extract-a-subset-of-a-colormap-as-a-new-colormap-in-matplotlib
    arr = np.linspace(0, 50, 100).reshape((10, 10))
    fig, ax = plt.subplots(ncols=2)
    cmap = plt.get_cmap('jet')
    new_cmap = truncate_colormap(cmap, 0.2, 0.8)
    ax[0].imshow(arr, interpolation='nearest', cmap=cmap)
    ax[1].imshow(arr, interpolation='nearest', cmap=new_cmap)
    plt.show()
    '''
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def hex_to_rgb(value):
    '''
    Converts hex to rgb colours
    value: string of 6 characters representing a hex colour.
    Returns: list length 3 of RGB values'''
    value = value.strip("#") # removes hash symbol if present
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def rgb_to_dec(value):
    '''
    Converts rgb to decimal colours (i.e. divides each value by 256)
    value: list (length 3) of RGB values
    Returns: list (length 3) of decimal values'''
    return [v/256 for v in value]

def get_continuous_cmap(hex_list, float_list=None):
    ''' creates and returns a color map that can be used in heat map figures.
        If float_list is not provided, colour map graduates linearly between each color in hex_list.
        If float_list is provided, each color in hex_list is mapped to the respective location in float_list.

        Parameters
        ----------
        hex_list: list of hex code strings
        float_list: list of floats between 0 and 1, same length as hex_list. Must start with 0 and end with 1.

        Returns
        ----------
        colour map'''
    rgb_list = [rgb_to_dec(hex_to_rgb(i)) for i in hex_list]
    if float_list:
        pass
    else:
        float_list = list(np.linspace(0,1,len(rgb_list)))

    cdict = dict()
    for num, col in enumerate(['red', 'green', 'blue']):
        col_list = [[float_list[i], rgb_list[i][num], rgb_list[i][num]] for i in range(len(float_list))]
        cdict[col] = col_list
    cmp = mcolors.LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=256)
    return cmp

############--------------------------------
##
## Polar stereo
##
############--------------------------------
#label 
#https://stackoverflow.com/questions/61778059/how-to-add-more-lines-of-longitude-latitude-in-northpolarstereo-projection

#https://stackoverflow.com/questions/73959687/placement-of-latitude-labels-in-cartopy-with-polar-stereographic-projection

'''
fig,axs= plt.subplots(subplot_kw={'projection': ccrs.NorthPolarStereo()},figsize=(14,16))

data,lons=add_cyclic_point(par[0],coord=par[0]['lon'])
dataSIG,lonsSIG=add_cyclic_point(abs(par[2][:,:]),coord=par[0]['lon'])

#data,lons=add_cyclic_point(data0,coord=par[0]['lon'])
#dataSIG,lonsSIG=add_cyclic_point(data0SIG,coord=par[0]['lon'])

lats=par[0].lat
#lons, lats = np.meshgrid(par[0].lon,par[0].lat)

CS1=axs.contourf(lons,lats,data,clevs,
            transform=ccrs.PlateCarree(),
            cmap=cmap,extent='both')
levels=[0,rs,1.0]
cs = axs.contourf(lonsSIG,lats,dataSIG, transform=ccrs.PlateCarree(),levels=levels,
            hatches=["", ".."], alpha=0.)

'''
'''
###
##Clim
##
# Value for r_extent is obtained by trial and error
# get it with `ax.get_ylim()` after running this code
r_extent = 4651194.319
r_extent *= 1.7       #increase a bit for better result

latS,latN,lonW,lonE,latlim,lonlim=dom.coord_domain(domain)

lonlat_proj = ccrs.PlateCarree()
use_proj = ccrs.NorthPolarStereo(central_longitude=0)

fig,axs= plt.subplots(subplot_kw={'projection': ccrs.NorthPolarStereo(central_longitude=0)},figsize=(14,16))

lat,lon=climb.latlon(climW)
dataClimW0,lons0=add_cyclic_point(climW,coord=climW['lon'])
dataClimW,lons=add_cyclic_point(dataClimW0,lons0)

lats=climW.lat
#dataSIG,lonsSIG=add_cyclic_point(abs(par[2][:,:]),coord=par[0]['lon'])
cmapClim.set_under('w')
#lons, lats = np.meshgrid(clim[lon] ,clim[lat])
CS1=axs.contourf(lons,lats, dataClimW,clevsClim,
            transform=ccrs.PlateCarree(),
            cmap=cmapClim,extend='both')
#ax.set_extent([-180, 180, 50, 90], lonlat_proj)
# Draw the coastines for each subplot
axs.coastlines(color='grey')
#axs.add_feature(cfeature.BORDERS, linestyle=':', alpha=1)

geom,lonX,latY=get_geom('SibH')
axs.add_geometries([geom], facecolor='None',edgecolor='black',crs=ccrs.PlateCarree(), lw=2., alpha=1.0)
#axs.text(40, 60, 'SibH',color='orange',horizontalalignment='right',transform=ccrs.PlateCarree())

geom,lonX,latY=get_geom('SCAND20120')
axs.add_geometries([geom], facecolor='None',edgecolor='black',crs=ccrs.PlateCarree(), lw=2., alpha=1.0)
#axs.text(120, 60, 'SCAND20120',color='fuchsia',horizontalalignment='right',transform=ccrs.PlateCarree())

#raw graticule (meridian and parallel lines)
gls = axs.gridlines(draw_labels=True, crs=lonlat_proj, lw=1, color="grey",linestyle='--',
        y_inline=True, xlocs=range(-180,180,30), ylocs=range(20,90,15))

# Adjust the location of the subplots on the page to make room for the colorbar
fig.subplots_adjust(bottom=0.35, top=0.8, left=0.20, right=0.80,
                wspace=0.05, hspace=0.5)
# Add a colorbar axis at the bottom of the graph
#([xmin,ymin,dx,dy])
cbar_ax = fig.add_axes([0.2, 0.30, 0.6, 0.02])

cbar=fig.colorbar(CS1, cax=cbar_ax,orientation='horizontal',label='$\mathrm{m \cdot s^{-1} \cdot K^{-1}}$')


# Prep circular boundary
circle_path = mpath.Path.unit_circle()
circle_path = mpath.Path(circle_path.vertices.copy() * r_extent,
                           circle_path.codes.copy())

#set circular boundary
#this method will not interfere with the gridlines' labels 
axs.set_boundary(circle_path)
axs.set_frame_on(False)  #hide the boundary frame

plt.draw()  # Enable the use of `gl._labels`

# Reposition the tick labels
for ea in gls.label_artists:
    pos = ea.get_position()
    if pos[0] == 150:
      # print("Position:", pos, ea.get_text() )
        ea.set_position([180, pos[1]])
        

## Add a big title at the top
ofileC='FIG5'+str(l)+'_clim_'+variable+'_'+data+'_'+domain+'_'+season+'_'+str(iyr)+'_'+str(fyr)
plt.suptitle('('+str(l)+')',y=0.85)
fig.savefig(plotsDir+ofileC+'.png',format='png',bbox_inches='tight')
#print(plotsDir+ofileC)
plt.show()
# In[21]:
'''

############-------------------------------
##
##Variables custom plot
##
###########--------------------------------

def get_var_infos(var):
    """
        Get default informations about variable as: label, units, levels, cmap,
        extend for regular plot, differences and bias.
        Parameters
        ----------
        var : str
            Variable name. Options are:
            - 'snc', 'frac_snow' (Snow Cover Extent)
            - 'tas', 't2m', 'tmp (Near-Surface Air Temperature)
            - 'pr' (Total Precipitation)
            - 'ta' (Air Temperature)
        Returns
        -------
        label : str
            Name of the variable.
        units : str
            Usual units of the variable.
        levels, levels_diff, levels_bias : ndarray
            Levels.
        cmap, cmap_diff, cmap_bias : str
            Usual colormap for the varibale.
        extend, extend, extend : {{'neither', 'min', 'max', 'both'}}
            Where to assign unique colors to out-of-bounds data and draw
            "extensions" (triangles, by default) on the colorbar.
        Example
        -------
        #>>> import sys
        #>>> sys.path.insert(1, '/home/mlalande/notebooks/utils')
        #>>> import utils as u
        >>>
       # >>> label, units, \
            levels, cmap, extend, \
            levels_diff, cmap_diff, extend_diff, \
            levels_bias, cmap_bias, extend_bias = u.get_var_infos('snc')
    """

    # Snow Cover Extent
    if var in ['snc', 'frac_snow','snow_cover_extent','SNC','SCE']:
        label = 'Snow Cover Extent'
        units = '%'

        levels = np.arange(5, 115, 10)
        cmap = 'YlGnBu'
        extend = 'neither'

        levels_diff = np.arange(-25, 30, 5)
        cmap_diff = 'bwr'
        extend_diff ='both'

        levels_std = np.arange(0, 40, 5)
        cmap_std = 'rainbow'
        extend_std = 'both'
## Snow depth
    elif var in ['snd','snod','SND','snow_depth']:
        label = 'Snow Depth' 
        units = 'cm'

        levels = np.arange(5, 115, 10)
        cmap = 'YlGnBu'
        extend = 'neither'

        levels_diff = np.arange(-25, 30, 5)
        cmap_diff = 'bwr'
        extend_diff ='both'

        levels_std = np.arange(0, 40, 5)
        cmap_std = 'rainbow'
        extend_std = 'both'
   # Sea ice concentration
    elif var in ['sic','SIC']:
        label = 'Sea Ice Concentration'
        units = '%'

        levels = np.arange(0, 110, 10)
        cmap = 'YlGnBu'
        extend = 'neither'

        levels_diff = np.arange(-25, 30, 5)
        cmap_diff = 'bwr'
        extend_diff ='both'

        levels_std = np.arange(0, 40, 5)
        cmap_std = 'rainbow'
        extend_std = 'both'

    # Near-Surface Air Temperature
    elif var in ['tas', 't2m', 'tmp','T2M','air']:
        label = 'Near-Surface Air Temperature'
        units = '°C'

        levels = np.arange(-30, 35, 5)
        cmap = 'RdYlBu_r'
        extend = 'neither'

        levels_diff = np.arange(-1.0, 1.1, 0.1)
        cmap_diff = 'RdBu_r'
        extend_diff = 'both'

        levels_std = np.arange(0, 5.25, 0.25)
        cmap_std = 'rainbow'
        extend_std = 'both'
   #Sea surface temperaure

    elif var in ['sst','SSTK','SST']:
        label = 'Sea surface temperature'
        units = '°C'

        levels = np.arange(-10, 45, 5)
        cmap = 'viridis'
        extend = 'neither'

        levels_diff = np.arange(-1.0, 1.1, 0.1)
        cmap_diff = 'bwr'
        extend_diff = 'both'

        levels_std = np.arange(0, 1.6, 0.1)
        cmap_std = 'rainbow'
        extend_std = 'both'

    elif var in ['thf','lhtfl','shtfl']:
        label = 'Total,latent or sensible heat flux'
        units = 'W / m^2'
        levels = np.arange(-30, 200, 10)
        cmap = 'bwr'
        extend = 'neither'

        levels_diff = np.arange(-5.0, 5.25, 0.25)
        cmap_diff = 'bwr'
        extend_diff = 'both'

        levels_std = np.arange(0, 65, 5)
        cmap_std = 'rainbow'
        extend_std = 'both'

         # Sea Level Pressure
    elif var in ['slp','msl', 'prmsl', 'psl','MSL']:
        label = 'Mean sea level pressure'
        units = 'hPa'

        levels = np.arange(1000, 1040,4 )
        cmap = 'RdBu_r'
        extend = 'neither'

        levels_diff = np.arange(-5, 6, 1)
        cmap_diff = 'seismic'
        extend_diff = 'both'

        levels_std = np.arange(0, 5.5, 0.5)
        cmap_std = 'rainbow'
        extend_std = 'both'

#HGT
    elif var in ['hgt','z','g300','z200']:
        label = 'Geopotential height'
        units = 'm'

        levels = np.arange(1000,60000,500 )
        cmap = 'RdBu_r'
        extend = 'neither'

        levels_diff = np.arange(-60, 65, 5)
        cmap_diff = 'seismic'
        extend_diff = 'both'

        levels_std = np.arange(0, 35, 5)
        cmap_std = 'rainbow'
        extend_std = 'both'

    # Total Precipitation
    elif var in ['pr', 'tp', 'precip']:
        label = 'Total Precipitation'
        units = 'mm/day'

        levels = np.arange(0, 5, 0.5)
        cmap = 'BrBG'
        extend = 'neither'

        levels_diff = np.arange(-1, 1, 0.2)
        cmap_diff = 'rainbow'
        extend_diff = 'both'

        levels_std = np.arange(-5, 5, 1)
        cmap_std = 'rainbow'
        extend_std = 'both'

    # Snowfall
    elif var == 'prsn':
        label = 'Snowfall'
        units = 'mm/day'
        cmap = 'DryWet'
        levels = np.arange(0, 5, 0.5)

    # Eastward Wind
    elif var in ['ua','u','uwnd','u925','v','vwnd','v925','U','V']:
        label = 'Wind'
        units = 'm/s'
        levels = np.arange(-7, 8, 1)
        cmap = 'bwr'
        extend = 'neither'

        levels_diff = np.arange(-0.3, 0.35, 0.05)
        cmap_diff = 'bwr'
        extend_diff = 'both'

        levels_std = np.arange(0, 3.25, 0.25)
        cmap_std = 'rainbow'
        extend_std = 'both'


    else:
        raise ValueError(
            f"""Invalid variable argument: '{var}'. 
             """
        )

    print('0-label','1- units', \
        '2-levels', '3-cmap', '4-extend', \
        '5-levels_diff', '6-cmap_diff', '7-extend_diff', \
        '8-levels_std', '9-cmap_std', '10-extend_std')
    info= {'label':label, 'units':units, \
            'levels':levels, 'cmap':cmap, 'extend':extend, \
            'levels_diff':levels_diff, 'cmap_diff':cmap_diff, 'extend_diff':extend_diff, \
            'levels_std':levels_std, 'cmap_std':cmap_std, 'extend_std':extend_std}
    return info


###
##NUMBER in CONTOUR
#https://scitools.org.uk/cartopy/docs/latest/gallery/scalar_data/contour_labels.html
#
import cartopy.crs as ccrs
import matplotlib.pyplot as plt


def sample_data(shape=(73, 145)):
    """Return ``lons``, ``lats`` and ``data`` of some fake data."""
    import numpy as np

    nlats, nlons = shape
    lats = np.linspace(-np.pi / 2, np.pi / 2, nlats)
    lons = np.linspace(0, 2 * np.pi, nlons)
    lons, lats = np.meshgrid(lons, lats)
    wave = 0.75 * (np.sin(2 * lats) ** 8) * np.cos(4 * lons)
    mean = 0.5 * np.cos(2 * lats) * ((np.sin(2 * lats)) ** 2 + 2)

    lats = np.rad2deg(lats)
    lons = np.rad2deg(lons)
    data = wave + mean

    return lons, lats, data


def main():
    fig = plt.figure()

    # Setup a global EckertIII map with faint coastlines.
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.EckertIII())
    ax.set_global()
    ax.coastlines('110m', alpha=0.1)

    # Use the same sample data as the waves example, but make it
    # more dependent on y for more interesting contours.
    x, y, z = sample_data((20, 40))
    z = z * -1.5 * y

    # Add colourful filled contours.
    filled_c = ax.contourf(x, y, z, transform=ccrs.PlateCarree())

    # And black line contours.
    line_c = ax.contour(x, y, z, levels=filled_c.levels,
                        colors=['black'],
                        transform=ccrs.PlateCarree())

    # Uncomment to make the line contours invisible.
    # plt.setp(line_c.collections, visible=False)

    # Add a colorbar for the filled contour.
    fig.colorbar(filled_c, orientation='horizontal')

    # Use the line contours to place contour labels.
    ax.clabel(
        line_c,  # Typically best results when labelling line contours.
        colors=['black'],
        manual=False,  # Automatic placement vs manual placement.
        inline=True,  # Cut the line where the label will be placed.
        fmt=' {:.0f} '.format,  # Labes as integers, with some extra space.
    )

    plt.show()


if __name__ == '__main__':
    main()


##to develop
##https://stackoverflow.com/questions/55303911/add-polygon-box-to-cartopy-python
#geom = geometry.box(minx=-79,maxx=-70,miny=-44.7,maxy=-49.3)
#axs.add_geometries([geom], facecolor='None',edgecolor='green',crs=cartopy.crs.PlateCarree(), alpha=0.3)

def get_geom(region):
    '''
    In figure:
    geom,lonX,latY=get_geom('HK')
    axs.add_geometries([geom], facecolor='None',edgecolor='black',crs=ccrs.PlateCarree(), alpha=0.3)
    axs.text(lonX+3, latY-2, 'HK',horizontalalignment='right',transform=ccrs.PlateCarree())

    '''
    latS,latN,lonW,lonE,latlim,lonlim=dom.coord_domain(region)
    geom = geometry.box(minx=lonW,maxx=lonE,miny=latS,maxy=latN)
    return geom,lonW,latN

###########
## WIND and ARROWs
########
#import matplotlib.pyplot as plt
#import cartopy.crs as ccrs
#from cartopy.examples.arrows import sample_data
#https://scitools.org.uk/cartopy/docs/v0.18/gallery/streamplot.html
#
#def main():
#    fig = plt.figure(figsize=(10, 5))
#    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
#    ax.set_extent([-90, 75, 10, 85], crs=ccrs.PlateCarree())
#    ax.coastlines()
#
#   x, y, u, v, vector_crs = sample_data(shape=(80, 100))
#    magnitude = (u ** 2 + v ** 2) ** 0.5
#    ax.streamplot(x, y, u, v, transform=vector_crs,
#                  linewidth=2, density=2, color=magnitude)
#    plt.show()
#
#
#if __name__ == '__main__':
#    main()

########
##
#Multiple single png into a pdf:
##
########
#import matplotlib.pyplot as plt
#https://stackoverflow.com/questions/19189488/use-a-loop-to-plot-n-charts-python
#nameList=['a','b','c','d']
#x=[[1,2,3,4],[1,2,3,4],[1,2,3,4],[1,2,3,4]]
#y=[[1,2,3,4],[1,2,3,4],[1,2,3,4],[1,2,3,4]]
#for i in range(len(x)):
    #plt.figure()
    #plt.plot(x[i],y[i])
    # Show/save figure as desired.
    #plt.savefig(nameList[i]+'.png',format='png',bbox_inches='tight')
    #plt.show()
# Can show all four figures at once by calling plt.show() here, outside the loop.
#plt.show()

#imagelist=['a.png','b.png','c.png','d.png']
#from fpdf import FPDF
#https://pyfpdf.readthedocs.io/en/latest/reference/image/index.html
#pdf = FPDF()
# imagelist is the list with all image filenames
#for image in imagelist:
#    pdf.add_page()
#    pdf.image(image,x = None, y = None, w = 0, h = 0, type = '', link = '')
#pdf.output("yourfile.pdf", "F")


###########3
###Multiplot
#####
#https://kpegion.github.io/Pangeo-at-AOES/examples/multi-panel-cartopy.html


###
#CUSTOM PLOT
##https://fabienmaussion.info/climate_system/projects/04_Getting_started_Antarctica.html
