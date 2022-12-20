'''
for season in DJF JFM;do
python climb_hgt.noaaV2c.py ${season}
done
'''
import sys
sys.path.insert(1, '/home/msantolaria/Documents/MyPythonLibrary/ClimAnag/')
import climbasis as climb
from climbasis import *
import domain as dom
import myplot
import glob
from eofs.xarray import Eof
from myplot import *
import obsinfo as obs
from obsinfo import *



source='/home/msantolaria/Documents/MyResearch/AnalysisObservations/'
resultsDir=source + 'Results/'
plotsDir=source + 'Plots/'

sourceData='/home/msantolaria/Documents/Data/'


iyr=1960
fyr=2014

domain='NHeq'
season='DJF'#sys.argv[1]
exp='obs'

model='noaaV2c'
variable='z_star200'
level=200
info=obs.get_obs('hgt',model)

fileName=info.get('filename')
dsall= xr.open_dataset(sourceData+model+'/'+fileName)['hgt']
ds=dsall.sel(level=level)
ds=dom.shifting_grid(ds)
lat,lon=climb.latlon(ds)
ylat=ds.coords[lat]
xlon=ds.coords[lon]
unit=ds.units
field=dom.field_dom(ds,domain)


##Computing z*=z-[z]
zonal=field.mean(dim=lon)
z_star=field-zonal


if season[0]=='m':
    rmon=int(season.split('mon')[1])
    vals,anoms=climb.monthly_selection(z_star,rmon,iyr,fyr)
else:
    vals,anoms=climb.seasonal_selection(z_star,season,iyr,fyr)


###
##Detrended anomalies
###
anoms_detrend=climb.detrend_dim(vals, 'time', deg=2)
###
##Clim and std dev
###
clim=vals.mean('time')
std=anoms.std('time')
std_det=anoms_detrend.std('time')
###
##Plotting
###
infoplot=myplot.get_var_infos('hgt')
cmapClim=infoplot.get('cmap')
units=infoplot.get('units')
clevsClim=np.arange(-200,220,20)#infoplot.get('levels')

cmapStd=infoplot.get('cmap_std')
clevsStd=np.arange(0,70,10)#infoplot.get('levels_std')

clevsT=np.arange(-14.0,16,2)#infoplot.get('levels_diff')
cmapT=infoplot.get('cmap_diff')

clevsRatio=np.arange(0,1.1,0.1)
cmapRatio='YlGn'

extend=infoplot.get('extend')
extentTF=False

###
##Clim
##
latS,latN,lonW,lonE,latlim,lonlim=dom.coord_domain(domain)

fig,axs= plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()},figsize=(8,10))
lat,lon=climb.latlon(clim)
lons, lats = np.meshgrid(clim[lon] ,clim[lat])
CS1=axs.contourf(lons,lats, clim,clevsClim,
            transform=ccrs.PlateCarree(),
            cmap=cmapClim,extend='both')
# Draw the coastines for each subplot
axs.coastlines()
axs.add_feature(cfeature.BORDERS, linestyle=':', alpha=1)
gl = axs.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
              linewidth=2, color='gray', alpha=0.0, linestyle='-')
gl.xlabels_top = False
gl.ylabels_right = False
gl.set_ylocator=mticker.FixedLocator(np.arange(latS,latN,10),0)
gl.set_xlocator = mticker.FixedLocator(np.arange(lonW,lonE,20),0)
gl.set_xformatter = LongitudeFormatter
gl.set_yformatter = LatitudeFormatter
gl.set_xlabel_style = {'color': 'black'}
gl.set_xlabel_style = {'color': 'black'}
# Adjust the location of the subplots on the page to make room for the colorbar
fig.subplots_adjust(bottom=0.35, top=0.7, left=0.20, right=0.80,
                wspace=0.05, hspace=0.5)
# Add a colorbar axis at the bottom of the graph
#([xmin,ymin,dx,dy])
cbar_ax = fig.add_axes([0.2, 0.35, 0.6, 0.02])
# Draw the colorbar
cbar=fig.colorbar(CS1, cax=cbar_ax,orientation='horizontal',label='%s'%(units))
## Add a big title at the top
ofileC='clim_'+variable+'_'+model+'_'+exp+'_'+domain+'_'+season+'_'+str(iyr)+'_'+str(fyr)
plt.suptitle(ofileC,y=0.68)
fig.savefig(plotsDir+ofileC+'.png',format='png',bbox_inches='tight')
print('Figure save at ',plotsDir, 'as',ofileC)
plt.show()
# In[21]:

###
##STD det
###
fig,axs= plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()},figsize=(8,10))

lat,lon=climb.latlon(clim)
lons, lats = np.meshgrid(std[lon] ,std[lat])
CS1=axs.contourf(lons,lats, std_det,clevsStd,
            transform=ccrs.PlateCarree(),
            cmap=cmapStd,extend='both')
# Draw the coastines for each subplot
axs.coastlines()
axs.add_feature(cfeature.BORDERS, linestyle=':', alpha=1)
gl = axs.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
              linewidth=2, color='gray', alpha=0.0, linestyle='-')
gl.xlabels_top = False
gl.ylabels_right = False
gl.set_ylocator=mticker.FixedLocator(np.arange(latS,latN,10),0)
gl.set_xlocator = mticker.FixedLocator(np.arange(lonW,lonE,20),0)
gl.set_xformatter = LongitudeFormatter
gl.set_yformatter = LatitudeFormatter
gl.set_xlabel_style = {'color': 'black'}
gl.set_xlabel_style = {'color': 'black'}
# Adjust the location of the subplots on the page to make room for the colorbar
fig.subplots_adjust(bottom=0.35, top=0.7, left=0.20, right=0.80,
                wspace=0.05, hspace=0.5)
# Add a colorbar axis at the bottom of the graph
#([xmin,ymin,dx,dy])
cbar_ax = fig.add_axes([0.2, 0.35, 0.6, 0.02])
# Draw the colorbar
cbar=fig.colorbar(CS1, cax=cbar_ax,orientation='horizontal',label='%s'%(units))
## Add a big title at the top
ofileSdet='anoms_det_std_'+variable+'_'+model+'_'+exp+'_'+domain+'_'+season+'_'+str(iyr)+'_'+str(fyr)
plt.suptitle(ofileSdet,y=0.68)
fig.savefig(plotsDir+ofileSdet+'.png',format='png',bbox_inches='tight')
print('Figure save at ',plotsDir, 'as',ofileSdet)
plt.show()

