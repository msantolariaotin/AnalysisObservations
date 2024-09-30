'''
for season in DJF JFM;do
python climb_prmsl.noaaV2c.py ${season}
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

domain='NAEext'
season='DJF'#sys.argv[1]
exp='obs'

model='noaaV2c'
variable='prmsl'
info=obs.get_obs(variable,model)

print(info)
fileName=info.get('filename')
ds= xr.open_dataset(sourceData+model+'/'+fileName)[variable]
# replace all values equal to -1000 with np.nan
#ds = ds_fill[variable].where(ds_fill[variable] != -1000.)
if info.get('gridlon')=='0_360':
    ds=dom.shifting_grid(ds)
lat,lon=climb.latlon(ds)
ylat=ds.coords[lat]
xlon=ds.coords[lon]
if ds.units=='Pa':
    ds=ds/100
    unit='hPa'
else:
    unit=ds.units
field=dom.field_dom(ds,domain)

if season[0]=='m':
    rmon=int(season.split('mon')[1])
    vals,anoms=climb.monthly_selection(field,rmon,iyr,fyr)
else:
    vals,anoms=climb.seasonal_selection(field,season,iyr,fyr)

###
##Detrended anomalies
###
anoms_detrend=climb.detrend_dim(vals, 'time', deg=1)
###
##Clim and std dev
###
clim=vals.mean('time')
##Plotting
###
infoplot=myplot.get_var_infos(variable)
cmapClim=infoplot.get('cmap')
units=infoplot.get('units')
clevsClim=infoplot.get('levels')

cmapStd=infoplot.get('cmap_std')
clevsStd=infoplot.get('levels_std')

clevsT=np.arange(-1.0,1.1,0.1)#infoplot.get('levels_diff')
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
#plt.show()
