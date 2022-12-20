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


iyr=1979
fyr=2019

domain='MA'
season='DJF'
exp='obs'

model='era5'
variable='tas'
#info=obs.get_obs(variable,model)
#print(info)
fileName='t925_era5_NH_mon_1979-2020.nc'#info.get('filename')
ds= xr.open_dataset(sourceData+model+'/'+fileName)['t'][:,0,:,:]
# replace all values equal to -1000 with np.nan
#ds = ds_fill[variable].where(ds_fill[variable] != -1000.)
#if info.get('gridlon')=='0_360':
##    ds=dom.shifting_grid(ds)
lat,lon=climb.latlon(ds)
ylat=ds.coords[lat]
xlon=ds.coords[lon]
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
std=vals.std('time')
std_det=anoms_detrend.std('time')

###
##Spatial trend
###
par=climb.trend_vect(anoms.time,anoms,'time')
trend=par[0]
intercept=par[1]
rvalue=par[2]
pvalue=par[3]
stderr=par[4]

###
##Ratio= 10*trend/std_dev(anom_notrend or anom)
###
ratio=10*abs(trend)/std_det

###
##U

variableU='u'
#info=obs.get_obs(variable,model)
#print(info)
fileNameU='u925_era5_NH_mon_1979-2020.nc'#info.get('filename')
dsU= xr.open_dataset(sourceData+model+'/'+fileNameU)['u'][:,0,:,:]
# replace all values equal to -1000 with np.nan
#ds = ds_fill[variable].where(ds_fill[variable] != -1000.)
#if info.get('gridlon')=='0_360':
 #   ds=dom.shifting_grid(ds)
latU,lonU=climb.latlon(dsU)
ylatU=dsU.coords[latU]
xlonU=dsU.coords[lonU]
unitU=dsU.units

fieldU=dom.field_dom(dsU,domain)

if season[0]=='m':
    rmon=int(season.split('mon')[1])
    valsU,anomsU=climb.monthly_selection(fieldU,rmon,iyr,fyr)
else:
    valsU,anomsU=climb.seasonal_selection(fieldU,season,iyr,fyr)

###
##Detrended anomalies
###
anomsU_detrend=climb.detrend_dim(valsU, 'time', deg=1)
###
##Clim and std dev
###
climU=valsU.mean('time')
stdU=anomsU.std('time')
stdU_det=anomsU_detrend.std('time')


variableV='v'
#info=obs.get_obs(variable,model)
#print(info)
fileNameV='v925_era5_NH_mon_1979-2020.nc'#info.get('filename')
dsV= xr.open_dataset(sourceData+model+'/'+fileNameV)['v'][:,0,:,:]
# replace all values equal to -1000 with np.nan
#ds = ds_fill[variable].where(ds_fill[variable] != -1000.)
#if info.get('gridlon')=='0_360':
 #   ds=dom.shifting_grid(ds)
latV,lonV=climb.latlon(dsV)
ylatV=dsV.coords[latV]
xlonV=dsV.coords[lonV]
c

fieldV=dom.field_dom(dsV,domain)

if season[0]=='m':
    rmon=int(season.split('mon')[1])
    valsV,anomsV=climb.monthly_selection(fieldV,rmon,iyr,fyr)
else:
    valsV,anomsV=climb.seasonal_selection(fieldV,season,iyr,fyr)

###
##Detrended anomalies
###
anomsV_detrend=climb.detrend_dim(valsV, 'time', deg=1)
###
##Clim and std dev
###
climV=valsV.mean('time')
stdV=anomsV.std('time')
stdV_det=anomsV_detrend.std('time')




###

infoplot=myplot.get_var_infos(variable)
cmapClim=infoplot.get('cmap')
units=infoplot.get('units')
clevsClim=infoplot.get('levels')

cmapStd=infoplot.get('cmap_std')
clevsStd=infoplot.get('levels_std')

clevsT=infoplot.get('levels_diff')
cmapT=infoplot.get('cmap_diff')

clevsRatio=np.arange(0,1.1,0.1)
cmapRatio='OrRd'

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
            cmap=cmapClim,extent='both')
##https://stackoverflow.com/questions/72529063/how-to-adjust-the-density-of-global-divergent-wind-vector-using-quiver
Q = axs.quiver(lons,lats,climUd,climVd,scale=10,color='dimgrey')
# make quiver key.
qk = axs.quiverkey(Q, 0.15, -0.20, 5, '5 m/s', labelpos='W')

# Draw the coastines for each subplot
axs.coastlines()
axs.add_feature(cfeature.BORDERS, linestyle=':', alpha=1)
#axs.add_feature(cfeature.NaturalEarthFeature('physical', 'ocean', '50m', edgecolor='None', facecolor='None'))
if extentTF==True:
    axs.set_extent([lonW, lonE, latS,latN])
    # Longitude and latitude labels
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
fig.savefig(plotsDir+ofileC+'.png',format='png')
print('Figure save at ',plotsDir, 'as',ofileC)
plt.show()
# In[21]:
'''
###
##STD
###
fig,axs= plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()},figsize=(8,10))
lat,lon=climb.latlon(clim)
lons, lats = np.meshgrid(std[lon] ,std[lat])
CS1=axs.contourf(lons,lats, std,clevsStd,
            transform=ccrs.PlateCarree(),
            cmap=cmapStd,extent='both')
# Draw the coastines for each subplot
axs.coastlines()
axs.add_feature(cfeature.BORDERS, linestyle=':', alpha=1)
#axs.add_feature(cfeature.NaturalEarthFeature('physical', 'ocean', '50m', edgecolor='None', facecolor='None'))
if extentTF==True:
    axs.set_extent([lonW, lonE, latS,latN])
    # Longitude and latitude labels
# Adjust the location of the subplots on the page to make room for the colorbar
fig.subplots_adjust(bottom=0.35, top=0.7, left=0.20, right=0.80,
                wspace=0.05, hspace=0.5)
# Add a colorbar axis at the bottom of the graph
#([xmin,ymin,dx,dy])
cbar_ax = fig.add_axes([0.2, 0.35, 0.6, 0.02])
# Draw the colorbar
cbar=fig.colorbar(CS1, cax=cbar_ax,orientation='horizontal',label='%s'%(units))
## Add a big title at the top
ofileS='anoms_std_'+variable+'_'+model+'_'+exp+'_'+domain+'_'+season+'_'+str(iyr)+'_'+str(fyr)
plt.suptitle(ofileS,y=0.68)
fig.savefig(plotsDir+ofileS+'.png',format='png')
print('Figure save at ',plotsDir, 'as',ofileS)
#plt.show()


###
##STD det
###
fig,axs= plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()},figsize=(8,10))

lat,lon=climb.latlon(clim)
lons, lats = np.meshgrid(std[lon] ,std[lat])
CS1=axs.contourf(lons,lats, std_det,clevsStd,
            transform=ccrs.PlateCarree(),
            cmap=cmapStd,extent='both')
# Draw the coastines for each subplot
axs.coastlines()
axs.add_feature(cfeature.BORDERS, linestyle=':', alpha=1)
#axs.add_feature(cfeature.NaturalEarthFeature('physical', 'ocean', '50m', edgecolor='None', facecolor='None'))
if extentTF==True:
    axs.set_extent([lonW, lonE, latS,latN])
    # Longitude and latitude labels
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
fig.savefig(plotsDir+ofileSdet+'.png',format='png')
print('Figure save at ',plotsDir, 'as',ofileSdet)
#plt.show()


###
##Spatial trend
###
fig,axs= plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()},figsize=(8,10))

lat,lon=climb.latlon(par[0])
lons, lats = np.meshgrid(par[0][lon] ,par[0][lat])
CS1=axs.contourf(lons,lats, 10*par[0][:,:],clevsT,
                transform=ccrs.PlateCarree(),
                cmap=cmapT,extent='both')
levels=[0,0.1,1.0]
cs = axs.contourf(lons,lats,par[3][:,:], transform=ccrs.PlateCarree(),levels=levels,
                hatches=["+", ""], alpha=0.)

# Draw the coastines for each subplot
axs.coastlines()
axs.add_feature(cfeature.BORDERS, linestyle=':', alpha=1)
#axs.add_feature(cfeature.NaturalEarthFeature('physical', 'ocean', '50m', edgecolor='None', facecolor='None'))
if extentTF==True:
    axs.set_extent([lonW, lonE, latS,latN])
    # Longitude and latitude labels
# Adjust the location of the subplots on the page to make room for the colorbar
fig.subplots_adjust(bottom=0.35, top=0.7, left=0.20, right=0.80,
                wspace=0.05, hspace=0.5)
# Add a colorbar axis at the bottom of the graph
#([xmin,ymin,dx,dy])
cbar_ax = fig.add_axes([0.2, 0.35, 0.6, 0.02])
# Draw the colorbar
cbar=fig.colorbar(CS1, cax=cbar_ax,orientation='horizontal',label='%s/dec'%(units))
## Add a big title at the top
ofileT='spatialtrend_'+variable+'_'+model+'_'+exp+'_'+domain+'_'+season+'_'+str(iyr)+'_'+str(fyr)
plt.suptitle(ofileT,y=0.68)
fig.savefig(plotsDir+ofileT+'.png',format='png')
print('Figure save at ',plotsDir, 'as',ofileT)
#plt.show()
###
##Ratio
###
fig,axs= plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()},figsize=(8,10))

lat,lon=climb.latlon(clim)
lons, lats = np.meshgrid(std[lon] ,std[lat])
CS1=axs.contourf(lons,lats, ratio,clevsRatio,
            transform=ccrs.PlateCarree(),
            cmap=cmapRatio,extent='both')
# Draw the coastines for each subplot
axs.coastlines()
axs.add_feature(cfeature.BORDERS, linestyle=':', alpha=1)
#axs.add_feature(cfeature.NaturalEarthFeature('physical', 'ocean', '50m', edgecolor='None', facecolor='None'))
if extentTF==True:
    axs.set_extent([lonW, lonE, latS,latN])
    # Longitude and latitude labels
# Adjust the location of the subplots on the page to make room for the colorbar
fig.subplots_adjust(bottom=0.35, top=0.7, left=0.20, right=0.80,
                wspace=0.05, hspace=0.5)
# Add a colorbar axis at the bottom of the graph
#([xmin,ymin,dx,dy])
cbar_ax = fig.add_axes([0.2, 0.35, 0.6, 0.02])
# Draw the colorbar
cbar=fig.colorbar(CS1, cax=cbar_ax,orientation='horizontal',label='%s'%(units))
## Add a big title at the top
ofileSdet='ratio_trend_anoms_det_'+variable+'_'+model+'_'+exp+'_'+domain+'_'+season+'_'+str(iyr)+'_'+str(fyr)
plt.suptitle(ofileSdet,y=0.68)
fig.savefig(plotsDir+ofileSdet+'.png',format='png')
print('Figure save at ',plotsDir, 'as',ofileSdet)
#plt.show()


##Computing time series of seasonal spatial average ----------------------------
print('Computing time series of seasonal spatial average')
ts_season=climb.spatial_average(anoms)
plotnameTs='timeseries_anoms_'+variable+'_'+model+'_'+exp+'_'+domain+'_'+season+'_'+str(iyr)+'_'+str(fyr)
np.savetxt(resultsDir+plotnameTs+'.txt',ts_season)
#print('saving .txt at',resultsDir+plotnameTs)

##Trends---------------------------------------
xd=np.arange(0,len(ts_season),1)
par=stats.linregress(xd,ts_season)
'''
'''
print('stats',par)
parmed=stats.theilslopes(ts_season,xd,0.90)
print('theilsen',parmed)
parmk=mk.original_test(ts_season)
print('trend', 'h', 'p', 'z', 'Tau', 's', 'var_s', 'slope', 'intercept')
print(parmk)
trendsOf=[par,parmed,parmk]
'''
'''
fig, ax = plt.subplots(figsize=(12, 8))
my_ticks=xd+iyr
print(xd,my_ticks)
ticks=np.arange(iyr,fyr+1,1)
#ax.plot(xd,ts_season,color='b', linewidth=2)
ax.plot(xd, ts_season, 'b:')
#ax.plot(xd, parmed[1] + parmed[0] * xd, 'r-',label="%s (%.2f)" % ('Theil-Sen', 10*parmed[0]))
#ax.plot(xd, parmed[1] + parmed[2] * xd, 'r--')
#ax.plot(xd, parmed[1] + parmed[3] * xd, 'r--')
ax.plot(xd, par[1] + par[0] * xd, 'b-',label="%s (%.2f);%s (%.2f)" % ('Lin reg', 10*par[0],'p-val',par[3]))
#ax.plot(xd, parmk[8] + parmk[7] * xd,'p-',label="%s (%.2f) sig %s" % ('MK',i 10*parmk[7],parmk[0]))

#ax.axhline(0, color='k')
#ax.set_ylim(0, 3)
ax.set_xlabel('Year')
ax.set_ylabel('%s'%(units))
ax.set_title(plotnameTs)
freq=5
plt.xticks(xd,my_ticks, rotation='vertical')
plt.legend()
plt.savefig(plotsDir+plotnameTs+'.png',format='png')
#plt.show()
print('save at', resultsDir,plotnameTs)
'''
'''
textfile = open(resultsDirTs+'trends_'+plotnameTs, "w")
for element in trendsOf:
    textfile.write(str(element) + "\n")
textfile.close()

'''
