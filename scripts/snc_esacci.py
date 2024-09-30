'''
for i in HK HM TP;do
python snc_esacci.py ${i}
done
'''
import sys
sys.path.insert(1, '/home/msantolaria/Documents/MyPythonLibrary/ClimAnag/')
import climbasis as climb
from climbasis import *
import domain as dom
import myplot
import glob


# In[3]:


source='/home/msantolaria/Documents/MyResearch/HMA-ClimAnalogs/Observations/'
resultsDir=source + 'Results/'
plotsDir=source + 'Plots/'
resultsDirT=source+ 'SpatialTrend/'
resultsDirTs=source+'TimeSeries/'
#resultsDir='/media/maria/MARIAEXT2/WinterTrendsObs/'
sourceData='/home/msantolaria/Documents/Data/'


# In[5]:


iyr=1982
imon=1
fmon=12
fyr=2014
domain=sys.argv[1]
season='DJFMA'
variable='snc'
units='%'
exp='obs'
decomp='original'
subdomain=domain
model='esacci'
fileName='ESACCI-L3C_SNOW-SCFG-AVHRR_MERGED-fv1.0_HKH_gapfilled_icefilled_montlhy_1982-2014.nc'
sourceData=sourceData+model+'/'
variable='snc'
ds = xr.open_dataset(sourceData+fileName)[variable]
units='%'
#print(ds0Y.units,'to',unitsY)

field=dom.field_dom(ds.sel(time=slice(str(iyr)+"-"+str(imon), str(fyr)+"-"+str(fmon))),domain)
ylat=field.coords['lat']
xlon=field.coords['lon']

print('field shape',field.shape)

if season[0]=='m':
    rmon=int(season.split('mon')[1])
    vals,anoms=climb.monthly_selection(field,rmon,iyr,fyr)
else:
    if season[0]=='D':
        print('Winter',iyr,'-',iyr+1,fyr-1,'-',fyr)
        vals,anoms=climb.seasonal_selection(field,season,6,iyr,6,fyr)
    else:
        vals,anoms=climb.seasonal_selection(field,season,1,iyr,12,fyr)

'''
anoms_detrend=climb.detrend_dim(vals, 'time', deg=1)

clim=vals.mean('time')
std=vals.std('time')
std_det=anoms_detrend.std('time')
#--------------------------------
subdomain=domain


##Clim-----------------
#clevs=np.arange(0,14.25,0.25)
clevs=np.arange(0,100,10)
color='GnBu'
exp='obs'
decomp='original'
lat,lon=climb.latlon(clim)
figclim=myplot.oneplot_ds_paper(variable=variable, decomp=decomp, exp=exp, model=model, clevs=clevs, color=color, units=units, subdomain=subdomain, ds=clim, xlon=clim[lon], ylat=clim[lat], season=season, iyr=iyr, fyr=fyr,extent=False,geom=True)
ofileC='clim_'+variable+'_'+model+'_'+exp+'_'+decomp+'_'+subdomain+'_'+season+'_'+str(iyr)+'_'+str(fyr)
figclim.savefig(plotsDir+ofileC+'.png',format='png')
print('Figure save at ',plotsDir, 'as',ofileC)
plt.show()
# In[21]:

##STD-----------------
clevs=np.arange(0,26,2)
color='rainbow'
exp='obs'
decomp='original'
lat,lon=climb.latlon(std)
figstd=myplot.oneplot_ds_paper(variable=variable, decomp=decomp, exp=exp, model=model, clevs=clevs, color=color, units=units, subdomain=subdomain, ds=std, xlon=clim[lon], ylat=clim[lat], season=season, iyr=iyr, fyr=fyr,extent=False,geom=True)
ofileS='anoms_std_'+variable+'_'+model+'_'+exp+'_'+decomp+'_'+subdomain+'_'+season+'_'+str(iyr)+'_'+str(fyr)
figstd.savefig(plotsDir+ofileS+'.png',format='png')
print('Figure save at ',plotsDir, 'as',ofileS)
plt.show()
###MEAN DET
clevs=np.arange(0,26,2)
color='rainbow'
exp='obs'
decomp='original-detrend'
lat,lon=climb.latlon(std_det)
figstd_det=myplot.oneplot_ds_paper(variable=variable, decomp=decomp, exp=exp, model=model, clevs=clevs, color=color, units=units, subdomain=subdomain, ds=std_det, xlon=clim[lon], ylat=clim[lat], season=season, iyr=iyr, fyr=fyr,extent=False,geom=True)
ofileS='anoms_det_mean_'+variable+'_'+model+'_'+exp+'_'+decomp+'_'+subdomain+'_'+season+'_'+str(iyr)+'_'+str(fyr)
figstd_det.savefig(plotsDir+ofileS+'.png',format='png')
print('Figure save at ',plotsDir, 'as',ofileS)
#-----------------------------------
print('Computing spatial trend')

coords2D={anoms.dims[1]: anoms.coords[anoms.dims[1]],anoms.dims[2]: anoms.coords[anoms.dims[2]]}
zero2D= xr.DataArray(np.zeros((anoms.shape[1],anoms.shape[2])),coords=coords2D,dims=[anoms.dims[1],anoms.dims[2]])
trend=xr.zeros_like(zero2D)
intercept=xr.zeros_like(zero2D)
rvalue=xr.zeros_like(zero2D)
pvalue=xr.zeros_like(zero2D)
stderr=xr.zeros_like(zero2D)

par=[]
for j in range(anoms.shape[1]):
    for i in range(anoms.shape[2]):
        if np.isnan(anoms[:,j,i].values).all()==True:
            trend[j,i]=np.nan
            intercept[j,i]=np.nan
            rvalue[j,i]=np.nan
            pvalue[j,i]=np.nan
            stderr[j,i]=np.nan
        else:
            Y=anoms[:,j,i]    
            #print(anoms[:,j,i])
            Y2=Y[np.logical_not(np.isnan(Y.values))]
            X=np.arange(0,len(vals[:,j,i]),1)
            X2=X[np.logical_not(np.isnan(Y.values))]
            par=stats.linregress(X2,Y2)
            trend[j,i]=par[0]
            intercept[j,i]=par[1]
            rvalue[j,i]=par[2]
            pvalue[j,i]=par[3]
            stderr[j,i]=par[4]
par=[trend,intercept,rvalue,pvalue,stderr]
print('0-trend','1-intercept','2-rvalue','3-pvalue','4-stderr')

ds_trend= xr.Dataset(
data_vars=dict(
           trend=([field.dims[1], field.dims[2]],par[0].values),
           intercept=([field.dims[1], field.dims[2]],par[1].values),
           rvalue=([field.dims[1], field.dims[2]],par[2].values),
           pvalue=([field.dims[1],field.dims[2]],par[3].values),
           stderr=([field.dims[1],field.dims[2]],par[4].values),
           ),
   coords=dict(
           latitude=([field.dims[1]],field.coords[field.dims[1]].values),
           longitude=([field.dims[2]],field.coords[field.dims[2]].values),),
attrs=dict(description="spatial trend"),
        )
ofileT='spatialtrend'+'_' +variable+ '.obs'+'_'+season+'_'+model+'_historical_'+domain+'_'+str(iyr)+'_'+str(fyr)
new_filenameT = resultsDirT+ofileT+'.nc'
print ('saving to ', new_filenameT)
ds_trend.to_netcdf(path=new_filenameT)
ds_trend.close()
#------------------------------------
clevsT=np.arange(-10.0,11.0,1)
color='RdBu_r'
exp='obs'
decomp='original'
lat,lon=climb.latlon(trend)
figtrend=myplot.oneplot_trend(variable=variable, decomp=decomp, exp=exp, model=model, clevs=clevsT, units=units, subdomain=subdomain, par=par, xlon=par[0][lon], ylat=par[0][lat], season=season, iyr=iyr, fyr=fyr,extent=False)
ofileT='spatialtrend_'+variable+'_'+model+'_'+exp+'_'+decomp+'_'+domain+'_'+season+'_'+str(iyr)+'_'+str(fyr)
figtrend.savefig(plotsDir+ofileT+'.png',format='png')
print('Figure save at ',plotsDir, 'as',ofileT)

ratio=10*abs(trend)/std_det


clevs=np.arange(0,1.5,0.01)
color='BuGn'
exp='obs'
decomp='original'
lat,lon=climb.latlon(trend)
figratio=myplot.oneplot_ds(variable=variable, decomp=decomp, exp=exp, model=model, clevs=clevs, color=color, units=units, subdomain=subdomain, ds=ratio, xlon=clim[lon], ylat=clim[lat], season=season, iyr=iyr, fyr=fyr,extent=False)

ofileR='ratio_trend-stddet_'+variable+'_'+model+'_'+exp+'_'+decomp+'_'+domain+'_'+season+'_'+str(iyr)+'_'+str(fyr)
figratio.savefig(plotsDir+ofileR+'.png',format='png')
print('Figure save at ',plotsDir, 'as',ofileR)
'''



##Computing time series of seasonal spatial average ----------------------------
print('Computing time series of seasonal spatial average')
#ts_season=climb.spatial_average(anoms)
wgts = np.cos(np.deg2rad(ylat))
tmp=np.ma.average(anoms,axis=1,weights=wgts)
ts_season=np.ma.average(tmp,axis=1)

plotnameTs='timeseries_anoms_'+variable+'_'+model+'_'+exp+'_'+decomp+'_'+subdomain+'_'+season+'_'+str(iyr)+'_'+str(fyr)
#plotnameTs='timeseries_vals'+variable+'_'+model4+'_'+exp+'_'+decomp+'_'+subdomain+'_'+season+'_'+str(iyr)+'_'+str(fyr)
np.savetxt(resultsDirTs+plotnameTs+'.txt',ts_season)
print(ts_season)
print('saving .txt at',resultsDirTs+plotnameTs)

##Trends---------------------------------------
xd=np.arange(0,len(ts_season),1)
par=stats.linregress(xd,ts_season)
print('stats',par)
parmed=stats.theilslopes(ts_season,xd,0.90)
print('theilsen',parmed)
parmk=mk.original_test(ts_season)
print('trend', 'h', 'p', 'z', 'Tau', 's', 'var_s', 'slope', 'intercept')
print(parmk)

trendsOf=[par,parmed,parmk]

fig, ax = plt.subplots(figsize=(12, 8))
my_ticks=xd+iyr
print(xd,my_ticks)
ticks=np.arange(iyr,fyr+1,1)
#ax.plot(xd,ts_season,color='b', linewidth=2)
ax.plot(xd, ts_season, 'b.')
ax.plot(xd, parmed[1] + parmed[0] * xd, 'r-',label="%s (%.2f)" % ('Theil-Sen', 10*parmed[0]))
ax.plot(xd, parmed[1] + parmed[2] * xd, 'r--')
ax.plot(xd, parmed[1] + parmed[3] * xd, 'r--')
ax.plot(xd, par[1] + par[0] * xd, 'g-',label="%s (%.2f)" % ('Lin reg', 10*par[0]))
ax.plot(xd, parmk[8] + parmk[7] * xd,'p-',label="%s (%.2f) sig %s" % ('MK', 10*parmk[7],parmk[0]))

#ax.axhline(0, color='k')
#ax.set_ylim(0, 3)
ax.set_xlabel('Year')
ax.set_ylabel('%s'%(units))
ax.set_title(plotnameTs)
freq=5
plt.xticks(xd,my_ticks, rotation='vertical')
plt.legend()
plt.savefig(resultsDirTs+plotnameTs+'.png',format='png')
#plt.show()
print('save at', resultsDirTs,plotnameTs)

textfile = open(resultsDirTs+'trends_'+plotnameTs, "w")
for element in trendsOf:
    textfile.write(str(element) + "\n")
textfile.close()
