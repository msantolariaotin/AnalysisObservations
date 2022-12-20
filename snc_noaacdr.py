'''
for i in DJFMA JASO;do
python snc_noaacdr.py ${i}
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


#resultsDir='/media/maria/MARIAEXT2/WinterTrendsObs/'
sourceData='/home/msantolaria/Documents/Data/'


# In[5]:


iyr=1982
imon=1
fmon=12
fyr=2014
domain='HMA'
season=sys.argv[1]
variable='snc'
units='%'

model='noaacdr'
fileName='nhsce_mon.v01r01_19661004_20210503_1.0x1.0.nc'
#Opening file
ds = xr.open_dataset(sourceData+model+'/'+fileName)['snow_cover_extent']
#print(ds2.units)
ds=100*ds
field=dom.field_dom(ds,domain)

if season[0]=='m':
    rmon=int(season.split('mon')[1])
    vals,anoms=climb.monthly_selection(field,rmon,iyr,fyr)
else:
    if season[0]=='D':
        print('Winter',iyr,'-',iyr+1,fyr-1,'-',fyr)
        vals,anoms=climb.seasonal_selection(field,season,6,iyr,6,fyr)
    else:
        vals,anoms=climb.seasonal_selection(field,season,1,iyr,12,fyr)


anoms_detrend=climb.detrend_dim(vals, 'time', deg=1)

clim=vals.mean('time')
std=vals.std('time')
std_det=anoms_detrend.std('time')
#--------------------------------
subdomain=domain


##Clim-----------------
#clevs=np.arange(0,14.25,0.25)
clevs=np.arange(0,8.25,0.25)
color='GnBu'
subdomain='HMA'
exp='obs'
decomp='original'
lat,lon=climb.latlon(clim)
figclim=myplot.oneplot_ds(variable=variable, decomp=decomp, exp=exp, model=model, clevs=clevs, color=color, units=units, subdomain=subdomain, ds=clim, xlon=clim[lon], ylat=clim[lat], season=season, iyr=iyr, fyr=fyr,extent=False)
ofileC='clim_'+variable+'_'+model+'_'+exp+'_'+decomp+'_'+subdomain+'_'+season+'_'+str(iyr)+'_'+str(fyr)
figclim.savefig(plotsDir+ofileC+'.png',format='png')
print('Figure save at ',plotsDir, 'as',ofileC)


# In[21]:


##STD-----------------
clevs=np.arange(0,3.5,0.1)
color='rainbow'
subdomain='HMA'
exp='obs'
decomp='original'
lat,lon=climb.latlon(std)
figstd=myplot.oneplot_ds(variable=variable, decomp=decomp, exp=exp, model=model, clevs=clevs, color=color, units=units, subdomain=subdomain, ds=std, xlon=clim[lon], ylat=clim[lat], season=season, iyr=iyr, fyr=fyr,extent=False)
ofileS='anoms_std_'+variable+'_'+model+'_'+exp+'_'+decomp+'_'+subdomain+'_'+season+'_'+str(iyr)+'_'+str(fyr)
figstd.savefig(plotsDir+ofileS+'.png',format='png')
print('Figure save at ',plotsDir, 'as',ofileS)

###MEAN DET
clevs=np.arange(0,3.5,0.1)
color='rainbow'
exp='obs'
decomp='original-detrend'
lat,lon=climb.latlon(std_det)
figstd_det=myplot.oneplot_ds(variable=variable, decomp=decomp, exp=exp, model=model, clevs=clevs, color=color, units=units, subdomain=subdomain, ds=std_det, xlon=clim[lon], ylat=clim[lat], season=season, iyr=iyr, fyr=fyr,extent=False)
ofileS='anoms_det_mean_'+variable+'_'+model+'_'+exp+'_'+decomp+'_'+subdomain+'_'+season+'_'+str(iyr)+'_'+str(fyr)
figstd_det.savefig(plotsDir+ofileS+'.png',format='png')
print('Figure save at ',plotsDir, 'as',ofileS)

#-----------------------------------
par=climb.trend_vect(anoms.time,anoms,'time')
trend=par[0]
intercept=par[1]
rvalue=par[2]
pvalue=par[3]
stderr=par[4]
#------------------------------------
clevs=np.arange(-1.5,1.6,0.01)
color='RdBu_r'
exp='obs'
decomp='original'
lat,lon=climb.latlon(trend)
figtrend=myplot.oneplot_trend(variable=variable, decomp=decomp, exp=exp, model=model, clevs=clevs, units=units, subdomain=subdomain, par=par, xlon=par[0][lon], ylat=par[0][lat], season=season, iyr=iyr, fyr=fyr,extent=False)
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
