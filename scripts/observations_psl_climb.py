#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
sys.path.insert(1, '/home/msantolaria/Documents/MyPythonLibrary/ClimAnag/')
import climbasis as climb
from climbasis import *
import domain as dom
import myplot
import glob


# In[2]:


from eofs.xarray import Eof


# In[3]:


source='/home/msantolaria/Documents/MyResearch/HMA-ClimAnalogs/Observations/'
resultsDir=source + 'Results/'
plotsDir=source + 'Plots/'
resultsDirTs=source+'TimeSeries/'

# In[4]:


#resultsDir='/media/maria/MARIAEXT2/WinterTrendsObs/'
sourceData='/home/msantolaria/Documents/Data/'


# In[5]:


#!ls /home/maria/Documents/Data/HadSLP2/


# In[6]:


variable='slp'
units='hPa'
domain='MA'
season='DJFMA'
iyr=1982
fyr=2014

exp='obs'
decomp='original'
subdomain='MA'
# In[23]:


model1='hadslp2'
filename1 = 'slp.mon_hadslp2_newtime_185001-201912.nc'
ds1 = xr.open_dataset(sourceData+model1+'/'+filename1)['slp']
print(ds1.units)
ds1 = ds1.assign_coords(longitude=(((ds1.longitude + 180) % 360) - 180))
ds1 = ds1.roll(longitude=int(len(ds1['longitude']) / 2), roll_coords=True)
field1=dom.field_dom(ds1,domain)


# In[24]:


model2=str('noaaV2c')
filename2 = 'prmsl_only.mon.mean.noaaV2c_185101-201412_2.0x2.0.nc'
ds2 = xr.open_dataset(sourceData+model2+'/'+filename2)['prmsl']
ds2 = ds2.assign_coords(lon=(((ds2.lon+ 180) % 360) - 180))
ds2 = ds2.roll(lon=int(len(ds2['lon']) / 2), roll_coords=True)
#print(ds1.units)
ds2=ds2/100
field2=dom.field_dom(ds2,domain)


# In[25]:


model3=str('noaaV3')
filename3 = 'prmsl_only.mon.mean.noaaV3_185101-201412_1.0x1.0.nc'
ds3 = xr.open_dataset(sourceData+model3+'/'+filename3)['prmsl']
ds3 = ds3.assign_coords(lon=(((ds3.lon+ 180) % 360) - 180))
ds3 = ds3.roll(lon=int(len(ds3['lon']) / 2), roll_coords=True)
#print(ds1.units)
ds3=ds3/100
field3=dom.field_dom(ds3,domain)


# In[26]:


model4=str('era5')
filename4='mslp_era5_NH_mon_1979-2020.nc'
ds4 = xr.open_dataset(sourceData+model4+'/'+filename4)['msl']
#print(ds1.units)
ds4=ds4/100
field4=dom.field_dom(ds4,domain)


# In[27]:


model5=str('eraint')
filename5='mslp.mon.eraint_197901_201512.nc'
ds5 = xr.open_dataset(sourceData+model5+'/'+filename5)['msl']
#print(ds1.units)
ds5 = ds5.assign_coords(longitude=(((ds5.longitude+ 180) % 360) - 180))
ds5 = ds5.roll(longitude=int(len(ds5['longitude']) / 2), roll_coords=True)
ds5=ds5/100
field5=dom.field_dom(ds5,domain)


# In[28]:


model6=str('ncep-ncar')
filename6='slp.mon.mean.ncep-ncar_194801-202203.nc'
ds6 = xr.open_dataset(sourceData+model6+'/'+filename6)['slp']
#print(ds1.units)
ds6 = ds6.assign_coords(lon=(((ds6.lon+ 180) % 360) - 180))
ds6 = ds6.roll(lon=int(len(ds6['lon']) / 2), roll_coords=True)
field6=dom.field_dom(ds6,domain)


# In[29]:


dsList=[field1,field2,field3,field4,field5,field6]
modelList=[model1,model2,model3,model4,model5,model6]


# In[30]:


valsList=[]
anomsList=[]

for elem in dsList:
    if season[0]=='D':
        print('Winter',iyr,'-',iyr+1,fyr-1,'-',fyr)
        vals,anoms=climb.seasonal_selection(elem,season,6,iyr,6,fyr)
    else:
        vals,anoms=climb.seasonal_selection(elem,season,1,iyr,12,fyr)
    valsList.append(vals)
    anomsList.append(anoms)


# In[31]:


climList=[]
stdList=[]
for elem in valsList:
    c=elem.mean('time')
    s=elem.std('time')
    climList.append(c)
    stdList.append(s)


# In[32]:
'''

parList=[]
trendList=[]
interceptList=[]
rvalueList=[]
pvalueList=[]
stderrList=[]
#---------------------------------------------------------
for elem in anomsList:
    par=climb.trend_vect(elem.time,elem,'time')
    parList.append(par)
    trendList.append(par[0])
    interceptList.append(par[1])
    rvalueList.append(par[2])
    pvalueList.append(par[3])
    stderrList.append(par[4])


# In[33]:


##Clim-----------------
clevs=np.arange(996,1034,2)
color='seismic'
subdomain='NAE'
exp='obs'
decomp='original'
for i in range(len(climList)):
    lat,lon=climb.latlon(climList[i])
    figclim=myplot.oneplot_ds(variable=variable, decomp=decomp, exp=exp, model=modelList[i], clevs=clevs, color=color, units=units, subdomain=subdomain, ds=climList[i], xlon=climList[i][lon], ylat=climList[i][lat], season=season, iyr=iyr, fyr=fyr)
    ofileC='clim_'+variable+'_'+modelList[i]+'_'+exp+'_'+decomp+'_'+domain+'_'+season+'_'+str(iyr)+'_'+str(fyr)
    figclim.savefig(plotsDir+ofileC+'.png',format='png')
    print('Figure save at ',plotsDir, 'as',ofileC)


# In[34]:


##STD-----------------
clevs=np.arange(0,4.2,0.2)
color='rainbow'
subdomain='NAE'
exp='obs'
decomp='original'
for i in range(len(stdList)):
    lat,lon=climb.latlon(stdList[i])
    figstd=myplot.oneplot_ds(variable=variable, decomp=decomp, exp=exp, model=modelList[i], clevs=clevs, color=color, units=units, subdomain=subdomain, ds=stdList[i], xlon=climList[i][lon], ylat=climList[i][lat], season=season, iyr=iyr, fyr=fyr)
    ofileS='anoms_std_'+variable+'_'+modelList[i]+'_'+exp+'_'+decomp+'_'+domain+'_'+season+'_'+str(iyr)+'_'+str(fyr)
    figstd.savefig(plotsDir+ofileS+'.png',format='png')
    print('Figure save at ',plotsDir, 'as',ofileS)


# In[35]:


clevs=np.arange(-1.5,1.6,0.01)
#color='RdBu_r'
subdomain='NAE'
exp='obs'
decomp='original'
for i in range(len(parList)):
    lat,lon=climb.latlon(trendList[i])
    figtrend=myplot.oneplot_trend(variable=variable, decomp=decomp, exp=exp, model=modelList[i], clevs=clevs, units=units, subdomain=subdomain, par=parList[i], xlon=parList[i][0][lon], ylat=parList[i][0][lat], season=season, iyr=iyr, fyr=fyr)
 #   figtrend=myplot.oneplot_trend(oneplot_trend(variable,decomp,exp,modelList[i],clevs,units,subdomain,parList[i],xlon,ylat,season,iyr,fyr))
    ofileT='spatialtrend_'+variable+'_'+modelList[i]+'_'+exp+'_'+decomp+'_'+domain+'_'+season+'_'+str(iyr)+'_'+str(fyr)
    figtrend.savefig(plotsDir+ofileT+'.png',format='png')
    print('Figure save at ',plotsDir, 'as',ofileT)

'''

##Computing time series of seasonal spatial average ----------------------------
print('Computing time series of seasonal spatial average')
tsList=[]
for i in np.arange(0,len(valsList),1):
    ts_season=climb.spatial_average(anomsList[i])
    #plotnameTs='timeseries_anoms_'+variable+'_'+modelList[i]+'_'+exp+'_'+decomp+'_'+subdomain+'_'+season+'_'+str(iyr)+'_'+str(fyr)
    plotnameTs='timeseries_anoms_'+variable+'_'+modelList[i]+'_'+exp+'_'+decomp+'_'+subdomain+'_'+season+'_'+str(iyr)+'_'+str(fyr)
    np.savetxt(resultsDirTs+plotnameTs+'.txt',ts_season)
    print(ts_season)
    tsList.append(ts_season)
    print('saving .txt at',resultsDirTs+plotnameTs)
# Plot the leading PC time series.
    fig, ax = plt.subplots(figsize=(12, 8))
    xd=np.arange(0,len(ts_season),1)
    my_ticks=xd+iyr
    print(xd,my_ticks)
#my_ticks=np.arange(iyr,fyr+1,1)
    ax.plot(xd,ts_season,color='b', linewidth=2)
    #ax.axhline(0, color='k')
    #ax.set_ylim(-0.5,0.5)
    ax.set_xlabel('Year')
    ax.set_ylabel('Vals %s'%(units))
    ax.set_title(plotnameTs)
    freq=5
    plt.xticks(xd,my_ticks, rotation='vertical')
    plt.savefig(resultsDirTs+plotnameTs+'.png',format='png')
    #plt.show()
    print('save at', resultsDir,plotnameTs)


