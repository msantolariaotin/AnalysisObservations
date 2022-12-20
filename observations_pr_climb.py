'''
for s in DJFMA JASO mon01 mon02 mon03 mon04 mon05 mon06 mon07 mon08 mon09 mon10 mon11 mon12;do
python3 observations_pr_climb.py ${s}
done
'''
# In[1]:


import sys
sys.path.insert(1, '/home/msantolaria/Documents/MyPythonLibrary/ClimAnag/')
import climbasis as climb
from climbasis import *
import domain as dom
import myplot
import glob


# In[2]:


source='/home/msantolaria/Documents/MyResearch/HMA-ClimAnalogs/Observations/'
resultsDir=source + 'Results/'
plotsDir=source + 'Plots/'


# In[6]:


sourceData='/home/msantolaria/Documents/Data/'
sourceData1=sourceData +'noaaV2c/'
sourceData2=sourceData +'era5/'
sourceData3=sourceData +'aphro/'
sourceData4=sourceData +'gpcc/'
sourceData5=sourceData +'cru/'
sourceData6=sourceData +'gpcp/'


# In[7]:


iyr=1982
fyr=2014
domain='HMA'
variable='pr'
units='mm/day'
season=sys.argv[1]


# In[8]:


model1='noaaV2c'
fileName1= 'prate.mon.mean.noaaV2c_185101-201412_2.0x2.0.nc'
ds1 = xr.open_dataset(sourceData1+fileName1)['prate']
print(model1,ds1.attrs)
print(ds1.units,'to')
ds1=3600*24*ds1
unit1='mm/day'
print(unit1)
field1=dom.field_dom(ds1,domain)


model3='aphro'
fileName3='APHRO_mon_MA_025deg_V1101_EXR1.1951-2015.nc'
ds3 = xr.open_dataset(sourceData3+fileName3)['precip']
print(model3,ds3.attrs)
field3=dom.field_dom(ds3,domain)
#print(ds3.attrs)


model4='gpcc'
fileName4='precip.mon.mean.nc'
ds4 = xr.open_dataset(sourceData4+fileName4)['precip']
print(model4,ds4.attrs)



model5='cru'
fileName5='cru_ts4.05.1901.2020.pre.dat.nc'
ds5 = xr.open_dataset(sourceData5+fileName5)['pre']
print(model5,ds5.attrs)
field05=dom.field_dom(ds5,domain)

month_length = field05.time.dt.days_in_month
#
field5=field05/month_length
units5='mm/day'
print(ds5.units,'to',units5)



model6='gpcp'
fileName6='precip_mon.gpcp_v2020_1891_2019_10.nc'
ds6 = xr.open_dataset(sourceData6+fileName6)['precip']
print(model6,ds6.attrs)
field06=dom.field_dom(ds6,domain)
month_length = field06.time.dt.days_in_month
#
field6=field06/month_length
units6='mm/day'
print(ds6.units,'to',units6)


# In[15]:


dsList=[field1,field3,field4,field5,field6]
modelList=[model1,model3,model4,model5,model6]


# In[16]:


valsList=[]
anomsList=[]

for elem in dsList:
    if season[0]=='m':
        rmon=int(season.split('mon')[1])
        vals,anoms=climb.monthly_selection(elem,rmon,iyr,fyr)
    else:
        if season[0]=='D':
            print('Winter',iyr,'-',iyr+1,fyr-1,'-',fyr)
            vals,anoms=climb.seasonal_selection(elem,season,6,iyr,6,fyr)
        else:
            vals,anoms=climb.seasonal_selection(elem,season,1,iyr,12,fyr)
    valsList.append(vals)
    anomsList.append(anoms)


# In[17]:


climList=[]
stdList=[]
for elem in valsList:
    c=elem.mean('time')
    s=elem.std('time')
    climList.append(c)
    stdList.append(s)



print(len(climList))# In[18]:

#climList[1].plot()


# In[19]:


parList=[]
trendList=[]
interceptList=[]
rvalueList=[]
pvalueList=[]
stderrList=[]
#---------------------------------------------------------
for elem in anomsList:
    par=climb.trend_vect(np.arange(0,elem.shape[0],1),elem,'time')
    parList.append(par)
    trendList.append(par[0])
    interceptList.append(par[1])
    rvalueList.append(par[2])
    pvalueList.append(par[3])
    stderrList.append(par[4])


# In[23]:


##Clim-----------------
#clevs=np.arange(0,14.25,0.25)
clevs=np.arange(0,8.25,0.25)
color='GnBu'
subdomain='HMA'
exp='obs'
decomp='original'
for i in range(len(climList)):
    lat,lon=climb.latlon(climList[i])
    figclim=myplot.oneplot_ds(variable=variable, decomp=decomp, exp=exp, model=modelList[i], clevs=clevs, color=color, units=units, subdomain=subdomain, ds=climList[i], xlon=climList[i][lon], ylat=climList[i][lat], season=season, iyr=iyr, fyr=fyr,extent=False)
    ofileC='clim_'+variable+'_'+modelList[i]+'_'+exp+'_'+decomp+'_'+subdomain+'_'+season+'_'+str(iyr)+'_'+str(fyr)
    figclim.savefig(plotsDir+ofileC+'.png',format='png')
    print('Figure save at ',plotsDir, 'as',ofileC)


# In[21]:


##STD-----------------
clevs=np.arange(0,3.5,0.1)
color='rainbow'
subdomain='HMA'
exp='obs'
decomp='original'
for i in range(len(stdList)):
    lat,lon=climb.latlon(stdList[i])
    figstd=myplot.oneplot_ds(variable=variable, decomp=decomp, exp=exp, model=modelList[i], clevs=clevs, color=color, units=units, subdomain=subdomain, ds=stdList[i], xlon=climList[i][lon], ylat=climList[i][lat], season=season, iyr=iyr, fyr=fyr,extent=False)
    ofileS='anoms_std_'+variable+'_'+modelList[i]+'_'+exp+'_'+decomp+'_'+subdomain+'_'+season+'_'+str(iyr)+'_'+str(fyr)
    figstd.savefig(plotsDir+ofileS+'.png',format='png')
    print('Figure save at ',plotsDir, 'as',ofileS)


# In[22]:


clevs=np.arange(-1.5,1.6,0.01)
#color='RdBu_r'
subdomain='HMA'
exp='obs'
decomp='original'
for i in range(len(parList)):
    lat,lon=climb.latlon(trendList[i])
    figtrend=myplot.oneplot_trend(variable=variable, decomp=decomp, exp=exp, model=modelList[i], clevs=clevs, units=units, subdomain=subdomain, par=parList[i], xlon=parList[i][0][lon], ylat=parList[i][0][lat], season=season, iyr=iyr, fyr=fyr,extent=False)
 #   figtrend=myplot.oneplot_trend(oneplot_trend(variable,decomp,exp,modelList[i],clevs,units,subdomain,parList[i],xlon,ylat,season,iyr,fyr))
    ofileT='spatialtrend_'+variable+'_'+modelList[i]+'_'+exp+'_'+decomp+'_'+subdomain+'_'+season+'_'+str(iyr)+'_'+str(fyr)
    figtrend.savefig(plotsDir+ofileT+'.png',format='png')
    print('Figure save at ',plotsDir, 'as',ofileT)


# In[ ]:





# In[ ]:




