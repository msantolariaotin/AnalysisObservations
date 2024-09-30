'''
for season in DJF MAM JJA SON;do
for subdomain in HK HM TP;do
python timeseries_tmp.cru.py ${decomp} ${season}
done
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


iyr=1982
fyr=2019

domain=sys.argv[2]
season=sys.argv[1]#'DJF'
exp='obs'

variable='tmp'
model='cru'
info=obs.get_obs(variable,model)
print(info)
fileName='cru_ts4_mon_1901.2020.tmp.dat.nc'
#fileName=info.get_obs(filename)
ds = xr.open_dataset(sourceData+model+'/'+fileName)[variable]
units=ds.units

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
print('stats',par)
parmed=stats.theilslopes(ts_season,xd,0.90)
print('theilsen',parmed)
parmk=mk.original_test(ts_season)
print('trend', 'h', 'p', 'z', 'Tau', 's', 'var_s', 'slope', 'intercept')
print(parmk)
trendsOf=[par,parmed,parmk]
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
#ax.plot(xd, parmk[8] + parmk[7] * xd,'p-',label="%s (%.2f) sig %s" % ('MK', 10*parmk[7],parmk[0]))
ax.plot(xd, par[1] + par[0] * xd, 'b-',label="%s (%.2f);%s (%.2f)" % ('Lin reg', 10*par[0],'p-val',par[3]))
#ax.axhline(0, color='k')
#ax.set_ylim(0, 3)
ax.set_xlabel('Year')
ax.set_ylabel('%s'%(units))
ax.set_title(plotnameTs)
freq=5
plt.xticks(xd,my_ticks, rotation='vertical')
plt.legend()
plt.savefig(plotsDir+plotnameTs+'.png',format='png')
plt.show()
print('save at', resultsDir,plotnameTs)
'''
textfile = open(resultsDirTs+'trends_'+plotnameTs, "w")
for element in trendsOf:
    textfile.write(str(element) + "\n")
textfile.close()

'''
