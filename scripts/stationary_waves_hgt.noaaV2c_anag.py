'''
for season in DJF MAM JJA SON;do
for decomp in dyntot residtot;do
python stationary_waves_hgt.noaaV2c_anag.py ${season} ${decomp}
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

###
expname='obs.1851.2014.prmsl.noaaV2c.NAE.80.40.100_period1960.2014_deg2_int5'
domainX='global'
outputDir='obs.ltdeg2.prmsl_noaaV2c.f3'
sourceData='/home/msantolaria/Documents/MyResearch/ClimAnag/'
sourceDataAn=sourceData+outputDir+'/ncAnalogs/'
resultsDir=sourceDataAn
plotsDir=sourceData+outputDir+'/Plots/'
sourceOriginal='/home/msantolaria/Documents/Data/'


model='obs'
iyrAn=1851
fyrAn=2014
iyrAnom=1960
iyrpiC=iyrAn
fyrAnom=2014
fyrpiC=fyrAn
deg=2
exp='obs'

variableX='prmsl'
dataX='noaaV2c'
domainXAn='NAE'
Na=str(80)
Ns=str(40)
Nr=str(100)


variableY='z200'
dataY='noaaV2c'
domainY='global'

iyr=1960
fyr=2014
season=sys.argv[1]#'DJF'
decomp=sys.argv[2]

variable=variableY
data=dataY
domain=domainY
subdomain='NHeq'#sys.argv[3]

###Analog
fileName=variable+'.anom.'+decomp+'_'+data+'_'+model+'_'+domain+'_'+str(iyrpiC)+'01_'+str(fyrpiC)+'12_Na'+Na+'_Ns'+Ns+'_Nr'+Nr+'_'+variableX+'.'+dataX+'.'+domainXAn

dsAn= xr.open_dataset(sourceDataAn+fileName+'.nc')[variable]
latAn,lonAn=climb.latlon(dsAn)
units=dsAn.units
fieldAn=dom.field_dom(dsAn,subdomain)
ylatAn=fieldAn.coords[latAn]
xlonAn=fieldAn.coords[lonAn]


##Computing z*=z-[z]
zonal=fieldAn.mean(dim=lonAn)
z_star=fieldAn-zonal


if season[0]=='m':
    rmon=int(season.split('mon')[1])
    anoms,anomsAn=climb.monthly_selection(z_star,rmon,iyr,fyr)
else:
    anoms,anomsAn=climb.seasonal_selection(z_star,season,iyr,fyr)


###
#clim=vals.mean('time')
std=anoms.std('time')
#std_det=anoms_detrend.std('time')
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

##STD det
###
fig,axs= plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()},figsize=(8,10))
latS,latN,lonW,lonE,latlim,lonlim=dom.coord_domain(subdomain)

lat,lon=climb.latlon(std)
lons, lats = np.meshgrid(std[lon] ,std[lat])
CS1=axs.contourf(lons,lats, std,clevsStd,
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
ofileSdet='anoms_det_std_zstar200.'+decomp+'_'+model+'_'+exp+'_'+domain+'_'+season+'_'+str(iyr)+'_'+str(fyr)
plt.suptitle(ofileSdet,y=0.68)
fig.savefig(plotsDir+ofileSdet+'.png',format='png',bbox_inches='tight')
print('Figure save at ',plotsDir, 'as',ofileSdet)
plt.show()

