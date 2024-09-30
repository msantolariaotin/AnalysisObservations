#!/usr/bin/env python
# coding: utf-8

import sys
sys.path.insert(1, '/home/msantolaria/Documents/MyPythonLibrary/ClimAnag/')
import climbasis as climb
import domain as dom
import dynad as dyn
from climbasis import *
from dynad import *


outputDir='/home/msantolaria/Documents/MyResearch/HMA-ClimAnalogs/'
outputExp='EXP0.1-ESACCI-nan/'

resultsDir=outputDir+outputExp+'MonAnalogs/'
plotsDir=outputDir+outputExp+'Plots/'

sourceData='/home/msantolaria/Documents/Data/'

imon=1
fmon=12
ouputVarName='dyntot'
iyr=1990
fyr=2010

rmon=int(sys.argv[1])
year=1994

months=['01','02','03','04','05','06','07','08','09','10','11','12']
monthName='mon'+months[rmon-1]

Na= 10 #80
Ns= 4#40
Nr= 2

member='obs'

##SLP----------------------------
dataX='hadslp2'
filenameX = 'slp.mon_hadslp2_newtime_185001-201912.nc'
sourceDataX=sourceData+'HadSLP2/'
variableX='slp'

dsX = xr.open_dataset(sourceDataX+filenameX)[variableX]
domainX='MA'
#dsX=ds0X.sel(time=~ds0X.get_index("time").duplicated()) ## HadSLP2 has double definition time in some dates, e.x 1975-07ds

unitsX='hPa'

fieldX=dom.field_dom(dsX.sel(time=slice(str(iyr)+"-"+str(imon), str(fyr)+"-"+str(fmon))),domainX)

ylatX=fieldX.coords['latitude']
xlonX=fieldX.coords['longitude']
#fieldX=dom.field_dom(dsX,domain)
print('fieldX shape',fieldX.shape)

#------------------------------
domainY='HK'
dataY='esacci'
fileNameY='ESACCI-L3C_SNOW-SCFG-AVHRR_MERGED-fv1.0_HKH_gapfilled_icefilled_montlhy_1982-2014.nc'
sourceDataY=sourceData+'ESACCI/'
variableY='snc'
dsY = xr.open_dataset(sourceDataY+fileNameY)[variableY]
unitsY='%'
#print(ds0Y.units,'to',unitsY)

fieldY=dom.field_dom(dsY.sel(time=slice(str(iyr)+"-"+str(imon), str(fyr)+"-"+str(fmon))),domainY)
ylatY=fieldY.coords['lat']
xlonY=fieldY.coords['lon']

print('fieldY shape',fieldY.shape)


##PROBLEM
fieldY.sel(time='1994').groupby('time.month').mean('time').plot(col='month',col_wrap=4)
fieldY.sel(time='1995').groupby('time.month').mean('time').plot(col='month',col_wrap=4)

#fieldXdrop = fieldX.sel(time=slice('1994-10', '1995-01'))
#print(fieldXdrop)
#fieldYdrop = fieldY.sel(time=slice('1994-10', '1995-01'))
#fieldYdrop.groupby('time.month').mean('time').plot(col='month',col_wrap=4)


### Monthly selection for historical fields----------------
valsX,anomsX=climb.monthly_selection(fieldX,rmon)
valsY,anomsY=climb.monthly_selection(fieldY,rmon)

print('anomsX shape',anomsX.shape,'anomsY',anomsY.shape)


monthlyValuesX=anomsX.values
monthlyValuesY=anomsY.values

def climAnalogs_2fields_esacci(monthlyValuesX,monthlyValuesY,year,Na,Ns,Nr,iyr,fyr,uncertainty=False):
#"FieldX and FieldY are monthly fields from which to built the analogsto improve and put the original ds"
print('Compute climate analogs of month',year,'over period',iyr,'-',fyr)
print(Na,Ns,Nr,iyr,fyr)
print(" I.- Preparing the data for computation : matrix form time x m (m=nlatxnlon) \n ")

datax=monthlyValuesX
#if np.count_nonzero(np.isnan(datax))!=0:
#    print('ERROR: nan values to be masked ')
dx=monthlyValuesX.shape
nlonx=dx[2]
nlatx=dx[1]
ntx=dx[0]
nsx=nlonx*nlatx
X=datax.reshape(ntx,nsx)
print('X shape',X.shape)
##Preparing the data for computation: associated field
datay=monthlyValuesY
#if np.count_nonzero(np.isnan(datay))!=0:
#    print('ERROR: nan values to be masked ')
dy=monthlyValuesY.shape
nlony=dy[2]
nlaty=dy[1]
nty=dy[0]
nsy=nlony*nlaty
Y=datay.reshape(nty,nsy)
#  print("Fields X and Y are in matrix form with dimension: time x m (m=nlatxnlon) \n",'X:',X.shape,'Y:',Y.shape,"\n")
    ###
print('Y shape',Y.shape)



##Target field (Xt SLP jan 1979, first element time=0)
# print("II.- Target field \n")
# print("Xt",variable1,season,iyr,'\n')
countYear=year-iyr
print('Year',year)
print('CountYear',countYear)
Xt=X[countYear,0:nsx]
Yt=Y[countYear,0:nsy]
nyr=fyr-iyr+1

print("III.- Computing the Euclidean distance between target SLP field",rmon,"/",iyr,"among the monthly fields of all nyr=",nyr,'\n ... \n')
distList=[] ## Contains the values of the euclidean distance (r=2) between jan 1979 and all jan (1980-1998)
orderList=[]

for i in range(nyr):
    tmp=np.linalg.norm(X[countYear,0:nsx]-X[i,0:nsx])
    distList.append(tmp)
    orderList.append(i)
f=list(zip(distList,orderList))
 #   print("List of Euclidean Distance and Year Analog (0=iyr) in ascending order \n(Note: first element of list is =0,distance between target field and itself) \n")
f_order=sorted(f, key=lambda x: x[0]) ## x[0] is distList elements: Ordering monthly fields by ascending order of their Euclidean distance to the target field
#    print('Show first five analogs: list(Euclidean distance, counter)',f_order[0:5],'\n \n')
#    print('Taking closest Na=',Na,'analogs of the list Euc.Dist (1st element removed) \n')
allyears= [x[1] for x in f_order]
years_analogs_counter=[]
if rmon==1 and 1995-iyr in allyears:
    index=allyears.index(1995-iyr)
    print('Removing 01-1995',len(allyears))
    print(index,rmon,'removing')
    allyears.pop(index)
    print(len(allyears),Na)
years_analogs_counter= allyears[1:Na+1]
print('Year analogs are selected')
print(np.asarray(years_analogs_counter)+iyr,'\n \n') ##Just to see which real years are the analogs

print('IV.- Construction column vector of selected Na analogs for X and Y with Na x 1 dimension \n \n')
    ##Constructing column vector of slected Na analogs for X and Y
Xc_colVector_Na=[]
Yc_colVector_Na=[]
for i in years_analogs_counter:
    tmp1=X[int(i),0:nsx]
    tmp2=Y[int(i),0:nsy]
    Xc_colVector_Na.append(tmp1)
    Yc_colVector_Na.append(tmp2)
 #   print(len(Xc_colVector_Na),len(Yc_colVector_Na),'\n \n')


# In[42]:


print('V.-  Resample randomly Ns analogs from Na (col vector Ns x 1 ). Do this process Nr times \n \n')
#Pass the list to the first argument and the number of elements you want to get to the second argument. A list is returned
Ns=Ns
Nr=Nr
# print('Ns:',Ns,'Nr:',Nr,'\n \n')
Xc_colVector_Ns_Nr=[]
Yc_colVector_Ns_Nr=[]
for i in range(Nr):
    Xc_colVector_Ns,Yc_colVector_Ns=zip(*random.sample(list(zip(Xc_colVector_Na, Yc_colVector_Na)), Ns))
    Xc=np.asarray(Xc_colVector_Ns)
    Yc=np.asarray(Yc_colVector_Ns)
    #Compute the (Moore-Penrose) pseudo-inverse of a matrix to get beta
    pseudoMP=np.linalg.pinv(Xc.T) # Dimension of matrix Xc are opposite (only by construction of code, algebra doesn't change)
    beta=np.dot(pseudoMP,Xt)
    Xca=np.dot(Xc.T,beta)
    ## Beta coefficients are applied to Ns Y patterns
    Yca=np.dot(Yc.T,beta)
    Xc_colVector_Ns_Nr.append(Xca)
    Yc_colVector_Ns_Nr.append(Yca)

Xc_colVector_Ns_Nr_mean=np.mean(np.asarray(Xc_colVector_Ns_Nr),axis=0) ##Xc_colVector_Ns_Nr.shape =(Nr,Ns,nsx)
Yc_colVector_Ns_Nr_mean=np.mean(np.asarray(Yc_colVector_Ns_Nr),axis=0) ##Xc_colVector_Ns_Nr.shape =(Nr,Ns,nsx)

Xca_avg=np.asarray(Xc_colVector_Ns_Nr_mean)
Yca_avg=np.asarray(Yc_colVector_Ns_Nr_mean)

Xca_Nr=np.asarray(Xc_colVector_Ns_Nr)
Yca_Nr=np.asarray(Yc_colVector_Ns_Nr)



# In[44]:


print('VI.- Finally, convert them to a nlatxnlon grid. \n \n Return: Xh, Yh (nlatxnlon) mon',rmon,'year',year,'analogs')
Xh=Xca_avg.reshape(nlatx,nlonx)
Yh=Yca_avg.reshape(nlaty,nlony)

Xh_Nr=Xca_Nr.reshape(Nr,nlatx,nlonx)
Yh_Nr=Yca_Nr.reshape(Nr,nlaty,nlony)
#if uncertainty==True:
#     return Xh,Yh,Xh_Nr,Yh_Nr
#else:
#     return Xh,Yh

print(Yh)
'''
# In[54]:


##Let's compare field
print('mon',rmon,'year',year)
plt.pcolormesh(Xh,cmap='RdBu')
plt.colorbar()


# In[53]:


anomsX.sel(time=str(rmon)+'-'+str(year)).plot()


# In[55]:


##Let's compare field
print('mon',rmon,'year',year)
plt.pcolormesh(Yh,cmap='YlGnBu')
plt.colorbar()
'''

# In[13]:
'''


XhList_years=[]
YhList_years=[]

XhList_Nr_years=[]
YhList_Nr_years=[]

for year in np.arange(iyr,fyr+1,1):
    Xh,Yh,Xh_Nr,Yh_Nr=dyn.climAnalogs_2fields(anomsX.values,anomsY.values,year,Na,Ns,Nr,iyr,fyr,uncertainty=True)
    XhList_years.append(Xh)
    YhList_years.append(Yh)
    XhList_Nr_years.append(Xh_Nr)
    YhList_Nr_years.append(Yh_Nr)
    


# In[41]:


Xh_year=np.asarray(XhList_years)
print(Xh_year.shape)

Yh_year=np.asarray(YhList_years)
print(Yh_year.shape)


Xh_Nr_years=np.asarray(XhList_Nr_years)
print(Xh_Nr_years.shape)

Yh_Nr_years=np.asarray(YhList_Nr_years)
print(Yh_Nr_years.shape)


# In[15]:


##Saving monthly field (mon01 1851-20149. NOTE: dims are define 'latitude' and 'longitude'  to facilitate the spatial regrid for model comparison

daX= xr.DataArray(data=Xh_year, dims=["time","latitude", "longitude"],
                  coords=[np.arange(iyr,fyr+1,1),
                fieldX.coords[fieldX.dims[1]].values,
                fieldX.coords[fieldX.dims[2]].values])
daX.name= variableX
daX.attrs['long_name'] = variableX+' analog field-selected domain'
daX.attrs['units'] = unitsX
ofileX='%s.anom.%s_%s_%s_%s_%s_%s_%i_%i_Na%i_Ns%i_Nr%i'%(variableX,ouputVarName,monthName,dataX,exp,member,domainX,iyr,fyr,Na,Ns,Nr)
new_filename_x = resultsDir+ofileX+'.nc'
print ('saving to ', new_filename_x)
daX.to_netcdf(path=new_filename_x)
daX.close()

print ('finished saving')




# In[16]:


Yh_year_filled = np.where(np.isnan(Yh_year), 9999, Yh_year)


# In[17]:


daY= xr.DataArray(data=Yh_year_filled, dims=["time","latitude", "longitude"],
                  coords=[np.arange(iyr,fyr+1,1),
                fieldY.coords[fieldY.dims[1]].values,
                fieldY.coords[fieldY.dims[2]].values])
daY.name= variableY
daY.attrs['long_name'] = variableY+'analog field'
daY.attrs['units'] = unitsY
ofileY='%s.anom.%s_%s_%s_%s_%s_%s_%i_%i_Na%i_Ns%i_Nr%i'%(variableY,ouputVarName,monthName,dataY,exp,member,domainY,iyr,fyr,Na,Ns,Nr)
new_filename_Y = resultsDir+ofileY+'.nc'
print ('saving to ', new_filename_Y)
daY.to_netcdf(path=new_filename_Y,compute=True)
daY.close()

print ('finished saving')


# In[18]:




##Residual----
fieldO=dom.field_dom(dsY.sel(time=slice(str(iyr)+"-"+str(imon), str(fyr)+"-"+str(fmon))),domainY)

valsO,anomsO=climb.monthly_selection(fieldO,rmon)


# In[19]:


anomsO.sel(time='1995')


# In[17]:


anomsO_filled=np.where(np.isnan(anomsO), -9999, anomsO)


# In[18]:


ana='/home/maria/Documents/ResearchProjects/HMA-ClimAnalogs/EXP2-ObsOnly/MonAnalogs/snc.anom.dyntot_mon01_esacci_hadslp2.slp-esacci.snc_obs_HMA_1990_2000_Na10_Ns7_Nr2.nc'


# In[19]:


dsD=xr.open_dataset(ana)['snc']


# In[20]:


np.max(dsD)


# In[21]:


dsR=anomsO_filled-dsD.values


# In[22]:


np.max(dsR)


# In[23]:


np.min(dsR)


# In[33]:


dsRnan=np.where(abs(dsR) > 1000, np.nan, dsR)


# In[34]:


np.isnan(dsRnan).all()


# In[20]:



substraction=np.zeros((anomsO.shape))


# In[ ]:



'''
