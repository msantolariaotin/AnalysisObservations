import sys
sys.path.insert(1, '/home/maria/Documents/MyPythonLibrary/')
import climbasis as climb
import domain as dom
import numpy as np
import random

def climAnalogs_2fields_ref(monthlyValuesX,monthlyValuesY,monthlyValuesXAn,domainXAn,year,Na,Ns,Nr,iyrAn,fyrAn,interval,uncertainty=False):
    '''
    Computed dynamical constructed analogs
    Parameters:
    -monthlyValuesX: dynamical field as target (SLP, geopotential...)
    -monthlyValuesY: associated field to be dynamically constructed ( temperature,snow cover, precipitation...)
    -monthlyValuesXref: dynamical field to pick/select the analogs (SLP piControl, SLP detrended, SLP original)
    -year: target year
    -Na,Ns,Nr: number of analogs, subsample, iterations for random subsample
    -iyrAn,fyrAn: period of choosing analogs, it can be different from field period anomalies( iyrAnom,fyrAnom)
    -interval: Interval of years around the target mon/year that are not considered
    -uncertainty: in case you want to save the random subsample to make a Montecarlo test for statistical significante (see O'Reilly et al. 2017 for details)
    
    >Application:
    Following Deser et al. 2016 methodology for CESM-LE large ensemble and NOAAv2c reanalysis:
    ##Dynamical adjustment
    Xh_year=[]
    Yh_year=[]
    for year in np.arange(iyrAnom,fyrAnom+1,1):
        Xh,Yh=climAnalogs_2fields_ref(monthlyValuesX=anomsX,monthlyValuesY=anomsYpiC,monthlyValuesXAn=anomsXpiC,
                                      domainXAn=domainXAn,year=year,Na=Na,Ns=Ns,Nr=Nr,
                                      iyrAn=iyrpiC,fyrAn=fyrpiC,interval=1,uncertainty=False)

        Xh_year.append(Xh)
        Yh_year.append(Yh)
    '''
    #print(Na,Ns,Nr,iyr,fyr)
    #print(" I.- Preparing the data for computation : matrix form time x m (m=nlatxnlon) \n ")
    datax=monthlyValuesX.values; dx=monthlyValuesX.shape;nlonx=dx[2];nlatx=dx[1];ntx=dx[0];nsx=nlonx*nlatx
    X=datax.reshape(ntx,nsx)
    
    datay=monthlyValuesY.values;dy=monthlyValuesY.shape;nlony=dy[2];nlaty=dy[1];nty=dy[0];nsy=nlony*nlaty
    Y=datay.reshape(nty,nsy)
    
    dataxAn=monthlyValuesXAn.values; dxAn=monthlyValuesXAn.shape;nlonxAn=dxAn[2];nlatxAn=dxAn[1];ntxAn=dxAn[0];nsxAn=nlonxAn*nlatxAn
    XAn=dataxAn.reshape(ntxAn,nsxAn)
  #  print(X.shape) 
    ##II.-Target field (Xt SLP jan 1979)
    #print(year)
    Xt_xarray=monthlyValuesX.sel(time=str(year))[0,:,:]
    Xt=Xt_xarray.values.reshape(nsx)
    
    Yt_xarray=monthlyValuesY.sel(time=str(year))[0,:,:]
    Yt=Yt_xarray.values.reshape(nsy)
    #
    ##III.- Euclidean distance
    distList=[] ## Contains the values of the euclidean distance (r=2) between jan 1979 and all jan (1980-1998)
    orderList=[]
    for i in np.arange(iyrAn,fyrAn+1,1):
        tmp= np.linalg.norm(dom.field_dom(monthlyValuesX.sel(time=str(year)),domainXAn).values-dom.field_dom(monthlyValuesXAn.sel(time=str(i)),domainXAn).values)
        #print(len(dom.field_dom(monthlyValuesX.sel(time=str(year)),domainXAn).values))
        distList.append(tmp)
        orderList.append(i)
    f=list(zip(distList,orderList))
    f_order=sorted(f, key=lambda x: x[0]) ## x[0] is distList elements: Ordering monthly fields by ascending order of their Euclidean distance to the target field
    ##Taking only the year analogs
    allyears= [x[1] for x in f_order]
    interval=interval
    forbiddenyears=(np.arange(year-interval,year+interval+1,1)).tolist()
    #print('Year',year,'Filter',forbiddenyears) # In this way we remove already the year itself
    allyears_filter=[x for x in allyears if x not in forbiddenyears]
    #print(allyears_filter)
    years_analogs_counter=[]
    years_analogs_counter= allyears_filter[:Na]
    #print(years_analogs_counter)
    #
    #IV.- Construction column vector of selected Na analogs for X and Y with Na x 1 dimension \n \n')
    ##Constructing column vector of slected Na analogs for X and Y
    Xc_colVector_Na=[]
    Yc_colVector_Na=[]
    for i in years_analogs_counter:
        tmp1=monthlyValuesXAn.sel(time=str(i))[0,:,:]#.values
        tmp2=monthlyValuesY.sel(time=str(i))[0,:,:] #extra dimension when selecting xarray
        Xc_colVector_Na.append(tmp1)
        Yc_colVector_Na.append(tmp2)
    #print(len(Xc_colVector_Na),len(Yc_colVector_Na),'\n \n')
    #V.-  Resample randomly Ns analogs from Na (col vector Ns x 1 ). Do this process Nr times \n \n')
    #Pass the list to the first argument and the number of elements you want to get to the second argument. A list is returned
    Ns=Ns;Nr=Nr
    #print('Ns:',Ns,'Nr:',Nr,'\n \n')
    Xc_colVector_Ns_Nr=[]
    Yc_colVector_Ns_Nr=[]
    for i in range(Nr):
        Xc_colVector_Ns,Yc_colVector_Ns=zip(*random.sample(list(zip(Xc_colVector_Na, Yc_colVector_Na)), Ns))
        Xc_array=np.asarray(Xc_colVector_Ns)
        Xc=Xc_array.reshape(Ns,nsx)
        Yc_array=np.asarray(Yc_colVector_Ns)
        Yc=Yc_array.reshape(Ns,nsy)
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

    #VI.- Finally, convert them to a nlatxnlon grid. \n \n Return: Xh, Yh (nlatxnlon) january 1979 analogs
    Xh=Xca_avg.reshape(nlatx,nlonx)
    Yh=Yca_avg.reshape(nlaty,nlony)
    Xh_Nr=Xca_Nr.reshape(Nr,nlatx,nlonx)
    Yh_Nr=Yca_Nr.reshape(Nr,nlaty,nlony)
    if uncertainty==True:
        return Xh,Yh,Xh_Nr,Yh_Nr
    else:
        return Xh,Yh

###### 
def climAnalogs_4fields_ref(monthlyValuesX,monthlyValuesY,monthlyValuesY2,monthlyValuesY3,monthlyValuesXAn,domainXAn,year,Na,Ns,Nr,iyrAn,fyrAn,interval,uncertainty=False):
    '''
    Computed dynamical constructed analogs
    Parameters:
    -monthlyValuesX: dynamical field as target (SLP, geopotential...)
    -monthlyValuesY: associated field to be dynamically constructed ( temperature,snow cover, precipitation...)
    -monthlyValuesXref: dynamical field to pick/select the analogs (SLP piControl, SLP detrended, SLP original)
    -year: target year
    -Na,Ns,Nr: number of analogs, subsample, iterations for random subsample
    -iyrAn,fyrAn: period of choosing analogs, it can be different from field period anomalies( iyrAnom,fyrAnom)
    -interval: Interval of years around the target mon/year that are not considered
    -uncertainty: in case you want to save the random subsample to make a Montecarlo test for statistical significante (see O'Reilly et al. 2017 for details)
    
    >Application:
    Following Deser et al. 2016 methodology for CESM-LE large ensemble and NOAAv2c reanalysis:
    ##Dynamical adjustment
    Xh_year=[]
    Yh_year=[]
    for year in np.arange(iyrAnom,fyrAnom+1,1):
        Xh,Yh=climAnalogs_2fields_ref(monthlyValuesX=anomsX,monthlyValuesY=anomsYpiC,monthlyValuesXAn=anomsXpiC,
                                      domainXAn=domainXAn,year=year,Na=Na,Ns=Ns,Nr=Nr,
                                      iyrAn=iyrpiC,fyrAn=fyrpiC,interval=1,uncertainty=False)

        Xh_year.append(Xh)
        Yh_year.append(Yh)
    '''
    #print(Na,Ns,Nr,iyr,fyr)
    #print(" I.- Preparing the data for computation : matrix form time x m (m=nlatxnlon) \n ")
    datax=monthlyValuesX.values; dx=monthlyValuesX.shape;nlonx=dx[2];nlatx=dx[1];ntx=dx[0];nsx=nlonx*nlatx
    X=datax.reshape(ntx,nsx)
    
    datay=monthlyValuesY.values;dy=monthlyValuesY.shape;nlony=dy[2];nlaty=dy[1];nty=dy[0];nsy=nlony*nlaty
    Y=datay.reshape(nty,nsy)
    
    datay2=monthlyValuesY2.values;dy2=monthlyValuesY2.shape;nlony2=dy2[2];nlaty2=dy2[1];nty2=dy2[0];nsy2=nlony2*nlaty2
    Y2=datay2.reshape(nty2,nsy2)
        
    datay3=monthlyValuesY3.values;dy3=monthlyValuesY3.shape;nlony3=dy3[2];nlaty3=dy3[1];nty3=dy3[0];nsy3=nlony3*nlaty3
    Y3=datay3.reshape(nty3,nsy3)
    
    dataxAn=monthlyValuesXAn.values; dxAn=monthlyValuesXAn.shape;nlonxAn=dxAn[2];nlatxAn=dxAn[1];ntxAn=dxAn[0];nsxAn=nlonxAn*nlatxAn
    XAn=dataxAn.reshape(ntxAn,nsxAn)
    print(X.shape)
     
    ##II.-Target field (Xt SLP jan 1979)
    #print(year)
    Xt_xarray=monthlyValuesX.sel(time=str(year))[0,:,:]
    Xt=Xt_xarray.values.reshape(nsx)
    
    Yt_xarray=monthlyValuesY.sel(time=str(year))[0,:,:]
    Yt=Yt_xarray.values.reshape(nsy)
    #
    ##III.- Euclidean distance
    distList=[] ## Contains the values of the euclidean distance (r=2) between jan 1979 and all jan (1980-1998)
    orderList=[]
    for i in np.arange(iyrAn,fyrAn+1,1):
        tmp= np.linalg.norm(dom.field_dom(monthlyValuesX.sel(time=str(year)),domainXAn).values-dom.field_dom(monthlyValuesXAn.sel(time=str(i)),domainXAn).values)
        #print(len(dom.field_dom(monthlyValuesX.sel(time=str(year)),domainXAn).values))
        distList.append(tmp)
        orderList.append(i)
    f=list(zip(distList,orderList))
    f_order=sorted(f, key=lambda x: x[0]) ## x[0] is distList elements: Ordering monthly fields by ascending order of their Euclidean distance to the target field
    ##Taking only the year analogs
    allyears= [x[1] for x in f_order]
    interval=interval
    forbiddenyears=(np.arange(year-interval,year+interval+1,1)).tolist()
    #print('Year',year,'Filter',forbiddenyears) # In this way we remove already the year itself
    allyears_filter=[x for x in allyears if x not in forbiddenyears]
    #print(allyears_filter)
    years_analogs_counter=[]
    years_analogs_counter= allyears_filter[:Na]
    #print(years_analogs_counter)
    #
    #IV.- Construction column vector of selected Na analogs for X and Y with Na x 1 dimension \n \n')
    ##Constructing column vector of slected Na analogs for X and Y
    Xc_colVector_Na=[]
    Yc_colVector_Na=[]
    Y2c_colVector_Na=[]
    Y3c_colVector_Na=[]
    for i in years_analogs_counter:
        tmp0=monthlyValuesXAn.sel(time=str(i))[0,:,:]#.values
        tmp1=monthlyValuesY.sel(time=str(i))[0,:,:] #extra dimension when selecting xarray
        tmp2=monthlyValuesY2.sel(time=str(i))[0,:,:] #extra dimension when selecting xarray
        tmp3=monthlyValuesY3.sel(time=str(i))[0,:,:] #extra dimension when selecting xarray
        Xc_colVector_Na.append(tmp0)
        Yc_colVector_Na.append(tmp1)
        Y2c_colVector_Na.append(tmp2)
        Y3c_colVector_Na.append(tmp3)
    #print(len(Xc_colVector_Na),len(Yc_colVector_Na),'\n \n')
    #V.-  Resample randomly Ns analogs from Na (col vector Ns x 1 ). Do this process Nr times \n \n')
    #Pass the list to the first argument and the number of elements you want to get to the second argument. A list is returned
    Ns=Ns;Nr=Nr
    #print('Ns:',Ns,'Nr:',Nr,'\n \n')
    Xc_colVector_Ns_Nr=[]
    Yc_colVector_Ns_Nr=[]
    Y2c_colVector_Ns_Nr=[]
    Y3c_colVector_Ns_Nr=[]
    for i in range(Nr):
        Xc_colVector_Ns,Yc_colVector_Ns,Y2c_colVector_Ns,Y3c_colVector_Ns=zip(*random.sample(list(zip(Xc_colVector_Na, Yc_colVector_Na,Y2c_colVector_Na,Y3c_colVector_Na)), Ns))
        Xc_array=np.asarray(Xc_colVector_Ns)
        Xc=Xc_array.reshape(Ns,nsx)
        Yc_array=np.asarray(Yc_colVector_Ns)
        Yc=Yc_array.reshape(Ns,nsy)
        Y2c_array=np.asarray(Y2c_colVector_Ns)
        Y2c=Y2c_array.reshape(Ns,nsy2)
        Y3c_array=np.asarray(Y3c_colVector_Ns)
        Y3c=Y3c_array.reshape(Ns,nsy3)
        #Compute the (Moore-Penrose) pseudo-inverse of a matrix to get beta
        pseudoMP=np.linalg.pinv(Xc.T) # Dimension of matrix Xc are opposite (only by construction of code, algebra doesn't change)
        beta=np.dot(pseudoMP,Xt)
        Xca=np.dot(Xc.T,beta)
        ## Beta coefficients are applied to Ns Y patterns
        Yca=np.dot(Yc.T,beta)
        Y2ca=np.dot(Y2c.T,beta)
        Y3ca=np.dot(Y3c.T,beta)
        Xc_colVector_Ns_Nr.append(Xca)
        Yc_colVector_Ns_Nr.append(Yca)
        Y2c_colVector_Ns_Nr.append(Y2ca)
        Y3c_colVector_Ns_Nr.append(Y3ca)

    Xc_colVector_Ns_Nr_mean=np.mean(np.asarray(Xc_colVector_Ns_Nr),axis=0) ##Xc_colVector_Ns_Nr.shape =(Nr,Ns,nsx)
    Yc_colVector_Ns_Nr_mean=np.mean(np.asarray(Yc_colVector_Ns_Nr),axis=0) ##Xc_colVector_Ns_Nr.shape =(Nr,Ns,nsx)
    Y2c_colVector_Ns_Nr_mean=np.mean(np.asarray(Y2c_colVector_Ns_Nr),axis=0) ##Xc_colVector_Ns_Nr.shape =(Nr,Ns,nsx)
    Y3c_colVector_Ns_Nr_mean=np.mean(np.asarray(Y3c_colVector_Ns_Nr),axis=0) ##Xc_colVector_Ns_Nr.shape =(Nr,Ns,nsx)

    Xca_avg=np.asarray(Xc_colVector_Ns_Nr_mean)
    Yca_avg=np.asarray(Yc_colVector_Ns_Nr_mean)
    Y2ca_avg=np.asarray(Y2c_colVector_Ns_Nr_mean)
    Y3ca_avg=np.asarray(Y3c_colVector_Ns_Nr_mean)

    Xca_Nr=np.asarray(Xc_colVector_Ns_Nr)
    Yca_Nr=np.asarray(Yc_colVector_Ns_Nr)
    Y2ca_Nr=np.asarray(Y2c_colVector_Ns_Nr)
    Y3ca_Nr=np.asarray(Y3c_colVector_Ns_Nr)

    #VI.- Finally, convert them to a nlatxnlon grid. \n \n Return: Xh, Yh (nlatxnlon) january 1979 analogs
    Xh=Xca_avg.reshape(nlatx,nlonx)
    Yh=Yca_avg.reshape(nlaty,nlony)
    Y2h=Y2ca_avg.reshape(nlaty2,nlony2)
    Y3h=Y3ca_avg.reshape(nlaty3,nlony3)
    
    Xh_Nr=Xca_Nr.reshape(Nr,nlatx,nlonx)
    Yh_Nr=Yca_Nr.reshape(Nr,nlaty,nlony)
    Y2h_Nr=Y2ca_Nr.reshape(Nr,nlaty2,nlony2)
    Y3h_Nr=Y3ca_Nr.reshape(Nr,nlaty3,nlony3)
    
    if uncertainty==True:
        return Xh,Yh,Y2h,Y3h,Xh_Nr,Yh_Nr,Y2h_Nr,Y3h_Nr
    else:
        return Xh,Yh,Y2h,Y3h

########
###ADAPTATIONS:
########
#
########
### Using 4 associated Y fields with exactly the same random subsample
########
def climAnalogs_4fields(monthlyValuesX,monthlyValuesY,monthlyValuesY2,monthlyValuesY3,year,Na,Ns,Nr,iyr,fyr,uncertainty=False):
    "FieldX and FieldY,Y2,Y3 are monthly fields from which to built the analogsto improve and put the original ds"
    #print('Compute climate analogs',year,'over period',iyr,'-',fyr)
    #print(Na,Ns,Nr,iyr,fyr)
 #   print(" I.- Preparing the data for computation : matrix form time x m (m=nlatxnlon) \n ")
    datax=monthlyValuesX
    dx=monthlyValuesX.shape
    nlonx=dx[2]
    nlatx=dx[1]
    ntx=dx[0]
    nsx=nlonx*nlatx
    X=datax.reshape(ntx,nsx)
    print(X.shape)
    ##Preparing the data for computation: associated field Y
    datay=monthlyValuesY
    dy=monthlyValuesY.shape
    nlony=dy[2]
    nlaty=dy[1]
    nty=dy[0]
    nsy=nlony*nlaty
    Y=datay.reshape(nty,nsy)
    ##Preparing the data for computation: associated field Y2
    datay2=monthlyValuesY2
    dy2=monthlyValuesY2.shape
    nlony2=dy2[2]
    nlaty2=dy2[1]
    nty2=dy2[0]
    nsy2=nlony2*nlaty2
    Y2=datay2.reshape(nty2,nsy2)
     ##Preparing the data for computation: associated field Y3
    datay3=monthlyValuesY3
    dy3=monthlyValuesY3.shape
    nlony3=dy3[2]
    nlaty3=dy3[1]
    nty3=dy3[0]
    nsy3=nlony3*nlaty3
    Y3=datay3.reshape(nty3,nsy3)

  #  print("Fields X and Y are in matrix form with dimension: time x m (m=nlatxnlon) \n",'X:',X.shape,'Y:',Y.shape,"\n")
    ###
    ##Target field (Xt SLP jan 1979, first element time=0)
   # print("II.- Target field \n")
   # print("Xt",variable1,season,iyr,'\n')
    countYear=year-iyr
    print('Year',year)
    print('CountYear',countYear)
    Xt=X[countYear,0:nsx]
    Yt=Y[countYear,0:nsy]
    Y2t=Y2[countYear,0:nsy2]
    Y3t=Y3[countYear,0:nsy3]
    distList=[] ## Contains the values of the euclidean distance (r=2) between jan 1979 and all jan (1980-1998)
    orderList=[]
  #  print("III.- Computing the Euclidean distance between target SLP field",season,"/",iyr,"among the monthly fields of all nyr=",nyr,'\n ... \n')
    nyr=fyr-iyr+1
    for i in range(nyr):
        tmp=np.linalg.norm(X[countYear,0:nsx]-X[i,0:nsx])
        distList.append(tmp)
        orderList.append(i)
    f=list(zip(distList,orderList))
 #   print("List of Euclidean Distance and Year Analog (0=iyr) in ascending order \n(Note: first element of list is =0,distance between target field and itself) \n")
    f_order=sorted(f, key=lambda x: x[0]) ## x[0] is distList elements: Ordering monthly fields by ascending order of their Euclidean distance to the target field
#    print('Show first five analogs: list(Euclidean distance, counter)',f_order[0:5],'\n \n')
#    print('Taking closest Na=',Na,'analogs of the list Euc.Dist (1st element removed) \n')
    Na=Na
    f_ca=f_order[1:Na+1] #Remove the first field as it is itself with dist=0.0, always in the first position (to improve by removing also the five years aroundm maybe)
    ## Taking only the year analogs
    years_analogs_counter= [x[1] for x in f_ca] ##x[1] taking the 2nd element of the list (counter), the year analogs
    #type(years_analogs)
 #   print('Year analogs are selected')
 #   print(np.asarray(years_analogs_counter)+iyr,'\n \n') ##Just to see which real years are the analogs
 #   print('IV.- Construction column vector of selected Na analogs for X and Y with Na x 1 dimension \n \n')
    ##Constructing column vector of slected Na analogs for X and Y
    Xc_colVector_Na=[]
    Yc_colVector_Na=[]
    Y2c_colVector_Na=[]
    Y3c_colVector_Na=[]
    for i in years_analogs_counter:
        tmp0=X[int(i),0:nsx]
        tmp1=Y[int(i),0:nsy]
        tmp2=Y2[int(i),0:nsy2]
        tmp3=Y3[int(i),0:nsy3]
        Xc_colVector_Na.append(tmp0)
        Yc_colVector_Na.append(tmp1)
        Y2c_colVector_Na.append(tmp2)
        Y3c_colVector_Na.append(tmp3)
#   print(len(Xc_colVector_Na),len(Yc_colVector_Na),'\n \n')
 #   print('V.-  Resample randomly Ns analogs from Na (col vector Ns x 1 ). Do this process Nr times \n \n')
    #Pass the list to the first argument and the number of elements you want to get to the second argument. A list is returned
    Ns=Ns
    Nr=Nr
   # print('Ns:',Ns,'Nr:',Nr,'\n \n')
    Xc_colVector_Ns_Nr=[]
    Yc_colVector_Ns_Nr=[]
    Y2c_colVector_Ns_Nr=[]
    Y3c_colVector_Ns_Nr=[]
    for i in range(Nr):
        Xc_colVector_Ns,Yc_colVector_Ns,Y2c_colVector_Ns,Y3c_colVector_Ns=zip(*random.sample(list(zip(Xc_colVector_Na, Yc_colVector_Na,Y2c_colVector_Na,Y3c_colVector_Na)), Ns))
        Xc=np.asarray(Xc_colVector_Ns)
        Yc=np.asarray(Yc_colVector_Ns)
        Y2c=np.asarray(Y2c_colVector_Ns)
        Y3c=np.asarray(Y3c_colVector_Ns)
        #Compute the (Moore-Penrose) pseudo-inverse of a matrix to get beta
        pseudoMP=np.linalg.pinv(Xc.T) # Dimension of matrix Xc are opposite (only by construction of code, algebra doesn't change)
        beta=np.dot(pseudoMP,Xt)
        Xca=np.dot(Xc.T,beta)
        ## Beta coefficients are applied to Ns Y patterns
        Yca=np.dot(Yc.T,beta)
        Y2ca=np.dot(Y2c.T,beta)
        Y3ca=np.dot(Y3c.T,beta)
        Xc_colVector_Ns_Nr.append(Xca)
        Yc_colVector_Ns_Nr.append(Yca)
        Y2c_colVector_Ns_Nr.append(Y2ca)
        Y3c_colVector_Ns_Nr.append(Y3ca)
    Xc_colVector_Ns_Nr_mean=np.mean(np.asarray(Xc_colVector_Ns_Nr),axis=0) ##Xc_colVector_Ns_Nr.shape =(Nr,Ns,nsx)
    Yc_colVector_Ns_Nr_mean=np.mean(np.asarray(Yc_colVector_Ns_Nr),axis=0) ##Xc_colVector_Ns_Nr.shape =(Nr,Ns,nsx)
    Y2c_colVector_Ns_Nr_mean=np.mean(np.asarray(Y2c_colVector_Ns_Nr),axis=0)
    Y3c_colVector_Ns_Nr_mean=np.mean(np.asarray(Y3c_colVector_Ns_Nr),axis=0)
    ##Average of "best" linear combinations
    Xca_avg=np.asarray(Xc_colVector_Ns_Nr_mean)
    Yca_avg=np.asarray(Yc_colVector_Ns_Nr_mean)
    Y2ca_avg=np.asarray(Y2c_colVector_Ns_Nr_mean)
    Y3ca_avg=np.asarray(Y3c_colVector_Ns_Nr_mean)
    
    ## Convert col vector of Nr into an array to reshape
    Xca_Nr=np.asarray(Xc_colVector_Ns_Nr)
    Yca_Nr=np.asarray(Yc_colVector_Ns_Nr)
    Y2ca_Nr=np.asarray(Y2c_colVector_Ns_Nr)
    Y3ca_Nr=np.asarray(Y3c_colVector_Ns_Nr)
    
    # print('VIII.- Finally, convert them to a nlatxnlon grid. \n \n Return: Xh, Yh (nlatxnlon) january 1979 analogs')
    Xh=Xca_avg.reshape(nlatx,nlonx)
    Yh=Yca_avg.reshape(nlaty,nlony)
    Y2h=Y2ca_avg.reshape(nlaty2,nlony2)
    Y3h=Y3ca_avg.reshape(nlaty3,nlony3)
    
    Xh_Nr=Xca_Nr.reshape(Nr,nlatx,nlonx)
    Yh_Nr=Yca_Nr.reshape(Nr,nlaty,nlony)
    Y2h_Nr=Y2ca_Nr.reshape(Nr,nlaty2,nlony2)
    Y3h_Nr=Y3ca_Nr.reshape(Nr,nlaty3,nlony3)
    
    if uncertainty==True:
        return Xh,Yh,Y2h,Y3h,Xh_Nr,Yh_Nr,Y2h_Nr,Y3h_Nr
    else:
        return Xh,Yh,Y2h,Y3h
 
####
##FIRST VERSION Deser et al. 2016: it is preferable to use the function climAnalogs_2fields_ref (at the top)
####
#Using piControl---------------------------------------------------------------------------------------------------------------------------------

def climAnalogs_piC_2fields(monthlyValuesXpic,monthlyValuesX,monthlyValuesYpic,year,Na,Ns,Nr,iyrpiC,fyrpiC):
 #  I.- Preparing the data for computation : matrix form time x m (m=nlatxnlon)
    print(Na,Ns,Nr)
    dataxpic=monthlyValuesXpic.values; dxpic=monthlyValuesXpic.shape
    nlonxpic=dxpic[2];nlatxpic=dxpic[1];ntxpic=dxpic[0];nsxpic=nlonxpic*nlatxpic
    Xpic=dataxpic.reshape(ntxpic,nsxpic)

    datax=monthlyValuesX.values; dx=monthlyValuesX.shape;nlonx=dx[2];nlatx=dx[1];ntx=dx[0];nsx=nlonx*nlatx
    X=datax.reshape(ntx,nsx)

    dataypic=monthlyValuesYpic.values;dypic=monthlyValuesYpic.shape
    nlonypic=dypic[2];nlatypic=dypic[1];ntypic=dypic[0];nsypic=nlonypic*nlatypic
    Ypic=dataypic.reshape(ntypic,nsypic)

   # "II.- Target field
    Xt_xarray=monthlyValuesX.sel(time=str(year))[0,:,:]
    Xt=Xt_xarray.values.reshape(nsx)

    ##III.- Euclidean distance
    distList=[] ## Contains the values of the euclidean distance 
    orderList=[] ## Contains the years associated of the euclidean distance values 
    for i in np.arange(iyrpiC,fyrpiC,1):
        tmp= np.linalg.norm(monthlyValuesX.sel(time=str(year)).values-monthlyValuesXpic.sel(time=str(i)).values)
        distList.append(tmp)
        orderList.append(i)
    f=list(zip(distList,orderList))
    f_order=sorted(f, key=lambda x: x[0]) ## x[0] is distList elements: Ordering monthly fields by ascending order of their Euclidean distance to the target field
    #print('Taking closest Na=',Na,'analogs of the list Euc.Dist (1st element removed) \n')
    Na=Na
    f_ca=f_order[1:Na+1] #Remove the first field as it is itself with dist=0.0, always in the first position (to improve by removing also the five years aroundm maybe)
    ##Taking only the year analogs
    years_analogs_counter= [x[1] for x in f_ca] ##x[1] taking the 2nd element of the list (counter), the year analogs
    #type(years_analogs)
    #print('Year analogs are selected')
    #print(np.asarray(years_analogs_counter),'\n \n') ##Just to see which real years are the analogs
    #IV.- Construction column vector of selected Na analogs for X and Y with Na x 1 dimension
    ##Constructing column vector of slected Na analogs for X and Y
    Xc_colVector_Na=[]
    Yc_colVector_Na=[]
    for i in years_analogs_counter:
        tmp1=monthlyValuesXpic.sel(time=str(i))[0,:,:]#.values
        tmp2=monthlyValuesYpic.sel(time=str(i))[0,:,:] #extra dimension when selecting xarray
        Xc_colVector_Na.append(tmp1)
        Yc_colVector_Na.append(tmp2)
    #print(len(Xc_colVector_Na),len(Yc_colVector_Na),'\n \n')
    #V.-Resample randomly Ns analogs from Na (col vector Ns x 1 ). Do this process Nr times
    Ns=Ns;Nr=Nr
    #print('Ns:',Ns,'Nr:',Nr,'\n \n')
    Xc_colVector_Ns_Nr=[]
    Yc_colVector_Ns_Nr=[]
    for i in range(Nr):
        Xc_colVector_Ns,Yc_colVector_Ns=zip(*random.sample(list(zip(Xc_colVector_Na, Yc_colVector_Na)), Ns))
        Xc_array=np.asarray(Xc_colVector_Ns)
        Xc=Xc_array.reshape(Ns,nsxpic)
        Yc_array=np.asarray(Yc_colVector_Ns)
        Yc=Yc_array.reshape(Ns,nsypic)
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
   #Xca_avg.shape
    Yca_avg=np.asarray(Yc_colVector_Ns_Nr_mean)
    Yca_avg.shape
   # print('VIII.- Finally, convert them to a nlatxnlon grid. \n \n Return: Xh, Yh (nlatxnlon) january 1979 analogs')
    Xh=Xca_avg.reshape(nlatxpic,nlonxpic)
    Yh=Yca_avg.reshape(nlatypic,nlonypic)
    return Xh,Yh


####
##TO USE in case the time period of piControl is so large that it cannot be read by xarray
####
def climAnalogs_piC_numpy_2fields(monthlyValuesXpic,monthlyValuesX,monthlyValuesYpic,year,Na,Ns,Nr,iyrAmon,fyrAmon,iyrpiC,fyrpiC):
 #  print(" I.- Preparing the data for computation : matrix form time x m (m=nlatxnlon) \n ")
    dataxpic=monthlyValuesXpic
    dxpic=monthlyValuesXpic.shape;nlonxpic=dxpic[2];nlatxpic=dxpic[1];ntxpic=dxpic[0];nsxpic=nlonxpic*nlatxpic
    Xpic=dataxpic.reshape(ntxpic,nsxpic)
    
    datax=monthlyValuesX
    dx=monthlyValuesX.shape;nlonx=dx[2];nlatx=dx[1];ntx=dx[0];nsx=nlonx*nlatx
    X=datax.reshape(ntx,nsx)
    
    dataypic=monthlyValuesYpic
    dypic=monthlyValuesYpic.shape;nlonypic=dypic[2];nlatypic=dypic[1];ntypic=dypic[0];nsypic=nlonypic*nlatypic
    Ypic=dataypic.reshape(ntypic,nsypic)
   # print("II.- Target field \n")
    countYear=year-iyrAmon
    #print('Year',year)
    #print('CountYear',countYear)
    Xt=X[countYear,0:nsx]
    distList=[] ## Contains the values of the euclidean distance (r=2) between jan 1979 and all jan (1980-1998)
    orderList=[]
  #  print("III.- Computing the Euclidean distance between target SLP field",season,"/",iyr,"among the monthly fields of all nyr=",nyr,'\n ... \n')
    nyrpiC=fyrpiC-iyrpiC+1
    for i in range(nyrpiC):
        tmp=np.linalg.norm(X[countYear,0:nsx]-Xpic[i,0:nsxpic])
        distList.append(tmp)
        orderList.append(i)
    f=list(zip(distList,orderList))
    f_order=sorted(f, key=lambda x: x[0]) ## x[0] is distList elements: Ordering monthly fields by ascending order of their Euclidean distance to the target field
    Na=Na
    f_ca=f_order[1:Na+1] #Remove the first field as it is itself with dist=0.0, always in the first position (to improve by removing also the five years aroundm maybe)
    years_analogs_counter= [x[1] for x in f_ca] ##x[1] taking the 2nd element of the list (counter), the year analogs
 #   print('IV.- Construction column vector of selected Na analogs for X and Y with Na x 1 dimension \n \n')
    ##Constructing column vector of slected Na analogs for X and Y
    Xc_colVector_Na=[]
    Yc_colVector_Na=[]
    for i in years_analogs_counter:
        tmp0=Xpic[int(i),0:nsx]
        tmp1=Ypic[int(i),0:nsypic]
        Xc_colVector_Na.append(tmp0)
        Yc_colVector_Na.append(tmp1)
 #   print('V.-  Resample randomly Ns analogs from Na (col vector Ns x 1 ). Do this process Nr times \n \n')
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
    ##Average of "best" linear combinations
    Xca_avg=np.asarray(Xc_colVector_Ns_Nr_mean)
    Yca_avg=np.asarray(Yc_colVector_Ns_Nr_mean)
   # print('VIII.- Finally, convert them to a nlatxnlon grid. \n \n Return: Xh, Yh (nlatxnlon) january 1979 analogs')
    Xh=Xca_avg.reshape(nlatx,nlonx)
    Yh=Yca_avg.reshape(nlatypic,nlonypic)
    return Xh,Yh


#######
### Adapted version to use ESA CCI snow cover product with full nans on 10-11-12/1994 and 01/1995
#######
##Using SNC-------

def climAnalogs_2fields_ref(monthlyValuesX,monthlyValuesY,monthlyValuesXAn,domainXAn,year,Na,Ns,Nr,iyrAn,fyrAn,interval,uncertainty=False):
    '''
    Computed dynamical constructed analogs
    Parameters:
    -monthlyValuesX: dynamical field as target (SLP, geopotential...)
    -monthlyValuesY: associated field to be dynamically constructed ( temperature,snow cover, precipitation...)
    -monthlyValuesXref: dynamical field to pick/select the analogs (SLP piControl, SLP detrended, SLP original)
    -year: target year
    -Na,Ns,Nr: number of analogs, subsample, iterations for random subsample
    -iyrAn,fyrAn: period of choosing analogs, it can be different from field period anomalies( iyrAnom,fyrAnom)
    -interval: Interval of years around the target mon/year that are not considered
    -uncertainty: in case you want to save the random subsample to make a Montecarlo test for statistical significante (see O'Reilly et al. 2017 for details)
    
    >Application:
    Following Deser et al. 2016 methodology for CESM-LE large ensemble and NOAAv2c reanalysis:
    ##Dynamical adjustment
    Xh_year=[]
    Yh_year=[]
    for year in np.arange(iyrAnom,fyrAnom+1,1):
        Xh,Yh=climAnalogs_2fields_ref(monthlyValuesX=anomsX,monthlyValuesY=anomsYpiC,monthlyValuesXAn=anomsXpiC,
                                      domainXAn=domainXAn,year=year,Na=Na,Ns=Ns,Nr=Nr,
                                      iyrAn=iyrpiC,fyrAn=fyrpiC,interval=1,uncertainty=False)

        Xh_year.append(Xh)
        Yh_year.append(Yh)
    '''
    #print(Na,Ns,Nr,iyr,fyr)
    #print(" I.- Preparing the data for computation : matrix form time x m (m=nlatxnlon) \n ")
    datax=monthlyValuesX.values; dx=monthlyValuesX.shape;nlonx=dx[2];nlatx=dx[1];ntx=dx[0];nsx=nlonx*nlatx
    X=datax.reshape(ntx,nsx)
    
    datay=monthlyValuesY.values;dy=monthlyValuesY.shape;nlony=dy[2];nlaty=dy[1];nty=dy[0];nsy=nlony*nlaty
    Y=datay.reshape(nty,nsy)
    
    dataxAn=monthlyValuesXAn.values; dxAn=monthlyValuesXAn.shape;nlonxAn=dxAn[2];nlatxAn=dxAn[1];ntxAn=dxAn[0];nsxAn=nlonxAn*nlatxAn
    XAn=dataxAn.reshape(ntxAn,nsxAn)
  #  print(X.shape) 
    ##II.-Target field (Xt SLP jan 1979)
    #print(year)
    Xt_xarray=monthlyValuesX.sel(time=str(year))[0,:,:]
    Xt=Xt_xarray.values.reshape(nsx)
    
    Yt_xarray=monthlyValuesY.sel(time=str(year))[0,:,:]
    Yt=Yt_xarray.values.reshape(nsy)
    #
    ##III.- Euclidean distance
    distList=[] ## Contains the values of the euclidean distance (r=2) between jan 1979 and all jan (1980-1998)
    orderList=[]
    range_years=np.arange(iyrAn,fyrAn+1,1).tolist()
    range_years.remove(1994)
    range_years.remove(1995)
    for i in range_years:
        tmp= np.linalg.norm(dom.field_dom(monthlyValuesX.sel(time=str(year)),domainXAn).values-dom.field_dom(monthlyValuesXAn.sel(time=str(i)),domainXAn).values)
        #print(len(dom.field_dom(monthlyValuesX.sel(time=str(year)),domainXAn).values))
        distList.append(tmp)
        orderList.append(i)
    f=list(zip(distList,orderList))
    f_order=sorted(f, key=lambda x: x[0]) ## x[0] is distList elements: Ordering monthly fields by ascending order of their Euclidean distance to the target field
    ##Taking only the year analogs
    allyears= [x[1] for x in f_order]
    interval=interval
    forbiddenyears=(np.arange(year-interval,year+interval+1,1)).tolist()
    #print('Year',year,'Filter',forbiddenyears) # In this way we remove already the year itself
    allyears_filter=[x for x in allyears if x not in forbiddenyears]
    #print(allyears_filter)
    years_analogs_counter=[]
    years_analogs_counter= allyears_filter[:Na]
    #print(years_analogs_counter)
    #
    #IV.- Construction column vector of selected Na analogs for X and Y with Na x 1 dimension \n \n')
    ##Constructing column vector of slected Na analogs for X and Y
    Xc_colVector_Na=[]
    Yc_colVector_Na=[]
    for i in years_analogs_counter:
        tmp1=monthlyValuesXAn.sel(time=str(i))[0,:,:]#.values
        tmp2=monthlyValuesY.sel(time=str(i))[0,:,:] #extra dimension when selecting xarray
        Xc_colVector_Na.append(tmp1)
        Yc_colVector_Na.append(tmp2)
    #print(len(Xc_colVector_Na),len(Yc_colVector_Na),'\n \n')
    #V.-  Resample randomly Ns analogs from Na (col vector Ns x 1 ). Do this process Nr times \n \n')
    #Pass the list to the first argument and the number of elements you want to get to the second argument. A list is returned
    Ns=Ns;Nr=Nr
    #print('Ns:',Ns,'Nr:',Nr,'\n \n')
    Xc_colVector_Ns_Nr=[]
    Yc_colVector_Ns_Nr=[]
    for i in range(Nr):
        Xc_colVector_Ns,Yc_colVector_Ns=zip(*random.sample(list(zip(Xc_colVector_Na, Yc_colVector_Na)), Ns))
        Xc_array=np.asarray(Xc_colVector_Ns)
        Xc=Xc_array.reshape(Ns,nsx)
        Yc_array=np.asarray(Yc_colVector_Ns)
        Yc=Yc_array.reshape(Ns,nsy)
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

    #VI.- Finally, convert them to a nlatxnlon grid. \n \n Return: Xh, Yh (nlatxnlon) january 1979 analogs
    Xh=Xca_avg.reshape(nlatx,nlonx)
    Yh=Yca_avg.reshape(nlaty,nlony)
    Xh_Nr=Xca_Nr.reshape(Nr,nlatx,nlonx)
    Yh_Nr=Yca_Nr.reshape(Nr,nlaty,nlony)
    if uncertainty==True:
        return Xh,Yh,Xh_Nr,Yh_Nr
    else:
        return Xh,Yh