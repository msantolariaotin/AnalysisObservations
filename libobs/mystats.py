##Basic libraries
import pandas as pd
from netCDF4 import Dataset
import numpy as np
import xarray as xr
import glob
import math


##Statistics
import pymannkendall as mk
from eofs.xarray import Eof
import statsmodels.api as sm
import scipy
from scipy import stats
from scipy.stats import t
from scipy.stats import f

#####
##
## Statistical significance trend/cor
##
####
def t_test_rs(alpha,N):
    confidence=1-alpha/2
    df=N-1
    t=stats.t.ppf(confidence,df)
    rs=math.sqrt((t*t)/(t*t + df))
    return rs
######---------------------------------------------
##
## Two-Sample t-Test for Equal Means
##
#####----------------------------------------------
#https://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm#:~:text=The%20two%2Dsample%20t%2Dtest,several%20variations%20on%20this%20test.
#The two-sample t-test for unpaired data is defined as: 
#H0: mu1=mu2
#H1: mu1!=mu2

#Test Statistic:  T = (Y1_mean - Y2_mean) / np.sqrt((s1**2/N1)+(s2**2+N2))
##where N1 and N2 are the sample sizes, Y1_mean and Y2_mean are the sample means, and s1**2 and s2**2 are the sample variances.
##Significance level: alphe
##Critical region: reject the null hypothesis that the two means are equal if  
## |T| > t_(1-alpha/2, df)
##where t_(1-alpha/2, df) is the critical value of the t distribution with df degrees of freedom
# if equal variances are assumed df=N1+N2-2

def t_test_difmeans(Y1,Y2,s1,s2,N1,N2,alpha):
    '''
    Ilustrative function, to apply at your own.
    '''
    Y1_mean=Y1.mean();Y2_mean=Y2.mean()
    s1=Y1.std();s2=Y2.std()
    s1_2=s1**2; s2_2=s2**2
    df1=N1;df2=N2
    frac1= s1_2/df1;frac2=s2_2/df2
    t_obs=(Y1_mean-Y2_mean)/np.sqrt(frac1+frac2)
    confidence=1-alpha/2
    df=df1+df2-2
    t_teo=stats.t.ppf(confidence,df)
    return t_teo,t_obs

#####----------------------------------------------
##                                        
## F Fisher test for standard deviation  
##
#####----------------------------------------------
##https://www.statology.org/f-test-python/
#An F-test is used to test whether two population variances are equal.
#The null and alternative hypotheses for the test are as follows:
#H0: s1**2 = s2**2 (the population variances are equal)
#H1: s1**2 ≠ s2**2 (the population variances are not equal)
#https://mgimond.github.io/Stats-in-R/F_test.html
#The null hypothesis is rejected i F is either too large ( s1**2/ s2**2 >>1 ) or too small (s1**2/ s2**2 <<1 ) based on the desired alpha.
def f_test(s1,s2,N1,N2,alpha):
    '''
Ilustrative function, to apply at your own.
    '''
    F_obs=s1**2/ s2**2 
    df1=N1-1 #numerator
    df2=N2-1  #denominator
    qmin=scipy.stats.f.ppf(q=alpha/2,dfn=df1,dfd=df2) 
    qmax=scipy.stats.f.ppf(q=1-alpha/2,dfn=df1,dfd=df2)
    return Fobs,qmin,qmax
'''
#Illustrative plot F-test:
x = np.linspace(0.2,3, 100) ##Fs
qmin=scipy.stats.f.ppf(q=0.025,dfn=26,dfd=26)
qmax=scipy.stats.f.ppf(q=1-0.025,dfn=26,dfd=26)
Ftheo=f.pdf(x, 26,26)
Fobs=1.2
plt.plot(x, Ftheo,'r-', lw=5, alpha=0.6, label='f pdf')
plt.axvline(x=qmin,color='blue')
plt.axvline(x=qmax,color='blue')
plt.axvline(x=Fobs,color='green')
plt.show()
'''

#####----------------------------------------------
##
##Autocorrelation series
##
####----------------------------------------------
'''From Javier Garcia-Serrano R codei ( CFU_Rfunc_test.txt):
#CFU_neff<-function(vector, fig=FALSE){
#
# Compute the effective degrees of freedom of a time series
#
# Description:
#
#       Number of effective degrees of freedom of "vector".
#       Based on: Zieba, A. (2010): Effective number of observations and unbiased estimators
#       of variance for autocorrelated data - an overview. Metrol. Meas. Syst. Vol XVII, No. 1,
#       pp. 3-16; index 330930, ISSN 0860-8229. www.metrology.pg.gda.pl
#
# Authors:
#
#      Created by Caio A. S. Coelho <caio.coelho@cptec.inpe.br>
#      Implemented in CFU_Rfunc by Javier García-Serrano <jgarcia@ic3.cat> (August 2011)


dof=length(vector)

a=acf(vector,lag.max=dof-1,plot=fig)$acf[2:dof,1,1]

s=0
for (k in 1:(dof-1)){
#    s=s+(((dof-k)/dof)*a[k])
    s=s+(((dof-k)/dof)*(a[k]**2))
}

#neff=dof/(1+(2*s))
neff=(dof/(1+(2*s)))-1

#if (neff>dof){neff=dof}

#outputs

#list(neff=neff)
neff

}
#https://scicoding.com/4-ways-of-calculating-autocorrelation-in-python/
'''
##Adapted to Python:

def autocor_neff(data):
    '''
Ilustrative function, to apply at your own.
data = [3, 16, 156, 47, 246, 176, 233, 140, 130, 
        101, 166, 201, 200, 116, 118, 247, 
        209, 52, 153, 232, 128, 27, 192, 168, 208, 
        187, 228, 86, 30, 151, 18, 254, 
        76, 112, 67, 244, 179, 150, 89, 49, 83, 147, 90, 
        33, 6, 158, 80, 35, 186, 127]
'''
    # Mean
    mean = np.mean(data)
    # Variance
    var = np.var(data)
    # Normalized data
    ndata = data - mean
    acorr = np.correlate(ndata, ndata, 'full')[len(ndata)-1:]
    acorr = acorr / var / len(ndata)
    
    dof=len(data)
    s=0
    for k in np.arange(0,dof4,1):
        s=s+(((dof-k)/dof)*(acorr[k]**2))

    neff=(dof/(1+(2*s)))-1
    return neff





