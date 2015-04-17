# -*- coding: utf-8 -*-
"""
Created on Fri Feb 13 15:33:02 2015

@author: brantd
"""
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import leastsq
fluxes = np.loadtxt("../agileanalysis/lcexperiments")
times2 = np.loadtxt("../agileanalysis/timesexperiments")
times=np.linspace(0,1,1e4)
def transit(t, depth, duration,midtransit,period,F0=1):
    transitflux=np.ones_like(t) #allocates an array of 1's
    #Make a boolean array where true is in transit
    secondhalf= (.5*duration)>((t-midtransit)%period)
    firsthalf= (period-.5*duration)<((t-midtransit)%period)
    intransit=firsthalf+secondhalf
    transitflux[intransit]-=depth #if in transit subtract flux by depth
    return F0*transitflux

def chi2(data,model,error):
    return np.sum((data-model)**2/error**2)
errors=np.zeros_like(times2)+0.01
def redchi2(data,model,error,dof):
    return chi2(data,model,error)/(len(data)-dof)
model=transit(times2,0.07,.2,.5,1)
print redchi2(fluxes,model,errors,3)

#for t in np.arange(.3,.5,.02):
#    model=transit(times2,0.07,.2,t,1)
#    print t,redchi2(fluxes,model,errors,3)
def lc(t,RpRs,aRs,P,i,t0):
    b=aRs*np.cos(np.radians(i))
    duration=P/(np.pi*aRs)*np.sqrt(1-b**2)
    depth=(RpRs)**2
    return transit(t,depth,duration,t0,P)
P=3
inclination=90
RpRs_steps=np.linspace(0,1,10)
aRs_steps=np.linspace(1,10,10)
chisquared=np.zeros((len(RpRs_steps),len(aRs_steps),len(times2)))
for i,RpRs_i in enumerate(RpRs_steps):
    for j,aRs_i in enumerate(aRs_steps):
        for k,t_i in enumerate(times2):
            model=lc(times2,RpRs_i,aRs_i,P,inclination,t_i)
            chi2=np.sum((fluxes-model)**2)
            chisquared[i,j,k]=chi2