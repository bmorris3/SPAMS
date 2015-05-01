# -*- coding: utf-8 -*-
"""
Created on Fri Feb 13 15:33:02 2015

@author: brantd
"""
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import leastsq
fluxes = np.loadtxt("lcexperiments")
times2 = np.loadtxt("timesexperiments")
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
def lc(RpRs,aRs,P,i,t0,t):
    b=aRs*np.cos(np.radians(i))
    duration=P/(np.pi*aRs)*np.sqrt(1-b**2)
    depth=(RpRs)**2
    return transit(t,depth,duration,t0,P)

###
lcguess=(0,1,2,3,4,5)
betterlc=leastsq(lc(1,1),lcguess,args=(times2))[0]           
###

plt.plot(times2,model)
plt.plot(times2,fluxes)
plt.show()

