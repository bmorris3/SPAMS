# -*- coding: utf-8 -*-
"""
Created on Fri Feb 13 15:33:02 2015

@author: brantd
"""
import numpy as np
from matplotlib import pyplot as plt

fluxes = np.loadtxt("fluxes_20150213.txt")
#times = np.loadtxt("times_20150213.txt")
times=np.linspace(0,1,1e4)
def transit(t, depth, duration,midtransit,period):
    transitflux=np.ones_like(t)
    bob= (.5*duration)>((t-midtransit)%period)
    steve= (period-.5*duration)<((t-midtransit)%period)
    intransit=bob+steve
    print intransit
    print ((t-midtransit)%period)
    #intransit= (midtransit-.5*duration<t)*(t<midtransit+.5*duration)
    transitflux[intransit]-=depth
    
    #transdrop=depth*np.ones(duration)
   # for i in range(0, len(transdrop)):
       # transitflux[i+offset]-=transdrop[i]
    return transitflux
plt.plot(times,transit(times,0.3,.1,100000,.35))
plt.show()