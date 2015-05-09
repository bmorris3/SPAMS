# -*- coding: utf-8 -*-
"""
Created on Fri Feb 13 15:33:02 2015

@author: brantd
"""
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from scipy.optimize import leastsq
d = np.loadtxt("../data/SDSS115219.99+024814.4.txt") #binary
e= np.loadtxt("../data/SDSS082346.00+201557.13.txt")
f=np.loadtxt("../data/SDSS160401.31+083109.01.txt")
def gridsearch(d):
    fluxes=d[:,1]
    times2=d[:,0]
    #fluxes = np.loadtxt("../data/SDSS082346.00+201557.13.txt", usecols=(1))
    #times2 = np.loadtxt("../data/SDSS082346.00+201557.13.txt", usecols=(0))
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
    #errors=np.zeros_like(times2)+0.01
    errors=d[:,2]
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
    P=.5
    inclination=90
    RpRs_steps=np.linspace(0,1,10)
    aRs_steps=np.linspace(1,100,10)
    chisquared=np.zeros((len(RpRs_steps),len(aRs_steps),len(times2)))
    for i,RpRs_i in enumerate(RpRs_steps):
        for j,aRs_i in enumerate(aRs_steps):
            for k,t_i in enumerate(times2):
                model=lc(times2,RpRs_i,aRs_i,P,inclination,t_i)
                chi2=np.sum((fluxes-model)**2)
                chisquared[i,j,k]=chi2
    
    
    def makeplot(i):
        fig = plt.figure()
        sheet = chisquared[:, :, i]
        ax = fig.add_subplot(111, projection='3d')
        XX, YY = np.meshgrid(RpRs_steps, aRs_steps)
        #ax.plot_surface(RpRs_steps, aRs_steps,sheet, rstride=1, cstride=1)
        ax.plot_surface(XX, YY, np.log(sheet), rstride=1, cstride=1)
        
        plt.show()
    bestRpRsarray=[] #creating an array of RpRs values based off Chi2
    for k,t_i in enumerate(times2):
        chisquared_sheet= chisquared[:,:,k]
        minindex= chisquared_sheet==np.min(chisquared_sheet)
        bestRpRs= RpRs_steps[np.sum(minindex, axis=1).astype(bool)][0]
        bestRpRsarray.append(bestRpRs)

    return times2, fluxes, bestRpRsarray
    
dtimes,dfluxes,dbestRpRsarray=gridsearch(d)
etimes,efluxes,ebestRpRsarray=gridsearch(e)
ftimes,ffluxes,fbestRpRsarray=gridsearch(f)

#plt.plot(times2,bestRpRsarray)
#plt.plot(times2,fluxes,'.')
f, (ax1, ax2, ax3)= plt.subplots(3,1)
ax1.plot(dtimes, dbestRpRsarray,'r')
ax1.plot(dtimes, dfluxes,'b.')
ax1.set_title('SDSS115219.99+024814.4')
#ax1.plt.grid(True)

ax2.plot(etimes, ebestRpRsarray,'r')
ax2.plot(etimes, efluxes,'b.')
ax2.set_title('SDSS082346.00+201557.13')
#ax2.plt.grid(True)
ax2.set_ylabel('Normalized flux')
ax3.plot(ftimes, fbestRpRsarray,'r')
ax3.plot(ftimes, ffluxes,'b.')
ax3.set_title('SDSS160401.31+083109.01')
ax3.set_xlabel('Time (Julian Date)')
#ax3.plt.grid(True)
plt.show()
