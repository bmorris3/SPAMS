# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 14:46:31 2014

@author: bmmorris
"""

import photometry
import trackSmooth
import pyfits
import numpy as np
from matplotlib import pyplot as plt
from glob import glob
from astropy.time import Time
from scipy import optimize
import datetime
import ephem

regionsfile  = open('stars.reg', 'r').read().splitlines() # MAKE A REGIONS FILE WITH YOUR STARS
images = sorted(glob('/astro/store/scratch/tmp/bmmorris/APO/UT150117/sdss08.*.fits'))[:-20] # Put in the expression to pick up the images that you want

starpositions = []
counter = 0
for line in regionsfile:
    if line.startswith('circle'):
        y_center, x_center = map(float,line.split('(')[1].split(')')[0].split(',')[:2])

        starpositions.append([y_center, x_center])
        #starnames.append('comp'+str(counter))
        counter += 1

makemasterdark = True
makemasterflat = True
if makemasterdark:
    #darkpaths = glob('/local/tmp/AGILE/D27apr2012M-4x4-10*.fits')
    darkpaths = ['/astro/store/scratch/tmp/bmmorris/APO/UT150117/dark-sdss08.{0:04d}.fits'.format(i) 
                 for i in range(1841, 1941)]
    testhdr = pyfits.getheader(darkpaths[0])
    fullx, fully = testhdr['FULLX'], testhdr['FULLY']
    alldarks = np.zeros((fullx, fully, len(darkpaths)))
    for i, darkpath in enumerate(darkpaths):
        alldarks[:,:,i] = pyfits.getdata(darkpath)[:fully, :fullx]
    masterdark = np.median(alldarks,axis=2)
    pyfits.writeto('masterdark.fits', masterdark, clobber=True)
else:
    masterdark = pyfits.getdata('masterdark.fits')    

if makemasterflat: 
    
    flatdarkpaths = glob('/astro/store/scratch/tmp/bmmorris/APO/UT150117/flatdarks.????.fits')
    testhdr = pyfits.getheader(flatdarkpaths[0])
    fullx, fully = testhdr['FULLX'], testhdr['FULLY']
    allflatdarks = np.zeros((fullx, fully, len(flatdarkpaths)))
    for i, darkpath in enumerate(flatdarkpaths):
        allflatdarks[:,:,i] = pyfits.getdata(darkpath)[:fully, :fullx]
    masterflatdark = np.median(allflatdarks,axis=2)
    
    flatpaths = glob('/astro/store/scratch/tmp/bmmorris/APO/UT150117/flats.????.fits')
    testhdr = pyfits.getheader(flatpaths[0])
    fullx, fully = testhdr['FULLX'], testhdr['FULLY']
    allflats = np.zeros((fullx, fully, len(flatpaths)))
    for i, flatpath in enumerate(flatpaths):
        allflats[:,:,i] = pyfits.getdata(flatpath)[:fully, :fullx] - masterflatdark
    masterflat = np.median(allflats,axis=2)
    masterflat /= np.median(masterflat)
    pyfits.writeto('masterflat.fits', masterflat, clobber=True)
else:
    masterflat = pyfits.getdata('masterflat.fits')    

apertureradii = [6.5]#np.arange(2, 15, 0.5) #[5.5]
stds = np.zeros_like(apertureradii,dtype=float)

times = np.zeros(len(images))
fluxes = np.zeros((len(images),len(starpositions), len(apertureradii)))
errors = np.zeros((len(images),len(starpositions), len(apertureradii)))
xcentroids = np.zeros((len(images),len(starpositions)))
ycentroids = np.zeros((len(images),len(starpositions)))
subtractedfluxes = np.zeros((len(images),len(starpositions), len(apertureradii)))
PSFsigmas = np.zeros((len(images),len(starpositions), len(apertureradii)))
airmass = np.zeros(len(images))
expdur = np.zeros(len(images))

means = np.zeros((len(images), len(apertureradii)))
medians = np.zeros((len(images), len(apertureradii)))

trackplots = False

for i in range(len(images)):
    print "%d of %d" % (i+1, len(images))
    imageheader = pyfits.getheader(images[i])
    imagedata = (pyfits.getdata(images[i])[:imageheader['FULLY'],:imageheader['FULLX']] - masterdark)/masterflat
    timeISO = imageheader['DATE-OBS']#imageheader['DATE-OBS'] 
    timeJD = Time(timeISO, format='isot', scale='utc').jd
    times[i] = timeJD
    zoomfactor = 25
    smoothconst = 2.5

    if trackplots and i == 0:
        fig = plt.figure(num=None, figsize=(18, 3), facecolor='w', edgecolor='k')
        fig.subplots_adjust(wspace=0.5)
        subplotsDimensions = 140
        photSubplotsOffset = 3
        statusSubplotOffset = 6
        plottingThings = [fig,subplotsDimensions,photSubplotsOffset]

    for j in range(len(starpositions)):
        deltax, deltay = starpositions[j]
        if i == 0:
            init_x = starpositions[j][0]
            init_y = starpositions[j][1]
        else:
            init_x = ycentroids[i-1][j]#starpositions[j][i-1]
            init_y = xcentroids[i-1][j]#starpositions[j][i-1]
        
        if trackplots:
            fig = plt.figure(num=None, figsize=(18, 3), facecolor='w', edgecolor='k')
            fig.subplots_adjust(wspace=0.5)
            subplotsDimensions = 140
            photSubplotsOffset = 3
            statusSubplotOffset = 6
            plottingThings = [fig,subplotsDimensions,photSubplotsOffset]
        else: 
            plottingThings = False


        [xCenter,yCenter,averageRadius, _] = trackSmooth.trackSmooth(imagedata,init_y, init_x, \
            smoothconst, plottingThings, preCropped=False, zoom=zoomfactor, plots=trackplots)
            
        for ap, apertureradius in enumerate(apertureradii):
            flux, error, photFlags, _ = photometry.phot(imagedata, xCenter, yCenter,
                                                            apertureradius,
                                                            plottingThings,
                                                            ccdGain=1.0,
                                                            plots=trackplots,\
                    annulusOuterRadiusFactor=2.0, annulusInnerRadiusFactor=1.5,\
                    sigmaclipping=True, returnsubtractedflux=True)

            fluxes[i, j, ap] = flux
            errors[i, j, ap] = error
        xcentroids[i, j] = xCenter
        ycentroids[i, j] = yCenter
    if trackplots:
        plt.show()

#testind = 0
#compind = 1
#lc = fluxes[:, testind, :]/fluxes[:, compind, :]
#
#lc_err = lc*np.sqrt( (errors[:, testind, :]/fluxes[:, testind, :])**2 + 
#         (errors[:, compind, :]/fluxes[:, compind, :])**2 )
#
#for ap in range(len(apertureradii)):
#    lcmedian = np.median(lc[:,ap])
#    lc[:,ap] /= lcmedian
#    lc_err[:,ap] /= lcmedian
#
#    m = np.mean(lc[:,ap])
#    s = np.std(lc[:,ap])
#    n = 4
#    nonoutliers = (lc[:,ap] < m + n*s)*(lc[:,ap] > m - n*s)
#    
#    stds[ap] = np.std(lc[nonoutliers, ap])
##    # remove airmass
##    def fitfunc(p):
##        return p[0]*airmass + p[1]
##    
##    def errfunc(p):
##        return (fitfunc(p) - lc)/lc_err
##    
##    initp = [1, 1]
##    bestp = optimize.leastsq(errfunc, initp)[0]
##    correctedlc = lc - fitfunc(bestp) + 1
#
#
##plt.plot(apertureradii, stds)
##np.savetxt('GD165lightcurve.txt',np.vstack([times, correctedlc, lc_err]))
#plt.errorbar(times, lc, yerr=lc_err, fmt='.')
#plt.show()

## Taken from git/research/koi351/variableaperture.ipynb
stds = np.zeros_like(apertureradii)
for j, ap in enumerate(apertureradii):
    target = fluxes[:, 0, j]
    compStars = fluxes[:,1:, j]

    numCompStars = np.shape(fluxes[:, 1:, j])[1]
    initP = np.zeros([numCompStars]) + 1./numCompStars
    #initP = np.ones((numCompStars, 1))/numCompStars
    ########################################################

    def errfunc(p,target):
        if all(p >=0.0): 
            #return np.dot(p,target.T) - compStars ## Find only positive coefficients
            return np.dot(p,compStars[:,:].T) - target ## Find only positive coefficients

    #return np.dot(p,compStarsOOT.T) - target
    bestFitP = optimize.leastsq(errfunc, initP[:] ,args=(target.astype(np.float64)),\
                                maxfev=10000000,epsfcn=np.finfo(np.float32).eps)[0]
    print '\nDefault weight:',1./numCompStars
    print 'Best fit regression coefficients:',bestFitP

    #self.comparisonStarWeights = np.vstack([compStarKeys,bestFitP])
    meanComparisonStar = np.dot(bestFitP, compStars.T)

    meanCompError = np.zeros_like(meanComparisonStar)
    for i in range(1,len(fluxes[0, :, j])):
        meanCompError += (bestFitP[i-1]*fluxes[:, i, j]/np.sum(bestFitP[i-1]*fluxes[:, i, j]))**2 * \
                         (errors[:, i, j]/fluxes[:, i, j])**2
    meanCompError = meanComparisonStar*np.sqrt(meanCompError)

    lc = fluxes[:, 0, j]/meanComparisonStar
    lc /= np.median(lc)
    lc_err = lc*np.sqrt((errors[:, 0, j]/fluxes[:, 0, j])**2 + (meanCompError/meanComparisonStar)**2)
    #reducedchi2 = np.sum(((lc-np.ones_like(lc))/lc_err)**2)/(len(lc) - 1)

    ####################################
    # Compute airmasses

    # SDSS 082346.00+201557.13
    starRA = ephem.hours('08:23:46.00')
    starDec = ephem.degrees('+20:15:57.13')
    observatory_minHorizon = '25:00:00'
    apo = ephem.Observer()
    apo.lat =  '32.0:46.0:49.0' 
    apo.long = '-105.0:49.0:13.0' 
    apo.elevation = 2788.0 # m
    apo.temp = 10.0  ## Celsius
    apo.horizon = observatory_minHorizon

    airmass = np.zeros_like(times)
    for i, t in enumerate(times):
        t_dt = Time(t, format='jd', scale='utc').datetime
        apo.date = t_dt.isoformat(' ')
        star = ephem.FixedBody()
        star._ra = ephem.hours(starRA)
        star._dec = ephem.degrees(starDec)
        star.compute(apo)
        z = np.pi/2 - float(star.alt)
        airmass[i] = 1./np.cos(z)

    ####################################

    # Fit airmass to the light curve
    def amtrend(p, t, airmass=airmass):
        '''p = [OOT flux, airmass coeff]'''
        return p[0]*(1 + p[1]*(airmass-1))

    def errfunc(p, t, y, airmass=airmass):
        '''
        p = [OOT flux, airmass coeff]
        y = light curve
        '''
        return amtrend(p, t) - y

    bestp = optimize.leastsq(errfunc, [1.0, 0.0], args=(times, lc))[0]
    correctedlc = lc/amtrend(bestp, times)
    correctedlc_err = lc_err/amtrend(bestp, times)



    # Make lc out of just comp stars
    comp0 = fluxes[:,1, j]
    comp1 = fluxes[:,2, j]
    bestp = optimize.leastsq(errfunc, [1.0, 0.0], args=(times, comp0/comp1))[0]
    correctedcomplc = comp0/comp1/amtrend(bestp, times)
    if False:
        fig, ax = plt.subplots(1, 2, sharex=True, sharey=True)
        ax[0].plot(times, comp0/comp1, '.')
        ax[0].plot(times, amtrend(bestp, times),'r')
        ax[1].plot(times, correctedcomplc,'.')
        plt.show()
    #correctedcomplc_err = complc_err/amtrend(bestp, times)

    stds[j] = np.std(correctedcomplc)
    

fig, ax = plt.subplots(1)
ax.errorbar(times, lc, yerr=lc_err, fmt='.', color='k')
ax.plot(times, amtrend(bestp, times), color='r')
plt.show()

import scipy
import scipy.fftpack

fig, ax = plt.subplots(2, 2, figsize=(14,14))
mininttime = int(np.min(times))
legendprops = {'numpoints':1, 'loc':'lower center'}

np.savetxt("lcexperiments",correctedlc)
np.savetxt("timesexperiments",times-mininttime)

ax[0, 0].plot(times - mininttime, target, '.', label='SDSS 08')
ax[0, 0].plot(times - mininttime, compStars[:,0], '.', label='Comp 1')
ax[0, 0].plot(times - mininttime, compStars[:,1], '.', label='Comp 2')
ax[0, 0].set_xlabel('JD - {0:d}'.format(mininttime))
ax[0, 0].set_ylabel('Counts')
ax[0, 0].set_title('Raw fluxes')
ax[0, 0].legend(**legendprops)

ax[0, 1].plot(times - mininttime, target, '.', label='SDSS 08')
ax[0, 1].plot(times - mininttime, meanComparisonStar, '.', label='Mean Comp')
ax[0, 1].set_xlabel('JD - {0:d}'.format(mininttime))
ax[0, 1].set_ylabel('Counts')
ax[0, 1].legend(**legendprops)
ax[0, 1].set_title('Mean comparison star')

ax[1, 0].errorbar(times - mininttime, correctedlc, yerr=correctedlc_err, 
                  fmt='.', color='k', ecolor='gray')
ax[1, 0].set_title('AM-corrected light curve')
ax[1, 0].set_xlabel('JD - {0:d}'.format(mininttime))
ax[1, 0].set_ylabel('Normalized flux')

FFT = np.abs(scipy.fft(correctedlc))
dt = times[1]-times[0] # days
dt *= 24*60
freqs = scipy.fftpack.fftfreq(correctedlc.size, dt)
#ax[1, 1].plot(freqs, 20*np.log10(FFT))
#ax[1, 1].set_title('FFT')
fig.subplots_adjust(hspace=0.3)
fig.savefig('quicklc.pdf', bbox_inches='tight')

np.savetxt('SDSS082346.00+201557.13.txt',np.vstack([times, correctedlc,correctedlc_err]).T, header='JD Flux Error')

from scipy.ndimage import filters
fig, ax = plt.subplots(1, figsize=(14,10))
ax.errorbar(times - mininttime, correctedlc, yerr=correctedlc_err, 
                  fmt='.', color='k', ecolor='gray')
ax.plot(times - mininttime, filters.gaussian_filter1d(correctedlc, 3), color='r', lw=2)
ax.set_title('Airmass-corrected light curve')
ax.set_xlabel('JD - {0:d}'.format(mininttime))
ax.set_ylabel('Normalized flux')
fig.savefig('lightcurve.pdf', bbox_inches='tight')
plt.show()  

