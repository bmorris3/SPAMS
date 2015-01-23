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
from astropy import wcs
import astropy.io.fits as fits
import ephem

secondnightimages = glob('/local/tmp/ARCSAT/20141107/WD2218+706_sdss_g_20141108*.fits')
thirdnightimages = glob('/local/tmp/ARCSAT/20141108/WD2218+706_sdss_g_20141109*.fits')
fourthnightimages = glob('/local/tmp/ARCSAT/20141109/WD2218+706_sdss_g_20141110*.fits')
images = sorted(glob('/local/tmp/ARCSAT/20141106/WD2218+706_sdss_g_20141107_06*.fits') +\
         glob('/local/tmp/ARCSAT/20141106/WD2218+706_sdss_g_20141107_07*.fits') +\
         secondnightimages + thirdnightimages + fourthnightimages)# Put in the expression to pick up the images that you want
#images = images[:-1]

regionsfile  = open('stars_fk5.reg', 'r').read().splitlines() # MAKE A REGIONS FILE WITH YOUR STARS
## Initial positions
starpositions = []
counter = 0
for line in regionsfile:
    if line.startswith('circle'):
        y_center, x_center = line.split('(')[1].split(')')[0].split(',')[:2]
        starpositions.append([y_center, x_center])
        #starnames.append('comp'+str(counter))
        counter += 1

makemasterdark = False
makemasterflat = False
if makemasterdark:
    darkpaths = glob('/local/tmp/ARCSAT/20141106/Dark_B2_20141107_0*.fits')
    testhdr = pyfits.getheader(darkpaths[0])
    fullx, fully = 512, 512#testhdr['FULLX'], testhdr['FULLY']
    alldarks = np.zeros((fullx, fully, len(darkpaths)))
    for i, darkpath in enumerate(darkpaths):
        if pyfits.getheader(darkpath)['EXPTIME'] == 45:
            alldarks[:,:,i] = pyfits.getdata(darkpath)[:fully, :fullx]
    masterdark = np.median(alldarks,axis=2)
    pyfits.writeto('masterdark.fits', masterdark, clobber=True)

    testhdr = pyfits.getheader(darkpaths[0])
    fullx, fully = 512, 512#testhdr['FULLX'], testhdr['FULLY']
    alldarks = np.zeros((fullx, fully, len(darkpaths)))
    for i, darkpath in enumerate(darkpaths):
        if pyfits.getheader(darkpath)['EXPTIME'] == 1:
            alldarks[:,:,i] = pyfits.getdata(darkpath)[:fully, :fullx]
    flatdark = np.median(alldarks,axis=2)
    pyfits.writeto('flatdark.fits', flatdark, clobber=True)

else:
    masterdark = pyfits.getdata('masterdark.fits')    

if makemasterflat: 
    flatpaths = glob('/local/tmp/ARCSAT/20141106/domeflat_sdss_g_*.fits')
    testhdr = pyfits.getheader(flatpaths[0])
    flatdark = pyfits.getdata('flatdark.fits')
    fullx, fully = 512, 512#testhdr['FULLX'], testhdr['FULLY']
    allflats = np.zeros((fullx, fully, len(flatpaths)))
    for i, flatpath in enumerate(flatpaths):
        allflats[:,:,i] = pyfits.getdata(flatpath)[:fully, :fullx] - flatdark
    masterflat = np.median(allflats,axis=2)
    masterflat /= np.median(masterflat)
    pyfits.writeto('masterflat.fits', masterflat, clobber=True)
else:
    masterflat = pyfits.getdata('masterflat.fits')    

apertureradii = [8]#np.arange(5,15)#[8]#
stds = np.zeros_like(apertureradii,dtype=float)
for ap, apertureradius in enumerate(apertureradii):
    times = np.zeros(len(images))
    fluxes = np.zeros((len(images),len(starpositions)))
    errors = np.zeros((len(images),len(starpositions)))
    xcentroids = np.zeros((len(images),len(starpositions)))
    ycentroids = np.zeros((len(images),len(starpositions)))
    subtractedfluxes = np.zeros((len(images),len(starpositions)))
    PSFsigmas = np.zeros((len(images),len(starpositions)))
    airmass = np.zeros(len(images))
    expdur = np.zeros(len(images))
    
    means = np.zeros(len(images))
    medians = np.zeros(len(images))
    trackplots = False#True
    
    for i in range(len(images)):
        print "%d of %d" % (i+1, len(images))
        imageheader = pyfits.getheader(images[i])
        imagedata = (pyfits.getdata(images[i]) - masterdark)/masterflat
        timeISO = imageheader['DATE-OBS']#imageheader['DATE-OBS'] 
        means[i] = np.mean(imagedata)
        medians[i] = np.median(imagedata)
        timeJD = Time(timeISO, format='isot', scale='utc').jd
        times[i] = timeJD
        airmass[i] = imageheader['AIRMASS']
        hdulist = fits.open(images[i])
        w = wcs.WCS(hdulist[0].header, hdulist)
        zoomfactor = 15#25#15#25      ##zoom=15, smooth=3.5
        smoothconst = 3.5# 2.5
        def trackandphot(init_x,init_y, plottingThings):
            [xCenter,yCenter,averageRadius, _] = trackSmooth.trackSmooth(imagedata,init_y, init_x, \
                smoothconst, plottingThings, preCropped=False, zoom=zoomfactor, plots=trackplots)
            flux, error, photFlags, _ = photometry.phot(imagedata, xCenter, yCenter,
                                                            apertureradius,
                                                            plottingThings,
                                                            ccdGain=1.0,
                                                            plots=trackplots,\
    #                annulusOuterRadiusFactor=4.0, annulusInnerRadiusFactor=1.5,\
                    annulusOuterRadiusFactor=2.5, annulusInnerRadiusFactor=1.1,\
                    sigmaclipping=True, returnsubtractedflux=True)
            return xCenter, yCenter, flux, error
    
        for j in range(len(starpositions)):
            deltax, deltay = starpositions[j]
            #init_x = starpositions[j][0]
            #init_y = starpositions[j][1]
            RAtodeg = float(ephem.hours(starpositions[j][0]))*180./np.pi
            Dectodeg = float(ephem.degrees(starpositions[j][1]))*180./np.pi
            init_x, init_y = w.wcs_world2pix(RAtodeg, Dectodeg, 0)
            
            if trackplots:
                fig = plt.figure(num=None, figsize=(18, 3), facecolor='w', edgecolor='k')
                fig.subplots_adjust(wspace=0.5)
                subplotsDimensions = 140
                photSubplotsOffset = 3
                statusSubplotOffset = 6
                plottingThings = [fig,subplotsDimensions,photSubplotsOffset]
            else: 
                plottingThings = False
            xCenter, yCenter, flux, error = trackandphot(init_x, init_y, plottingThings)
            
            fluxes[i][j] = flux
            errors[i][j] = error
            xcentroids[i][j] = xCenter
            ycentroids[i][j] = yCenter
            if trackplots:
                plt.show()
    
    lc = fluxes[:,0]/fluxes[:,1]
    lc_err = lc*np.sqrt( (errors[:,0]/fluxes[:,0])**2 + (errors[:,1]/fluxes[:,1])**2 )
    lcmedian = np.median(lc)
    lc /= lcmedian
    lc_err /= lcmedian
    
    # remove airmass
    def fitfunc(p):
        return p[0]*airmass + p[1]
    
    def errfunc(p):
        return (fitfunc(p) - lc)/lc_err
    
    initp = [1, 1]
    bestp = optimize.leastsq(errfunc, initp)[0]
    correctedlc = lc - fitfunc(bestp) + 1
    stds[ap] = np.std(correctedlc)

#plt.plot(apertureradii, stds)
#np.savetxt('GD165lightcurve.txt',np.vstack([times, correctedlc]))
fig, ax = plt.subplots(figsize=(12,12))

alltimes = np.copy(times)
goodpoints = np.abs(correctedlc - np.median(correctedlc)) < 0.4#*np.std(correctedlc)
#times = times[goodpoints]
#fluxes = fluxes[goodpoints]
#correctedlc = correctedlc[goodpoints]
#lc_err = lc_err[goodpoints]

## Taken from git/research/koi351/variableaperture.ipynb
oot = np.ones_like(times).astype(bool)
target = fluxes[:,0]
targetOOT = fluxes[oot,0]
compStars = fluxes[:,1:]
compStarsOOT = fluxes[oot,1:]

numCompStars = np.shape(fluxes[:,1:])[1]
initP = np.zeros([numCompStars])+ 1./numCompStars
########################################################

def errfunc(p,target):
    if all(p >=0.0): 
        return np.dot(p,compStarsOOT.T) - targetOOT ## Find only positive coefficients

#return np.dot(p,compStarsOOT.T) - target
bestFitP = optimize.leastsq(errfunc,initP[:],args=(targetOOT.astype(np.float64)),\
                            maxfev=10000000,epsfcn=np.finfo(np.float32).eps)[0]
print '\nDefault weight:',1./numCompStars
print 'Best fit regression coefficients:',bestFitP

#self.comparisonStarWeights = np.vstack([compStarKeys,bestFitP])
meanComparisonStar = np.dot(bestFitP,compStars.T)

meanCompError = np.zeros_like(meanComparisonStar)
for i in range(1,len(fluxes[0,:])):
    meanCompError += (bestFitP[i-1]*fluxes[:,i]/np.sum(bestFitP[i-1]*fluxes[:,i]))**2 *(errors[:,i]/fluxes[:,i])**2
meanCompError = meanComparisonStar*np.sqrt(meanCompError)

lc = fluxes[:,0]/meanComparisonStar
lc /= np.median(lc)
lc_err = lc*np.sqrt((errors[:,0]/fluxes[:,0])**2 + (meanCompError/meanComparisonStar)**2)
reducedchi2 = np.sum(((lc-np.ones_like(lc))/lc_err)**2)/(len(lc) - 1)

sigmaclipstd = np.std(lc[np.abs(lc - np.median(lc)) < 4*np.std(lc)])

#clcmean = np.mean(correctedlc)
#ax.plot(times, lc,'.')#, yerr=lc_err, fmt='.')
#import matplotlib
#matplotlib.rcParams['font.size'] = 15
minint = int(np.min(times))
times_offset = times - minint
ax.errorbar(times_offset, lc, yerr=lc_err, fmt='.', color='k')
ax.set_xlabel('JD - %d' % minint)
ax.set_ylabel('Flux')
ax.axhline(1,ls='--',color='k')
ax.set_title('WD 2218+706: $\chi^2'+(' = %.3f' % reducedchi2)+r',\sigma = '+('%.3f$' % sigmaclipstd))
fig.savefig('plots/WD2218+706_lc.pdf')
#ax.set_ylim([0.8,1.2])
#plt.show()


difftimes = np.diff(times)
bigjumps = np.abs(difftimes - np.median(difftimes)) > 2*np.std(difftimes)
bigjumpinds = np.arange(len(difftimes))[bigjumps]
bigjumpinds = np.concatenate([[0],bigjumpinds,[len(times)]])
#fig, ax = plt.subplots(2, 2, figsize=(14,10))
fig = plt.figure(figsize=(14,10))
for i in range(len(bigjumpinds)-1):
    ax = fig.add_subplot(2,3,i+1)
    plttimes = times[1+bigjumpinds[i]:bigjumpinds[i+1]]
    minint = int(np.min(plttimes))
    plttimes = plttimes - minint
    pltlc = lc[1+bigjumpinds[i]:bigjumpinds[i+1]]
    pltlc_err = lc_err[1+bigjumpinds[i]:bigjumpinds[i+1]]
    ax.errorbar(plttimes, pltlc, yerr=pltlc_err, fmt='k.')
    ax.set_xlabel('JD - %d' % minint)
    ax.set_ylabel('Flux')
fig.tight_layout()
plt.show()




    