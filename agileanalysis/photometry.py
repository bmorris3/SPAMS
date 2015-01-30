'''oscaar v2.0 
   Module for differential photometry
   Developed by Brett Morris, 2011-2013'''
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from scipy import optimize

def sigmaclip(data, Nsigma=4):
    passthese = (np.abs(data - np.median(data)) < Nsigma*np.std(data))
    return data[passthese]

def fitgaussian(x, y, initparams):
    '''
    Return the best fit parameters `amplitude`, `x0` and `sigma` to a fit
    to a guassian for one dimensional data set `x` and `y`, with initial 
    parameters in list `initparams`.
    '''
    gaussian = lambda x, amplitude, x0, sigma: amplitude*np.exp(-(x-x0)**2/(2*sigma**2))
    errfunc = lambda p, x: gaussian(x, p[0], p[1], p[2]) - y
    bestFitP = optimize.leastsq(errfunc, initparams[:], args=(x),\
                                maxfev=10000000, epsfcn=np.finfo(np.float32).eps)[0]
    return bestFitP

#def genericPrayerBead(fitfunction, errfunction, bestFitP, x, bestFitModel, residuals, errors, plots=False):
#	"""Was meant to be way more generic than it is because of the leastsq 'args' parameter """
#	PBparameterTraces = np.zeros([len(bestFitModel),len(bestFitP)])
#	if plots: plt.ion()
#	for i in range(len(bestFitModel)):
#		if i == 0: 
#			modelPlusResiduals = bestFitModel + residuals
#			shiftedSigmas = errors
#		else: 
#			modelPlusResiduals = bestFitModel + shiftedResiduals
#	
#		data = np.copy(modelPlusResiduals) 	## Add the shifted residuals to the best fit model
#		bestFitP =  (1.0*np.array(bestFitP,dtype=np.float64)).tolist()
#		
#		args = (np.array(modelPlusResiduals),errors)
#			
#		PBiterationBestFitPs = optimize.leastsq(errfunction,bestFitP[:],args=args,\
#				epsfcn=1*np.finfo(np.float64).eps,xtol=1*np.finfo(np.float64).eps,maxfev=100*100*(len(data)+1))[0]
#		PBparameterTraces[i,:] = PBiterationBestFitPs	## record the best fit parameters
#		shiftedResiduals = np.roll(residuals,i)		## shift the residuals over one, repeat
#		shiftedSigmas = np.roll(errors,i)
#		if plots:
#			plt.clf()
#			plt.plot(x,data,'o')
#			plt.plot(x,fitfunction(PBiterationBestFitPs,x))
#			#plt.plot()
#			plt.draw()
#	if plots: plt.ioff()
#	uncertainties_std = np.std(PBparameterTraces,axis=0)	## Std of the best fits for each param is ~ the uncertainty on each param
#	uncertainties_env = np.max(PBparameterTraces,axis=0) - np.min(PBparameterTraces,axis=0)
#	## measure PB uncertainties using the range of values and the std, and take the larger of the two
#	uncertainties = np.max(np.vstack([uncertainties_std,uncertainties_env]),axis=0)
#	
#	return PBparameterTraces, uncertainties


def phot(image, xCentroid, yCentroid, apertureRadius, plottingThings, annulusOuterRadiusFactor=2.8, annulusInnerRadiusFactor=1.40, ccdGain=1, sigmaclipping=False, plots=False, returnsubtractedflux=False):
    '''
    Method for aperture photometry. 
    
    Parameters
    ----------
    image : numpy.ndarray
        FITS image opened with PyFITS
    
    xCentroid : float
        Stellar centroid along the x-axis (determined by trackSmooth or equivalent)
                
    yCentroid : float
        Stellar centroid along the y-axis (determined by trackSmooth or equivalent)
                
    apertureRadius : float
        Radius in pixels from centroid to use for source aperture
                     
    annulusInnerRadiusFactor : float
        Measure the background for sky background subtraction fron an annulus from a factor of 
        `annulusInnerRadiusFactor` bigger than the `apertureRadius` to one a factor `annulusOuterRadiusFactor` bigger.
    
    annulusOuterRadiusFactor : float
        Measure the background for sky background subtraction fron an annulus a factor of 
        `annulusInnerRadiusFactor` bigger than the `apertureRadius` to one a factor `annulusOuterRadiusFactor` bigger.
                          
    ccdGain : float
        Gain of your detector, used to calculate the photon noise
    
    plots : bool
            If `plots`=True, display plots showing the aperture radius and 
            annulus radii overplotted on the image of the star
                   
    Returns
    -------
    rawFlux : float
        The background-subtracted flux measured within the aperture
    
    rawError : float
        The photon noise (limiting statistical) Poisson uncertainty on the measurement of `rawFlux`
    
    errorFlag : bool
        Boolean corresponding to whether or not any error occured when running oscaar.phot(). If an error occured, the flag is
        True; otherwise False.
               
     Core developer: Brett Morris (NASA-GSFC)
    '''
    if plots:
        [fig,subplotsDimensions,photSubplotsOffset] = plottingThings
        if photSubplotsOffset == 0: plt.clf()
    annulusRadiusInner = annulusInnerRadiusFactor*apertureRadius 
    annulusRadiusOuter = annulusOuterRadiusFactor*apertureRadius

    ## From the full image, cut out just the bit around the star that we're interested in
    imageCrop = image[xCentroid-annulusRadiusOuter+1:xCentroid+annulusRadiusOuter+2,yCentroid-annulusRadiusOuter+1:yCentroid+annulusRadiusOuter+2]
    [dimy,dimx] = imageCrop.shape
    XX, YY = np.meshgrid(np.arange(dimx),np.arange(dimy))    
    x = (XX - annulusRadiusOuter)**2
    y = (YY - annulusRadiusOuter)**2
    ## Assemble arrays marking the pixels marked as either source or background pixels
    sourceIndices = x + y <= apertureRadius**2
    skyIndices = (x + y <= annulusRadiusOuter**2)*(x + y >= annulusRadiusInner**2)

    if sigmaclipping:    
        clippedbackground = np.median(sigmaclip(imageCrop[skyIndices]))

        rawFlux = np.sum(imageCrop[sourceIndices] - clippedbackground)*ccdGain
        rawError = np.sqrt(np.sum(imageCrop[sourceIndices]*ccdGain) + ccdGain*clippedbackground) ## Poisson-uncertainty
    else:
        rawFlux = np.sum(imageCrop[sourceIndices] - np.median(imageCrop[skyIndices]))*ccdGain
        rawError = np.sqrt(np.sum(imageCrop[sourceIndices]*ccdGain) + np.median(ccdGain*imageCrop[skyIndices])) ## Poisson-uncertainty

#    fig = plt.figure()
#    cutoff = 4*np.std(imageCrop[skyIndices])
#    clipped = sigmaclip(imageCrop[skyIndices])
#    plt.hist(imageCrop[skyIndices],1000,facecolor='r',histtype='stepfilled',alpha=0.2)
#    plt.hist(clipped,1000,facecolor='w',histtype='stepfilled')
#    plt.axvline(ymin=0,ymax=1,x=np.median(imageCrop[skyIndices])+cutoff)
#    plt.axvline(ymin=0,ymax=1,x=np.median(imageCrop[skyIndices])-cutoff)
#    plt.axvline(ymin=0,ymax=1,x=np.median(imageCrop[skyIndices]),color='g')
#    plt.axvline(ymin=0,ymax=1,x=np.median(clipped),color='m')
#    plt.title('Difference in medians: %f' % ((np.median(imageCrop[skyIndices]) - np.median(clipped))/np.median(imageCrop[skyIndices])))
#    plt.show()

    if plots:
        def format_coord(x, y):
            ''' Function to also give data value on mouse over with imshow. '''
            col = int(x+0.5)
            row = int(y+0.5)
            try:
                return 'x=%i, y=%i, Flux=%1.1f' % (x, y, imageCrop[row,col])
            except:
                return 'x=%i, y=%i' % (x, y)
       
        med = np.median(imageCrop)
        dsig = np.std(imageCrop)
        
        ax = fig.add_subplot(subplotsDimensions+photSubplotsOffset+1)
        ax.imshow(imageCrop, cmap=cm.gray, interpolation="nearest",vmin = med-0.5*dsig, vmax =med+2*dsig)
       
        theta = np.arange(0,360)*(np.pi/180)
        rcos = lambda r, theta: annulusRadiusOuter + r*np.cos(theta)
        rsin = lambda r, theta: annulusRadiusOuter + r*np.sin(theta)
        ax.plot(rcos(apertureRadius,theta),rsin(apertureRadius,theta),'m',linewidth=4)
        ax.plot(rcos(annulusRadiusInner,theta),rsin(annulusRadiusInner,theta),'r',linewidth=4)
        ax.plot(rcos(annulusRadiusOuter,theta),rsin(annulusRadiusOuter,theta),'r',linewidth=4)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_title('Aperture')
        ax.set_xlim([-.5,dimx-.5])
        ax.set_ylim([-.5,dimy-.5])
        ax.format_coord = format_coord 
        plt.draw()
    # new feature for diagnostics
    if not returnsubtractedflux: 
        return [rawFlux, rawError, False]
    else: 
        if sigmaclipping: 
            return [rawFlux, rawError, False, clippedbackground]
        else: 
            return [rawFlux, rawError, False, np.median(imageCrop[skyIndices])]

def multirad(image, xCentroid, yCentroid, apertureRadii, plottingThings, annulusOuterRadiusFactor=2.8, annulusInnerRadiusFactor=1.40, ccdGain=1, sigmaclipping=False, plots=False):
    '''
    Method for aperture photometry. 
    
    Parameters
    ----------
    image : numpy.ndarray
        FITS image opened with PyFITS
    
    xCentroid : float
        Stellar centroid along the x-axis (determined by trackSmooth or equivalent)
                
    yCentroid : float
        Stellar centroid along the y-axis (determined by trackSmooth or equivalent)
                
    apertureRadii : list
        List of aperture radii (floats) to feed to phot().
                     
    annulusInnerRadiusFactor : float
        Measure the background for sky background subtraction fron an annulus from a factor of 
        `annulusInnerRadiusFactor` bigger than the `apertureRadius` to one a factor `annulusOuterRadiusFactor` bigger.
    
    annulusOuterRadiusFactor : float
        Measure the background for sky background subtraction fron an annulus a factor of 
        `annulusInnerRadiusFactor` bigger than the `apertureRadius` to one a factor `annulusOuterRadiusFactor` bigger.
                          
    ccdGain : float
        Gain of your detector, used to calculate the photon noise
    
    plots : bool
            If `plots`=True, display plots showing the aperture radius and 
            annulus radii overplotted on the image of the star
                   
    Returns
    -------
    rawFlux : float
        The background-subtracted flux measured within the aperture
    
    rawError : float
        The photon noise (limiting statistical) Poisson uncertainty on the measurement of `rawFlux`
    
    errorFlag : bool
        Boolean corresponding to whether or not any error occured when running oscaar.phot(). If an error occured, the flag is
        True; otherwise False.
               
     Core developer: Brett Morris (NASA-GSFC)
    '''

    #[apertureRadiusMin, apertureRadiusMax, apertureRadiusStep] = apertureRadiusSettings
    #apertureRadii = np.arange(apertureRadiusMin, apertureRadiusMax, apertureRadiusStep)

    fluxes = []
    errors = []
    photFlags = []
    for apertureRadius in apertureRadii:
        flux, error, photFlag = phot(image, xCentroid, yCentroid, apertureRadius, \
                                plottingThings, annulusOuterRadiusFactor=annulusOuterRadiusFactor, \
                                annulusInnerRadiusFactor=annulusInnerRadiusFactor, \
                                ccdGain=ccdGain, plots=False, sigmaclipping=sigmaclipping)
        fluxes.append(flux)
        errors.append(error)
        photFlags.append(photFlag)
    annulusRadiusOuter = annulusOuterRadiusFactor*np.max(apertureRadii)
    imageCrop = image[xCentroid-annulusRadiusOuter+1:xCentroid+annulusRadiusOuter+2,yCentroid-annulusRadiusOuter+1:yCentroid+annulusRadiusOuter+2]
    [dimy,dimx] = imageCrop.shape

    if plots:
        [fig,subplotsDimensions,photSubplotsOffset] = plottingThings
        if photSubplotsOffset == 0: plt.clf()
        def format_coord(x, y):
            ''' Function to also give data value on mouse over with imshow. '''
            col = int(x+0.5)
            row = int(y+0.5)
            try:
                return 'x=%i, y=%i, Flux=%1.1f' % (x, y, imageCrop[row,col])
            except:
                return 'x=%i, y=%i' % (x, y)
       
        med = np.median(imageCrop)
        dsig = np.std(imageCrop)
        
        ax = fig.add_subplot(subplotsDimensions+photSubplotsOffset+1)
        ax.imshow(imageCrop, cmap=cm.gray, interpolation="nearest",vmin = med-0.5*dsig, vmax =med+2*dsig)
       
        theta = np.arange(0,360)*(np.pi/180)
        rcos = lambda r, theta: annulusRadiusOuter + r*np.cos(theta)
        rsin = lambda r, theta: annulusRadiusOuter + r*np.sin(theta)
        for apertureRadius in apertureRadii:
            ax.plot(rcos(apertureRadius,theta),rsin(apertureRadius,theta),linewidth=4)
        #ax.plot(rcos(annulusRadiusInner,theta),rsin(annulusRadiusInner,theta),'r',linewidth=4)
        #ax.plot(rcos(annulusRadiusOuter,theta),rsin(annulusRadiusOuter,theta),'r',linewidth=4)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_title('Aperture')
        ax.set_xlim([-.5,dimx-.5])
        ax.set_ylim([-.5,dimy-.5])
        ax.format_coord = format_coord 
        plt.draw()            
    return fluxes, errors, photFlags

def phot_nobgsubtraction(image, xCentroid, yCentroid, apertureRadius, plottingThings, ccdGain=1, plots=False):
        if plots:
            [fig,subplotsDimensions,photSubplotsOffset] = plottingThings
            if photSubplotsOffset == 0: plt.clf()
        #annulusRadiusInner = annulusInnerRadiusFactor*apertureRadius 
        #annulusRadiusOuter = annulusOuterRadiusFactor*apertureRadius
    
        ## From the full image, cut out just the bit around the star that we're interested in
        annulusRadiusOuter = apertureRadius
        imageCrop = image[xCentroid-annulusRadiusOuter+1:xCentroid+annulusRadiusOuter+2,yCentroid-annulusRadiusOuter+1:yCentroid+annulusRadiusOuter+2]
        [dimy,dimx] = imageCrop.shape
        XX, YY = np.meshgrid(np.arange(dimx),np.arange(dimy))    
        x = (XX - annulusRadiusOuter)**2
        y = (YY - annulusRadiusOuter)**2
        ## Assemble arrays marking the pixels marked as either source or background pixels
        sourceIndices = x + y <= apertureRadius**2
        #print sum(sum(sourceIndices))
        Npixels = sum(sum(sourceIndices))
        #skyIndices = (x + y <= annulusRadiusOuter**2)*(x + y >= annulusRadiusInner**2)
    
    #    if sigmaclipping:    
    #        clippedbackground = np.median(sigmaclip(imageCrop[skyIndices]))
    #
    #        rawFlux = np.sum(imageCrop[sourceIndices] - clippedbackground)*ccdGain
    #        rawError = np.sqrt(np.sum(imageCrop[sourceIndices]*ccdGain) + ccdGain*clippedbackground) ## Poisson-uncertainty
    #    else:
    #    rawFlux = np.sum(imageCrop[sourceIndices] - np.median(imageCrop[skyIndices]))*ccdGain
    #    rawError = np.sqrt(np.sum(imageCrop[sourceIndices]*ccdGain) + np.median(ccdGain*imageCrop[skyIndices])) ## Poisson-uncertainty
        rawFlux = np.sum(imageCrop[sourceIndices])*ccdGain
        rawError = np.sqrt(np.sum(imageCrop[sourceIndices]*ccdGain)) ## Poisson-uncertainty
    
        return [rawFlux, rawError, False, Npixels]#, np.median(imageCrop[skyIndices])]

def phot_tmp(image, xCentroid, yCentroid, apertureRadius, backgroundarray, plottingThings, annulusOuterRadiusFactor=2.8, annulusInnerRadiusFactor=1.40, ccdGain=1, sigmaclipping=False, plots=False, returnsubtractedflux=False):
    '''
    Method for aperture photometry. 
    
    phot_tmp specific features
    --------------------------
    Extra tricks for removing background, i.e. removing
    the mean, median, gaussian fit to the histogram of the pixels in the background
    annulus, etc. 
    
    Also does gaussian fit to the PSF and returns the width (sigma)    
    
    Takes a boolean array the same size as the image, which indicates whether
    each pixel is part of the source or part of the background, determined by
    a cross correlation with the image
    
    Parameters
    ----------
    image : numpy.ndarray
        FITS image opened with PyFITS
    
    xCentroid : float
        Stellar centroid along the x-axis (determined by trackSmooth or equivalent)
                
    yCentroid : float
        Stellar centroid along the y-axis (determined by trackSmooth or equivalent)
                
    apertureRadius : float
        Radius in pixels from centroid to use for source aperture
                     
    annulusInnerRadiusFactor : float
        Measure the background for sky background subtraction fron an annulus from a factor of 
        `annulusInnerRadiusFactor` bigger than the `apertureRadius` to one a factor `annulusOuterRadiusFactor` bigger.
    
    annulusOuterRadiusFactor : float
        Measure the background for sky background subtraction fron an annulus a factor of 
        `annulusInnerRadiusFactor` bigger than the `apertureRadius` to one a factor `annulusOuterRadiusFactor` bigger.
                          
    ccdGain : float
        Gain of your detector, used to calculate the photon noise
    
    plots : bool
            If `plots`=True, display plots showing the aperture radius and 
            annulus radii overplotted on the image of the star
                   
    Returns
    -------
    rawFlux : float
        The background-subtracted flux measured within the aperture
    
    rawError : float
        The photon noise (limiting statistical) Poisson uncertainty on the measurement of `rawFlux`
    
    errorFlag : bool
        Boolean corresponding to whether or not any error occured when running oscaar.phot(). If an error occured, the flag is
        True; otherwise False.
               
     Core developer: Brett Morris (NASA-GSFC)
    '''
    if plots:
        [fig,subplotsDimensions,photSubplotsOffset] = plottingThings
        if photSubplotsOffset == 0: plt.clf()
    annulusRadiusInner = annulusInnerRadiusFactor*apertureRadius 
    annulusRadiusOuter = annulusOuterRadiusFactor*apertureRadius

    ## From the full image, cut out just the bit around the star that we're interested in
    imageCrop = image[xCentroid-annulusRadiusOuter+1:xCentroid+annulusRadiusOuter+2,yCentroid-annulusRadiusOuter+1:yCentroid+annulusRadiusOuter+2]
    [dimy,dimx] = imageCrop.shape
    XX, YY = np.meshgrid(np.arange(dimx),np.arange(dimy))    
    x = (XX - annulusRadiusOuter)**2
    y = (YY - annulusRadiusOuter)**2
    ## Assemble arrays marking the pixels marked as either source or background pixels
    sourceIndices = x + y <= apertureRadius**2
    skyIndices = (x + y <= annulusRadiusOuter**2)*(x + y >= annulusRadiusInner**2)

#    backgroundsubtractionfunction = lambda x: np.median(x)
#    def backgroundsubtractionfunction(rawdata):
#        n, xedges = np.histogram(rawdata, 1000)
#        smallbins = [np.mean([xedges[i], xedges[i+1]]) for i in range(len(xedges)-1)]
#        gaussian = lambda x, amplitude, x0, sigma: amplitude*np.exp(-(x-x0)**2/(2*sigma**2))
#        errfunc = lambda p, x: gaussian(x, p[0], p[1], p[2]) - n
#        initP = [np.max(n), 0.03, 0.01]
#        bestFitP = optimize.leastsq(errfunc,initP[:], args=(smallbins),maxfev=10000000,epsfcn=np.finfo(np.float32).eps)[0]
#        #print bestFitP
#        x0 = bestFitP[1]
#        return x0

    ## Use cross-correlation's identification of background pixels to make
    ## sure that no sources are accidentally incorporated
    backgroundcrop = backgroundarray[xCentroid-annulusRadiusOuter+1:xCentroid+annulusRadiusOuter+2,yCentroid-annulusRadiusOuter+1:yCentroid+annulusRadiusOuter+2]
    truebackgroundpixels = backgroundcrop * skyIndices
    background = imageCrop[truebackgroundpixels]
    #plt.show()

    clippedbackground = np.median(background)

    #def backgroundsubtractionfunction(background):
    #clippedbackground = backgroundsubtractionfunction(sigmaclip(imageCrop[skyIndices], Nsigma=3))


    rawFlux = np.sum(imageCrop[sourceIndices] - clippedbackground)*ccdGain
    rawError = np.sqrt(np.sum(imageCrop[sourceIndices]*ccdGain) + ccdGain*clippedbackground) ## Poisson-uncertainty

    #############################################
    # Measure the PSF width in each dimension
    tinyimageCrop = image[xCentroid-apertureRadius+1:xCentroid+apertureRadius+2, yCentroid-apertureRadius+1:yCentroid+apertureRadius+2]
#    plt.imshow(tinyimageCrop)
    xsum_tinyimage = np.sum(tinyimageCrop, axis=0)
    ysum_tinyimage = np.sum(tinyimageCrop, axis=1)
    
    amplitude_x, x0_x, PSFsigma_x = fitgaussian(range(len(xsum_tinyimage)), xsum_tinyimage, [np.max(xsum_tinyimage), len(xsum_tinyimage)/2., 1])
    amplitude_y, x0_y, PSFsigma_y = fitgaussian(range(len(ysum_tinyimage)), ysum_tinyimage, [np.max(ysum_tinyimage), len(ysum_tinyimage)/2., 1])
    #amplitude = 0.5*(amplitude_x + amplitude_y)
    #x0 = 0.5*(amplitude_x + amplitude_y)
    PSFsigma = 0.5*(PSFsigma_x + PSFsigma_y)
    #print PSFsigma,'= sigma'
    #print PSFsigma*2.354,'= FWHM'
    #plt.plot(xsum_tinyimage)
    #plt.plot(ysum_tinyimage)
    #plt.show()


#    fig = plt.figure()
#    cutoff = 3*np.std(imageCrop[skyIndices])
#    clipped = sigmaclip(imageCrop[skyIndices], Nsigma=3)
#    plt.hist(imageCrop[skyIndices],1000,facecolor='r',histtype='stepfilled',alpha=0.2)
#    n, xedges, _ = plt.hist(clipped,1000,facecolor='w',histtype='stepfilled')
#
#    smallbins = [np.mean([xedges[i], xedges[i+1]]) for i in range(len(xedges)-1)]
#    gaussian = lambda x, amplitude, x0, sigma: amplitude*np.exp(-(x-x0)**2/(2*sigma**2))
#    errfunc = lambda p, x: gaussian(x, p[0], p[1], p[2]) - n
#    initP = [np.max(n), 0.08, 0.02]
#    bestFitP = optimize.leastsq(errfunc,initP[:], args=(smallbins),maxfev=10000000,epsfcn=np.finfo(np.float32).eps)[0]
#    x0 = bestFitP[1]
#    plt.plot(smallbins, gaussian(smallbins, bestFitP[0], bestFitP[1], bestFitP[2]), 'r', lw=2)
#
#    plt.axvline(ymin=0,ymax=1,x=np.median(imageCrop[skyIndices])+cutoff)
#    plt.axvline(ymin=0,ymax=1,x=np.median(imageCrop[skyIndices])-cutoff)
#   # plt.axvline(ymin=0,ymax=1,x=np.median(imageCrop[skyIndices]),color='g')
#   # plt.axvline(ymin=0,ymax=1,x=np.median(clipped),color='m')
##    plt.title('Difference in medians: %f' % ((np.median(imageCrop[skyIndices]) - np.median(clipped))/np.median(imageCrop[skyIndices])))
#    plt.axvline(ymin=0,ymax=1,x=np.median(clipped),color='b')
#    plt.axvline(ymin=0,ymax=1,x=np.mean(clipped),color='m')
#    #plt.title('Difference medians, mean: %f' % ((np.mean(clipped) - np.median(clipped))/np.median(imageCrop[skyIndices])))
#    plt.title('Difference medians, fit: %f' % ((x0 - np.median(clipped))/np.median(clipped)))
#    plt.show()

#    plt.clf()

    #plt.show()

    if plots:
        def format_coord(x, y):
            ''' Function to also give data value on mouse over with imshow. '''
            col = int(x+0.5)
            row = int(y+0.5)
            try:
                return 'x=%i, y=%i, Flux=%1.1f' % (x, y, imageCrop[row,col])
            except:
                return 'x=%i, y=%i' % (x, y)
       
        med = np.median(imageCrop)
        dsig = np.std(imageCrop)
        
        ax = fig.add_subplot(subplotsDimensions+photSubplotsOffset+1)
        ax.imshow(imageCrop, cmap=cm.gray, interpolation="nearest",vmin = med-0.5*dsig, vmax =med+2*dsig)
       
        theta = np.arange(0,360)*(np.pi/180)
        rcos = lambda r, theta: annulusRadiusOuter + r*np.cos(theta)
        rsin = lambda r, theta: annulusRadiusOuter + r*np.sin(theta)
        ax.plot(rcos(apertureRadius,theta),rsin(apertureRadius,theta),'m',linewidth=4)
        ax.plot(rcos(annulusRadiusInner,theta),rsin(annulusRadiusInner,theta),'r',linewidth=4)
        ax.plot(rcos(annulusRadiusOuter,theta),rsin(annulusRadiusOuter,theta),'r',linewidth=4)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_title('Aperture')
        ax.set_xlim([-.5,dimx-.5])
        ax.set_ylim([-.5,dimy-.5])
        ax.format_coord = format_coord 
        plt.draw()
    # new feature for diagnostics
    if not returnsubtractedflux: 
        return [rawFlux, rawError, False]
    else: 
        if sigmaclipping: 
            return [rawFlux, rawError, False, clippedbackground, PSFsigma]
        else: 
            return [rawFlux, rawError, False, np.median(imageCrop[skyIndices]), PSFsigma]

