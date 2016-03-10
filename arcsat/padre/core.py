# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os

import astropy.units as u
from astropy.io import fits
from astropy.time import Time
from astropy import wcs
from astropy.coordinates import Angle
from astropy.utils.console import ProgressBar

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

from .cals import master_dark as calculate_master_dark
from .cals import master_flat as calculate_master_flat
from .photometry import phot, track_smooth
from .lightcurve import LightCurve

__all__ = ["Observation"]


class Observation(object):
    """
    Container for one observation
    """
    def __init__(self, dark_paths=None, flat_paths=None, data_image_paths=None,
                 output_dir_path=None, target_name=None, first_obs_date=None,
                 regions_file_path=None):
        """
        Parameters
        ----------
        dark_paths : string
            Paths to dark frames            
        flat_paths : string
            Paths to flat fields
        data_images_paths : list
            Paths to the data images
        output_dir_path : string
            Path to the data output directory
        target_name : string
            Name of the target
        first_obs_date : string
            string representation of the first date of observations
        regions_file_path : string
            Path to regions file denoting positions of stars to track
        """
        self.dark_paths = dark_paths
        self.flat_paths = flat_paths
        self.data_image_paths = data_image_paths
        self.target_name = target_name
        self.first_obs_date = first_obs_date
        self.regions_file_path = regions_file_path
        self.master_dark = None
        self.master_flat = None
        self.initial_centroids = None
        self.cached_light_curves = dict()
        
        # Make output directory if one does not already exist there
        self.output_dir = os.path.join(output_dir_path, "{0}_{1}"
                                       .format(target_name, first_obs_date))
        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)
        
    def get_calibration_frames(self, force=False, master_flat_path=None,
                               master_dark_path=None, master_dark_scale_factor=1):
        """
        Get the master flat field, master dark frame.
        
        Parameters
        ----------
        force : bool
            Load existing master calibration frames if those files exist or
            ``force``==`False`, otherwise calculate calibration frames.

        master_flat_path : string or `None` (optional)
            Path to the master flat you'd like to use. If `None`, generate a
            master flat.

        master_dark_path : string or `None` (optional)
            Path to the master dark you'd like to use. If `None`, generate a
            master dark.

        master_dark_scale_factor : float
            Scale the master_dark frame by this factor. If, for example, you're
            using another target's master dark which had a different exposure
            length, you can set ``master_dark_scale_factor`` to the ratio of
            exposure times to scale the fluxes appropriately
        """
        print(master_flat_path)
        if master_flat_path is None:
            master_flat_path = os.path.join(self.output_dir, 'masterflat.fits')

        if master_dark_path is None:
            master_dark_path = os.path.join(self.output_dir, 'masterdark.fits')

        print(master_flat_path, 'exists={0}'.format(os.path.exists(master_flat_path)))

        if not os.path.exists(master_dark_path) or force:
            master_dark = calculate_master_dark(self.dark_paths,
                                                self.data_image_paths[-1])
            fits.writeto(master_dark_path, master_dark, clobber=True)
        else:
            master_dark = fits.getdata(master_dark_path)

        if not os.path.exists(master_flat_path) or force:
            print('make master flat')
            master_flat = calculate_master_flat(self.flat_paths,
                                                self.dark_paths)
            fits.writeto(master_flat_path, master_flat, clobber=True)
        else:
            master_flat = fits.getdata(master_flat_path)

        self.master_dark = master_dark*master_dark_scale_factor
        self.master_flat = master_flat

    def parse_regions_file(self):
        """
        Parse the contents of the regions file to get initial centroids of 
        stars.
        """
        regions_file_contents = open(self.regions_file_path, 
                                     'r').read().splitlines()
        self.initial_centroids = []
        for line in regions_file_contents:
            if line.startswith('circle'):
                y_center, x_center = line.split('(')[1].split(')')[0].split(',')[:2]
                self.initial_centroids.append([y_center, x_center])

    def generate_calibrated_images(self, output_dir, output_suffix='corrected'):
        """
        Dark subtract and flat field all data images, save them to `output_dir`.

        This is useful for producing images that can be passed to astrometry.net
        """
        n_images = len(self.data_image_paths)
        print("Generating calibrated images...")
        with ProgressBar(n_images) as bar:
            for i, image_path in enumerate(self.data_image_paths):
                bar.update(i)
                image_header = fits.getheader(image_path)
                image_data = ((fits.getdata(image_path) - self.master_dark) /
                              self.master_flat)

                # TODO: Should this line be in our out?
                #image_data[image_data < 0] = np.median(image_data)

                file_name = image_path.split(os.sep)[-1].split('.fits')[0]
                fits.writeto(os.path.join(output_dir,
                                          file_name+output_suffix+'.fits'),
                             image_data, clobber=True, header=image_header)

    def calculate_photometry(self, aperture_radii, calibrated=True,
                             track_plots=False):
        """
        Calculate photometry for all ``aperture_radii``.
        """
        if self.initial_centroids is None:
            self.parse_regions_file()

        n_images = len(self.data_image_paths)
        n_comparison_stars = len(self.initial_centroids)
        n_aperture_radii = len(aperture_radii)
        times = np.zeros(n_images)
        exp_times = u.Quantity(np.zeros(n_images), unit=u.min)
        fluxes = np.zeros((n_images, n_comparison_stars, n_aperture_radii))
        errors = np.zeros((n_images, n_comparison_stars, n_aperture_radii))
        xcentroids = np.zeros((n_images, n_comparison_stars))
        ycentroids = np.zeros((n_images, n_comparison_stars))
        airmass = np.zeros(n_images)
        means = np.zeros(n_images)
        medians = np.zeros(n_images)
        
        trackplots = track_plots
        
        # Loop over images
        print("Calculating photometry...")
        with ProgressBar(n_images) as bar:
            for i, image_path in enumerate(self.data_image_paths):
                bar.update(i)
                image_header = fits.getheader(image_path)
                if not calibrated:
                    image_data = ((fits.getdata(image_path) - self.master_dark) /
                                  self.master_flat)
                else:
                    image_data = fits.getdata(image_path)
                exp_times[i] = image_header['EXPTIME']*u.s
                jd_start_exposure = Time(image_header['JD'], format='jd')
                jd_mid_exposure = (jd_start_exposure + 0.5*exp_times[i]).jd
                times[i] = jd_mid_exposure
                means[i] = np.mean(image_data)
                medians[i] = np.median(image_data)
                airmass[i] = image_header['AIRMASS']
                hdulist = fits.open(image_path)
                w = wcs.WCS(hdulist[0].header, hdulist)

                # Loop over aperture radii
                for ap, apertureradius in enumerate(aperture_radii):

                    zoomfactor = 10#20
                    smoothconst = 3#2.5

                    for j in range(n_comparison_stars):
                        # If coordinates stored in sexigesimal notation, use WCS
                        # to find the target, otherwise use pixel positions.
                        if ':' in self.initial_centroids[j][0]:
                                RA = Angle(self.initial_centroids[j][0],
                                           unit=u.hourangle)
                                Dec = Angle(self.initial_centroids[j][1],
                                            unit=u.deg)
                                init_x, init_y = w.wcs_world2pix(RA.degree,
                                                                 Dec.degree,
                                                                 0)
                        elif i == 0:
                            init_x = self.initial_centroids[j][0]
                            init_y = self.initial_centroids[j][1]
                        else:
                            init_x = xcentroids[i-1, j]
                            init_y = ycentroids[i-1, j]

                        if trackplots:
                            fig = plt.figure(num=None, figsize=(18, 3), facecolor='w', edgecolor='k')
                            fig.subplots_adjust(wspace=0.5)
                            subplotsDimensions = 140
                            photSubplotsOffset = 3
                            statusSubplotOffset = 6
                            plottingThings = [fig,subplotsDimensions,photSubplotsOffset]
                        else:
                            plottingThings = False
                        try:
                            [xCenter,yCenter,averageRadius, _] = track_smooth(image_data,init_y, init_x, \
                                smoothconst, plottingThings, preCropped=False, zoom=zoomfactor, plots=trackplots)
                            flux, error, photFlags, _ = phot(image_data, xCenter, yCenter,
                                                                            apertureradius,
                                                                            plottingThings,
                                                                            ccdGain=1.0,
                                                                            plots=trackplots,\
                    #                annulusOuterRadiusFactor=4.0, annulusInnerRadiusFactor=1.5,\
                                    annulusOuterRadiusFactor=3, annulusInnerRadiusFactor=1.5,\
                                    sigmaclipping=True, returnsubtractedflux=True)
                        except (ValueError, IndexError):
                            no_value = np.nan
                            flux = no_value
                            error = no_value
                            xCenter = no_value
                            yCenter = no_value
                            print('ValueError raised at {0}'.format(times[i]))

                        fluxes[i, j, ap] = flux
                        errors[i, j, ap] = error
                        xcentroids[i, j] = xCenter
                        ycentroids[i, j] = yCenter
                        if trackplots:
                            plt.show()

        self.fluxes = fluxes
        self.errors = errors
        self.times = times
        self.xcentroids = xcentroids
        self.ycentroids = ycentroids
        self.aperture_radii = aperture_radii
        
    def light_curve(self, star_index=0):
        """
        Calculate light curve
        
        star_index : int
            Choose the target star index
        """
        light_curves = []
        for aperture_radius_index in range(self.fluxes.shape[2]):
            light_curve_key = (star_index, aperture_radius_index)

            if light_curve_key in self.cached_light_curves:
                return self.cached_light_curves[light_curve_key]
            else:
                comp_stars = np.arange(self.fluxes.shape[1]) != star_index
                target_fluxes = self.fluxes[:, star_index, aperture_radius_index]
                comp_star_fluxes = self.fluxes[:, comp_stars, aperture_radius_index]

                n_comp_stars = comp_star_fluxes.shape[1]
                initP = np.zeros([n_comp_stars])+ 1./n_comp_stars
                ########################################################

                def errfunc(p,target):
                    if all(p >=0.0):
                        return np.dot(p, comp_star_fluxes.T) - target_fluxes ## Find only positive coefficients

                #return np.dot(p,compStarsOOT.T) - target
                best_p = optimize.leastsq(errfunc, initP[:],
                                          args=(target_fluxes.astype(np.float64)),
                                                maxfev=10000000,
                                                epsfcn=np.finfo(np.float32).eps)[0]
                #print '\nDefault weight:',1./numCompStars
                #print 'Best fit regression coefficients:',bestFitP

                #self.comparisonStarWeights = np.vstack([compStarKeys,bestFitP])
                meanComparisonStar = np.dot(best_p, comp_star_fluxes.T)

                meanCompError = np.zeros_like(meanComparisonStar)
                comp_star_indices = np.arange(self.fluxes.shape[1])[comp_stars]
                for i in comp_star_indices:
                    meanCompError += ((best_p[i-1]*self.fluxes[:, i, aperture_radius_index] /
                                       np.nansum(best_p[i-1]*self.fluxes[:, i, aperture_radius_index]))**2 *
                                       (self.errors[:, i, aperture_radius_index] /
                                       self.fluxes[:, i, aperture_radius_index])**2)
                meanCompError = meanComparisonStar*np.sqrt(meanCompError)

                lc = self.fluxes[:, star_index, aperture_radius_index] / meanComparisonStar
                lc /= np.nanmedian(lc)
                lc_err = (lc*np.sqrt((self.errors[:, star_index, aperture_radius_index] /
                          self.fluxes[:, star_index, aperture_radius_index])**2 +
                          (meanCompError/meanComparisonStar)**2))
                self.cached_light_curves[light_curve_key] = (lc, lc_err)

                lc_name = ("{0}\n target={1}, aprad={2}"
                           .format(self.target_name,
                                   star_index,
                                   self.aperture_radii[aperture_radius_index]))

                # metadata = {}
                # metadata_attrs = ['xcentroids', 'ycentroids', 'fluxes', 'errors']
                # for attr in metadata_attrs:
                #     metadata[attr] = getattr(self, attr)
                light_curves.append(LightCurve(times=self.times, fluxes=lc,
                                               errors=lc_err, name=lc_name,
                                               raw_fluxes=self.fluxes,
                                               raw_errors=self.errors,
                                               xcentroids=self.xcentroids,
                                               ycentroids=self.ycentroids))
        return light_curves


