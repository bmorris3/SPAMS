# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from astropy.time import Time
import astropy.units as u
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy import optimize

__all__ = ['LightCurve', 'split_light_curve_at_gaps', 
           'concatenate_light_curves']

class LightCurve(object):
    """
    Container object for light curves
    """
    def __init__(self, times=None, fluxes=None, errors=None, name=None,
                 meta=None, raw_fluxes=None, raw_errors=None, xcentroids=None,
                 ycentroids=None, weights=None):

        if isinstance(times[0], Time) and isinstance(times, np.ndarray):
            times = Time(times)
        elif not isinstance(times, Time):
            times = Time(times, format='jd')

        self.times = times
        self.fluxes = fluxes
        self.errors = errors
        self.name = name
        self.meta = meta
        self.raw_fluxes = raw_fluxes
        self.raw_errors = raw_errors
        self.xcentroids = xcentroids
        self.ycentroids = ycentroids
        self.weights = weights
        self.cached_light_curves = {}
        self.meta = meta

    def plot(self, ax=None, show=True,
             plot_kwargs={'color':'b', 'marker':'o', 'lw':0}):
        """
        Plot light curve
        """
        if ax is None:
            ax = plt.gca()

        ax.plot(self.times.jd, self.fluxes, **plot_kwargs)
        ax.set(xlabel='Time', ylabel='Flux', title=self.name)

        if show:
            plt.show()

    def megaplot(self, show=True,
             plot_kwargs={'color':'b', 'marker':'o', 'lw':0}):
        """
        Plot light curve and other data
        """
        plt.figure(figsize=(15, 15))
        gs = GridSpec(3, 3)
        ax_lc = plt.subplot(gs[0, :])
        ax_rawflux = plt.subplot(gs[1, :], sharex=ax_lc)
        ax_centroids = plt.subplot(gs[2, 0])
        ax_x = plt.subplot(gs[2, 1])
        ax_y = plt.subplot(gs[2, 2])

        ax_lc.plot(self.times.jd, self.fluxes, **plot_kwargs)
        ax_lc.set(ylabel='Flux', title=self.name)

        ax_rawflux.plot(self.times.jd, self.raw_fluxes[:, 1:, 0])
        ax_rawflux.plot(self.times.jd, self.raw_fluxes[:, 0, 0], 'k.')
        ax_rawflux.set(xlabel='Time', ylabel='Flux')

        ax_centroids.plot(self.xcentroids[:, 1:], self.ycentroids[:, 1:], '.')
        ax_centroids.plot(self.xcentroids[:, 0], self.ycentroids[:, 0], 'ks')
        ax_centroids.set(xlabel='x', ylabel='y')

        ax_x.plot(self.times.jd, self.xcentroids[:, 1:])
        ax_x.plot(self.times.jd, self.xcentroids[:, 0], 'k.')
        ax_x.set(xlabel='Time', ylabel='x')

        ax_y.plot(self.times.jd, self.ycentroids[:, 1:])
        ax_y.plot(self.times.jd, self.ycentroids[:, 0], 'k.')
        ax_y.set(xlabel='Time', ylabel='y')

        if show:
            plt.show()

    def __repr__(self):
        """Repr string should be 'LightCurve: name' """
        return "<{0}: {1}>".format(self.__class__.__name__, self.name)

    @property
    def times_jd(self):
        return self.times.jd

    def save_to(self, path, overwrite=True):
        """
        Save times, fluxes, errors to ``path``.txt
        """
        if not path.endswith('.npz'):
            path += '.npz'
        
        if not os.path.exists(path) or overwrite:
            # attrs = ['times_jd', 'fluxes', 'errors']
            # output_array = np.zeros((len(self.fluxes), len(attrs)), dtype=float)
            # for i, attr in enumerate(attrs):
            #     output_array[:, i] = getattr(self, attr)
            # np.savetxt(os.path.join(path), output_array)
            save_attrs = ['times_jd', 'fluxes', 'errors', 'raw_fluxes',
                          'raw_errors', 'xcentroids', 'ycentroids']
            save_dict = {}
            for attr in save_attrs:
                save_dict[attr] = getattr(self, attr)
            np.savez_compressed(os.path.join(path), **save_dict)

    @classmethod
    def load_from(cls, path):
        """
        Load times, fluxes, errors to ``path``.txt
        """
        if not path.endswith('.npz'):
            path += '.npz'

        in_file = np.load(os.path.join(path))

        # times, fluxes, errors = np.loadtxt(os.path.join(path), unpack=True)
        # name = path.split(os.sep)[-1].split('.')[0]
        # return cls(times=times, fluxes=fluxes, errors=errors, name=name)
        return cls(times=in_file['times_jd'], fluxes=in_file['fluxes'],
                   errors=in_file['errors'], raw_fluxes=in_file['raw_fluxes'],
                   raw_errors=in_file['raw_errors'],
                   xcentroids=in_file['xcentroids'],
                   ycentroids=in_file['ycentroids'])

    def fit_linear_baseline(self, plots=False, sigma=5):
        """
        Find linear baseline trend with `~numpy.polyfit`. Do sigma clipping
        outlier rejection before fitting for the trend.
        """
        # Remove linear baseline trend
        if len(self.fluxes) > 1:
            order = 1
            print()
            non_outliers = (np.abs(self.fluxes - self.fluxes.mean()) < 
                            sigma*np.std(self.fluxes) + (self.fluxes < 
                            self.fluxes.mean()))
            
            linear_baseline = np.polyfit(self.times.jd[non_outliers],
                                         self.fluxes[non_outliers], order)
            linear_baseline_fit = np.polyval(linear_baseline, self.times.jd)
    
            if plots:
                fig, ax = plt.subplots(1, figsize=(8, 6))
                ax.axhline(1, ls='--', color='k')
                ax.plot(self.times.jd[non_outliers], self.fluxes[non_outliers], 'bo', label='data')
                ax.plot(self.times.jd[-non_outliers], self.fluxes[-non_outliers],'ro', label='outliers')
                ax.plot(self.times.jd, linear_baseline_fit, 'r', label='fit')
                ax.set(title='Linear trend fit', xlabel='JD', ylabel='Flux')
                ax.legend()
                plt.show()
        else:
            linear_baseline = None

        return linear_baseline


    def remove_linear_baseline(self, plots=False):
        """
        Find and remove linear baseline trend
        """
        linear_baseline = self.fit_linear_baseline(plots=plots)
        if linear_baseline is not None:
            linear_baseline_fit = np.polyval(linear_baseline, self.times.jd)
            self.fluxes =  self.fluxes/linear_baseline_fit
            self.errors = self.errors/linear_baseline_fit

    # def generate_light_curve(self, star_index=0):
    #
    #     for aperture_radius_index in range(self.fluxes.shape[2]):
    #         light_curve_key = (star_index, aperture_radius_index)
    #
    #         if light_curve_key in self.cached_light_curves:
    #             return self.cached_light_curves[light_curve_key]
    #         else:
    #             comp_stars = np.arange(self.fluxes.shape[1]) != star_index
    #             target_fluxes = self.fluxes[:, star_index, aperture_radius_index]
    #             comp_star_fluxes = self.fluxes[:, comp_stars, aperture_radius_index]
    #
    #             n_comp_stars = comp_star_fluxes.shape[1]
    #             initP = np.zeros([n_comp_stars])+ 1./n_comp_stars
    #             ########################################################
    #
    #             def errfunc(p,target):
    #                 if all(p >=0.0):
    #                     return np.dot(p, comp_star_fluxes.T) - target_fluxes ## Find only positive coefficients
    #
    #             #return np.dot(p,compStarsOOT.T) - target
    #             best_p = optimize.leastsq(errfunc, initP[:],
    #                                       args=(target_fluxes.astype(np.float64)),
    #                                             maxfev=10000000,
    #                                             epsfcn=np.finfo(np.float32).eps)[0]
    #             #print '\nDefault weight:',1./numCompStars
    #             #print 'Best fit regression coefficients:',bestFitP
    #
    #             #self.comparisonStarWeights = np.vstack([compStarKeys,bestFitP])
    #             meanComparisonStar = np.dot(best_p, comp_star_fluxes.T)
    #
    #             meanCompError = np.zeros_like(meanComparisonStar)
    #             comp_star_indices = np.arange(self.fluxes.shape[1])[comp_stars]
    #             for i in comp_star_indices:
    #                 meanCompError += ((best_p[i-1]*self.fluxes[:, i, aperture_radius_index] /
    #                                    np.sum(best_p[i-1]*self.fluxes[:, i, aperture_radius_index]))**2 *
    #                                    (self.errors[:, i, aperture_radius_index] /
    #                                    self.fluxes[:, i, aperture_radius_index])**2)
    #             meanCompError = meanComparisonStar*np.sqrt(meanCompError)
    #
    #             lc = self.fluxes[:, star_index, aperture_radius_index] / meanComparisonStar
    #             lc /= np.median(lc)
    #             lc_err = (lc*np.sqrt((self.errors[:, star_index, aperture_radius_index] /
    #                       self.fluxes[:, star_index, aperture_radius_index])**2 +
    #                       (meanCompError/meanComparisonStar)**2))
    #             self.cached_light_curves[light_curve_key] = (lc, lc_err)
    #
    #             lc_name = ("{0}\n target={1}, aprad={2}"
    #                        .format(self.target_name,
    #                                star_index,
    #                                self.aperture_radii[aperture_radius_index]))
    #
    #             # metadata = {}
    #             # metadata_attrs = ['xcentroids', 'ycentroids', 'fluxes', 'errors']
    #             # for attr in metadata_attrs:
    #             #     metadata[attr] = getattr(self, attr)
    #             return self.__class__.__init__(times=self.times, fluxes=lc,
    #                                            errors=lc_err, name=lc_name,
    #                                            raw_fluxes=self.fluxes,
    #                                            raw_errors=self.errors,
    #                                            xcentroids=self.xcentroids,
    #                                            ycentroids=self.ycentroids,
    #                                            weights=best_p)


def split_light_curve_at_gaps(light_curve, gap_median_factor=3):
    """
    Split a light curve by time gaps, into smaller light curves.
    """
    sort_by_time = np.argsort(light_curve.times.jd)
    dt = np.diff(light_curve.times.jd[sort_by_time])
    split_indices = np.argwhere(dt > gap_median_factor*np.median(dt)).T[0]
    split_times = [np.mean(light_curve.times.jd[sort_by_time][[split_index-1,
                                                               split_index+1]])
                   for split_index in split_indices]
    split_times = np.concatenate([[light_curve.times.jd.min()],
                                  split_times,
                                  [light_curve.times.jd.max()]])

    # # Visualize breaks in time
    # plt.plot(light_curve.times.jd)
    # [plt.axhline(splittime) for splittime in split_times]
    # plt.show()

    chunks = []

    for i in range(len(split_times) - 1):
        start_jd, end_jd = split_times[i], split_times[i+1]
        chunk = ((start_jd < light_curve.times.jd) &
                 (light_curve.times.jd < end_jd))
        chunks.append(LightCurve(times=light_curve.times.jd[chunk],
                                 fluxes=light_curve.fluxes[chunk],
                                 errors=light_curve.errors[chunk],
                                 name=light_curve.name+'_chunk{0}'.format(i)))
    return chunks


def concatenate_light_curves(light_curve_list, name=None):
    """
    Combine multiple light curves into one `LightCurve` object
    """
    times = []
    fluxes = []
    errors = []
    for light_curve in light_curve_list:
        times.append(light_curve.times.jd)
        fluxes.append(light_curve.fluxes)
        errors.append(light_curve.errors)
    times, fluxes, errors = [np.concatenate(i)
                             for i in [times, fluxes, errors]]

    times = Time(times, format='jd')
    return LightCurve(times=times, fluxes=fluxes, errors=errors,
                      name=name)

