# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from astropy.time import Time
import astropy.units as u
import os
import numpy as np
import matplotlib.pyplot as plt

__all__ = ['LightCurve', 'split_light_curve_at_gaps', 
           'concatenate_light_curves']

class LightCurve(object):
    """
    Container object for light curves
    """
    def __init__(self, times=None, fluxes=None, errors=None, name=None,
                 meta=None):

        if isinstance(times[0], Time) and isinstance(times, np.ndarray):
            times = Time(times)
        elif not isinstance(times, Time):
            times = Time(times, format='jd')

        self.times = times
        self.fluxes = fluxes
        self.errors = errors
        self.name = name
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

    def __repr__(self):
        """Repr string should be 'LightCurve: name' """
        return "{0}: {1}".format(self.__class__.__name__, self.name)

    @property
    def times_jd(self):
        return self.times.jd

    def save_to(self, path, overwrite=True):
        """
        Save times, fluxes, errors to ``path``.txt
        """
        if not path.endswith('.txt'):
            path += '.txt'
        
        if not os.path.exists(path) or overwrite:
            attrs = ['times_jd', 'fluxes', 'errors']
            output_array = np.zeros((len(self.fluxes), len(attrs)), dtype=float)
            for i, attr in enumerate(attrs):
                output_array[:, i] = getattr(self, attr)
            np.savetxt(os.path.join(path), output_array)

    @classmethod
    def load_from(cls, path):
        """
        Load times, fluxes, errors to ``path``.txt
        """
        if not path.endswith('.txt'):
            path += '.txt'
        
        times, fluxes, errors = np.loadtxt(os.path.join(path), unpack=True)
        name = path.split(os.sep)[-1].split('.')[0]
        return cls(times=times, fluxes=fluxes, errors=errors, name=name)

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

