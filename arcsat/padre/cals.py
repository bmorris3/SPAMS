# -*- coding: utf-8 -*-
"""
Calibration tools for ARCSAT.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from astropy.io import fits

__all__ = ['master_dark', 'master_flat']

def master_dark(dark_images_paths, exp_same_as_image_at_path):
    """
    Median-combine darks from list ``dark_images_paths`` with the same exposure 
    length as the file at path ``exp_same_as_image_at_path``.
        
    Parameters
    ----------
    dark_images_paths : list
        Paths to dark frames
        
    exp_same_as_image_at_path : string
        Only use dark frames with exposure length equivalent to the exposure
        length of the file at this path
        
    Returns
    -------
    median_dark : `~numpy.ndarray`
        Pixelwise median of all satisfactory darks
    """
    ref_header = fits.getheader(exp_same_as_image_at_path)
    ref_exposure_length = ref_header['EXPTIME']
    ref_temp = ref_header['SET-TEMP']

    last_header = fits.getheader(dark_images_paths[-1])
    rows, cols = last_header['NAXIS1'], last_header['NAXIS2']

    exposure_lengths = [fits.getheader(dark_path)['EXPTIME']
                        for dark_path in dark_images_paths]
    set_temps = [fits.getheader(dark_path).get('SET-TEMP', 100)
                 for dark_path in dark_images_paths]                            
                            
    darks_correct_exptime = [dark_path for dark_path, exp_len, set_temp in
                             zip(dark_images_paths, exposure_lengths, set_temps)
                             if exp_len == ref_exposure_length and set_temp == ref_temp]

    if len(darks_correct_exptime) == 0:
        raise ValueError("No darks with proper exposure length found. Needed "
                         "duration {0} s and temp {1} C, got durations {2} s " 
                         "and temp {3}".format(ref_exposure_length, ref_temp,
                                               sorted(set(exposure_lengths)), 
                                               sorted(set(set_temps))))

    darks_cube = np.zeros((rows, cols, len(darks_correct_exptime)))
    for i, dark_path in enumerate(darks_correct_exptime):
        darks_cube[:, :, i] = fits.getdata(dark_path)
    median_dark = np.median(darks_cube, axis=2)

    return median_dark

def master_flat(flat_images_paths, dark_images_paths=None, master_flat_dark=None):
    """
    Median-combine flats from list ``flat_images_paths`` after dark correcting
    images with the same exposure length.
        
    Parameters
    ----------
    flat_images_paths : list
        Paths to flat field images    
    
    dark_images_paths : list (optional)
        Paths to dark frames
        
    master_flat_dark : `~numpy.ndarray` (optional)
        Master dark frame of same exposure length as flats
        
    Returns
    -------
    flat : `~numpy.ndarray`
        Flat field (normalized)
    """
    # Get dark frame for flats
    if master_flat_dark is None:
        flat_dark = master_dark(dark_images_paths, flat_images_paths[-1])
    else:
        flat_dark = master_flat_dark

    first_header = fits.getheader(flat_images_paths[0])
    rows, cols = first_header['NAXIS1'], first_header['NAXIS2']

    flats_cube = np.zeros((rows, cols, len(flat_images_paths)))
    for i, flat_path in enumerate(flat_images_paths):
        flats_cube[:, :, i] = fits.getdata(flat_path)
    
    median_flat = np.median(flats_cube, axis=2) - flat_dark
    normalized_flat = median_flat/np.median(median_flat)
    return normalized_flat




    
