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
    exposure_length = fits.getheader(exp_same_as_image_at_path)['EXPTIME']

    first_header = fits.getheader(dark_images_paths[0])
    rows, cols = first_header['NAXIS1'], first_header['NAXIS2']
    
    darks_correct_exptime = [dark_path for dark_path in dark_images_paths
                             if fits.getheader(dark_path)['EXPTIME'] == exposure_length]
        
    darks_cube = np.zeros((rows, cols, len(darks_correct_exptime)))
    for i, dark_path in enumerate(darks_correct_exptime):
        darks_cube[:, :, i] = fits.getdata(dark_path)
    
    median_dark = np.median(darks_cube, axis=2)
    return median_dark

def master_flat(flat_images_paths, dark_images_paths):
    """
    Median-combine flats from list ``flat_images_paths`` after dark correcting
    images with the same exposure length.
        
    Parameters
    ----------
    flat_images_paths : list
        Paths to flat field images    
    
    dark_images_paths : list
        Paths to dark frames
        
    Returns
    -------
    flat : `~numpy.ndarray`
        Flat field (normalized)
    """
    # Get dark frame for flats
    flat_dark = master_dark(dark_images_paths, flat_images_paths[0])    

    first_header = fits.getheader(flat_images_paths[0])
    rows, cols = first_header['NAXIS1'], first_header['NAXIS2']

    flats_cube = np.zeros((rows, cols, len(flat_images_paths)))
    for i, flat_path in enumerate(flat_images_paths):
        flats_cube[:, :, i] = fits.getdata(flat_path)
    
    median_flat = np.median(flats_cube, axis=2) - flat_dark
    normalized_flat = median_flat/np.median(median_flat)
    return normalized_flat




    
