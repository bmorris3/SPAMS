from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np
import matplotlib.pyplot as plt
import h5py
from astropy.io import fits
from glob import glob

paths = sorted(glob('/Users/bmmorris/data/20160509/sdss1519_sdss_g_*.fits'))
first_fits_header = fits.getheader(paths[0])
first_fits_data = fits.getdata(paths[0])[:first_fits_header['NAXIS1'],
                                         :first_fits_header['NAXIS2']]

f = h5py.File('20160509.hdf5', 'w')
#f = h5py.File('20160509.hdf5', 'r+')

# Create image cube
data = f.create_dataset('images', dtype=np.float,
                        shape=(len(paths), first_fits_data.shape[0],
                               first_fits_data.shape[1]))
# Assemble dictionary of metadata
keys = set([i for i in first_fits_header.keys()
            if len(i) > 0 and i != 'HISTORY'])
meta = {key: [] for key in keys}

# Collect data, metadata from each FITS file
for i in range(len(paths)):
    f = fits.open(paths[i])
    data[i, ...] = f[0].data
    header = f[0].header

    for key in keys:
        meta[key].append(header[key])

# Store metadata in HDF5 attributes
for key in keys:
    print(key)
    if type(meta[key][0]) == str:
        dtype = "S"
    else:
        dtype = None
    data.attrs.create(key, np.array(meta[key], dtype=dtype))

f.close()