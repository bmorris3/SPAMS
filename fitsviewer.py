import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
import sys
import os
from glob import glob

# Open image from command line argument
#path = sys.argv[1]

path = sorted(glob('*.fits'), key=os.path.getmtime)[-1]
print("Opening: {0}".format(path))
img = fits.open(path)[0].data

def format_coord(x, y):
    '''
    Function to give data value on mouse over with imshow.
    '''
    try:
        return 'x={0:.1f}, y={1:.1f}, Flux={2:.0f}'.format(x, y, img[y,x])
    except:
        return 'x={0:d}, y={1:d}'.format(x, y)

# Get decent image scaling
m = np.median(img)
N = 4
s = np.std(img)

fig, ax = plt.subplots()
ax.imshow(img, interpolation='nearest', vmin=m-N*s/4, vmax=m+N*s,
          origin='lower', cmap=plt.cm.binary_r)
ax.format_coord = format_coord
plt.show()
