"""
convert -delay 5 -loop 0 *.png field.gif

"""


import numpy as np
from matplotlib import pyplot as plt
import pyfits
from glob import glob

path="/astro/store/scratch/tmp/bmmorris/APO/UT150108/Q1UW02/UT150117/"

imagepath=glob(path+"sdss*")
for i in range(len(imagepath)):
    
    image=pyfits.getdata(imagepath[i])[:,:256]
    header=pyfits.getheader(imagepath[i])
    n=1.25
    m=np.mean(image)
    s=np.std(image)
    import matplotlib.cm as cm
    flat=pyfits.getdata(path+'quickanalysis/flat.fits')[:,:256]
    
    plt.imshow(image/flat,vmin=m-n*s,vmax=m+n*s,cmap=cm.bone,origin='lower',
               interpolation='nearest')
    plt.title(header["DATE-OBS"].split("T")[1])
    plt.yticks([])
    plt.xticks([])
    plt.savefig("plots/"+str(i).zfill(4)+".png",bbox_inches="tight")
    plt.close()
    
    
   # plt.show()
