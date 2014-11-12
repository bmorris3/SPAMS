import numpy as np
from matplotlib import pyplot as plt

rawdata = np.loadtxt('WD2218+706lightcurve.txt')

times = rawdata[0,:]
lc = rawdata[1,:]
lc_err = rawdata[2,:]
R_wd = 1.3 #earth radii
R_p = np.linspace(0,1,100) #earth radius
R_moon = .272
R_mars = .53
R_mercury = .38
R_Kepler32f = .789
R_Ceres = .075

depth = (R_p/R_wd)**2
Nsigma = depth/np.std(lc)
#fig, ax = plt.subplots()
#ax.plot(times, lc, 'o')
#plt.show()

radii = [R_moon, R_mars, R_mercury, R_Kepler32f, R_Ceres]
radiinames = ['Moon', 'Mars', 'Mercury', 'Kepler-32f', 'Ceres']

fig, ax = plt.subplots()
#ax.axhline(5,ls='--', lw = 2, color = 'm')
for i, thing in enumerate(radii):
    ax.axvline(thing,ls='--', lw = 1, color = 'k')
    ax.annotate(radiinames[i], xy=(thing-.02, 15), textcoords = 'data', ha = 'center', rotation = 90, va = 'center')
ax.plot(R_p, Nsigma, lw = 2)
ax.fill_between(R_p, 0, 5, color = 'k', alpha = 0.5)
ax.set_xlabel('Planet Radius [$R_\oplus$]')
ax.set_ylabel('$\sigma$')
ax.set_title('ARCSAT Detection Confidence')
fig.savefig('plots/ARCSATDetectionConfidence.pdf', bbox_inches='tight')
plt.show()




