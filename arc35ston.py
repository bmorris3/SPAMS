import numpy as np
from matplotlib import pyplot as plt

rawdatawd = np.loadtxt('WD2218+706lightcurve.txt')
rawdataGD = np.loadtxt('GD165lightcurve.txt')
timeswd = rawdatawd[0,:]
lcwd = rawdatawd[1,:]
lc_errwd = rawdatawd[2,:]
timesGD = rawdataGD[0,:]
lcGD = rawdataGD[1,:]
lc_errGD = rawdataGD[2,:]


R_wd = 1.3 #earth radii
R_p = np.linspace(0,1,100) #earth radius
R_moon = .272
R_mars = .53
R_mercury = .38
R_Kepler32f = .789
R_Ceres = .075

depth = (R_p/R_wd)**2
Nsigmawd = depth/np.std(lcwd)
#fig, ax = plt.subplots()
NsigmaGD = depth/np.std(lcGD)
#ax.plot(times, lc, 'o')
#plt.show()

radii = [R_moon, R_mars, R_mercury, R_Kepler32f, R_Ceres]
radiinames = ['Moon', 'Mars', 'Mercury', 'Kepler-32f', 'Ceres']

fig, ax = plt.subplots()
#ax.axhline(5,ls='--', lw = 2, color = 'm')
for i, thing in enumerate(radii):
    ax.axvline(thing,ls='--', lw = 1, color = 'k')
    ax.annotate(radiinames[i], xy=(thing-.02, 70), textcoords = 'data', ha = 'center', rotation = 90, va = 'center')
ax.fill_between(R_p, 0, 5, color = 'k', alpha = 0.7)
ax.plot(R_p, Nsigmawd, lw = 2, label='ARCSAT',color='#3385FF')
ax.plot(R_p, NsigmaGD, lw = 2, label='3.5m',color='r')
ax.annotate('Marginal Detection', xy=(1, 2.5), textcoords = 'data', ha = 'right', rotation = 0, va = 'center',color='w')
ax.legend(loc='center left')
ax.set_xlabel('Planet Radius [$R_\oplus$]')
ax.set_ylabel('$\sigma$')
ax.set_title('ARCSAT Detection Confidence')
fig.savefig('plots/ARC35DetectionConfidence.pdf', bbox_inches='tight')
plt.show()




