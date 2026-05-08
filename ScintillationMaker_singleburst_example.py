import numpy as np
from numpy import newaxis as na
import matplotlib.style as mplstyle
mplstyle.use('fast')
import matplotlib.pyplot as plt
from scipy.signal import correlate #,correlation_lags 

from lib_ScintillationMaker import simulate_singleburst_spectrum

mHz = 1.0e-3
kHz = 1.0e+3
MHz = 1.0e+6
GHz = 1.0e+9
minute = 60.
hour = 3600.
day = 24.*hour

#user input
N_nu = 1000000
bandwidth = 1.
nu_s = bandwidth/N_nu*10 #10 data points per scintillation bandwidth

#derived quantities
nu = np.linspace(0.,bandwidth,num=N_nu)

# The higher N_im, the closer the output gets to the ensemble average, but at computational cost
E = simulate_singleburst_spectrum(nu,nu_s,N_im=10000,th_lim=4)

I = np.abs(E)**2
I = I/np.mean(I) 

ACF_out = correlate(I-np.mean(I),I-np.mean(I),mode='same')/correlate(np.ones_like(I),np.ones_like(I),mode='same')
dnu = np.mean(np.diff(nu))

# Theoretical expectation in ensemble average
def GaussianACF1D(Dnu,nu_s,modindex,offset):
    model = np.empty(len(Dnu),dtype=float)
    model = (modindex**2 /( 1 + Dnu**2/nu_s**2 ) + 1)*(1+offset) - 1
    return model
Dnu = np.linspace(-bandwidth/2.,bandwidth/2.,num=len(nu))
ACF_in = GaussianACF1D(Dnu,nu_s,1,0)

figure = plt.figure(figsize=(8, 4), dpi=200)
ax = figure.add_subplot(1,1,1)
N_sets = 2
for i in range(N_sets):
    N_nu_set = (N_nu//N_sets)*(i+1)
    data = I[:N_nu_set]
    ACF_out = correlate(data-np.mean(data),data-np.mean(data),mode='same')/correlate(np.ones_like(data),np.ones_like(data),mode='same')
    Dnu_out = nu[:N_nu_set]-nu[N_nu_set//2]
    N_scint = (nu[N_nu_set-1]-nu[0])/nu_s
    ax.plot(Dnu_out/nu_s,ACF_out,label=r"{0:.0f} BW/$\nu_s$".format(N_scint))

ax.plot(Dnu/nu_s,ACF_in,label="theory",color="black")
ax.legend()
ax.set_xlim([-8,8])
ax.set_xlabel(r"$\delta\nu/nu_s$ ")
plt.show()
