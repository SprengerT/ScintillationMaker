 
import numpy as np
from numpy import newaxis as na
import matplotlib.style as mplstyle
mplstyle.use('fast')
import matplotlib.pyplot as plt
from scipy.signal import correlate,correlation_lags 

from lib_ScintillationMaker import simulate_scintillation

mHz = 1.0e-3
kHz = 1.0e+3
MHz = 1.0e+6
GHz = 1.0e+9
minute = 60.
hour = 3600.
day = 24.*hour

#user input
N_t = 2000
N_nu = 2000
bandwidth = 250.*MHz
nu0 = 1360.*MHz
duration = 120.*minute
nu_s = 1.*MHz
t_s = 1.*minute

#derived quantities
nu = np.linspace(0.,bandwidth,num=N_nu)
nu = nu-np.mean(nu)+nu0
t = np.linspace(0.,duration,num=N_t)
t = t-np.mean(t)

# The higher N_im, the closer the output gets to the ensemble average, but at computational cost
E = simulate_scintillation(t,nu,t_s,nu_s,N_im=1000)

I = np.abs(E)**2
I = I/np.mean(I)

##################
#### The following are diagnostics and can be removed ####

figure = plt.figure(figsize=(8, 4), dpi=300)
ax = figure.add_subplot(1,1,1)
plot = ax.pcolormesh(t,nu,np.swapaxes(I,0,1))
ax.set_xlabel(r"time $t$ [s]")
ax.set_ylabel(r"frequency $\nu$ [Hz]")
figure.colorbar(plot, ax=ax)
plt.show()

ACF_out = correlate(I-np.mean(I),I-np.mean(I),mode='same')/correlate(np.ones_like(I),np.ones_like(I),mode='same')
dt = np.mean(np.diff(t))
dnu = np.mean(np.diff(nu))
Dt = dt*correlation_lags(N_t,N_t,mode='same')
Dnu = dnu*correlation_lags(N_nu,N_nu,mode='same')

# Theoretical expectation in ensemble average
def GaussianACF2D(Dt,Dnu,nu_s,t_s,modindex,offset):
    model = np.empty((len(Dt),len(Dnu)),dtype=float)
    model[:,:] = (modindex**2 /( 1 + Dnu[na,:]**2/nu_s**2 ) * np.exp( -Dt[:,na]**2/(2*t_s**2*(1 + Dnu[na,:]**2/nu_s**2)) ) + 1)*(1+offset) - 1
    return model
ACF_in = GaussianACF2D(Dt,Dnu,nu_s,t_s,1,0)

figure = plt.figure(figsize=(8, 4), dpi=300)
ax = figure.add_subplot(1,2,1)
plot = ax.pcolormesh(Dt,Dnu,np.swapaxes(ACF_in,0,1))
ax.set_ylabel(r"$\delta\nu$ [Hz]")
ax.set_xlabel(r"$\delta t$ [s]")
ax.set_title("input ACF (analytical)")
figure.colorbar(plot, ax=ax)
ax = figure.add_subplot(1,2,2)
plot = ax.pcolormesh(Dt,Dnu,np.swapaxes(ACF_out,0,1))
ax.set_xlabel(r"$\delta t$ [s]")
ax.set_title("output ACF (numerical)")
figure.colorbar(plot, ax=ax)
plt.show()

