# ScintillationMaker
Quick tool to simulate scintillation dynamic spectra in the image approximation with user specified scintillation bandwidth and timescale.

## Theoretical background
We start from a random dispersive phase $\phi$ and neglect any frequency dependence of it. We also assume that the effective velocity is constant throughout the observation. Then, the electric field is given by

${E(t,\nu) \propto \sum_n \mu_n \exp\left( i \phi_n + i \Phi_n(t,\nu) \right) }$

where $\mu$ is the amplitude of each image and $\Phi$ is its geometrical phase due to the length of its scattered path. Using canonical notations, the geometrical phase is given by

${\Phi_n = \pi\frac{\nu}{c}D_\text{eff}\left( \vec{\theta}\_n - \frac{\vec{V}\_\text{eff}}{D_\text{eff}}t \right)^2}$

We now assume that the images' amplitudes are distributed like a Gaussian that is constant over frequency, which yields the often observed Lorentzian spectral ACF (although in reality there is a rather strong spectral evolution of the width of the scattering disk!):

${\mu_n \propto \exp\left( -\frac{\vec{\theta}_n^2}{2\theta_s^2} \right)}$

$\theta_s$ is the characteristic scattering angle defining the size of the scattering disk.
For this model, the scintillation bandwidth $\nu_s$ and timescale $t_s$ are known:

${\nu_s = \frac{c}{\pi D_\text{eff} \theta_s^2}}$

${t_s = \frac{c}{2\pi\nu_0\vert \vec{V}_\text{eff} \vert \theta_s}}$

where $\nu_0$ is the central frequency.

In contrast to the effective distance, effective velocity, and the characteristic scattering angle, these scales are directly accessible to the observer and fully define the observed scintillation in this model. Thus, we absorb physical constants into our variables. Furthermore, we also absorb the central time and central frequency because the first is an arbitrary choice and the second is only important relative to the full bandwidth and the scintillation bandwidth. The redefined variables are:

${\tilde{\theta} = \theta/\theta_s}$

${\tilde{\nu} = \frac{\nu-\nu_0}{\nu_s}}$

${\tilde{t} = \frac{t-t_0}{t_{s}}}$

We choose a coordinate system where the image position $\tilde{\theta}$ is separated into a part parallel to the effective velocity -- pointing in positive direction of the velocity -- and a part perpendicular to it. Then we obtain the geometric phase in observational parameters:

$\Phi_n = \tilde{\nu}\times\left( {\tilde{\theta}}\_{\shortparallel,n}^2 + \tilde{\theta}_{\perp,n}^2 \right) - \tilde{t} {\tilde{\theta}}\_{\shortparallel,n} + {\cal O}\left(\frac{\nu_s}{\nu_0}\right) + \text{const.}$

We ignore terms containing ${\nu_s/\nu_0}$ because we are not interested in spectral evolution here. As a result, the input scintillation bandwidth is correct at all frequencies rather than in a narrow region around the input central frequency. Terms that do not contain $\tilde{\theta}$ only contribute to the phase which is lost when forming the intensity. Terms that do not contain $\tilde{t}$ and $\tilde{\nu}$ contribute but cannot be distuingished from the random phase $\phi_n$. Thus, they can too be ignored by absorbing them into this term, which leaves us with a very simple formula.

The obtained formula can be used to simulate scintillation with desired scales without specifying any physical distances or velocities.

## Practical implementation

Images are drawn from a uniform random distribution over $\tilde{\theta}\_{\shortparallel}$ and $\tilde{\theta}_{\perp}$ where both range over a circle of radius 3 around the origin by default. The number of images can be set (the default is 1000). The circular window on a theoretically infinite screen is an approximation done to save computation time. The range can be increased beyond the default of 3. All proportionality relations are replaced by equations. The overall magnification of the source is not known and the result needs a normalisation specified by the user.
