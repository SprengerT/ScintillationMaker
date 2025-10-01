# ScintillationMaker
Quick tool to simulate scintillation dynamic spectra in the image approximation with user specified scintillation bandwidth and timescale.

## Theoretical background
We start from a random dispersive phase $\phi$ and neglect any frequency dependence of it. We also assume that the effective velocity is constant throughout the observation. Then, the electric field is given by

$E(t,\nu) \propto \sum_n \mu_n \exp\left( i \phi_n + i \Phi_n(t,\nu) \right) $

where $\mu$ is the amplitude of each image and $\Phi$ is its geometrical phase due to the length of its scattered path. Using canonical notations, the geometrical phase is given by

$\Phi = \pi\frac{\nu}{c}D_\text{eff}\left( \vec{\theta}_n - \frac{\vec{V}_\text{eff}}{D_\text{eff}}t \right)^2$

We now assume that the images' amplitudes are distributed like a Gaussian that is constant over frequency, which yields the often observed Lorentzian spectral ACF (although in reality there is a rather strong spectral evolution of the width of the scattering disk!):

$\mu_n \propto \exp\left( -\frac{\vec{\theta}_n^2}{2\theta_s^2} \right)$

$\theta_s$ is the characteristic scattering angle defining the size of the scattering disk.
For this model, the scintillation bandwidth $\nu_s$ and timescale $t_s$ are known:

$\nu_s = \frac{c}{\pi D_\text{eff} \theta_s^2}$

$t_s = \frac{c}{2\pi\nu_0\vert \vec{V}_\text{eff} \vert \theta_s}$

where $\nu_0$ is the central frequency.

In contrast to the effective distance, effective velocity, and the charcteristic scattering angle, these scales are directly accessible to the observer and fully define the observed scintillation in this model. Thus, we absorb physical constants into our variables:

$\tilde{\theta} = \theta/\theta_s$

$\tilde{\nu} = \nu/\nu_s$

$\tilde{t} = \tfrac{t}{t_{s}}$

We choose a coordinate system where the image position $\tilde{\theta}$ is separated into a part parallel to the effective velocity and a part perpendicular to it. Then we obtain the geometric phase in observational parameters:

$\Phi = \tilde{\nu}\left{\left( \tilde{\theta}_\parallel - \frac{\nu_s}{2\nu_0}\tilde{t} \right)^2 + \tilde{\theta}_\perp^2 \right}$

The minus sign here is not following directly from this transformation. However, we can define the coordinate systems of $\vec{\theta}$ and $\vec{V}_\text{eff}$ such that the sign is correct, which is the usual choice of coordinates.

The obtained formula can be used to simulate scintillation with desired scales without specifying any  physical distances or velocities.

## Practical implementation
TODO
