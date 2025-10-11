import poincare_maps as pm
import numpy as np
import matplotlib.pyplot as plt

s = 0.9
c = "blue"
dpi = 120
figsize = (9, 6)


def plot_orbit(particle: pm.Particle, percentage: float = 100):
    if percentage < 0 or percentage > 100:
        raise ValueError("Percentage must be between 0 and 100.")

    points = int(np.floor(particle.t.shape[0] * percentage / 100) - 1)

    t = particle.t[:points]
    theta = particle.theta[:points]
    psip = particle.psip[:points]
    rho = particle.rho[:points]
    zeta = particle.zeta[:points]
    pzeta = particle.pzeta[:points]
    ptheta = particle.ptheta[:points]

    fig = plt.figure(**{"figsize": figsize, "layout": "constrained", "dpi": dpi})
    ax = fig.subplots(6, 1, sharex=True)
    ax[0].scatter(t, theta, s, c)
    ax[1].scatter(t, psip, s, c)
    ax[2].scatter(t, rho, s, c)
    ax[3].scatter(t, zeta, s, c)
    ax[4].scatter(t, ptheta, s, c)
    ax[5].scatter(t, pzeta, s, c)
    # Zoom out Pzeta plot
    current_ylim = np.array(ax[5].get_ylim())
    ax[5].set_ylim(np.sort([current_ylim[0] / 3, current_ylim[1] * 3]))

    ax[0].set_title(r"$\theta$")
    ax[1].set_title(r"$\psi_p$")
    ax[2].set_title(r"$\rho_{||}$")
    ax[3].set_title(r"$\zeta$")
    ax[4].set_title(r"$P_\theta$")
    ax[5].set_title(r"$P_\zeta$")

    plt.show()
    plt.close()
