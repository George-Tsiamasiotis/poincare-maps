import numpy as np
import matplotlib.pyplot as plt
from pyncare import Particle

plt.rcParams["text.usetex"] = True

s = 0.3
c = "blue"
dpi = 120
figsize = (9, 6)


def orbit_plot(particle: Particle, percentage: float = 100):
    if percentage < 0 or percentage > 100:
        raise ValueError("Percentage must be between 0 and 100.")

    points = int(np.floor(particle.evolution.time.shape[0] * percentage / 100) - 1)

    time = particle.evolution.time[:points]
    theta = particle.evolution.theta[:points]
    psip = particle.evolution.psip[:points]
    rho = particle.evolution.rho[:points]
    zeta = particle.evolution.zeta[:points]
    pzeta = particle.evolution.pzeta[:points]
    ptheta = particle.evolution.ptheta[:points]

    fig = plt.figure(figsize=figsize, layout="constrained", dpi=dpi)
    ax = fig.subplots(6, 1, sharex=True)
    ax[0].scatter(time, theta, s, c)
    ax[1].scatter(time, psip, s, c)
    ax[2].scatter(time, rho, s, c)
    ax[3].scatter(time, zeta, s, c)
    ax[4].scatter(time, ptheta, s, c)
    ax[5].scatter(time, pzeta, s, c)
    # Zoom out Pzeta plot
    if abs(np.nanmax(np.diff(pzeta))) < 1e-6:
        current_ylim = np.array(ax[5].get_ylim())
        ax[5].set_ylim(np.sort([current_ylim[0] / 3, current_ylim[1] * 3]))

    ax[0].set_xlabel(r"$\theta$")
    ax[1].set_xlabel(r"$\psi_p$")
    ax[2].set_xlabel(r"$\rho_{||}$")
    ax[3].set_xlabel(r"$\zeta$")
    ax[4].set_xlabel(r"$P_\theta$")
    ax[5].set_xlabel(r"$P_\zeta$")

    plt.show()
    plt.close()
