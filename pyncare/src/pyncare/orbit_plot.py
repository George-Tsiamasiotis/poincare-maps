from warnings import warn
import numpy as np
import matplotlib.pyplot as plt
from math import floor, log10
from pyncare import Particle

plt.rcParams["text.usetex"] = True

s = 0.3
c = "blue"
dpi = 120
figsize = (9, 6)
target_points = 50_000


def orbit_plot(particle: Particle, percentage: float = 100, downsample: bool = True):
    if percentage < 0 or percentage > 100:
        raise ValueError("Percentage must be between 0 and 100.")

    step = 1
    if downsample:
        length = len(particle.evolution.time)
        oom = floor(log10(length))
        target_oom = floor(log10(target_points))
        if oom > target_oom:
            step = 10 ** (oom - target_oom)

    points = int(np.floor(particle.evolution.time.shape[0] * percentage / 100) - 1)

    time = particle.evolution.time[:points][::step]
    theta = particle.evolution.theta[:points][::step]
    psip = particle.evolution.psip[:points][::step]
    rho = particle.evolution.rho[:points][::step]
    zeta = particle.evolution.zeta[:points][::step]
    pzeta = particle.evolution.pzeta[:points][::step]
    ptheta = particle.evolution.ptheta[:points][::step]
    energy = particle.evolution.energy[:points][::step]

    if downsample and len(time) > target_points * 10:
        warn("Downsampling did not work..")

    fig = plt.figure(figsize=figsize, layout="constrained", dpi=dpi)
    ax = fig.subplots(7, 1, sharex=True)
    ax[0].scatter(time, theta, s, c)
    ax[1].scatter(time, psip, s, c)
    ax[2].scatter(time, rho, s, c)
    ax[3].scatter(time, zeta, s, c)
    ax[4].scatter(time, ptheta, s, c)
    ax[5].scatter(time, pzeta, s, c)
    ax[6].scatter(time, energy, s, c)
    # Zoom out Pzeta plot
    if abs(np.std(pzeta)) < 1e-6:
        current_ylim = np.array(ax[5].get_ylim())
        ax[5].set_ylim(np.sort([current_ylim[0] / 3, current_ylim[1] * 3]))
    # Zoom out Energy plot
    if abs(np.std(energy)) < 1e-6:
        current_ylim = np.array(ax[6].get_ylim())
        ax[6].set_ylim(np.sort([current_ylim[0] / 2, current_ylim[1] * 2]))

    ax[0].set_xlabel(r"$\theta$")
    ax[1].set_xlabel(r"$\psi_p$")
    ax[2].set_xlabel(r"$\rho_{||}$")
    ax[3].set_xlabel(r"$\zeta$")
    ax[4].set_xlabel(r"$P_\theta$")
    ax[5].set_xlabel(r"$P_\zeta$")
    ax[6].set_xlabel(r"Energy")

    plt.show()
    plt.close()
