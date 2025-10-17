import poincare_maps as pm
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["text.usetex"] = True


def poincare_plot(ax, p: pm.Poincare, walls=None):

    angles = p.get_angles()
    fluxes = p.get_fluxes()

    s, c = 0.05, "blue"

    for i in range(len(angles)):

        ax.scatter(pi_mod(angles[i][1:]), fluxes[i][1:], s, c)
        # angle0, flux0 = angles[i, 0], fluxes[i, 1]
        # ax.scatter(pi_mod(angle0), flux0, marker="_", s=10, color="k")

    if p.angle == "theta":
        ax.set_xlabel(r"$\zeta$")
        ax.set_ylabel(r"$\psi_p$", rotation=0)
        ax.set_title(rf"$\zeta-\psi_p,$ cross section at $\theta={p.intersection:.4g}$")
        wall = walls[0]
    elif p.angle == "zeta":
        ax.set_xlabel(r"$\theta$")
        ax.set_ylabel(r"$\psi$", rotation=0)
        ax.set_title(rf"$\theta-\psi,$ cross section at $\zeta={p.intersection:.4g}$")
        wall = walls[1]

    ylim = ax.get_ylim()
    ax.axhline(y=wall, color="red")

    ax.set_xlim(-np.pi, np.pi)
    ax.set_ylim(ylim)
    ax.set_xticks(
        np.linspace(-np.pi, np.pi, 5),
        [r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"],
    )


def pi_mod(arr: np.ndarray):
    a = np.mod(arr, 2 * np.pi)
    a = a - 2 * np.pi * (a > np.pi)
    return a
