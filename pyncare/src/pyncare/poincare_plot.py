from matplotlib.axes import Axes
import pyncare as pc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

plt.rcParams["text.usetex"] = True
matplotlib.use("gtk3agg")


s = 0.3
c = "blue"
marker = "o"


def poincare_plot(ax: Axes, p: pc.Poincare, wall: float = np.nan):
    # TODO: add walls

    angles = p.angles
    fluxes = p.fluxes

    for i in range(len(angles)):
        ax.scatter(pi_mod(angles[i]), fluxes[i], s, c, marker=marker)

    if p.section == "ConstTheta":
        ax.set_xlabel(r"$\zeta$")
        ax.set_ylabel(r"$\psi_p$", rotation=0)
        ax.set_title(rf"$\zeta-\psi_p,$ cross section at $\theta={p.alpha:.4g}$")
    elif p.section == "ConstZeta":
        ax.set_xlabel(r"$\theta$")
        ax.set_ylabel(r"$\psi$", rotation=0)
        ax.set_title(rf"$\theta-\psi,$ cross section at $\zeta={p.alpha:.4g}$")

    ax.set_xlim(-np.pi, np.pi)
    ax.set_xticks(
        np.linspace(-np.pi, np.pi, 5),
        [r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"],
    )

    if not np.isnan(wall):
        ax.axhline(y=wall, c="r")

    ax.set_ylim(np.clip(ax.get_ylim(), a_min=0, a_max=None).tolist())

    plt.show()
    plt.close()


def pi_mod(arr: np.ndarray):
    a = np.mod(arr, 2 * np.pi)
    a = a - 2 * np.pi * (a > np.pi)
    return a
