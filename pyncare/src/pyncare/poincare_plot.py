import pyncare as pc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

plt.rcParams["text.usetex"] = True
matplotlib.use("gtk3agg")


s = 0.05
c = "blue"
dpi = 150
figsize = (10, 7)


def poincare_plot(p: pc.Poincare):
    # TODO: add walls

    angles = p.angles
    fluxes = p.fluxes

    fig = plt.figure(figsize=figsize, layout="constrained", dpi=dpi)
    ax = fig.subplots()

    for i in range(len(angles)):
        ax.scatter(pi_mod(angles[i]), fluxes[i], s, c)

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

    plt.show()
    plt.close()


def pi_mod(arr: np.ndarray):
    a = np.mod(arr, 2 * np.pi)
    a = a - 2 * np.pi * (a > np.pi)
    return a
