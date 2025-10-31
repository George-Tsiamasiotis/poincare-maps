import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from pyncare import Current

plt.rcParams["text.usetex"] = True


def g_plot(ax: Axes, current: Current):
    psip_data = current.psip_data
    g_data = current.g_data
    # Smooth derivative curve
    psips = np.linspace(current.psip_data[0], current.psip_data[-1], 1000)
    dg_dpsip = np.zeros((len(psips)))
    for i in range(len(dg_dpsip)):
        dg_dpsip[i] = current.dg_dpsip(psips[i])

    dax = ax.twinx()
    ax.scatter(psip_data, g_data, c="k", s=2, zorder=2, alpha=0.8, label="data points")
    ax.plot(psip_data, g_data, c="b", label=r"$g(\psi_p)$")
    dax.plot(psips, dg_dpsip, c="r", label=r"$\partial g(\psi_p)\partial \psi_p$")

    ax.set_xlabel(r"$\psi_p$")
    ax.set_ylabel(r"$g(\psi_p)$")
    dax.set_ylabel(r"$\partial g(\psi_p)/\partial \psi_p$")
    ax.grid(True)
    ax.margins(0)
    ax.legend()


def i_plot(ax: Axes, current: Current):
    psip_data = current.psip_data
    i_data = current.i_data
    psips = np.linspace(current.psip_data[0], current.psip_data[-1], 1000)
    di_dpsip = np.zeros((len(psips)))
    # Smooth derivative curve
    for i in range(len(di_dpsip)):
        di_dpsip[i] = current.di_dpsip(psips[i])

    dax = ax.twinx()
    ax.scatter(psip_data, i_data, c="k", s=2, zorder=2, alpha=0.8, label="data points")
    ax.plot(psip_data, i_data, c="b", label=r"$I(\psi_p)$")
    dax.plot(psips, di_dpsip, c="r", label=r"$\partial I(\psi_p)\partial \psi_p$")

    ax.set_ylabel(r"$I(\psi_p)$")
    dax.set_ylabel(r"$\partial I(\psi_p)/\partial \psi_p$")
    ax.set_xlabel(r"$\psi_p$")
    ax.grid(True)
    ax.margins(0)
    ax.legend()
