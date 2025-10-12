import poincare_maps as pm
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["text.usetex"] = True

levels = 30
cmap = "managua"


def b_plot(ax, bfield: pm.Bfield):

    theta_data = bfield.theta_data()
    psip_data = bfield.psip_data()
    r_data = bfield.r_data()
    z_data = bfield.z_data()
    b_data = bfield.b_data()
    _, psi_grid = np.meshgrid(theta_data, psip_data)

    contour = ax.contourf(r_data, z_data, b_data, **{"levels": levels, "cmap": cmap})
    ax.axis("equal")
    ax.set_title("Magnetic flux surfaces")
    ax.set_xlabel(r"$R[m]$")
    ax.set_ylabel(r"$Z[m]$")
    plt.colorbar(contour, ax=ax)


def db_plots(axx, axy, bfield: pm.Bfield):

    theta_data = bfield.theta_data()
    psip_data = bfield.psip_data()
    r_data = bfield.r_data()
    z_data = bfield.z_data()
    db_dpsip_grid = bfield.db_dpsip_data()
    db_dtheta_grid = bfield.db_dtheta_data()
    _, psi_grid = np.meshgrid(theta_data, psip_data)

    contour_kw = {"levels": levels, "cmap": cmap}

    contourx = axx.contourf(r_data, z_data, db_dpsip_grid, **contour_kw)
    contoury = axy.contourf(r_data, z_data, db_dtheta_grid, **contour_kw)
    plt.colorbar(contourx, ax=axx)
    plt.colorbar(contoury, ax=axy)
    axx.axis("equal")
    axy.axis("equal")
    axx.set_xlabel(r"$R[m]$")
    axx.set_ylabel(r"$Z[m]$")
    axy.set_xlabel(r"$R[m]$")
    axy.set_ylabel(r"$Z[m]$")
    axx.set_title(r"$\partial B/\partial\theta$")
    axy.set_title(r"$\partial B/\partial\psi_p$")


def q_plot(ax, qfactor: pm.Qfactor):
    psip_data = qfactor.psip_data()
    q_data = qfactor.q_data()

    ax.scatter(psip_data, q_data, c="k", s=2, zorder=2, alpha=0.8)
    ax.plot(psip_data, q_data, c="b")

    ax.set_xlabel(r"$\psi_p$")
    ax.set_ylabel(r"$q(\psi_p)$")
    ax.set_title(r"$q(\psi_p)$")
    ax.grid(True)
    ax.margins(0)


def psi_plot(ax, qfactor: pm.Qfactor):
    psip_data = qfactor.psip_data()
    psi_data = qfactor.psi_data()

    ax.scatter(psip_data, psi_data, c="k", s=2, zorder=2, alpha=0.8)
    ax.plot(psip_data, psi_data, c="b")

    ax.set_xlabel(r"$\psi_p$")
    ax.set_ylabel(r"$\psi(\psi_p)$")
    ax.set_title(r"$\psi(\psi_p)$")
    ax.grid(True)
    ax.margins(0)


def g_plot(ax, current: pm.Current):
    psip_data = current.psip_data()
    g_data = current.g_data()
    dg_dpsip_data = current.di_dpsip_data()
    dax = ax.twinx()

    ax.scatter(psip_data, g_data, c="k", s=2, zorder=2, alpha=0.8)
    ax.plot(psip_data, g_data, c="b")
    dax.plot(psip_data, dg_dpsip_data, c="r")

    ax.set_xlabel(r"$\psi_p$")
    ax.set_ylabel(r"$g(\psi_p)$")
    dax.set_ylabel(r"$\partial g(\psi_p)\partial \psi_p$")
    ax.set_title("Toroidal current")
    ax.grid(True)
    ax.margins(0)


def i_plot(ax, current: pm.Current):
    psip_data = current.psip_data()
    i_data = current.i_data()
    di_dpsip_data = current.di_dpsip_data()
    dax = ax.twinx()

    ax.scatter(psip_data, i_data, c="k", s=2, zorder=2, alpha=0.8)
    ax.plot(psip_data, i_data, c="b")
    dax.plot(psip_data, di_dpsip_data, c="r")

    ax.set_xlabel(r"$\psi_p$")
    ax.set_ylabel(r"$I(\psi_p)$")
    dax.set_ylabel(r"$\partial I(\psi_p)\partial \psi_p$")
    ax.set_title("Poloidal current")
    ax.grid(True)
    ax.margins(0)
