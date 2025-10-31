import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from pyncare import Bfield

plt.rcParams["text.usetex"] = True

levels = 30
cmap = "managua"


def b_plot(ax: Axes, bfield: Bfield):

    r_data = bfield.r_data
    z_data = bfield.z_data
    b_data = bfield.b_data

    contour = ax.contourf(r_data, z_data, b_data, **{"levels": levels, "cmap": cmap})
    ax.axis("equal")
    ax.set_title("Magnetic field strength")
    ax.set_xlabel(r"$R[m]$")
    ax.set_ylabel(r"$Z[m]$")
    plt.colorbar(contour, ax=ax)


def db_plot(axx: Axes, axy: Axes, bfield: Bfield):

    r_data = bfield.r_data
    z_data = bfield.z_data
    db_dpsip_grid = bfield.db_dpsip_data
    db_dtheta_grid = bfield.db_dtheta_data

    contour_kw = {"levels": levels, "cmap": cmap}

    contourx = axx.contourf(r_data, z_data, db_dpsip_grid, **contour_kw)
    contoury = axy.contourf(r_data, z_data, db_dtheta_grid, **contour_kw)
    axx.axis("equal")
    axy.axis("equal")
    axx.set_xlabel(r"$R[m]$")
    axx.set_ylabel(r"$Z[m]$")
    axy.set_xlabel(r"$R[m]$")
    axy.set_ylabel(r"$Z[m]$")
    axx.set_title(r"$\partial B/\partial\theta$")
    axy.set_title(r"$\partial B/\partial\psi_p$")
    plt.colorbar(contourx, ax=axx)
    plt.colorbar(contoury, ax=axy)
