import poincare_maps as pm
import numpy as np
import matplotlib.pyplot as plt

levels = 30
cmap = "managua"


def b_plot(ax, bfield: pm.Bfield):

    theta_data = bfield.theta_data()
    psip_data = bfield.psip_data()
    rgrid, zgrid = bfield.rz_grid()
    bgrid = bfield.b_grid()
    _, psi_grid = np.meshgrid(theta_data, psip_data)

    contour = ax.contourf(rgrid, zgrid, bgrid, **{"levels": levels, "cmap": cmap})
    ax.axis("equal")
    plt.colorbar(contour, ax=ax)


def db_plots(axx, axy, bfield: pm.Bfield):

    theta_data = bfield.theta_data()
    psip_data = bfield.psip_data()
    rgrid, zgrid = bfield.rz_grid()
    db_dpsip_grid, db_dtheta_grid = bfield.db_grids()
    _, psi_grid = np.meshgrid(theta_data, psip_data)

    contour_kw = {"levels": levels, "cmap": cmap}

    contourx = axx.contourf(rgrid, zgrid, db_dpsip_grid, **contour_kw)
    contoury = axy.contourf(rgrid, zgrid, db_dtheta_grid, **contour_kw)
    axx.axis("equal")
    axy.axis("equal")
    plt.colorbar(contourx, ax=axx)
    plt.colorbar(contoury, ax=axy)
