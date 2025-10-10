import poincare_maps as pm
import numpy as np
import matplotlib.pyplot as plt

levels = 30
cmap = "managua"


def b_plots(ax, bfield: pm.Bfield):

    theta_data = bfield.theta_data()
    psip_data = bfield.psip_data()
    rgrid, zgrid, bgrid = bfield.rzb_grids()
    _, psi_grid = np.meshgrid(theta_data, psip_data)

    contour = ax.contourf(rgrid, zgrid, bgrid, **{"levels": levels, "cmap": cmap})
    ax.axis("equal")
    plt.colorbar(contour, ax=ax)
