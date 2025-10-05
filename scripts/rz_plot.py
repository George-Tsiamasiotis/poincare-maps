"""Magnetic Filed Profile Plot."""

import poincare_maps as pm
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("gtk3agg")

bfield = pm.Bfield("./data.nc", "bicubic")

theta_data = bfield.theta_data()
psip_data = bfield.psip_data()
rgrid, zgrid, bgrid = bfield.rzb_grids()
_, psi_grid = np.meshgrid(theta_data, psip_data)

fig = plt.figure(**{"figsize": (6, 6), "layout": "constrained"})
fig.suptitle("Magnetic Field Profile")

ax = fig.subplots(1, 1)

contour = ax.contourf(rgrid, zgrid, bgrid, **{"levels": 30, "cmap": "managua"})
ax.axis("equal")

plt.show()
