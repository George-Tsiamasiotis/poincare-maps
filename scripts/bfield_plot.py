import poincare_maps as pm
import matplotlib
import matplotlib.pyplot as plt


from plots import b_plot, db_plots

matplotlib.use("gtk3agg")

bfield = pm.Bfield("./data.nc", "bicubic")

fig = plt.figure(**{"figsize": (15, 5), "layout": "constrained"})
fig.suptitle("Magnetic Field Profile")

ax = fig.subplots(1, 3)
b_plot(ax[0], bfield)
db_plots(ax[1], ax[2], bfield)

plt.show()
