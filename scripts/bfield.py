import poincare_maps as pm
import matplotlib
import matplotlib.pyplot as plt


from plots import b_plots

matplotlib.use("gtk3agg")

bfield = pm.Bfield("./data.nc", "bicubic")

fig = plt.figure(**{"figsize": (6, 6), "layout": "constrained"})
fig.suptitle("Magnetic Field Profile")

ax = fig.subplots(1, 1)
b_plots(ax, bfield)


plt.show()
