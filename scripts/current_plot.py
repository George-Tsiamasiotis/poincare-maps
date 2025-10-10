import poincare_maps as pm
import matplotlib
import matplotlib.pyplot as plt


from plots import g_plot, i_plot

matplotlib.use("gtk3agg")

current = pm.Current("./data.nc", "akima")

fig = plt.figure(**{"figsize": (15, 5), "layout": "constrained"})
fig.suptitle("Plasma Current")

ax = fig.subplots(1, 2)
g_plot(ax[0], current)
i_plot(ax[1], current)

plt.show()
