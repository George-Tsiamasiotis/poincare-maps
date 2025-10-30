#!/usr/bin/env -S uv run --script --active
# /// script
# requires-python = ">=3.14"
# dependencies = [
#     "matplotlib",
#     "pyncare",
# ]
# ///
import matplotlib
import matplotlib.pyplot as plt
from pyncare import Bfield, b_plot, db_plot


matplotlib.use("gtk3agg")

bfield = Bfield("./data.nc", "bicubic")

fig = plt.figure(**{"figsize": (15, 5), "layout": "constrained"})
fig.suptitle("Magnetic Field Profile")

ax = fig.subplots(1, 3)
b_plot(ax[0], bfield)
db_plot(ax[1], ax[2], bfield)

plt.show()
