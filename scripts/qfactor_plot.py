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
from pyncare import Qfactor, q_plot, psi_plot


matplotlib.use("gtk3agg")

qfactor = Qfactor("./data.nc", "akima")

fig = plt.figure(figsize=(11, 5), layout="constrained")
fig.suptitle("q-factor Profile")

ax = fig.subplots(1, 2)
q_plot(ax[0], qfactor)
psi_plot(ax[1], qfactor)

plt.show()
