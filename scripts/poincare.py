#!/usr/bin/env -S uv run --script --active
# /// script
# requires-python = ">=3.14"
# dependencies = [
#   "matplotlib",
#   "numpy",
#   "pyncare"
# ]
# ///
import matplotlib.pyplot as plt
import pyncare as pc
import numpy as np

qfactor = pc.Qfactor("./data.nc", "akima")
currents = pc.Currents("./data.nc", "akima")
bfield = pc.Bfield("./data.nc", "bicubic")
perturbation = pc.Perturbation(
    [
        pc.Harmonic("./data.nc", "akima", m=0, n=1, phase=0),
        pc.Harmonic("./data.nc", "akima", m=1, n=9, phase=0),
        pc.Harmonic("./data.nc", "akima", m=1, n=7, phase=0),
    ]
)

num = 100
init = pc.PoincareInit(
    thetas=np.zeros(num),
    psips=np.linspace(0.0, qfactor.psip_wall, num),
    rhos=0.01 * np.ones(num),
    zetas=-0 * np.ones(num),
    mus=np.zeros(num),
)

params = pc.MappingParameters(section="theta", alpha=np.pi, intersections=1000)
poincare = pc.Poincare(init=init, params=params)
poincare.run(
    qfactor=qfactor,
    currents=currents,
    bfield=bfield,
    perturbation=perturbation,
)
print(poincare)

fig = plt.figure(figsize=(10, 5), layout="constrained", dpi=120)
ax = fig.add_subplot()
pc.poincare_plot(ax, poincare, wall=qfactor.psip_wall)
