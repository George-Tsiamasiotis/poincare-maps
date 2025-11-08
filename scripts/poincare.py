#!/usr/bin/env -S uv run --script --active
# /// script
# requires-python = ">=3.14"
# dependencies = [
#   "matplotlib",
#   "numpy",
#   "pyncare"
# ]
# ///
import pyncare as pc
import numpy as np

qfactor = pc.Qfactor("./data.nc", "akima")
currents = pc.Currents("./data.nc", "akima")
bfield = pc.Bfield("./data.nc", "bicubic")
perturbation = pc.Perturbation(
    [
        pc.Harmonic("./data.nc", "akima", m=1, n=10, phase=0),
        pc.Harmonic("./data.nc", "akima", m=1, n=8, phase=0),
    ]
)

num = 100
init = pc.PoincareInit(
    thetas=np.linspace(0, np.pi, num),
    psips=np.linspace(0.0, qfactor.psip_wall, num),
    rhos=0.01 * np.ones(num),
    zetas=np.zeros(num),
    mus=np.zeros(num),
)

params = pc.MappingParameters(section="theta", alpha=np.pi, intersections=300)
poincare = pc.Poincare(init=init, params=params)
poincare.run(
    qfactor=qfactor,
    currents=currents,
    bfield=bfield,
    perturbation=perturbation,
)

pc.poincare_plot(poincare, wall=qfactor.psip_wall)
