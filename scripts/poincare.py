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
current = pc.Current("./data.nc", "akima")
bfield = pc.Bfield("./data.nc", "bicubic")
per = pc.Perturbation(
    [
        pc.Harmonic("./data.nc", "akima", m=1, n=2, phase=0),
        pc.Harmonic("./data.nc", "akima", m=3, n=2, phase=0),
    ]
)

num = 20
init = pc.PoincareInit(
    thetas=np.linspace(0, np.pi, num),
    psips=np.linspace(0.0, 0.8, num) * qfactor.psip_wall,
    rhos=0.001 * np.ones(num),
    zetas=np.zeros(num),
    mus=np.zeros(num),
)

mapping = pc.Mapping(section="theta", alpha=3.14, intersections=200)
poincare = pc.Poincare(init=init, mapping=mapping)
poincare.run(
    qfactor=qfactor,
    current=current,
    bfield=bfield,
    per=per,
)

pc.poincare_plot(poincare)
