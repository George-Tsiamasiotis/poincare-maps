#!/usr/bin/env -S uv run --script --active
# /// script
# requires-python = ">=3.14"
# dependencies = [
#   "pyncare"
#   "matplotlib"
# ]
# ///
import matplotlib
import pyncare as pc

matplotlib.use("gtk3agg")

qfactor = pc.Qfactor("./data.nc", "akima")
current = pc.Current("./data.nc", "akima")
bfield = pc.Bfield("./data.nc", "bicubic")
per = pc.Perturbation(
    [
        pc.Harmonic("./data.nc", "akima", m=1, n=2, phase=0),
        pc.Harmonic("./data.nc", "akima", m=3, n=2, phase=0),
    ]
)

initial = pc.InitialConditions(
    t0=0,
    theta0=3.14,
    psip0=0.5 * qfactor.psip_wall,
    rho0=0.001,
    zeta0=0.0,
    mu=0,
)

particle = pc.Particle(initial)

particle.integrate(
    qfactor=qfactor,
    current=current,
    bfield=bfield,
    per=per,
    t_eval=[0, 20000],
)
print(particle)

pc.orbit_plot(particle)
