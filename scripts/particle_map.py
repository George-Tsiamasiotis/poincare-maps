#!/usr/bin/env -S uv run --script --active
# /// script
# requires-python = ">=3.14"
# dependencies = [
#   "matplotlib",
#   "pyncare"
# ]
# ///
import matplotlib
import pyncare as pc

matplotlib.use("gtk3agg")

qfactor = pc.Qfactor("./data.nc", "akima")
currents = pc.Currents("./data.nc", "akima")
bfield = pc.Bfield("./data.nc", "bicubic")
perturbation = pc.Perturbation(
    [
        pc.Harmonic("./data.nc", "akima", m=1, n=7),
        pc.Harmonic("./data.nc", "akima", m=1, n=9),
    ]
)

initial = pc.InitialConditions(
    time0=0,
    theta0=3.14,
    psip0=0.5 * qfactor.psip_wall,
    rho0=0.001,
    zeta0=0.0,
    mu=0,
)

params = pc.MappingParameters(section="theta", alpha=3.14, intersections=100)

particle = pc.Particle(initial)

particle.map(
    qfactor=qfactor,
    currents=currents,
    bfield=bfield,
    perturbation=perturbation,
    params=params,
)
print(particle)

pc.orbit_plot(particle)
