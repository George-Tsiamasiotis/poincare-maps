#!/usr/bin/env -S uv run --script --active
# /// script
# requires-python = ">=3.14"
# dependencies = [
#   "pyncare"
# ]
# ///
import pyncare as pm

qfactor = pm.Qfactor("./data.nc", "akima")
current = pm.Current("./data.nc", "akima")
bfield = pm.Bfield("./data.nc", "bicubic")
per = pm.Perturbation(
    [
        pm.Harmonic("./data.nc", "akima", m=1, n=2, phase=0),
        pm.Harmonic("./data.nc", "akima", m=3, n=2, phase=0),
    ]
)

initial = pm.InitialConditions(
    t0=0,
    theta0=3.14,
    psip0=0.5 * qfactor.psip_wall,
    rho0=0.001,
    zeta0=0.0,
    mu=0,
)

particle = pm.Particle(initial)

particle.integrate(
    qfactor=qfactor, current=current, bfield=bfield, per=per, t_eval=[0, 100000]
)

print(particle)
