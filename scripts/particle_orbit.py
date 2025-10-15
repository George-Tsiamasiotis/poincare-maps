import poincare_maps as pm
import matplotlib

from plots import plot_orbit

matplotlib.use("gtk3agg")

bfield = pm.Bfield("./data.nc", "bicubic")
qfactor = pm.Qfactor("./data.nc", "cubic")
current = pm.Current("./data.nc", "steffen")
harmonics = [
    pm.Harmonic("./data.nc", "akima", 1, 7, 0),
    pm.Harmonic("./data.nc", "akima", 1, 8, 0),
]
per = pm.Perturbation(harmonics)

init = pm.InitialConditions(
    t0=0,
    theta0=0,
    psip0=0.025,
    rho0=0.01,
    zeta0=0.0,
    mu=0,
)

particle = pm.Particle(init)
particle.run_ode(
    bfield=bfield,
    current=current,
    qfactor=qfactor,
    per=per,
    t_eval=(0.0, 15000),
    steps=0,
)
print(particle)

plot_orbit(particle)
