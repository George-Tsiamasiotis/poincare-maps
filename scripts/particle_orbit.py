import poincare_maps as pm
import matplotlib

from plots import plot_orbit

matplotlib.use("gtk3agg")

bfield = pm.Bfield("./data.nc", "bicubic")
qfactor = pm.Qfactor("./data.nc", "cubic")
current = pm.Current("./data.nc", "steffen")

init = pm.InitialConditions(
    t0=0,
    theta0=0,
    psip0=0.5 * qfactor.psip_wall,
    rho0=0.01,
    zeta0=0.0,
    mu=0,
)

particle = pm.Particle(init)
particle.run_ode(
    bfield=bfield,
    current=current,
    qfactor=qfactor,
    t_eval=(0.0, 2000),
    steps=0,
)
print(particle)

plot_orbit(particle)
