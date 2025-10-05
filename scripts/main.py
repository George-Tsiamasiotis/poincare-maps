import poincare_maps as pm
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("gtk3agg")

bfield = pm.Bfield("./data.nc", "bicubic")
qfactor = pm.Qfactor("./data.nc", "akima")
current = pm.Current("./data.nc", "akima")

init = pm.InitialConditions(
    t0=0,
    theta0=0,
    psip0=0.05,
    rho0=0.004,
    zeta0=0,
    mu=1e-4,
    pzeta=-0.04,
)

particle = pm.Particle(init)
particle.run(
    bfield=bfield,
    current=current,
    qfactor=qfactor,
    step_size=1e-3,
    steps=20000,
)

points = np.asarray(particle.get_points()).T
time = points[0]
thetas = points[1]
psips = points[2]
rhos = points[3]
zetas = points[4]

fig = plt.figure(**{"figsize": (13, 7), "layout": "constrained"})
ax = fig.subplots(3, 2)
ax[0, 0].plot(time, thetas)
ax[0, 1].plot(time, psips)
ax[1, 0].plot(time, rhos)
ax[1, 1].plot(time, zetas)
ax[2, 0].plot(thetas, psips)

ax[0, 0].set_title("θ")
ax[0, 1].set_title("ψp")
ax[1, 0].set_title("ρ")
ax[1, 1].set_title("ζ")
ax[2, 0].set_title("θ-ψp")

plt.show()
