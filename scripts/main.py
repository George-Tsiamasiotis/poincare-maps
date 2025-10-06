import poincare_maps as pm
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("gtk3agg")

bfield = pm.Bfield("./data.nc", "bicubic")
qfactor = pm.Qfactor("./data.nc", "cubic")
current = pm.Current("./data.nc", "steffen")

init = pm.InitialConditions(
    t0=0,
    theta0=3.14,
    psip0=0.05,
    rho0=0.05,
    zeta0=0.1,
    mu=1e-4,
)

particle = pm.Particle(init)
particle.run(
    bfield=bfield,
    current=current,
    qfactor=qfactor,
    step_size=1e-1,
    steps=5000,
)
print(particle)

t = particle.t
theta = particle.theta
psip = particle.psip
rho = particle.rho
zeta = particle.zeta
pzeta = particle.pzeta
ptheta = particle.ptheta

fig = plt.figure(**{"figsize": (9, 6), "layout": "constrained", "dpi": 120})
ax = fig.subplots(6, 1, sharex=True)
s, c = 0.08, "blue"
ax[0].scatter(t, theta, s, c)
ax[1].scatter(t, psip, s, c)
ax[2].scatter(t, rho, s, c)
ax[3].scatter(t, zeta, s, c)
ax[4].scatter(t, ptheta, s, c)
ax[5].scatter(t, pzeta, s, c)
ax[5].set_ylim([-0.05, 0.0])

ax[0].set_title("θ")
ax[1].set_title("ψp")
ax[2].set_title("ρ")
ax[3].set_title("ζ")
ax[4].set_title("Pθ")
ax[5].set_title("Pζ")

plt.show()
plt.close()

fig = plt.figure(**{"figsize": (6, 4), "layout": "constrained"})
ax = fig.subplots()
ax.scatter(np.mod(theta, 2 * np.pi), ptheta, s, c)
ax.set_title("θ-Pθ")

plt.show()
plt.close()
