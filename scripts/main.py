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
    theta0=0,
    psip0=0.05,
    rho0=0.2,
    zeta0=0.0,
    mu=1e-6,
)

init = pm.InitialConditions(
    t0=0,
    theta0=0,
    psip0=0.02,
    rho0=0.05,
    zeta0=0,
    mu=1e-6,
)
particle = pm.Particle(init)
particle.run_ode(
    bfield=bfield,
    current=current,
    qfactor=qfactor,
    t_eval=(0.0, 200),
    steps=50000,
)
print(particle)

point_jump = 1
t = particle.t[::point_jump]
theta = particle.theta[::point_jump]
psip = particle.psip[::point_jump]
rho = particle.rho[::point_jump]
zeta = particle.zeta[::point_jump]
pzeta = particle.pzeta[::point_jump]
ptheta = particle.ptheta[::point_jump]

fig = plt.figure(**{"figsize": (9, 6), "layout": "constrained", "dpi": 120})
ax = fig.subplots(6, 1, sharex=True)
s, c = 0.9, "blue"
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
# ax.scatter(np.mod(theta, 2 * np.pi), ptheta, s, c)
ax.scatter(np.array(particle.t), np.array(zeta) / np.array(theta), s, c)
# ax.scatter(np.array(particle.theta), np.array(zeta), s, c)
# ax.set_title("θ-Pθ")

plt.show()
plt.close()
