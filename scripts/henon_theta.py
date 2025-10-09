import poincare_maps as pm
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("gtk3agg")

bfield = pm.Bfield("./data.nc", "bicubic")
qfactor = pm.Qfactor("./data.nc", "cubic")
current = pm.Current("./data.nc", "steffen")

psip_wall = qfactor.psip_wall
psi_wall = qfactor.psi_wall

init = pm.InitialConditions(
    t0=0,
    theta0=0,
    psip0=0.04,
    rho0=0.05,
    zeta0=0,
    mu=1e-6,
)

particle = pm.Particle(init)
particle.run_henon(
    bfield=bfield,
    current=current,
    qfactor=qfactor,
    angle="theta",
    intersection=np.pi / 2,
    turns=100,
)
print(particle)

t = particle.t
theta = particle.theta
psip = particle.psip
rho = particle.rho
zeta = particle.zeta
pzeta = particle.pzeta
ptheta = particle.ptheta
psi = particle.psi

fig = plt.figure(**{"figsize": (9, 6), "layout": "constrained", "dpi": 120})
ax = fig.subplots(6, 1, sharex=True)
s, c = 0.9, "blue"
ax[0].scatter(t, theta, s, c)
ax[1].scatter(t, psip, s, c)
ax[2].scatter(t, rho, s, c)
ax[3].scatter(t, zeta, s, c)
ax[4].scatter(t, ptheta, s, c)
ax[5].scatter(t, pzeta, s, c)
ax[5].set_ylim([-0.05, 0.05])

ax[0].set_title("θ")
ax[1].set_title("ψp")
ax[2].set_title("ρ")
ax[3].set_title("ζ")
ax[4].set_title("Pθ")
ax[5].set_title("Pζ")

plt.show()
plt.close()


# ==========================================================================


def pi_mod(arr: np.ndarray):
    a = np.mod(arr, 2 * np.pi)
    a = a - 2 * np.pi * (a > np.pi)
    return a


fig = plt.figure(**{"figsize": (10, 7), "layout": "constrained", "dpi": 100})
ax = fig.subplots()
s, c, marker = 3, "black", "."
ax.scatter(pi_mod(zeta), psip, s, c, marker=marker)

ax.set_title("ζ-ψp")
ax.set_xlim(-np.pi, np.pi)
ax.set_ylim(0, psip_wall)
ax.set_xticks(
    np.linspace(-np.pi, np.pi, 5),
    ["-π", "-π/2", "0", "π/2", "π"],
)

plt.show()
plt.close()
