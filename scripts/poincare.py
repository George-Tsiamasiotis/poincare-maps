import poincare_maps as pm
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("gtk3agg")

bfield = pm.Bfield("./data.nc", "bicubic")
qfactor = pm.Qfactor("./data.nc", "cubic")
current = pm.Current("./data.nc", "steffen")

# psip_wall = qfactor.psip_wall
# psi_wall = qfactor.psi_wall
psip_wall = 0.1
psi_wall = 0.1

points = 60
psips = np.linspace(0.01, 0.9 * psip_wall, points)
thetas = np.array([0] * points)
zetas = np.array([0] * points)

poincare = pm.Poincare()

for i in range(points):

    init = pm.InitialConditions(
        t0=0,
        theta0=thetas[i],
        psip0=psips[i],
        rho0=0.01,
        zeta0=0,
        mu=0,
    )

    particle = pm.Particle(init)
    poincare.add_particle(particle)

poincare.run(
    bfield=bfield,
    current=current,
    qfactor=qfactor,
    angle="zeta",
    intersection=np.pi / 2,
    turns=200,
)

angles = poincare.get_angles()
fluxes = poincare.get_fluxes()

# ==========================================================================


def pi_mod(arr: np.ndarray):
    a = np.mod(arr, 2 * np.pi)
    a = a - 2 * np.pi * (a > np.pi)
    return a


fig = plt.figure(**{"figsize": (11, 6), "layout": "constrained", "dpi": 150})
ax = fig.subplots()
s, c, marker = 0.5, "black", "."
for i in range(len(angles)):

    ax.scatter(pi_mod(angles[i]), fluxes[i], s, c, marker=marker)

ax.set_title("angle-flux")
ax.set_xlim(-np.pi, np.pi)
ax.set_ylim(0, psi_wall)
ax.set_xticks(
    np.linspace(-np.pi, np.pi, 5),
    ["-π", "-π/2", "0", "π/2", "π"],
)

plt.show()
plt.close()
