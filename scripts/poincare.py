import poincare_maps as pm
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from plots import poincare_plot

matplotlib.use("gtk3agg")

bfield = pm.Bfield("./data.nc", "bicubic")
qfactor = pm.Qfactor("./data.nc", "akima")
current = pm.Current("./data.nc", "akima")

psip_wall = qfactor.psip_wall
psi_wall = qfactor.psi_wall
# For `lar` feature, since q=1
# psip_wall = 0.1
# psi_wall = 0.1

points = 60
psips = np.linspace(0.01 * psip_wall, 0.99 * psip_wall, points)
thetas = np.array([0] * points)
zetas = np.array([0] * points)

poincare = pm.Poincare()

for i in range(points):

    init = pm.InitialConditions(
        t0=0,
        theta0=thetas[i],
        psip0=psips[i],
        rho0=0.01,
        zeta0=zetas[i],
        mu=0,
    )

    particle = pm.Particle(init)
    poincare.add_particle(particle)

poincare.run(
    bfield=bfield,
    current=current,
    qfactor=qfactor,
    angle="theta",
    intersection=np.pi,
    turns=200,
)

fig = plt.figure(**{"figsize": (10, 6), "layout": "constrained"})
ax = fig.subplots()

poincare_plot(ax, poincare, walls=(psip_wall, psi_wall))

plt.show()
plt.close()
