import poincare_maps as pm
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from plots import poincare_plot

matplotlib.use("gtk3agg")

bfield = pm.Bfield("./data.nc", "bicubic")
qfactor = pm.Qfactor("./data.nc", "akima")
current = pm.Current("./data.nc", "akima")
harmonics = [
    pm.Harmonic("./data.nc", "akima", 1, 7, 0),
    pm.Harmonic("./data.nc", "akima", 1, 9, 0),
]
per = pm.Perturbation(harmonics)

psip_wall = qfactor.psip_wall
psi_wall = qfactor.psi_wall

points = 60
psips = np.linspace(0, psip_wall, points)
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
    per=per,
    angle="theta",
    intersection=np.pi,
    turns=600,
)
print(poincare)

fig = plt.figure(**{"figsize": (10, 6), "layout": "constrained"})
ax = fig.subplots()

poincare_plot(ax, poincare, walls=(psip_wall, psi_wall))

plt.show()
plt.close()
