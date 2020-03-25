import pybert as pb
import pygimli as pg
import numpy as np
import matplotlib.pyplot as plt

mesh3d = pg.load('mesh.bms')

rho = []
for cell in mesh3d.cells():
    if cell.center().x() < 0.53 / 3:
        rho.append(50)
    elif cell.center().x() > (2 * 0.53 / 3):
        rho.append(150)
    else:
        rho.append(100)
rho = np.array(rho)

scheme  = pb.load('dd2.shm')

sim = pb.ERTManager()
sim_result = sim.simulate(mesh = mesh3d, res = rho, scheme = scheme)

sim_result.set('valid', np.ones(3904))
# invert
inv = pb.ERTManager()
inv.setMesh(mesh3d)
inv.setData(sim_result)
inv_result = inv.invert(err = 0.01, lam = 50)

# interpolate the 3dmesh to the 2dmesh for plotting
mesh2d = pg.load('2d.bms')
model2d = pg.interpolate(mesh3d, inv_result, mesh2d.cellCenters())
ax, c, gci = pg.show(mesh2d, model2d, logScale = False, showMesh = True, hold = True, cMap = 'jet', return_gci = True, colorBar = False, cMin = 25, cMax = 175)
cbar = plt.colorbar(gci, label=r"Resistivity ($\Omega$m)", orientation = 'horizontal')
plt.tight_layout()
plt.savefig('model.png')
plt.show()
