# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 13:31:39 2020

@author: Benjamin
"""

import os
import numpy as np
import pygimli as pg
import pygimli.meshtools as mt
from pygimli.physics import ert
from pygimli.physics.ert import simulate as simulateERT
pyvista = pg.optImport("pyvista")

import matplotlib.pyplot as plt


def _check_area(area):
    """
    Check that the area argument is valid.
    For example, the west limit should not be greater than the east limit.
    """
    x1, x2, y1, y2 = area
    assert x1 <= x2, \
        "Invalid area dimensions {}, {}. x1 must be < x2.".format(x1, x2)
    assert y1 <= y2, \
        "Invalid area dimensions {}, {}. y1 must be < y2.".format(y1, y2)

def regular(area, shape, z=None):
    nx, ny = shape
    x1, x2, y1, y2 = area
    _check_area(area)
    xs = np.linspace(x1, x2, nx)
    ys = np.linspace(y1, y2, ny)
    # Must pass ys, xs in this order because meshgrid uses the first argument
    # for the columns
    arrays = np.meshgrid(ys, xs)[::-1]
    if z is not None:
        arrays.append(z*np.ones(nx*ny, dtype=np.float))
    return [i.ravel() for i in arrays]


#%% create model geometry

depth = 60  # mesh depth
width = 350 # mesh width

widthlandfill = 50 # 5, 15, 25, 35
Hlandfill = 10 # 5, 15, 25, 35

# -------------------
tickness_liner = 1

HAno = 1.5 
depthAno = -5.0 # 3.75, 13.75 or 23.75
widthAno = tickness_liner/2 # 5, 15, 25, 35
ticknessAno = tickness_liner/2
shift = (widthlandfill + tickness_liner)/2 - ticknessAno/2
#shift = 0

# world = mt.createWorld(start=[-width/2, -width/2, 0], end=[width/2, width/2, -depth],
#                         worldMarker=True) # ,area=10 markerPosition=[width,width]


def create_plc():

    world = mt.createCube(size=[width, width, depth], 
                          pos=[width/2, width/2, -depth/2],
                          marker=1,
                          area=75,
                          markerPosition=[width/5, width/5, -depth/2]
                          )

    for i, b in enumerate(world.boundaries()):
        if b.norm()[2] == 1.0:
            b.setMarker(pg.core.MARKER_BOUND_HOMOGEN_NEUMANN)
        else:
            b.setMarker(pg.core.MARKER_BOUND_MIXED)


    liner = mt.createCube(size=[widthlandfill + tickness_liner, widthlandfill + tickness_liner, Hlandfill], 
            pos=[width/2, width/2, -Hlandfill/2], 
            marker=2, 
            area=1
            )
              
    landfill = mt.createCube(size=[widthlandfill, widthlandfill, Hlandfill], 
                        pos=[width/2, width/2, -Hlandfill/2], 
                        marker=3, 
                        )


    ano = mt.createCube(size=[widthAno, ticknessAno, HAno], 
                        pos=[width/2+shift, width/2, depthAno], 
                        marker=4, 
                        area=1e-3
                        )

    # plc = mt.mergePLC([world, liner, landfill, ano])
    plc = mt.mergePLC([world,liner,landfill, ano])

    return world,landfill,ano,plc

#%%

def create_scheme(world,landfill,ano,plc):

    n_sensors_surf = 16*16
    sensors_surf = np.zeros((n_sensors_surf, 3))

    shiftx = 51
    shifty = 0

    area = [landfill.xmin()+shiftx,
            landfill.xmax()+shiftx,
            landfill.ymin()+shifty,
            landfill.ymax()++shifty,
            ]

    shape=(16,16)
    xy= np.vstack(regular(area,shape))
    id_sensors = np.arange(0,len(xy[0,:]),2)
    sensors_surf = np.zeros([len(id_sensors),3])
    sensors_surf[:,0:2]=xy[0:2,id_sensors].T

    np.shape(sensors_surf)
    # add MALM electrodes 
    EAp= [ano.xmax()/2+ano.xmin()/2,
          ano.ymax()/2+ano.ymin()/2,
          depthAno]
    EBp = [world.xmax()/1.1, world.ymax()/2, 0]
    ENp = [world.xmin()+10, world.ymax()/2, 0]
        

    elecs= np.vstack([sensors_surf,EAp,EBp,ENp])

    np.savetxt('coords_elecs.txt', elecs)

    plt.figure()

    plt.axes()
    rectangle = plt.Rectangle((width/2-widthlandfill/2,width/2-widthlandfill/2), widthlandfill, widthlandfill,ec="red")
    plt.gca().add_patch(rectangle)

    plt.plot(ano.xmin(), ano.ymin(), '*r',label='Ano position', markersize=12)


    plt.plot(elecs[:,0], elecs[:,1], '.g',label='Elec positions', markersize=2)
    plt.plot(elecs[-3,0], elecs[-3,1], '+k',label='Elec A (Masse)', markersize=12)
    plt.plot(elecs[-2,0], elecs[-2,1], '.y',label='Elec B remote', markersize=12)
    plt.plot(elecs[-1,0], elecs[-1,1], '.b',label='Elec N remote', markersize=12)

    plt.axis('scaled')
    plt.xlim([world.xmin(),world.xmax()])
    plt.ylim([world.xmin(),world.xmax()])

    plt.savefig('./elecs.png', dpi=450)


    malmScheme = pg.DataContainerERT()
    for ei in elecs:
        malmScheme.createSensor([ei[0], ei[1], ei[2]])

    EA = len(elecs) -3
    EB = len(elecs) -2
    EN = len(elecs) -1
    meas = []
    for j in np.arange(0,len(elecs)-3):
        M=j
        meas.append([EA, EB,M, EN])

    meas = np.array(meas)
    malmScheme.resize(len(meas))    

    for i, j in enumerate("abmn"):
        malmScheme.set(j, meas[:, i])
        
    malmScheme.save('./malmScheme.shm')
    print(malmScheme)

    # %% put the electrodes as node into the mesh
    print(len(malmScheme.sensorPositions()))
    for i, pos in enumerate(malmScheme.sensorPositions()):
        # if i >len(malmScheme.sensorPositions())-150: 
            # print(pos)
            # print(i)
            plc.createNode(pos, marker=-99)  # electrode
            plc.createNode(pos + [0, 0, -0.01])  # refinement

    mesh = mt.createMesh(plc)
    mesh.exportVTK('./fwd_mesh_landfill.vtk')

    return malmScheme, mesh


def run_fwd_model(malmScheme, mesh):
    #%%
    SoilR = 100
    linerR = 1e6
    AnoR = SoilR
    noiseLevel = 5 #0.01 #1 
    noiseAbs= 1e-6 #1e-6 # 1e-6

    u_elecs = simulateERT(mesh, scheme=malmScheme, res= [[1,SoilR],[2,linerR],[3,SoilR],[4,AnoR]],
                            calcOnly=False, verbose=True, sr=False,
                            noiseLevel=noiseLevel, noiseAbs=noiseAbs, seed=1337)


    u_elecs.save('solution.data')


if __name__ == '__main__':

    world,landfill,ano,plc = create_plc()
    malmScheme, mesh = create_scheme(world,landfill,ano,plc)
    run_fwd_model(malmScheme, mesh)


