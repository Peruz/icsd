import pybert as pb
import pygimli as pg
import pygimli.meshtools as mt
import numpy as np
from functools import reduce

check_meshes = True # plot 2d mesh

def mk_mesh(width = 53, height = 52, thickness = 0.02, layers = 4, refinment = 0.002, elec = None, VRTe = None, source = (0.265, 0.26, 0.01)):
    """
    make 2D mesh and associated 3D transdimensional mesh
    
    it seems there are some problems when VRTe and elec have some common x and y,
    even if z is different
    """

    all2Dcoordinates = np.vstack((elec[:, range(2)], VRTe[:, range(2)], source[0:2]))
    
    unq, count = np.unique(all2Dcoordinates, axis = 0, return_counts = True)
    if len(unq) != len(all2Dcoordinates):
        print('the following electrodes have the same x y coord: ', unq[count>1])
    unique2Dcoordinates =  np.vstack(set(tuple(row) for row in all2Dcoordinates))
    #unique2Dcoordinates = np.unique(all2Dcoordinates, axis = 0)
    zeros = np.zeros(len(unique2Dcoordinates))
    unique2Dcoordinates = np.column_stack((unique2Dcoordinates, zeros))
    
    print('unique', len(all2Dcoordinates))
    print('unique', len(unique2Dcoordinates))
    #for row in unique2Dcoordinates:
    #    print(row)
    #print(unique2Dcoordinates)
    #print('total nodes', len(all2Dcoordinates))
    #print('total unique 2d nodes', len(unique2Dcoordinates))

    # geom = mt.createRectangle(start = [0, 0], end = [width, height], boundaryMarker = 0, marker = 2)
    geom = mt.createRectangle(start = [0, 0], end = [width, height])

    for pos in unique2Dcoordinates:
        geom.createNode(pg.RVector3(pos))
        #geom.createNode(pos + pg.RVector3(0, refinment))
        #geom.createNode(pos - pg.RVector3(0, refinment))
        #geom.createNode(pos + pg.RVector3(refinment, 0))
        #geom.createNode(pos - pg.RVector3(refinment, 0))

    # Make 2D triangular mesh
    mesh2d = mt.createMesh(geom, quality = 33, area = 0.0003)

    # if check_meshes:
        #ax, c = pg.show(mesh2d)
        #ax.plot(pg.x(scheme_sim.sensorPositions()), pg.y(scheme_sim.sensorPositions()), "ko", markersize = 3)
        #pg.wait()

    mesh2d.exportVTK('2d.vtk')
    mesh2d.save('2d.bms')

    # 2D to 3D, give each cell individual marker (will be automatically translated in 3D)
    for cell in mesh2d.cells():
        cell.setMarker(cell.id())

    # 3D
    mesh3d = pg.createMesh3D(mesh2d, np.linspace(0, thickness, layers + 1))

    # PRINT SOME INFO
    string = "3D Mesh has %d cells but only %d inversion parameters."
    print("-" * 80)
    print(string % (mesh3d.cellCount(), len(np.unique(mesh3d.cellMarkers()))))
    print('mesh nodes:  ', mesh3d.nodeCount(), 'node markers:  ', set(mesh3d.nodeMarkers()))
    print("-" * 80)

    # Set refernce nodes in corners (necessary for closed geometries),
    # using existing nodes to keep consistent with the 2D mesh
    lower_left_node = mesh3d.findNearestNode([mesh3d.xmin(), mesh3d.ymin(), 0.0])
    mesh3d.node(lower_left_node).setMarker(pg.MARKER_NODE_CALIBRATION)
    
    lower_right_node = mesh3d.findNearestNode([mesh3d.xmax(), mesh3d.ymin(), 0.0])
    mesh3d.node(lower_right_node).setMarker(pg.MARKER_NODE_REFERENCEELECTRODE)

    # mark source as electrode
    source_node = mesh3d.findNearestNode(source)
    mesh3d.node(source_node).setMarker(pg.MARKER_NODE_ELECTRODE)
    
    # mark VRTe and elec as elec
    for VRTe_row in VRTe:
        nearest = mesh3d.findNearestNode([VRTe_row[0], VRTe_row[1], VRTe_row[2]])
        mesh3d.node(nearest).setMarker(pg.MARKER_NODE_ELECTRODE)

    for elec_row in elec:
        nearest = mesh3d.findNearestNode([elec_row[0], elec_row[1], elec_row[2]])
        mesh3d.node(nearest).setMarker(pg.MARKER_NODE_ELECTRODE)
    
    print("-" * 80)
    print('\nafter adding electrodes:')
    print('mesh nodes:  ', mesh3d.nodeCount(), 'node markers:  ', set(mesh3d.nodeMarkers()))
    print("-" * 80)
  
    # print the coordnates of the electrodes
    elec_99 = 0
    for n in mesh3d.nodeMarkers():
        if n == -99:
            elec_99 += 1
    print('number of electrodes', elec_99)
  
    xn = [n.pos()[0] for n in mesh3d.nodes()]
    yn = [n.pos()[1] for n in mesh3d.nodes()]
    zn = [n.pos()[2] for n in mesh3d.nodes()]
    mn = mesh3d.nodeMarkers()
    #for i,n in enumerate(mn):
    #    if n == -99:
    #        print(xn[i], yn[i], zn[i])

    mesh3d.exportVTK("3d_transdimensional.vtk")
    mesh3d.save("3d_transdimensional.bms")
    return(mesh2d, mesh3d)