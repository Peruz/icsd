import matplotlib.pyplot as plt
import numpy as np
import pybert as pb
import pygimli as pg
import pygimli.meshtools as mt
import sys
import matplotlib

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 12}
matplotlib.rc('font', **font)

make_meshes = True # make or load
make_sequences = True # make or load bert schemes
invert = True # run inversion or load rho to plot
calc = True

# Rhizotron geometry
height = 0.53
width = 0.52
thickness = 0.02
# Mesh parameters
layers = 4
refinement = 0.001
#
VRTe_nums = range(65, 354)

plt.figure()
# Elec coordinates
elec_X = np.linspace(0.09, 0.44, num = 8)
elec_Y = np.linspace(0.43, 0.08, num = 8)
elec_z = 0.00
elec_xyz = []
for y in elec_Y:
    for x in elec_X:
            elec_xyz.append([x, y, elec_z])
            plt.scatter(x, y)
elec_xyz = np.array(elec_xyz)

# VRTe coordinates
VRTe_X = np.linspace(0.025, 0.505, num = 17)
VRTe_Y = np.linspace(0.02, 0.5, num = 17)
VRTe_z = 0.01
VRTe_xyz = []
for x in VRTe_X:
        for y in VRTe_Y:
            VRTe_xyz.append([x, y, VRTe_z])
            plt.scatter(x, y)
VRTe_xyz = np.array(VRTe_xyz)

# Source 
source_xyz = np.array([0.25, 0.27, 0.01])

# Get unique 2d coordinates
all_coord = np.vstack((elec_xyz, VRTe_xyz, source_xyz)) 

if calc:
    # Data containers
    scheme_source = pb.DataContainerERT()
    scheme_VRTe = pb.DataContainerERT()
    scheme_source.setSensorPositions(all_coord) 
    scheme_VRTe.setSensorPositions(all_coord) 


    ### MESH
    if make_meshes:
        try:
            from mesh import mk_mesh
        except:
            sys.exit('no mesh.py')
        else:
            mesh2d, mesh3d = mk_mesh(width = width, height = height, thickness = thickness,
                                     layers = layers, refinment = refinement,
                                     elec = elec_xyz, VRTe = VRTe_xyz, source = source_xyz)
        pg.show(mesh2d)
        pg.wait()

    elif not make_meshes:  
        mesh3d = pg.load('3d_transdimensional.bms')
        mesh3d = pg.load('mesh.bms')
        mesh2d = pg.load('2d.bms')

    ### SEQUENCES
    if make_sequences:
        try:
            from Rhi_204 import *
        except:
            sys.exit('there should be the short module with the functions to make the sequences')
        seq = mk_seq() # seq with 1
        VRTe_seq = mk_full_malm(seq, VRTe = VRTe_nums) # seq for all VRTe
        source_seq = mk_full_malm(seq, VRTe = 354) # seq with source

    test_seq = np.array([[353, 63, 1, 2], [353, 63, 2, 3]])

    ### BERT SCHEMES
    scheme_VRTe.resize(len(VRTe_seq))
    scheme_source.resize(len(source_seq))

    scheme_test = pb.DataContainerERT()
    scheme_test.resize(2)

    for i, j in enumerate("abmn"):
        scheme_VRTe.set(j, VRTe_seq[:, i] -1)
        scheme_source.set(j, source_seq[:,i] -1)
        scheme_test.set(j, test_seq[:, i] -1)

    scheme_VRTe.set("valid", np.ones(len(VRTe_seq)))
    scheme_source.set("valid", np.ones(len(source_seq)))
    scheme_test.set("valid", np.ones(len(test_seq)))

    #pb.exportData(scheme_source, 'source.txt')
    scheme_VRTe.setSensorPositions(all_coord) 
    scheme_source.setSensorPositions(all_coord) 
    scheme_test.setSensorPositions(all_coord) 

    ### SIMULATION
    # simulate both VRTe and source

    sim_source = pb.ERTManager()
    VRTe_source = pb.ERTManager()

    rho = []
    for cell in mesh3d.cells():
        if cell.center().x() < height / 3:
            rho.append(1)
        elif cell.center().x() > (2 * height / 3):
            rho.append(1)
        else:
            rho.append(1)
    rho = np.array(rho)

    scheme_source.save('source.txt')
    scheme_VRTe.save('VRTe.txt')
    source_sim_results = sim_source.simulate(mesh = mesh3d, res = rho, scheme = scheme_source)
    VRTe_sim_results = VRTe_source.simulate(mesh = mesh3d, res = rho, scheme = scheme_VRTe)
    source_sim_results.save('source_sim_results.ohm')
    VRTe_sim_results.save('VRTe_sim_results.ohm')


### Extract and save vectors for iCSD
def get_r(file_ohm):
    """
    read file with simulation results and calculate the resistances from k and rhoa

    parameters:
    ===
    file_ohm: name of the file to read (t: str, ext: ohm)

    returns:
    ===
    data_r: list with the data resistances (t: list)
    """
    # check args
    if type(file_ohm) != str:
        raise TypeError('### the name of the file has to be a string')
    if not file_ohm.endswith('.ohm'):
        raise ValueError('### file extension has to be .ohm')
    # funtion
    data = pb.load(file_ohm)
    data_k = data['k'].array()
    data_rhoa = data['rhoa'].array()
    # rhoa = R * k 
    data_r = data_rhoa / data_k
    return(data_r, data.sensorPositions())

source_r, _ = get_r('source_sim_results.ohm')
VRTe_r, sensPos = get_r('VRTe_sim_results.ohm')
source_r_np = np.array(source_r)
VRTe_r_np = np.array(VRTe_r)

VRTe_coord = []
for i, s in enumerate(sensPos):
    if i+1 in VRTe_nums: # add 1 to keep correct
        VRTe_coord.append([s[0], s[1], s[2]])
VRTe_coord_np = np.array(VRTe_coord)

np.savetxt('ObsData.txt', source_r_np)
np.savetxt('VRTeSim.txt', VRTe_r_np)
np.savetxt('VRTeCoord.txt', VRTe_coord_np, fmt = '%1.3f')

print('number of rows in ObsData (source simulation): ', len(source_r_np))
print('number of rows in VRTeSim: ', len(VRTe_r))
print('taking number of VRTe coord: ', len(VRTe_coord_np))

