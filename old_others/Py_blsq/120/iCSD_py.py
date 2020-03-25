import numpy as np
from scipy.optimize import lsq_linear
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import argparse

ap = argparse.ArgumentParser(description = 'Discrete current source density inversion')
ap.add_argument(
    '-VRTeSim', metavar = 'BERT_VRTe_calculation', type = str, help = 'file containing the forward calculation for each virtual source, going into the A matrix', default = 'VRTeSim.txt')
ap.add_argument(
    '-Obs', metavar = 'Observed_MALM_data', type = str, help = 'measured MALM data, b column of the linear system', default = 'ObsData.txt')
ap.add_argument(
    '-VRTeCoord', metavar = 'VRTe_coordinates', type = str, help = 'coordinates of the VRTe for plotting', default = 'VRTeCoord.txt')
ap.add_argument(
    '-wr', metavar = 'weight regularization', type = float, help = 'first derivative regularization weight', default = 0.1, nargs = '*')
ap.add_argument(
    '-wc', metavar = 'weight constrain', type = float, help = 'current conservation constrain, sum current fractions = 1', default = 100)
ap.add_argument(
    '--sc', metavar = 'source coords', type = str, help = 'coordinates of the sources, format = x1,y1 x2,y2', nargs = '*')
ap.add_argument(
    '--pareto', metavar = 'pareto mode', type = str, help = 'activate the pareto mode and gets lower, upper and number of regularization weights', nargs = '*')
ap.add_argument(
    '-nx', metavar = 'number of VRTe on x', type = int, help = 'number of VRTe for each row in x direction', default = 10)
ap.add_argument(
    '-ny', metavar = 'number of VRTe on y', type = int, help = 'number of VRTe for each row in y direction', default = 12)

""" functions """

def mkRegMatrix(nVRTe, nx, ny):
    reg = []
    vrte = range(1, nVRTe + 1)
    vrte = np.reshape(vrte,(ny, nx))
    for y in range(ny):
        for x in range(nx):
            minus = vrte[y, x]
            if x + 1 in range(nx):
                plus = vrte[y , x + 1]
                row = np.zeros(nVRTe, int)
                row[minus -1] = - 1
                row[plus -1] = + 1
                reg.append(row)
            if y + 1 in range(ny):
                plus = vrte[y + 1, x]
                row = np.zeros(nVRTe)
                row[minus -1] = - 1
                row[plus -1] = + 1
                reg.append(row)
    reg = np.array(reg)
    return reg


def grid_xy(x,y):
    Xm = np.linspace(min(x), max(x), 50)
    Ym = np.linspace(min(y), max(y), 50)
    XI,YI = np.meshgrid(Xm,Ym)
    print('returning meshgrid...')
    return XI, YI


def plotCSD(VRTeCoord, VRTeSolution, RealSources = None):

    x, y = VRTeCoord[:, 0], VRTeCoord[:, 1]
    z = np.copy(VRTeSolution)

    if RealSources != None:
        print(RealSources)
        RSs = [(float(rc.split(',')[0]), float(rc.split(',')[1])) for rc in RealSources]
        numRealSources = int(len(RSs)  / 2)
    
    XI, YI = grid_xy(x, y)
    points = np.column_stack((x,y))
    grid = griddata(points, z, (XI, YI), method = 'nearest') # grid is a np array

    print('plotting surface...')
    plt.figure('surface')
    dict_figures = {1:"surface"} # make a dict to keep tract of the figures, add first fig "Surface"
    plt.imshow(grid, extent = (min(x),max(x),min(y),max(y)),aspect='auto', origin='lower',cmap='jet')

    cbar = plt.colorbar()
    cbar.set_label('Fraction of Current Source', labelpad = 10, fontsize=14)

    # set colors for edge and face
    norm_z = (z - min(z)) / (max(z) - min(z))
    edgecolor_norm_z = plt.cm.Greys(norm_z)
    facecolor_norm_z = plt.cm.jet(norm_z)

    # go back to Surface fig
    plt.figure('surface')
    gcf_num=plt.gcf().number
    print("back to: ", dict_figures[gcf_num])
    # overlap data and surface
    if RealSources != None:
        for s in RSs: # add the real sources as red dots
            Rx,Ry = s[0], s[1]
            plt.plot(Rx, Ry,'ow')
    plt.scatter(x, y, facecolor=facecolor_norm_z, edgecolor=edgecolor_norm_z, cmap = 'jet')
    # set labels and fonts
    plt.ylabel('y [m]',fontsize=15)
    plt.xlabel('x [m]',fontsize=15)
    axes = plt.gca()
    axes.set_xlim([0,0.53])
    axes.set_ylim([0,0.52])
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.tight_layout()
    plt.savefig('iCSD',dpi=350)
    plt.show()

def iCSD(A, b, nVRTe, constrain_weight, regularization_weight, nx, ny):
    A = A.reshape((-1, nVRTe), order = 'F') # reshape by columns (Fortran-like)
    w = np.ones(A.shape[0])
    
    # add constrain sum x = 1
    cA = np.ones(A.shape[1])
    cb = np.array([1])
    cw = np.ones(1) * constrain_weight # weight of sum VRTe = 1
    
    A_c = np.vstack((A, cA))
    b_c = np.concatenate((b, cb))
    w_c = np.concatenate((w, cw))
    
    # make regularization der x = 0
    
    # make regularization matrix
    rA = mkRegMatrix(nVRTe = nVRTe, nx = nx, ny = ny)
    # add 0 to b as regularization target
    rb = np.zeros(rA.shape[0])
    # regulaztion weights vector
    rw = np.ones(rA.shape[0]) * regularization_weight
    
    # add regularization
    A_cr = np.vstack((A_c, rA))
    b_cr = np.concatenate((b_c, rb))
    w_cr = np.concatenate((w_c, rw))
    
    # make weight matrix for A
    W = np.zeros((A_cr.shape[0], A_cr.shape[0]))
    np.fill_diagonal(W, w_cr)
    
    # applied weights
    A_crw = np.dot(W, A_cr)
    b_crw = w_cr * b_cr # these remain 0 in this case 
    
    # solve with lsq_linear 
    x = lsq_linear(A_crw, b_crw, bounds = (0, 1)) 
    
    return x

def ResidualAnalysis(x, constrain_weight, regularizatio_weight):
    fitting_res = x.fun[:nData]
    constrain_res = x.fun[nData + 1] / constrain_weight
    regularization_res = x.fun[nData + 2:] / regularization_weight # this constrain is not included in the reg function
    reg_sum = np.sum(np.square(regularization_res)) # ||Ax - b||2
    fit_sum = np.sum(np.square(fitting_res)) # ||x - x0||2 
 
    return reg_sum, fit_sum
 
""" main """

# parse parameters
args = ap.parse_args()
VRTeSim = args.VRTeSim
Obs = args.Obs
VRTeCoord = args.VRTeCoord
regularization_weight = args.wr
constrain_weight = args.wc
SC = args.sc
pareto = args.pareto
nx = args.nx
ny =args.ny

# load files
VRTeCoord = np.loadtxt(VRTeCoord)
b = np.loadtxt(Obs)
A = np.loadtxt(VRTeSim)

# check file dimensions
nVRTe = VRTeCoord.shape[0]
nData = b.shape[0]
nSim = A.shape[0]
b_expected = nSim / nVRTe
if b_expected != nData:
    sys.exit('lenghts of the input files do not agree')

# current source inversion
if pareto is not None:
    pareto_fit_res = []
    pareto_reg_res = []
    # run iCSD multiple times while changing the weights to explore the L-curve
    w_pareto = np.linspace(float(pareto[0]), float(pareto[1]), float(pareto[2]))
    print('pareto weights are\n', w_pareto)
    for w in w_pareto:
        x = iCSD(A, b, nVRTe, constrain_weight, w, nx, ny)
        fr, rr = ResidualAnalysis(x, constrain_weight, regularization_weight)
        pareto_fit_res.append(fr)
        pareto_reg_res.append(rr)
    plt.loglog(pareto_fit_res, pareto_reg_res, 'or')
    plt.xlabel('Residual')
    plt.ylabel('Roughness')
    plt.show()
    print(pareto_fit_res, pareto_reg_res)

else:
    #  run single regularization current source density inversion
    x = iCSD(A, b, nVRTe, constrain_weight, regularization_weight, nx, ny)
    ResidualAnalysis(x, constrain_weight, regularization_weight)
    # plot solution CSD distribution
    plotCSD(VRTeCoord, x.x, SC)
