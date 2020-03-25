from scipy.interpolate import griddata
from matplotlib.backends.backend_pdf import PdfPages
import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import natsort as ns

def grid_xy(x,y):
    Xm = np.linspace(min(x), max(x), 50)
    Ym = np.linspace(min(y), max(y), 50)
    XI,YI = np.meshgrid(Xm,Ym)
    print('returning meshgrid...')
    return XI, YI


def plotCSD(coordinates,PestOutPar,RealSources):
    param = np.loadtxt(PestOutPar,skiprows=1,usecols=1)
    coord = np.loadtxt(coordinates)
    x, y = coord[:,0], coord[:,1]
    z = np.copy(param)
    if RealSources != None:
        print(RealSources)
        RSs = [(float(rc.split(',')[0]), float(rc.split(',')[1])) for rc in RealSources]
        numRealSources = int(len(RSs)  / 2)
    XI, YI = grid_xy(x, y)
    points = np.column_stack((x,y))
    grid = griddata(points, z, (XI, YI), method = 'linear') # grid is a np array
    print('plotting surface...')
    f = plt.figure('surface', figsize = (7,7))
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
    plt.scatter(x, y, facecolor=facecolor_norm_z, edgecolor=edgecolor_norm_z, cmap = 'jet')
    if RealSources != None:
        for s in RSs: # add the real sources as red dots
            Rx,Ry = s[0], s[1]
            plt.plot(Rx,Ry,'or')
    # set labels and fonts
    plt.ylabel('y [m]',fontsize=15)
    plt.xlabel('x [m]',fontsize=15)
    axes = plt.gca()
    axes.set_xlim([0,0.53])
    axes.set_ylim([0,0.52])
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.title(PestOutPar)
    return f


allpar = [p for p in os.listdir() if ('.par.') in p]

allpar = ns.natsorted(allpar, key=lambda y: y.lower())

ap = argparse.ArgumentParser(description = 'plot the Pareto front from pest.pod file')
ap.add_argument('--sc', metavar = 'source coords', type = str, help = 'coordinates of the sources, format = x1,y1 x2,y2', nargs = '*')
ap.add_argument('-vec', metavar = 'VRTe coords', type = str, help = 'file containing the coordinates of the virtual electrodes', default = 'VRTeCoord.txt')
args = ap.parse_args()
SC = args.sc
VRTeC = args.vec

# extract iteration to print
# find pst file and pareto section
pst = [f for f in os.listdir() if f.endswith('.pst')]
print(pst) 
pls = []
paretoSect =False
if len(pst) == 1:
    with open(pst[0]) as fid:
        for i,l in enumerate(fid):
            if '* pareto' in l:
                paretoSect = True
            if paretoSect and not '#' in l:
                pls.append(l)

pls = [l.replace('\n','') for l in pls if l != '\n']
regline = pls[2].split()
minreg, Maxreg, nIter = float(regline[0]), float(regline[1]), int(regline[2])
Iterline = pls[3].split()
Iter, allIter, lastIter = int(Iterline[0]), int(Iterline[1]), int(Iterline[2])

# calculate iterations to print

i2p = list(np.linspace(Iter, allIter*(nIter - 1) + Iter + lastIter, nIter+1, dtype = int))
print(i2p)
print(allpar)
with PdfPages('iCSD_Pareto.pdf') as pdf:
    for i in i2p:
        print(i)
        fname = allpar[i]

        f = plotCSD(VRTeC, fname, SC)
        pdf.savefig()
        plt.close(f)


# print pareto curve
pod = pst[0].replace('pst','pod')
data = pd.read_csv(pod, skiprows = 0, delim_whitespace = True )
plt.plot(data.iloc[i2p,0], data.iloc[i2p,1], 'ro-')
plt.xlabel('MALM data objective function')
plt.ylabel('Regularization objective function')
#
plt.ylim((0, 1))
#plt.xscale('log')
##plt.yscale('log')
plt.savefig('ParetoFront.png', dpi = 300)
plt.show()
