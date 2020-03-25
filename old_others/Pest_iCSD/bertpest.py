import numpy as np
import re
#import pygimli as pg
import pybert as pb
import os
import sys
from subprocess import call
import matplotlib.pyplot as plt
from pygimli.meshtools import readGmsh
from scipy.interpolate import griddata
import matplotlib.patheffects as PathEffects

""" Use module importlib, function reload: importlib.reload(bp) """

def mkShm(coordinates,sequence):
    """
    load electrodes coordinates, e.g.,"coordinates.txt" and sequence, combines all in a BERT format file .shm, ready for DCMOD
    """
    OutNameSeq = os.path.splitext(sequence)[0] + '.shm' 
    # clean a bit
    if os.path.isfile(OutNameSeq):
        os.remove(OutNameSeq)
    # load files
    coord = np.loadtxt(coordinates)
    seq = np.loadtxt(sequence)
    # get length of the files
    num_coord = np.asarray(coord.shape[0], dtype = 'i')
    num_meas = np.asarray(seq.shape[0], dtype = 'i')
    # write file.shm
    with open (OutNameSeq,'ab') as f_handle:
        np.savetxt (f_handle,num_coord[None], fmt = '%i')
        np.savetxt (f_handle,np.array(['# x y z']), fmt = '%s')
        np.savetxt (f_handle,coord, fmt = '%3.3f')
        np.savetxt (f_handle,num_meas[None], fmt = '%i')
        np.savetxt (f_handle,np.array(['# a b m n']), fmt = '%s')
        np.savetxt (f_handle,seq, fmt = '%i %i %i %i')

def mkMALMskip0(A,B):
    """
    make a ''skip0 MALM'' sequence for Rhizotron configuration
    take care: using transpose for verticals, flip for diagonals
    """
    
    A,B = np.array(A), np.array(B)
    # define matrix
    electrodes = np.arange(1,65).reshape((8,8))
    electrodes_flip = np.fliplr(electrodes)
    electrodes_T = np.transpose(electrodes)
    # horizontal and vertical
    hor = electrodes.repeat(2, axis = 1)[:,1:-1].reshape(-1,2)
    print(hor)
    ver = electrodes_T.repeat(2, axis = 1)[:,1:-1].reshape(-1,2)
    print(ver)
    # diagonals in the two directions
    diag1 = np.vstack(electrodes.diagonal(i).repeat(2)[1:-1].reshape(-1,2) for i in range(-6,7))
    print(diag1.shape[0])
    diag2 = np.vstack(electrodes_flip.diagonal(i).repeat(2)[1:-1].reshape(-1,2) for i in range(-6,7)) 
    print(diag2.shape[0])
    # combine and remove rows containing A or B
    seqMN_all = np.vstack((hor,ver,diag1,diag2))
    print('seqMN all: ',seqMN_all)
    mask = np.any(np.equal(seqMN_all,A) | np.equal(seqMN_all,B), axis = 1)
    seqMN = seqMN_all[~mask]
    # current electrode columns
    seqAB = np.ones_like(seqMN)
    seqAB[:,0] *= A
    seqAB[:,1] *= B
    # combine columns
    seq = np.column_stack((seqAB,seqMN))
    np.savetxt('sequence.txt',seq, fmt = '%i %i %i %i')
    print('after removing potential dipoles that contain current electrodes: ',seq.shape[0])
    return(seq)


def mkMALMseqs(MALMelec, VRTelec):
    """ 
    takes the couple of lab electrodes (A,B) and makes the suitable MALM and VRTeSeq sequences substituting A with each VRTe
    MALMelec = (A,B); among the real electrodes, they are indeed the electrodes used for current injection during MALM
    VRTelec = (VRTeFirst, VRTeLast) first and last of VRTe numbers, as in the BERT sequence; e.g., VRTe = [65,66,67,....201]
    """ 
    OutNameVRTeSeq = 'VRTeSeq.txt'
    OutNameLabSynSeq = 'MALMSeq.txt'
    # clean a bit
    if os.path.isfile(OutNameVRTeSeq):
        os.remove(OutNameVRTeSeq)
    if os.path.isfile(OutNameLabSynSeq):
        os.remove(OutNameLabSynSeq)
    A, B = MALMelec[0], MALMelec[1]
    VRTeFirst, VRTeLast = VRTelec[0], VRTelec[1]
    VRTe = np.arange(VRTeFirst, VRTeLast+1) # define number of the VRTe electrodes in the sequence
    
    MALM = mkMALMskip0(A,B)
    # LABseqSyn similar to LAB but modify number of one to match the BERT sequence and synthetic coordinates
    # this is done because we prefare to use a different node rather than moving the node (different mesh and blabla... mess)
    # i.e., in the lab we move the electrode, in the synthetic simulation we change node.
    # Note, this is done only if a pure synthetic acquisition and iCSD is desidered. Otherwise just use lab data and VRTeseq
    
    # VRTe sequences, create long sequence containing all the acquisitions for each VRTe
    NumVRTe = VRTe.shape[0]
    seqLen = MALM.shape[0]
    VRTeseq = np.tile(MALM,(NumVRTe,1))
    VRTeseqLen = VRTeseq.shape[0]
    # substitute of the electrode to change
    for i in range(NumVRTe):
        VRTeseq[i*seqLen:(i+1)*seqLen,0] = VRTe[i] 
    np.savetxt('VRTeSeq.txt',VRTeseq,'%i %i %i %i')
    np.savetxt('MALMSeq.txt',MALM,'%i %i %i %i')

def readBERT(infile):
    coordinates = []
    data = []
    read_coord = 0
    read_data = 0
    fid = open (infile)
    for line in fid:
        if line[0] == '#':
            if line[2] == 'x':
                read_coord = 1 # start coordinates
            if line[2] == 'a':
                read_coord = 0 # stop coordinates
                read_data = 1 # start data
        else:
            if read_coord == 1:
                line_coord = re.findall('[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?',line)
                coordinates.append(line_coord)
            if read_data == 1:
                line_data = re.findall('[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?',line)
                data.append(line_data)
    del coordinates[-1]
    del data[-1]
    coord = np.array(coordinates, dtype = 'f' )
    num_elec = coord.shape[0] 
    data = np.array(data, dtype = 'f')
    num_data = data.shape[0] 
    return(coord,data,num_elec,num_data)

def plotData(infile,interval):
    """ reads and plot data (reads BERT format using "readBERT") """
    coord,data,num_elec,num_data = readBERT(infile)
    num_VRTe = num_data/interval
    checkLenghts = num_VRTe.is_integer()
    if not checkLenghts:
        sys.exit('checkLenghts failled, data length/VRTe is not an integer')
    else:
        num_VRTe_Int = np.array(num_VRTe, dtype ='i')
        vector_data_num = np.arange(interval).reshape(-1,1)
    for i in range(num_VRTe_Int):
        first = i * interval
        last = interval * (i+1)
        plt.plot(vector_data_num,data[first:last,4],'o')
    plt.show()

def mkBMS(mshfile):
    # if present, remove previous file
    if os.path.isfile('./mesh.bms'):
        os.remove('./mesh.bms')
    # take input name
    print ('convert: ',fname)
    # read and convert file
    mesh=readGmsh(fname, verbose=True)
    # save coverted file
    mesh.save('mesh.bms')

def plotCSD(coordinates,PestOuttPar,RealSourceXY):
    param = np.loadtxt(PestOuttPar,skiprows=1,usecols=1)
    coord = np.loadtxt(coordinates)
    x, y = coord[:,0], coord[:,1]
    z = np.copy(param)
    Rx,Ry = RealSourceXY[0], RealSourceXY[1]
    print(Rx)
    XI, YI = grid_xy(x, y)
    points = np.column_stack((x,y))
    grid = griddata(points, z, (XI, YI), method = 'linear') # grid is a np array
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
    plt.scatter(x, y, facecolor=facecolor_norm_z, edgecolor=edgecolor_norm_z, cmap = 'jet')
    plt.plot(Rx,Ry,'or')
    # set labels and fonts
    plt.ylabel('y [m]',fontsize=15)
    plt.xlabel('x [m]',fontsize=15)
    axes = plt.gca()
    axes.set_xlim([0,0.53])
    axes.set_ylim([0,0.52])
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.savefig('iCSD',dpi=350)
    plt.show()
def plotCSD5(coordinates,PestOuttPar,RealSourceXY):
    param = np.loadtxt(PestOuttPar,skiprows=1,usecols=1)
    coord = np.loadtxt(coordinates)
    x, y = coord[:,0], coord[:,1]
    z = np.copy(param)
    Rx1,Ry1 = RealSourceXY[0], RealSourceXY[1]
    Rx2,Ry2 = RealSourceXY[2], RealSourceXY[3]
    Rx3,Ry3 = RealSourceXY[4], RealSourceXY[5]
    Rx4,Ry4 = RealSourceXY[6], RealSourceXY[7]
    Rx5,Ry5 = RealSourceXY[8], RealSourceXY[9]
    XI, YI = grid_xy(x, y)
    points = np.column_stack((x,y))
    grid = griddata(points, z, (XI, YI), method = 'linear') # grid is a np array
    print('plotting surface...')
    plt.figure('surface')
    dict_figures = {1:"surface"} # make a dict to keep tract of the figures, add first fig "Surface"
    plt.imshow(grid, extent = (min(x),max(x),min(y),max(y)),aspect='auto', origin='lower',cmap='viridis')

    cbar = plt.colorbar()
    cbar.set_label('Fraction of Current Source', labelpad = 10, fontsize=14)

    # set colors for edge and face
    norm_z = (z - min(z)) / (max(z) - min(z))
    edgecolor_norm_z = plt.cm.Greys(norm_z)
    facecolor_norm_z = plt.cm.viridis(norm_z)

    # go back to Surface fig
    plt.figure('surface')
    gcf_num=plt.gcf().number
    print("back to: ", dict_figures[gcf_num])
    # overlap data and surface
    plt.scatter(x,y,facecolor=facecolor_norm_z,edgecolor=edgecolor_norm_z)
    plt.plot(Rx1,Ry1,'or')
    plt.plot(Rx2,Ry2,'or')
    plt.plot(Rx3,Ry3,'or')
    plt.plot(Rx4,Ry4,'or')
    plt.plot(Rx5,Ry5,'or')
    # set labels and fonts
    plt.ylabel('y [m]',fontsize=15)
    plt.xlabel('x [m]',fontsize=15)
    axes = plt.gca()
    axes.set_xlim([0,0.53])
    axes.set_ylim([0,0.52])
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.savefig('iCSD',dpi=350)
    plt.show()

def grid_xy(x,y):
    Xm = np.linspace(min(x), max(x), 50)
    Ym = np.linspace(min(y), max(y), 50)
    XI,YI = np.meshgrid(Xm,Ym)
    print('returning meshgrid...')
    return XI, YI

