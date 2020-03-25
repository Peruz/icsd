import numpy as np
import os

def mkMALMskip0AB(A,B):
    A,B = np.array(A), np.array(B)
    # define matrix
    electrodes = np.arange(1,65).reshape((8,8))
    electrodes_flip = np.fliplr(electrodes)
    electrodes_T = np.transpose(electrodes)
    # horizontal and vertical
    hor = electrodes.repeat(2, axis = 1)[:,1:-1].reshape(-1,2)
    #print('horizontal: ', hor)
    ver = electrodes_T.repeat(2, axis = 1)[:,1:-1].reshape(-1,2)
    #print('vertical: ', ver)
    # diagonals in the two directions
    diag1 = np.vstack(electrodes.diagonal(i).repeat(2)[1:-1].reshape(-1,2) for i in range(-6,7))
    #print('diag1', diag1)
    diag2 = np.vstack(electrodes_flip.diagonal(i).repeat(2)[1:-1].reshape(-1,2) for i in range(-6,7)) 
    #print('diag2', diag2)
    # combine and remove rows containing A or B
    seqMN_all = np.vstack((hor,ver,diag1,diag2))
    #print('seqMN all: ',seqMN_all)
    mask = np.any(np.equal(seqMN_all,1) | np.equal(seqMN_all,64), axis = 1)
    seqMN = seqMN_all[~mask]
    # current electrode columns
    seqAB = np.ones_like(seqMN)
    seqAB[:,0] *= A
    seqAB[:,1] *= B
    # combine columns
    seq = np.column_stack((seqAB,seqMN))
    #print('seq', seq)
    print('a: ', A,'b: ', B,'\nafter removing potential dipoles that contain current electrodes: ',seq.shape[0])
    return(seq)


seqComplete = np.array([0,0,0,0])
for e in range(65, 371):
    seqAB = mkMALMskip0AB(e, 64)
    seqComplete = np.vstack((seqComplete, seqAB))
seqComplete = seqComplete[1:,:]

print(seqComplete)
np.savetxt('completeMALM.txt', seqComplete, '%i %i %i %i')
