# -*- coding: utf-8 -*-
"""
Created on Mon May 11 17:29:01 2020
@author: Benjamin
"""

import numpy as np
from scipy.sparse import diags
import matplotlib.pyplot as plt


#%%
def nx_ny(coord):
    """find number of nodes in each direction, has to be a regular grid"""
    nx = np.unique(np.round(coord[:, 0], 3)).shape[0]
    ny = np.unique(np.round(coord[:, 1], 3)).shape[0]
    
    return nx,ny
    
def nx_ny_nz(coord):
    """find number of nodes in each direction, has to be a regular grid"""
    nx = np.unique(np.round(coord[:, 0], 3)).shape[0]
    ny = np.unique(np.round(coord[:, 1], 3)).shape[0]
    nz = np.unique(np.round(coord[:, 2], 3)).shape[0]
    
    return nx,ny,nz

#%% Smoothing methods for different mesh types

#%% 2d

def regularize_A(coord,nVRTe):
    """create and append rows for to A, 
    for spatial regularization (simple model smoothing). 
    Working only on 2d regular meshes
    """
    reg = []
    vrte = range(1, nVRTe + 1)
    
    nx,ny= nx_ny(coord)

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
    reg_A = np.array(reg)
    
    return reg_A
        
def regularize_A_x_y(coord,alphaSx,alphaSy):
    """create and append rows for spatial regularization to A, 
    second derivative is applied in both direction x and y
    math:: Dx = ?? Dy=
    We used a ponderate diagonal matrix with coeffcient (1,-2, 1)
    """
    nx,ny= nx_ny(coord)

    ncells = nx*ny;       
    Dx=diags([1, -2, 1], [-1, 0, 1], shape=(ncells, ncells)).todense()
    idx=np.arange(0,len(Dx),nx)
    Dx[idx,:] = 0
    idx2=np.arange(nx,len(Dx),nx)
    Dx[idx2,:] = 0
    
    Dy=diags([1, -2, 1], [-nx, 0, nx], shape=(ncells, ncells)).todense()
    idy=np.arange(0,nx)
    Dy[idy,:] = 0
    idy2=np.arange(((ny-1)*nx+1),(nx)*(ny))
    Dy[idy2,:] = 0
    
    reg_Ax = alphaSx*np.array(Dx.transpose()*Dx)
    reg_Ay = alphaSy*np.array(Dy.transpose()*Dy)
    
    return reg_Ax, reg_Ay
    
#%% 3d
def regularize_A_3d(nVRTe,coord):
    """model smoothing consisting in creating and 
    appending rows for spatial regularization to A
    """
    nx,ny,nz= nx_ny_nz(coord)
    reg = []
    vrte = range(1, nVRTe + 1)
    vrte = np.reshape(vrte,(ny, nx, nz))
    for z in range(nz):
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
    reg_A = np.array(reg)
    
    return reg_A
        
def regularize_A_UnstructuredMesh2d(coord,nVRTe,k_neighbors=2): 
    """model smoothing consisting in creating and appending rows for spatial regularization to A. 
    Adapted for unstructured mesh since it uses the k_neighbors method, default k=2. Also working on regular grid 2d"""
    reg = []
    for VRTEnb in range(nVRTe):
        dist =  np.linalg.norm(coord[VRTEnb]-coord, axis=1)
        closest = np.argsort(dist)
        k = k_neighbors  # For each point, find the k closest current sources
        Ind = closest[1:k+1]
        row = np.zeros(nVRTe) # add a line to the regularisation A with k non-null coefficients
        knorm = dist[closest[1:k+1]]/dist[closest[1:k+1]].sum(axis=0,keepdims=1)
        row[Ind]= -knorm
        row[VRTEnb]= 1 # = one for the actual current source
        reg.append(row)
        test=[1]
        mask = np.in1d(test, VRTEnb)
        if mask.any()==True: 
            print('nll')
            # self.fc = plt.figure('TEST regularisation')
            # ax = self.fc.add_subplot(111, projection='3d')
            # ax.scatter(self.coord[VRTEnb,0], self.coord[VRTEnb,1], self.coord[VRTEnb,2], linewidths=12,
            #            facecolor = 'green', edgecolor = 'green')
            # ax.scatter(self.coord[Ind,0], self.coord[Ind,1], self.coord[Ind,2], linewidths=12,
            #            facecolor = 'red', edgecolor = 'red')
            # ax.set_xlim([min(self.coord_x),max(self.coord_x)])
            # ax.set_ylim([min(self.coord_y),max(self.coord_y)])
            # ax.set_zlim([min(self.coord_z),max(self.coord_z)])
            # self.fc.savefig(self.path2save+ 'TEST regularisation', dpi = 600)
            # #plt.show()
            # #plt.close()
        reg_A = np.array(reg)
        
    return reg_A


def regularize_A_UnstructuredMesh3d(coord,nVRTe,k_neighbors=9): 
    """model smoothing consisting in creating and appending rows for spatial regularization to A. 
    Adapted for unstructured mesh since it uses the k_neighbors method, default k=4. Also working on regular grid 2d"""
    reg = []
    for VRTEnb in range(nVRTe):
        dist =  np.linalg.norm(coord[VRTEnb]-coord, axis=1)
        closest = np.argsort(dist)
        k = k_neighbors  # For each point, find the k closest current sources
        # print('kkkkkkkkkkkkkkkkkkkkkkkk')
        # print(k_neighbors)
        Ind = closest[1:k+1]
        row = np.zeros(nVRTe) # add a line to the regularisation A with k non-null coefficients
        knorm = dist[closest[1:k+1]]/dist[closest[1:k+1]].sum(axis=0,keepdims=1)
        row[Ind]= -knorm
        row[VRTEnb]= 1 # = one for the actual current source
        reg.append(row)
        test=[1]
        mask = np.in1d(test, VRTEnb)
        if mask.any()==True: 
            print('nll')
            # self.fc = plt.figure('TEST regularisation')
            # ax = self.fc.add_subplot(111, projection='3d')
            # ax.scatter(self.coord[VRTEnb,0], self.coord[VRTEnb,1], self.coord[VRTEnb,2], linewidths=12,
            #             facecolor = 'green', edgecolor = 'green')
            # ax.scatter(self.coord[Ind,0], self.coord[Ind,1], self.coord[Ind,2], linewidths=12,
            #             facecolor = 'red', edgecolor = 'red')
            # ax.set_xlim([min(self.coord_x),max(self.coord_x)])
            # ax.set_ylim([min(self.coord_y),max(self.coord_y)])
            # ax.set_zlim([min(self.coord_z),max(self.coord_z)])
            # self.fc.savefig(self.path2save+ 'TEST regularisation', dpi = 600)
            #plt.show()
            #plt.close()
        reg_A = np.array(reg)
        
    return reg_A


def regularize_A_x_y_z(coord):
    """Model smoothing in 3d, not tested not working"""
    # self.estimateM0()
    nx,ny,nz= nx_ny_nz(coord)

    reg_Axz = []
    reg_Ayz = []
    for z in range(nz):
        ncells = nx*ny;       
        Dx=diags([1, -2, 1], [-1, 0, 1], shape=(ncells, ncells)).todense()
        idx=np.arange(0,len(Dx),nx)
        Dx[idx,:] = 0
        idx2=np.arange(nx,len(Dx),nx)
        Dx[idx2,:] = 0
        
        Dy=diags([1, -2, 1], [-nx, 0, nx], shape=(ncells, ncells)).todense()
        idy=np.arange(0,nx)
        Dy[idy,:] = 0
        idy2=np.arange(((ny-1)*nx+1),(nx)*(ny))
        Dy[idy2,:] = 0
        
        reg_Ax = alphaSx*np.array(Dx.transpose()*Dx)
        reg_Ay = alphaSy*np.array(Dy.transpose()*Dy)

        reg_Axz.append(reg_Ax)
        reg_Ayz.append(reg_Ay)

    reg_Ax= np.reshape(reg_Axz,[160,20])
    reg_Ay= np.reshape(reg_Ayz,[160,20])
    
    return reg_Ax, reg_Ay

#%% Initiate vectors to build regularisation matrice for A, b

def regularize_b(reg_A):
    """initiate vector b with zeros, the length is determined by the number of regul rows in A"""
    reg_b = np.zeros(reg_A.shape[0])
    
    return reg_b

def regularize_w(reg_A,wr,x0_prior,**kwargs):
    """create vector with weights, the length is determined by the number of regul rows in A such as
    .. math :: A = (G'*Wd*G + lambda*Wm)  
               b = G'*Wd*d + lambda*Wm*m0;   
     """
    if x0_prior==True:

        print('reg Wm (smallness + spatial reg) * lambda=' + str(wr))
        reg_w_0_b = np.ones(reg_A.shape[0]) * kwargs.get('x0') * wr
        reg_w_0_A = np.ones(reg_A.shape[0])* wr
        # reg_w_0_A = np.ones(reg_A.shape[0]) * kwargs.get('x0') * wr

        return reg_w_0_b, reg_w_0_A
        
    else:
        
        reg_w = np.ones(reg_A.shape[0]) * wr
        
        return reg_w

            
#%% RELATIVE SMALLNESS conditions (m-m0)
        
def ponderate_smallnessX0(alphaSxy,alphax0,**kwargs):
    """ Create relative smallness instance and applied smallness coefficient (\alpha_{x_{0}}) weight 
    
    .. math :: X_{0} = A*\alpha_{x_{0}}

    Parameters
    ------------
    self
    """
    # Smoothing matrices size are different from normal regularisation to anisotropic (alphaSxy)
    # Need to differentiate the two cases
    # In both cases the smoothing matrices are ponderated by alphax0
    if alphaSxy==True:
        reg_smallx0 = np.ones(kwargs.get('reg_Ax').shape)*alphax0
    else:
        reg_smallx0 = np.ones(kwargs.get('reg_A').shape)*alphax0 
    
    return reg_smallx0
    
def sum_smallness_smoothness(alphaSxy,x0_prior,**kwargs):
    """sum smallness and spatial regularisation

    .. math:: W_{m}=\alpha_{s}I+{D_{x}}^{T}D_{x} + D_{z}}^{T}D_{z}

    Parameters
    ------------
    self
    """
    # Test all the 4 possible cases i.e. 
    # with/without smallness regularisation (x0_prior)
    # with/without anisotropic smoothing (alphaSxy)
    # sum all contributions
    if (alphaSxy==True and x0_prior==False):
        # sum reg Ax, reg Ay
        reg_A_ss= kwargs.get('reg_Ax') + kwargs.get('reg_Ay')
    
    elif (alphaSxy==True and x0_prior==True):
        # sum small x0, reg Ax, reg Ay 
        reg_A_ss= kwargs.get('reg_Ax') + kwargs.get('reg_Ay')+ kwargs.get('reg_smallx0')

    elif (alphaSxy==False and x0_prior==True):
        # reg_A= reg_A + small x0
        reg_A_ss= kwargs.get('reg_A') + kwargs.get('reg_smallx0')
        
    elif (alphaSxy==False and x0_prior==False):
        reg_A_ss= kwargs.get('reg_A')
                
    return reg_A_ss
        