import re
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.interpolate import griddata
from scipy.linalg import lu
from scipy.optimize import lsq_linear
from matplotlib.backends.backend_pdf import PdfPages
import sys
import os

class iCSD2d_Class():
    """ Create a icsd inversion object.
    
    Parameters
    ----------
    coord_file : str, mandatory
        coordinates of the VRTe for plotting
    wr : float, optional
        Weight regularization
    wc : float, optional
        The height of the instrument above the ground. Can be specified per coil
        in the .csv.
    """
    def __init__(self,dirName):
        """ these functions are called only once, even for pareto,
        as they do not depend on the regularization weight"""
        # load
        self.dirName = dirName
        self.sim='VRTeSim.txt'
        self.obs='ObsData.txt'
        self.coord_file='VRTeCoord.txt'
        self.wr=25 #weight regularization
        self.wc=10000 #current conservation constrain, sum current fractions = 1
        self.sc=[] #coordinates of the sources, format = x1,y1 x2,y2'
        self.retElec=None #coordinates of the return electrode, format = x1,y1')
        self.pareto=False #if True run many icsd to explore the Pareto front
        self.errRmin=1 #min R to which the err is applied before passing to constant error
        self.pareto_MinErr=0.001
        self.pareto_MaxErr=1
        self.pareto_nSteps=10
        self.obs_err='const'
        
    def icsd_init(self):
        """ these functions are called only once, even for pareto,
        as they do not depend on the regularization weight"""
        self.createdirs()
        # load
        self.load_coord()
        self.load_obs()
        self.load_sim()
        # check vector sizes
        self.check_nVRTe()
        # # reshape VRTe vector into matrix A 
        self.reshape_A()
        # # define mode to weights the data        
        self.obs_w_f()
        # # constrain (curent conservation)
        self.con_A_f()
        self.con_b_f()
        self.con_w_f()
        # get VRTe grid size and make grid
        self.nx_ny()
        self.mkGrid_XI_YI()
        # # append spatial regularization
        self.regularize_A()
        self.regularize_b()
        # # stack data, constrain, and regularization 
        self.stack_A()
        self.stack_b()

    ### mkdirs 
    def createdirs(self):
        self.path2load = self.dirName
        
        self.path2save= self.dirName + 'fig/'
        try:
            # Create target Directory
            os.mkdir(self.path2save)
            print("Directory " , self.path2save ,  " Created ") 
        except FileExistsError:
            print("Directory " , self.path2save ,  " already exists")
            
    ### LOAD 
    def load_coord(self):
        print('Load coordinates')
        self.coord = np.loadtxt(self.path2load + self.coord_file)
        self.coord_x, self.coord_y =  self.coord[:, 0],  self.coord[:, 1]
        return self.coord

    def load_obs(self):
        print('Load observations')
        self.b = np.loadtxt(self.path2load+ self.obs)

    def load_sim(self):
        print('Load simulations')
        self.A = np.loadtxt(self.path2load+ self.sim)
        
    ### MAKE LINEAR SYSTEM 

    def check_nVRTe(self):
        if self.A.shape[0] / self.b.shape[0] ==  self.coord.shape[0]:
            self.nVRTe = self.coord.shape[0]
        else:
            raise ValueError('### dimensions of the files do not agree')

    def reshape_A(self):
        self.A = self.A.reshape((-1, self.nVRTe), order = 'F')

    def obs_w_f(self):
        """weight the observations, can also ignore observations by setting w = 0"""
        if self.obs_err == 'const':
            self.obs_w = np.ones(self.b.shape[0])
        elif self.obs_err == 'sqrt':
            self.obs_w = 1 / np.sqrt(np.abs(self.b))
            self.obs_w[self.obs_w >= self.errRmin] = 1    
    
    ### CONSTRAIN 

    def con_A_f(self):
        self.con_A = np.ones(self.A.shape[1])

    def con_b_f(self):
        self.con_b = np.ones(1)

    def con_w_f(self):
        self.con_w = np.ones(1) * self.wc
        
    ### REGULARIZATION 

    def nx_ny(self):
        """find number of nodes in each direction, has to be a regular grid"""
        self.nx = np.unique(np.round(self.coord[:, 0], 3)).shape[0]
        self.ny = np.unique(np.round(self.coord[:, 1], 3)).shape[0]
        print(self.nx, self.ny)

    def regularize_A(self):
        """create and append rows for spacial regularization to A"""
        reg = []
        vrte = range(1, self.nVRTe + 1)
        vrte = np.reshape(vrte,(self.ny, self.nx))
        for y in range(self.ny):
            for x in range(self.nx):
                minus = vrte[y, x]
                if x + 1 in range(self.nx):
                    plus = vrte[y , x + 1]
                    row = np.zeros(self.nVRTe, int)
                    row[minus -1] = - 1
                    row[plus -1] = + 1
                    reg.append(row)
                    if y + 1 in range(self.ny):
                        plus = vrte[y + 1, x]
                        row = np.zeros(self.nVRTe)
                        row[minus -1] = - 1
                        row[plus -1] = + 1
                        reg.append(row)
        self.reg_A = np.array(reg)
        print(len(self.reg_A))

    def regularize_b(self):
        """append 0's to b"""
        self.reg_b = np.zeros(self.reg_A.shape[0])

    def regularize_w(self):
        """create vector with weights, the length is determined by the number of regul rows in A"""
        self.reg_w = np.ones(self.reg_A.shape[0]) * self.wr

    ### VERTICAL STACK EQUATIONS

    def stack_A(self):
        self.A_s = np.vstack((self.A, self.con_A, self.reg_A))

    def stack_b(self):
        self.b_s = np.concatenate((self.b, self.con_b, self.reg_b))

    def stack_w(self):
        """create vector with weights for observation, constrain, and regularization
        then use it as diagonal for the weight matrix"""
        w = np.concatenate((self.obs_w, self.con_w, self.reg_w))
        W = np.zeros((w.shape[0], w.shape[0]))
        np.fill_diagonal(W, w)
        self.W_s = W

    ### APPLY WEIGHTS 

    def weight_A(self):
        """Apply the weights to A"""
        self.A_w = np.matmul(self.W_s, self.A_s)

    def weight_b(self):
        """Apply the weights to b"""
        self.b_w = np.matmul(self.b_s, self.W_s)

    ### PREPARE FOR ICSD

    def prepare4iCSD(self):
        """ this fucntion are called for each weight, keep them separated for pareto"""
        # create regularization part of the weight matrix
        self.regularize_w()
        self.stack_w()
        # apply weights with matrix multiplication
        self.weight_A()
        self.weight_b()

    ### LSQ 

    def iCSD(self):
        """solve linear system, given A matrix (VRTe, constrain, regul) and b (observations)"""
        self.x = lsq_linear(self.A_w, self.b_w, bounds = (0, 1))

    def ResidualAnalysis(self):
        fitting_res = self.x.fun[0 : self.b.shape[0]]
        # constrain_res = self.x.fun[self.b.shape[0] + 1] / self.wc
        regularization_res = self.x.fun[self.b.shape[0] + 2 :] / self.wr # constrain not included in the reg function
        self.reg_sum = np.sum(np.square(regularization_res)) # ||Ax - b||2
        self.fit_sum = np.sum(np.square(fitting_res)) # ||x - x0||2
        self.pareto_list_FitRes.append(self.fit_sum)
        self.pareto_list_RegRes.append(self.reg_sum)

    def run_pareto(self):
        """
        run iCSD multiple times while changing the weights to explore the L-curve
        """
        self.pareto_weights = np.linspace(self.pareto_MinErr, self.pareto_MaxErr, self.pareto_nSteps)
        print('pareto weights are\n', self.pareto_weights)

        self.pareto_list_FitRes = []
        self.pareto_list_RegRes = []
        with PdfPages('iCSD_pareto.pdf') as pdf:
            for self.wr in self.pareto_weights:
                print('-' * 80)
                print('Pareto mode, running new iCSDwith regularizaiton weight: ', self.wr)
                self.prepare4iCSD()
                self.iCSD()
                print(min(self.x.x), max(self.x.x))
                self.plotCSD()
                pdf.savefig(self.path2save + self.f)
                plt.close(self.f)

                self.ResidualAnalysis()

            self._plotPareto_()
            pdf.savefig(self.path2save+self.p)
            plt.close(self.p)
    ### PLOT 


    def _fig_RealSources_(self):
        """ add known point sources if present """
        if self.sc == None:
            return
        print('reading real sources: ', self.sc)
        for s in self.sc:
            sx = float(s.split(',')[0])
            sy = float(s.split(',')[1])
            plt.plot(sx, sy,'ow', markersize = 10, markeredgecolor = 'k')

    def _fig_ReturnElec_(self):
        """ plot the return electrode """
        if self.retElec == None:
            return
        print('reading return electrode: ', self.retElec)
        retElecx = float(self.retElec.split(',')[0])
        retElecy = float(self.retElec.split(',')[1])
        plt.plot(retElecx, retElecy,'sw', markersize = 10)

    def _fig_VRTe_(self):
        """ plot the VRTe current franctions """
        norm_z = (self.x.x - min(self.x.x)) / (max(self.x.x) - min(self.x.x))
        grey_cm = plt.cm.get_cmap('Greys')
        edgecolor_norm_z = grey_cm(norm_z)
        jet_cm = plt.cm.get_cmap('jet')
        facecolor_norm_z = jet_cm(norm_z)
        plt.scatter(self.coord_x, self.coord_y, facecolor = facecolor_norm_z, edgecolor = edgecolor_norm_z, cmap = 'jet')

    def mkGrid_XI_YI(self):
        """ grid for interpolation """
        Xm = np.linspace(min(self.coord_x), max(self.coord_x), 500)
        Ym = np.linspace(min(self.coord_y), max(self.coord_y), 500)
        self.XI, self.YI = np.meshgrid(Xm, Ym)

    def _fig_Interpolation_(self):
        """ plot the interpolation of the VRTe current fractions """
        points = np.column_stack((self.coord_x, self.coord_y))
        grid = griddata(points, self.x.x, (self.XI, self.YI), method = 'linear') # grid is a np array
        plt.imshow(grid,
        extent = (min (self.coord_x), max(self.coord_x), min(self.coord_y), max(self.coord_y)),
        aspect = 'auto', origin = 'lower', cmap= 'jet')
        cbar = plt.colorbar()
        # plt.clim(0,0.12)
        cbar.set_label('Fraction of Current Source', labelpad = 10)

    def _fig_Axis_Labels_(self):
        plt.ylabel('y [m]',fontsize=12)
        plt.xlabel('x [m]',fontsize=12)
        axes = plt.gca()
        axes.set_xlim([0,0.53])
        axes.set_ylim([0,0.52])
        plt.tick_params(axis='both', which='major')
        plt.tight_layout()

    def plotCSD(self):
        self.f = plt.figure('surface')
        self._fig_Interpolation_()
        self._fig_VRTe_()
        self._fig_RealSources_()
        self._fig_ReturnElec_()
        self._fig_Axis_Labels_()
        if not self.pareto:
            self._plotFIT_()

    def _plotPareto_(self):
        self.p = plt.figure('pareto', figsize=(10,4))
        plt.plot(self.pareto_list_FitRes, self.pareto_list_RegRes, 'or')
        ax = plt.gca()
        ax.tick_params(axis = 'both', which = 'both', direction = 'out')
        ax.grid()
        plt.xlabel('Residual')
        plt.ylabel('Roughness')
        plt.tight_layout()
        plt.savefig(self.path2save+'ParetoFront.png', dpi = 600)

    def _plotFIT_(self):
        plt.figure()
        plt.subplot(121)
        plt.plot(self.x.fun[:204] + self.b_w[:204], 'or', label = 'Inverted CSD')
        plt.plot(self.b_w[:204], 'ob', label = 'True model')
        plt.xlabel('Measurement number')
        plt.ylabel('R [Ohm]')
        plt.legend()
        plt.subplot(122)
        plt.plot(self.x.fun[:204] + self.b_w[:204], self.b_w[:204], 'or')
        plt.xlabel('Inverted CSD, R [Ohm]')
        plt.ylabel('True model, R [Ohm]')
        plt.tight_layout()
        plt.savefig(self.path2save+'Fit.png', dpi = 600)
        plt.show()

    def writeFIT(self):
        np.savetxt('inverted.txt', self.x.fun[:204] + self.b_w[:204])
        np.savetxt('truemodel.txt', self.b_w[:204])
        
        
    def run_single(self):
        self.prepare4iCSD()
        self.iCSD()
        self.plotCSD()
        self.writeFIT()
        self.f.savefig(self.path2save+ 'iCSD', dpi = 600)
        plt.show()
        