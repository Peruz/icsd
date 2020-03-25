import re
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.interpolate import griddata
from scipy.linalg import lu
from scipy.optimize import lsq_linear
from matplotlib.backends.backend_pdf import PdfPages

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 12}
matplotlib.rc('font', **font)

class iCSD_class():
    def get_args(self):
        ap = argparse.ArgumentParser(description = ' current source density inversion')
        ap.add_argument('-sim', metavar = 'BERT_VRTe_calculation', type = str,
            help = 'file containing the forward calculation for each virtual source, going into the A matrix',
            default = 'VRTeSim.txt')
        ap.add_argument('-obs', metavar = 'Observed_MALM_data', type = str,
            help = 'measured MALM data, b column of the linear system',
            default = 'ObsData.txt')
        ap.add_argument('-coord_file', metavar = 'VRTe_coordinates', type = str,
            help = 'coordinates of the VRTe for plotting',
            default = 'VRTeCoord.txt')
        ap.add_argument('-wr', metavar = 'weight regularization', type = float,
            help = 'first derivative regularization weight', default = 25, nargs = '*')
        ap.add_argument('-wc', metavar = 'weight constrain', type = float,
            help = 'current conservation constrain, sum current fractions = 1', default = 10000)
        ap.add_argument('-sc', metavar = 'source coords', type = str,
            help = 'coordinates of the sources, format = x1,y1 x2,y2', nargs = '*')
        ap.add_argument('-retElec', metavar = 'Return Elec', type = str,
            help = 'coordinates of the return electrode, format = x1,y1') # rhi elec 64 = 0.44,0.8
        ap.add_argument('-pareto', action = 'store_true',
            help = 'if True run many icsd to explore the Pareto front', default = False)
        ap.add_argument('-errRmin', metavar = 'min R in err', type = float,
            help = 'min R to which the err is applied before passing to constant error', default = 1)
        ap.add_argument('-pareto_MinErr', metavar = 'min Pareto error', type = float,
            help = 'starting error for pareto mode', default = 1)
        ap.add_argument('-pareto_MaxErr', metavar = 'max Pareto error', type = float,
            help = 'maximum error for pareto mode', default = 500)
        ap.add_argument('-pareto_nSteps', metavar = 'number steps Pareto', type = int,
            help = 'number of error steps for the pareto mode', default = 21)
        ap.add_argument('-obs_err', metavar = 'obs weight mode', type = str,
            help = 'choose between constant weight and w = 1/sqrt(abs(obs))', default = 'const',
            choices = ['const', 'sqrt'])
        self.args = ap.parse_args()

    ### PLOT 

    def _fig_RealSources_(self):
        """ add known point sources if present """
        if self.args.sc == None:
            return
        print('reading real sources: ', self.args.sc)
        for s in self.args.sc:
            sx = float(s.split(',')[0])
            sy = float(s.split(',')[1])
            plt.plot(sx, sy,'ow', markersize = 10, markeredgecolor = 'k')

    def _fig_ReturnElec_(self):
        """ plot the return electrode """
        if self.args.retElec == None:
            return
        print('reading return electrode: ', self.args.retElec)
        retElecx = float(self.args.retElec.split(',')[0])
        retElecy = float(self.args.retElec.split(',')[1])
        plt.plot(retElecx, retElecy,'sw', markersize = 10)

    def _fig_VRTe_(self):
        """ plot the VRTe current franctions """
        norm_z = (self.x.x - min(self.x.x)) / (max(self.x.x) - min(self.x.x))
        grey_cm = plt.cm.get_cmap('Greys')
        edgecolor_norm_z = grey_cm(norm_z)
        jet_cm = plt.cm.get_cmap('jet')
        facecolor_norm_z = jet_cm(norm_z)
        plt.scatter(self.coord_x, self.coord_y, facecolor = facecolor_norm_z, edgecolor = edgecolor_norm_z, cmap = 'jet')

    def _mkGrid_XI_YI_(self):
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

    def _plotCSD_(self):
        self.f = plt.figure('surface')
        self._fig_Interpolation_()
        self._fig_VRTe_()
        self._fig_RealSources_()
        self._fig_ReturnElec_()
        self._fig_Axis_Labels_()
        if not self.args.pareto:
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
        plt.savefig('ParetoFront.png', dpi = 600)

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
        plt.savefig('Fit.png', dpi = 600)
        plt.show()

    def _writeFIT_(self):
        np.savetxt('inverted.txt', self.x.fun[:204] + self.b_w[:204])
        np.savetxt('truemodel.txt', self.b_w[:204])
    ### LOAD 

    def _load_coord_(self):
        self.coord = np.loadtxt(self.args.coord_file)
        self.coord_x, self.coord_y = self.coord[:, 0], self.coord[:, 1]

    def _load_obs_(self):
        self.b = np.loadtxt(self.args.obs)

    def _load_sim_(self):
        self.A = np.loadtxt(self.args.sim)

    ### MAKE LINEAR SYSTEM 

    def _check_nVRTe_(self):
        if self.A.shape[0] / self.b.shape[0] ==  self.coord.shape[0]:
            self.nVRTe = self.coord.shape[0]
        else:
            raise ValueError('### dimensions of the files do not agree')

    def _reshape_A_(self):
        self.A = self.A.reshape((-1, self.nVRTe), order = 'F')

    def _obs_w_(self):
        """weight the observations, can also ignore observations by setting w = 0"""
        if self.args.obs_err == 'const':
            self.obs_w = np.ones(self.b.shape[0])
        elif self.args.obs_err == 'sqrt':
            self.obs_w = 1 / np.sqrt(np.abs(self.b))
            self.obs_w[self.obs_w >= self.args.errRmin] = 1

    ### CONSTRAIN 

    def _con_A_(self):
        self.con_A = np.ones(self.A.shape[1])

    def _con_b_(self):
        self.con_b = np.ones(1)

    def _con_w_(self):
        self.con_w = np.ones(1) * self.args.wc

    ### REGULARIZATION 

    def _nx_ny_(self):
        """find number of nodes in each direction, has to be a regular grid"""
        self.nx = np.unique(np.round(self.coord[:, 0], 3)).shape[0]
        self.ny = np.unique(np.round(self.coord[:, 1], 3)).shape[0]
        print(self.nx, self.ny)

    def _regularize_A_(self):
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

    def _regularize_b_(self):
        """append 0's to b"""
        self.reg_b = np.zeros(self.reg_A.shape[0])

    def _regularize_w_(self):
        """create vector with weights, the length is determined by the number of regul rows in A"""
        self.reg_w = np.ones(self.reg_A.shape[0]) * self.args.wr

    ### VERTICAL STACK EQUATIONS

    def _stack_A_(self):
        self.A_s = np.vstack((self.A, self.con_A, self.reg_A))

    def _stack_b_(self):
        self.b_s = np.concatenate((self.b, self.con_b, self.reg_b))

    def _stack_w_(self):
        """create vector with weights for observation, constrain, and regularization
        then use it as diagonal for the weight matrix"""
        w = np.concatenate((self.obs_w, self.con_w, self.reg_w))
        W = np.zeros((w.shape[0], w.shape[0]))
        np.fill_diagonal(W, w)
        self.W_s = W

    ### APPLY WEIGHTS 

    def _weight_A_(self):
        """Apply the weights to A"""
        self.A_w = np.matmul(self.W_s, self.A_s)

    def _weight_b_(self):
        """Apply the weights to b"""
        self.b_w = np.matmul(self.b_s, self.W_s)

    ### PREPARE FOR ICSD

    def _prepare4iCSD_(self):
        """ this fucntion are called for each weight, keep them separated for pareto"""
        # create regularization part of the weight matrix
        self._regularize_w_()
        self._stack_w_()
        # apply weights with matrix multiplication
        self._weight_A_()
        self._weight_b_()

    ### LSQ 

    def _iCSD_(self):
        """solve linear system, given A matrix (VRTe, constrain, regul) and b (observations)"""
        self.x = lsq_linear(self.A_w, self.b_w, bounds = (0, 1))

    def _ResidualAnalysis_(self):
        fitting_res = self.x.fun[0 : self.b.shape[0]]
        # constrain_res = self.x.fun[self.b.shape[0] + 1] / self.args.wc
        regularization_res = self.x.fun[self.b.shape[0] + 2 :] / self.args.wr # constrain not included in the reg function
        self.reg_sum = np.sum(np.square(regularization_res)) # ||Ax - b||2
        self.fit_sum = np.sum(np.square(fitting_res)) # ||x - x0||2
        self.pareto_list_FitRes.append(self.fit_sum)
        self.pareto_list_RegRes.append(self.reg_sum)

    ### PUBLIC METHODS

    def icsd_init(self):
        """ these functions are called only once, even for pareto,
        as they do not depend on the regularization weight"""
        # load
        self._load_coord_()
        self._load_obs_()
        self._load_sim_()
        # check vector sizes
        self._check_nVRTe_()
        # reshape VRTe vector into matrix A 
        self._reshape_A_()
        # define mode to weights the data        
        self._obs_w_()
        # constrain (curent conservation)
        self._con_A_()
        self._con_b_()
        self._con_w_()
        # get VRTe grid size and make grid
        self._nx_ny_()
        self._mkGrid_XI_YI_()
        # append spatial regularization
        self._regularize_A_()
        self._regularize_b_()
        # stack data, constrain, and regularization 
        self._stack_A_()
        self._stack_b_()

    def run_pareto(self):
        """
        run iCSD multiple times while changing the weights to explore the L-curve
        """
        self.pareto_weights = np.linspace(self.args.pareto_MinErr, self.args.pareto_MaxErr, self.args.pareto_nSteps)
        print('pareto weights are\n', self.pareto_weights)

        self.pareto_list_FitRes = []
        self.pareto_list_RegRes = []
        with PdfPages('iCSD_pareto.pdf') as pdf:
            for self.args.wr in self.pareto_weights:
                print('-' * 80)
                print('Pareto mode, running new iCSDwith regularizaiton weight: ', self.args.wr)
                self._prepare4iCSD_()
                self._iCSD_()
                print(min(self.x.x), max(self.x.x))
                self._plotCSD_()
                pdf.savefig(self.f)
                plt.close(self.f)

                self._ResidualAnalysis_()

            self._plotPareto_()
            pdf.savefig(self.p)
            plt.close(self.p)

    def run_single(self):
        self._prepare4iCSD_()
        self._iCSD_()
        self._plotCSD_()
        self._writeFIT_()
        self.f.savefig('iCSD', dpi = 600)
        plt.show()

if __name__ == "__main__":
    icsd = iCSD_class()
    icsd.get_args()
    icsd.icsd_init()
    if icsd.args.pareto == True:
        print('run pareto')
        icsd.run_pareto()
    elif icsd.args.pareto == False:
        print('run single inversion')
        icsd.run_single()
