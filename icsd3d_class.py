import os
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.interpolate import griddata
from scipy.linalg import lu
from scipy.optimize import lsq_linear, curve_fit, least_squares, leastsq
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm
from matplotlib.ticker import MultipleLocator
from kneed import KneeLocator
from scipy.interpolate import griddata as gd
from mpl_toolkits.mplot3d import Axes3D
import pyvista as pv
from scipy.sparse import diags

import sys
from numpy import linalg as LA

class iCSD3d_Class():
    """ Create a icsd inversion object.
    
    Parameters
    ----------
    coord_file : str, mandatory
        coordinates of the VRTe for plotting
    wr : float, optional
        Weight regularization
    wc : float, optional
    """
    def __init__(self,dirName):
        self.dirName = dirName
        self.clim = []
       # load
        self.type='3d'
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
        self.obs_err='const' #const or sqrt - choose between constant weight and w = 1/sqrt(abs(obs))
        # IMPLEMENT obs_err based on reciprocal analysis i.e. estimated standard deviation of the data errors;                         % estimated standard deviation of the traveltime data errors
        self.k=4  # For each point, find the k closest current sources
        self.TL=False # Time lapse inversion
        self.x0_prior=False #  relative smallness regularization as a prior criterion for the inversion; i.ethe algorithm minimizes ||m−m0||2
        self.x0_ini_guess=False # initial guess
        self.knee=False # L-curve knee automatic detection
        self.KneeWr=[]
        self.regMesh='strc' # strc or unstrc
        self.plotElecs=False # strc or unstrc
        self.alphax0=1 # weight on model smallness relative to m0
        self.alphaSxy=False # weight on model smoothness in z-direction [TO IMPLEMENT]
        self.alphaSx=1 # weight on model smoothness in x-direction [TO IMPLEMENT]
        self.alphaSy=1 # weight on model smoothness in y-direction [TO IMPLEMENT]
        self.alphaSz=1 # weight on model smoothness in z-direction [TO IMPLEMENT]
        self.mesh=None # weight on model smoothness in z-direction [TO IMPLEMENT]

    def icsd_init(self):
        """ these functions are called only once, even for pareto,
        as they do not depend on the regularization weight"""
        # load
        self.createdirs()
        self.load_coord()
        self.load_obs()
        self.load_sim()
        if self.plotElecs==True:
            self.load_geom() # geometry file containing electrodes position includinf remotes 
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
        if self.type=='2d':
            self.nx_ny()
        else:
            self.nx_ny_nz()
        self.mkGrid_XI_YI()
        # # append spatial regularization
        if self.type=='2d': # 2D CASE -----------------------------------------
            if self.regMesh=='strc':
                if self.alphaSxy==True:
                    self.regularize_A_x_y()
                    if self.x0_prior==True:
                        self.regularize_smallnessX0() #add smallness regularisation
                    self.regularize_sum_AX0()
                else:
                    if self.x0_prior==True:
                        raise ValueError('#### dimensions of matrices do not agree - change regularisation types')
                    else:
                        self.regularize_A()
            else:
                self.regularize_A_UnstructuredMesh3d()
        else:              # 3D CASE -----------------------------------------
            if self.regMesh=='strc':
                if self.x0_prior==True:
                        raise ValueError('### dimensions of matrices do not agree - change regularisation type')
                else:
                    self.regularize_A_3d() # working for structured mesh
            elif self.regMesh=='unstrc':
                self.regularize_A_UnstructuredMesh3d() # working for unstructured mesh
                if self.x0_prior==True:
                    self.regularize_smallnessX0() #add smallness regularisation
            elif self.alphaSxy==True:
                if self.x0_prior==True:
                    self.regularize_smallnessX0() #add smallness regularisation
                self.regularize_sum_AX0()
        self.regularize_b()
        # # stack data, constrain, and regularization 
        self.stack_A()
        self.stack_b()
        
        
# NEW FUNCTIONS to introduce 
# function penalized some virtual sources (split rhizotron)
# load_names for tl analysis
# Introduce error weigting using reciprocal error instead of constant or w = 1/sqrt(abs(obs))
# use the ERT mesh as mesh for virtual position of current sources
        
# NEW FUNCTIONS COMPARE TO the 2d case    

    def load_geom(self):
        geom_files = [f for f in os.listdir(self.path2load) if f.endswith('.geom')]
        if len(geom_files) != 1:
            raise ValueError('should be only one geom file in the current directory')

        fileNameElec = geom_files[0]
        #fileNameElec= '*.geom'
        # RemLineNb= search_string_in_file(fileNameElec, '#Remote')

        line_number = 0
        line_of_results = []
        # Open the file in read only mode
        with open(self.path2load + fileNameElec, 'r') as read_obj:
            # Read all lines in the file one by one
            for line in read_obj:
                # For each line, check if line contains the string
                line_number += 1
                if ('#Remote') in line:
                    # If yes, then add the line number & line as a tuple in the list
                    line_of_results.append((line_number))
        self.RemLineNb= int(line_of_results[0])-1

        self.coordE = np.loadtxt(self.path2load+ fileNameElec)
        self.pointsE= np.vstack(self.coordE[:self.RemLineNb,1:4])

    def plotCSD3d_pyvista(self):

        filename = self.path2load + 'ExportSol.dat'
        data_2_plot =np.genfromtxt(filename)
        coord_x= data_2_plot[:,0]
        coord_y= data_2_plot[:,1]
        coord_z= data_2_plot[:,2]
        step=(max(coord_x)-min(coord_x))/10
        num =10
        opacity = [0, 0, 0.1, 0.3, 0.6, 0.9, 1]

        grid = pv.UniformGrid()
        spc=(max(coord_x)-min(coord_x))/10
        # spc=1
        xdim = int(round((max(coord_x)-min(coord_x))/spc))
        ydim = int(round((max(coord_y)-min(coord_y))/spc))
        zdim = int(round((max(coord_z)-min(coord_z))/spc))
        grid.dimensions = (xdim, ydim, zdim)
        grid.dimensions = np.array(grid.dimensions) +1
        grid.origin = (min(coord_x), min(coord_y), min(coord_z)) # The bottom left corner of the data set
        grid.spacing = (spc, spc,spc) # These are the cell sizes along each axis
        

        coord= data_2_plot[:,:-1]
        pv.set_plot_theme('document')
        poly = pv.PolyData(coord)
        pvfig = pv.Plotter(notebook=False,window_size=[600, 600])
        
        if self.mesh!=None:
            print('plot mesh.vtk')
            ModelVtk = pv.read(self.path2load + self.mesh)
            cmap = plt.cm.get_cmap('viridis', 2)
            pvfig.add_bounding_box()
            pvfig.add_mesh(cmap=cmap,mesh=ModelVtk,scalars='rhomap', opacity=0.2)    # add a dataset to the scene
        pvfig.add_mesh(poly, point_size=15.0, scalars=data_2_plot[:,3], opacity=opacity, render_points_as_spheres=True,cmap='jet')
        print('interpolation spacing='+ str(spc))
        interpolated = grid.interpolate(poly, radius=spc)
        cmap = plt.cm.get_cmap('jet',10)
        contours = interpolated.contour()
        # pvfig.add_mesh(interpolated, show_scalar_bar=False, cmap=cmap,opacity=0.3)
        pvfig.add_mesh(contours, show_scalar_bar=False, opacity= opacity,cmap='jet')
        pvfig.show_bounds(bounds=[min(coord_x), max(coord_x), 
                                        min(coord_y),max(coord_y),
                                        min(coord_z), 0],font_size=16)
        pvfig.add_axes()
        pvfig.show_axes()
        pvfig.add_scalar_bar('Normalized Current density',width=0.25,vertical=False,position_x=0.3)
        # pvfig.update_scalar_bar_range([minTs, cbarmax])  
        if self.plotElecs==True:
            pvfig.add_points(self.pointsE)
            # p.add_point_labels(pointsE,coordE[:RemLineNb,0].astype(int), point_size=15, font_size=35,
            #     shape_opacity=0.01, margin=4.)
            pvfig.add_point_labels(self.pointsE,self.coordE[:self.RemLineNb,0].astype(int), point_size=15, font_size=15,
                shape_opacity=0.01, margin=4.)
        pvfig.show(auto_close=True)  
        if self.knee==True:
           if self.wr==self.KneeWr:
               # pvfig.screenshot('Pts_iCSD_knee'+ str(self.ObsName) + '.png')
                pvfig.screenshot(self.path2save+ 'Pts_iCSD_knee_wr'+ self.obs +  str(self.KneeWr) + '.png')
        else:
                pvfig.screenshot(self.path2save+ 'Pts_iCSD_wr'+ str(self.wr) + '.png')

        
    def plotCSD3d(self):
        self.f = plt.figure('volume')

        step=(max(self.coord_x)-min(self.coord_x))/10

        xlin=np.arange(min(self.coord_x),max(self.coord_x),step)
        ylin=np.arange(min(self.coord_y),max(self.coord_y),step)
        zlin=np.arange(min(self.coord_z),max(self.coord_z),step)
        #generate new grid
        X,Y,Z=np.meshgrid(xlin,ylin,zlin)
        data_2_plot= self.x.x
        #interpolate "data.v" on new grid "inter_mesh"
        # V = gd((self.coord_x,self.coord_y,self.coord_z), data_2_plot, (X,Y,Z), method='linear')
        
        ax=self.f.gca(projection='3d')
        sc=ax.scatter(self.coord_x, self.coord_y, self.coord_z, c=data_2_plot, cmap ='coolwarm', s=data_2_plot*1e4,
                 )
        cbar = plt.colorbar(sc)
        cbar.set_label('# current density')
        # ax.view_init(azim=-101, elev=35)
        title= 'Scattered current sources density' 
        plt.title(title)
        plt.savefig(self.path2save+ 'icsd_scatter', dpi=550,bbox_inches='tight',pad_inches = 0)

        if self.knee==True:
           if self.wr==self.KneeWr:
               plt.savefig(self.path2save+ 'icsd_knee_scatter'+ str(self.KneeWr) + '.png',dpi=550,bbox_inches='tight',pad_inches = 0)
        plt.show()

      
    def DetectKneePt(self):
        self.kn = KneeLocator(self.pareto_list_FitRes,self.pareto_list_RegRes, 
                         curve='convex', direction='decreasing')
#        print('knee xloc=' + str(self.kn.knee))
        self.IdPtkneew= np.where(self.kn.knee==self.pareto_list_FitRes)[0]
#        print('Id pt knee=' + str(self.IdPtkneew))
        self.pareto_weights[self.IdPtkneew]
#        print('Wr knee=' + str(self.pareto_weights[self.IdPtkneew]))
        if len(self.IdPtkneew)<1:
             self.IdPtkneew=1
             print('No knee detection possible, put 1 as default')
#            raise ValueError('No knee detection possible')
             
    ### Individual Misfit  
    def normF1(self):
        self.F1=[]
        for i in range(np.shape(self.A)[1]):
            F1i = LA.norm((self.b-self.A[:,i]))
            self.F1.append(F1i)
        self.norm_F1 = (self.F1 - min(self.F1)) / (max(self.F1) - min(self.F1))

    def misfit_2_initialX0(self):
        self.x0F1=1./((self.norm_F1+1)*(self.norm_F1+1))
        self.x0F1_sum= self.x0F1/sum(self.x0F1) # normlize such as sum equal to 1
        # self.x0F1_sum=np.ones(self.x0F1_sum.shape)*0.1


    def plotmisfitF1(self):
        fig, axs = plt.subplots(nrows=1, ncols=2, sharex=True)
        ax = axs[0]
        if self.type=='2d':
            points = np.column_stack((self.coord_x, self.coord_y))
            grid = griddata(points, self.norm_F1, (self.XI, self.YI), method = 'linear') # grid is a np array
            im1 = ax.imshow(grid,norm=LogNorm(vmin=0.3, vmax=0.7), extent = (min (self.coord_x), max(self.coord_x), min(self.coord_y), max(self.coord_y)),
            aspect = 'auto', origin = 'lower', cmap= 'jet')
            ax.set_ylabel('y [m]',fontsize=15)
            ax.set_xlabel('x [m]',fontsize=15)
            ax.set_title('F1 misfit',fontsize=15)

            #axes = plt.gca()
            #axes.set_xlim([0,0.53])
            #axes.set_ylim([0,0.52])
            ax.tick_params(axis='both', which='major', labelsize=15)
            #ax.set_tight_layout()
            ax.set_aspect(1.0)
            cbar1 = plt.colorbar(im1,ax=ax, format="%.2f",fraction=0.046, pad=0.04)
            cbar1.set_label('Normalised misfit', labelpad = 5, fontsize=14)
    
            ax = axs[1]
            grid = griddata(points, self.x0F1, (self.XI, self.YI), method = 'linear') # grid is a np array
            im2 = ax.imshow(grid,norm=LogNorm(vmin=0.3, vmax=0.7), extent = (min (self.coord_x), max(self.coord_x), min(self.coord_y), max(self.coord_y)),
            aspect = 'auto', origin = 'lower', cmap= 'jet')
            ax.set_ylabel('y [m]',fontsize=15)
            ax.set_xlabel('x [m]',fontsize=15)
            ax.set_title('x0 solution',fontsize=15)
            #axes = plt.gca()
            #axes.set_xlim([0,0.53])
            #axes.set_ylim([0,0.52])
            ax.tick_params(axis='both', which='major', labelsize=15)
            ax.set_aspect(1.0)
            cbar2 = plt.colorbar(im2,ax=ax, format="%.3f",fraction=0.046, pad=0.04)
            cbar2.set_label('x0 solution', labelpad = 5, fontsize=14)
        else:
            print('3d case to write')
            
        # fig.suptitle(self.ObsName, y=0.80)
        fig.tight_layout()
        # fig.savefig('F1_x0F1'+ self.ObsName, dpi = 450, bbox_inches='tight')


    def run_misfitF1(self):
        self.normF1()
        self.misfit_2_initialX0()
        self.plotmisfitF1()
        
    def regularize_A_UnstructuredMesh3d(self): # should also work for the 2d case
        print('regularize_A_UnstructuredMesh3d')
        reg = []
        for VRTEnb in range(self.nVRTe):
            dist =  np.linalg.norm(self.coord[VRTEnb]-self.coord, axis=1)
            closest = np.argsort(dist)
            k = self.k  # For each point, find the k closest current sources
            #print('k=' + str(k))
            Ind = closest[1:k+1]
            row = np.zeros(self.nVRTe) # add a line to the regularisation A with k non-null coefficients
            #row[Ind]= -1/k # ponderate coeff for k closest current sources
            knorm = dist[closest[1:k+1]]/dist[closest[1:k+1]].sum(axis=0,keepdims=1)
            row[Ind]= -knorm
            row[VRTEnb]= 1 # = one for the actual current source
            reg.append(row)
            test=[1]
            mask = np.in1d(test, VRTEnb)
            #print(mask)
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
        self.reg_A = np.array(reg)

        
    def regularize_A_3d(self):
        """create and append rows for spacial regularization to A"""
        reg = []
        vrte = range(1, self.nVRTe + 1)
        vrte = np.reshape(vrte,(self.ny, self.nx, self.nz))
        for z in range(self.nz):
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

    def nx_ny_nz(self):
        """find number of nodes in each direction, has to be a regular grid"""
        self.nx = np.unique(np.round(self.coord[:, 0], 3)).shape[0]
        self.ny = np.unique(np.round(self.coord[:, 1], 3)).shape[0]
        self.nz = np.unique(np.round(self.coord[:, 2], 3)).shape[0]

    def Export_sol(self):
        fileName= self.path2load+ 'ExportSol.dat'
        if self.type=='2d':
            ExpSOL = np.vstack([self.coord_x,self.coord_y,self.x.x])
            ExpSOL= ExpSOL.T
            f = open(fileName,'w')
            np.savetxt(f, ExpSOL, fmt='%1.2f %1.2f %1.6f', delimiter='\t',header='X Y i')   # X is an array
            f.close()
        else:
            ExpSOL = np.vstack([self.coord_x,self.coord_y,self.coord_z,self.x.x])
            ExpSOL= ExpSOL.T
            f = open(fileName,'w')
            np.savetxt(f, ExpSOL, fmt='%1.2f %1.2f %1.2f %1.6f', delimiter='\t',header='X Y Z i')   # X is an array
            f.close()

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
        if self.type=='2d':
            self.coord_x, self.coord_y =  self.coord[:, 0],  self.coord[:, 1]
        else:
            self.coord_x, self.coord_y, self.coord_z = self.coord[:, 0], self.coord[:, 1], self.coord[:, 2]
        return self.coord

    def load_obs(self):
        print('Load observations')
        self.b = np.loadtxt(self.path2load+ self.obs)

    def load_sim(self):
        print('Load simulations')
        self.A = np.loadtxt(self.path2load+ self.sim)
        print('*'*36)

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
        elif self.obs_err == 'reciprocals': #[TO IMPLEMENT]
            self.obs_w = 1 / np.sqrt(np.abs(self.sd_rec))

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

    def regularize_A(self):
        """create and append rows for spatial regularization to A"""
        print('Reg A (Luca''s implementation)')
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
        
    def regularize_A_x_y(self):
        """create and append rows for spatial regularization to A"""
        print('regularize_A_x_y')
        ncells = self.nx*self.ny;       
        Dx=diags([1, -2, 1], [-1, 0, 1], shape=(ncells, ncells)).todense()
        idx=np.arange(0,len(Dx)-1,self.nx)
        Dx[idx,:] = 0;
        idx2=np.arange(self.nx,len(Dx)-1,self.nx)
        Dx[idx2,:] = 0;
        
        Dy=diags([1, -2, 1], [-self.nx, 0, self.nx], shape=(ncells, ncells)).todense()
        idy=np.arange(0,self.nx-1)
        Dy[idy,:] = 0;
        idy2=np.arange(((self.ny-2)*self.nx),(self.nx-1)*(self.ny-1))
        Dy[idy2,:] = 0;
        
        self.reg_Ax = self.alphaSx*np.array(Dx.transpose()*Dx)
        self.reg_Ay = self.alphaSy*np.array(Dy.transpose()*Dy)

    ### RELATIVE SMALLNESS conditions (m-m0)
        
    def regularize_smallnessX0(self):
        print('reg smallnessX0')
        if self.alphaSxy==True:
            self.reg_smallx0 = np.ones(self.reg_Ax.shape)*self.alphax0
        else:
            self.reg_smallx0 = np.ones(self.reg_A.shape)*self.alphax0
        # self.reg_smallx0 = np.ones(self.x0F1_sum.shape)*self.alphax0
        
    def regularize_sum_AX0(self):
        """sum smallness and spatial regularisation"""
        if (self.alphaSxy==True and self.x0_prior==True):
            print("""sum small x0, reg Ax, reg Ay""") 
            self.reg_A= self.reg_smallx0 + self.reg_Ax + self.reg_Ay
        elif (self.alphaSxy==True and self.x0_prior==False):
            print("""sum reg Ax, reg Ay""") 
            self.reg_A= self.reg_Ax + self.reg_Ay
        elif (self.alphaSxy==False and self.x0_prior==True):
            print("""reg_A= reg_A + small x0""") 
            self.reg_A= self.reg_A + self.reg_smallx0
            
    def regularize_b(self):
        self.reg_b = np.zeros(self.reg_A.shape[0])

    def regularize_w(self):
        """create vector with weights, the length is determined by the number of regul rows in A"""
        self.reg_w = np.ones(self.reg_A.shape[0]) * self.wr
        if self.x0_prior==True:
            print('reg Wm (smallness + spatial reg) * lambda=' + str(self.wr))
            self.reg_w_0_b = np.ones(self.reg_A.shape[0]) * self.x0F1_sum * self.wr
            self.reg_w_0_A = np.ones(self.reg_A.shape[0])* self.wr
            
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
        
        if self.x0_prior==True:
            #add lines for relative smallness 
            wa = np.concatenate((self.obs_w, self.con_w, self.reg_w_0_A))
            wb = np.concatenate((self.obs_w, self.con_w, self.reg_w_0_b))
            W = np.zeros((wa.shape[0], wa.shape[0]))
            np.fill_diagonal(W, wa)
            self.W_s_A = W
            np.fill_diagonal(W, wb)
            self.W_s_b = W

    ### APPLY WEIGHTS 

    def weight_A(self):
        """Apply the weights to A"""
        if self.x0_prior==True:
            self.A_w = np.matmul(self.W_s_A, self.A_s)
        else:
            self.A_w = np.matmul(self.W_s, self.A_s)
            
    def weight_b(self):
        """Apply the weights to b"""
        if self.x0_prior==True:
            self.b_w = np.matmul(self.b_s,  self.W_s_b)
        else:
            self.b_w = np.matmul(self.b_s,  self.W_s)
    
    ### PREPARE FOR ICSD

    def prepare4iCSD(self):
        """ this fucntion are called for each weight, keep them separated for pareto"""
        # create regularization part of the weight matrix
        if self.x0_ini_guess==True:
            self.run_misfitF1()
        elif self.x0_prior==True:
            self.run_misfitF1()
        self.regularize_w()
        self.stack_w()
        # apply weights with matrix multiplication
        self.weight_A()
        self.weight_b()

    ### LSQ 

    def iCSD(self):
        """solve linear system, given A matrix (VRTe, constrain, regul) and b (observations)"""
        if self.x0_ini_guess==False:
            print('No initial guess')
            self.x = lsq_linear(self.A_w, self.b_w, bounds = (0, 1))
            print('*' * 20)
            print('CURRENT Sum=' + str(np.sum(self.x.x)))
        else:
            print('Initial guess x0')
            a = self.A_w 
            b = self.b_w
            def func(x, a, b):
                return (b - np.dot(a, x))
            self.x = least_squares(func, x0=self.x0F1_sum, bounds = (0, 1), args=(a, b)) # Add initial guess
            print('CURRENT Sum=' + str(np.sum(self.x.x)))
        self.Export_sol()
        
        
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
        with PdfPages(self.path2save+ 'iCSD_pareto.pdf') as pdf:
            for self.wr in self.pareto_weights:              
                if self.wr== self.pareto_weights[0]:
                    if self.x0_prior==True:
                        self.run_misfitF1()
                    if self.x0_ini_guess==True:
                        print('initial x0 guess')
                        self.run_misfitF1()
            
                print('-' * 80)
                print('Pareto mode, running new iCSDwith regularizaiton weight: ', self.wr)
                self.prepare4iCSD()
                self.iCSD()
                if self.type=="2d":
                    self.plotCSD()
                else:
                    self.plotCSD3d()
                    # self.plotCSD3d_pyvista()

                pdf.savefig(self.f)
                plt.close(self.f)

                self.ResidualAnalysis()

            
            self.DetectKneePt()
            self.wr=float(self.pareto_weights[self.IdPtkneew])
            print('Knee detected for wr=' + str(self.wr))
            self._plotPareto_()
            pdf.savefig(self.p)
            plt.close(self.p)

            if self.knee==True:
               self.plot_knee_icsd()
               
    ### PLOT
    def plot_knee_icsd(self):
        self.KneeWr=self.wr
        # self.wr=float(self.pareto_weights[self.IdPtkneew])
        self.kn.plot_knee_normalized()
        # self.kn.savefig('Lc', dpi = 450)
        # pdf.savefig(self.p)
        # plt.close(self.p)

        self.run_single()
        #self.f.savefig('iK'+ str("{:02d}".format(int(ObsName.translate({ord(i): None for i in 'OW'})))), dpi = 450,
        #    bbox_inches='tight',pad_inches = 0)
        plt.show()
            

    def _fig_RealSources_(self):
        """ add known point sources if present """
        if self.sc == None:
            return
        # print('reading real sources: ', self.sc)
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
        if self.type=='2d':
            points = np.column_stack((self.coord_x, self.coord_y))
            grid = griddata(points, self.x.x, (self.XI, self.YI), method = 'linear') # grid is a np array
            plt.imshow(grid,
                       extent = (min (self.coord_x), max(self.coord_x), min(self.coord_y), max(self.coord_y)),
                       aspect = 'auto', origin = 'lower', cmap= 'jet')
            cbar = plt.colorbar()
            # plt.clim(0,np.percentile(self.x.x, 80))
            if self.clim:
                plt.clim(self.clim[0],self.clim[1])
            cbar.set_label('Fraction of Current Source', labelpad = 10)
        # else:
        #     print('3d case to write')   
        if self.mesh!=None:
            print('to write')
            # print('plot mesh.vtk')
            # ModelVtk = pv.read(self.path2load + self.mesh)
            # mesh=pg.load(path2files+ icsd3d_TL_RWU.mesh)
            # pg.show(mesh,data=rhomap,label='rhomap')
            
    def _fig_Axis_Labels_(self):
        plt.ylabel('y [m]',fontsize=12)
        plt.xlabel('x [m]',fontsize=12)
        axes = plt.gca()
        # axes.set_xlim([0,0.53])
        # axes.set_ylim([0,0.52])
        plt.tick_params(axis='both', which='major')
        plt.tight_layout()
        axes.set_aspect('equal')
        
    def plotCSD(self):
        self.f = plt.figure('surface')
        if self.mesh!=None:
            self.f, ax = plt.subplots('surface',nrows=2)

        self._fig_Interpolation_()
        self._fig_VRTe_()
        self._fig_RealSources_()
        self._fig_ReturnElec_()
        self._fig_Axis_Labels_()
        if not self.pareto:
            self._plotFIT_()

    def _plotPareto_(self):
        self.p, self.ax = plt.subplots()
        # self.p = plt.figure('pareto', figsize=(10,4))
        self.ax.annotate('Wr=' + str(int(self.wr)), xy=(float(np.asarray(self.pareto_list_FitRes)[self.IdPtkneew]), 
                                         float(np.asarray(self.pareto_list_RegRes)[self.IdPtkneew])), 
                                      xytext=(float(np.asarray(self.pareto_list_FitRes)[self.IdPtkneew])+max(self.pareto_list_FitRes)/3, 
                                         float(np.asarray(self.pareto_list_RegRes)[self.IdPtkneew])+max(self.pareto_list_RegRes)/3),
                                      arrowprops=dict(facecolor='black', shrink=0.05))
        plt.plot(float(np.asarray(self.pareto_list_FitRes)[self.IdPtkneew]), float(np.asarray(self.pareto_list_RegRes)[self.IdPtkneew]), 'og')
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
        print('write FIT')
        np.savetxt(self.path2load+'invertedR.txt', self.x.fun[:] + self.b_w[:])
        np.savetxt(self.path2load+'truemodelR.txt', self.b_w[:])
        # np.savetxt(self.path2load+'Ap.txt', self.A_w[:])
        np.savetxt(self.path2load+'b_w.txt', self.b_w[:])
        np.savetxt(self.path2load+'b_s.txt', self.b_s[:])
        # np.savetxt(self.path2load+'reg_b.txt', self.reg_b[:])
        # np.savetxt(self.path2load+'b_s.txt', self.b_s[:])
        if self.x0_prior==True:
            np.savetxt(self.path2load+'x0F1_sum.wtxt', self.x0F1_sum)
        np.savetxt(self.path2load+'reg_A.txt',self.reg_A)


        
    def run_single(self):
        self.prepare4iCSD()            
        self.iCSD()
        if self.type=='2d':
            self.plotCSD()
        else:
            print('3d case to plot using pyvista')
            self.plotCSD3d()
            self.plotCSD3d_pyvista()
        self.writeFIT()
        self.f.savefig(self.path2save+'iCSD', dpi = 600)
        plt.show()
        
    

        

             

