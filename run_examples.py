import os

maindir='E:/Padova/Software/SourceInversion/icsd_dev/'
os.chdir(maindir)

# -----------------------------------#
# Exemple cotton
# -----------------------------------#
path2files="examples/Cotton_Soil/Nail1/"

from icsd2d_class import iCSD2d_Class as i2d
icsd=i2d(dirName=path2files)   
icsd.icsd_init()


icsd.run_single()
# icsd.run_pareto()

# -----------------------------------#
# apply here a smallness contrainst to the regularisation
from icsd3d_class import iCSD3d_Class as i3d

icsd=i3d(dirName=path2files)   
icsd.type='2d'
icsd.obs_err='sqrt' # choose between constant weight and w = 1/sqrt(abs(obs))
icsd.x0_prior=False
icsd.x0_ini_guess=False # initial guess
icsd.icsd_init()

icsd.run_single()

# run pareto
icsd.pareto_MinErr=0.001
icsd.pareto_MaxErr=100
icsd.knee=True
icsd.run_pareto()

# -----------------------------------#
# Exemple with a 3d landfill geometry --> ARTICLE comp. Geosciences ?
# -----------------------------------#
path2files="examples/Ano0_Test89_lambda_200/"

from icsd3d_class import iCSD3d_Class as i3d
icsd3d=i3d(dirName=path2files)   
icsd3d.sim="SNoAno.txt"
icsd3d.obs="ONoAno_synt.txt"
icsd3d.coord_file="VRTeCoord.txt"
icsd3d.RegMesh=='strc'
icsd3d.x0_prior=False
icsd3d.x0_ini_guess=False # initial guess
icsd3d.wr=200
icsd3d.obs_err='sqrt' # choose between constant weight and w = 1/sqrt(abs(obs))

icsd3d.icsd_init()        
icsd3d.run_single()


icsd3d.pareto_MinErr=0.001
icsd3d.pareto_MaxErr=100
icsd3d.knee=True
icsd3d.run_pareto()


    # filename = 'ExportSol.dat'
    # data_2_plot =np.genfromtxt(filename)
    # coord_x= data_2_plot[:,0]
    # coord_y= data_2_plot[:,1]
    # coord_z= data_2_plot[:,2]
    # step=(max(coord_x)-min(coord_x))/10
    # num =10
    
    # xlin=np.arange(min(coord_x),max(coord_x),step)
    # ylin=np.arange(min(coord_y),max(coord_y),step)
    # zlin=np.arange(min(coord_z),max(coord_z),step)
    # #generate new grid
    # X,Y,Z=np.meshgrid(xlin,ylin,zlin)
    
    # from scipy.interpolate import griddata as gd
    # #interpolate "data.v" on new grid "inter_mesh"
    # V = gd((coord_x,coord_y,coord_z), data_2_plot[:,3], (X,Y,Z), method='linear')
    
    # from mpl_toolkits.mplot3d import Axes3D
    
    # fig = plt.figure()
    # ax=fig.gca(projection='3d')
    # # Elecs, Elec99, Elec991, ax = MR.SensorsPos(mesh3d,pltRem=True,pltVRTE=False,savefig=False, axs=ax)
    # sc=ax.scatter(coord_x, coord_y, coord_z, c=data_2_plot[:,3], cmap ='coolwarm', s=data_2_plot[:,3]*1e4,
    #               vmin=0, vmax=0.01)
    # ax.plot(coords[:,0],coords[:,1],'k--')
    # cbar = plt.colorbar(sc)
    # cbar.set_label('# current density')
    # ax.view_init(azim=-101, elev=35)
    # # plt.legend(sc,'Test1')
    # title= 'Scattered current sources density -' + ssName
    # # plt.title(title)
    # # plt.savefig(path2file + ssName, dpi=550,bbox_inches='tight',pad_inches = 0)
    # plt.show()
    
# -----------------------------------#
# Example time-lapse data SOIL article
#-> Importance of the smallness minimisation
# -----------------------------------#


