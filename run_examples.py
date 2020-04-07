import os
import numpy as np

maindir='E:/Padova/Software/SourceInversion/icsd_dev/'
os.chdir(maindir)

# -----------------------------------#
# Exemple cotton
# -----------------------------------#
# path2files="examples/Cotton_Soil/Nail1/"
path2files="examples/Cotton_Soil/Nail1/"

from icsd2d_class import iCSD2d_Class as i2d
icsd=i2d(dirName=path2files)   
icsd.icsd_init()

# 
        
icsd.run_single()
icsd.run_pareto()

# -----------------------------------#
# apply here a smallness contrainst to the regularisation
from icsd3d_class import iCSD3d_Class as i3d

icsd=i3d(dirName=path2files)   
icsd.type='2d'
icsd.obs_err='sqrt' # choose between constant weight and w = 1/sqrt(abs(obs))
icsd.wr=10 #weight regularization
icsd.x0_prior=True
icsd.x0_ini_guess=False # initial guess
icsd.icsd_init()

icsd.run_single()

# run pareto
icsd.pareto_MinErr=0.001
icsd.pareto_MaxErr=1
icsd.knee=True
icsd.run_pareto()

# -----------------------------------#
# Exemple with a 3d landfill geometry --> ARTICLE comp. Geosciences ?
# -----------------------------------#
from icsd3d_class import iCSD3d_Class as i3d
path2files="examples/Landfill_3d/"

icsd3d_landfill=i3d(dirName=path2files)   
icsd3d_landfill.type='3d'
icsd3d_landfill.sim="SNoAno.txt"
icsd3d_landfill.obs="ONoAno_synt.txt"
icsd3d_landfill.coord_file="VRTeCoord.txt"
icsd3d_landfill.regMesh=='strc'
icsd3d_landfill.x0_prior=True
icsd3d_landfill.x0_ini_guess=True # initial guess
# icsd3d_landfill.wr=1
icsd3d_landfill.obs_err='const' # choose between constant weight and w = 1/sqrt(abs(obs))
icsd3d_landfill.mesh='mesh3d_rhomap.vtk'
icsd3d_landfill.plotElecs=False     
icsd3d_landfill.icsd_init()  
icsd3d_landfill.run_single()


icsd3d_landfill.pareto_MinErr=0.01
icsd3d_landfill.pareto_MaxErr=100
icsd3d_landfill.knee=True
icsd3d_landfill.run_pareto()


# -----------------------------------#
# Example time-lapse salt tracer
# show the importance of an anisotropic spatial regularisation
# -----------------------------------#
from icsd3d_class import iCSD3d_Class as i3d
import pyvista as pv

tstep=2
path2files="examples/Salt_tracer/t"+str(tstep) +'/'

icsd3d_TL_RWU=i3d(dirName=path2files)   
icsd3d_TL_RWU.type='2d'
icsd3d_TL_RWU.sim="SMALMtimeReg"+str(tstep)+".txt"
icsd3d_TL_RWU.obs="OMALMtimeReg"+str(tstep)+"_synt.txt"
icsd3d_TL_RWU.coord_file="VRTeCoord.txt"
icsd3d_TL_RWU.regMesh='strc'
icsd3d_TL_RWU.x0_prior=False
icsd3d_TL_RWU.x0_ini_guess=True # initial guess
# icsd3d_TL_RWU.wr=60
icsd3d_TL_RWU.plotElecs=False

icsd3d_TL_RWU.alphaSxy=True
icsd3d_TL_RWU.alphaSx=1 # weight on model smoothness in x-direction [TO IMPLEMENT]
icsd3d_TL_RWU.alphaSy=1 # weight on model smoothness in x-direction [TO IMPLEMENT]
# icsd3d_TL_RWU.mesh='Figures/ERTmodel_reg_'+str(tstep)+'.vtk'      

icsd3d_TL_RWU.icsd_init() 
icsd3d_TL_RWU.run_single()

icsd3d_TL_RWU.pareto_MinErr=0.1
icsd3d_TL_RWU.pareto_MaxErr=200
icsd3d_TL_RWU.knee=True
icsd3d_TL_RWU.run_pareto()

# -----------------------------------#
# Example time-lapse data SOIL article
#-> Importance of the smallness minimisation
# -----------------------------------#
from icsd3d_class import iCSD3d_Class as i3d
import pyvista as pv
import matplotlib.pyplot as plt

path2files="examples/RWU_TimeLapse/ADAM/"

icsd3d_TL_RWU=i3d(dirName=path2files)   
icsd3d_TL_RWU.type='3d'
icsd3d_TL_RWU.sim="VRTeSimADAMM0.txt"
icsd3d_TL_RWU.obs="ObsData_ADAMM0.txt"
icsd3d_TL_RWU.coord_file="VRTeCoord.txt"
icsd3d_TL_RWU.regMesh='unstrc'
icsd3d_TL_RWU.x0_prior=True
icsd3d_TL_RWU.x0_ini_guess=True # initial guess
icsd3d_TL_RWU.wr=60
icsd3d_TL_RWU.plotElecs=True
icsd3d_TL_RWU.icsd_init() 
# icsd3d_TL_RWU.mesh='mesh3d_rhomap.vtk'      
icsd3d_TL_RWU.run_single()


icsd3d_TL_RWU.pareto_MinErr=0.1
icsd3d_TL_RWU.pareto_MaxErr=200
icsd3d_TL_RWU.knee=True
icsd3d_TL_RWU.run_pareto()
