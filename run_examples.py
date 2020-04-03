import os
import numpy as np

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
icsd.run_pareto()

# -----------------------------------#
# apply here a smallness contrainst to the regularisation
from icsd3d_class import iCSD3d_Class as i3d

icsd=i3d(dirName=path2files)   
icsd.type='2d'
icsd.obs_err='sqrt' # choose between constant weight and w = 1/sqrt(abs(obs))
icsd.wr=100 #weight regularization
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
path2files="examples/Landfill_3d/"

icsd3d_landfill=i3d(dirName=path2files)   
icsd3d_landfill.type='3d'
icsd3d_landfill.sim="SNoAno.txt"
icsd3d_landfill.obs="ONoAno_synt.txt"
icsd3d_landfill.coord_file="VRTeCoord.txt"
icsd3d_landfill.regMesh=='strc'
icsd3d_landfill.x0_prior=False
icsd3d_landfill.x0_ini_guess=False # initial guess
icsd3d_landfill.wr=200
icsd3d_landfill.obs_err='sqrt' # choose between constant weight and w = 1/sqrt(abs(obs))

icsd3d_landfill.icsd_init()        
icsd3d_landfill.run_single()


icsd3d_landfill.pareto_MinErr=0.001
icsd3d_landfill.pareto_MaxErr=100
icsd3d_landfill.knee=True
icsd3d_landfill.run_pareto()

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
icsd3d_TL_RWU.x0_prior=False
icsd3d_TL_RWU.x0_ini_guess=False # initial guess
# icsd3d_TL_RWU.wr=200
icsd3d_TL_RWU.plotElecs=True
icsd3d_TL_RWU.icsd_init()        
# icsd3d_TL_RWU.run_single()

icsd3d_TL_RWU.pareto_MinErr=0.001
icsd3d_TL_RWU.pareto_MaxErr=200
icsd3d_TL_RWU.knee=True
icsd3d_TL_RWU.run_pareto()
