"""
Inversion of current source density apply to a roots imaging
------------------------------------------------------------
"""
import os
import matplotlib.pyplot as plt

# -----------------------------------#
# Exemple Landfill
# -----------------------------------#
path2files="./Landfill_3d/"

from icsd.icsd3d import iCSD3d as i3d 
from icsd.plotters import mpl_plot
import numpy as np


coords_elecs = np.loadtxt(path2files+'coords_elecs.txt')[:-3,:]
len(coords_elecs)
mpl_plot.showObs2d(path2files,coords_elecs,
                   filename='solution.data')

icsd3d_landfill=i3d(dirName=path2files)   
icsd3d_landfill.regMesh='strc'
icsd3d_landfill.type='3d'
icsd3d_landfill.obs_err='sqrt' # choose between constant weight and w = 1/sqrt(abs(obs))
icsd3d_landfill.wr=1 #weight regularization
icsd3d_landfill.alphaSxy=False
icsd3d_landfill.x0_prior=False
icsd3d_landfill.x0_ini_guess=False # initial guess

    
icsd3d_landfill.createSurvey(fname_obs='solution.data',fname_sim='Green.data')

icsd3d_landfill.invert(show=True)

# Unconstrainsted current source densities inversion
icsd3d_landfill.invert(wr=1,x0_prior=False)

icsd3d_landfill.invert()
# Estimate initial model
icsd3d_landfill.estimateM0(method_m0='F1', show=True)

# Constrainsted current source densities inversion
icsd3d_landfill.invert(regMesh='strc',wr=1,x0_prior=False)
# icsd3d_landfill.invertregMesh='strc',x0_prior=True)


import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D


fig = plt.figure()
ax=fig.gca(projection='3d')
icsd3d_landfill.showResults(ax=ax,
                            mesh='fwd_mesh_landfill.vtk')
plt.show()







