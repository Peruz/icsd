"""
Inversion of current source density apply to a roots imaging
------------------------------------------------------------
"""
import os
import matplotlib.pyplot as plt

# -----------------------------------#
# Exemple cotton
# -----------------------------------#
path2files="./Cotton_Soil/nail1/"

from icsd.icsd3d import iCSD3d as i3d 
from icsd.plotters import mpl_plot

icsd=i3d(dirName=path2files)   
icsd.type='2d'
icsd.obs_err='sqrt' # choose between constant weight and w = 1/sqrt(abs(obs))
icsd.wr=1 #weight regularization
icsd.alphaSxy=False
icsd.x0_prior=False
icsd.x0_ini_guess=False # initial guess
#icsd.icsd_init()

icsd.createSurvey(fname_obs='ObsData.txt',fname_sim='VRTeSim.txt')
# icsd.plotElecs=False

m0 = icsd.estimateM0(method_m0='F1',show=True)
m0 = icsd.estimateM0(method_m0='Pearson',show=True)

# icsd.clim=[0,0.1]
# icsd.run_single()

sol= icsd.invert(x0_prior=False,wr=1)


 
fig, ax = plt.subplots()
ax.scatter([0.2,0.4],[0.1,0.3],color='r',marker='v')
icsd.showResults(ax=ax)
# fig, ax = plt.subplots()
plt.show()


icsd.showResults()

# icsd.ModelResolution(jacMi=305)

# run pareto
icsd.pareto_MinErr=0.01
icsd.pareto_MaxErr=10
icsd.pareto_nSteps=10

icsd.knee=False
icsd.run_pareto()



# %% 