"""
ICSD: calibration in water
==========================

This example illustrate the use of the ICSD code for a calibration experiment consisting in emerging an isolate cable in a box filled with water. 
The tip of the cable release the current into the water. The aim is to find the location of the cable tip after measuring the potentials at different location of the box.
For further information see Peruzzo et al. (2019)
"""


# Import packages
import matplotlib.pyplot as plt
from icsd.icsd3d import iCSD3d as i3d 
from icsd.plotters import mpl_plot


path2files="./Cotton_Soil/nail1/"

icsd=i3d(dirName=path2files)   
icsd.type='2d'
icsd.obs_err='sqrt' # choose between constant weight and w = 1/sqrt(abs(obs))
icsd.wr=1 #weight regularization
icsd.alphaSxy=False
icsd.x0_prior=False
icsd.x0_ini_guess=False # initial guess

# create a survey container
icsd.createSurvey(fname_obs='ObsData.txt',fname_sim='VRTeSim.txt')

# First method to estimate m0
m0 = icsd.estimateM0(method_m0='F1',show=True)

# Second method to estimate the initial model based on Pearson analysis
m0 = icsd.estimateM0(method_m0='Pearson',show=True)

# Invert the data
sol= icsd.invert(x0_prior=False,wr=1)


# Plot the results 
fig, ax = plt.subplots()
icsd.showResults(ax=ax)
plt.show()


# Another way to plot the result
icsd.showResults()


# Run pareto analysis between 0.01 to 10 of regularisation weight
icsd.pareto_MinErr=0.01
icsd.pareto_MaxErr=10
icsd.pareto_nSteps=10
icsd.run_pareto()

