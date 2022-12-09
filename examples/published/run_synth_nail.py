"""
Multiple sources inversion of CSD
-----------------------------------
The procedure to generate the data is describe in Peruzzo et al. (2020). 
The synthetic test evaluates the iCSD capability to locate a current point source (white dot). 
The white circles represent the VRTe and the white square the return electrode. 
The source of current was centered between four VRTe to highlight the distribution of the current density.

*Peruzzo, L., Chou, C., Wu, Y. et al. Imaging of plant current pathways for non-invasive root Phenotyping using a newly developed electrical current source density approach. Plant Soil 450, 567–584 (2020). https://doi.org/10.1007/s11104-020-04529-w*
"""



# %%
# **Import the icsd package**
import matplotlib.pyplot as plt
from icsd.icsd3d import iCSD3d as i3d 

# %%
# **Change inversion parameters** 
path2files="./iCSD_fig7/"
icsd=i3d(dirName=path2files)   

from icsd.icsd3d import iCSD3d as i3d 

icsd=i3d(dirName=path2files)   
icsd.type='2d'
icsd.obs_err='sqrt' # choose between constant weight and w = 1/sqrt(abs(obs))
icsd.wr=1 #weight regularization
icsd.alphaSxy=False
icsd.x0_prior=False
icsd.x0_ini_guess=False # initial guess

icsd.createSurvey(fname_obs='ObsData.txt',fname_sim='VRTeSim.txt')
sol= icsd.invert(x0_prior=False,wr=1,show=True)



