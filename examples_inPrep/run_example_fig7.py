"""
Inversion of current source density fig7 manuscript Peruzzo et al. 2020
-----------------------------------------------------------------------
Peruzzo, L., Chou, C., Wu, Y. et al. Imaging of plant current pathways for non-invasive root Phenotyping using a newly developed electrical current source density approach. Plant Soil 450, 567–584 (2020). https://doi.org/10.1007/s11104-020-04529-w
"""
import os
import matplotlib.pyplot as plt

path2files="./iCSD_fig7/"
from icsd3d_class import iCSD3d as i3d
from plotters import mpl_plot

# mpl_plot.showObs2d(path2files)

icsd=i3d(dirName=path2files)   
icsd.type='2d'
icsd.obs_err='sqrt' # choose between constant weight and w = 1/sqrt(abs(obs))
icsd.wr=1 #weight regularization
icsd.alphaSxy=True
icsd.x0_prior=False
icsd.x0_ini_guess=False # initial guess
icsd.method_m0='Pearson'
icsd.logTrans=True # initial guess

new = icsd.createSurvey(fname_obs='ObsData.txt',fname_sim='VRTeSim.txt')

# icsd.plotElecs=False

m0 = icsd.estimateM0(method_m0='F1',show=True)

# # icsd.clim=[0,0.1]
# # icsd.run_single()

# sol= icsd.invert(x0_prior=False,wr=0.1)

# fig, ax = plt.subplots()
# icsd.showResults(ax=ax)
# plt.show()

icsd.invert(pareto=True, pareto_MinErr=0.01, pareto_MaxErr=1)

# icsd=i3d(dirName=path2files)   
# #icsd.createSurvey()
# sol= icsd.invert(pareto=False,
#                  alphaSxy= False, 
#                  x0_ini_guess=False, 
#                  x0_prior=False,wr=33.3)
# icsd.showResults()

# icsd.ModelResolution(jacMi=305)

#%%
# vrtecoord = icsd.surveys[0].coord
# Robs = icsd.surveys[0].b
# #Rsim = np.reshape(icsd.surveys[0].A,[len(vrtecoord)*len(Robs)])
# Rsim = icsd.surveys[0].A
# Rsim.shape

# # buiding matrix A
# A = np.zeros((Rsim.shape[0] + 1  + len(vrtecoord), len(vrtecoord)))
# A[:Rsim.shape[0],:] = Rsim
# A[Rsim.shape[0] + 1,:] = 1 # charge conservation
# L1 = np.eye(len(vrtecoord)) + np.eye(len(vrtecoord), k=-1)*-1
# L1[0,1] = -1
# A[Rsim.shape[0]+1:,:] = L1

# # building matrix b
# b = np.zeros((len(vrtecoord), 1))
# b[:Rsim.shape[0],0] = Robs
# b[Rsim.shape[1],0] = 1 # charge conservation

# # solving
# x, residues, rank, singularValues = lstsq(Rsim, Robs)
# print(x)
# fig, ax = plt.subplots()
# ax.semilogy(x, 'o')

# # view results
# fig, ax = plt.subplots()
# cax = ax.scatter(vrtecoord[:,0], vrtecoord[:,1], 50, x)
# # cax = ax.tricontourf(vrtecoord[:,0], vrtecoord[:,1], x)
# fig.colorbar(cax, label='Source')
# #ax.scatter(vrtecoord[:,0],vrtecoord[:,1],color='b',marker='o')

