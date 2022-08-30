from icsd.icsd3d import iCSD3d

from icsd.importers import read
from icsd.inversion import priorM0, regularize, solve
from icsd.plotters import mpl_plot

import os 
from pathlib import Path
import numpy as np
import numpy.testing as npt
import pytest

import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

MODULE_DIR = Path(os.path.dirname(__file__))


def test_plot2d(file="/data/nail1/"):
	''' Plot estimate m0 and results of ISCD'''

	icsd=iCSD3d(dirName=str(MODULE_DIR) + file)   
	icsd.type='2d'
	icsd.obs_err='sqrt' # choose between constant weight and w = 1/sqrt(abs(obs))
	icsd.wr=1 #weight regularization
	icsd.alphaSxy=False
	icsd.x0_prior=False
	icsd.x0_ini_guess=False # initial guess 

	icsd.createSurvey(fname_obs='ObsData.txt',fname_sim='VRTeSim.txt')
	m0 = icsd.estimateM0(method_m0='F1',show=True)
	m0 = icsd.estimateM0(method_m0='Pearson',show=True)
	sol= icsd.invert(wr=1)
 
	fig, ax = plt.subplots()
	ax.scatter([0.2,0.4],[0.1,0.3],color='r',marker='v')
	icsd.showResults(ax=ax)
	# fig, ax = plt.subplots()
	plt.show()
	plt.close('all')


def test_pareto_plot(file="/data/nail1/"):
	''' Save results to pdf, plot L-curve'''

	icsd=iCSD3d(dirName=str(MODULE_DIR) + file)  
	icsd.type='2d'
	icsd.obs_err='sqrt' # choose between constant weight and w = 1/sqrt(abs(obs))
	icsd.wr=1 #weight regularization
	icsd.alphaSxy=False
	icsd.x0_prior=False
	icsd.x0_ini_guess=False # initial guess 


	icsd.createSurvey(fname_obs='ObsData.txt',fname_sim='VRTeSim.txt')
	m0 = icsd.estimateM0(method_m0='F1',show=True)
	m0 = icsd.estimateM0(method_m0='Pearson',show=True)
	plt.close('all')
	icsd.run_pareto()


def test_plot3d(file="/data/Landfill_3d/"):

	coords_elecs = np.loadtxt(str(MODULE_DIR) + file +'coords_elecs.txt')[:-3,:]
	mpl_plot.showObs2d(str(MODULE_DIR) + file,coords_elecs,
	                   filename='solution.data')

	icsd3d_landfill=iCSD3d(dirName=str(MODULE_DIR) + file)   
	icsd3d_landfill.regMesh='strc'
	icsd3d_landfill.type='3d'
	icsd3d_landfill.obs_err='sqrt' # choose between constant weight and w = 1/sqrt(abs(obs))
	icsd3d_landfill.wr=1 #weight regularization
	icsd3d_landfill.alphaSxy=False
	icsd3d_landfill.x0_prior=False
	icsd3d_landfill.x0_ini_guess=False # initial guess

	    
	icsd3d_landfill.createSurvey(fname_obs='solution.data',fname_sim='Green.data')

	# Unconstrainsted current source densities inversion
	icsd3d_landfill.invert(wr=1,x0_prior=False)

	# Estimate initial model
	icsd3d_landfill.estimateM0(method_m0='F1', show=True)

	# Constrainsted current source densities inversion
	icsd3d_landfill.invert(regMesh='strc',wr=1,x0_prior=False)
	# icsd3d_landfill.invertregMesh='strc',x0_prior=True)




	fig = plt.figure(projection='3d')
	ax=fig.gca()
	icsd3d_landfill.showResults(ax=ax,
	                            mesh='fwd_mesh_landfill.vtk')
	plt.show()


	fig = plt.figure(projection='3d')
	ax=fig.gca()
	icsd3d_landfill.showResults(ax=ax,
	                            mesh='fwd_mesh_landfill.vtk',
                                gif3d=True)
	plt.show()



if __name__ == '__main__':

    test_plot3d()




    



