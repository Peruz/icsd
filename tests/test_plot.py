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


MODULE_DIR = Path(os.path.dirname(__file__))
EXAMPLE_FILE_FOLDER = str(MODULE_DIR) + "/data/nail1/"


def test_plot2d():
	''' Plot estimate m0 and results of ISCD'''

	icsd=iCSD3d(dirName=EXAMPLE_FILE_FOLDER)   
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


def test_pareto_plot():
	''' Save results to pdf, plot L-curve'''

	icsd=iCSD3d(dirName=EXAMPLE_FILE_FOLDER)   
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

