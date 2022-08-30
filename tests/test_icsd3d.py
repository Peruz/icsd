from icsd.icsd3d import iCSD3d

from icsd.importers import read
from icsd.inversion import priorM0, regularize, solve
from icsd.exporters import save

import os 
from pathlib import Path
import numpy as np
import numpy.testing as npt
import pytest



# ---
# TESTING IMPORTER MODULE

MODULE_DIR = Path(os.path.dirname(__file__))


def import_example_2d_regular():

	EXAMPLE_FILE_FOLDER = str(MODULE_DIR) + "/data/nail1/"

	ex = {}
	_, _, ex['coord'] = read.load_coord(EXAMPLE_FILE_FOLDER)
	ex['obs'] = read.load_obs(EXAMPLE_FILE_FOLDER)
	ex['sim'] = read.load_sim(EXAMPLE_FILE_FOLDER)

	ex['nVRTe'] = solve.check_nVRTe(ex['sim'], ex['obs'], ex['coord'])
	ex['sim_reshaped'] = solve.reshape_A(ex['sim'] , ex['nVRTe'])
	ex['nx'], ex['ny'] = regularize.nx_ny(ex['coord'])

	return ex


def import_example_3d_regular():

	EXAMPLE_FILE_FOLDER = str(MODULE_DIR) + "/data/Landfill_3d/"

	ex = {}
	_, _, ex['coord'] = read.load_coord(EXAMPLE_FILE_FOLDER)
	ex['obs'] = read.load_obs(EXAMPLE_FILE_FOLDER, filename='solution.data')
	ex['sim'] = read.load_sim(EXAMPLE_FILE_FOLDER, data_sim='Green.data')

	ex['nVRTe'] = solve.check_nVRTe(ex['sim'], ex['obs'], ex['coord'])
	ex['sim_reshaped'] = solve.reshape_A(ex['sim'] , ex['nVRTe'])
	ex['nx'], ex['ny'] = regularize.nx_ny(ex['coord'])

	return ex


# test support in numpy
# https://numpy.org/doc/stable/reference/generated/numpy.testing.assert_allclose.html


@pytest.fixture(name="coord")
def fixture_coord():
	"""
	Return a sample coord
	"""
	EXAMPLE_FILE_FOLDER = str(MODULE_DIR) + "/data/nail1/"
	coord_x, coord_y, coord = read.load_coord(EXAMPLE_FILE_FOLDER)
	return coord

@pytest.fixture(name="obs")
def fixture_obs():
	"""
	Return a sample obs
	"""
	EXAMPLE_FILE_FOLDER = str(MODULE_DIR) + "/data/nail1/"
	obs = read.load_obs(EXAMPLE_FILE_FOLDER)
	return obs

@pytest.fixture(name="sim")
def fixture_sim():
	"""
	Return a sample obs
	"""
	EXAMPLE_FILE_FOLDER = str(MODULE_DIR) + "/data/nail1/"
	sim = read.load_sim(EXAMPLE_FILE_FOLDER)
	return sim


def test_load_coord():
	# test geometry of the electrodes
	EXAMPLE_FILE_FOLDER = str(MODULE_DIR) + "/data/nail1/"
	coord_x, coord_y, coord = read.load_coord(EXAMPLE_FILE_FOLDER)
	assert ~np.isnan(coord).any()


def test_load_obs():
	# test for NaN values
	EXAMPLE_FILE_FOLDER = str(MODULE_DIR) + "/data/nail1/"
	obs = read.load_obs(EXAMPLE_FILE_FOLDER)
	assert ~np.isnan(obs).any()

def test_load_sim():
	# test load obs 
	# test for NaN values
	EXAMPLE_FILE_FOLDER = str(MODULE_DIR) + "/data/nail1/"
	sim = read.load_sim(EXAMPLE_FILE_FOLDER)
	assert ~np.isnan(sim).any()


# ---
# TESTING INVERSION MODULE

# PRIOR M0 TESTS
# --------------

def test_sum_product_moment():
	# test if sum equal to 1 and bounds between 0 and 1 
	ex = import_example_2d_regular()
	M0 = priorM0.product_moment(ex['sim_reshaped'],ex['obs'])
	sumM0 = sum(M0)
	assert min(M0)>=0


def test_sum_misfit_norm_F1():
	# test if sum equal to 1 and bounds between 0 and 1 
	ex = import_example_2d_regular()
	_, M0 = priorM0.misfitF1_2_initialX0(ex['sim_reshaped'],ex['obs'])
	sumM0 = sum(M0)
	npt.assert_almost_equal(sumM0,1,decimal=3)
	assert min(M0)>=0
	assert max(M0)<=1


# REGULARISATION
# --------------

def test_regular_grid():
	ex = import_example_2d_regular()
	# test if the grid is regular based on nb of points in each directions
	nx, ny = regularize.nx_ny(ex['coord'])
	npt.assert_allclose(ex['coord'][1:nx-1,0], ex['coord'][0,0], rtol=1e-5, atol=0)

# def test_smoothing_2d():
# 	ex = import_example_2d_regular()
# 	# test if regularisation of green functions matrice A has the correct size (==VRTE)
# 	reg_A = regularize.regularize_A(ex['coord'], ex['nVRTe'])
# 	assert reg_A.shape == (ex['nx'], ex['ny'])


# def test_smoothing_3d():
# 	# test if regularisation of green functions matrice A has the correct size (==VRTE)
# 	ex = import_example_3d_regular()
# 	reg_A = regularize.regularize_A_3d(ex['coord'], ex['nVRTe'])
# 	assert reg_A.shape == (ex['nx'], ex['ny'], ex['nz'])


# def test_smoothing_unstructured_mesh2d():
# 	# test if regularisation of green functions matrice A has the correct size (==VRTE)
# 	ex = import_example_2d_regular()
# 	reg_A = regularize.regularize_A_UnstructuredMesh2d(ex['coord'], ex['nVRTe'])
# 	nx, ny = regularize.nx_ny(ex['coord'])
# 	assert reg_A.shape == (ex['nx'], ex['ny'])

# def test_smoothing_unstructured_mesh3d():
# 	# test if regularisation of green functions matrice A has the correct size (==VRTE)
# 	ex = import_example_3d_regular()
# 	reg_A = regularize.regularize_A_UnstructuredMesh2d(ex['coord'], ex['nVRTe'])
# 	nx, ny = regularize.nx_ny(ex['coord'])
# 	assert reg_A.shape == (ex['nx'], ex['ny'], ex['nz'])

# SOLVER
# ------

# def test_current_conservation():
# 	# the sum of the current in the solution must be equal to 1
# 	x = solve.iCSD(A_w, b_w, dim, coord, path)
# 	assert np.sum(x.x) == 1


def test_export_sol2d():

	ex = import_example_2d_regular()

	icsd=iCSD3d(dirName=str(MODULE_DIR) + "/data/nail1/")   
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

	save.export_sol(ex['coord'], sol, dim=icsd.type, path=None, filename_root="ICSDsol.dat")


def test_export_sol3d():

	ex = import_example_3d_regular()

	icsd=iCSD3d(dirName=str(MODULE_DIR) + "/data/Landfill_3d/")   
	icsd.type='3d'
	icsd.obs_err='sqrt' # choose between constant weight and w = 1/sqrt(abs(obs))
	icsd.wr=1 #weight regularization
	icsd.alphaSxy=False
	icsd.x0_prior=False
	icsd.x0_ini_guess=False # initial guess 
	icsd.createSurvey(fname_obs='solution.data',fname_sim='Green.data')
	sol= icsd.invert(wr=1)

	save.export_sol(ex['coord'], sol, dim=icsd.type, path=None, filename_root="ICSDsol.dat")

    

# ---
# TESTING MAIN ICSD3 CLASS


if __name__ == '__main__':

    #test_sum_misfit_norm_F1()
    test_regular_grid()
