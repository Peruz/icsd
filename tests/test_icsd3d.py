from icsd.icsd3d import iCSD3d

from icsd.importers import read
from icsd.inversion import priorM0, regularize, solve

import os 
from pathlib import Path
import numpy as np
import numpy.testing as npt
import pytest



# ---
# TESTING IMPORTER MODULE

MODULE_DIR = Path(os.path.dirname(__file__))
EXAMPLE_FILE_FOLDER = str(MODULE_DIR) + "/data/nail1/"


def import_example_2d_regular():
	ex = {}
	_, _, ex['coord'] = read.load_coord(EXAMPLE_FILE_FOLDER)
	ex['obs'] = read.load_obs(EXAMPLE_FILE_FOLDER)
	ex['sim'] = read.load_sim(EXAMPLE_FILE_FOLDER)

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
	coord_x, coord_y, coord = read.load_coord(EXAMPLE_FILE_FOLDER)
	return coord

@pytest.fixture(name="obs")
def fixture_obs():
	"""
	Return a sample obs
	"""
	obs = read.load_obs(EXAMPLE_FILE_FOLDER)
	return obs

@pytest.fixture(name="sim")
def fixture_sim():
	"""
	Return a sample obs
	"""
	sim = read.load_sim(EXAMPLE_FILE_FOLDER)
	return sim



def test_load_coord():
	# test geometry of the electrodes
	coord_x, coord_y, coord = read.load_coord(EXAMPLE_FILE_FOLDER)
	assert ~np.isnan(coord).any()


def test_load_obs():
	# test for NaN values
	obs = read.load_obs(EXAMPLE_FILE_FOLDER)
	assert ~np.isnan(obs).any()

def test_load_sim():
	# test load obs 
	# test for NaN values
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
	assert npt.assert_allclose(sumM0,1,rtol=1e-3)
	assert min(M0)>=0
	assert max(M0)<=0


def test_sum_misfit_norm_F1():
	# test if sum equal to 1 and bounds between 0 and 1 
	ex = import_example_2d_regular()
	_, M0 = priorM0.misfitF1_2_initialX0(ex['sim_reshaped'],ex['obs'])
	sumM0 = sum(M0)
	assert npt.assert_allclose(sumM0,1,rtol=1e-3)
	assert min(M0)>=0
	assert max(M0)<=0


# REGULARISATION
# --------------

def test_regular_grid(coord):
	ex = import_example_2d_regular()
	# test if the grid is regular based on nb of points in each directions
	if coord.ndim==2:
		nx, ny = regularize.nx_ny(coord)
		assert nx == ny
		return nx, ny
	else:
		nx_ny_nz = regularize.nx_ny_nz(coord)
		assert nx == ny == nz

def test_smoothing_2d(coord, nVRTe, nx, ny):
	ex = import_example_2d_regular()
	# test if regularisation of green functions matrice A has the correct size (==VRTE)
	reg_A = regularize.regularize_A(ex['coord'], ex['nVRTe'])
	assert reg_A.shape == (nx, ny)


def test_smoothing_3d(coord, nVRTe, nx, ny, nz):
	# test if regularisation of green functions matrice A has the correct size (==VRTE)
	reg_A = regularize.regularize_A_3d(coord, nVRTe)
	assert reg_A.shape == (nx, ny, nz)


def test_smoothing_unstructured_mesh(coord, nVRTe, k_neighbors):
	# test if regularisation of green functions matrice A has the correct size (==VRTE)
	if coord.ndims==2:
		reg_A = regularize.regularize_A_UnstructuredMesh2d(coord, nVRTe)
		nx, ny = regularize.nx_ny(coord)
		assert reg_A.shape == (nx, ny, nz)
	else:
		reg_A = regularize.regularize_A_UnstructuredMesh3d(coord, nVRTe)
		nx_ny_nz = regularize.nx_ny_nz(coord)
		assert reg_A.shape == (nx, ny, nz)


# SOLVER
# ------

def test_current_conservation():
	# the sum of the current in the solution must be equal to 1
	x = solve.iCSD(A_w, b_w, dim, coord, path)
	assert np.sum(x.x) == 1



# ---
# TESTING MAIN ICSD3 CLASS

def test_init_exemple():
    # ex = import_example_2d_regular()

    # survey.obs_w = obs_w_f(self.obs_err, survey.b, self.errRmin, sd_rec=None)

    # # set constrain (curent conservation)
    # survey.con_A = con_A_f(survey.A)
    # survey.con_b = con_b_f(survey.b)
    # survey.con_w = con_w_f(self.wc)

    # # append spatial regularization (add lines to the matrice)
    # survey.reg_A = self._parseModelReg(survey)
    # survey.reg_b = regularize_b(self.reg_A)

    # # stack data, constrain, and regularization
    # survey.A_s = stack_A(survey.A, survey.con_A, survey.reg_A)
    # survey.b_s = stack_b(survey.b, survey.con_b, survey.reg_b)
    pass


def test_prepare4icsd_exemple():
	pass
	

