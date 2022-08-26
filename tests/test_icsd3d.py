from icsd.icsd3d import iCSD3d
from icsd.inversion import priorM0, regularize, solve

import numpy as np


# PRIOR M0 TESTS
# --------------
A = np.ones([3,3])
b = np.ones([1,3])

def test_sum_product_moment():
	# test if sum equal to 1 and bounds between 0 and 1 
	M0 = priorM0.product_moment(A,b)
	assert sum(M0) == 1
	assert min(M0)>=0
	assert max(M0)<=0


def test_sum_misfit_norm_F1():
	# test if sum equal to 1 and bounds between 0 and 1 
	_, M0 = priorM0.misfitF1_2_initialX0(A,b)
	assert sum(M0) == 1
	assert min(M0)>=0
	assert max(M0)<=0


# REGULARISATION
# --------------

def test_regular_grid(coord):
	# test if the grid is regular based on nb of points in each directions
	if coord.ndims==2:
		nx, ny = regularize.nx_ny(coord)
		assert nx == ny
	else:
		nx_ny_nz = regularize.nx_ny_nz(coord)
		assert nx == ny == nz

def test_smoothing_2d(coord, nVRTe, nx, ny):
	# test if regularisation of green functions matrice A has the correct size (==VRTE)
	reg_A = regularize.regularize_A(coord, nVRTe)
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