# Copyright (c) 2018 The ICSD Developers.
# https://github.com/Peruz/icsd/graphs/contributors
# Distributed under the terms of the BSD 3-Clause License.
# SPDX-License-Identifier: BSD-3-Clause
#
"""
Solver
"""
import numpy as np
from scipy.optimize import least_squares, lsq_linear

#%% Solve linear system Ax=b


def iCSD(A_w, b_w, dim, coord, path, **kwargs):
    """
    Solve linear system, given weigted A matrix (VRTe, constrain, regul) and weigted b (observations).

    Parameters
    ----------
    * x0_ini_guess : 1D-arrays
        Initial guess
    * A_w : 1D-arrays
        Kernel of green functions
    * b_w : 1D-array
        Weigted observations
    * dim : int
        Survey dimension i.e 2d or 3d
    * coord : 1D-arrays
        Coordinates of the virtual sources
    Returns
    -------
    x : 1D-arrays
        Solution
    """

    if kwargs.get("x0") is None:
        # No initial guess use lsq_linear solver
        x = lsq_linear(A_w, b_w, bounds=(0, 1), verbose=0)
        print("*" * 20)
        print("CURRENT Sum=" + str(np.sum(x.x)))
    else:
        # Initial guess x0 use least_squares solver
        a = A_w
        b = b_w
        def func(x, a, b):
            return b - np.dot(a, x)

        x = least_squares(
            func, x0=kwargs.get("x0"), bounds=(0, 1), args=(a, b)
        )  # Add initial guess
    return x


#%% MAKE LINEAR SYSTEM


def check_nVRTe(A, b, coord):
    if A.shape[0] / b.shape[0] == coord.shape[0]:
        nVRTe = coord.shape[0]
    else:
        raise ValueError("### dimensions of the files do not agree")
    return nVRTe


def reshape_A(A, nVRTe):
    A = A.reshape((-1, nVRTe), order="F")
    return A


def obs_w_f(obs_err, b, errRmin, sd_rec=None):
    """weight the observations, can also ignore observations by setting w = 0"""
    if obs_err == "const":
        obs_w = np.ones(b.shape[0])
    elif obs_err == "sqrt":
        obs_w = 1 / np.sqrt(np.abs(b))
        if (b == 0).any():
            print("b = 0 could be a problem, check presence of 0 and filter if needed")
        if (obs_w >= 10 * errRmin).any():
            print("errRmin not correctly set, adjust")
        obs_w[obs_w >= errRmin] = 1
    return obs_w


#%% CONSTRAIN


def con_A_f(A):
    """Set current conservation constrainst on A (rows of ones)"""
    con_A = np.ones(A.shape[1])
    return con_A


def con_b_f(b):
    """Set current conservation constrainst on b"""
    con_b = np.ones(1)
    return con_b


def con_w_f(wc):
    """Set current conservation constrainst weight; default is wc=1e6"""
    con_w = np.ones(1) * wc
    return con_w


#%% VERTICAL STACK EQUATIONS


def stack_A(A, con_A, reg_A):
    """Stack A (green fcts), constrainsts and regularisation"""
    # con_A = _con_A_f(A)
    A_s = np.vstack((A, con_A, reg_A))
    return A_s


def stack_b(b, con_b, reg_b):
    """Stack b, constrainsts and regularisation"""
    # con_b = _con_b_f(b)
    b_s = np.concatenate((b, con_b, reg_b))
    return b_s


def stack_w(obs_w, con_w, x0_prior, **kwargs):
    """create vector with weights for observation, constrain, and regularization
    then use it as diagonal for the weight matrix"""
    # con_w = _con_w_f(wc)
    reg_w = kwargs.get("reg_w")
    reg_w_0_A = kwargs.get("reg_w_0_A")
    reg_w_0_b = kwargs.get("reg_w_0_b")

    if x0_prior:  # if relative smallness
        wa = np.concatenate((obs_w, con_w, reg_w_0_A))
        wb = np.concatenate((obs_w, con_w, reg_w_0_b))
        W = np.zeros((wa.shape[0], wa.shape[0]))
        np.fill_diagonal(W, wa)
        W_s_A = W
        np.fill_diagonal(W, wb)
        W_s_b = W
        return W_s_A, W_s_b

    else:
        w = np.concatenate((obs_w, con_w, reg_w))
        W = np.zeros((w.shape[0], w.shape[0]))
        np.fill_diagonal(W, w)
        W_s = W
        return W_s


#%% APPLY WEIGHTS


def weight_A(x0_prior, A_s, **kwargs):
    """Apply the weights to A"""

    if x0_prior:
        A_w = np.matmul(kwargs.get("W_s_A"), A_s)
    else:
        A_w = np.matmul(kwargs.get("W_s"), A_s)
    return A_w


def weight_b(x0_prior, b_s, **kwargs):
    """Apply the weights to b"""
    if x0_prior:
        b_w = np.matmul(b_s, kwargs.get("W_s_b"))
    else:
        b_w = np.matmul(b_s, kwargs.get("W_s"))
    return b_w
