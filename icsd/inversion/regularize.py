# Copyright (c) 2018 The ICSD Developers.
# https://github.com/Peruz/icsd/graphs/contributors
# Distributed under the terms of the BSD 3-Clause License.
# SPDX-License-Identifier: BSD-3-Clause
#
"""
Regularisation functions
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import diags


#%%
def nx_ny(coord):
    """find number of nodes in each direction, has to be a regular grid"""
    nx = np.unique(np.round(coord[:, 0], 3)).shape[0]
    ny = np.unique(np.round(coord[:, 1], 3)).shape[0]

    return nx, ny


def nx_ny_nz(coord):
    """find number of nodes in each direction, has to be a regular grid"""
    nx = np.unique(np.round(coord[:, 0], 3)).shape[0]
    ny = np.unique(np.round(coord[:, 1], 3)).shape[0]
    nz = np.unique(np.round(coord[:, 2], 3)).shape[0]

    return nx, ny, nz


#%% Smoothing methods for different mesh types

#%% 2d


def regularize_A(coord, nVRTe):
    """create and append rows for to A,
    for spatial regularization (simple model smoothing).
    Working only on 2d regular meshes
    """
    reg = []
    vrte = range(1, nVRTe + 1)

    nx, ny = nx_ny(coord)

    vrte = np.reshape(vrte, (ny, nx))
    for y in range(ny):
        for x in range(nx):
            minus = vrte[y, x]
            if x + 1 in range(nx):
                plus = vrte[y, x + 1]
                row = np.zeros(nVRTe, int)
                row[minus - 1] = -1
                row[plus - 1] = +1
                reg.append(row)
                if y + 1 in range(ny):
                    plus = vrte[y + 1, x]
                    row = np.zeros(nVRTe)
                    row[minus - 1] = -1
                    row[plus - 1] = +1
                    reg.append(row)
    reg_A = np.array(reg)

    return reg_A

#%% 3d
def regularize_A_3d(nVRTe, coord):
    """model smoothing consisting in creating and
    appending rows for spatial regularization to A
    """
    nx, ny, nz = nx_ny_nz(coord)
    reg = []
    vrte = range(1, nVRTe + 1)
    vrte = np.reshape(vrte, (ny, nx, nz))
    for z in range(nz):
        for y in range(ny):
            for x in range(nx):
                minus = vrte[y, x]
                if x + 1 in range(nx):
                    plus = vrte[y, x + 1]
                    row = np.zeros(nVRTe, int)
                    row[minus - 1] = -1
                    row[plus - 1] = +1
                    reg.append(row)
                    if y + 1 in range(ny):
                        plus = vrte[y + 1, x]
                        row = np.zeros(nVRTe)
                        row[minus - 1] = -1
                        row[plus - 1] = +1
                        reg.append(row)
    reg_A = np.array(reg)

    return reg_A


def regularize_A_UnstructuredMesh2d(coord, nVRTe, k_neighbors=2):
    """model smoothing consisting in creating and appending rows for spatial regularization to A.
    Adapted for unstructured mesh since it uses the k_neighbors method, default k=2. Also working on regular grid 2d"""
    reg = []
    for VRTEnb in range(nVRTe):
        dist = np.linalg.norm(coord[VRTEnb] - coord, axis=1)
        closest = np.argsort(dist)
        k = k_neighbors  # For each point, find the k closest current sources
        Ind = closest[1 : k + 1]
        row = np.zeros(
            nVRTe
        )  # add a line to the regularisation A with k non-null coefficients
        knorm = dist[closest[1 : k + 1]] / dist[closest[1 : k + 1]].sum(
            axis=0, keepdims=1
        )
        row[Ind] = -knorm
        row[VRTEnb] = 1  # = one for the actual current source
        reg.append(row)
        test = [1]
        mask = np.in1d(test, VRTEnb)
        reg_A = np.array(reg)

    return reg_A


def regularize_A_UnstructuredMesh3d(coord, nVRTe, k_neighbors=9):
    """model smoothing consisting in creating and appending rows for spatial regularization to A.
    Adapted for unstructured mesh since it uses the k_neighbors method, default k=4. Also working on regular grid 2d"""
    reg = []
    for VRTEnb in range(nVRTe):
        dist = np.linalg.norm(coord[VRTEnb] - coord, axis=1)
        closest = np.argsort(dist)
        k = k_neighbors  # For each point, find the k closest current sources
        Ind = closest[1 : k + 1]
        row = np.zeros(
            nVRTe
        )  # add a line to the regularisation A with k non-null coefficients
        knorm = dist[closest[1 : k + 1]] / dist[closest[1 : k + 1]].sum(
            axis=0, keepdims=1
        )
        row[Ind] = -knorm
        row[VRTEnb] = 1  # = one for the actual current source
        reg.append(row)
        test = [1]
        mask = np.in1d(test, VRTEnb)
        reg_A = np.array(reg)

    return reg_A


#%% Initiate vectors to build regularisation matrice for A, b

def regularize_b(reg_A):
    """initiate vector b with zeros, the length is determined by the number of regul rows in A"""
    reg_b = np.zeros(reg_A.shape[0])

    return reg_b


def regularize_w(reg_A, wr, x0_prior_flag, **kwargs):
    """create vector with weights, the length is determined by the number of regul rows in A such as
    .. math :: A = (G'*Wd*G + lambda*Wm)
               b = G'*Wd*d + lambda*Wm*m0;
    """
    if x0_prior_flag:
        # reg Wm (smallness + spatial reg) * lambda=" + str(wr))
        reg_w_0_b = np.ones(reg_A.shape[0]) * kwargs.get("x0") * wr
        reg_w_0_A = np.ones(reg_A.shape[0]) * wr
        return reg_w_0_b, reg_w_0_A
    else:
        reg_w = np.ones(reg_A.shape[0]) * wr

        return reg_w
