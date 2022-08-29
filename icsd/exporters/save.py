# -*- coding: utf-8 -*-
"""
Created on Mon May 11 15:19:59 2020

@author: Benjamin
"""

import os

import numpy as np

#%% SAVE FCTS


def export_sol(coord, vec_sol, dim, path=None, filename_root="Solution.dat"):
    """Export solution x

    Parameters
    ----------
    outputdir : str
        Path where the .csv files will be saved.
    """

    if path is None:
        cwd = os.getcwd()
        path = cwd

    filename = path + filename_root
    if dim == "2d":
        ExpSOL = np.vstack([coord[:, 0], coord[:, 1], vec_sol])
        ExpSOL = ExpSOL.T
        f = open(filename, "w")
        np.savetxt(
            f, ExpSOL, fmt="%1.2f %1.2f %1.6f", delimiter="\t", header="X Y i"
        )  # X is an array
        f.close()
    else:
        ExpSOL = np.vstack([coord[:, 0], coord[:, 1], coord[:, 2], vec_sol])
        ExpSOL = ExpSOL.T
        f = open(filename, "w")
        np.savetxt(
            f, ExpSOL, fmt="%1.2f %1.2f %1.2f %1.6f", delimiter="\t", header="X Y Z i"
        )  # X is an array
        f.close()
