# Copyright (c) 2018 The ICSD Developers.
# https://github.com/Peruz/icsd/graphs/contributors
# Distributed under the terms of the BSD 3-Clause License.
# SPDX-License-Identifier: BSD-3-Clause
#
"""
Exporter
"""
import os
import numpy as np

#%% SAVE FCTS


def export_sol(coord, vec_sol, dim, path=None, filename_root="ICSDsol.dat"):
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
        ExpSOL = np.vstack([coord[:, 0], coord[:, 1], vec_sol.x])
        f = open(filename, "w")
        np.savetxt(
            f, ExpSOL.T, fmt="%1.2f %1.2f %1.6f", delimiter="\t", header="X Y i"
        )
        f.close()
    else:
        ExpSOL = np.vstack([coord[:, 0], coord[:, 1], coord[:, 2], vec_sol.x])
        f = open(filename, "w")
        np.savetxt(
            f, ExpSOL.T, fmt="%1.2f %1.2f %1.2f %1.6f", delimiter="\t", header="X Y Z i"
        )  # X is an array
        f.close()
