# Copyright (c) 2018 The ICSD Developers.
# https://github.com/Peruz/icsd/graphs/contributors
# Distributed under the terms of the BSD 3-Clause License.
# SPDX-License-Identifier: BSD-3-Clause
#
"""
Gridder
"""

import numpy as np


def mkGrid_XI_YI(coord_x, coord_y, nbe=500):
    """
    Grid for interpolation

    Parameters
    ----------
    * coord_x : 1D-arrays
        The x coordinates of the grid points
    * coord_y : float or 1D-array
        The y coordinate of the grid points
    Returns
    -------
    XI, YI : 1D-arrays
        Meshgrid points
    """
    Xm = np.linspace(min(coord_x), max(coord_x), nbe)
    Ym = np.linspace(min(coord_y), max(coord_y), nbe)
    (XI, YI) = np.meshgrid(Xm, Ym)

    return XI, YI
