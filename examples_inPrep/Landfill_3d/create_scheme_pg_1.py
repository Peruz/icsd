# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 13:31:39 2020

@author: Benjamin
"""

import os
import numpy as np
import pygimli as pg
import pygimli.meshtools as mt
from pygimli.physics import ert
from pygimli.physics.ert import simulate as simulateERT
pyvista = pg.optImport("pyvista")

import matplotlib.pyplot as plt


def _check_area(area):
    """
    Check that the area argument is valid.
    For example, the west limit should not be greater than the east limit.
    """
    x1, x2, y1, y2 = area
    assert x1 <= x2, \
        "Invalid area dimensions {}, {}. x1 must be < x2.".format(x1, x2)
    assert y1 <= y2, \
        "Invalid area dimensions {}, {}. y1 must be < y2.".format(y1, y2)

def regular(area, shape, z=None):
    nx, ny = shape
    x1, x2, y1, y2 = area
    _check_area(area)
    xs = np.linspace(x1, x2, nx)
    ys = np.linspace(y1, y2, ny)
    # Must pass ys, xs in this order because meshgrid uses the first argument
    # for the columns
    arrays = np.meshgrid(ys, xs)[::-1]
    if z is not None:
        arrays.append(z*np.ones(nx*ny, dtype=np.float))
    return [i.ravel() for i in arrays]


#%% create model geometry

depth = 60  # mesh depth
width = 350 # mesh width

widthlandfill = 50 # 5, 15, 25, 35
Hlandfill = 5 # 5, 15, 25, 35

# -------------------
tickness_liner = 10

HAno = 2.5 
depthAno = -3.0 # 3.75, 13.75 or 23.75
widthAno = 5 # 5, 15, 25, 35
ticknessAno = tickness_liner/2
shift = (widthlandfill + tickness_liner)/2 - ticknessAno/2


# Define electrode position

mesh = pg.load('./fwd_mesh_landfill.vtk')





