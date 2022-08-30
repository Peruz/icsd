# Copyright (c) 2018 The ICSD Developers.
# https://github.com/Peruz/icsd/graphs/contributors
# Distributed under the terms of the BSD 3-Clause License.
# SPDX-License-Identifier: BSD-3-Clause
#
"""
Importers
"""

import os
import numpy as np
import pandas as pd

def load_coord(path, filename='VRTeCoord.txt', dim=2):
    """load coordinates of the virtual current sources

    Parameters
    ----------

    """
    coord = np.loadtxt(path + filename)
    if dim == 2:
        coord_x, coord_y = coord[:, 0], coord[:, 1]
        return coord_x, coord_y, coord
    else:
        coord_x, coord_y, coord_z = coord[:, 0], coord[:, 1], coord[:, 2]
        return coord_x, coord_y, coord_z, coord


def _read_pygimli(filename):
    ''' Read pygimli data format '''
    
    with open(filename) as f:
        lines = f.readlines()
        
    nb_of_sensors = int(lines[0])
    nb_of_quad = int(lines[nb_of_sensors+2])
    header = lines[nb_of_sensors+3]
    header = header.split(" ")
    header = header[1:-1]
    
    data = pd.read_csv(filename,sep='\t',
                      skiprows=nb_of_sensors+4,
                      header=None,
                      nrows=nb_of_quad)
    data.columns = header
    data["r"] =  data["rhoa"]/data["k"]
    
    return data["r"].to_numpy()
        

def load_obs(path, filename='ObsData.txt', index=None):
    """load the observations file (normalised voltages)

    Parameters
    ----------

    """
    if filename.endswith(".data"):  
        b = _read_pygimli(path + filename)
    else:
        b = np.loadtxt(path + filename)

    return b


def load_sim(path, data_sim='VRTeSim.txt'):
    """load the simulated green functions file
    """

    if isinstance(data_sim, str):
        if data_sim.endswith(".data"):  
            A = _read_pygimli(path + data_sim)
        else:
            A = np.loadtxt(path + data_sim)
    else:
        A = data_sim

    return A
