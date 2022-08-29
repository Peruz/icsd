# -*- coding: utf-8 -*-
"""
Created on Mon May 11 15:18:31 2020

@author: Benjamin
"""

import os
import numpy as np

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


def load_obs(path, filename='ObsData.txt', index=None):
    """load the observations file (normalised voltages)

    Parameters
    ----------

    """
    if filename.endswith(".data"):  # TDIP data importer
        bfile = pg.load(filename)
        b = (bfile["m" + str(index)] / bfile["k"]).array()

    else:
        b = np.loadtxt(path + filename)

    return b


def load_sim(path, data_sim='VRTeSim.txt'):
    """load the simulated green functions file
    """
    if isinstance(data_sim, str):
        A = np.loadtxt(path + data_sim)
    else:
        A = data_sim

    return A


def load_geom(path):
    """load the geometry of the acquisition (*geom file custum for Mise-à-la-masse data)

    Parameters
    ----------

    """
    geom_files = [f for f in os.listdir(path) if f.endswith(".geom")]
    if len(geom_files) != 1:
        raise ValueError("should be only one geom file in the current directory")

    fileNameElec = geom_files[0]
    line_number = 0
    line_of_injection = []
    line_of_remotes = []
    # Open the file in read only mode
    with open(path + fileNameElec, "r") as read_obj:
        # Read all lines in the file one by one
        for line in read_obj:
            # For each line, check if line contains the string
            # print(line)
            line_number += 1
            if ("#Remote") in line:
                # If yes, then add the line number & line as a tuple in the list
                line_of_remotes.append((line_number))
            if ("#Injection") in line:
                line_of_injection.append((line_number))

    # print(line_of_injection)
    # print(line_of_remotes)
    RemLineNb = int(line_of_remotes[0]) - 1
    Injection = int(line_of_injection[0]) - 1

    coordE = np.loadtxt(path + fileNameElec)
    pointsE = np.vstack(coordE[:RemLineNb, 1:4])

    return RemLineNb, Injection, coordE, pointsE

