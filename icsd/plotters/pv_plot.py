# Copyright (c) 2018 The ICSD Developers.
# https://github.com/Peruz/icsd/graphs/contributors
# Distributed under the terms of the BSD 3-Clause License.
# SPDX-License-Identifier: BSD-3-Clause
#
"""
3D plots using pyvista
"""

import os
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv


#%% Plot CSD in 3d using pyvista package

def plotCSD3d_pv(solution,
    coords, path=None, **kwargs
):  # add kwards or *args for pyvista plot

    if path is None:
        cwd = os.getcwd()
        path = cwd

    coord_x = coords[:, 0]
    coord_y = coords[:, 1]
    coord_z = coords[:, 2]

    opacity = [0, 0, 0.1, 0.3, 0.6, 0.9, 1]
    grid = pv.UniformGrid()
    spc = (max(coord_x) - min(coord_x)) / 10
    xdim = int(round((max(coord_x) - min(coord_x)) / spc))
    ydim = int(round((max(coord_y) - min(coord_y)) / spc))
    zdim = int(round((max(coord_z) - min(coord_z)) / spc))
    grid.dimensions = (xdim, ydim, zdim)
    grid.dimensions = np.array(grid.dimensions) + 1
    grid.origin = (
        min(coord_x),
        min(coord_y),
        min(coord_z),
    )  # The bottom left corner of the data set
    grid.spacing = (spc, spc, spc)  # These are the cell sizes along each axis

    pv.set_plot_theme("document")
    poly = pv.PolyData(coords)
    pvfig = pv.Plotter(notebook=False, window_size=[600, 600])

    if kwargs.get("mesh") is not None:
        ModelVtk = pv.read(path + kwargs.get("mesh"))
        cmap = plt.cm.get_cmap("viridis", 2)
        pvfig.add_bounding_box()
        pvfig.add_mesh(
            cmap=cmap, mesh=ModelVtk, scalars="Marker", opacity=0.2
        )  # add a dataset to the scene

    pvfig.add_mesh(
        poly,
        point_size=15.0,
        scalars=solution,
        opacity=opacity,
        render_points_as_spheres=True,
        cmap="jet",
        show_scalar_bar=False,
    )
    print("interpolation spacing=" + str(spc))
    interpolated = grid.interpolate(poly, radius=spc * 2)
    cmap = plt.cm.get_cmap("jet", 10)
    contours = interpolated.contour()
    pvfig.add_mesh(contours, show_scalar_bar=False, opacity=opacity, cmap="jet")
    pvfig.show_bounds(
        bounds=[
            min(coord_x),
            max(coord_x),
            min(coord_y),
            max(coord_y),
            min(coord_z),
            0,
        ],
        font_size=16,
    )
    pvfig.add_axes()
    pvfig.show_axes()
    pvfig.add_scalar_bar(
        "Normalized Current density", width=0.25, vertical=False, position_x=0.3
    )
    pvfig.show(auto_close=True)

    knee = kwargs.get("knee")
    wr = kwargs.get("wr")
    if knee is True:
        KneeWr = kwargs.get("KneeWr")
    gif3d = kwargs.get("gif3d")
    if knee == True:
        if wr == KneeWr:
            pvfig.screenshot(path + "Pts_iCSD_knee_wr" + str(KneeWr) + ".png")
            if gif3d == True:
                viewup = [0.5, 0.5, 1]
                path = pvfig.generate_orbital_path(
                    factor=2.0, shift=poly.length, viewup=viewup, n_points=36
                )
                pvfig.open_movie(path + "orbit.gif", framerate=4)
                pvfig.orbit_on_path(path, write_frames=True, viewup=[0, 0, 1])
        else:
            pvfig.screenshot(path + "Pts_iCSD_wr" + str(wr) + ".png")
            if gif3d == True:
                viewup = [0.5, 0.5, 1]
                path = pvfig.generate_orbital_path(
                    factor=2.0, shift=poly.length, viewup=viewup, n_points=36
                )
                pvfig.open_movie("orbit.gif", framerate=4)
                pvfig.orbit_on_path(path, write_frames=True, viewup=[0, 0, 1])

    pvfig.close()
    pass
