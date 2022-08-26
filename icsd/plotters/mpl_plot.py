# -*- coding: utf-8 -*-
"""
Created on Mon May 11 14:44:26 2020

@author: Benjamin
"""

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata

from icsd.gridder.mkgrid import mkGrid_XI_YI
from icsd.importers.read import *


def _fig_Interpolation_(ax, f, coord, data, **kwargs):
    """plot the interpolation of the VRTe current fractions"""

    coord_x, coord_y = parseCoord(coord, dim="2d")
    XI, YI = mkGrid_XI_YI(coord_x, coord_y)
    points = np.column_stack((coord_x, coord_y))

    grid = griddata(points, data, (XI, YI), method="linear")  # grid is a np array
    img = ax.imshow(
        grid,
        interpolation="spline16",
        extent=(min(coord_x), max(coord_x), min(coord_y), max(coord_y)),
        aspect="auto",
        origin="lower",
        cmap="jet",
    )

    # img = ax.matshow(grid,
    #            extent = (min (coord_x), max(coord_x), min(coord_y), max(coord_y)),
    #            aspect = 'auto', origin = 'lower', cmap= 'jet')

    # cbar = f.colorbar(img,ax=ax,
    #                     orientation='vertical',
    #                     shrink=0.6)
    cbar = plt.colorbar(img, ax=ax, orientation="vertical", shrink=0.6)

    if kwargs.get("clim") is not None:
        clim = kwargs.get("clim")
    else:
        clim = [min(data), max(data)]
    img.set_clim(vmin=clim[0], vmax=clim[1])

    if kwargs.get("lgd_label") is not None:
        cbar.set_label(kwargs.get("lgd_label"), labelpad=10)
    else:
        cbar.set_label("Fraction of Current Source", labelpad=10)


def _fig_RealSources_(ax, sc):
    """add known point sources if present"""
    if sc == None:
        return
    for s in sc:
        sx = float(s.split(",")[0])
        sy = float(s.split(",")[1])
        ax.plot(sx, sy, "ow", markersize=10, markeredgecolor="k")


def _fig_ReturnElec_(ax, retElec):
    """plot the return electrode"""
    if retElec == None:
        return
    print("reading return electrode: ", retElec)
    retElecx = float(retElec.split(",")[0])
    retElecy = float(retElec.split(",")[1])
    plt.plot(retElecx, retElecy, "sw", markersize=10)


def _fig_VRTe_(ax, coord, data_sol, **kwargs):
    """plot the VRTe current fractions"""
    coord_x, coord_y = parseCoord(coord, dim="2d")

    if kwargs.get("clim") is not None:
        clim = kwargs.get("clim")
        norm_z = (data_sol - clim[0]) / (clim[1] - clim[0])
    else:
        norm_z = (data_sol - min(data_sol)) / (max(data_sol) - min(data_sol))
    grey_cm = plt.cm.get_cmap("Greys")
    edgecolor_norm_z = grey_cm(norm_z)
    jet_cm = plt.cm.get_cmap("jet")
    facecolor_norm_z = jet_cm(norm_z)
    ax.scatter(
        coord_x,
        coord_y,
        facecolor=facecolor_norm_z,
        edgecolor=edgecolor_norm_z,
        cmap="jet",
    )


def _fig_Axis_Labels_(ax, title):
    ax.set_title(title)
    ax.set_ylabel("y [m]", fontsize=12)
    ax.set_xlabel("x [m]", fontsize=12)
    # axes = plt.gca()
    # axes.set_xlim([0,0.53])
    # axes.set_ylim([0,0.52])
    ax.tick_params(axis="both", which="major")
    plt.tight_layout()
    ax.set_aspect("equal")


def _fig_ModelParameter_mi_(coord, mpi):
    plt.plot(coord[mpi, 0], coord[mpi, 1], "sr", markersize=10)


def parseCoord(coord, dim):

    coord_x = coord[:, 0]
    coord_y = coord[:, 1]

    if dim == "3d":
        coord_z = coord[:, 2]
        return coord_x, coord_y, coord_z
    else:
        return coord_x, coord_y


def plotRemotes(path, dim, pltRemotes=False):

    if dim == "2d":
        print("not yet implemented")
        return
    else:
        if pltRemotes == True:
            RemLineNb, Injection, coordE, pointsE = load_geom(
                path
            )  # geometry file containing electrodes position includinf remotes
            ax.scatter(
                coordE[RemLineNb : RemLineNb + 2, 1],
                coordE[RemLineNb : RemLineNb + 2, 2],
                coordE[RemLineNb : RemLineNb + 2, 3],
                marker="v",
                color="black",
                s=60,
                label="Remotes",
            )
            ax.scatter(
                coordE[Injection, 1],
                coordE[Injection, 2],
                coordE[Injection, 3],
                marker="*",
                color="black",
                s=60,
                label="A",
            )
            ax.scatter(
                coordE[:RemLineNb, 1],
                coordE[:RemLineNb, 2],
                coordE[:RemLineNb, 3],
                marker="s",
                color="black",
                s=60,
                label="Velecs",
            )
            ax.view_init(azim=-101, elev=35)


#%% Specific plot functions for ICSD outputs


def plotPareto(wr, pareto_list_FitRes, pareto_list_RegRes, IdPtkneew, path):
    p, ax = plt.subplots()  #'L-curve'
    # ax.annotate('Wr=' + str(int(wr)), xy=(float(np.asarray(pareto_list_FitRes)[IdPtkneew]),
    #                                  float(np.asarray(pareto_list_RegRes)[IdPtkneew])),
    #                               xytext=(float(np.asarray(pareto_list_FitRes)[IdPtkneew])+max(pareto_list_FitRes)/3,
    #                                  float(np.asarray(pareto_list_RegRes)[IdPtkneew])+max(pareto_list_RegRes)/3),
    #                               arrowprops=dict(facecolor='black', shrink=0.05))
    plt.plot(
        float(np.asarray(pareto_list_FitRes)[IdPtkneew]),
        float(np.asarray(pareto_list_RegRes)[IdPtkneew]),
        "og",
    )
    plt.plot(pareto_list_FitRes, pareto_list_RegRes, "or")
    ax = plt.gca()
    ax.tick_params(axis="both", which="both", direction="out")
    ax.grid()
    plt.xlabel("Residual")
    plt.ylabel("Roughness")
    plt.tight_layout()
    plt.show()
    plt.savefig(path + "ParetoFront.png", dpi=600)

    return p, ax


def plot_knee_icsd(wr, kn):
    """Plot CSD for the best regularisation parameter after L-curve automatic analysis using a knee-locator

    Parameters
    ------------
    self
    """
    KneeWr = wr
    kn.plot_knee_normalized()
    # plt.show(kn)
    # run_single()


def plotFIT(b, b_w, xfun, path):
    stopat = len(b)  # 204
    plt.figure()
    plt.subplot(121)
    plt.plot(xfun[:stopat] + b_w[:stopat], "or", label="Inverted CSD")
    plt.plot(b_w[:stopat], "ob", label="True model")
    plt.xlabel("Measurement number")
    plt.ylabel("R [Ohm]")
    plt.legend()
    plt.subplot(122)
    plt.plot(xfun[:stopat] + b_w[:stopat], b_w[:stopat], "or")
    plt.xlabel("Inverted CSD, R [Ohm]")
    plt.ylabel("True model, R [Ohm]")
    plt.tight_layout()
    plt.savefig(path + "Fit.png", dpi=600)
    plt.show()


def plotCSD2d(
    coord, data_sol, b, b_w, path, pareto, retElec=None, sc=None, ax=None, **kwargs
):
    """Plot CSD in 2d, using matplotlib and scipy interpolation

    Parameters
    ------------
    self
    """

    clim = None
    if kwargs.get("clim") is not None:
        clim = kwargs.get("clim")

    if kwargs.get("index") is not None:
        fig_name = "CSD 2d T" + str(kwargs.get("index"))
    else:
        fig_name = "CSD 2d"

    if ax == None:
        f = plt.figure("CSD 2d")
        ax = plt.subplot()
        # ax = plt.gca()
    else:
        f = plt.figure("CSD 2d")
        ax = plt.subplot()

    _fig_Interpolation_(ax, f, coord, data_sol, clim=clim)
    _fig_VRTe_(ax, coord, data_sol, clim=clim)
    _fig_RealSources_(ax, sc)
    _fig_ReturnElec_(ax, retElec)

    if kwargs.get("title_wr") is not None:
        if kwargs.get("index") is not None:
            title = fig_name + r"  $\lambda$=" + str(kwargs.get("title_wr"))
        else:
            title = r"$\lambda$=" + str(kwargs.get("title_wr"))
        _fig_Axis_Labels_(ax, title)

    return f, ax

    if kwargs.get("xfun") is not None:
        xfun = kwargs.get("xfun")

    if not pareto:
        plotFIT(b, b_w, xfun, path)


def plotCSD3d(
    wr,
    coord,
    data,
    path,
    filename,
    knee,
    KneeWr,
    ax=None,
    title=None,
    pltRemotes=False,
    **kwargs
):
    """plot scattered 3d current sources density for a given regularisation weight wr
    (can be the knee location if pareto-curve mode is run)

    Parameters
    ----------
    sc: sources coordinates

    kwargs (to add) : 'sc' (plot source position (for a synthetic case experiment)

    """

    coord_x, coord_y, coord_z = parseCoord(coord, dim="3d")

    if ax == None:
        print("no Axis")
        # ax = plt.gca()
        f = plt.figure("volume")
        ax = f.gca(projection="3d")

    step = (max(coord_x) - min(coord_x)) / 10

    xlin = np.arange(min(coord_x), max(coord_x), step)
    ylin = np.arange(min(coord_y), max(coord_y), step)
    zlin = np.arange(min(coord_z), max(coord_z), step)
    # generate new grid
    X, Y, Z = np.meshgrid(xlin, ylin, zlin)

    data_2_plot = data
    # fig = plt.figure()
    # ax = fig.add_subplot(projection='3d')

    # fig = plt.figure()
    # ax = fig.add_subplot(projection='3d')

    # sc=ax.scatter(coord_x, coord_y, coord_z, c=data_2_plot, cmap ='coolwarm', s=data_2_plot*1e4,
    #           )
    sc = ax.scatter(
        coord_x,
        coord_y,
        coord_z,
        c=data_2_plot,
        cmap="coolwarm",
        s=10,
    )
    ax.view_init(azim=0, elev=90)

    # sc=ax.scatter(coord_x, coord_y, coord_z, c=data_2_plot, cmap ='coolwarm', s=data_2_plot*1e4,
    #           )
    # if self.clim is None:
    #     print('none clim')
    # sc.set_clim(self.clim)
    cbar = plt.colorbar(sc, ax=ax)
    # self.labels()
    # cbar.set_label(self.physLabel)
    # ax.set_zlim([-10,0])

    # label=kwargs.get('label_nor', 'normal')

    for key, value in kwargs.items():
        print("{0} = {1}".format(key, value))

        if key == "zlim":
            ax.set_zlim([kwargs.get(key)[0], kwargs.get(key)[1]])

    plotRemotes(path, dim="3d", pltRemotes=False)  # plot remotes and injection position

    if title == None:
        title = "Scattered current sources density, wr=" + str(wr)
    else:
        title = title + ", wr=" + str(wr)
    # plt.legend()
    ax.set_ylabel("y [m]", fontsize=12)
    ax.set_xlabel("x [m]", fontsize=12)
    ax.set_title(title)
    plt.savefig(
        path + filename + "_icsd_scatter.png",
        dpi=550,
        bbox_inches="tight",
        pad_inches=0,
    )

    if knee == True:
        if wr == KneeWr:
            plt.savefig(
                path + filename + "icsd_knee_scatter" + str(KneeWr) + ".png",
                dpi=550,
                bbox_inches="tight",
                pad_inches=0,
            )

    plt.show()

    return f, ax


#%% Generic plot functions


def scatter2d(coord, data, label, path, filename, pltRemotes=False, ax=None, **kwargs):

    coord_x, coord_y = parseCoord(coord, dim="2d")

    # f = plt.figure('volume')

    if ax == None:
        f = plt.figure("volume")
        print("ax = None")

    step = (max(coord_x) - min(coord_x)) / 10
    xlin = np.arange(min(coord_x), max(coord_x), step)
    ylin = np.arange(min(coord_y), max(coord_y), step)

    X, Y = np.meshgrid(xlin, ylin)
    ax = f.gca()
    sc = ax.scatter(coord_x, coord_y, c=data, cmap="coolwarm")
    cbar = plt.colorbar(sc)
    cbar.set_label(label)
    ax.set_ylabel("y [m]", fontsize=15)
    ax.set_xlabel("x [m]", fontsize=15)

    plotRemotes(path, dim="2d", pltRemotes=False)  # plot remotes and injection position

    plt.legend()
    # plt.title(title)
    plt.savefig(
        path + filename + "_icsd_scatter.png",
        dpi=550,
        bbox_inches="tight",
        pad_inches=0,
    )
    plt.show()

    return ax


def scatter3d(coord, data, label, path, filename, pltRemotes=False, ax=None, **kwargs):

    coord_x, coord_y, coord_z = parseCoord(coord, dim="3d")

    f = plt.figure("volume")

    if ax == None:
        print("no Axis")
        ax = plt.gca()
    # ax=f.gca(projection='3d')
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")

    step = (max(coord_x) - min(coord_x)) / 10
    xlin = np.arange(min(coord_x), max(coord_x), step)
    ylin = np.arange(min(coord_y), max(coord_y), step)
    zlin = np.arange(min(coord_z), max(coord_z), step)
    X, Y, Z = np.meshgrid(xlin, ylin, zlin)
    sc = ax.scatter(coord_x, coord_y, coord_z, c=data, cmap="coolwarm")
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label(label)
    ax.set_ylabel("y [m]", fontsize=15)
    ax.set_xlabel("x [m]", fontsize=15)

    plotRemotes(path, dim="3d", pltRemotes=False)  # plot remotes and injection position

    # plt.legend()
    # plt.title(title)
    plt.savefig(
        path + filename + "_icsd_scatter.png",
        dpi=550,
        bbox_inches="tight",
        pad_inches=0,
    )
    plt.show()

    return f


def labels(method):
    """Parse graphical labels to plot"""
    if method == "F1":
        physLabel = "normed misfit F1"
    if method == "Pearson":
        physLabel = "Pearson r coeff"

    return physLabel


def plotContour2d(
    coord, data_sol, physLabel, path, retElec=None, sc=None, ax=None, **kwargs
):
    """Plot contour in 2d, using matplotlib and scipy interpolation

    Parameters
    ------------
    self
    """

    if kwargs.get("fig_name") is not None:
        fig_name = kwargs.get("fig_name")
    else:
        fig_name = "2d"

    if kwargs.get("index") is not None:
        fig_name += "  T" + str(kwargs.get("index"))
        # print(fig_name)

    if ax == None:
        f = plt.figure(fig_name)
        ax = plt.gca()

    _fig_Interpolation_(ax, f, coord, data_sol, lgd_label=physLabel)
    _fig_VRTe_(ax, coord, data_sol)
    _fig_RealSources_(ax, sc)
    # _fig_ReturnElec_(ax,retElec)
    # ax.set_xlabel('x [m]')
    # ax.set_ylabel('y [m]')
    # ax.set_aspect('equal', adjustable='box')
    _fig_Axis_Labels_(ax, title=fig_name)

    if kwargs.get("jac") is not None:
        _fig_ModelParameter_mi_(coord, kwargs.get("jac"))

    if kwargs.get("title_wr") is not None:
        title = r"$\lambda$=" + str(title_wr)
        _fig_Axis_Labels_(title)

    return ax, f


def showObs2d(path, ax=None, **kwargs):
    """Plot contour in 2d, using matplotlib and scipy interpolation.
    Required surface and borehole electrode to make the 2d interpolation possible
    Parameters
    ------------
    self
    """
    filename = "ObsData.txt"
    clim = None
    if ax == None:
        f = plt.figure("ObsData")
        ax = plt.gca()

    if kwargs.get("filename") is not None:
        filename = kwargs.get("filename")
    if kwargs.get("clim") is not None:
        clim = kwargs.get("clim")
    if kwargs.get("index") is not None:
        index = kwargs.get("index")
    else:
        index = 0

    RemLineNb, Injection, coordE, pointsE = load_geom(
        path
    )  # geometry file containing electrodes position includinf remotes
    data_obs = load_obs(path, filename, index)

    _fig_Interpolation_(ax, pointsE, data_obs, lgd_label="U/I", clim=clim)
    _fig_VRTe_(ax, pointsE, data_obs)
    # _fig_RealSources_(ax,sc=None)
    # _fig_ReturnElec_(ax,retElec=None)
    _fig_Axis_Labels_(ax, title="Obs T=" + str(index))

    return ax


def current_streamlines(path, Res=1, mesh=None, **kwargs):
    """
    Current streamlines
    -------
    Only suitable for MALM with fix potential electrodes.
    Correction of resistivity Res

    """
    filename = "ObsData.txt"
    f = plt.figure("ObsData")
    ax = plt.gca()

    if kwargs.get("filename") is not None:
        filename = kwargs.get("filename")

    RemLineNb, Injection, coordE, pointsE = load_geom(
        path
    )  # geometry file containing electrodes position includinf remotes
    # print(path)
    # print(path)

    data_obs = load_obs(path, filename)

    xn = 30
    yn = 30

    xx = np.linspace(min(coordE[:, 1]), max(coordE[:, 1]), xn)
    yy = np.linspace(min(coordE[:, 2]), max(coordE[:, 2]), yn)

    xx, yy = np.meshgrid(xx, yy)
    points = np.transpose(np.vstack((coordE[:, 1], coordE[:, 2])))
    u_interp = interpolate.griddata(points, data_obs, (xx, yy), method="cubic")
    uu = np.reshape(u_interp, [xn * yn])
    data = uu / Res

    if mesh is None:
        mesh = pg.createGrid(
            x=np.linspace(min(coordE[:, 1]), max(coordE[:, 1]), xn),
            y=np.linspace(min(coordE[:, 2]), max(coordE[:, 2]), yn),
        )

    if kwargs.get("vmin"):
        vmin = kwargs.get("vmin")
    else:
        vmin = min(Obs)

    if kwargs.get("vmax"):
        vmax = kwargs.get("vmax")
    else:
        vmax = max(Obs)

    if kwargs.get("ax"):
        ax = kwargs.get("ax")
    else:
        fig, ax = plt.subplots()

    sc = ax.scatter(
        coordE[:, 1], coordE[:, 2], c=Obs, cmap="coolwarm", s=5e2, vmin=vmin, vmax=vmax
    )  # norm=matplotlib.colors.Normalize()
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label("V")
    ax.set_ylabel("y [m]", fontsize=15)
    ax.set_xlabel("x [m]", fontsize=15)

    # if len(kwargs.get('sensors'))>0:
    #     sensors = kwargs.get('sensors')
    #     ax.scatter(sensors[:,0],sensors[:,1],color='k',marker='.',label='pot. elec')
    #     for i in range(len(sensors[:,0])):
    #            ax.annotate(str(i+1), (sensors[i,0],sensors[i,1]))
    #     if kwargs.get('A'):
    #         A = kwargs.get('A')
    #         ax.scatter(sensors[A,0],sensors[A,1],color='r',marker='v',label='A. elec')
    #     if kwargs.get('B'):
    #         B = kwargs.get('B')
    #         ax.scatter(sensors[B,0],sensors[B,1],color='b',marker='v',label='B. elec')
    #     if kwargs.get('Nfix'):
    #         Nfix = kwargs.get('Nfix')
    #         ax.scatter(sensors[Nfix,0],sensors[Nfix,1],color='g',marker='v',label='Nfix. elec')

    # if kwargs.get('gridCoarse'):
    #     gridCoarse = pg.createGrid(x=np.linspace(min(sensors[:,0]), max(sensors[:,0]),xn/2),
    #                      y=np.linspace(min(sensors[:,1]), max(sensors[:,1]),yn/2))

    drawStreams(
        ax, mesh, -pg.solver.grad(mesh, data), color="green", quiver=True, linewidth=3.0
    )
