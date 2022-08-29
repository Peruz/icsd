import os
from copy import copy, deepcopy

import matplotlib.pyplot as plt
import numpy as np
#from kneed import KneeLocator
from matplotlib.backends.backend_pdf import PdfPages
from scipy.linalg import block_diag
from scipy.optimize import least_squares, lsq_linear

from .exporters.save import *
from .gridder.mkgrid import mkGrid_XI_YI
from .importers.read import *
from .inversion.priorM0 import *
from .inversion.regularize import *
from .inversion.solve import *
from .plotters.mpl_plot import (
    labels,
    plot_knee_icsd,
    plotContour2d,
    plotCSD2d,
    plotCSD3d,
    plotPareto,
    scatter2d,
    scatter3d,
)
from .plotters.pv_plot import plotCSD3d_pv


class iCSD3d(object):

    """
    Create a icsd inversion object.
    """

    def __init__(self, dirName):

        self.surveys = []
        self.df = (
            None  # main dataframe containing resistance values (TDIP data in progress)
        )

        # input files directory
        self.dirName = dirName

        # survey files to import (default names)
        self.type = "2d"  # or '3d'
        self.sim = "VRTeSim.txt"
        self.obs = "ObsData.txt"
        self.coord_file = "VRTeCoord.txt"  # coordinates of the virtual sources

        # inversion default parameters
        self.wr = 25  # weight regularization
        self.wc = 1e5  # 1e5 #current conservation constrain, sum current fractions = 1
        self.pareto = False  # if True run many icsd to explore the Pareto front
        self.pareto_MinErr = 0.001
        self.pareto_MaxErr = 1
        self.pareto_nSteps = 10
        self.k = 4  # For each point, find the k closest current sources
        self.x0_prior = False  #  relative smallness regularization as a prior criterion for the inversion; i.e the algorithm minimizes ||m−m0||2
        self.x0_ini_guess = False  # initial guess
        self.knee = False  # L-curve knee automatic detection
        self.KneeWr = []

        self.logTrans = False
        # time lapse option
        # self.iTimeLapse = False # to enable timelapse inversion [Not implemented yet]
        self.TL = False  # Time lapse inversion (see slides VIU Venice to implement)
        self.TL_gamma = 1  # Time lapse inversion weight

        # tdip option
        self.TDIP_flag = False  # import tdip data
        # self.sequential=False #if True run a sequential inversion using previsous gate as reference model m0

        # mesh type
        self.regMesh = "strc"  # strc or unstrc

        # Initial model
        self.alphax0 = 1  # weight on model smallness relative to m0
        self.inix0 = (
            "cst"  #  intial physical model: None or 'cst' if constant vector *0.1
        )
        self.method_m0 = "F1"  # method to estimate intial physical model (F1= misfit from a single current source OR ProductMoment method)

        # spatial regularization
        self.alphaSxy = False  # weight on model smoothness in z-direction
        self.alphaSx = 1  # weight on model smoothness in x-direction
        self.alphaSy = 1  # weight on model smoothness in y-direction
        self.alphaSz = 1  # weight on model smoothness in z-direction [TO IMPLEMENT]

        # data regularisation
        self.obs_err = "sqrt"  # const or sqrt or reciprocal [not yet implemented] - choose between constant weight and w = 1/sqrt(abs(obs))
        self.errRmin = (
            1  # min R to which the err is applied before passing to constant error
        )

        # graphic options
        self.sc = []  # coordinates of the sources, format = x1,y1 x2,y2'
        self.retElec = None  # coordinates of the return electrode, format = x1,y1')
        self.clim = []
        self.plotElecs = False
        self.gif3d = False  # create gif orbit
        self.title = None  #
        self.mesh = None  # mesh3d .vtk to plot with the results of icsd

        # IMPLEMENT obs_err based on reciprocal analysis i.e. estimated standard deviation of the data errors;                         % estimated standard deviation of the traveltime data errors
        # self.solution=None
        # processing outputs [To write]
        # self.models = []
        # self.rmses = []

        if self.surveys:
            print("Existing survey")

    def icsd_init(self, survey):
        """these functions are called only once, even for pareto,
        as they do not depend on the regularization weight"""
        # create directories

        print("initiation ICSD")
        self.createdirs(survey)

        # load virtual sources coordinates
        if self.type == "2d":
            self.coord_x, self.coord_y, survey.coord = load_coord(
                survey.path2load, self.coord_file, dim=2,
            )
        else:
            self.coord_x, self.coord_y, self.coord_z, survey.coord = load_coord(
                survey.path2load, self.coord_file, dim=3,
            )

        if self.TDIP_flag == False:
            # load observations resistances b
            survey.b = load_obs(survey.path2load, survey.obs)
            # load simulated resistances A (i.e. Green function)
            survey.A = load_sim(survey.path2load, survey.sim)

        print("log transformation: " + str(self.logTrans))
        # Log transformation before inversion
        if self.logTrans == True:
            # check the min

            # IF data and model are in the same range
            # TranslateMIN= np.min(np.array([min(survey.b),min(survey.A)]))
            # survey.A = np.log(survey.A + 1 - TranslateMIN) # translate, then transform */
            # survey.b = np.log(survey.b + 1 - TranslateMIN) # translate, then transform */

            # IF data is very small compare to the model

            TranslateMIN_A = min(survey.A)
            survey.A = np.log(
                survey.A + 1 - TranslateMIN_A
            )  # translate, then transform */
            TranslateMIN_B = min(survey.b)
            survey.b = np.log(
                survey.b + 1 - TranslateMIN_B
            )  # translate, then transform */

            print("TO DO need to log-trans data errors as well")
            # TO DO need to log-trans data errors as well

        # load observations electrode coordinates
        if self.plotElecs == True:
            (
                survey.RemLineNb,
                survey.Injection,
                survey.coordE,
                survey.pointsE,
            ) = load_geom(
                self.path2load
            )  # geometry file containing electrodes position includinf remotes

        # check vector sizes
        survey.nVRTe = check_nVRTe(survey.A, survey.b, survey.coord)

        # reshape A vector into matrix A of nVRTE collumns
        survey.A = reshape_A(survey.A, survey.nVRTe)

        print("obs_err: " + str(self.obs_err))
        # define mode to weights the data (const, sqrt or reciprocal)
        survey.obs_w = obs_w_f(self.obs_err, survey.b, self.errRmin, sd_rec=None)

        # set constrain (curent conservation)
        survey.con_A = con_A_f(survey.A)
        survey.con_b = con_b_f(survey.b)
        survey.con_w = con_w_f(self.wc)

        # append spatial regularization (add lines to the matrice)
        survey.reg_A = self._parseModelReg(survey)
        survey.reg_b = regularize_b(self.reg_A)

        # stack data, constrain, and regularization
        survey.A_s = stack_A(survey.A, survey.con_A, survey.reg_A)
        survey.b_s = stack_b(survey.b, survey.con_b, survey.reg_b)

    ### mkdirs
    def createdirs(self, survey):
        survey.path2load = self.dirName
        survey.path2save = survey.path2load + "figs"
        # self.path2save= self.dirName + 'fig/'
        # print(self.path2save)
        # cwd = os.getcwd()
        # print(cwd+self.path2save)
        # try:
        #     # Create target Directory
        #     os.mkdir(self.path2save)
        #     print("Directory " , self.path2save ,  " Created ")
        # except FileExistsError:
        #     print("Directory " , self.path2save ,  " already exists")

    ### PREPARE FOR ICSD

    def prepare4iCSD(self, index=0):
        """this function is called for each weight, keep them separated for pareto
        # 1 - Estimate M0 initial if smallness guided inversion or initial guess
        # 2-  Create vector with weights-->
        (data weigth, constrainsts weight and regularisation weigth)

        """
        if (
            self.surveys[index].x0_ini_guess == True
            or self.surveys[index].x0_prior == True
        ):
            self._estimateM0_(index=index)

        # Create vector with weight (data weigth, constrainsts weight and regularisation weigth)
        if self.surveys[index].x0_prior == True:  # if relative smallness
            self.reg_w_0_b, self.reg_w_0_A = regularize_w(
                self.surveys[index].reg_A,
                self.wr,
                self.surveys[index].x0_prior,
                x0=self.x0,
            )

            # stack weight matrix
            self.W_s_A, self.W_s_b = stack_w(
                self.surveys[index].obs_w,
                self.surveys[index].con_w,
                self.surveys[index].x0_prior,
                reg_w_0_A=self.reg_w_0_A,
                reg_w_0_b=self.reg_w_0_b,
            )

            # apply weight to A and b (data weigth, constrainsts weight and regularisation weigth)
            self.A_w = weight_A(
                self.surveys[index].x0_prior, self.surveys[index].A_s, W_s_A=self.W_s_A
            )
            self.b_w = weight_b(
                self.surveys[index].x0_prior, self.surveys[index].b_s, W_s_b=self.W_s_b
            )

        else:
            self.reg_w = regularize_w(self.surveys[index].reg_A, self.wr, self.x0_prior)

            self.W_s = stack_w(
                self.surveys[index].obs_w,
                self.surveys[index].con_w,
                self.x0_prior,
                reg_w=self.reg_w,
            )
            self.A_w = weight_A(self.x0_prior, self.surveys[index].A_s, W_s=self.W_s)
            self.b_w = weight_b(self.x0_prior, self.surveys[index].b_s, W_s=self.W_s)

    def _parseModelReg(self, survey):
        """Parse regularisation smoothing and
        prior constrainst parameters before inversion
        """

        # 2D CASE -----------------------------------------
        
        print(self.type == "2d")
        if self.type == "2d":
            if self.regMesh == "strc":  # structured (grid) mesh of virtual sources
                if (
                    self.alphaSxy == True
                ):  # spatial smoothing with coeff  alphaSx and alphaSy
                    self.reg_Ax, self.reg_Ay = regularize_A_x_y(
                        survey.coord, self.alphaSx, self.alphaSy
                    )
                else:
                    if self.x0_prior == True:
                        raise ValueError(
                            "#### dimensions of matrices do not agree - change regularisation types"
                        )
                    else:
                        print('Im here')
                        self.reg_A = regularize_A(survey.coord, survey.nVRTe)
            else:
                self.reg_A = regularize_A_UnstructuredMesh3d(
                    survey.coord, self.nVRTe, self.k
                )

        # 3D CASE -----------------------------------------
        else:
            if self.regMesh == "strc":
                if self.x0_prior == True:
                    # self.regularize_A_x_y_z()
                    raise ValueError(
                        "### Not yet ready to work with m0 and 3d grid mesh, regularize_A_x_y_z need to be tested"
                    )
                else:
                    self.reg_A = regularize_A_3d(
                        self.nVRTe, survey.coord
                    )  # working for structured mesh

            elif self.regMesh == "unstrc":
                self.reg_A = regularize_A_UnstructuredMesh3d(
                    survey.coord, self.nVRTe, self.k
                )

            if self.alphaSxy == True:  # anisotropic smoothing
                # print(self.alphax0)
                self.reg_smallx0 = ponderate_smallnessX0(
                    self.alphaSxy, self.alphax0, reg_Ax=self.reg_Ax
                )
                self.reg_A = sum_smallness_smoothness(
                    self.alphaSxy,
                    self.x0_prior,
                    reg_smallx0=self.reg_smallx0,
                    reg_Ax=self.reg_Ax,
                    reg_Ay=self.reg_Ay,
                )
            else:
                self.reg_smallx0 = ponderate_smallnessX0(
                    self.alphaSxy, self.alphax0, reg_A=self.reg_A
                )
                self.reg_A = sum_smallness_smoothness(
                    self.alphaSxy,
                    self.x0_prior,
                    reg_smallx0=self.reg_smallx0,
                    reg_A=self.reg_A,
                )

        return self.reg_A

    def run_all_gates(self, showfig=False):
        """stack all the IP gates and invert simultaneously # Not tested"""
        # stack A and b

        self.prepare4iCSD(index=0)
        A_all_gates = self.A_w
        b_all_gates = self.b_w

        idx = np.hstack([0, np.arange(1, len(self.surveys), 2)])

        # for i in range(1,len(self.surveys),2):
        for i in idx:
            # gateid = [1,2,3,8,16,20]
            # gateid = [1,2]
            # for i in gateid:
            # print(i)
            self.prepare4iCSD(index=i)
            A_all_gates = block_diag(A_all_gates, self.A_w)
            b_all_gates = np.hstack([b_all_gates, self.b_w])

        # print(np.shape(A_all_gates))
        # print(np.shape(b_all_gates))

        # constrainsted inversion
        # if (self.x0_ini_guess == True or self.x0_prior == True):
        #     self.x = iCSD(self.A_w,self.b_w,
        #                   self.type,
        #                   self.surveys[index].coord,
        #                   self.surveys[index].path2load,
        #                   x0=self.x0)

        # UNconstrainsted inversion
        # else:
        xsol = iCSD(
            A_all_gates,
            b_all_gates,
            self.type,
            self.surveys[0].coord,
            self.surveys[0].path2load,
        )

        # split system to recover the solution for each gate
        i_s = 0
        i_end = np.shape(self.A_w)[1]

        # xsol_new = {}

        # for i in range(0,3):

        #     xsol_new[i] = xsol
        #     xsol_new[i].x_split = xsol.x[i_s:i_end]
        #     i_s = np.shape(self.A_w)[1]*(i+1)
        #     i_end = np.shape(self.A_w)[1]*(i+2)

        # sol_backup = deepcopy(self.surveys[0])
        xsol_backup = deepcopy(xsol)
        # sol_backup.tmp = xsol
        # sol_backup.tmp.xsol = xsol.x
        # print(xsol)
        # print(self.surveys[i].solution)

        # for i in range(1,len(self.surveys),2):
        for i in idx:
            print(i)
            xsol_backup = deepcopy(xsol)
            xsol_backup.x = deepcopy(xsol.x)
            # print(len(xsol.x))
            # print(len(xsol_backup.x))
            # print(self.surveys[i])
            setattr(self.surveys[i], "solution", xsol_backup)
            setattr(self.surveys[i].solution, "x", xsol_backup.x[i_s:i_end])
            # print(xsol_backup.x)
            # print(len(xsol_backup.x))
            # print(i_s,i_end)
            # print(xsol_backup.x[i_s:i_end])
            # print(self.surveys[i].solution.xtest)
            i_s = np.shape(self.A_w)[1] * (i + 1)
            i_end = np.shape(self.A_w)[1] * (i + 2)
            # del sol_backup
            print(self.surveys[i].solution.x)

        # for i in range(0,3):
        #     self.surveys[i].solution.x = sol_backup[i].solution.x
        #     # self.surveys[i].solution.x=sol_backup.x[:i_end]
        #     self.surveys[i].solution.x = sol_backup[i].solution.x[i_s:i_end]
        #     print(self.surveys[i].solution.x)
        #     print(i)
        #     print(sol_backup.x)

        # i_s = np.shape(self.A_w)[1]*(i+1)
        # i_end = np.shape(self.A_w)[1]*(i+2)
        # del sol_backup

    def run_single(self, index=0, show=False):
        """Run a single inversion (unique regularisation weight)
        Equivalent to several steps::
            self.prepare4iCSD()
            self.plotCSD()
            self.misfit()

        """
        print("run_single i=" + str(index))
        if not hasattr(self.surveys[index], "A"):  # if A matrix does not exist
            self.icsd_init(self.surveys[index])

        self.prepare4iCSD(index=index)

        # constrainsted inversion
        if (
            self.surveys[index].x0_ini_guess == True
            or self.surveys[index].x0_prior == True
        ):
            print("constrainsted inversion")
            self.x = iCSD(
                self.A_w,
                self.b_w,
                self.type,
                self.surveys[index].coord,
                self.surveys[index].path2load,
                x0=self.x0,
            )

        # UNconstrainsted inversion
        else:
            print("UNconstrainsted inversion")
            self.x = iCSD(
                self.A_w,
                self.b_w,
                self.type,
                self.surveys[index].coord,
                self.surveys[index].path2load,
            )

        self.surveys[index].solution = self.x

        if show == True:
            ax, f = self.showResults(index=index)
            # self.RMSAnalysis()
            self.misfit()
            # f.savefig(self.path2save+'iCSD', dpi = 600)
            plt.tight_layout()
            plt.show()
            plt.close(f)
            # add result to container
            # _add_to_container(self.x)

            # plt.close(f)
            # self.ModelResolution()
            return ax, f
        else:
            pass

    def run_pareto(self):
        """
        run iCSD multiple times while changing the weights to explore the L-curve
        """
        # self.icsd_init(self.surveys[0])
        if hasattr(self, "pareto_weights") is False:
            self.pareto_weights = np.linspace(
                self.pareto_MinErr, self.pareto_MaxErr, self.pareto_nSteps
            )

        print("pareto weights are\n", self.pareto_weights)

        self.pareto_list_FitRes = []
        self.pareto_list_RegRes = []

        savename = self.surveys[0].path2save + "iCSD_pareto.pdf"
        if self.x0_prior:
            savename = self.surveys[0].path2save + "_xprior_iCSD_pareto.pdf"

        with PdfPages(savename) as pdf:

            # LOOP ON REG WEIGHTS
            for self.wr in self.pareto_weights:

                # RUN ICSD
                self.prepare4iCSD()
                if self.x0_ini_guess == True or self.x0_prior == True:
                    self.x = iCSD(
                        self.A_w,
                        self.b_w,
                        self.type,
                        self.surveys[0].coord,
                        self.surveys[0].path2load,
                        x0=self.x0,
                    )
                else:
                    self.x = iCSD(
                        self.A_w,
                        self.b_w,
                        self.type,
                        self.surveys[0].coord,
                        self.surveys[0].path2load,
                    )

                # PLOT ICSD
                if self.type == "2d":
                    self.f, _ = plotCSD2d(
                        self.surveys[0].coord,
                        self.x.x,
                        self.surveys[0].b,
                        self.b_w,
                        self.surveys[0].path2load,
                        self.pareto,
                        retElec=None,
                        sc=None,
                        title_wr=self.wr,
                    )
                else:
                    # pass
                    # plotCSD3d(wr,coord,data,path,filename,knee,KneeWr,ax=None,title=None,pltRemotes=False,**kwargs)
                    self.f, _ = plotCSD3d(
                        self.wr,
                        self.surveys[0].coord,
                        self.x.x,
                        self.surveys[0].path2load,
                        self.obs,
                        self.knee,
                        self.KneeWr,
                        title=None,
                    )

                pdf.savefig(self.f)
                # print(self.f)
                # plt.show()
                # plt.close(self.f)

                self.residualAnalysis()
                self.misfit()

            try:
                # Detect knee point wr of the L-curve
                self.detectKneePt()
            except:
                print("Can" "t find a knee")

            self.wr = float(self.pareto_weights[self.IdPtkneew])
            print("Knee detected for wr=" + str(self.wr))

            # Plot the L-curve
            self.p, self.ax = plotPareto(
                self.wr,
                self.pareto_list_FitRes,
                self.pareto_list_RegRes,
                self.IdPtkneew,
                self.surveys[0].path2save,
            )
            pdf.savefig(self.p)
            # plt.close(self.p)

            # Plot knee curve
            if self.knee:
                plot_knee_icsd(self.wr, self.kn)
                self.run_single()

        pass

    def detectKneePt(self):
        self.kn = KneeLocator(
            self.pareto_list_FitRes,
            self.pareto_list_RegRes,
            curve="convex",
            direction="decreasing",
        )
        self.IdPtkneew = np.where(self.kn.knee == self.pareto_list_FitRes)[0]
        self.pareto_weights[self.IdPtkneew]
        if len(self.IdPtkneew) < 1:
            self.IdPtkneew = 1
            print("No knee detection possible, put 1 as default")
        #            raise ValueError('No knee detection possible')
        else:
            self.knee = True

    def invert(self, show=False, **kwargs):
        """Invert voltage to current densities.

        Parameters
        ----------
        * x0_prior:
            If True, relative smallness x0_prior model
        * regMesh:'unstrc' or 'strc'
            Specify mesh type
        * wr:
            Specify regularisation weight
        * pareto :
            If True, run L-curve analysis.
        * pareto_MinErr, pareto_MaxErr, pareto_nSteps:
            Specify pareto paramters

        Returns:

        * x :
            pd.dataframe
            Data container
        """
        for (k, v) in kwargs.items():
            setattr(self, k, v)
            print("{0} = {1}".format(k, v))

        if self.pareto == False:
            if len(self.surveys) < 2:
                if show:
                    ax, f = self.run_single(show=show)
                    return self.x, ax, f
                else:
                    self.run_single(show=show)
                    return self.x

            else:
                if self.TL == True:  # flag for time-lapse inversion
                    if self.TL_gamma == 1:
                        for i, survey in enumerate(self.surveys):
                            self.run_single(index=i)
                    else:
                        raise ValueError("weigted TL inversion not yet implemented")
                        # self.run_single_TL()
                else:  # TDIP survey
                    if "index" in kwargs:
                        self.run_single(index=kwargs["index"])
                    elif (
                        "simultaneous" in kwargs
                    ):  # stack all the gates and invert simultaneously
                        self.run_all_gates()
                    elif (
                        "sequential" in kwargs
                    ):  # use previous time steps/ IP gates as initial model m0
                        for i, survey in enumerate(self.surveys):
                            if i == 0:
                                self.run_single(index=i)
                            else:
                                # self.sequential = True
                                self.surveys[i].x0_ini_guess = True
                                self.surveys[i].x0_prior = True
                                self.x0 = self.surveys[0].solution.x
                                self.run_single(index=i, show=show)
                    else:
                        for i, survey in enumerate(self.surveys):
                            self.run_single(index=i, show=show)

        else:
            self.run_pareto()

        # return self.surveys

    #%% DEFINE SURVEY container for observations and simulations

    def createSurvey(
        self, fname_obs="ObsData.txt", fname_sim="VRTeSim.txt", append=False
    ):
        """
        Parameters
        ----------
        fname_obs : TYPE
            DESCRIPTION.
        fname_sim : TYPE
            DESCRIPTION.
        append : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        None.

        """
        # set attribute according to the first survey
        if len(self.surveys) == 0:
            print("no existing survey")
            # survey = iCSD3d(dirName=self.dirName)
            survey = deepcopy(self)
            survey.obs = fname_obs
            survey.sim = fname_sim
            survey.icsd_init(survey)
            print(survey.path2load)
            self.surveys.append(survey)
        else:  # check we have the same configuration than other survey
            if append:
                print("append new survey")
                # survey = iCSD3d(dirName=self.dirName)
                survey = deepcopy(self)
                survey.obs = fname_obs
                survey.sim = fname_sim
                survey.icsd_init(survey)
                # check = [a == b for a,b, in zip(self.nVRTe, self.surveys[len(surveys-1)].nVRTe)]
                # if all(check) is True:
                self.surveys.append(survey)
            else:
                print("Replace existing survey")
                survey = deepcopy(self)
                survey.obs = fname_obs
                survey.sim = fname_sim
                survey.icsd_init(survey)
                self.surveys[0] = survey
        return self.surveys

    def createTimeLapseSurvey(self, fnames_obs, fnames_sim):
        """Import multiple surveys and return
        a surveys object.

        Parameters
        ----------
        fnames_obs : list of str
            List of file to be parsed or directory where the files are.
        """
        self.surveys = []
        if isinstance(fnames_obs, list):  # it's a list of filename
            if len(fnames_obs) < 2:
                raise ValueError("at least two files needed for timelapse inversion")
        else:  # it's a directory and we import all the files inside
            if os.path.isdir(fnames_obs):
                fnames_obs = [
                    os.path.join(fnames_obs, f)
                    for f in np.sort(os.listdir(fnames_obs))
                    if f[0] != "."
                ]
                # this filter out hidden file as well
            else:
                raise ValueError(
                    "dirname should be a directory path or a list of filenames"
                )
        for fnames_obs, fnames_sim in zip(fnames_obs, fnames_sim):
            self.createSurvey(fnames_obs, fnames_sim, append=True)

        return self.surveys

    def createTDIPSurvey(self, fname_obs, fname_sim, **kwargs):
        """create TDIP survey and return a survey object.

        Parameters
        ----------
        fname : *.bin file containing TDIP infos
        """

        Vp_norm = True  # Apply K geometric factor
        for key, value in kwargs.items():
            if key == "knorm":
                Vp_norm = value

        print("createTDIPSurvey")
        # read TDIP file
        # normalisation of all Vs by K factor
        self.Vs, A = loadTDIPSurvey(
            self.dirName + fname_obs, self.dirName + fname_sim, Vp_norm=Vp_norm
        )  # *.data (pygimli format)
        # tdip_sim = load_sim(self.dirName, fname_sim) #

        # self.Vs = np.vstack(tdip_obs)
        # self.Vs = np.transpose(self.Vs)

        # For now, transform each gates to several surveys as it is done for timelapse icsd
        self.surveys = []
        print("relaxing CC constrainst for time domain IP")
        for i in range(self.Vs.shape[1]):
            # survey = iCSD3d(dirName=self.dirName)
            if i > 0:
                self.wc = 1e-5  # relaxing CC constrainst for time domain IP
                # self.errRmin=1e2
            survey = copy(self)
            print(survey.wc)

            if A.ndim > 1:
                survey.A = A[:, i]
            else:
                survey.A = A
            survey.TDIP_flag = True  # activate tdip flag
            survey.b = self.Vs[:, i]
            survey.icsd_init(survey)
            self.surveys.append(survey)

        return self.surveys

    #%% INITIAL MODEL

    def estimateM0(self, method_m0="F1", show=True, ax=None, index=None):
        """Estimate initial M0 model for constrainsted inversion

        Parameters
        ----------
        * method_m0: 'F1' OR 'Pearson'
            Specify the method to estimate the initial M0
        * show:
            Show M0 misfit
        * ax:
            Specify plot axis
        Returns:

        """
        if index is not None:
            m0 = self._parseM0_(
                method_m0, index=index
            )  # define the method to estimate M0
            self.surveys[index].physLabel = labels(
                method_m0
            )  # lgd labeling for the plot
        else:
            for i, survey in enumerate(self.surveys):
                m0 = self._parseM0_(
                    method_m0, index=i
                )  # define the method to estimate M0
                self.surveys[i].physLabel = labels(
                    method_m0
                )  # lgd labeling for the plot
                if show == True:
                    ax, f = self.showEstimateM0(index=i, ax=ax)
                    return m0, ax, f

        return m0

    def _estimateM0_(self, show=True, index=None):
        self.estimateM0(method_m0=self.method_m0, index=index)

    def _parseM0_(self, method_m0, index=0):
        """Parse M0 parameters"""
        # for i, survey in enumerate(self.surveys):
        #     print(i)
        survey = self.surveys[index]
        # print(survey.b[0])
        if method_m0 is not None:
            if method_m0 == "F1":
                survey.path2load = self.dirName
                (
                    self.surveys[index].norm_F1,
                    self.surveys[index].x0,
                ) = misfitF1_2_initialX0(survey.A, survey.b)
            elif method_m0 == "Pearson":
                self.surveys[index].x0 = product_moment(survey.A, survey.b)
        elif self.inix0 is not None:
            if survey.inix0 == "cst":
                self.surveys[index].x0 = np.ones(survey.b.shape) * 0.1
        self.x0 = np.copy(self.surveys[index].x0)
        return self.surveys[index].x0

    #%% PLOT fcts

    def showEstimateM0(self, index=0, ax=None):
        """Show initial model m0 estimation.

        Parameters
        ----------
        index : int, optional
            Index of the survey to plot.
        ax : Matplotlib.Axes, optional
            If specified, the graph will be plotted against this axis.
        """

        # fig, ax = plt.subplots()
        if self.type == "2d":
            ax, f = plotContour2d(
                self.surveys[index].coord,
                self.surveys[index].x0,
                self.surveys[index].physLabel,
                path=self.surveys[index].path2load,
                retElec=None,
                sc=None,
                index=index,
                ax=ax,
                fig_name="Estimated m0",
            )
        else:
            ax, f = scatter3d(
                self.surveys[index].coord,
                self.surveys[index].x0,
                self.surveys[index].physLabel,
                self.surveys[index].path2load,
                self.surveys[index].obs,
                ax=ax,
            )
        plt.tight_layout()
        plt.show()
        # plt.close(f)

        return ax, f

    def showResults(
        self,
        ax=None,
        data=None,
        clim=None,
        cmap="viridis_r",
        plotElecs=False,
        sc=False,
        retElec=None,
        mesh=None,
        gif3d=False,
        title=None,
        cut=False,
        **kwargs
    ):
        """Show inverted model.

        Parameters
        ----------
        ax : Matplotlib.Axes, optional
            If specified, the graph will be plotted against this axis.
        clim : array, optional
            Minimum and Maximum value of the colorbar.
        cmap : str, optional
            Name of the Matplotlib colormap to use.
        plotElecs : bool, optional
            If `True` add to the ICSD plot measuring electrodes as points
        sc : bool, optional
            If `True` add to the ICSD plot remotes and injection electrodes as points
        retElec : array, optional
            Coordinates of the return electrode, format = x1,y1')
            If Not None add to the ICSD plot the return B electrode
        mesh :  str, optional
            Specify name of the vtk file
            If Not None add mesh3d.vtk to plot with the results of icsd (for 3d using pyvista)
        gif3d :  bool, optional
            If `True` record a gif using orbital function of pyvista
          title :  str, optional
            Specify inversion titlename to be add to the plot
        """
        self.clim = clim
        self.plotElecs = plotElecs
        self.sc = sc
        self.retElec = retElec  # return electrode
        self.mesh_over = mesh  # plot mesh over
        self.gif3d = gif3d  # animated
        self.title = title  # plot title

        index = 0
        for key, value in kwargs.items():
            if key == "index":
                index = value

        if data is None:
            data = self.surveys[index].solution.x

        if self.type == "2d":
            f, ax = plotCSD2d(
                self.surveys[index].coord,
                data,
                self.surveys[index].b,
                self.b_w,
                self.surveys[index].path2load,
                self.pareto,
                retElec=None,
                sc=None,
                ax=ax,
                title_wr=self.wr,
                index=index,
                clim=self.clim,
            )
        else:
            f = plotCSD3d(
                self.wr,
                self.surveys[index].coord,
                data,
                self.surveys[index].path2load,
                self.surveys[index].obs,
                self.surveys[index].knee,
                self.surveys[index].KneeWr,
                ax=ax,
                title=None,
                **kwargs
            )
            if cut == True:
                plotCSD3d_pyvistaSlice()
            else:
                plotCSD3d_pv(
                    self.surveys[index].coord,
                    path=self.surveys[index].path2load,
                    filename_root="Solution.dat",
                    knee=self.surveys[index].knee,
                    wr=self.surveys[index].wr,
                    KneeWr=self.surveys[index].KneeWr,
                    mesh=self.mesh_over,
                    plotElecs=plotElecs,
                    gif3d=self.surveys[index].gif3d,
                    **kwargs
                )
        return ax, f

    def showCurrentStreamlines(self, index=0, ax=None):
        """Show Current Streamlines.

        Parameters
        ----------
        """
        if self.type == "2d":
            print("not yet implemented in 3d")
            # f = current_streamlines(path)
        else:
            print("not yet implemented in 3d")

    #%% POST inversion analysis

    def residualAnalysis(self):
        fitting_res = self.x.fun[0 : self.surveys[0].b.shape[0]]
        # constrain_res = self.x.fun[self.b.shape[0] + 1] / self.wc
        regularization_res = (
            self.x.fun[self.surveys[0].b.shape[0] + 2 :] / self.wr
        )  # constrain not included in the reg function
        self.reg_sum = np.sum(np.square(regularization_res))
        self.fit_sum = np.sum(np.square(fitting_res))
        self.pareto_list_FitRes.append(self.fit_sum)
        self.pareto_list_RegRes.append(self.reg_sum)

    def RMSAnalysis(self):
        self.rms = np.sum(np.square(self.x.fun)) / len(self.x.fun)
        print("RMS error=" + str(self.rms))

    def misfit(self):
        """

        Parameters
        ----------

        """
        residuals = self.x.fun[0 : self.surveys[0].b.shape[0]] - self.surveys[0].b
        val = np.linalg.norm(self.surveys[0].obs_w * residuals)

        # val = np.sum(self.obs_w*(residuals**2))

        misfit = val * self.wr
        # print(self.b_w[0 : self.b.shape[0]])
        # print('RMS error='+ str(misfit))

        return misfit

    def modelResolution(self, sol=None, coord=None, jacMi=None):

        if sol is None:
            sol = self.x

        if coord is None:
            coord = self.coord

        J = sol.jac[self.b.shape[0] + 1 :, jacMi]

        # writeFIT(self.path2load,self.x0_prior,self.x.fun,self.b_w,self.b_s,self.x0,self.reg_A,J,saveall=False)

        if self.type == "2d":
            f = plotContour2d(
                coord,
                J,
                "Jac matrice",
                path=self.path2load,
                retElec=None,
                sc=None,
                jac=jacMi,
            )
        else:
            f = scatter3d(coord, J, "Jac matrice", self.path2load, self.obs)
        plt.tight_layout()
        plt.show()
        plt.close(f)


#%% TO DO

# NEW FUNCTIONS to introduce
# function penalized some virtual sources (split rhizotron)
# load_names for tl analysis
# Introduce error weigting using reciprocal error instead of constant or w = 1/sqrt(abs(obs))
# use the ERT mesh as mesh for virtual position of current sources
