# Copyright (c) 2018 The ICSD Developers.
# https://github.com/Peruz/icsd/graphs/contributors
# Distributed under the terms of the BSD 3-Clause License.
# SPDX-License-Identifier: BSD-3-Clause
#
"""
Initial model estimations functions
"""

import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as LA
from scipy.stats import pearsonr


### Individual Misfit
def _normF1(A, b, **kwargs):
    """compute the norm between observation data and individual green functions"""
    F1 = []
    for i in range(np.shape(A)[1]):
        F1i = LA.norm(
            (b - A[:, i])
        )  # norm between observation and simulated current source i
        F1.append(F1i)
    return np.array(F1)


def misfitF1_2_initialX0(A, b, **kwargs):
    """Transform the misfit F1 (punctual source inversion) to an initial solution M0
    using a 1/x^2 transformation"""

    norm_F1 = _normF1(A, b)
    x0F1 = 1 / norm_F1  # Inverse misfit using a 1/x^2 transformation
    x0F1_sum = x0F1 / sum(x0F1)  # normalize such as sum equal to 1
    M0 = x0F1_sum
    return norm_F1, M0


def product_moment(A, b):
    """Compute the product moment correlation after Binley et al. 1999
    .. math::

        r_{k}= \frac{\sum_{i}(D_{I}-\overline{D})(F_{i}(I_{k})-\overline{F}(I_{k}))}{\sqrt{\sum_{i}(D_{I}-\overline{D})^{2}}\sum_{i}(F_{i}(I_{k})-\overline{F}(I_{k}))^{2}}
    where $D_{i}$ is the $i^{th}$ measured transfer resistance and $F_{i}(I_{k})$ is the $i^{th}$  transfer resistance computed to unit current at location k.
    """
    # Estimate a covariance matrix, given data observation and weights and tranfert resitances measured.
    rpearson = []
    for i in range(np.shape(A)[1]):
        corr, _ = pearsonr(b, A[:, i])
        rpearson.append(corr)
        if i == 1:
            print(A[:, i])
    M0 = rpearson / np.sum(rpearson)  # normalize such as sum equal to 1
    TranslateMIN = np.min(M0)
    M0 = M0 - TranslateMIN  # translate, then transform */
    return M0
