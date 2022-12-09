"""
Single source Estimation of CSD
--------------------------------
The procedure to generate the data is describe in Peruzzo et al. (2020). 
The synthetic test evaluates the iCSD capability to locate a current point source (white dot). 
The white circles represent the VRTe and the white square the return electrode. 
The source of current was centered between four VRTe to highlight the distribution of the current density.

Here we compared the estimation using the F1 method and the Pearson method.

.. math::

    r_{k}= \frac{\sum_{i}(D_{I}-\overline{D})(F_{i}(I_{k})-\overline{F}(I_{k}))}{\sqrt{\sum_{i}(D_{I}-\overline{D})^{2}}\sum_{i}(F_{i}(I_{k})-\overline{F}(I_{k}))^{2}}
where $D_{i}$ is the $i^{th}$ measured transfer resistance and $F_{i}(I_{k})$ is the $i^{th}$  transfer resistance computed to unit current at location k.

*Peruzzo, L., Chou, C., Wu, Y. et al. Imaging of plant current pathways for non-invasive root Phenotyping using a newly developed electrical current source density approach. Plant Soil 450, 567–584 (2020). https://doi.org/10.1007/s11104-020-04529-w*
"""


# %%
# **Import the icsd package**
import matplotlib.pyplot as plt
from icsd.icsd3d import iCSD3d as i3d

# %%
# **Change inversion parameters**
# The default dimension type is 2d
path2files = "./iCSD_fig7/"
icsd = i3d(dirName=path2files)

icsd.createSurvey(fname_obs="ObsData.txt", fname_sim="VRTeSim.txt")
m0 = icsd.estimateM0(methodM0="F1", show=True)

# %%
# **Change inversion parameters**
# Change the estimate to Pearson

m0 = icsd.estimateM0(methodM0="Pearson", show=False)

fig, ax = plt.subplots(1)
icsd.showEstimateM0(ax=ax)
