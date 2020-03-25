Current source density inversion for imaging of plant roots
===

## icsd_class

actual code for the iCSD

## synthetic_complete is an example where:

in main_pybert_calc: 
* a 2d inv - 3d fwd routine is used for the forward calculation of VRTe responses
* the same mesh is used to calculated the response of a synthetic source
* the obtained vectors are saved for the iCSD

in icsd:
* a pareto inversion is used to invert the position of the synthetic source, calling icsd_class

## 3dinv_2dplot

example of 3D ert inversion interpolated to plot in 2D with matplotlib.
The gmsh rhizotron mesh is used, for the 2dinv-3dfwd see other examples.

## examples

some examples of laboratory acquisitions and iCSD (based on cotton experiments)

## malm_sequences

python code to generate the 204 malm sequence that is used for the rhizotron work
