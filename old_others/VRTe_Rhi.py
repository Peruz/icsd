# this represents the model #
import numpy as np
import os
import sys

#remove old modeloutput
if os.path.isfile('PESTinput.txt'):
    os.remove('PESTinput.txt')
#load new poarameter values
weights = np.loadtxt('Weight.wgt', skiprows = 0) # skip #data
#print(weights)
num_weights = weights.shape[0]
#load all VRTe simulations
VRTeData = np.loadtxt('VRTe.sim').reshape(-1,1)
num_VRTeData = VRTeData.shape[0]
num_eachVRTeData = num_VRTeData/num_weights
check = num_eachVRTeData.is_integer()
if not check:
    sys.exit('check failled, error with number of weights and VRTe lenght')
else:
    num_eachVRTeData_Int = np.array(num_eachVRTeData, dtype = 'i')

out = np.zeros((num_eachVRTeData_Int, 1))
#print('first 3 data in initial globalVRTe are: ',out[0:3])
for i in range(0,num_weights):
    w = weights[i]
    first = i * num_eachVRTeData_Int
    last = (i + 1) * num_eachVRTeData_Int
    iVRTe = VRTeData[first:last]
    iVRTeWeighted = iVRTe * weights[i]
    #print('VRTe num: ',i+1,'\nfirst 3 data are: ',iVRTe[0:3],'\nweighted data are: ',iVRTeWeighted[0:3]) 
    out = out + iVRTeWeighted
    #print('added contribution VRTe',i+1,' with value',w,' to globalVRTe; first 3 data now are:\n ',out[0:3])

np.savetxt('PESTinput.txt', out, header='data',fmt='%16.8f') # save and add #data (pest marker)
