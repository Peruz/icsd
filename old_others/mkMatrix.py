import re
import numpy as np
import sys
nd = 61 # number of data per VRTe

def readBERT(infile):
    coordinates = []
    data = []
    read_coord = 0
    read_data = 0
    fid = open (infile)
    for line in fid:
        if line[0] == '#':
            if line[2] == 'x':
                read_coord = 1 # start coordinates
            if line[2] == 'a':
                read_coord = 0 # stop coordinates
                read_data = 1 # start data
        else:
            if read_coord == 1:
                line_coord = re.findall('[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?',line)
                coordinates.append(line_coord)
            if read_data == 1:
                line_data = re.findall('[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?',line)
                data.append(line_data)
    del coordinates[-1]
    del data[-1]
    coord = np.array(coordinates, dtype = 'f' )
    num_elec = coord.shape[0] 
    data = np.array(data, dtype = 'f')
    num_data = data.shape[0] 
    return(coord,data,num_elec,num_data)


coord, data, num_elec, num_data = readBERT(sys.argv[1])
num_VRTe = 120
matrix = np.ones((nd,num_VRTe))

for i in range(num_VRTe):
    col = data[i * nd : (i + 1) * nd, 4]
    matrix[:, i] = col

np.savetxt('matrix.txt', matrix)
np.savetxt('coordinates.txt', coord)




