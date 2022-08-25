# -*- coding: utf-8 -*-
"""
Created on Mon May 11 15:18:31 2020

@author: Benjamin
"""

import numpy as np
import os 

#import pybert as pb
#from pybert import tdip
#import pygimli as pg

#%% LOAD from ICSD3d format: .txt file tab separated         

def load_coord(path, filename, dim):
    """ load coordinates of the virtual current sources
    
    Parameters
    ----------

    """
    coord = np.loadtxt(path + filename)
    if dim == '2d':
        coord_x, coord_y =  coord[:, 0],  coord[:, 1]
        return coord_x, coord_y, coord
    else:
        coord_x, coord_y, coord_z = coord[:, 0], coord[:, 1], coord[:, 2]
        return coord_x, coord_y, coord_z, coord

def load_obs(path, filename,index=None):
    """ load the observations file (normalised voltages)
    
    Parameters
    ----------

    """       
    if filename.endswith('.data'): # TDIP data importer
        bfile = pg.load(filename)
        # b = (bfile['m'+str(index)]).array()
        # # print(str(index)+str(b[0]))
        # import matplotlib.pyplot as plt 
        # fig, ax = plt.subplots()
        # ax.plot(b)
        b = (bfile['m'+str(index)]/bfile['k']).array()
        # fig, ax = plt.subplots()
        # ax.plot(b)

    else:
        b = np.loadtxt(path+ filename)
                
    return b

def load_sim(path, data_sim, knorm=True):
    """ load the simulated green functions file
    
    Parameters
    ----------

    """
    if isinstance(data_sim, str):
        if '.data' in data_sim:
            tdip_data_sim = pg.load(path + data_sim)
            A = []
            i=0
            if knorm == True:
                while True:
                    try: 
                      A.append(((tdip_data_sim['m'+str(i)])/tdip_data_sim['k']).array())
                      i += 1
                      # print(i)
                    except:
                     break
            else:
                while True:
                    try: 
                      A.append((tdip_data_sim['m'+str(i)]).array())
                      i += 1
                    except:
                     break
                 
            A = np.vstack(A) 
            A = np.transpose(A) 
            
        else:
            A = np.loadtxt(path+ data_sim)
    else:
        A = data_sim
    print(np.shape(A))
    print('*'*36)
    
    return A

def load_geom(path):
    """ load the geometry of the acquisition (*geom file custum for Mise-à-la-masse data)
    
    Parameters
    ----------

    """
    geom_files = [f for f in os.listdir(path) if f.endswith('.geom')]
    if len(geom_files) != 1:
        raise ValueError('should be only one geom file in the current directory')
    
    fileNameElec = geom_files[0]  
    line_number = 0
    line_of_injection = []
    line_of_remotes = []
    # Open the file in read only mode
    with open(path + fileNameElec, 'r') as read_obj:
        # Read all lines in the file one by one
        for line in read_obj:
            # For each line, check if line contains the string
            # print(line)
            line_number += 1
            if ('#Remote') in line:
                # If yes, then add the line number & line as a tuple in the list
                line_of_remotes.append((line_number))
            if ('#Injection') in line:
                line_of_injection.append((line_number))
    
    # print(line_of_injection)
    # print(line_of_remotes)
    RemLineNb= int(line_of_remotes[0])-1
    Injection= int(line_of_injection[0])-1
    
    coordE = np.loadtxt(path+ fileNameElec)
    pointsE= np.vstack(coordE[:RemLineNb,1:4])
    
    return RemLineNb, Injection, coordE, pointsE
            
        
#%% UTILS: LOAD from gimli or resipy format         

    
def dataImport(SimFile=None,ObsFile=None):
    """Data importer for common data files (Resipy and Gimli)
    Import and parse observation files, simulated file and geometry file
    
    Parameters
    ----------

    """
    
    if fileExt=='*.data':
        print('pygimli format import')
        A, b = load_pg_data(SimFile=None,ObsFile=None)
    if fileExt=='*.dat':
        print('resipy format import') 
        A, b = load_resipy_data(SimFile=None,ObsFile=None)

        
def loadTDIPSurvey(fname_obs,fname_sim, Vp_norm=True):
    """Data importer for common data files (Gimli)
    Import and parse observation files, simulated green fct file
    
    Parameters
    ----------
    - Vp_norm: True: normalise chargeability by multiplying with Vp
    Return
    ----------
    - Vs containing all gates + utlimate gate is Vp
    - Vs_green containing all gates + utlimate gate is Vp_green
    """     

    Vs = []
    Vs_green = []
    if isinstance(fname_obs, str):
        if '.dat' in fname_obs:
            i=0 
            tdip_obs = pg.load(fname_obs)
            Vs.append(((tdip_obs['r'])).array())
            # print(np.shape(Vs))
            if Vp_norm == True:
                # print('import obs')
                while True:
                    try: 
                      try:
                          Vs.append(((tdip_obs['M'+str(i)])*tdip_obs['r']).array())
                      except:
                          Vs.append(((tdip_obs['m'+str(i)])*tdip_obs['r']).array())
                      i += 1
                      # print(i)
                      # print(np.shape(Vs))
                    except:
                     break
            else:
                while True:
                    try: 
                      Vs.append((tdip_obs['m'+str(i)]).array())
                      i += 1
                    except:
                        break
            Vs = np.vstack(Vs) 
            Vs = np.transpose(Vs)

        if '.dat' in fname_sim:
            i=0
            tdip_sim = pg.load(fname_sim)
            Vs_green.append(((tdip_sim['r'])).array())
            if Vp_norm == True:
                while True:
                    try: 
                      Vs_green.append(((tdip_sim['m'+str(i)])*tdip_sim['r']).array())
                      i += 1
                    except:
                     break
            else:
                while True:
                    try: 
                      Vs_green.append((tdip_sim['m'+str(i)]).array())
                      i += 1
                    except:
                        break
            Vs_green = np.vstack(Vs_green) 
            Vs_green = np.transpose(Vs_green) 
            # print(np.shape(Vs_green))
        else:
            Vs_green = np.loadtxt(fname_sim)
            Vs_green = Vs_green*np.ones([i+1,np.shape(Vs_green)[0]])
            Vs_green = np.transpose(Vs_green) 
            # print(np.shape(Vs_green))
    else:   
        print('Can''t import numpy array data')


       
    return Vs, Vs_green

