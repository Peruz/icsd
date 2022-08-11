# -*- coding: utf-8 -*-
"""
Created on Mon May 11 16:22:08 2020
@author: Benjamin
Estimation of initial model based on the physical assumption that a single source current
can describe the pattern of the masse anomaly
"""
from numpy import linalg as LA
import numpy as np
from scipy.stats import pearsonr
import matplotlib.pyplot as plt


### Individual Misfit  
def _normF1(A,b,**kwargs):
    """compute the norm between observation data and individual green functions"""
    F1=[]
    for i in range(np.shape(A)[1]):
        F1i = LA.norm((b-A[:,i]))# norm between observation and simulated current source i
        F1.append(F1i)
        
        for key, value in kwargs.items():
            if key == 'int_plot':
                TEST = value
                
            # if i in np.arange(5,len(TEST[1][:,0]),50):
            #     fig, axs = plt.subplots(1, 2, figsize=(20,5))
            #     # axs[0].plot(A[:,i])
            #     # axs[0].set_title('A[:,i]')
            #     # axs[0].legend(['b:'+ str(b[0])+ '  A: ' + str(A[1,i])])
            #     # axs[1].plot(b)
            #     # axs[1].set_title('b' + str(i))
    
            #     sc=axs[0].scatter(TEST[0][:,1], TEST[0][:,2], c=b,
            #                   cmap ='coolwarm',s=200) 
            #     axs[0].set_ylabel('y [m]',fontsize=15)
            #     axs[0].set_xlabel('x [m]',fontsize=15)
            #     axs[0].set_title('b: ' + str(b[0]))
            #     cbar = plt.colorbar(sc,ax=axs[0])
                
            #     sc=axs[1].scatter(TEST[0][:,1], TEST[0][:,2], c=A[:,i],
            #                   cmap ='coolwarm',s=200)    #vmin=-200,vmax=200)            
            #     axs[1].scatter(TEST[1][i,0], TEST[1][i,1], c='g')
            #     axs[1].set_ylabel('y [m]',fontsize=15)
            #     axs[1].set_xlabel('x [m]',fontsize=15)
            #     axs[1].set_title('A[:,i], i=' + str(i) + '  A0:' + str(A[0,i]))
            #     cbar = plt.colorbar(sc,ax=axs[1])

         
    # normalise such as the sum is equal to 1 for current conservation
    # norm_F1 = (F1 - min(F1)) / (max(F1) - min(F1)) 
    # plt.figure()
    # plt.plot(F1)  
    

    
    # print(type(norm_F1))
    # print(type(F1))
        
    # return norm_F1
    return np.array(F1)

def misfitF1_2_initialX0(A,b,**kwargs):
    """Transform the misfit F1 (punctual source inversion) to an initial solution M0  
    using a 1/x^2 transformation"""

    # diff_Ab = []
    # plt.figure()  
    # for i in range(np.shape(A)[1]):
    # # for i in list(np.arange(0, 90, 20)):
    #     plt.plot(A[:,i],'*')
    #     # diff_Ab.append(abs(A[:,i])-abs(b))
    #     diff_Ab.append(A[:,i]-b)
    # plt.plot(b,'-.')
    
    
    # # diff_Ab = np.vstack(diff_Ab)
    # # np.argmin(np.sum(abs(diff_Ab),axis=1))
    # # np.argmin(np.mean(diff_Ab,axis=0))
    
    # diff_Ab[40]
    # # diff_Ab[0]

    # plt.figure()  
    # # plt.plot(A[:,np.argmin(np.sum(diff_Ab,axis=0))],'*')
    # plt.plot(A[:,40],'*')
    # plt.plot(b,'-.')
  
    norm_F1 = _normF1(A,b)

    # x0F1_2=1./((norm_F1+1)*(norm_F1+1)) # Inverse misfit using a 1/x^2 transformation
    # x0F1=1./((norm_F1+1)) # Inverse misfit using a 1/x^2 transformation
    x0F1=1/norm_F1 # Inverse misfit using a 1/x^2 transformation
    x0F1_sum= x0F1/sum(x0F1) # normalize such as sum equal to 1
    M0=x0F1_sum

    # plt.figure()
    # plt.plot(x0F1_sum)
    # plt.title('sum=1')

    return norm_F1, M0


def product_moment(A,b):
    """ Compute the product moment correlation after Binley et al. 1999
    .. math:: 

        r_{k}= \frac{\sum_{i}(D_{I}-\overline{D})(F_{i}(I_{k})-\overline{F}(I_{k}))}{\sqrt{\sum_{i}(D_{I}-\overline{D})^{2}}\sum_{i}(F_{i}(I_{k})-\overline{F}(I_{k}))^{2}}
    where $D_{i}$ is the $i^{th}$ measured transfer resistance and $F_{i}(I_{k})$ is the $i^{th}$  transfer resistance computed to unit current at location k. 
    """
    # self.icsd_init()            
    # Estimate a covariance matrix, given data observation and weights and tranfert resitances measured.
    rpearson=[]
    for i in range(np.shape(A)[1]):            
        corr, _ = pearsonr(b, A[:,i])
        rpearson.append(corr)
        if i ==1:
            print(A[:,i])
    # M0=rpearson
    # print(np.isnan(rpearson))
    M0= rpearson/np.sum(rpearson) # normalize such as sum equal to 1
    TranslateMIN= np.min(M0) 
    M0 = M0 - TranslateMIN # translate, then transform */
    # check if sum equal to 1 and bounds between 0 and 1
    # print(np.sum(rpearson))
    return M0
        
        
