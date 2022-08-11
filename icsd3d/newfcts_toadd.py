# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 12:49:57 2020

@author: Benjamin
"""

# data misfit
def dataMisfit(p, obs, ini0):
    misfit = fmodel(p, ini0) - obs
    if forwardModel == 'Q':
        misfit = misfit # to help the solver with small Q
    if forwardModel == 'QP':
        misfit[len(misfit)//2:] *= 1e5 # TODO to be tested out
    return misfit 


# model misfit only for conductivities not depths
def modelMisfit(p):
    cond = p[:np.sum(vc)] # smoothing only for parameters elements (depth or cond)
    return cond[:-1] - cond[1:]

# set up regularisation
# p : parameter, app : ECa,
# pn : consecutive previous profile (for lateral smoothing)
# spn : profile from other survey (for time-lapse)
if regularization  == 'l1':
    def objfunc(p, app, pn, spn, alpha, beta, gamma, ini0):
        return np.sqrt(np.sum(np.abs(dataMisfit(p, app, ini0)))/len(app)
                       + alpha*np.sum(np.abs(modelMisfit(p)))/np.sum(vc)
                       + beta*np.sum(np.abs(p - pn))/len(p)
                       + gamma*np.sum(np.abs(p - spn))/len(p))
elif regularization == 'l2':
    def objfunc(p, app, pn, spn, alpha, beta, gamma, ini0):
        return np.sqrt(np.sum(dataMisfit(p, app, ini0)**2)/len(app)
                       + alpha*np.sum(modelMisfit(p)**2)/np.sum(vc)
                       + beta*np.sum((p - pn)**2)/len(p)
                       + gamma*np.sum((p - spn)**2)/len(p))

rmse[j] = np.sqrt(np.sum(dataMisfit(out, obs, ini0)**2)/np.sum(obs**2)/len(obs))*100

