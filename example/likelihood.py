#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.ndimage.filters import gaussian_filter as gf
import scipy.stats as stats
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C

def gp_train(index_bin, params_train, obs_train):
    
    """
    
    Training of Gaussian processes per bin.
    It takes as input the bin indexes, training parameters and observables.
    It gives as output a list of GP and scaling values to avoid over-regularization.
    
    Parameters
    ----------
    
    index_bin: int
               bin index of the multipole for PS or S/N for Peaks
    params_train: numpy.ndarray of shape (#cosmo, #parameters)
                  training parameters from simulations
    obs_train: numpy.ndarray of shape (#cosmo, #realisations, #bins)
               training observables from simulations (e.g. the summary statistics)
   
    Returns
    -------
    gp: object
     GP for a given bin
    scaling: float
         scaling for the data values to avoid over-regularization
    """
    
    
    #means over the realisations for each bin looping over each cosmology
    obs_means=np.array([np.mean(obs_temp,axis=0)[index_bin] for obs_temp in obs_train])
    
    #sem over the realisations for each bin looping over each cosmology
    alpha_diag=np.array([stats.sem(obs_temp,axis=0)[index_bin] for obs_temp in obs_train])
    
    #data scaling factor
    scaling = (np.mean(obs_means))
    
    #scale the data values to avoid over-regularization
    obs_means /= scaling
    
    #scale the standard errors of the mean values to avoid over-regularization
    alpha_diag /= scaling
    
    #define the kernel
    kernel = C(5.0, (1e-4, 1e4)) * RBF([3, 0.3, 5], (1e-4, 1e4))
    
    #instantiate the gaussian process
    gp = GaussianProcessRegressor(kernel=kernel,alpha=(alpha_diag)**2,n_restarts_optimizer=50,normalize_y=True)
    gp.fit(params_train,obs_means)
    
    return gp, scaling
    


def GP_pred(params_pred, gp_list, scaling):
    
    """
    
    Computes the prediction for a new point in parameter space.
    It is called by the likelihoood to get thee theoretical predictions
    for the observables.
    
    Parameters
    ----------
    
    params_pred: list of floats
                  new point in parameter space to get the prediction for
               
    gp_list: list of objects
            list of GP for a given bin
    scaling: numpy.ndarray of shape (#cosmo, #realisations, #bins)
               training observables from simulations (e.g. PS or Peaks)
   
    Returns
    -------
    gp: list of objects
       list of GP per bin
    scaling: numpy.ndarray
         scaling for the data values to avoid over-regularization
    """
     
    pred_list=[]
    sigma_list=[]
    
    for gp in gp_list:
        
        pred,sigma=gp.predict(np.array(params_pred).reshape(1,len(params_pred)),return_std=True)
        
        pred_list.append(pred[0])
        
        sigma_list.append(sigma[0])
    
    return np.array(pred_list * scaling),\
                np.array(sigma_list * scaling)

def lnlike(theta, data, icov, gp_list, scaling, norm):
    
    """
    
    Likelihood function (loglike of a Gaussian)
    
    Parameters
    ----------
    
    theta: list
           point in parameter space
    data: numpy.ndarray 
          observed data
    icov: numpy.ndarray of shape (#bins, #bins)
          inverse of the covariance matrix
    gp_list: list of objects
             list of GP per bin
    scaling: numpy.ndarray
             scaling for the data values to avoid over-regularization
    norm: float
          the Hartlap factor
   
    Returns
    -------
    lnlikelihood: float
               value for the likelihood function
    """
 
    M_nu, Omega_m, A_s=theta
    
    obs_pred=GP_pred(theta,gp_list,scaling)[0] 

    lnlikelihood=-0.5*norm*np.dot(np.dot((data.reshape(len(data),1)-obs_pred.reshape(len(data),1)).T,icov),(data.reshape(len(data),1)-obs_pred.reshape(len(data),1)))
    
    
    if M_nu<=0:
        return -np.inf
    
    return lnlikelihood


def lnprior(theta, M_nu_min, M_nu_max, Omega_m_min, Omega_m_max, A_s_min, A_s_max):
    """
    
    Prior (log, flat). If the parameter values are within the minimum and maximum values
    returns 0.0. Otherwise returns -inf.
    
    Parameters
    ----------
    
    theta: list
           point in parameter space
    M_nu_min: float 
              minimum value for the neutrino mass
    M_nu_max: float 
              maximum value for the neutrino mass
    Omega_m_min: float 
              minimum value for the matter density parameter
    Omega_m_max: float 
              maximum value for the matter density parameter
    A_s_min: float 
              minimum value for the amplitude of the primordial power spectrum
    A_s_max: float 
              maximum value for the amplitude of the primordial power spectrum
   
    Returns
    -------
    0.0 if the parameter values are within the minimum and maximum values
    otherwise -np.inf
    """
    
    M_nu, Omega_m, A_s = theta
    
    if (M_nu_min < M_nu < M_nu_max and
            Omega_m_min < Omega_m < Omega_m_max and
            A_s_min < A_s < A_s_max):
        
        return 0.0
    
    return -np.inf

def lnpost(theta, data,icov, gp_list, scaling, norm, M_nu_min, M_nu_max, Omega_m_min, Omega_m_max, A_s_min, A_s_max):
    
    """
    
    Posterior (log). It is computed as lnprior + lnlike.
    
    Parameters
    ----------
    
    theta: list
           point in parameter space
    data: numpy.ndarray 
          observed data
    icov: numpy.ndarray of shape (#bins, #bins)
          inverse of the covariance matrix
    gp_list: list of objects
             list of GP per bin
    scaling: numpy.ndarray
             scaling for the data values to avoid over-regularization
    norm: float
          the Hartlap factor
    M_nu_min: float 
              minimum value for the neutrino mass
    M_nu_max: float 
              maximum value for the neutrino mass
    Omega_m_min: float 
              minimum value for the matter density parameter
    Omega_m_max: float 
              maximum value for the matter density parameter
    A_s_min: float 
              minimum value for the amplitude of the primordial power spectrum
    A_s_max: float 
              maximum value for the amplitude of the primordial power spectrum
   
    Returns
    -------
    lp+lnlike: float
               value for the posterior distribution
    """
        
    lp = lnprior(theta, M_nu_min, M_nu_max, Omega_m_min, Omega_m_max, A_s_min, A_s_max)
    if not np.isfinite(lp):
        return -np.inf
    
    return lp+lnlike(theta, data, icov, gp_list, scaling, norm)


