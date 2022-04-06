#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.ndimage.filters import gaussian_filter as gf




def kappa_eff(kappa_z1, kappa_z2, lambda_1, lambda_2):
    """
    
    Computes the effective convergence map at :math:`z_{eff}` interpolating
    between two maps at known redshifts :math:`z_1` and :math:`z_2`.
    
    Parameters
    ----------
    kappa_z1: numpy.ndarray of size (npix, npix)
              convergence map at :math:`z_1`
    kappa_z2: numpy.ndarray of size (npix, npix)
              convergence map at :math:`z_2`
    lambda_1: float
              first weight interpolation parameter
    lambda_2: float
              second weight interpolation parameter
              
    Returns
    -------
    kappa_eff: numpy.ndarray of size (npix, npix)
            Effective convergence map at redshift :math:`z_{eff}`
    """
    
    
    kappa_eff = kappa_z1*(lambda_1*lambda_2) + kappa_z2*(1 - lambda_1*lambda_2)
    
    return kappa_eff
    

def power_spectrum(map_data, resolution, npix):
    """
    
    Computes the power of a 2D convergence map of size (npix, npix)
    
    Parameters
    ----------
    
    map_data: numpy.ndarray of size (npix, npix)
              convergence map
    resolution: float
                resolution of the map in arcmin
    npix: int
          number of pixels per side  
    
    Returns
    -------
    ell: numpy.ndarray 
       multipoles   
    power_spectrum_1d: numpy.ndarray
                    1D power spectrum of the map
    
    """
    
    pixel_size = np.pi * resolution / 180. / 60. #rad/pixel    
    data_ft = np.abs(np.fft.fft2(map_data))
    data_ft_shifted = np.fft.fftshift(data_ft)
    power_spectrum_2d = np.abs(data_ft_shifted * np.conjugate(data_ft_shifted)) / npix**2
    nyquist = np.int(data_ft_shifted.shape[0] / 2)
    center = power_spectrum_2d.shape[0]/2
    v, u = np.indices((power_spectrum_2d.shape))
    k = np.sqrt((u - center)**2 + (v - center)**2)
    k = k.astype('int32')

    tbin = np.bincount(k.ravel(), power_spectrum_2d.ravel())
    nr = np.bincount(k.ravel())
    radialprofile = tbin / nr
    power_spectrum_1d = radialprofile[:nyquist] * (pixel_size)**2
    k = np.arange(power_spectrum_1d.shape[0])
    ell = 2. * np.pi * k / pixel_size / 360

    return ell, power_spectrum_1d



def takes_params(fn):
    
    """
    
    Extracts cosmological parameters from the file name 
    
    Parameters
    ----------
    
    fn: str
        file name of observable array
   
    Returns
    -------
    mnv: float 
       neutrino mass   
    Om: float
       matter density parameter
    As: float
       amplite of the primordial power spectrum   
    
    """
    
    mnv=float(fn.split('convergence_gal_mnv')[1].split('_')[0])
    Om=float(fn.split('om')[1].split('_')[0])
    As=float(fn.split('As')[1].split('_')[0])

    return mnv,Om,As
    

