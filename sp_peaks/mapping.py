"""MAPPING.A

:Name: mapping.py
 
:Description: This package contains methods to create kappa maps and snr maps from shear catalogs                      
                                                                               
:Authors: Lucie Baumont <lucie.baumont@cea.fr> Andreas Tersenov <atersenov@physics.uoc.gr> Martin Kilbinger <martin.kilbinger@cea.fr> 
"""
import numpy as np
from lenspack.geometry.projections.gnom import radec2xy
from lenspack.utils import bin2d
from lenspack.image.inversion import ks93
from lenspack.image.transforms import starlet2d
from lenspack.starlet_l1norm import noise_coeff, get_l1norm_noisy
from astropy.stats import mad_std
from scipy import ndimage as ndi

def bin_shear_field(ra, dec, g1_sim, g2_sim, size_x_deg=10, size_y_deg=10, pixel_size_emap_amin=0.4):
    """
    Bin a shear field into a 2D map.

    Parameters
    ----------
    ra : array-like
        Right ascension coordinates in degrees.
    dec : array-like
        Declination coordinates in degrees.
    g1_sim : array-like
        Simulated shear component g1.
    g2_sim : array-like
        Simulated shear component g2.
    size_x_deg : float, optional
        Size of the map in degrees along the x-axis (default is 10).
    size_y_deg : float, optional
        Size of the map in degrees along the y-axis (default is 10).
    pixel_size_emap_amin : float, optional
        Pixel size of the output 2D map in arcminutes (default is 0.4).

    Returns
    -------
    numpy.ndarray, numpy.ndarray
        Two 2D maps representing the binned shear field components (e1map, e2map).

    """
    x, y = radec2xy(np.mean(ra), np.mean(dec), ra, dec) # Project (ra,dec) -> (x,y)

    Nx = int(size_x_deg / pixel_size_emap_amin * 60)
    Ny = int(size_y_deg / pixel_size_emap_amin * 60)

    e1map, e2map = bin2d(x, y, npix=(Nx, Ny), v=(g1_sim, g2_sim)) # bin the shear field into a 2D map
    return e1map, e2map

def add_noise_to_kappa_map(kappa_map, shape_noise, n_gal, pix_arcmin):
    """
    Add shape noise to a kappa map.

    Parameters
    ----------
    kappa_map : numpy.ndarray
        The input kappa map.
    shape_noise : float
        Shape noise level.
    n_gal : int
        Number of galaxies.
    pix_arcmin : float
        Pixel size in arcminutes.

    Returns
    -------
    numpy.ndarray, numpy.ndarray
        Two 2D maps representing the noisy kappa map and the shape noise map.

    """    
    sigma_noise_CFIS = shape_noise / (np.sqrt(2 * n_gal * pix_arcmin**2))
    noise_map_CFIS_z05 = sigma_noise_CFIS * np.random.randn(kappa_map.shape[0], kappa_map.shape[1]) # generate noise map
    kappa_map_noisy = kappa_map + noise_map_CFIS_z05 # Add noise to the mass map
    return kappa_map_noisy, noise_map_CFIS_z05

def smooth_kappa_map(kappa_map, pixel_size_emap_amin):
    """
    Smooth a kappa map using Gaussian filtering.

    Parameters
    ----------
    kappa_map : numpy.ndarray
        The input kappa map.
    pixel_size_emap_amin : float
        Pixel size of the kappa map in arcminutes.

    Returns
    -------
    numpy.ndarray
        The smoothed kappa map.

    """    
    # Set the standard deviation of the Gaussian filter based on the pixel size of the kappa map
    precision_Peaks = 2 / pixel_size_emap_amin # pixel_size_emap_amin is the pixel size of the kappa map in arcminutes
    kappa_map_smoothed = ndi.gaussian_filter(kappa_map, precision_Peaks)
    return kappa_map_smoothed

def convert_to_snr_map(kappa_map_smoothed, noise_map_smoothed):
    """
    Convert a smoothed kappa map and a noise map to a signal-to-noise ratio (SNR) map.

    Parameters
    ----------
    kappa_map_smoothed : numpy.ndarray
        The smoothed kappa map.
    noise_map_smoothed : numpy.ndarray
        The smoothed noise map.

    Returns
    -------
    numpy.ndarray
        The SNR map.

    """    
    snr_map = kappa_map_smoothed / np.std(noise_map_smoothed)
    return snr_map

def compute_multiscale_snr_maps(image, noise, nscales):
    """
    Compute SNR maps for each wavelet scale of a noisy image.
    
    Parameters:
        image (numpy.ndarray): The noiseless image.
        noise (numpy.ndarray): The noise to be added to the image.
        nscales (int): Number of wavelet scales for starlet decomposition.
        
    Returns:
        snr_maps (list of numpy.ndarray): List of SNR maps for each scale.
    """
    # Add noise to the noiseless image
    image_noisy = image + noise
    
    # Perform starlet decomposition
    image_starlet = starlet2d(image_noisy, nscales)
    
    # Estimate the noise level
    noise_estimate = mad_std(image_noisy)
    coeff_j = noise_coeff(image, nscales)
    
    snr_maps = []
    for image_j, std_co in zip(image_starlet, coeff_j):
        sigma_j = std_co * noise_estimate
        
        # Compute SNR map
        snr_map = image_j / sigma_j
        snr_maps.append(snr_map)
    
    return snr_maps
