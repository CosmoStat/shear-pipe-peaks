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

def create_kappa_map(ra, dec, g1_sim, g2_sim, size_x_deg=10, size_y_deg=10, pixel_size_emap_amin=0.4):
    x, y = radec2xy(np.mean(ra), np.mean(dec), ra, dec) # Project (ra,dec) -> (x,y)

    Nx = int(size_x_deg / pixel_size_emap_amin * 60)
    Ny = int(size_y_deg / pixel_size_emap_amin * 60)

    e1map, e2map = bin2d(x, y, npix=(Nx, Ny), v=(g1_sim, g2_sim)) # bin the shear field into a 2D map
    emap = np.array([e1map,e2map]) # stack the two components into a single array

    kappaE, kappaB = ks93(e1map, -e2map) # make kappa map (the minus sign has to be here for our data conventions)
    return kappaE, kappaB


def add_noise_to_kappa_map(kappa_map, shape_noise, n_gal, pix_arcmin):
    sigma_noise_CFIS = shape_noise / (np.sqrt(2 * n_gal * pix_arcmin**2))
    noise_map_CFIS_z05 = sigma_noise_CFIS * np.random.randn(kappa_map.shape[0], kappa_map.shape[1]) # generate noise map
    kappa_map_noisy = kappa_map + noise_map_CFIS_z05 # Add noise to the mass map
    return kappa_map_noisy, noise_map_CFIS_z05

def smooth_kappa_map(kappa_map, pixel_size_emap_amin):
    # Set the standard deviation of the Gaussian filter based on the pixel size of the kappa map
    precision_Peaks = 2 / pixel_size_emap_amin # pixel_size_emap_amin is the pixel size of the kappa map in arcminutes
    kappa_map_smoothed = ndi.gaussian_filter(kappa_map, precision_Peaks)
    return kappa_map_smoothed

def convert_to_snr_map(kappa_map_smoothed, noise_map_smoothed):
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
