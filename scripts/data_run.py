import os
import numpy as np
import matplotlib.pyplot as plt
from pycs.astro.wl.mass_mapping import *
from pycs.sparsity.sparse2d.starlet import *
from pycs.misc.cosmostat_init import *
from pycs.astro.wl.hos_peaks_l1 import *
import sp_peaks
from sp_peaks import slics
from sp_peaks import mapping
from sp_peaks import summary_statistics
from sp_peaks import plotting

# BOOKEEPING

# Read the file paths from master_file.txt
filename = ".././input/master_file.txt"
with open(filename, 'r') as file:
    file_paths = file.readlines()
    file_paths = [path.strip() for path in file_paths]
filename_data = slics.parse_SLICS_filenames(file_paths)

# Path to the .dat file
dat_file_path = "/home/tersenov/shear-pipe-peaks/example/CosmoTable.dat" # for CANDIDE
# dat_file_path = "/Users/atersenov/Software/shear-pipe-peaks/example/CosmoTable.dat" # for local machine

# Read the cosmological parameters from the .dat file
cosmo_params = slics.read_SLICS_cosmo_params(dat_file_path)
mapped_params = slics.map_cosmo_params_to_data(filename_data, dat_file_path)


# CONSTANTS AND PARAMETERS
N_GAL = 7 
SIZE_X_DEG = 10.
SIZE_Y_DEG = 10.
PIX_ARCMIN = 1.
SHAPE_NOISE = 0.44
NSCALES = 5
MIN_SNR = -2
MAX_SNR = 6
NBINS=31
NBINS_L1 = 40

# LOAD & READ CATALOGUE DATA
# CATALOG_FILE = "/n17data/tersenov/SLICS/Cosmo_DES/16_a/LOS4/DES_MocksCat_16_a_4_Bin3_LOS4_R4.dat" # for CANDIDE
CATALOG_FILE = "/Users/atersenov/Software/test_cosmostat/data/DES_MocksCat_16_a_4_Bin3_LOS4_R4.dat" # for local machine
catalog_data = slics.read_catalogue_pd(CATALOG_FILE)
ra = catalog_data['RA']
dec = catalog_data['Dec']
g1_sim = catalog_data['gamma1_sim']
g2_sim = catalog_data['gamma2_sim']
kappa_sim = catalog_data['kappa_sim']

# DATA PROCESSING
x, y = radec2xy(np.mean(ra), np.mean(dec), ra, dec)
Nx, Ny = int(SIZE_X_DEG / PIX_ARCMIN * 60), int(SIZE_Y_DEG / PIX_ARCMIN * 60)
galmap = bin2d(x, y, npix=(Nx,Ny))
mask = (galmap > 0).astype(int)

sigma_noise = np.zeros_like(galmap)
sigma_noise[mask != 0] = SHAPE_NOISE / np.sqrt(2 * galmap[mask != 0])
sigma_noise[mask == 0] = np.max(sigma_noise[mask != 0]) # set the noise to the maximum value in the map where there are galaxies
noise_map_CFIS_z05 = sigma_noise * np.random.randn(sigma_noise.shape[0], sigma_noise.shape[1]) # generate noise map

e1map, e2map = bin2d(x, y, npix=(Nx, Ny), v=(g1_sim, g2_sim)) 
# Shear data class initialization
d = shear_data()
d.g1 = e1map
d.g2 = -e2map
(nx,ny) = e1map.shape
d.mask = mask
# Shear noise covariance matrix
Ncov = np.zeros((nx,ny))
Ncov[mask > 0] = 2. * sigma_noise[mask > 0]**2
Ncov[mask == 0] = 1e9 # set the noise to the maximum value in the map where there are galaxies
d.Ncov = Ncov
d.nx = nx
d.ny = ny  
# Mass mapping class initialization
M = massmap2d(name='mass')
M.init_massmap(d.nx,d.ny)
M.DEF_niter = 50
Inpaint = True
M.niter_debias = 30
M.Verbose = False 
# KS
ks =  M.gamma_to_cf_kappa(e1map,-e2map) 
ks = ks.real
ks_noisy = ks + noise_map_CFIS_z05

WT = starlet2d(gen2=False,l2norm=False, verb=False)
WT.init_starlet(nx, ny, nscale=NSCALES)
H = HOS_starlet_l1norm_peaks(WT)
H.set_bins(Min=MIN_SNR, Max=MAX_SNR, nbins=NBINS)
H.set_data(ks_noisy, SigmaMap=sigma_noise, Mask=mask)
H.get_mono_scale_peaks(ks_noisy, sigma_noise, mask=mask)
H.get_wtpeaks(Mask=mask)
pc = H.Peaks_Count
H.get_wtl1(NBINS_L1*2, Mask=mask)