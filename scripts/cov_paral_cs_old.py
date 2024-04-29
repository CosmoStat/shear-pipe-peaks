import random
import os
import numpy as np
import matplotlib.pyplot as plt

from pycs.astro.wl.mass_mapping import *
from pycs.sparsity.sparse2d.starlet import *
from pycs.misc.cosmostat_init import *
from pycs.astro.wl.hos_peaks_l1 import *

from sp_peaks import slics
from sp_peaks import mapping
from sp_peaks import summary_statistics
from sp_peaks import plotting
from multiprocessing import Pool

def compute_statistics(args):
    tile_file, method = args
    
    catalog_data = slics.read_catalogue_pd(tile_file)
    
    ra = catalog_data['RA']
    dec = catalog_data['Dec']
    g1_sim = catalog_data['gamma1_sim']
    g2_sim = catalog_data['gamma2_sim']

    x, y = radec2xy(np.mean(ra), np.mean(dec), ra, dec)
    Nx, Ny = int(SIZE_X_DEG / PIX_ARCMIN * 60), int(SIZE_Y_DEG / PIX_ARCMIN * 60)
    galmap = bin2d(x, y, npix=(Nx,Ny))
    mask = (galmap > 0).astype(int)
    
    sigma_noise = np.zeros_like(galmap)
    sigma_noise[mask != 0] = SHAPE_NOISE / np.sqrt(2 * galmap[mask != 0])
    sigma_noise[mask == 0] = np.max(sigma_noise[mask != 0]) # set the noise to the maximum value in the map where there are galaxies        noise_map_CFIS_z05 = sigma_noise * np.random.randn(sigma_noise.shape[0], sigma_noise.shape[1]) # generate noise map
    # noise_map_CFIS_z05 = sigma_noise * np.random.randn(sigma_noise.shape[0], sigma_noise.shape[1]) # generate noise map
    
    # # constant noise level
    # sigma_noise = np.ones_like(galmap) * SHAPE_NOISE / np.sqrt(2 * N_GAL)
    
    e1map, e2map = bin2d(x, y, npix=(Nx, Ny), v=(g1_sim, g2_sim)) 
    noise_e1 = np.random.randn(*e1map.shape) * sigma_noise
    noise_e2 = np.random.randn(*e2map.shape) * sigma_noise
    
    e1map = e1map + noise_e1 * mask # add noise to the map
    e2map = e2map + noise_e2 * mask # add noise to the map

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

    if method == 'ks':
        # Mass mapping class initialization
        M = massmap2d(name='mass')
        M.init_massmap(d.nx,d.ny)
        M.Verbose = False 
        ks = M.gamma_to_cf_kappa(e1map,-e2map) 
        ks = ks.real
        mass_map = ks
        
    if method == 'ksi':
        M = massmap2d(name='mass')
        M.init_massmap(d.nx, d.ny)
        M.DEF_niter = 50
        M.niter_debias = 30
        M.Verbose = False
        ksi =  M.iks(d.g1, d.g2, mask) 
        ksi = ksi.real    
        mass_map = ksi

    WT = starlet2d(gen2=False,l2norm=False, verb=False)
    WT.init_starlet(nx, ny, nscale=NSCALES)

    H = HOS_starlet_l1norm_peaks(WT)
    H.set_bins(Min=MIN_SNR, Max=MAX_SNR, nbins=NBINS)
    H.set_data(mass_map, SigmaMap=sigma_noise, Mask=mask)
    H.get_mono_scale_peaks(mass_map, sigma_noise, smoothing_sigma=6, mask=mask)
    peak_counts_single = H.Mono_Peaks_Count
    H.get_wtpeaks(Mask=mask)
    peak_counts_multi = H.Peaks_Count
    H.get_wtl1(NBINS_L1*2, Mask=mask, min_snr=-6, max_snr=6)
    l1norm_histogram = H.l1norm

    return peak_counts_single, peak_counts_multi, l1norm_histogram

#change directory to the directory where the script is located
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Define the path to the "master_file_cov.txt"
master_file_path = ".././input/master_file_cov.txt"

# Read the file paths from the "master_file_cov.txt"
with open(master_file_path, "r") as file:
    file_paths = file.readlines()
    file_paths = [path.strip() for path in file_paths]

# Parse these file paths
parsed_cov_data = slics.parse_cov_SLICS_filenames(file_paths)

los_numbers = np.unique(parsed_cov_data['LOS']) # List of all LOS numbers
num_realizations = 124 # Number of realizations
num_tiles_per_realization = 19 # Number of tiles to select for each realization

num_bins = 4
bin_number = 1

# Reconstruct 124 realisations of the survey by picking each tile from a random LOS, ensuring that each file is only included once.
collections_of_files = slics.survey_realizations_reconstruction(num_realizations, num_tiles_per_realization, bin_number, parsed_cov_data['LOS'], file_paths)

# Constants and Parameters
N_GAL = 7 
SIZE_X_DEG = 10.
SIZE_Y_DEG = 10.
PIX_ARCMIN = 1.
SHAPE_NOISE = 0.1
NSCALES = 5
# Histogram parameters
MIN_SNR = -2
MAX_SNR = 6
NBINS = 31
NBINS_L1 = 40

NUM_REALIZATIONS = 124  # Number of realizations


SS_PC_data_vectors = []
MS_PC_data_vectors = []
l1_norm_data_vectors = []

# Loop over realizations
for realization_files in collections_of_files:
    # Create a pool of worker processes
    with Pool(processes=num_tiles_per_realization) as pool:
        # Compute statistics in parallel for each tile in this realization
        args = [(file, 'ksi') for file in realization_files]
        results = pool.map(compute_statistics, args)
        
    # Now, process the results to separate them into individual lists
    peak_counts_single_realization = [result[0] for result in results]
    peak_counts_multi_realization = [result[1] for result in results]
    l1_norm_histogram_realization = [result[2] for result in results]

    # Compute the average vectors for this realization
    average_peak_counts_single = np.mean(peak_counts_single_realization, axis=0)
    average_peak_counts_multi = np.mean(peak_counts_multi_realization, axis=0)
    average_l1_norm_histogram = np.mean(l1_norm_histogram_realization, axis=0)

    # Append the average vectors for this realization to the lists of data vectors
    SS_PC_data_vectors.append(average_peak_counts_single)
    MS_PC_data_vectors.append(average_peak_counts_multi)
    l1_norm_data_vectors.append(average_l1_norm_histogram)


# Convert the list of data vectors into a NumPy array
SS_PC_data_vectors = np.array(SS_PC_data_vectors)
MS_PC_data_vectors = np.array(MS_PC_data_vectors)
l1_norm_data_vectors = np.array(l1_norm_data_vectors)
    
# save the data vectors as a .npy file
np.save('.././output/covariances/datavectors/data_vector_SS_PC_ksi_pc_smoothing_6_noise01.npy', SS_PC_data_vectors)
np.save('.././output/covariances/datavectors/data_vector_MS_PC_cs_ksi_pc_smoothing_6_noise01.npy', MS_PC_data_vectors)
np.save('.././output/covariances/datavectors/L1_norm_data_vector_cs_ksi_pc_smoothing_6_noise01.npy', l1_norm_data_vectors)

# ######################
# actual_NUM_REALIZATIONS = NUM_REALIZATIONS
# ######################

# # Reshape the data vectors 
# MS_PC_data_vectors_reshaped = MS_PC_data_vectors.reshape(actual_NUM_REALIZATIONS, -1)
# l1_norm_data_vectors_reshaped = l1_norm_data_vectors.reshape(actual_NUM_REALIZATIONS, -1)

# # Compute the average histogram vector across all realizations
# mean_SS_PC_over_realizations = np.mean(SS_PC_data_vectors, axis=0)
# mean_MS_PC_over_realizations = np.mean(MS_PC_data_vectors_reshaped, axis=0)
# mean_l1_norm_over_realizations = np.mean(l1_norm_data_vectors_reshaped, axis=0)

# # Compute the deviations of histograms in each realization from the average vector
# deviations_SS_PC = SS_PC_data_vectors - mean_SS_PC_over_realizations
# deviations_MS_PC = MS_PC_data_vectors_reshaped - mean_MS_PC_over_realizations
# deviations_l1_norm = l1_norm_data_vectors_reshaped - mean_l1_norm_over_realizations

# # Compute the covariance matrix 
# num_realizations_SS_PC, num_bins_SS_PC = SS_PC_data_vectors.shape
# num_realizations_MS_PC, num_bins_MS_PC = MS_PC_data_vectors_reshaped.shape
# num_realizations_l1_norm, num_bins_l1_norm = l1_norm_data_vectors_reshaped.shape

# covariance_matrix_SS_PC = np.dot(deviations_SS_PC.T, deviations_SS_PC) / (num_realizations_SS_PC - 1)
# covariance_matrix_MS_PC = np.dot(deviations_MS_PC.T, deviations_MS_PC) / (num_realizations_MS_PC - 1)
# covariance_matrix_l1_norm = np.dot(deviations_l1_norm.T, deviations_l1_norm) / (num_realizations_l1_norm - 1)

# # save the covariance matrix as a .npy file
# np.save('.././output/covariances/cov/covariance_matrix_SS_PC__ks_high_noise_less_bins.npy', covariance_matrix_SS_PC)
# np.save('.././output/covariances/cov/covariance_matrix_MS_PC_ks_high_noise_less_bins.npy', covariance_matrix_MS_PC)
# np.save('.././output/covariances/cov/covariance_matrix_l1_norm_ks_high_noise_less_bins.npy', covariance_matrix_l1_norm)

# # Calculate the diagonal of the covariance matrix
# diagonal_SS_PC = np.sqrt(np.diag(covariance_matrix_SS_PC))
# diagonal_MS_PC = np.sqrt(np.diag(covariance_matrix_MS_PC))
# diagonal_l1_norm = np.sqrt(np.diag(covariance_matrix_l1_norm))

# # Check for zero values and replace them with a small positive value
# diagonal_SS_PC[diagonal_SS_PC == 0] = 1e-10
# diagonal_MS_PC[diagonal_MS_PC == 0] = 1e-10
# diagonal_l1_norm[diagonal_l1_norm == 0] = 1e-10

# # Calculate the correlation coefficients
# correlation_matrix_SS_PC = covariance_matrix_SS_PC / np.outer(diagonal_SS_PC, diagonal_SS_PC)
# correlation_matrix_MS_PC = covariance_matrix_MS_PC / np.outer(diagonal_MS_PC, diagonal_MS_PC)
# correlation_matrix_l1_norm = covariance_matrix_l1_norm / np.outer(diagonal_l1_norm, diagonal_l1_norm)


