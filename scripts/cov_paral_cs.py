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
bin_number = 2

# Reconstruct 124 realisations of the survey by picking each tile from a random LOS, ensuring that each file is only included once.
collections_of_files = slics.survey_realizations_reconstruction(num_realizations, num_tiles_per_realization, bin_number, parsed_cov_data['LOS'], file_paths)

# Constants and Parameters
N_GAL = 7 
SIZE_X_DEG = 10.
SIZE_Y_DEG = 10.
PIX_ARCMIN = 1.
SHAPE_NOISE = 0.44
NSCALES = 5

# Histogram parameters
MIN_SNR = -2
MAX_SNR = 6
NBINS = 31
NBINS_L1 = 40

NUM_REALIZATIONS = 124  # Number of realizations

def compute_statistics(tile_file):
    catalog_data = slics.read_catalogue_pd(tile_file)
    
    ra = catalog_data['RA']
    dec = catalog_data['Dec']
    g1_sim = catalog_data['gamma1_sim']
    g2_sim = catalog_data['gamma2_sim']

    x, y = radec2xy(np.mean(ra), np.mean(dec), ra, dec)
    Nx, Ny = int(SIZE_X_DEG / PIX_ARCMIN * 60), int(SIZE_Y_DEG / PIX_ARCMIN * 60)
    galmap = bin2d(x, y, npix=(Nx,Ny))
    mask = (galmap > 0).astype(int)
    sigma_noise = np.where(mask != 0, SHAPE_NOISE / np.sqrt(2 * galmap), np.max(sigma_noise[mask != 0]))
    noise_map_CFIS_z05 = sigma_noise * np.random.randn(sigma_noise.shape[0], sigma_noise.shape[1]) # generate noise map
    e1map, e2map = bin2d(x, y, npix=(Nx, Ny), v=(g1_sim, g2_sim)) 

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
    H.get_mono_scale_peaks(ks_noisy, sigma_noise, mask=mask)
    peak_counts_single = H.Mono_Peaks_Count

    H.set_data(ks_noisy, SigmaMap=sigma_noise, Mask=mask)
    H.get_wtpeaks(Mask=mask)
    peak_counts_multi = H.Peaks_Count

    H.get_wtl1(NBINS_L1*2, Mask=mask)
    l1norm_histogram = H.l1norm

    return peak_counts_single, peak_counts_multi, l1norm_histogram


SS_PC_data_vectors = []
MS_PC_data_vectors = []
l1_norm_data_vectors = []

# Loop over realizations
for realization_files in collections_of_files:
    # Create a pool of worker processes
    with Pool(processes=num_tiles_per_realization) as pool:
        # Compute L1-norm histograms in parallel for each tile in this realization
        peak_counts_single_realization, peak_counts_multi_realization, l1_norm_histogram_realization = pool.map(compute_statistics, realization_files)

    # Compute the average vectors for this realization
    average_peak_counts_single = np.mean(peak_counts_single_realization, axis=0)
    average_peak_counts_multi = np.mean(peak_counts_multi_realization, axis=0)
    average_l1_norm_histogram = np.mean(l1_norm_histogram_realization, axis=0)

    # Append the average L1-norm histogram vector for this realization to the list of data vectors
    SS_PC_data_vectors.append(average_peak_counts_single)
    MS_PC_data_vectors.append(average_peak_counts_multi)
    l1_norm_data_vectors.append(average_l1_norm_histogram)



# Convert the list of data vectors into a NumPy array
SS_PC_data_vectors = np.array(SS_PC_data_vectors)
MS_PC_data_vectors = np.array(MS_PC_data_vectors)
l1_norm_data_vectors = np.array(l1_norm_data_vectors)
    
# save the data vectors as a .npy file
np.save('.././output/data_vector_SS_PC.npy', SS_PC_data_vectors)
np.save('.././output/data_vector_MS_PC.npy', MS_PC_data_vectors)
np.save('.././output/L1_norm_data_vector.npy', l1_norm_data_vectors)

# Reshape the data vectors 
MS_PC_data_vectors_reshaped = MS_PC_data_vectors.reshape(NUM_REALIZATIONS, -1)
l1_norm_data_vectors_reshaped = l1_norm_data_vectors.reshape(NUM_REALIZATIONS, -1)

# Compute the average histogram vector across all realizations
mean_SS_PC_over_realizations = np.mean(SS_PC_data_vectors, axis=0)
mean_MS_PC_over_realizations = np.mean(MS_PC_data_vectors_reshaped, axis=0)
mean_l1_norm_over_realizations = np.mean(l1_norm_data_vectors_reshaped, axis=0)

# Compute the deviations of histograms in each realization from the average vector
deviations_SS_PC = SS_PC_data_vectors - mean_SS_PC_over_realizations
deviations_MS_PC = MS_PC_data_vectors_reshaped - mean_MS_PC_over_realizations
deviations_l1_norm = l1_norm_data_vectors_reshaped - mean_l1_norm_over_realizations

# Compute the covariance matrix 
num_realizations_SS_PC, num_bins_SS_PC = SS_PC_data_vectors.shape[0]
num_realizations_MS_PC, num_bins_MS_PC = MS_PC_data_vectors_reshaped.shape
num_realizations_l1_norm, num_bins_l1_norm = l1_norm_data_vectors_reshaped.shape

covariance_matrix_SS_PC = np.dot(deviations_SS_PC.T, deviations_SS_PC) / (num_realizations_SS_PC - 1)
covariance_matrix_MS_PC = np.dot(deviations_MS_PC.T, deviations_MS_PC) / (num_realizations_MS_PC - 1)
covariance_matrix_l1_norm = np.dot(deviations_l1_norm.T, deviations_l1_norm) / (num_realizations_l1_norm - 1)

# save the covariance matrix as a .npy file
np.save('.././output/covariance_matrix_SS_PC.npy', covariance_matrix_SS_PC)
np.save('.././output/covariance_matrix_MS_PC.npy', covariance_matrix_MS_PC)
np.save('.././output/covariance_matrix_l1_norm.npy', covariance_matrix_l1_norm)

# Calculate the diagonal of the covariance matrix
diagonal_SS_PC = np.sqrt(np.diag(covariance_matrix_SS_PC))
diagonal_MS_PC = np.sqrt(np.diag(covariance_matrix_MS_PC))
diagonal_l1_norm = np.sqrt(np.diag(covariance_matrix_l1_norm))

# Check for zero values and replace them with a small positive value
diagonal_SS_PC[diagonal_SS_PC == 0] = 1e-10
diagonal_MS_PC[diagonal_MS_PC == 0] = 1e-10
diagonal_l1_norm[diagonal_l1_norm == 0] = 1e-10

# Calculate the correlation coefficients
correlation_matrix_SS_PC = covariance_matrix_SS_PC / np.outer(diagonal_SS_PC, diagonal_SS_PC)
correlation_matrix_MS_PC = covariance_matrix_MS_PC / np.outer(diagonal_MS_PC, diagonal_MS_PC)
correlation_matrix_l1_norm = covariance_matrix_l1_norm / np.outer(diagonal_l1_norm, diagonal_l1_norm)

# Plotting
plotting.plot_map(covariance_matrix_SS_PC, title='SS-PC', cmap='viridis', vmin=None, vmax=None)
plotting.plot_map(covariance_matrix_MS_PC, title='MS-PC', cmap='viridis', vmin=None, vmax=None)
plotting.plot_map(covariance_matrix_l1_norm, title='l1-norm', cmap='viridis', vmin=None, vmax=None)

# Try different plotting
correlation_matrix_SS_PC = np.corrcoef(covariance_matrix_SS_PC.T)
correlation_matrix_MS_PC = np.corrcoef(covariance_matrix_MS_PC.T)
correlation_matrix_l1_norm = np.corrcoef(covariance_matrix_l1_norm.T)

plotting.plot_map(correlation_matrix_SS_PC, title='SS-PC', cmap='viridis', vmin=None, vmax=None)
plotting.plot_map(correlation_matrix_MS_PC, title='MS-PC', cmap='viridis', vmin=None, vmax=None)
plotting.plot_map(correlation_matrix_l1_norm, title='l1-norm', cmap='viridis', vmin=None, vmax=None)

# save the figure as pdf
plt.savefig('.././output/covariance_matrix_SS_PC.pdf')
plt.savefig('.././output/covariance_matrix_MS_PC.pdf')
plt.savefig('.././output/covariance_matrix_l1_norm.pdf')
