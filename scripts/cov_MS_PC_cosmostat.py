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
NBINS=31

NUM_REALIZATIONS = 124  # Number of realizations

# Initialize an empty list to store data vectors for each realization
data_vectors = []

# Initialize an empty list to store data the Multi-Scale Peak Counts vectors for each realization
MS_PC_data_vectors = []

# Loop over realizations
for realization_files in collections_of_files:  
    # Initialize an empty list to store peak counts vectors for each scale in this realization
    peak_counts_realization = []

    # Loop over files (tiles) in this realization
    for tile_file in realization_files:
        # Load the catalog data for this tile (similar to your previous code)
        # catalog_data = slics.read_catalogue(tile_file)
        catalog_data = slics.read_catalogue_pd(tile_file)

        # Extract data from the catalog (similar to your previous code)
        ra = catalog_data['RA']
        dec = catalog_data['Dec']
        g1_sim = catalog_data['gamma1_sim']
        g2_sim = catalog_data['gamma2_sim']

                # Calculate peak counts for this tile
        x, y = radec2xy(np.mean(ra), np.mean(dec), ra, dec)
        Nx, Ny = int(SIZE_X_DEG / PIX_ARCMIN * 60), int(SIZE_Y_DEG / PIX_ARCMIN * 60)
        galmap = bin2d(x, y, npix=(Nx,Ny))
        mask = (galmap > 0).astype(int)
        sigma_noise = np.where(mask != 0, SHAPE_NOISE / np.sqrt(2 * galmap), np.max(sigma_noise[mask != 0]))
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
        H.get_wtpeaks(Mask=mask)
        peak_counts_multi = H.Peaks_Count

        # Append peak counts for this tile to the list for this realization
        peak_counts_realization.append(peak_counts_multi)

    # Compute the average peak counts vector for this realization
    average_peak_counts = np.mean(peak_counts_realization, axis=0)

    # Append the average peak counts vector for this realization to the list of data vectors
    MS_PC_data_vectors.append(average_peak_counts)

# Convert the list of data vectors into a NumPy array
MS_PC_data_vectors = np.array(MS_PC_data_vectors)
# save the data vectors as a .npy file
np.save('.././output/data_vector_MS_PC.npy', MS_PC_data_vectors)




data_vectors_reshaped = MS_PC_data_vectors.reshape(NUM_REALIZATIONS, -1)  

# Compute the average peak counts vector across all realizations for multiscale peak counts
mean_PC_over_realizations_multi = np.mean(data_vectors_reshaped, axis=0)

# Compute the deviations of peak counts in each realization from the average vector
deviations_multi = data_vectors_reshaped - mean_PC_over_realizations_multi

# Compute the covariance matrix for multiscale peak counts
num_realizations_multi, num_bins_multi = data_vectors_reshaped.shape
covariance_matrix_multi = np.dot(deviations_multi.T, deviations_multi) / (num_realizations_multi - 1)
# save the covariance matrix as a .npy file
np.save('.././output/covariance_matrix_MS_PC.npy', covariance_matrix_multi)

# Normalize the covariance matrix for multiscale peak counts
diagonal_multi = np.sqrt(np.diag(covariance_matrix_multi))
# Add a small epsilon to the diagonal to avoid division by zero
epsilon = 1e-10
covariance_matrix_normalized_multi = covariance_matrix_multi / (np.outer(diagonal_multi, diagonal_multi)+epsilon)


plotting.plot_map(covariance_matrix_normalized_multi, title='Multi-scale Peak Counts', cmap='viridis', vmin=None, vmax=None)
# save the figure as pdf
plt.savefig('.././output/covariance_matrix_MS_PC.pdf')
