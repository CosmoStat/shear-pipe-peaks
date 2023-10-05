import os
import random
import sp_peaks
from sp_peaks import slics
from sp_peaks import mapping
from sp_peaks import summary_statistics
from sp_peaks import plotting
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
from lenspack.geometry.projections.gnom import radec2xy
from lenspack.utils import bin2d
from lenspack.image.inversion import ks93
import lenspack.peaks as peaks
from lenspack.starlet_l1norm import noise_coeff, get_l1norm_noisy
from lenspack.image.transforms import starlet2d
from astropy.stats import mad_std

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
PIX_ARCMIN = 0.4
SHAPE_NOISE = 0.44
NSCALES = 5
NBINS = 40 
KAPPA_SNR = np.linspace(-2, 6, 31)
NUM_REALIZATIONS = 124  # Number of realizations

# Initialize an empty list to store data vectors for each realization
data_vectors = []

# Loop over realizations
for realization_files in collections_of_files: 
    # Initialize an empty list to store peak counts vectors for each tile in this realization
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

        # Calculate peak counts for this tile (similar to your previous code)
        e1map, e2map = mapping.bin_shear_field(ra, dec, g1_sim, g2_sim)
        kappaE, _ = ks93(e1map, -e2map)
        kappaE_noisy, noise_map_CFIS_z05 = mapping.add_noise_to_kappa_map(kappaE, SHAPE_NOISE, N_GAL, PIX_ARCMIN)
        kappaE_noisy_smoothed = mapping.smooth_kappa_map(kappaE_noisy, PIX_ARCMIN) 
        snr = mapping.convert_to_snr_map(kappaE_noisy_smoothed, kappaE_noisy_smoothed)
        kappa_th_center_snr, peak_counts_single = summary_statistics.compute_single_scale_peak_counts(snr, KAPPA_SNR)

        # Append peak counts for this tile to the list
        peak_counts_realization.append(peak_counts_single)

    # Compute the average peak counts vector for this realization
    average_peak_counts = np.mean(peak_counts_realization, axis=0)

    # Append the average peak counts vector for this realization to the list of data vectors
    data_vectors.append(average_peak_counts)

# Convert the list of data vectors into a NumPy array
data_vectors = np.array(data_vectors)
# save the data vectors as a .npy file
np.save('.././output/data_vector_SS_PC.npy', data_vectors)

# Compute the average peak counts vector across all realizations
mean_PC_over_realizations = np.mean(data_vectors, axis=0)

# Compute the deviations of peak counts in each realization from the average vector
deviations = data_vectors - mean_PC_over_realizations

# Compute the covariance matrix
# The covariance_matrix will be a square matrix of shape (num_bins, num_bins).
num_realizations = data_vectors.shape[0]
covariance_matrix = np.dot(deviations.T, deviations) / (num_realizations - 1)
# save the covariance matrix as a .npy file
np.save('.././output/covariance_matrix_SS_PC.npy', covariance_matrix)

# Normalize the covariance matrix to the unity of the diagonal
diagonal_sqrt = np.sqrt(np.diag(covariance_matrix))
covariance_matrix_normalized = covariance_matrix / np.outer(diagonal_sqrt, diagonal_sqrt)

# # Normalize the covariance matrix
# max_cov_value = np.max(covariance_matrix)
# covariance_matrix_normalized = covariance_matrix / max_cov_value

# Plot the covariance matrix
plotting.plot_map(covariance_matrix_normalized, title='Covariance Matrix of Peak Counts', cmap='viridis', vmin=None, vmax=None)
# save the figure as pdf
plt.savefig('.././output/covariance_matrix_SS_PC.pdf', bbox_inches='tight')
