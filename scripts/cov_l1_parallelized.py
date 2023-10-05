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
PIX_ARCMIN = 0.4
SHAPE_NOISE = 0.44
NSCALES = 5
NBINS = 40 
KAPPA_SNR = np.linspace(-2, 6, 31)
NUM_REALIZATIONS = 124  # Number of realizations

# Initialize an empty list to store data vectors for the L1-norm histogram for each realization
L1_norm_data_vectors = []


# function to compute L1-norm histograms for a single tile
def compute_l1_norm_histogram(tile_file):
    catalog_data = slics.read_catalogue_pd(tile_file)

    # Extract data from the catalog 
    ra = catalog_data['RA']
    dec = catalog_data['Dec']
    g1_sim = catalog_data['gamma1_sim']
    g2_sim = catalog_data['gamma2_sim']

    # Calculate shear and kappa maps
    e1map, e2map = mapping.bin_shear_field(ra, dec, g1_sim, g2_sim)
    kappaE, _ = ks93(e1map, -e2map)

    # Add noise to the kappa map
    kappaE_noisy, noise_map_CFIS_z05 = mapping.add_noise_to_kappa_map(kappaE, SHAPE_NOISE, N_GAL, PIX_ARCMIN)

    # Compute the L1-norm histogram for each scale
    bins_l1, l1norm_histogram = get_l1norm_noisy(kappaE, noise_map_CFIS_z05, NSCALES, NBINS*2)

    return l1norm_histogram

# Loop over realizations
for realization_files in collections_of_files:
    # Create a pool of worker processes
    with Pool(processes=num_tiles_per_realization) as pool:
        # Compute L1-norm histograms in parallel for each tile in this realization
        l1_norm_histogram_realization = pool.map(compute_l1_norm_histogram, realization_files)

    # Compute the average L1-norm histogram vector for this realization
    average_l1_norm_histogram = np.mean(l1_norm_histogram_realization, axis=0)

    # Append the average L1-norm histogram vector for this realization to the list of data vectors
    L1_norm_data_vectors.append(average_l1_norm_histogram)


# Convert the list of data vectors into a NumPy array
L1_norm_data_vectors = np.array(L1_norm_data_vectors)
# save the data vectors as a .npy file
np.save('.././output/L1_norm_data_vector.npy', L1_norm_data_vectors)

# Reshape the L1-norm data vectors 
L1_norm_data_vectors_reshaped = L1_norm_data_vectors.reshape(NUM_REALIZATIONS, -1)

# Compute the average L1-norm histogram vector across all realizations for the L1-norm
mean_L1_norm_over_realizations = np.mean(L1_norm_data_vectors_reshaped, axis=0)

# Compute the deviations of L1-norm histograms in each realization from the average vector
deviations_L1_norm = L1_norm_data_vectors_reshaped - mean_L1_norm_over_realizations

# Compute the covariance matrix for the L1-norm histograms
num_realizations_L1_norm, num_bins_L1_norm = L1_norm_data_vectors_reshaped.shape
covariance_matrix_L1_norm = np.dot(deviations_L1_norm.T, deviations_L1_norm) / (num_realizations_L1_norm - 1)
# save the covariance matrix as a .npy file
np.save('.././output/covariance_matrix_l1_norm.npy', covariance_matrix_L1_norm)

# Calculate the diagonal of the covariance matrix
diagonal_L1_norm = np.sqrt(np.diag(covariance_matrix_L1_norm))

# Check for zero values and replace them with a small positive value
diagonal_L1_norm[diagonal_L1_norm == 0] = 1e-10

# Calculate the correlation coefficients
correlation_matrix_L1_norm = covariance_matrix_L1_norm / np.outer(diagonal_L1_norm, diagonal_L1_norm)

plotting.plot_map(correlation_matrix_L1_norm, title='l1-norm', cmap='viridis', vmin=None, vmax=None)
# save the figure as pdf
plt.savefig('.././output/covariance_matrix_l1_norm.pdf')
