import random
import sys  
import os
import numpy as np
from numpy import linalg as LA
from scipy import ndimage as ndi
from scipy.special import erf
import matplotlib.pyplot as plt
from lenspack.utils import bin2d
from lenspack.geometry.projections.gnom import radec2xy
from lenspack.image.inversion import ks93
from lenspack.geometry.projections import gnom
import lenspack.peaks as peaks
from lenspack.image.transforms import starlet2d
from pycs.astro.wl.mass_mapping import *
from pycs.misc.im_isospec import *
from pycs.sparsity.sparse2d.starlet import *
from pycs.misc.cosmostat_init import *
from pycs.misc.mr_prog import *
from pycs.misc.utilHSS import *
from pycs.misc.im1d_tend import *
from pycs.misc.stats import *
from pycs.sparsity.sparse2d.dct import dct2d, idct2d
from pycs.sparsity.sparse2d.dct_inpainting import dct_inpainting
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
NBINS=31

NUM_REALIZATIONS = 124  # Number of realizations

# Initialize an empty list to store data vectors for each realization
data_vectors = []

# Function to compute peak counts for a single tile
def compute_peak_counts(tile_file):

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
    H.get_mono_scale_peaks(ks_noisy, sigma_noise, mask=mask)

    peak_counts_single = H.Mono_Peaks_Count

    return peak_counts_single



# Loop over realizations
for realization_files in collections_of_files: 
    # Create a pool of worker processes
    with Pool(processes=num_tiles_per_realization) as pool:
        # Compute peak counts in parallel for each tile in this realization
        peak_counts_realization = pool.map(compute_peak_counts, realization_files)

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
