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
import argparse
import pandas as pd

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

def find_corresponding_files_for_tile(tile_file):
    """
    Find the files corresponding to the other bins of the same LOS and tile.
    
    Parameters
    ----------
    tile_file : str
        The path to the file for the tile of interest.
        
    Returns
    -------
    list of str
        A list of paths to the files corresponding to the other bins of the same LOS and tile.
    """
        
    # Extract the relevant information from the file path
    file_name = tile_file.split("/")[-1]
    file_parts = file_name.split("_")
    bin_number = int(file_parts[4][3:])
    LOS_number = int(file_parts[5][3:])
    tile_number = int(file_parts[6][1:-4])
    
    num_bins = 4
    
    # Find the files corresponding to the other bins of the same LOS and tile
    corresponding_files = [tile_file]
    for bin_num in range(1, num_bins + 1):
        if bin_num != bin_number:
            # find the corresponding file
            corresponding_file = f"/n17data/tersenov/Cov_DES_SLICS/LOS{LOS_number:03d}/DES_MocksCat_SLICS_4_Bin{bin_number:01d}_LOS{LOS_number:03d}_R{tile_number:02d}.dat"
            corresponding_files.append(corresponding_file)
            
    return corresponding_files


def compute_statistics(args):
    
    tile_file, seed, method = args
    
    corresponding_files = find_corresponding_files_for_tile(tile_file)
    
    # combine the corresponding files into a single catalog
    catalog_data = pd.DataFrame()
    for file_path in corresponding_files:
        df = slics.read_catalogue_pd(file_path, all_col=False)
        catalog_data = pd.concat([catalog_data, df], ignore_index=True)

    
    # if the seed is not 0, shuffle the catalog
    if seed != 0:
        # Shuffle the last two columns
        last_two_cols = catalog_data.iloc[:, 2:].sample(frac=1, random_state=seed).reset_index(drop=True)
        shuffled_catalog = catalog_data.copy()
        shuffled_catalog['gamma1_sim'] = last_two_cols['gamma1_sim']
        shuffled_catalog['gamma2_sim'] = last_two_cols['gamma2_sim']
        catalog_data = shuffled_catalog
    
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
    sigma_noise[mask == 0] = np.max(sigma_noise[mask != 0]) # set the noise to the maximum value in the map where there are galaxies       
    
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



def compute_cov_datavectors(realization, master_file_path, output_dir, num_tiles, mass_mapping_method, run):
    
    #change directory to the directory where the script is located
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    
    # Read the file paths from the "master_file_cov.txt"
    with open(master_file_path, "r") as file:
        file_paths = file.readlines()
        file_paths = [path.strip() for path in file_paths]
        
    # Parse these file paths
    parsed_cov_data = slics.parse_cov_SLICS_filenames(file_paths)

    los_numbers = np.unique(parsed_cov_data['LOS']) # List of all LOS numbers
    num_realizations = 124 # Number of realizations
    num_tiles_per_realization = num_tiles # Number of tiles to select for each realization

    num_bins = 4
    bin_number = 1

    # Reconstruct 124 realisations of the survey by picking each tile from a random LOS, ensuring that each file is only included once.
    collections_of_files = slics.survey_realizations_reconstruction(num_realizations, num_tiles_per_realization, bin_number, parsed_cov_data['LOS'], file_paths)
        
    SS_PC_data_vectors = []
    MS_PC_data_vectors = []
    l1_norm_data_vectors = []    

    # Loop over realizations
    for realization_files in collections_of_files:
        # Create a pool of worker processes
        with Pool(processes=num_tiles_per_realization) as pool:
            args = [(file, realization, mass_mapping_method) for file in realization_files]
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
    
    #crete the data vector directories if they do not exist
    if not os.path.exists(os.path.join(output_dir, 'SS_PC')):
        os.makedirs(os.path.join(output_dir, 'SS_PC'))
    if not os.path.exists(os.path.join(output_dir, 'MS_PC')):
        os.makedirs(os.path.join(output_dir, 'MS_PC'))
    if not os.path.exists(os.path.join(output_dir, 'L1_norm')):
        os.makedirs(os.path.join(output_dir, 'L1_norm'))
        
    # # save the data vectors as a .npy file
    # np.save('.././output/covariances/datavectors/data_vector_SS_PC_ksi_merged_bins_sm6_noise044.npy', SS_PC_data_vectors)
    # np.save('.././output/covariances/datavectors/data_vector_MS_PC_cs_ksi_merged_bins_sm6_noise044.npy', MS_PC_data_vectors)
    # np.save('.././output/covariances/datavectors/L1_norm_data_vector_cs_ksi_merged_bins_sm6_noise044.npy', l1_norm_data_vectors)

    np.save(os.path.join(output_dir, f'SS_PC/merged_bins_SS_PC_{mass_mapping_method}_noise_seed{realization}_run{run}.npy'), SS_PC_data_vectors)    
    np.save(os.path.join(output_dir, f'MS_PC/merged_bins_MS_PC_{mass_mapping_method}_noise_seed{realization}_run{run}.npy'), MS_PC_data_vectors)
    np.save(os.path.join(output_dir, f'L1_norm/merged_bins_L1_norm_{mass_mapping_method}_noise_seed{realization}_run{run}.npy'), l1_norm_data_vectors)
    



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process data from covariance dataset for a given noise realization and mass mapping method')
    parser.add_argument('realization', type=str, help='Identifier for the noise realization seed being processed.')
    parser.add_argument('master_file_path', type=str, help='Path to the master file listing all catalog files.')
    parser.add_argument('output_dir', type=str, help='Directory where summary statistics should be saved.')
    parser.add_argument('num_tiles', type=int, help='Number of tiles in the footprint.')
    parser.add_argument('mass_mapping_method', type=str, help='Mass mapping method to use.')
    parser.add_argument('run_number', type=int, help='Run identifier.')
    args = parser.parse_args()

    # output_directory = '/n17data/tersenov/SLICS/Cosmo_DES/summary_stats'  # Adjust this path as necessary
    master_file_path = '/home/tersenov/shear-pipe-peaks/input/master_file_cov.txt'
    
    SEED = int(args.realization)
    # set the seed for the random number generator to be used in the script
    np.random.seed(SEED)


    compute_cov_datavectors(
        SEED,
        master_file_path,
        args.output_dir,
        args.num_tiles,
        args.mass_mapping_method,
        args.run_number
    )