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
from multiprocessing import Pool
from glob import glob

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

def make_shear_map(CATALOG_FILE, add_noise=True):
    catalog_data = slics.read_catalogue_pd(CATALOG_FILE)
    ra = catalog_data['RA']
    dec = catalog_data['Dec']
    g1_sim = catalog_data['gamma1_sim']
    g2_sim = catalog_data['gamma2_sim']
    
    x, y = radec2xy(np.mean(ra), np.mean(dec), ra, dec)
    Nx, Ny = int(SIZE_X_DEG / PIX_ARCMIN * 60), int(SIZE_Y_DEG / PIX_ARCMIN * 60)
    galmap = bin2d(x, y, npix=(Nx, Ny))
    mask = (galmap > 0).astype(int)

    sigma_noise = np.zeros_like(galmap)
    sigma_noise[mask != 0] = SHAPE_NOISE / np.sqrt(2 * galmap[mask != 0])
    sigma_noise[mask == 0] = np.max(sigma_noise[mask != 0])
    # noise_map = sigma_noise * np.random.randn(sigma_noise.shape[0], sigma_noise.shape[1])

    e1map, e2map = bin2d(x, y, npix=(Nx, Ny), v=(g1_sim, g2_sim))

    if add_noise:
        # Add noise only if requested
        noise_e1 = np.random.randn(*e1map.shape) * sigma_noise
        noise_e2 = np.random.randn(*e2map.shape) * sigma_noise
        e1map_noisy = e1map + noise_e1 * mask
        e2map_noisy = e2map + noise_e2 * mask
        return e1map_noisy, e2map_noisy, mask, sigma_noise
    else:
        # Return the maps without added noise
        return e1map, e2map, mask, sigma_noise


def make_mass_map(e1map, e2map, mask, sigma_noise, method='ks'):
    d = shear_data()
    d.g1 = e1map
    d.g2 = -e2map
    (nx, ny) = e1map.shape
    d.mask = mask
    Ncov = np.zeros((nx, ny))
    Ncov[mask > 0] = 2. * sigma_noise[mask > 0]**2
    Ncov[mask == 0] = 1e9
    d.Ncov = Ncov
    d.nx = nx
    d.ny = ny

    if method == 'ks':
        M = massmap2d(name='mass')
        M.init_massmap(d.nx, d.ny)
        M.DEF_niter = 50
        M.niter_debias = 30
        M.Verbose = False
        ks = M.gamma_to_cf_kappa(e1map, -e2map)
        ks = ks.real
        return ks


def summary_statistics(ks_noisy, sigma_noise, mask, nscales=NSCALES, min_snr=MIN_SNR, max_snr=MAX_SNR, nbins=NBINS, nbins_l1=NBINS_L1):
    nx, ny = ks_noisy.shape
    WT = starlet2d(gen2=False, l2norm=False, verb=False)
    WT.init_starlet(nx, ny, nscale=nscales)
    H = HOS_starlet_l1norm_peaks(WT)
    H.set_bins(Min=min_snr, Max=max_snr, nbins=nbins)
    H.set_data(ks_noisy, SigmaMap=sigma_noise, Mask=mask)
    H.get_mono_scale_peaks(ks_noisy, sigma_noise, mask=mask)
    H.get_wtpeaks(Mask=mask)
    pc = H.Peaks_Count
    H.get_wtl1(nbins_l1*2, Mask=mask)

    H.plot_mono_peaks_histogram()
    H.plot_peaks_histo(log_scale=True)
    H.plot_l1norm()
    
    return H.Mono_Peaks_Count, H.Peaks_Count, H.l1norm

def process_tile(filename, mass_mapping_method='ks', add_noise=False, save_mass_map=False, mass_map_output_file=None):
    e1map, e2map, mask, sigma_noise = make_shear_map(filename, add_noise=add_noise)
    ks = make_mass_map(e1map, e2map, mask, sigma_noise, method=mass_mapping_method)
    ks_noisy = ks
    peaks_mono, peaks_multi, l1norm = summary_statistics(ks_noisy, sigma_noise, mask)

    if save_mass_map:
        np.save(mass_map_output_file, ks*mask)

    return peaks_mono, peaks_multi, l1norm


def worker(args):
    filename, mass_mapping_method, add_noise, save_mass_map, output_dir = args
    # Construct output filename for the mass map
    base_name, file_ext = os.path.splitext(os.path.basename(filename))
    new_file_ext = '.npy'
    mass_map_output_file = os.path.join(output_dir, f"{base_name}_{mass_mapping_method}{new_file_ext}") if save_mass_map else None
    
    # Simulate processing the tile
    summary_statistics = process_tile(filename, mass_mapping_method, add_noise, save_mass_map, mass_map_output_file)
    return summary_statistics


def process_footprint(file_list_path, output_dir=None, mass_mapping_method='ks', add_noise=True, save_mass_map=False, num_processes=19):
    """
    Processes a footprint by parallelizing over the number of tiles in a footprint
    and returns a single averaged vector for each summary statistic across all tiles.
    """
    if save_mass_map and output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    with open(file_list_path, 'r') as file:
        filenames = file.read().splitlines()
    
    args = [(filename, mass_mapping_method, add_noise, save_mass_map, output_dir) for filename in filenames]
    
    with Pool(num_processes) as pool:
        results = pool.map(worker, args)
    
    # Initialize lists to hold each type of summary statistic separately
    SS_PC_data = []
    MS_PC_data = []
    l1_norm_data = []
    
    # Iterate over the results to collect the statistics
    for SS_PC, MS_PC, l1_norm in results:
        SS_PC_data.append(SS_PC)
        MS_PC_data.append(MS_PC)
        l1_norm_data.append(l1_norm)
    
    # Calculate the mean of each summary statistic across all tiles
    SS_PC_mean = np.mean(SS_PC_data, axis=0)
    MS_PC_mean = np.mean(MS_PC_data, axis=0)
    l1_norm_mean = np.mean(l1_norm_data, axis=0)
    
    # Return the averaged statistics
    return SS_PC_mean, MS_PC_mean, l1_norm_mean

def process_cosmo(cosmology, cosmo_dir, output_dir, mass_mapping_method='ks', add_noise=True, save_mass_map=False, num_processes=19):
    """
    Processes files for a specific cosmology and bin by parallelizing over the tiles in a footprint
    and saves a 10-row table (2 seeds x 5 LOS) for each bin.
    """
    if save_mass_map and output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    # Define bins, seeds, and LOS based on your setup
    bins = ['Bin1', 'Bin2', 'Bin3', 'Bin4']
    seeds = ['a', 'f']  # Example seeds
    LOS_range = range(1, 6)  # Example LOS 1 through 5
    
    for bin in bins:
        all_results = []  # Reset for each bin
        
        for seed in seeds:
            for LOS in LOS_range:
                # Find the specific file for each combination of cosmology, bin, seed, LOS
                file_name_pattern = f"{cosmology}_{seed}_LOS{LOS}_{bin}.txt"
                file_list_path = glob(os.path.join(cosmo_dir, file_name_pattern))
                
                # Skip if no file matches the pattern
                if not file_list_path:
                    print(f"File not found for pattern: {file_name_pattern}")
                    continue
                
                # Assuming only one file matches, use the first one found
                SS_PC_mean, MS_PC_mean, l1_norm_mean = process_footprint(
                    file_list_path[0],
                    output_dir,
                    mass_mapping_method,
                    add_noise,
                    save_mass_map,
                    num_processes
                )
                
                # Append the results with full vectors preserved
                all_results.append((seed, LOS, SS_PC_mean, MS_PC_mean, l1_norm_mean))
        
        # Convert the results to a structured numpy array
        dtype = [('seed', 'U10'), ('LOS', int),
                 ('SS_PC_mean', np.object_), ('MS_PC_mean', np.object_), ('l1_norm_mean', np.object_)]
        results_array = np.array(all_results, dtype=dtype)
        
        # Save the results array for the current bin
        save_path = os.path.join(output_dir, f"{cosmology}_{bin}_{mass_mapping_method}_run.npy")
        np.save(save_path, results_array)
        print(f"Saved: {save_path}")


# Example usage
cosmology = '19'  # Specify the cosmology you're processing
cosmo_dir = "../output/tiles_lists/"  # Directory containing the .txt files
output_directory = '/n17data/tersenov/SLICS/Cosmo_DES/summary_stats'  # Directory to save the summary stats

# Process the specified cosmology and save the results for each bin
process_cosmo(
    cosmology,
    cosmo_dir,
    output_directory,
    'ks',  # Example mass mapping method
    True,  # Add noise
    True,  # Save mass maps
    19  # Number of processes
)





