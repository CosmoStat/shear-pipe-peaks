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
from collections import defaultdict
from multiprocessing import Pool
import argparse
import pandas as pd


# CONSTANTS AND PARAMETERS
N_GAL = 7 
SIZE_X_DEG = 10.
SIZE_Y_DEG = 10.
PIX_ARCMIN = 1.
SHAPE_NOISE = 0.44
NSCALES = 6
MIN_SNR = -2
MAX_SNR = 6
NBINS = 21
NBINS_L1 = 20


def merge_catalogs(file_paths, all_col=True):
    """
    Merges catalogs from given file paths into a single DataFrame, using a custom reading function.

    Parameters:
    file_paths : list[str]
        The file paths of the catalogs to be merged.
    all_col : bool
        Flag indicating whether to read all columns or just essential ones.

    Returns:
    pandas.DataFrame
        The combined catalog.
    """
    combined_df = pd.DataFrame()
    for file_path in file_paths:
        df = slics.read_catalogue_pd(file_path, all_col=all_col)
        combined_df = pd.concat([combined_df, df], ignore_index=True)
    
    return combined_df

def generate_file_lists(master_file_path, cosmology):
    """
    Generates lists of file paths for each unique combination of seed, LOS, and bin
    for the specified cosmology by reading from a master file.

    Parameters:
    master_file_path : str
        Path to the master file containing all filenames.
    cosmology : str
        The specific cosmology to process.

    Returns:
    dict
        A dictionary where each key is a tuple (seed, LOS, bin) and the value
        is the list of file paths for that combination.
    """
    groups = defaultdict(list)
    with open(master_file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue

            parts = line.split('/')
            filename = parts[-1]
            file_parts = filename.split('_')
            
            # Adjust these indices according to your file naming conventions
            file_cosmology = file_parts[2]
            seed = file_parts[3]
            LOS = file_parts[5]
            bin_part = file_parts[6]
            
            if file_cosmology == cosmology:
                key = (seed, LOS, bin_part)
                groups[key].append(line)

    return groups


def make_shear_map(CATALOG_FILE, add_noise=True, random_seed=None):
    """
    Generates shear maps from a catalog file, with an option to add Gaussian noise.

    Parameters:
    CATALOG_FILE : str
        Path to the catalog file containing galaxy data.
    add_noise : bool, optional
        Determines if noise should be added to the shear maps. Default is True.

    Returns:
    Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
        A tuple containing e1map, e2map, mask, and sigma_noise arrays.
    """
    catalog_data = CATALOG_FILE
    ra = catalog_data['RA']
    dec = catalog_data['Dec']
    g1_sim = catalog_data['gamma1_sim']
    g2_sim = catalog_data['gamma2_sim']
    
    x, y = radec2xy(np.mean(ra), np.mean(dec), ra, dec)
    # convert x,y from rad to arcmin
    x, y = x * 60 * 180 / np.pi, y * 60 * 180 / np.pi
    
    Nx, Ny = int(SIZE_X_DEG / PIX_ARCMIN * 60), int(SIZE_Y_DEG / PIX_ARCMIN * 60)
    galmap = bin2d(x, y, npix=(Nx, Ny))
    mask = (galmap > 0).astype(int)

    # variable noise level
    sigma_noise = np.zeros_like(galmap)
    sigma_noise[mask != 0] = SHAPE_NOISE / np.sqrt(2 * galmap[mask != 0])
    sigma_noise[mask == 0] = np.max(sigma_noise[mask != 0])

    e1map, e2map = bin2d(x, y, npix=(Nx, Ny), v=(g1_sim, g2_sim))
    
    # Set random seed if provided
    if random_seed is not None:
        np.random.seed(random_seed)

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
    """
    Creates a mass map from shear maps using the specified mass mapping method.

    Parameters:
    e1map, e2map : np.ndarray
        Arrays of the first and second shear components.
    mask : np.ndarray
        A binary mask indicating the presence of galaxies.
    sigma_noise : np.ndarray
        The noise level in the shear measurements.
    method : str, optional
        The mass mapping method to use. Default is 'ks' (Kaiser-Squires).

    Returns:
    np.ndarray
        The generated mass map.
    """
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
        M.Verbose = False
        ks = M.gamma_to_cf_kappa(e1map, -e2map)
        ks = ks.real
        return ks
    
    if method == 'ksi':
        M = massmap2d(name='mass')
        M.init_massmap(d.nx, d.ny)
        M.DEF_niter = 50
        M.niter_debias = 30
        M.Verbose = False
        ksi =  M.iks(d.g1, d.g2, mask) 
        ksi = ksi.real
        return ksi

    if method == 'wiener':
        M = massmap2d(name='mass')
        M.init_massmap(d.nx, d.ny)
        M.DEF_niter = 50
        M.niter_debias = 30
        M.Verbose = False
        pn = readfits('/home/tersenov/shear-pipe-peaks/input/exp_wiener_miceDSV_noise_powspec.fits')
        ps1d = readfits('/home/tersenov/shear-pipe-peaks/input/exp_wiener_miceDSV_signal_powspec.fits')
        d.ps1d = ps1d
        ke_inp_pwiener, kb_winp = M.prox_wiener_filtering(InshearData=d, PowSpecSignal=d.ps1d, Pn=pn, niter=M.DEF_niter, Inpaint=True) 
        return ke_inp_pwiener
    
    if method == 'mca':
        M = massmap2d(name='mass')
        M.init_massmap(d.nx, d.ny)
        M.DEF_niter = 30
        M.Verbose = False
        ps1d = readfits('/home/tersenov/shear-pipe-peaks/input/exp_wiener_miceDSV_signal_powspec.fits')
        d.ps1d = ps1d        
        mcalens, _, _, _ = M.sparse_wiener_filtering(d, d.ps1d, Nsigma=3, niter=M.DEF_niter, Inpaint=True, Bmode=True)
        return mcalens
    
    

def summary_statistics(mass_map, sigma_noise, mask, nscales=NSCALES, min_snr=MIN_SNR, max_snr=MAX_SNR, nbins=NBINS, nbins_l1=NBINS_L1):
    """
    Computes summary statistics from a noisy kappa map using wavelet transforms.

    Parameters:
    mass_map : np.ndarray
        (Noisy) kappa map.
    sigma_noise : np.ndarray
        Noise level in the kappa map.
    mask : np.ndarray
        Binary mask indicating the observational field.
    nscales : int
        Number of wavelet scales to use.
    min_snr, max_snr : float
        Minimum and maximum signal-to-noise ratios for peak detection.
    nbins : int
        Number of bins for histogramming peaks.
    nbins_l1 : int
        Number of bins for the L1-norm histogram.

    Returns:
    Tuple[np.ndarray, np.ndarray, np.ndarray]
        Arrays of the computed Mono_Peaks_Count, Peaks_Count, and l1norm.
    """
    nx, ny = mass_map.shape
    WT = starlet2d(gen2=False, l2norm=False, verb=False)
    WT.init_starlet(nx, ny, nscale=nscales)
    H = HOS_starlet_l1norm_peaks(WT)
    H.set_bins(Min=min_snr, Max=max_snr, nbins=nbins)
    H.set_data(mass_map, SigmaMap=sigma_noise, Mask=mask)
    H.get_mono_scale_peaks(mass_map, sigma_noise, smoothing_sigma=12, mask=mask)
    H.get_wtpeaks(Mask=mask)
    pc = H.Peaks_Count
    H.get_wtl1(nbins_l1*2, Mask=mask, min_snr=-6, max_snr=6)

    return H.Mono_Peaks_Count, H.Peaks_Count, H.l1norm

def process_tile(catalog, mass_mapping_method='ks', add_noise=True, save_mass_map=False, mass_map_output_file=None, random_seed=None):
    
    e1map, e2map, mask, sigma_noise = make_shear_map(catalog, add_noise=add_noise, random_seed=random_seed)
    mass_map = make_mass_map(e1map, e2map, mask, sigma_noise, method=mass_mapping_method)
    peaks_mono, peaks_multi, l1norm = summary_statistics(mass_map, sigma_noise, mask)

    if save_mass_map:
        np.save(mass_map_output_file, mass_map*mask)

    return peaks_mono, peaks_multi, l1norm

def worker(args):
    
    catalog, mass_mapping_method, add_noise, save_mass_map, run_number, random_seed = args
    
    # Construct output filename for the mass map
    base_name = '/n17data/tersenov/SLICS/Cosmo_DES/summary_stats/combined_bins/mass_maps/'
    new_file_ext = '.npy'
    mass_map_output_file = f"{base_name}merged_catalog_{mass_mapping_method}_run{run_number}_rs{random_seed}{new_file_ext}" if save_mass_map else None

    # Simulate processing the tile
    summary_statistics = process_tile(catalog, mass_mapping_method, add_noise, save_mass_map=False, mass_map_output_file=mass_map_output_file, random_seed=random_seed)
    
    return summary_statistics

def process_footprint(catalog_list, output_dir=None, mass_mapping_method='ks', add_noise=True, save_mass_map=False, num_processes=19, run_number=0, random_seed=None):

    if save_mass_map and output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    args = [(catalog, mass_mapping_method, add_noise, save_mass_map, run_number, random_seed) for catalog in catalog_list]
    
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

def process_cosmo(cosmology, master_file_path, output_dir, mass_mapping_method='ks', add_noise=True, save_mass_map=False, num_processes=19, run_number=0, random_seed=None):

    if save_mass_map and output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    file_lists = generate_file_lists(master_file_path, cosmology)

    bins = ['Bin1', 'Bin2', 'Bin3', 'Bin4']
    seeds = ['a', 'f'] 
    LOSs = ['LOS1','LOS2','LOS3','LOS4','LOS5'] 

    all_results = []  # Store results for each bin here

    for seed in seeds:
        for LOS in LOSs:
            file_paths_for_every_bin = []
            footprint_catalogs = []
            
            for bin in bins:
                
                # Generate the key for looking up in the file lists
                key = (seed, bin, LOS)
                file_paths = file_lists.get(key, [])
                
                if not file_paths:
                    print(f"File not found for pattern: {key}")
                    continue  # Skip if no files for this combination
                
                # append the results for the current bin
                file_paths_for_every_bin.append(file_paths)

            
            for i in range(19):
                # merge the catalogs
                combined_catalog = merge_catalogs(np.array(file_paths_for_every_bin)[:,i][:])
                # add the combined catalog to the list of all catalogs for the footprint
                footprint_catalogs.append(combined_catalog)
                
            SS_PC_mean, MS_PC_mean, l1_norm_mean = process_footprint(
                # combined_catalog,
                footprint_catalogs,
                output_dir,
                mass_mapping_method,
                add_noise,
                save_mass_map,
                num_processes,
                run_number,
                random_seed
            )
            
            all_results.append((seed, LOS, SS_PC_mean, MS_PC_mean, l1_norm_mean))
                

    # Convert the results to a structured numpy array
    dtype = [('seed', 'U10'), ('LOS', 'U10'),
                ('SS_PC_mean', np.object_), ('MS_PC_mean', np.object_), ('l1_norm_mean', np.object_)]
    results_array = np.array(all_results, dtype=dtype)
    
    # Save the results array for the current bin
    save_path = os.path.join(output_dir, f"{cosmology}_{mass_mapping_method}_run{run_number}_rs{random_seed}.npy")
    np.save(save_path, results_array)
    print(f"Saved: {save_path}")



if __name__ == "__main__":
    cosmology = 'fid'
    master_file_path = '/home/tersenov/shear-pipe-peaks/input/master_file.txt'  # Adjust this path as necessary
    output_directory = '/n17data/tersenov/SLICS/Cosmo_DES/summary_stats/combined_bins'  # Adjust this path as necessary

    process_cosmo(
        cosmology,
        master_file_path,
        output_directory,
        'ksi',
        True,  # Consider if you want these as command-line arguments too
        True,
        19,
        3,
        0
    )
