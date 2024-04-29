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
from scipy import ndimage

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
    catalog_data = slics.read_catalogue_pd(CATALOG_FILE)
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
    
    # # constant noise level
    # sigma_noise = np.ones_like(galmap) * SHAPE_NOISE / np.sqrt(2 * N_GAL)
    # noise_map = sigma_noise * np.random.randn(sigma_noise.shape[0], sigma_noise.shape[1])

    e1map, e2map = bin2d(x, y, npix=(Nx, Ny), v=(g1_sim, g2_sim)) 

    if add_noise:
        # Add noise only if requested
        noise_e1 = np.random.randn(*e1map.shape) * sigma_noise
        noise_e2 = np.random.randn(*e2map.shape) * sigma_noise
        e1map_noisy = e1map + noise_e1 * mask
        e2map_noisy = e2map + noise_e2 * mask
        
        # # smooth the noisy maps with a gaussian filter 
        # sm_e1map_noisy = ndimage.gaussian_filter(e1map_noisy, 2)
        # sm_e2map_noisy = ndimage.gaussian_filter(e2map_noisy, 2)

        # # the final e1map will be the e1map_noisy where the mask is 1 and the smoothed e1map where the mask is 0
        # e1map_noisy = e1map_noisy * mask + sm_e1map_noisy * (1 - mask)
        # e2map_noisy = e2map_noisy * mask + sm_e2map_noisy * (1 - mask)
        
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
        ke_inp_pwiener, kb_winp = M.prox_wiener_filtering(d.g1, d.g2, d.ps1d, d.Ncov, Pn=pn, Inpaint=True) #, Pn=Pn) # ,ktr=InShearData.ktr)
        return ke_inp_pwiener

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
    H.get_mono_scale_peaks(mass_map, sigma_noise, smoothing_sigma=6, mask=mask)
    H.get_wtpeaks(Mask=mask)
    pc = H.Peaks_Count
    H.get_wtl1(nbins_l1*2, Mask=mask, min_snr=-6, max_snr=6)

    return H.Mono_Peaks_Count, H.Peaks_Count, H.l1norm

def process_tile(filename, mass_mapping_method='ks', add_noise=False, save_mass_map=False, mass_map_output_file=None):
    """
    Processes a single tile: generates shear and mass maps and computes summary statistics.

    Parameters:
    filename : str
        Path to the catalog file for the tile.
    mass_mapping_method : str
        Method used for mass mapping.
    add_noise : bool
        If True, noise is added to the shear maps.
    save_mass_map : bool
        If True, the generated mass map is saved to disk.
    mass_map_output_file : str, optional
        Path where the mass map should be saved, if applicable.

    Returns:
    Tuple[np.ndarray, np.ndarray, np.ndarray]
        Tuple of summary statistics: Mono_Peaks_Count, Peaks_Count, l1norm.
    """
    e1map, e2map, mask, sigma_noise = make_shear_map(filename, add_noise=add_noise)
    mass_map = make_mass_map(e1map, e2map, mask, sigma_noise, method=mass_mapping_method)
    peaks_mono, peaks_multi, l1norm = summary_statistics(mass_map, sigma_noise, mask)

    if save_mass_map:
        np.save(mass_map_output_file, mass_map*mask)

    return peaks_mono, peaks_multi, l1norm


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

def worker(args):
    """
    Worker function for processing a single tile in parallel.

    Parameters:
    args : tuple
        Arguments to pass to process_tile, including filename, mass mapping method,
        add_noise flag, save_mass_map flag, and mass_map_output_file path.

    Returns:
    Tuple[np.ndarray, np.ndarray, np.ndarray]
        The summary statistics from processing the tile.
    """
    filename, mass_mapping_method, add_noise, save_mass_map, run_number = args
    
    # Construct output filename for the mass map
    base_name, file_ext = os.path.splitext(filename)
    new_file_ext = '.npy'
    mass_map_output_file = f"{base_name}_{mass_mapping_method}_run{run_number}{new_file_ext}" if save_mass_map else None
    
    # Simulate processing the tile
    summary_statistics = process_tile(filename, mass_mapping_method, add_noise, save_mass_map, mass_map_output_file)
    return summary_statistics

def process_footprint(file_list_path, output_dir=None, mass_mapping_method='ks', add_noise=True, save_mass_map=False, num_processes=19, run_number=0):
    """
    Processes a set of tiles (a footprint) in parallel, averaging summary statistics.

    Parameters:
    file_paths : list[str]
        List of paths to catalog files to process.
    output_dir : str, optional
        Directory where output files should be saved.
    mass_mapping_method : str
        Method used for mass mapping.
    add_noise : bool
        If True, noise is added to the shear maps.
    save_mass_map : bool
        If True, mass maps are saved.
    num_processes : int
        Number of parallel processes to use.

    Returns:
    Tuple[np.ndarray, np.ndarray, np.ndarray]
        Averaged summary statistics across all processed tiles.
    """
    if save_mass_map and output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    args = [(filename, mass_mapping_method, add_noise, save_mass_map, run_number) for filename in file_list_path]
    print(args)
    
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

def process_cosmo(cosmology, master_file_path, output_dir, mass_mapping_method='ks', add_noise=True, save_mass_map=False, num_processes=19, run_number=0):
    """
    Processes all tiles for a specific cosmology, generating and saving summary statistics
    for each bin in a specified directory.

    Parameters:
    cosmology : str
        Identifier for the cosmology being processed.
    master_file_path : str
        Path to the master file listing all catalog files.
    output_dir : str
        Directory where summary statistics should be saved.
    mass_mapping_method : str
        Mass mapping method to use.
    add_noise : bool
        If True, noise is added to shear maps.
    save_mass_map : bool
        If True, saves the mass maps alongside the catalog files.
    num_processes : int
        Number of processes for parallel execution.
    run_number : int
        Run identifier.

    Returns:
    None
    """
    if save_mass_map and output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    file_lists = generate_file_lists(master_file_path, cosmology)

    bins = ['Bin1', 'Bin2', 'Bin3', 'Bin4']
    seeds = ['a', 'f'] 
    LOSs = ['LOS1','LOS2','LOS3','LOS4','LOS5'] 

    for bin in bins:
        all_results = []  # Store results for each bin here
        
        for seed in seeds:
            for LOS in LOSs:
                # Generate the key for looking up in the file lists
                key = (seed, bin, LOS)
                file_paths = file_lists.get(key, [])
                
                if not file_paths:
                    print(f"File not found for pattern: {key}")
                    continue  # Skip if no files for this combination
                
                SS_PC_mean, MS_PC_mean, l1_norm_mean = process_footprint(
                    file_paths,
                    output_dir,
                    mass_mapping_method,
                    add_noise,
                    save_mass_map,
                    num_processes,
                    run_number
                )
                
                all_results.append((seed, LOS, SS_PC_mean, MS_PC_mean, l1_norm_mean))

        # Convert the results to a structured numpy array
        dtype = [('seed', 'U10'), ('LOS', 'U10'),
                 ('SS_PC_mean', np.object_), ('MS_PC_mean', np.object_), ('l1_norm_mean', np.object_)]
        results_array = np.array(all_results, dtype=dtype)
        
        # Save the results array for the current bin
        save_path = os.path.join(output_dir, f"{cosmology}_{bin}_{mass_mapping_method}_run{run_number}.npy")
        np.save(save_path, results_array)
        print(f"Saved: {save_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process data for a given cosmology using a specific mass mapping method.')
    parser.add_argument('cosmology', type=str, help='Identifier for the cosmology being processed.')
    parser.add_argument('master_file_path', type=str, help='Path to the master file listing all catalog files.')
    parser.add_argument('output_dir', type=str, help='Directory where summary statistics should be saved.')
    parser.add_argument('num_tiles', type=int, help='Number of tiles in the footprint.')
    parser.add_argument('mass_mapping_method', type=str, help='Mass mapping method to use.')
    parser.add_argument('run_number', type=int, help='Run identifier.')
    args = parser.parse_args()

    # Assuming master_file_path and output_directory are constants or can be derived from other parameters
    master_file_path = '/home/tersenov/shear-pipe-peaks/input/master_file.txt'  # Adjust this path as necessary
    output_directory = '/n17data/tersenov/SLICS/Cosmo_DES/summary_stats'  # Adjust this path as necessary

    process_cosmo(
        args.cosmology,
        master_file_path,
        output_directory,
        args.mass_mapping_method,
        True,  # Consider if you want these as command-line arguments too
        True,
        args.num_tiles,
        args.run_number
    )
