"""SLICS.

:Name: slics.py

:Description: This package contains methods to interface with
    the (Cosmo-)SLICS simulations (Hernois-Deraps et al).

:Authors: Martin Kilbinger <martin.kilbinger@cea.fr> Andreas Tersenov <atersenov@physics.uoc.gr>
"""
import os
from astropy.io import ascii
import numpy as np



col_names = [
    "RA",
    "Dec",
    "e1_data",
    "e2_data",
    "w",
    "redshift_true_sim",
    "gamma1_sim",
    "gamma2_sim",
    "kappa_sim",
    "S_metacal_data",
]

col_names_essential = [
    "RA",
    "Dec",
    "e1_data",
    "e2_data",
    "redshift_true_sim",
]

def read_catalogue(file_path, all_col=True):
    """Read Catalogue.

    Read SLICS catalogue and return content.

    Parameters
    ----------
    file_path : str
        input file path
    all_col : bool, optional
        if True returns all columns; if False (default) only
        essential columns

    Raises
    ------
    IOError
        if input file not found         

    Returns
    -------
    dict
        catalogue content

    """
    if not os.path.exists(file_path):
        raise IOError(f"SLICS catlogue file {file_path} does not exist")

    if all_col:
        include_names = None
    else:
        include_names = col_names_essential

    dat = ascii.read(file_path, names=col_names, include_names=include_names)

    return dat

def parse_SLICS_filenames(file_paths):
    """Reads the information in the filename.

    Parameters
    ----------
    file_paths : list of str (.txt file)
        List containing the paths to the files to be processed

    Returns
    -------
    data : recarray
        Numpy recarray containing the information extracted from the file names

    """

    # Make empty recarray to store the data
    data = np.recarray(len(file_paths), dtype=[('id', int), ('seed', 'U1'), ('bin', int), ('LOS', int), ('tile', int)])

    # Iterate over the file paths and process each file
    for i, file_path in enumerate(file_paths):
        # Extract the file name from the file path
        file_name = file_path.split("/")[-1]
        
        # Split file name into parts
        file_parts = file_name.split("_")
        
        id = int(file_parts[2])
        seed = file_parts[3]
        unknown_number = int(file_parts[4])
        bin = int(file_parts[5][3:])  # Extract the number after "Bin"
        LOS = int(file_parts[6][3:])  # Extract the number after "LOS"
        tile = int(file_parts[7][1:-4])  # Extract the number after "R"

        # Assign the extracted data to the corresponding fields in the recarray
        data[i]['id'] = id
        data[i]['seed'] = seed
        data[i]['bin'] = bin
        data[i]['LOS'] = LOS
        data[i]['tile'] = tile
        # print(data[i])

    return data 

def read_SLICS_cosmo_params(file_path):
    """Reads the cosmological parameters that correspond to each ID, from the .dat file.

    Parameters
    ----------
    file_path : str
        Path to the .dat file containing the cosmological parameters

    Returns
    -------
    cosmo_params : dict
        Dictionary mapping each ID to its corresponding cosmological parameters

    """
    if not os.path.exists(file_path):
        raise IOError(f"SLICS cosmological parameters file {file_path} does not exist")

    dat = ascii.read(file_path, names=['ID', 'Om', 'h', 'w_0', 'sigma_8', 'Oc'])

    cosmo_params = {}
    
    for row in dat:
        id = int(row['ID'])
        params = {
            'Om': float(row['Om']),
            'h': float(row['h']),
            'w_0': float(row['w_0']),
            'sigma_8': float(row['sigma_8']),
            'Oc': float(row['Oc'])
        }
        cosmo_params[id] = params

    return cosmo_params

def map_cosmo_params_to_data(data, dat_file_path):
    """Maps cosmological parameters to data based on IDs.

    Parameters
    ----------
    data : numpy recarray
        The data containing IDs to be mapped to cosmological parameters.
    dat_file_path : str
        Path to the .dat file containing the cosmological parameters.

    Returns
    -------
    mapped_params : list of dict
        List of dictionaries containing cosmological parameters
        mapped to the data based on IDs.
    """
    # Read the cosmological parameters from the .dat file
    cosmo_params = read_SLICS_cosmo_params(dat_file_path)

    # Map the IDs in the recarray to the corresponding cosmological parameters
    mapped_params = []
    for row in data:
        id = row['id']
        params = cosmo_params.get(id)
        if params:
            mapped_params.append(params)
        else:
            print(f"No parameters found for ID {id}")

    return mapped_params