"""SLICS.

:Name: slics.py

:Description: This package contains methods to interface with
    the (Cosmo-)SLICS simulations (Hernois-Deraps et al).

:Authors: Martin Kilbinger <martin.kilbinger@cea.fr>
"""

import os
import numpy as np

from astropy.io import ascii
from astropy.table import vstack

from cs_util import cat as cs_cat


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
    "gamma1_sim",
    "gamma2_sim",
    "redshift_true_sim",
]

max_los = 4
max_zbin = 4
max_tile = 19


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
        raise IOError(f"SLICS catalogue file {file_path} does not exist")

    if all_col:
        include_names = None
    else:
        include_names = col_names_essential

    dat = ascii.read(file_path, names=col_names, include_names=include_names)

    return dat


def get_cat_name(cosmo_id, zbin, los, tile):
    """Get Cat Name.

    Return catalogue file name.

    Parameters
    ----------
    cosmo_id : str
        cosmology model ID
    zbin : int
        redshift bin number
    los : int
        line-of-sight number
    tile : int
        tile number (realisation)

    Returns
    -------
    str
        catalogue name

    """
    cat_name = (
        f"{cosmo_id}/LOS{los}/DES_MocksCat_{cosmo_id}_4_Bin{zbin}_LOS{los}_"
        + f"R{tile}.dat"
    )

    return cat_name


def read_multiple_catalogues(
    root_directory,
    cosmo_id="fid_f",
    zbins=None,
    lsos=None,
    tiles=None,
    combine="add",
    verbose=False,
):
    """Read Multiple Catalogues.

    Read and combine a number of SLICS catalogues.

    Parameters
    ----------
    root_directory : str
        root directory of input catalogues
    cosmo_id : str
        cosmology model ID
    zbins : list
        list of redshift bins to combine, default is None (combine all)
    lsos : list
        list of lines of sight numbers to combine, default is None
        (combine all)
    tiles : list
        list of tile numbers to combine, default is None (combine all)

    Returns
    -------
    dict
        catalogue content

    """
    if not os.path.exists(root_directory):
        raise IOError(f"Root directory {root_directory} does not exist")

    if not zbins:
        zbins = np.arange(max_zbin) + 1
    if not lsos:
        lsos = np.arange(max_los) + 1
    if not tiles:
        tiles = np.arrange(max_tile) + 1

    dat_comb = None

    # Loop over input catalogue identifiers
    for zbin in zbins:
        for los in lsos:
            for tile in tiles:

                # Get single catalogue
                cat_name = get_cat_name(cosmo_id, zbin, los, tile)
                cat_path = f"{root_directory}/{cat_name}"

                if verbose:
                    print(f"Reading catalogue {cat_name}...")
                dat = read_catalogue(cat_path, all_col=False)

                if dat_comb is None:
                    # First time: Copy catalogue
                    dat_comb = dat.copy()
                else:
                    # Other times: Combine with previous
                    if combine == "add":
                        dat_comb = vstack([dat_comb, dat])
                    else:
                        pass

    return dat_comb


def resample_z(dat, dndz_path, z_max=None):

    # Read external dndz file
    z_centers_ext, dndz_ext, z_edges_ext = cs_cat.read_dndz(dndz_path)

    # Create redshift histogram of SLICS catalogue, using external zbins
    dndz_slics, z_edges_slics = np.histogram(
        dat["redshift_true_sim"],
        bins=z_edges_ext
    )

    # Cut redshifts if desired
    if z_max is not None:
        idx_z_max = np.where(z_edges_ext < z_max)
        dndz_ext = dndz_ext[idx_z_max]
        dndz_slics = dndz_slics[idx_z_max]

    # Normalise redshift distributions
    dndz_ext = dndz_ext / sum(dndz_ext)
    dndz_slics = dndz_slics / sum(dndz_slics)

    probability = dndz_ext / dndz_slics

    

