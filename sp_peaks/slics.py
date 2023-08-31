"""SLICS.

:Name: slics.py

:Description: This package contains methods to interface with
    the (Cosmo-)SLICS simulations (Hernois-Deraps et al).

:Authors: Martin Kilbinger <martin.kilbinger@cea.fr>
"""

import os
import random
import numpy as np

from astropy.io import ascii
from astropy.table import vstack

from lenspack.geometry import measures

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


def drop_marked(dat):
    """Drop Marked.

    Remove rows with marked (np.inf) redshifts.

    Parameters
    ----------
    dat : dict
        input SLICS catalogue

    """
    idx_drop = np.isnan(dat["redshift_true_sim"])
    dat_out = dat[~idx_drop]


def resample_z(dat, dndz_path, n_goal, z_max=None, verbose=False):
    """Resample Z.

    Resample catalogue according to external redshift distribution from file.

    Parameters
    ----------
    dat : dict
        input SLICS catalogue
    dndz_path : str
        path to external redshift distribution file
    n_goal : int
        number of galaxies to be reached after resampling
    z_max : float, optional
        maximum redshift, default is ``None`` (use all redshifts)
    verbose : bool, optional
        verbose output if ``True``, default is ``False``
    """
    # Read external dndz file
    z_centers_ext, dndz_ext, z_edges_ext = cs_cat.read_dndz(dndz_path)

    # Create redshift histogram of SLICS catalogue, using external zbins
    dndz_slics, _ = np.histogram(
        dat["redshift_true_sim"],
        bins=z_edges_ext
    )

    # Cut redshifts if desired
    if z_max is not None:
        idx_z_max = np.where(z_edges_ext < z_max)
        z_centers_ext = z_centers_ext[idx_z_max]
        dndz_ext = dndz_ext[idx_z_max]
        dndz_slics = dndz_slics[idx_z_max]

    # Normalise external redshift distributions to desired number of resampled
    # objects
    dndz_ext = dndz_ext / sum(dndz_ext) * n_goal

    # Ratio of external to SLICS histogram = fraction of galaxies in each
    # bin to resample
    ratio = dndz_ext
    w_nonz = dndz_slics > 0

    # Set ratio where SLICS histogram is non-zero
    ratio[w_nonz] = ratio[w_nonz] / dndz_slics[w_nonz]

    # Set to 0 where SLICS histogram is zero -> resample zero galaxies
    # in this bin
    ratio[~w_nonz] = 0

    idx_zero = np.where(~w_nonz)[0]
    if len(idx_zero) > 0:
        print(
            f"Warning: in {len(idx_zero)} z-bins the number of resampled galaxies"
            + " will be set to zero."
        )

    idx_over = np.where(ratio > 1)[0]
    if len(idx_over) > 0:
        print(
            f"Warning: in {len(idx_over)} z-bins the number of resampled"
            + " galaxies will be set to the original number."
        )

        # Truncate ratio to one
        ratio[idx_over] = 1

    # Get indices in external redshift histogram of SLICS input catalogue
    idx_z = np.digitize(dat["redshift_true_sim"], z_edges_ext)

    n_tot = 0
    for idx in range(len(dndz_slics)):

        if ratio[idx] == 1:
            # No galaxy is removed in this z-bin
            continue

        # Get index list for this z-bin                                 
        w = (np.where(idx_z == idx + 1))[0]

        if ratio[idx] > 0:
            # Number of objects to remove
            n_drop = dndz_slics[idx] - int(ratio[idx] * dndz_slics[idx])

            # Create sample of Indices of objects to be dropped for this bin
            i_drop = np.array(random.sample(list(w), n_drop))

            # Mark objects to be dropped with invalid redshift              
            dat["redshift_true_sim"][i_drop] = np.inf

        else:
            # Remove all objects in this z-bin
            n_drop = len(w)
            dat["redshift_true_sim"][w] = np.inf

        n_tot = n_tot + n_drop

    if verbose:
        print(f'dropping {n_tot} to match nofz'.format(n_tot))

    drop_marked(dat)


def get_extend(dat):
    """Get extend.

    Return extend of catalogue.

    Parameters
    ----------
    dat : dict
        input SLICS catalogue

    Returns
    -------
    list
        extend ([ra_min, ra_max, dec_min, dec_max])

    """
    ra_min = np.amin(dat["RA"])
    ra_max = np.amax(dat["RA"])
    dec_min = np.amin(dat["Dec"])
    dec_max = np.amax(dat["Dec"])

    return [ra_min, ra_max, dec_min, dec_max]


def get_number_density(dat):
    """Get Number Density

    Return galaxy number density, not accounting for masking

    Parameters
    ----------
    dat : dict
        input SLICS catalogue

    Returns
    -------
    float
        number of galaxies [arcmin^{-2}]

    """
    extend = get_extend(dat)

    solid_angle_deg = measures.solid_angle(extend)
    solid_angle_amin = solid_angle_deg * 60 ** 2

    ngal = len(dat) / solid_angle_amin

    return ngal
