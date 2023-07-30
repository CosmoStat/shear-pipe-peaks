# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.7
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # CosmoSLICS simulation testing
# 07/2023                                                                       

# %matplotlib inline                                                            
# %load_ext autoreload                                                          
# %autoreload 2                                                                                                                                                 

import matplotlib.pylab as plt
import os
from cs_util import cat as cs_cat

# +
import sp_peaks           
print(f"sp_peaks version = {sp_peaks.__version__}")

from sp_peaks import slics
# -

root_directory = "/n17data/tersenov/SLICS/Cosmo_DES"

cat_path = f"{root_directory}/06_f/LOS3/DES_MocksCat_06_f_4_Bin3_LOS3_R19.dat"

# Load catalogue, all columns
dat = slics.read_catalogue(cat_path)

# DEBUG
dat_comb = dat

# Print column names
print(dat.dtype.names)

# Print first line
print(dat[0])

# Load only essential columns
dat_ess = slics.read_catalogue(cat_path, all_col=False)

from astropy.table import vstack
x = vstack([dat, dat])

len(dat)

print(dat_ess[0])

# Combine all four redshift bins for given cosmo ID, line of sight, and tile number
dat_comb = slics.read_multiple_catalogues(
    root_directory,
    cosmo_id="fid_f",
    zbins=None,
    lsos=[2],
    tiles=[5],
    combine="add",
    verbose=True,
)

len(dat_comb)

# +
# Read CFIS redshift distribution (blind version "A", ShapePipe)
dndz_CFIS_path = "/n17data/mkilbing/astro/data/CFIS/v1.0/nz/dndz_SP_A.txt"
res = slics.resample_z(dat_comb, dndz_CFIS_path, len(dat_comb) / 4, z_max=1.8)
# -


