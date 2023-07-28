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
z_centers_ext, dndz_ext, z_edges_ext = cs_cat.read_dndz(dndz_CFIS_path)

nz = len(z_edges_ext)
print(len(z_centers_ext), len(dndz_ext), len(z_edges_ext))
# -

dndz_slics, z_edges_slics = np.histogram(dat_comb["redshift_true_sim"], bins=z_edges_CFIS)

# CHeck that z
max(z_edges_slics - z_edges_ext)

plt.plot(res[1][:-1], res[0] / sum(res[0]), '-')
plt.plot(z_centers_ext, dndz_ext / sum(dndz_ext), '-')

res[1][np.where(res[0] == 0)]

z_max = 1.8

idx_z_max = np.where(z_edges_ext < z_max)                               
dndz_ext = dndz_ext[idx_z_max]                                          
dndz_slics = dndz_slics[idx_z_max] 

probability = dndz_ext / dndz_slics 

probability


