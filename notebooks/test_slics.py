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
import copy
import numpy as np
from cs_util import cat as cs_cat

# +
import sp_peaks           
print(f"sp_peaks version = {sp_peaks.__version__}")

from sp_peaks import slics

# +
# Set input directories

# SLICS simulations
root_directory = "/n17data/tersenov/SLICS/Cosmo_DES"
#root_directory = "."

# CFIS redshift distribution
dndz_CFIS_directory = "/n17data/mkilbing/astro/data/CFIS/v1.0/nz"
# -

# ## Read SLICS catalogue

cat_path = f"{root_directory}/06_f/LOS3/DES_MocksCat_06_f_4_Bin3_LOS3_R19.dat"

# Load catalogue, all columns
dat = slics.read_catalogue(cat_path)

# Print column names
print(dat.dtype.names)

# Print first line
print(dat[0])

# Load only essential columns
dat_ess = slics.read_catalogue(cat_path, all_col=False)

print(dat_ess[0])

# ## Read and combine multiple SLICS catalogues

# Combine all four redshift bins for given cosmo ID, line of sight, and tile number
dat_comb = slics.read_multiple_catalogues(
    root_directory,
    cosmo_id="06_f",
    zbins=None,
    lsos=[2],
    tiles=[5],
    combine="add",
    verbose=True,
)

print(f"Number of galaxies = {len(dat_comb)}")
n_gal = slics.get_number_density(dat_comb)
print(f"Number density = {n_gal:.2f} arcmin^{{-2}}")

# ### Resample SLICS to match exteral dn/dz

# Set CFIS redshift distribution (blind version "A", ShapePipe)
dndz_CFIS_path = f"{dndz_CFIS_directory}/dndz_SP_A.txt"

# +
# Test resampling. Compare redshift distributions

# External (CFIS) redshift histogram
z_centers_ext, dndz_ext, z_edges_ext = cs_cat.read_dndz(dndz_CFIS_path)

# Original SLICS redshift histogram
dndz_slics, _ = np.histogram(dat_comb["redshift_true_sim"], bins=z_edges_ext)

# Original SLICS normalised redshift histogram
dndz_slics_norm, _ = np.histogram(dat_comb["redshift_true_sim"], bins=z_edges_ext, density=True)
# -

# Set the number of objects to resample.
# Has to be smaller than number of input objects.
n_goal = len(dat_comb) / 2

# Resample
slics.resample_z(dat_comb, dndz_CFIS_path, n_goal, z_max=1.8, verbose=True)

# +
# Testing

# Resampled SLICS redshift histogram
dndz_resampled, _ = np.histogram(dat_comb["redshift_true_sim"], bins=z_edges_ext)

# +
# Testing: resampled numbers are never higher than original ones
fig, ax = plt.subplots(figsize=(8, 8))

ax.plot(
    z_centers_ext, dndz_slics, '-',
    z_centers_ext, dndz_resampled, '-.',
)
ax.set_xlim([0, 2])
plt.savefig("dndz_slics_res.pdf")

# +
# Testing: ratio of resampled to original numbers are never larger than unity

fig, ax = plt.subplots(figsize=(8, 8))

ax.plot(
    z_centers_ext, dndz_resampled / dndz_slics, '-',
)
ax.set_xlim([0, 2])
plt.savefig("dndz_slics_res_ratio.pdf")

# +
# Resampled SLICS redshift histogram
dndz_resampled_norm, _ = np.histogram(dat_comb["redshift_true_sim"], bins=z_edges_ext, density=True)

# Testing: resampled dndz follows external dndz

fig, ax = plt.subplots(figsize=(8, 8))

ax.plot(z_centers_ext, dndz_ext, '-', label='CFIS')
ax.plot(z_centers_ext, dndz_slics_norm, '-', label='SLICS')
ax.plot(z_centers_ext, dndz_resampled_norm, '-', label='resampled')
ax.set_xlim([0, 2])
_ = ax.legend()
plt.savefig("dndz_CFIS_slics_res.pdf")
# -


