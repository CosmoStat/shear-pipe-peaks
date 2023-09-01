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

# +
import sp_peaks           
print(f"sp_peaks version = {sp_peaks.__version__}")

from sp_peaks import slics
# -

cat_path = "/n17data/tersenov/SLICS/Cosmo_DES/06_f/LOS3/DES_MocksCat_06_f_4_Bin3_LOS3_R19.dat"

# Load catalogue, all columns
dat = slics.read_catalogue(cat_path)

# Print column names
print(dat.dtype.names)

# Print first line
print(dat[0])

# Load only essential columns
dat_ess = slics.read_catalogue(cat_path, all_col=False)

print(dat_ess[0])


