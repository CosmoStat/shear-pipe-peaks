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
