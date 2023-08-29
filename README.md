# Weak-lensing peak counts with ShapePipe processed data

This python package contains methods for cosmological parameter inference from weak-lensing peak counts.
It works with galaxy catalogue data processed with ShapePipe, but is general to be used with other data.

To use this code, download this repository by clicking on the `Clone` button above.

The easiest way to install the required software is using `conda`.
You can follow the instructions to install Anaconda [here](https://docs.anaconda.com/anaconda/install/index.html) or miniconda [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

Included is code to reproduce results and plots from Ayçoberry et al. (2022); see below.

## Content

1. [The python library](#the-python-library)
   1. [Installation](#installation)
   1. [Usage](#usage)
   1. [Examples](#examples)
1. [The impact of systematic errors on weak-lensing peak counts for UNIONS](#the-impact-of-systematic-errors-on-weak-lensing-peak-counts-for-unions)
   1. [Installation](#installation)
   1. [Description of the input files](#description-of-the-input-files)
   1. [Example](#example)

## The python library

### Installation

To create the conda environment type
```bash
conda create -f environment.yml
```

### Usage

Once installed, any library files can be used via `import` in a python script or notebook, e.g.:
```python
import sp_peaks
from sp_peaks import slics
```

### Examples

Example notebooks and associated python scripts can be run in the `notebooks` folder.


## The impact of systematic errors on weak-lensing peak counts for UNIONS 

This part reproduces the results of the paper "UNIONS: The impact of systematic errors on weak-lensing peak counts".

Authors: Emma Ayçoberry, Virginia Ajani, Axel Guinot, Martin Kilbinger, Valeria Pettorino, Samuel Farrens, Jean-Luc Starck, Raphaël Gavazzi, Michael J. Hudson

Abstract:
The Ultraviolet Near-Infrared Optical Northern Survey (UNIONS) is an ongoing deep photometric multi-band survey of the Northern sky. As part of UNIONS, the Canada-France Imaging Survey (CFIS) provides r-band data with a median seeing of 0.65 arcsec, which we use to study weak-lensing peak counts for cosmological inference.
This work aims to assess systematics effects for weak-lensing peak counts and their impact on cosmological parameters for the UNIONS survey. In particular, we present results on local calibration, metacalibration shear bias, baryonic feedback, the source galaxy redshift estimate, intrinsic alignment, and cluster member dilution. We expect constraints to become more reliable with future (larger) data catalogues, for which the current pipeline will provide a starting point. This paper investigates for the first time with UNIONS weak-lensing data and peak counts the impact of the local calibration, residual multiplicative shear bias, redshift uncertainty, baryonic feedback, intrinsic alignment and cluster member dilution. The value of matter density parameter is the most impacted and can shift up to ~ 0.03 which corresponds to 0.5 sigma depending on the choices for each systematics. We expect constraints to become more reliable with future (larger) data catalogues, for which the current code provides a starting point.



### Installation

To create the conda environment type
```bash
conda create -f environment_3.6.yml
```
Note that python version 3.6.13 is required.


#### Description of the input files

**Note about redshifts:** as we found in the paper that the change in the peaks ditribution (and on the corresponding constraints) between redshift z=0.65 and z=0.68 is not significant, we provide here only the arrays at z=0.65. If interested in the distributions at z=0.68 please contact emma.aycoberry@iap.fr or vajani@phys.ethz.ch. 


- List of the peaks from the simulations at redshift z=0.65:
list_cosmo_peaks_z065.txt 

- Peaks distribution coming from simulations: 
these need to be downloaded from [here](https://zenodo.org/record/6344515#.Yk2j6S0QOqA), [DOI: 10.5281/zenodo.6344515](https://zenodo.org/record/6344515#.Yk2k3C0QOqA). This step is already present in the example notebook where you create a folder called peaks_z065/ in the repo where the peaks distribution coming from simulations have to be stored and download the peaks as:

  `zenodo_get 10.5281/zenodo.6344515`

  in the folder `input/peaks_z065`

- Baryonic correction:
Fid_correction.npy, HighAGN_correction.npy and LowAGN_correction.npy

- Peaks from CFIS-P3 data with different parameters:
peaks_mean_global.npy = peaks obtain with the global calibration, peaks_mean_Xdeg.npy = peaks obtain with the calibration on X square degree, peaks_mean_dm_1deg.npy = peaks obtain with the calibration on 1 square degree, and the multiplicative shear bias delta m = 0.007

- Peaks for the covariance:
convergence_gal_mnv0.00000_om0.30000_As2.1000_peaks_2arcmin_0.65_b030_snr_min_max_ngal_7.npy

##### **Description of the parameters that can be used in the notebook**
- param_z: the redshift of the desired simulations
Default is ’065’ (if interested in the distributions at other redshifts, please contact emma.aycoberry@iap.fr or vajani@phys.ethz.ch.)

- param_z_cov: the redshift of the desired simulations, different of the previous file to respect the name given by the simulations
Default is ’0.65’

- param_cal: the way the calibration is done
Can be ‘global’ (global calibration), ’05’ (calibration on 0.5 square degree), 1’, ‘2’, ‘4’, ‘dm_1’ (calibration on 1 square degree, with the multiplicative shear bias delta m = 0.007)

- param_baryonic_correction: which baryonic correction we want to apply
Can be ‘no’ if we don’t want to correct the simulations, or ‘Fid’, ‘HighAGN’, ‘LowAGN’ to correct the simulations

- param_cut: up to which index of SNR we are going
Can be 19 if we want to cut the last bins or 30 if we want to keep all the bins.


______________________________________
##### _Configuration needed to obtain the different figures of the paper_
- Fig. 7: test the calibration size
param_z = '065' 
param_z_cov = '0.65'
param_cal = ‘global’ or ’05’ or ‘1’ or ‘2’ or ‘4’ for the different cases
param_baryonic_correction = 0
param_cut = 30

- Fig. 8: test the residual shear bias
param_z = '065' 
param_z_cov = '0.65'
param_cal = ‘global’ or ‘1’ or ‘dm_1’ for the different cases
param_baryonic_correction = 0
param_cut = 30


- Fig. 11: test the baryonic correction
param_z = '065' 
param_z_cov = '0.65'
param_cal = ‘global’ 
param_baryonic_correction = 0 or ‘LowAGN’ or ‘Fid’ or ‘HighAGN’ for the different cases
param_cut = 30

- Fig. 12: test the intrinsic alignment and boost factor
param_z = '065' 
param_z_cov = '0.65' 
param_cal = ‘global’ 
param_baryonic_correction = 0
param_cut = 30 or 19 for the different cases

- Fig. 13: ideal model vs conservative model
Ideal model:
param_z = '065' 
param_z_cov = '0.65' 
param_cal = ‘global’ 
param_baryonic_correction = 0
param_cut = 30

*Conservative model:*
param_z = '065',

param_z_cov = '0.65', 

param_cal = ‘dm_1’,  

param_baryonic_correction = 'Fid', 

param_cut = 19, 

#### Example

After installing the package, see [Dependencies](#dependencies), go to the example folder and launch the jupyter notebook [constraints_CFIS-P3.ipynb](https://github.com/CosmoStat/shear-pipe-peaks/blob/main/example/constraints_CFIS-P3.ipynb) by typing:

`cd example`

`jupyter notebook`

and run it following the instructions in the description of the cells.

Follow the [Configuration needed to obtain the different figures of the paper](#configuration-needed-to-obtain-the-different-figures-of-the-paper) to get the results of the paper or play with the different [input parameters](#description-of-the-parameters-that-can-be-used-in-the-notebook) to get different configurations. 

Contact vajani@phys.ethz.ch or emma.aycoberry@iap.fr for questions concerning the repo.
