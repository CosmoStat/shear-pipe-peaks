## The impact of systematic errors on weak-lensing peak counts for UNIONS 

This repository containts the results of the paper: UNIONS: The impact of systematic errors on weak-lensing peak counts

Authors: Emma Ayçoberry, Virginia Ajani, Axel Guinot, Martin Kilbinger, Valeria Pettorino, Samuel Farrens, Jean-Luc Starck, Raphaël Gavazzi


The Ultraviolet Near-Infrared Optical Northern Survey (UNIONS) is an ongoing deep photometric multi-band survey of the Northern sky. As part of UNIONS, the Canada-France Imaging Survey (CFIS) provides r-band data with a median seeing of 0.65 arcsec, which we use to study weak-lensing peak counts for cosmological inference.
This work aims to assess systematics effects for weak-lensing peak counts and their impact on cosmological parameters for the UNIONS survey. In particular, we present results on local calibration, metacalibration shear bias, baryonic feedback, the source galaxy redshift estimate, intrinsic alignment, and the cluster member dilution. We expect constraints to become more reliable with future (larger) data catalogues, for which the current pipeline will provide a starting point.


### Content

1. [Dependencies](#dependencies)
2. [Description of the input files](#description-of-the-input-files)
3. [Example](#example)


#### Dependencies

The following python packages should be installed with their specific dependencies to be able to run the example:

- [numpy](https://numpy.org/install/)
- [astropy](https://www.astropy.org)
- [lenspack](https://github.com/CosmoStat/lenspack.git)
- [scipy](https://scipy.org/install/)
- [scikit-learn](https://scikit-learn.org/stable/install.html)>=0.20
- [matplotlib](https://matplotlib.org/stable/users/installing/index.html)
- [emcee](https://emcee.readthedocs.io/en/v2.2.1/user/install/)==3.1.1
- [chainconsumer](https://samreay.github.io/ChainConsumer/)
- [getdist](https://getdist.readthedocs.io/en/latest/intro.html#getting-started)
- [zenodo_get](https://github.com/dvolgyes/zenodo_get) 
- [jupyter](https://jupyter.org/install)






#### Description of the input files


- List of the peaks from the simulations at the 2 redshifts: 
list_cosmo_peaks_z065.txt and list_cosmo_peaks_z068.txt

- Peaks distribution coming from simulations: 
these need to be downloaded from [here](https://zenodo.org/record/6344515#.Yk2j6S0QOqA), [DOI: 10.5281/zenodo.6344515](https://zenodo.org/record/6344515#.Yk2k3C0QOqA). This step is already present in the example notebook where you create a folder called peaks_z065/ in the repo where the peaks distribution coming from simulations have to be stored and download the peaks as:

  `zenodo_get 10.5281/zenodo.6344515`

  in the folder `input/peaks_z065`

- Baryonic correction:
Fid_correction.npy, HighAGN_correction.npy and LowAGN_correction.npy

- Peaks from CFIS-P3 data with different parameters:
peaks_mean_global.npy = peaks obtain with the global calibration, peaks_mean_Xdeg.npy = peaks obtain with the calibration on X square degree, peaks_mean_dm_1deg.npy = peaks obtain with the calibration on 1 square degree, and the multiplicative shear bias delta m = 0.007

- Peaks for the covariance (at the 2 redshift):
convergence_gal_mnv0.00000_om0.30000_As2.1000_peaks_2arcmin_0.65_b030_snr_min_max_ngal_7.npy
convergence_gal_mnv0.00000_om0.30000_As2.1000_peaks_2arcmin_0.68_b030_snr_min_max_ngal_7.npy

##### **Description of the parameters that can be used in the notebook**
- param_z: the redshift of the wanted simulations
Can be ’065’ or ‘068’

- param_z_cov: the redshift of the wanted simulations, different of the previous file to respect the name given by the simulations
Can be ‘0.65’ or ‘0.68’

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

- Fig. 9: test the redshift
param_z = '065' or ’068’ for the different cases
param_z_cov = '0.65' or ‘0.68’  for the different cases
param_cal = ‘global’ 
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

To run the example you need to git clone this repository as 

`git clone https://github.com/CosmoStat/shear-pipe-peaks.git`

ad then go to the example folder and launch the jupyter notebook [constraints_CFIS-P3.ipynb](https://github.com/CosmoStat/shear-pipe-peaks/blob/main/example/constraints_CFIS-P3.ipynb) by doing:

`cd example`

`jupyter notebook`

and run it following the instructions in the description of the cells.

Follow the [Configuration needed to obtain the different figures of the paper](#configuration-needed-to-obtain-the-different-figures-of-the-paper) to get the results of the paper or play with the different [input parameters](#description-of-the-parameters-that-can-be-used-in-the-notebook) to get different configurations. 

Contact vajani@phys.ethz.ch or emmaaycoberry@gmail.com for questions concerning the repo.
