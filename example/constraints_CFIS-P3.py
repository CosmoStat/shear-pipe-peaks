import argparse
import numpy as np
import emcee
import numpy.linalg as la
import matplotlib.pyplot as plt
from sklearn.externals.joblib import Parallel, delayed, cpu_count
from multiprocessing import cpu_count, Pool
import time
from chainconsumer import ChainConsumer
from utils import *
from likelihood import *
import os
import sys

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""                                                                
 
SCRIPT NAME:                                                                   
constraints_CFIS-P3.py

EXAMPLE:
python constraints_CFIS-P3.py '/home/baumont/software/shear-pipe-peaks' --nproc 4  
"""
)
parser.add_argument('maindir', help='work directory containing input and output directories')
parser.add_argument('--nproc', help='number of processes', type=int, default=1)

args = parser.parse_args()

#for tex in ChainConsumer                                     
pref = os.environ['CONDA_PREFIX']
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'

#Specify the free parameters for redshift, shear calibration, baryonic correction and cut in S/N                                                    
# Check the README to know the possibility for the different cases and reproduce the plots                                                          
param_z = '065'
param_z_cov = '0.65'
param_cal = 'dm_1deg'
param_baryonic_correction = 'Fid'
param_cut = 19
n_patches = 13
maindir = args.maindir
nproc = args.nproc
#Load the file names for peaks file                                                                                                                 
file_name_peaks_Maps = (np.array(np.loadtxt(maindir+'/input/list_cosmo_peaks_z{}.txt'.format(param_z), usecols=0,dtype=str)))

#Extract cosmological parameters values from the simulations associated to eachpeak counts file                                                    
params = np.array([takes_params(f) for f in file_name_peaks_Maps])

#identifies fiducial cosmology index
index_fiducial = 15
params_fiducial = params[index_fiducial]

#load the peaks distributions for the theoretical prediction                                                                                        
Peaks_Maps_DM = np.array([np.load(maindir+'/input/peaks_z{}/%s'.format(param_z)%(fn), mmap_mode='r') for fn in file_name_peaks_Maps])[:,:,:param_cut]

# load baryon corrections                                 
bar_corr = np.load(maindir+'/input/{}_correction.npy'.format(param_baryonic_correction))[:param_cut]

#Apply the baryonic correction                                                                                                                      
Peaks_Maps = np.copy(Peaks_Maps_DM)

if(param_baryonic_correction == 'no'):
    # no baryonic correction
    for i in range(Peaks_Maps_DM.shape[0]):
        for j in range(Peaks_Maps_DM.shape[1]):
            Peaks_Maps[i,j,:] = Peaks_Maps_DM[i,j,:] * 1
else:
    # apply the choosen baryonic correction
    for i in range(Peaks_Maps_DM.shape[0]):
        for j in range(Peaks_Maps_DM.shape[1]):
            Peaks_Maps[i,j,:] = Peaks_Maps_DM[i,j,:] * bar_corr

# take the mean over realizations                            
Peaks_Maps_mean=np.mean(Peaks_Maps, axis=1)

# create array for snr histogram                                    
snr_array=np.linspace(-2,6,31)
snr_centers=0.5*(snr_array[1:]+snr_array[:-1])

#Load peaks distribution for data                                     
peaks_data = np.load(maindir+'/input/peaks_mean_{}.npy'.format(param_cal), mmap_mode='r')[:param_cut]

#Train Gaussian Processes regressor                                    
#this returns a tuple consisting of (list of GP, scaling)                                                                                           
ncpu = cpu_count()
gp_scaling=np.array([Parallel(n_jobs = ncpu, verbose = 5)(delayed(gp_train)(index_bin, params_train = params, obs_train =  Peaks_Maps) for index_bin in range(Peaks_Maps.shape[2]))]).reshape(Peaks_Maps.shape[2], 2)

gp_list=gp_scaling[:,0]
scaling=gp_scaling[:,1]

#Covariance matrix                                  
#Load peaks to compute covariance matrix
cov_peaks_DM = np.load(maindir+'/input/convergence_gal_mnv0.00000_om0.30000_As2.1000_peaks_2arcmin_{}_b030_snr_min_max_ngal_7.npy'.format(param_z_cov),\
 mmap_mode='r')[:,:param_cut]

#Apply Baryonic Correction                                                                                                                          
cov_peaks = np.copy(cov_peaks_DM)
for i in range(Peaks_Maps_DM.shape[1]):
    cov_peaks[i,:] = cov_peaks_DM[i,:] * bar_corr

#Compute covariance matrix and scale it for sky coverage                                                                                            
cov=(1/n_patches)*np.cov(cov_peaks.T)

#compute the inverse of the covariance                                             
icov=la.inv(cov)

# get constraints with MCMC                                                                                                                         
# compute hartlap factor              
n_real=Peaks_Maps_DM.shape[1]
n_bins=len(peaks_data)
norm=(n_real-n_bins-2)/(n_real-1)

#Define priors (range of sims)              
M_nu_min = 0.06  # minimum from oscillation experiments               
M_nu_max = 0.62
Omega_m_min = 0.18
Omega_m_max = 0.42
A_s_min = 1.29
A_s_max = 2.91

#Specify number of dimensions for parameter space, number of walkers, initial position                                                              
ndim, nwalkers = 3,250
pos = [params_fiducial +  1e-3*np.random.randn(ndim) for i in range(nwalkers)]

print("{0} CPUs".format(ncpu))

with Pool(processes=nproc) as pool:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnpost, pool=pool, args=[peaks_data, icov,gp_list, scaling, norm,M_nu_min, M_nu_max, Omega_m_min, Omega_m_max, A_s_min, A_s_max])
    start = time.time()
    sampler.run_mcmc(pos, 6500, progress=True)
    end = time.time()
    multi_time = end - start
    print("Multiprocessing took {0:.1f} seconds".format(multi_time))
# remove burn-in                                                        
samples = sampler.chain[:,200:, :].reshape((-1, ndim))
#save results
np.save(maindir+'/output/constraints_z{}_{}_{}corr_{}snr.npy'.format(param_z,param_cal,param_baryonic_correction,param_cut),samples)
