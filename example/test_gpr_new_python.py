import argparse
import numpy as np
import emcee
import numpy.linalg as la
import matplotlib.pyplot as plt
# from sklearn.externals.joblib import Parallel, delayed, cpu_count
import joblib
from joblib import Parallel, delayed, cpu_count
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
test_gpr_new_python.py

EXAMPLE:
python test_gpr_new_python.py '/home/tersenov/shear-pipe-peaks' --nproc 5  
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
index_cosmo = 2
params_train = np.delete(params, index_cosmo, axis=0)
obs_train = np.delete(Peaks_Maps, index_cosmo, axis=0)

ncpu = cpu_count()
gp_scaling=np.array([Parallel(n_jobs = ncpu, verbose = 5)(delayed(gp_train)(index_bin, params_train = params_train, obs_train =  obs_train) for index_bin in range(Peaks_Maps.shape[2]))]).reshape(Peaks_Maps.shape[2], 2)

# Load the parameters and observables for the left-out cosmology (the first cosmology)
left_out_params = params[index_cosmo]  
left_out_obs = Peaks_Maps[index_cosmo]  

gp_list=gp_scaling[:,0]
scaling=gp_scaling[:,1]

test = GP_pred(left_out_params,gp_list,scaling)
peak_counts_pred = test[0]

# Calculate the absolute percentage error for each bin
abs_percentage_error = np.abs((Peaks_Maps_mean[0] - peak_counts_pred) / Peaks_Maps_mean[0]) * 100
# Calculate the mean MAPE across all bins
mean_mape = np.mean(abs_percentage_error)
# save the mean MAPE
np.savetxt(maindir+'/output/mean_mape.txt', [mean_mape])

plt.figure()
plt.plot(test[0] - Peaks_Maps_mean[0])
plt.show()
#save the plot
plt.savefig(maindir+'/output/test_gpr.pdf')
