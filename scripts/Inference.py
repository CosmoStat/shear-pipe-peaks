import os
import numpy as np
import matplotlib.pyplot as plt
import emcee
import numpy.linalg as la
from getdist import plots, MCSamples, parampriors
from joblib import Parallel, delayed, cpu_count
from multiprocessing import cpu_count, Pool
import time
from chainconsumer import ChainConsumer
from utils import *
from scipy.ndimage.filters import gaussian_filter as gf
import scipy.stats as stats
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C
        
## Creating arrays for cosmological parameters and peak counts data vectors
# Read the CosmoTable.dat and parse into a dictionary
cosmo_table_path = '/home/tersenov/shear-pipe-peaks/example/CosmoTable.dat'
cosmo_params = {}

with open(cosmo_table_path, 'r') as file:
    lines = file.readlines()[1:]  # Skip the header

for line in lines:
    parts = line.split()
    cosmo_id = int(parts[0])
    parameters = [float(param) for param in parts[1:]]
    cosmo_params[cosmo_id] = parameters

# Initialize the params array based on the cosmo_params dictionary
params = np.array([cosmo_params[cosmo_id] for cosmo_id in range(25)])


n_cosmologies = 25
n_realizations = 10  # Each file has 10 realizations (2 seeds x 5 LOS)
n_peaks = 30  # Number of peak counts bins per realization

# Initialize the Peaks_Maps array with zeros
Peaks_Maps = np.zeros((n_cosmologies, n_realizations, n_peaks))

# Assume we're working with Bin1 
bin_number = 1

for cosmo_id in range(n_cosmologies):
    file_name = f"{cosmo_id:02d}_Bin{bin_number}_ks_run666.npy"
    file_path = f"/n17data/tersenov/SLICS/Cosmo_DES/summary_stats/{file_name}"
    data = np.load(file_path, allow_pickle=True)
    for realization_id in range(n_realizations):
        peak_counts = data[realization_id][2]
        Peaks_Maps[cosmo_id, realization_id, :] = peak_counts
        
param_cut = 19
Peaks_Maps = Peaks_Maps[:,:,:param_cut]
Peaks_Maps.shape

#for tex in ChainConsumer
pref = os.environ['CONDA_PREFIX']
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'

initial_length_scales = [0.5 * (0.5482 - 0.1019), # Om
                         0.5 * (0.8129 - 0.6034), # h
                         0.5 * (-0.5223 + 1.9866), # w_0 (note the sign inversion for positive value)
                         0.5 * (1.3428 - 0.4716), # sigma_8
                         0.5 * (0.5009 - 0.0546)] # Oc

# Adjust bounds for each parameter's length scale
length_scale_bounds = [(0.1 * ls, 10 * ls) for ls in initial_length_scales]  # Example bounds


def gp_train_new(index_bin, params_train, obs_train):
    
    """
    Training of Gaussian processes per bin.
    It takes as input the bin indexes, training parameters and observables.
    It gives as output a list of GP and scaling values to avoid over-regularization.
    
    Parameters
    ----------
    
    index_bin: int
               bin index of the multipole for PS or S/N for Peaks
    params_train: numpy.ndarray of shape (#cosmo, #parameters)
                  training parameters from simulations
    obs_train: numpy.ndarray of shape (#cosmo, #realisations, #bins)
               training observables from simulations (e.g. the summary statistics)
   
    Returns
    -------
    gp: object
     GP for a given bin
    scaling: float
         scaling for the data values to avoid over-regularization
    """
    #means over the realisations for each bin looping over each cosmology
    obs_means=np.array([np.mean(obs_temp,axis=0)[index_bin] for obs_temp in obs_train])
    
    #sem over the realisations for each bin looping over each cosmology
    alpha_diag=np.array([stats.sem(obs_temp,axis=0)[index_bin] for obs_temp in obs_train])
    
    #data scaling factor
    scaling = (np.mean(obs_means))
    print(scaling)
    if scaling == 0 or np.isnan(scaling):
        raise ValueError("Scaling factor is zero or NaN.")
    
    #scale the data values to avoid over-regularization
    obs_means /= scaling
    
    #scale the standard errors of the mean values to avoid over-regularization
    alpha_diag /= scaling
    
    #define the kernel
    # kernel = C(5.0, (1e-4, 1e4)) * RBF([2, 0.3, 5, 4, 2], (1e-4, 1e4))
    kernel = C(5.0, (1e-4, 1e4)) * RBF([0.22315, 0.10475, 0.73215, 0.4356, 0.22315],[(0.022315, 2.2315), (0.010475, 1.0475), (0.073215, 7.3215), (0.04356, 4.356), (0.022315, 2.2315)])
    
    #instantiate the gaussian process
    gp = GaussianProcessRegressor(kernel=kernel,alpha=(alpha_diag)**2,n_restarts_optimizer=50,normalize_y=True)
    gp.fit(params_train,obs_means)
    
    return gp, scaling
    
def GP_pred(params_pred, gp_list, scaling):
    
    """
    
    Computes the prediction for a new point in parameter space.
    It is called by the likelihoood to get thee theoretical predictions
    for the observables.
    
    Parameters
    ----------
    
    params_pred: list of floats
                  new point in parameter space to get the prediction for
               
    gp_list: list of objects
            list of GP for a given bin
    scaling: numpy.ndarray of shape (#cosmo, #realisations, #bins)
               training observables from simulations (e.g. PS or Peaks)
   
    Returns
    -------
    gp: list of objects
       list of GP per bin
    scaling: numpy.ndarray
         scaling for the data values to avoid over-regularization
    """
     
    pred_list=[]
    sigma_list=[]
    
    for gp in gp_list:
        
        pred,sigma=gp.predict(np.array(params_pred).reshape(1,len(params_pred)),return_std=True)
        
        pred_list.append(pred[0])
        
        sigma_list.append(sigma[0])
    
    return np.array(pred_list * scaling),\
                np.array(sigma_list * scaling)    
    
#this returns a tuple consisting of (list of GP, scaling)
ncpu = cpu_count()
gp_scaling=np.array([Parallel(n_jobs = ncpu, verbose = 5)(delayed(gp_train_new)(index_bin, params_train = params, obs_train = Peaks_Maps) for index_bin in range(Peaks_Maps.shape[2]))]).reshape(Peaks_Maps.shape[2], 2)

gp_list=gp_scaling[:,0]
scaling=gp_scaling[:,1]

params_fiducial = np.array([0.2905, 0.6898, -1., 0.8364,0.2432])
test_fid = GP_pred(params_fiducial,gp_list,scaling)


bins = np.linspace(-2, 6, 31)
binscenter = 0.5 * (bins[:-1] + bins[1:])
plt.plot(binscenter[:19], test_fid[0])
plt.show()
plt.savefig('/home/tersenov/shear-pipe-peaks/output/TEST_GP_HIGH_NOISE.pdf')

### Get constraints with MCMC
cov_SS_PC = np.load('/home/tersenov/shear-pipe-peaks/output/covariance_matrix_SS_PC_cs.npy', allow_pickle=True)[:param_cut,:param_cut]

n_patches = 19
cov = (1/n_patches)*cov_SS_PC
icov = la.inv(cov)


def lnlike_new(theta, data, icov, gp_list, scaling, norm):
    """
    Likelihood function (loglike of a Gaussian)
    
    Parameters
    ----------
    theta: list
           point in parameter space
    data: numpy.ndarray 
          observed data
    icov: numpy.ndarray of shape (#bins, #bins)
          inverse of the covariance matrix
    gp_list: list of objects
             list of GP per bin
    scaling: numpy.ndarray
             scaling for the data values to avoid over-regularization
    norm: float
          the Hartlap factor
    Returns
    -------
    lnlikelihood: float
               value for the likelihood function
    """
    Om, h, w_0, sigma_8, Oc = theta
    obs_pred = GP_pred(theta,gp_list,scaling)[0] 
    lnlikelihood = -0.5 * norm * np.dot(np.dot((data.reshape(len(data),1)-obs_pred.reshape(len(data),1)).T, icov),(data.reshape(len(data),1)-obs_pred.reshape(len(data),1)))

    return lnlikelihood

def lnprior_new(theta, Om_min, Om_max, h_min, h_max, w_0_min, w_0_max, sigma_8_min, sigma_8_max, Oc_min, Oc_max):
    """
    Prior (log, flat). If the parameter values are within the minimum and maximum values
    returns 0.0. Otherwise returns -inf.
    
    Parameters
    ----------
    
    theta: list
           point in parameter space
   
    Returns
    -------
    0.0 if the parameter values are within the minimum and maximum values
    otherwise -np.inf
    """
    Om, h, w_0, sigma_8, Oc = theta
    
    if (Om_min < Om < Om_max and
        h_min < h < h_max and
        w_0_min < w_0 < w_0_max and
        sigma_8_min < sigma_8 < sigma_8_max and
        Oc_min < Oc < Oc_max):
        return 0.0
    
    return -np.inf

def lnpost_new(theta, data, icov, gp_list, scaling, norm, Om_min, Om_max, h_min, h_max, w_0_min, w_0_max, sigma_8_min, sigma_8_max, Oc_min, Oc_max):
    """
    Posterior (log). It is computed as lnprior + lnlike.
    Parameters
    ----------
    
    theta: list
           point in parameter space
    data: numpy.ndarray 
          observed data
    icov: numpy.ndarray of shape (#bins, #bins)
          inverse of the covariance matrix
    gp_list: list of objects
             list of GP per bin
    scaling: numpy.ndarray
             scaling for the data values to avoid over-regularization
    norm: float
          the Hartlap factor
    Returns
    -------
    lp+lnlike: float
               value for the posterior distribution
    """
    lp = lnprior_new(theta, Om_min, Om_max, h_min, h_max, w_0_min, w_0_max, sigma_8_min, sigma_8_max, Oc_min, Oc_max)
    if not np.isfinite(lp):
        return -np.inf
    
    return lp+lnlike_new(theta, data, icov, gp_list, scaling, norm)

#### Compute the Hartlap factor
peaks_data = test_fid[0]
n_real = Peaks_Maps.shape[1]
n_bins = len(peaks_data)
norm = (n_real-n_bins-2)/(n_real-1)

norm=1

#### Define values for the prior for the parameters
Om_min = 0.1019
Om_max = 0.5482
h_min = 0.6034
h_max = 0.8129
w_0_min = -1.9866
w_0_max = -0.5223
sigma_8_min = 0.4716
sigma_8_max = 1.3428
Oc_min = 0.0546
Oc_max = 0.5009

#### Specify number of dimensions for parameter space, number of walkers, initial position
ndim, nwalkers = 5, 250
test_init_params = np.array([ 0.33,  0.65, -1.5 ,  1.1,  0.16])
pos = [test_init_params +  1e-3*np.random.randn(ndim) for i in range(nwalkers)]

#### Run MCMC
print("{0} CPUs".format(ncpu))

with Pool(processes=ncpu//2) as pool:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnpost_new, pool=pool, args=[peaks_data, icov, gp_list, scaling, norm, Om_min, Om_max, h_min, h_max, w_0_min, w_0_max, sigma_8_min, sigma_8_max, Oc_min, Oc_max])
    start = time.time()
    sampler.run_mcmc(pos, 12000, progress=True)
    end = time.time()
    multi_time = end - start
    print("Multiprocessing took {0:.1f} seconds".format(multi_time))

samples = sampler.chain[:,200:, :].reshape((-1, ndim))

np.save('/home/tersenov/shear-pipe-peaks/output/NEW_HIGH_NOISE_constraints.npy',samples)

