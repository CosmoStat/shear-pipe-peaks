"""
PLOTTING - A module for creating various types of plots.

:Name: summary_statistics.py 

:Description: This module provides functions for creating plots, including:
                - Plotting 2D maps using colormaps.
                - Plotting histograms of peak counts.
                - Plotting L1-norm histograms for different scales.

:Authors: Andreas Tersenov <atersenov@physics.uoc.gr>
"""

import matplotlib.pyplot as plt
import numpy as np

def plot_map(map_data, title='', cmap='inferno', vmin=None, vmax=None):
    """
    Plot a 2D map using a colormap.

    Parameters:
        map_data (numpy.ndarray): The 2D map data.
        title (str): Title of the plot.
        cmap (str): Colormap name.
        vmin (float): Minimum value for colormap scaling.
        vmax (float): Maximum value for colormap scaling.
    """
    plt.figure()
    img = plt.imshow(map_data, cmap=cmap, vmin=vmin, vmax=vmax, origin='lower')
    plt.title(title)
    plt.colorbar(img)
    plt.show()    

def plot_peak_count_histograms(kappa_th_center, peak_counts, title, xlabel, ylabel, log_scale=False):
    """
    Plot histograms of peak counts.

    Parameters:
        kappa_th_center (numpy.ndarray): Array of kappa threshold centers.
        peak_counts (numpy.ndarray or list of numpy.ndarray): Peak counts for each scale.
        title (str): Title of the plot.
        xlabel (str): Label for the x-axis.
        ylabel (str): Label for the y-axis.
        log_scale (bool, optional): Whether to use a logarithmic scale for the y-axis. Default is False.
    """
    plt.figure()
    
    if isinstance(peak_counts, list):  # Multiscale case
        for scale, peak_count in enumerate(peak_counts):
            plt.plot(kappa_th_center, peak_count, label=f'Scale {scale + 1}')
    else:  # Single-scale case
        plt.plot(kappa_th_center, peak_counts)
    
    plt.legend()
    if log_scale:
        plt.yscale('log')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True)
    plt.title(title)
    plt.show()

def plot_l1norm_histograms(bins, l1norm_histogram, title, xlabel, ylabel, xlim=None, log_scale=False):
    """
    Plot L1-norm histograms for each scale.

    Parameters:
        bins (list of numpy.ndarray): List of bin edges for each scale.
        l1norm_histogram (list of numpy.ndarray): List of L1-norm histograms for each scale.
        title (str): Title of the plot.
        xlabel (str): Label for the x-axis.
        ylabel (str): Label for the y-axis.
        xlim (tuple or list): Limits for the x-axis (e.g., xlim=(0, 6)).
        log_scale (bool): Whether to use a logarithmic scale for the y-axis.
    """
    plt.figure()

    for scale, l1norm_hist in enumerate(l1norm_histogram):
        plt.plot(bins[scale], l1norm_hist, label=f'Scale {scale}')

    plt.legend()
    plt.xticks()
    plt.yticks()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True)
    plt.title(title)
    
    if xlim:
        plt.xlim(xlim)
        
    if log_scale:
        plt.yscale('log')
        
    plt.show()