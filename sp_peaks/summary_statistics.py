"""SUMMARY STATISTICS.

:Name: summary_statistics.py
 
:Description: This package contains methods to compute summary statistics from mass maps                      
                                                                               
:Authors: Lucie Baumont <lucie.baumont@cea.fr> Martin Kilbinger <martin.kilbinger@cea.fr> 
"""

def compute_single_scale_peak_counts(snr_map, kappa_snr):
    """
    Compute peak counts for a single SNR map.

    Parameters:
        snr_map (numpy.ndarray): SNR map.
        kappa_snr (numpy.ndarray): Array of kappa values corresponding to the SNR map.

    Returns:
        kappa_th_center_snr (numpy.ndarray): Array of kappa threshold centers for peak counts.
        peak_counts (numpy.ndarray): Peak counts for the given SNR map.
    """
    kappa_th_center_snr = 0.5 * (kappa_snr[:-1] + kappa_snr[1:])
    peak_counts = peaks.peaks_histogram(snr_map, kappa_snr)[0]
    return kappa_th_center_snr, peak_counts

def compute_multiscale_peak_counts(snr_maps, kappa_snr):
    """
    Compute peak counts for each wavelet scale of SNR maps.

    Parameters:
        snr_maps (list of numpy.ndarray): List of SNR maps for each scale.
        kappa_snr (numpy.ndarray): Array of kappa values corresponding to SNR maps.

    Returns:
        kappa_th_center_snr (numpy.ndarray): Array of kappa threshold centers for peak counts.
        peak_counts (list of numpy.ndarray): List of peak counts for each scale.
    """
    kappa_th_center_snr = 0.5 * (kappa_snr[:-1] + kappa_snr[1:])
    
    peak_counts = [peaks.peaks_histogram(snr_map, kappa_snr)[0] for snr_map in snr_maps]

    return kappa_th_center_snr, peak_counts

