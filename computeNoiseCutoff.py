# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 15:53:41 2022

@author: Noam Roth
"""


def noise_cutoff(amps, quantile_length=.25, n_bins=100, nc_threshold=5, percent_threshold=0.10):
    """
    A new metric to determine whether a unit's amplitude distribution is cut off
    (at floor), without assuming a Gaussian distribution.
    This metric takes the amplitude distribution, computes the mean and std
    of an upper quartile of the distribution, and determines how many standard
    deviations away from that mean a lower quartile lies.
    Parameters
    ----------
    amps : ndarray_like
        The amplitudes (in uV) of the spikes.
    quantile_length : float
        The size of the upper quartile of the amplitude distribution.
    n_bins : int
        The number of bins used to compute a histogram of the amplitude
        distribution.
    n_low_bins : int
        The number of bins used in the lower part of the distribution (where
        cutoff is determined).
     nc_threshold: float
        the noise cutoff result has to be lower than this for a neuron to fail
    percent_threshold: float
        the first bin has to be greater than percent_threshold for neuron the to fail
    Returns
    -------
    cutoff : float
        Number of standard deviations that the lower mean is outside of the
        mean of the upper quartile.
    See Also
    --------
    missed_spikes_est
    Examples
    --------
    1) Compute whether a unit's amplitude distribution is cut off
        >>> amps = spks_b['amps'][unit_idxs]
        >>> cutoff = bb.metrics.noise_cutoff(amps, quantile_length=.25, n_bins=100)
    """
    cutoff = np.float64(np.nan)
    first_low_quantile = np.float64(np.nan)
    fail_criteria = np.ones(1).astype(bool)[0]

    if amps.size > 1:  # ensure there are amplitudes available to analyze
        bins_list = np.linspace(0, np.max(amps), n_bins)  # list of bins to compute the amplitude histogram
        n, bins = np.histogram(amps, bins=bins_list)  # construct amplitude histogram
        idx_peak = np.argmax(n)  # peak of amplitude distribution
        # don't count zeros #len(n) - idx_peak, compute the length of the top half of the distribution -- ignoring zero bins
        length_top_half = len(np.where(n[idx_peak:-1] > 0)[0])
        # the remaining part of the distribution, which we will compare the low quantile to
        high_quantile = 2 * quantile_length
        # the first bin (index) of the high quantile part of the distribution
        high_quantile_start_ind = int(np.ceil(high_quantile * length_top_half + idx_peak))
        # bins to consider in the high quantile (of all non-zero bins)
        indices_bins_high_quantile = np.arange(high_quantile_start_ind, len(n))
        idx_use = np.where(n[indices_bins_high_quantile] >= 1)[0]

        if len(n[indices_bins_high_quantile]) > 0:  # ensure there are amplitudes in these bins
            # mean of all amp values in high quantile bins
            mean_high_quantile = np.mean(n[indices_bins_high_quantile][idx_use])
            std_high_quantile = np.std(n[indices_bins_high_quantile][idx_use])
            if std_high_quantile > 0:
                first_low_quantile = n[(n != 0)][1]  # take the second bin
                cutoff = (first_low_quantile - mean_high_quantile) / std_high_quantile
                peak_bin_height = np.max(n)
                percent_of_peak = percent_threshold * peak_bin_height

                fail_criteria = (cutoff > nc_threshold) & (first_low_quantile > percent_of_peak)

    nc_pass = ~fail_criteria
    return nc_pass, cutoff, first_low_quantile