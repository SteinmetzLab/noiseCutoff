function [nc_pass, cutoff, first_low_quantile] = noise_cutoff(amps, varargin)
    %     A new metric to determine whether a unit's amplitude distribution is cut off
    %     (at floor), without assuming a Gaussian distribution.
    %     This metric takes the amplitude distribution, computes the mean and std
    %     of an upper quartile of the distribution, and determines how many standard
    %     deviations away from that mean a lower quartile lies.
    %     Parameters
    %     ----------
    %     amps : ndarray_like
    %         The amplitudes (in uV) of the spikes.
    %     quantile_length : float
    %         The size of the upper quartile of the amplitude distribution.
    %     n_bins : int
    %         The number of bins used to compute a histogram of the amplitude
    %         distribution.
    %      nc_threshold: float
    %         the noise cutoff result has to be lower than this for a neuron to fail
    %     percent_threshold: float
    %         the first bin has to be greater than percent_threshold for neuron the to fail
    %     Returns
    %     -------
    %     cutoff : float
    %         Number of standard deviations that the lower mean is outside of the
    %         mean of the upper quartile.
    %     Examples
    %     --------
    %     1) Compute whether a unit's amplitude distribution is cut off
    %         >>> amps = spks_b['amps'][unit_idxs]
    %         >>> cutoff = bb.metrics.noise_cutoff(amps, quantile_length=.25, n_bins=100)

    nin = nargin;
    if nin < 4
        %set default params
        quantile_length=.25;
        n_bins=100;
        nc_threshold=5;
        percent_threshold=0.1;
    else
        quantile_length = varargin{1};
        n_bins = varargin{2};
        nc_threshold = varargin{3};
        percent_threshold = varargin{4};
    end

    cutoff = nan;
    first_low_quantile = nan;
    fail_criteria = ones(1, 'logical');

    if length(amps) > 1
        bins_list = linspace(0,max(amps),n_bins);
        [n, ~] = histcounts(amps, bins_list);
        idx_nz = find(n>0);  % indices of nonzeros
        [~, idx_peak_nz] = max(n(idx_nz));
        [~, idx_peak] = max(n);
        length_top_half = length(idx_nz) - idx_peak_nz;
        % the remaining part of the distribution, which we will compare the low quantile to
        high_quantile = 2 * quantile_length;
        % the first bin (index) of the high quantile part of the distribution
        high_quantile_start_ind = ceil(high_quantile * length_top_half + idx_peak);
        % bins to consider in the high qunatile (of all non-zero bins)
        indices_bins_high_quantile = high_quantile_start_ind:length(n);
        idx_use = find(n(indices_bins_high_quantile)>=1);

        if ~isempty(n(indices_bins_high_quantile))
            % mean of all amp values in high quantile bins
            mean_high_quantile = mean(n(indices_bins_high_quantile(idx_use)));
            std_high_quantile = std(n(indices_bins_high_quantile(idx_use)));
            if std_high_quantile > 0
                nonZero_n = n(n~=0); 
                first_low_quantile = nonZero_n(2); % taking the second bin
                cutoff = (first_low_quantile - mean_high_quantile) / std_high_quantile;
                peak_bin_height = max(n);
                percent_of_peak = percent_threshold * peak_bin_height;

                fail_criteria = (cutoff > nc_threshold) && (first_low_quantile > percent_of_peak);
                nc_pass = ~fail_criteria;
            else
                cutoff = nan;
                nc_pass = false;
            end
        else
            cutoff = nan;
            nc_pass = false;
        end
    else
        cutoff = nan;
        nc_pass = false;
    end

end