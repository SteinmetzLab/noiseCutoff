function plotNoiseCutoff(amps, spikeTimes)
%     '''
%     Parameters
%     ----------
%     amps : ndarray_like
%         The amplitudes (in uV) of the spikes. (spikes.amps)
%     spikeTimes : ndarray_like
%         The spike times. (spikes.times)
%     '''
%


[nc_pass, nc_value, first_bin_height] = noise_cutoff(amps);

n_bins=100;
percent_threshold=0.10;

f = figure; f.Color = 'w';
fp = f.Position;
f.Position = [fp(1) fp(2) 1300 369];

subplot(1,2,1);
scatter(spikeTimes,amps,'k.')
ylim([0 max(amps)]);
xlabel('Time (s)');
ylabel('Template amplitude (KS)');
% if ~isempty(params.cidx)
%     t1 = title(sprintf('Cluster #%d: FR=%.2f', params.cidx, firingRate));
% else
%     t1 = title(sprintf('FR=%.2f', firingRate));
% end
hold on;
box off;

if nc_pass
    title('Cutoff metric value: ' + string(round(nc_value, 2)),'Color', 'green')
else
    title('Cutoff metric value: ' + string(round(nc_value, 2)),'Color', 'red')
    
end



subplot(1,2,2)
h  = histogram(amps,n_bins, 'orientation', 'horizontal');
peak_bin_height = max(h.Values);
percent_label =  round(first_bin_height / peak_bin_height, 2) * 100;
xline(percent_threshold*peak_bin_height);
ylim([0,max(amps)]);
xlabel('Count');
hold on;
box off;


if nc_pass
    title('Low bin: ' + string(round(percent_label, 2)) + '{}% of peak ', 'Color', 'green');
else
    title('Low bin: ' + string(round(percent_label, 2)) + '{}% of peak ','Color', 'red');
    
end

end
