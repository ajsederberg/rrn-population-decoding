function [xc_vs_time, xc_vs_time_hist] = simdataNoiseCorr(proc_simdata)
% Computes the noise correlation by subtracting the trial-average response
% to each stimulus, then computing cross-correlations across cells at each
% point in time. 

has_stim = squeeze(any(proc_simdata.stimulus_info_array ~= 0, 1));
stim_ids = find(has_stim);

stim_vals = 0;
for ii = 1:length(stim_ids)
    stim_vals = stim_vals + proc_simdata.stimulus_info_array(:, stim_ids(ii))*ii;
end

% compute responses with stimulus-average subtracted
stim_resp = 0*proc_simdata.dr_out_peristim;
for ii = 1:length(stim_ids)
    stim_resp(stim_vals == ii, :, :) = bsxfun(@minus, proc_simdata.dr_out_peristim(stim_vals == ii, :, :), ... 
        mean(proc_simdata.dr_out_peristim(stim_vals == ii, :, :), 1));
    
end

%% compute noise correlations at each point in time
num_neur = size(stim_resp, 2);
xc_vs_time = zeros(num_neur, num_neur, size(stim_resp, 3));
xc_bins = linspace(-1, 1, 51);
xc_vs_time_hist = zeros(length(xc_bins)-1, size(xc_vs_time, 3));
for i_t = 1:size(stim_resp, 3)
    t_mat = corrcoef(stim_resp(:, :, i_t));
    xc_vs_time(:, :, i_t) = t_mat;

    xc_vs_time_hist(:, i_t) = histcounts(offDiag(t_mat), xc_bins);

end

