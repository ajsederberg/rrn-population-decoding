function processed_data = computePeristimNetworkRout(sim_data_table, input_params, network_param)
% [r_out_peristim, dr_out_peristim, peri_event_resps_cell, detrended_peri_event_resps_cell] = computePeristimNetworkRout(sim_data_table, input_params, network_param)
% This function computes the stimulus-sorted responses from the full
% timeseries in sim_data_table using the stimulus onset times in
% input_params. 

% First initialize variables. In units of network tau (=1, est. 20 ms)
dt = 0.2;   % .2*20 ms = 4 ms res. 
% previously taken to 2*trial_length; this isn't really necessary. 
pst_length = dt/2 + input_params.trial_length + max(1./input_params.freq_range);
peri_stim_time = dt:dt:pst_length;
num_t = length(peri_stim_time);
num_neurons = network_param.numNeur;
num_reps = input_params.num_reps;
num_freqs = length(input_params.freq_range);
num_amps = length(input_params.amp_range);

if strcmp(input_params.ampfreq_combo, 'sweep')
    num_stims = num_freqs*num_amps*num_reps;
elseif strcmp(input_params.ampfreq_combo, 'paired')
    num_stims = num_freqs*num_reps;
else
    disp('invalid ampfreq_combo value ( must be "sweep" or "paired"')
end
    
r_out_peristim = zeros(num_stims, num_neurons, num_t);

% stimulus ID info: has a "true" in the amplitdue, frequency position
% corresponding to the i-th stimulus
stimulus_info_array = zeros(num_stims, num_amps, num_freqs);
% look up the precise input times for each trial
precise_impulse_times = cell(num_stims, 1);
precise_impuse_peritimes = cell(num_stims, 1);
input_freq_vals = cell(num_stims, 1);
input_amp_vals = cell(num_stims, 1);
% all input times
all_input_times = cell2mat(input_params.input_times);

% isolate the baseline time: pre-stimulus activity after 100tau burn-in. 
pre_stim_time = sim_data_table.t_out < input_params.burn_in_time;
baseline_time = sim_data_table.t_out > 50 & pre_stim_time;

% compute the static fixed point firing rates 
r_out_baseline = network_param.fun_io(sim_data_table.x_out(baseline_time, :));
baseline_activity = mean(r_out_baseline, 1);


% plot responses (single trial) for each stim
for i_sp = 1:num_stims
    
    start_time = input_params.trial_start_times(i_sp);
    stim_range = sim_data_table.t_out >= start_time-dt & sim_data_table.t_out <= (start_time + pst_length + dt);
    
    t_vals = sim_data_table.t_out(stim_range);
    r_vals = network_param.fun_io(sim_data_table.x_out(stim_range, :));
    
    % setting linear/extrap prevents the edge t-vals (0 or
    % peri_stim_time(end)) from being "out of bounds"
    r_interp = interp1(t_vals - start_time, r_vals, peri_stim_time, 'linear', 'extrap');
    
    
    r_out_peristim(i_sp, :, :) = permute(r_interp, [3 2 1]);
    
    amp_val = input_params.input_ampfreq(i_sp, 1);
    freq_val = input_params.input_ampfreq(i_sp, 2);
    
    input_freq_vals{i_sp} = freq_val;
    input_amp_vals{i_sp} = amp_val;
    
    i_amp = input_params.amp_range == amp_val;
    j_freq = input_params.freq_range == freq_val;
    stimulus_info_array(i_sp, i_amp, j_freq) = true;
    
    ait_inds = all_input_times - start_time < input_params.trial_length & ...
        all_input_times - start_time > 0;
    precise_impulse_times{i_sp} = all_input_times(ait_inds);
    precise_impuse_peritimes{i_sp} = num2cell(all_input_times(ait_inds) - start_time)';
end

dr_out_peristim = bsxfun(@minus, r_out_peristim, reshape(baseline_activity, [1 num_neurons, 1]));
%% compute the peri-event responses

% this has the # of time points for 1 period of the frequency of each input
peri_event_time = cellfun(@(x) length(dt:dt:1/x), input_freq_vals, ...
    'UniformOutput', false);
% this has the start times of each stimulus
% stim_start_times = (input_params.stim_start_times');
% this puts the peristimulus responses into cell array format, with each
% cell representing 1 trial
r_out_peristim_cell = mat2cell(dr_out_peristim, ones(num_stims, 1), num_neurons, num_t);
% % these are the precise input times of each impulse
% input_times = cellfun(@(x) num2cell(x), input_params.input_times(:), 'uniformoutput', false);

event_fun = @(r_out, freq_T, peri_inp_times) ...
    cellfun(@(inp_t) r_out(:, :, sum(peri_stim_time < inp_t) + (1:freq_T)), ...
    peri_inp_times, 'uniformoutput', false);
    
peri_event_resps = cellfun(@(rr, f_T, i_t) cell2mat(event_fun(rr, f_T, i_t)), ...
    r_out_peristim_cell, peri_event_time, precise_impuse_peritimes , ...
    'UniformOutput', false);

% detrend
peri_event_resps_cell = cellfun(@(x) mat2cell(x, size(x, 1), ones(size(x, 2),1), size(x, 3)), ...
    peri_event_resps, 'UniformOutput', false);
detrended_peri_event_resps_cell = cellfun(@(x) cell2mat(cellfun(@(y) ...
    permute(detrend(permute(squeeze(y), [2 1])), [2 3 1]), x, 'UniformOutput', ...
    false)), peri_event_resps_cell, 'UniformOutput', false);


% Put in struct:
% processed_data.r_out_peristim = r_out_peristim; % rout_peristim can be
                                                  % computed from the baseline
processed_data.baseline_activity = baseline_activity;
processed_data.dr_out_peristim = dr_out_peristim;
processed_data.peri_stim_time = peri_stim_time;
processed_data.stimulus_info_array = stimulus_info_array;
processed_data.precise_impulse_times = precise_impulse_times;
processed_data.precise_impuse_peritimes = precise_impuse_peritimes;

% These were never used; don't waste space to save it until needed
% processed_data.peri_event_resps_cell = peri_event_resps_cell; 
% processed_data.detrended_peri_event_resps_cell = detrended_peri_event_resps_cell;
% processed_data.peri_event_time = peri_event_time;


%% if multiple trials present, compute the trial-average
stim_ave_dr_out = zeros(num_amps, num_freqs, size(dr_out_peristim, 2), size(dr_out_peristim, 3));
for i_a = 1:num_amps
    for i_f = 1:num_freqs
        stim_trial_inds = stimulus_info_array(:, i_a, i_f) == 1;
        if any(stim_trial_inds)
            stim_ave = mean(dr_out_peristim(stim_trial_inds, :, :), 1);
        
            stim_ave_dr_out(i_a, i_f, :, :) = permute(stim_ave, [4 1 2 3]);
        end
    end
end
processed_data.stim_ave_dr_out = stim_ave_dr_out;

% mean_peri_event_resp = cellfun(@(x) squeeze(mean((x), 1)), ...
%     peri_event_resps, 'UniformOutput', false);
% 
% mean_peri_event_resp_detrended = cellfun(@(x) squeeze(mean((x), 1)), ...
%     detrended_peri_event_resps_cell, 'UniformOutput', false);

% compute modulation indices?