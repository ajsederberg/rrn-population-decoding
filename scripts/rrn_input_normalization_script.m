function sim_name = rrn_input_normalization_script(network_code, pulse_width, rng_index, efficient_save)
% For a given network connectivity (class), simulate 5x5 combinations of
% input amplitude and frequency and find the combination of
% (amplitude,frequency) pairs A and B such that freq(A) = 2*freq(B). 

% num_neurons = 100;
% inh_netwt_frac = 0.55;
% network_gain = 1.1;
% network_index = 3;    
% network_code = encodeNetworkMetaparams(num_neurons, network_gain, network_index, inh_netwt_frac)

% rng_index = 1;

network_param = initializeNetworkParams(network_code, rng_index);

%% preliminary: find the basic range of amplitudes that puts the network in a dynamic region
% network tau = 1, so frequencies in freq_range are in units of that.
% Supposing a network tau of 0.020 s, then 8 Hz (low end of freq. range in
% Licata and Raposo papers) is 0.16 per tau and 16 Hz (high end of freq
% range) is 0.32 per tau. We'll expand slightly from this range; 
stim_param.freq_range = linspace(0.16, 0.32, 2);      % prelim: only look at low/high end of range
% amplitude range from .1 to 1.5. 
% % try a modified amp_range if it's a sparse matrix? 
%  on 2/25/19, changed upper bound to log10(5) and lowerbound to log10(.1)
stim_param.amp_range = logspace(log10(0.1), log10(5), 15);


if isfield(network_param, 'sparsity') && strcmp(network_param.iofun_name, 'tanh')
    if network_param.sparsity < .5
        stim_param.amp_range = stim_param.amp_range / network_param.sparsity / (2*pulse_width);
    end
end
% if strcmp(network_param.network_tag

stim_param.num_reps = 5; 
% ISI of 5 s is 250 tau. 
stim_param.isi = 150;

% pulse_width represents the timescale of upstream processing. This
% shouldbe set to sthat the input-triggered responses are about as
% modulated as observed in experiment. (cf Fig. 8, Licata et al)
% pulse_width = 0.5;
stim_param.pulse_width = pulse_width;
stim_param.ampfreq_combo = 'sweep';

[sim_data_table_prelim, input_params_prelim] = simulateRRNwInputs(network_param, stim_param);

proc_simdata_prelim = computePeristimNetworkRout(sim_data_table_prelim, input_params_prelim, network_param);
%%
tr_L = input_params_prelim.trial_length;
ps_t = proc_simdata_prelim.peri_stim_time;
halftime_stim_ave = mean(proc_simdata_prelim.stim_ave_dr_out(:, :, :, ps_t > 0 & ps_t < tr_L/2), 4);
stim_ave_by_cell = mean(proc_simdata_prelim.stim_ave_dr_out(:, :, :, ps_t > tr_L/2 & ...
    ps_t < tr_L), 4);

% calculate the fraction of cells active over 0.5 ("highly active" units)
frac_active = mean(stim_ave_by_cell > 0.5 | halftime_stim_ave > 0.5, 3);
%%
% find the amplitude at which the higher-frequency stimulus has less than
% 10% of units highly active (i.e. frac_active < 0.1)
amp_low = max(0.01, round(100*input_params_prelim.amp_range(max(1, sum(frac_active(:, 2) < 0.1)))/sqrt(5))/100);
amp_high = min(10*amp_low, input_params_prelim.amp_range(end));



%% then repeat trials to make a smooth average evoke FR surface. 
% select points int he frequency range that we will eventually study


if ~contains(network_param.network_tag, 'spEnoise')
    stim_param.freq_range = [0.12 0.14 0.16 0.2 0.24 0.28 0.32 0.40]; % linspace(0.1, 0.5, 5);  % now sample five points
    stim_param.amp_range = linspace(amp_low, amp_high, 6); % now sample 6 points

    stim_param.num_reps = 20;
else
    % the noisy simulations are much more memory-consuming. Only take 5
    % reps. 
    stim_param.freq_range = [0.14 0.16 0.2 0.28 0.32 0.40]; % linspace(0.1, 0.5, 5);  % now sample five points
    stim_param.amp_range = linspace(amp_low, amp_high, 5); % now sample 5 points

    stim_param.num_reps = 5;
end
tic
[sim_data_table, input_params] = simulateRRNwInputs(network_param, stim_param);
toc
%%
% This is a parameter sweep simulation. 
sim_name = [network_param.network_tag '_AFsweep_pw' num2str(round(100*pulse_width))];


%% basic processing
proc_simdata = computePeristimNetworkRout(sim_data_table, input_params, network_param);

% if efficient_save
    % don't save the full sim_data_table 
    saveInParfor([network_param.matsdir 'sim_data_' sim_name], ...
        input_params, network_param, proc_simdata);    
% else
%     saveInParfor([network_param.matsdir 'sim_data_' sim_name], sim_data_table, ...
%         input_params, network_param, proc_simdata, sim_name);
% end
%% Plot basic overview of simuluation : one panel per stimulus combo 

nC = length(input_params.freq_range);
nR = length(input_params.amp_range);

y_range = [min(proc_simdata.dr_out_peristim(:)) max(proc_simdata.dr_out_peristim(:))];
y_range = min(1, y_range);
y_range = max(0, y_range);
% evoke_fr = zeros(num_stims, num_neurons);
t_vals = proc_simdata.peri_stim_time;
makeMyFigure(10*nC, 10*nR);
for i_sp = 1:nC*nR
    j_freq = mod(i_sp-1, nC) + 1;
    i_amp = floor((i_sp-1)/nC) + 1;
    dr_vals = squeeze(proc_simdata.stim_ave_dr_out(i_amp, j_freq, :, :));

%     [i_amp, j_freq] = find(squeeze(proc_simdata.stimulus_info_array(i_sp, :, :)));
%     ij_sp = j_freq + (i_amp-1)*nC;
    
    subplot(nR, nC, i_sp)
    plot(t_vals, dr_vals)
    title({['amp = ' num2str(input_params.amp_range(i_amp), '%1.2g')]; ...
        ['freq = ' num2str(input_params.freq_range(j_freq))]})
    ylim(y_range)
end
print(gcf, '-dpdf', [network_param.plotdir 'overview_' sim_name]);

%% spline the average network firing rate surface to find matching freq/amp inputs
splineParamSweepSim(proc_simdata, input_params, network_param, sim_name);

% %% uncomment to debug the stimulus generation. 
% figure();
% for i_sp = 1:num_stims
% % i_sp = 1;
%     start_time = input_params.trial_start_times(i_sp);
%     stim_range = sim_data_table.t_out > start_time & sim_data_table.t_out < (start_time + input_params.trial_length);
%     
%     t_vals = start_time + linspace(0, 50, 1000)'; %sim_data_table.t_out(stim_range);
%     stim_vals = 0*t_vals;
%     for i_t = 1:length(t_vals)
%         stim_vals(i_t) = input_params.input_fun(t_vals(i_t));
%     end
%     subplot(nR, nC, i_sp)
%     plot(t_vals, stim_vals)
%     title({['amp = ' num2str(input_params.input_ampfreq(i_sp, 1), '%1.1g')]; ...
%         ['freq = ' num2str(input_params.input_ampfreq(i_sp, 2))]})
% end
    
