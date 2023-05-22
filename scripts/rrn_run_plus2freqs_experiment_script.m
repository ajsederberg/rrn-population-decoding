function sim_name = rrn_run_plus2freqs_experiment_script(network_code, pulse_width, rng_index, efficient_save, freq_str, f2f_type)
% For a given network specified by network_code and initialization seed
% rng_index, run the frequency discrimination test. 

% rng_index = 1;

network_param = generateNetworkParam(network_code, rng_index);

% this will look up the combination of frequencies/amplitudes of stimulus
% that gave roughly balanced inputs. 
ampfreq_simname = [network_param.network_tag ...
    '_AFsweep_pw' num2str(round(100*pulse_width))];
ampfreq_lookup = load([network_param.matsdir '/amplitude_lookup_' ampfreq_simname], ...
    'amp1_vals', 'amp2_vals', 'target_freqs', 'v_fr');
%%
% network tau = 1, so frequencies in freq_range are in units of that.
% Supposing a network tau of 0.020 s (characteristic frequency of 50 Hz) 
% Then 8 Hz (low end of freq. range in
% Licata and Raposo papers) is 0.16 per tau and 16 Hz (high end of freq
% range) is 0.32 per tau. 
[~, freq_ind16] = min(abs(ampfreq_lookup.target_freqs - 0.16));

% This part does not change, because we want to be on the same contour as
% for the 16 vs 32 experiment. 
% amplitude range from .15 to 1.5.  
amp_vals =  [ampfreq_lookup.amp1_vals(:, freq_ind16) ampfreq_lookup.amp2_vals(:, freq_ind16)];
ok_inds = ~any(isnan(amp_vals) | amp_vals > 10, 2);
if strcmp(f2f_type, 'low')
    amp_ind = find(ok_inds, 1, 'first');
elseif strcmp(f2f_type, 'med')
    a_inds = find(ok_inds);
    amp_ind = a_inds(round(end/2));
else
    amp_ind = find(ok_inds, 1, 'last' ); % take the last one (highest fr)
end


% Now we want to look up different frequencies, but on the same curve as
% these values [7, 14] or [10, 20]
if strcmp(freq_str, '714')
    % 7 Hz means 0.14 
    [~, freq_ind] = min(abs(ampfreq_lookup.target_freqs - 0.14));

elseif strcmp(freq_str, '1020')
    % 10 Hz means 0.20
      [~, freq_ind] = min(abs(ampfreq_lookup.target_freqs - 0.20));

end
% check for NaNs
amp_vals =  [ampfreq_lookup.amp1_vals(amp_ind, freq_ind) ampfreq_lookup.amp2_vals(amp_ind, freq_ind)];
ok_inds = ~any(isnan(amp_vals) | amp_vals > 10, 2);
if ~ok_inds
    disp(['amplitude look-up failed for ' freq_str])
    return;
end
    
% set frequencies
stim_param.freq_range = ampfreq_lookup.target_freqs(freq_ind)*[1 2];


stim_param.amp_range = [ampfreq_lookup.amp1_vals(amp_ind, freq_ind) ampfreq_lookup.amp2_vals(amp_ind, freq_ind)];
stim_param.num_reps = 400;
% ISI of 3 s is 150 tau. 
stim_param.isi = 150;

% pulse_width represents the timescale of upstream processing. This
% shouldbe set so that the input-triggered responses are about as
% modulated as observed in experiment. (cf Fig. 8, Licata et al)
% pulse_width = 0.5;
stim_param.pulse_width = pulse_width;
stim_param.ampfreq_combo = 'paired';
[sim_data_table, input_params] = simulateRRNwInputs(network_param, stim_param);
%%
% This is the f2fdiscExp : freq vs. 2*freq discrim experiment
sim_name = [network_param.network_tag ...
    '_f' freq_str 'f' f2f_type 'discExp_pw' num2str(round(100*pulse_width))];


%% basic processing
proc_simdata = computePeristimNetworkRout(sim_data_table, input_params, network_param);

if efficient_save
    % don't save the full sim_data_table [*CAUTION* ... ???]
    saveInParfor([network_param.matsdir 'sim_data_' sim_name], ...
        input_params, network_param, proc_simdata);    
else
    saveInParfor([network_param.matsdir 'sim_data_' sim_name], sim_data_table, ...
        input_params, network_param, proc_simdata, sim_name);
end
%% Plot basic overview of simuluation : one panel per stimulus combo 

nC = length(input_params.freq_range);
nR = length(input_params.amp_range);

y_range = [min(proc_simdata.dr_out_peristim(:)) max(proc_simdata.dr_out_peristim(:))];

% evoke_fr = zeros(num_stims, num_neurons);
t_vals = proc_simdata.peri_stim_time;
makeMyFigure(10*nC, 10*nR);
for i_sp = 1:nC*nR
    j_freq = mod(i_sp-1, nC) + 1;
    i_amp = floor((i_sp-1)/nR) + 1;
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
    
