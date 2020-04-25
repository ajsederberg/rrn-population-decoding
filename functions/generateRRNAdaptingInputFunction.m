function [input_fun, input_times, input_amplitudes, t0tf, stim_param] = generateRRNAdaptingInputFunction(stim_param)
% This function generates a handle to a function, input_fun, that provides
% the modulation and timing of input to the network at time t. 

% Start time of inputs is 100 (network tau = 1). If this is in the
% neighborhood of 20 ms, then 1s = 50 tau, and a reasonable ISI is 250 tau.
% The means that a frequency of 10-20 Hz is 0.2 - 0.4. 
freq_range = stim_param.freq_range;
amp_range = stim_param.amp_range;
isi = stim_param.isi;
num_reps = stim_param.num_reps;

num_freqs = length(freq_range);
num_amps = length(amp_range);

% trial length is 1 second, or 50 tau  
trial_length = 50;
stim_param.trial_length = trial_length;

% determine the time bin for generating input pulses;  these should be
% smaller than 0.1/max(freq_range)
dt_resolution = 0.0001;
dt = min(0.01, ceil(0.1/dt_resolution/max(freq_range))*dt_resolution);
stim_param.pulse_dt = dt;

if strcmp(stim_param.ampfreq_combo, 'sweep')
    num_stims = num_freqs*num_amps;
elseif strcmp(stim_param.ampfreq_combo, 'paired')
    num_stims = num_freqs;
elseif strcmp(stim_param.ampfreq_combo, 'adaptsweep')
    num_stims = 1;
    %%%%%%%%%%%%%% Adaptation options
    tau_adapt = stim_param.tau_adapt;   % input amplitude is a*exp(-dt/tau_adapt)
                                        % where dt is the time since the last
                                        % stim pulse
    %%%%%%%%%%%%%%%
else
    disp('Sweep across frequency and amplitude combos, or take freq/amp as pairs?')
    disp('Please set stim_param.ampfreq_combo to "sweep" or "paired"')
    input_fun = [];
    return;
end   


%%%%%%%%%%%%%%%
input_dt_ISI = cell(1, num_stims);
%%%%%%%%%%%%%%%
input_times = cell(1, num_stims);
stim_start_times = cell(1, num_stims);
input_amplitudes = cell(1, num_stims);
input_frequency = cell(1, num_stims);
% define trial start times, starting a t = 100 and spaced by isi (in units
% of network tau)
stim_param.burn_in_time = 100;
trial_start_times = stim_param.burn_in_time + (isi : isi : (num_stims*num_reps*isi));
input_ampfreq = zeros(length(trial_start_times), 2);

stim_param.trial_start_times = trial_start_times;
t0tf = [0 trial_start_times(end)+isi];

remaining_trial_start_times = trial_start_times;

%%
for i_stim = 1:num_stims
    if strcmp(stim_param.ampfreq_combo, 'sweep')
        i_freq = mod(i_stim-1, num_freqs) + 1;
        i_amp = ceil(i_stim/num_freqs);
    else
        i_freq = i_stim;
        i_amp = i_stim;
    end

    % select amplitude of this stimulus 
    input_amplitudes{i_stim} = amp_range(i_amp);
    input_frequency{i_stim} = freq_range(i_freq);
    % randomly select a set of input times
    inp_times = remaining_trial_start_times( sort(randperm(length(...
        remaining_trial_start_times), num_reps)));
    
    % remove these trial start times from the list of remaining trial start
    % times. 
    remaining_trial_start_times = setdiff(remaining_trial_start_times, ...
        inp_times);
    
    % find the indices of those times to save the stimulus parameters
    [~, trial_inds] = intersect(trial_start_times, inp_times);
    input_ampfreq(trial_inds, :) = repmat([amp_range(i_amp) freq_range(i_freq)], length(trial_inds), 1);
    
    % At each trial time: generate a poisson process with rate freq_range(i_freq). 
    
    % old way:
    % repeat trials (add background noise?) 
    % inp1_times = sort(randperm(length(trial_start_times), round(length(trial_start_times)/2)));
    % inp2_times = setdiff(1:length(trial_start_times), inp1_times);

    % % change input scheme. 
    % At each trial time: generate a poisson process with rate freq_range(i_freq). 
    stim_rate = freq_range(i_freq);
    num_bins = round(trial_length/(2*dt));
    impulse_window_times = 2*dt:2*dt:trial_length;
    stim_impulse_windows = rand(length(inp_times), num_bins) < stim_rate*trial_length/num_bins;

    stim_generation_probability = rand(length(inp_times), num_bins);
    [~, sgp_ord] = sort(stim_generation_probability, 2);
    s_impulse_time_ind = sort(sgp_ord(:, 1:round(stim_rate*trial_length)), 2);
    s_trial_ind = repmat((1:length(inp_times))', 1, size(s_impulse_time_ind, 2));
    
    s_impulse_time_ind = s_impulse_time_ind(:);
    s_trial_ind = s_trial_ind(:);
    
%     % get the trial index and the time index for each event - POISSON
%     % METHOD
%     [s_trial_ind, s_impulse_time_ind] = find(stim_impulse_windows);

    s_trial_times = inp_times(s_trial_ind) + impulse_window_times(s_impulse_time_ind);
    % s2_trial_times = trial_start_times(inp2_times(s2_trial_ind)) + impulse_window_times(s2_impulse_time_ind);
    input_times{i_stim} = sort(s_trial_times);
    
    %%%%%%%%%%%%%% Adaptation options
    if strcmp(stim_param.ampfreq_combo, 'adaptsweep')

        input_dt_ISI{i_stim} = [0 diff(input_times{i_stim})];   % dt_ISI is the time since the last
                                        % stim pulse
    else
        input_dt_ISI{i_stim} = 0*input_times{i_stim};
    end
    %%%%%%%%%%%%%%%
    stim_start_times{i_stim} = inp_times;
end
pulse_width = stim_param.pulse_width;  % this represents the filtering by upstream networks; setting to half of THIS network's time
                    % constant?  
% stim_param.pulse_width = pulse_width;
stim_param.input_frequency = input_frequency;
stim_param.input_ampfreq = input_ampfreq; 
stim_param.stim_start_times = stim_start_times;
% for Gaussian pulses, use this. 
% pulse_fun = @(t, a) exp(-t.^2/(2*a^2));
pulse_fun = @(t, a) (t.^2/a^2).*exp((t>0).*(1-t/a)).*(t > 0);

input_fun = @(t) sum(cell2mat(cellfun(@(in_time, in_pat) sum(pulse_fun(t - in_time, pulse_width), 2).*in_pat, ...
    input_times, input_amplitudes, 'uniformoutput', false)), 2);

stim_param.pulse_fun = pulse_fun;
% % if you want to use an exponential input function, use this
% pulse_fun = @(t, a) (t.^2/a^2).*exp(1-t/a).*(t > 0);
% input_fun = @(t) inefficientInputFun(t, input_times, input_amplitudes, pulse_fun, pulse_width);
