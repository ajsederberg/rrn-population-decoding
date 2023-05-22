% Script to demonstrate use of code set
% The main script ('rrn_main_script.m') will loop through many combinations
% of inputs characteristics, network characteristis, and so on. This is a
% demonstration script to show how to use the basic functionality.


% These were the parameters used to generate figures in paper.
% Parameters for spnorm_g400_n500_ei100_rng## are (Fig. 1-3):
% num_neurons = [500];
% inh_netwt_frac = [0.5];
% network_gain = [4 ];
% pulse_width_vals = [ .5 ];
% rng_vals = [1 2 3 4 5 6 7 8 9 10 11 12 13 17]; % note your random number
% generator may seed
% different simulations,
% but it should be pretty
% robust across seeds.

% Parameters for spnormE_g200_n500_ei100_rng## are (Fig. 4-5):
% num_neurons = [500];
% inh_netwt_frac = [0.5];
% network_gain = [2 ];
% pulse_width_vals = [ .5 ];
% rng_vals = 1:6;

%% 1. Define network meta-parameters of interest.
% Set the number of neurons and the fraction of inputs that are inhibitory
% (0.5 means balanced network). (80% of neurons are excitatory.)
% Network gain scales the strength of inputs.
num_neurons = 500;
inh_netwt_frac = 0.5;
network_gain = 4;

rng_val = 42; % set the random number generator seed for replicability

% These parameters are encoded in a key that is used to generate parameter
% structs for simulations below.
% The function for encoding parameters is:
% network_code = encodeNetworkMetaparams(num_neurons, network_gain, network_index, inh_netwt_frac)

%  "network_index" controls the connectivity patters
%   connect codes    1          2           3       4           5           6           7           8         9        10         11           12
% connect_tags = {'logNorm', 'splogNorm', 'norm', 'spnorm', 'spnormSL', 'spnorm2P', 'norm2PSL', 'spnormB', 'spfix', 'spnormE', 'spEnoise', 'spnormEL'};
% spnormE is same as spnorm, but only excitatory units receive inputs.
% Note that 'spnorm' takes the input cells out of the exc pool, but doesn't
% renormalized the e-i balance in the rest of the network. 'spnormB'
% rectifies this, by scaling up the excitatory connections.

% Set network_index to 10 for 'spnormE':
network_index = 10;
network_code = encodeNetworkMetaparams(num_neurons, network_gain, network_index, inh_netwt_frac);

% you can check what parameters will be generated for this simulation here
check_pars = generateNetworkParam(network_code, rng_val);
check_codes = check_pars.network_tag;
pulse_width = .5;
%% 2. Find input current values to normalize low/high frequencies
% Takes about 2 minutes when running on ~8 CPUs each with ~1.5 GB memory
eff_save = true;    % efficient save - does not save all files - might not be used anymore
tic
% Find the combinations of input amp/freq that generate the same overall firing rate
rrn_input_normalization_script(network_code, pulse_width, rng_val, eff_save);
toc

%% 3. Run the "freq2freq" experiment: (about 65 seconds for 500 neuron network)
% Options:

% input_current_level : this allowed selection of different overall firing
% rate. Options are 'low' or 'med' which select the lowest/median contour
% found in 'rrn_input_normalization_script'
% Default: '' (empty string) selects the highest-FR contour.
input_current_level = '';

tic
% Simulate experiments with default average input pulse rates (sometimes
% called frequency; note, however, inputs are not periodic)
sim_name = rrn_run_freq2freq_experiment_script(network_code, ...
    pulse_width, rng_val, eff_save, input_current_level);
% This function saves everything in the locations indicated by the network
% code.
toc
%% takes a long time
% this string is used to name files - so classifier results can be
% associated to the correct simulation files.
tic
input_type = ['f2f' input_current_level];
run_fitMultiClassifiersModelData(network_code, pulse_width, rng_val, input_type);
% outputs are saved to the results directory
toc
%% you can also change the input "frequencies", from the default of 8 vs. 16
% to 7 vs 14 (freq_str = '714';) or 10 vs 20 (freq_str = 1020);
freq_str = '714'; %options: '714', '1020'
input_type = ['f' freq_str 'f' input_current_level];

sim_name = rrn_run_plus2freqs_experiment_script(network_code, pulse_width, ...
    rng_val, eff_save, freq_str, input_current_level);
run_fitMultiClassifiersModelData(network_code, pulse_width, rng_index, input_type);

disp([freq_str ' run, t = ' num2str(toc, '%1.0f')])


