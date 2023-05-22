% Script to manage all simulations. 
% Parameters for spnorm_g400_n500_ei100_rng## are: 
% num_neurons = [500];
% inh_netwt_frac = [0.5];
% network_gain = [4 ];
% pulse_width_vals = [ .5 1.5 5];

%% Define network meta-parameters of interest. 
% num_neurons = 100;
% inh_netwt_frac = 0.55;
% network_gain = 1.1;
% network_index = 3; 

num_neurons = [500];
inh_netwt_frac = [0.5];
network_gain = [1 2];

[nn_g, in_g, ng_g] = meshgrid(num_neurons, inh_netwt_frac, network_gain);
%  network_gain, network_index, inh_netwt_frac
%   connect codes    1          2           3       4           5           6           7           8         9        10         11           12 
% connect_tags = {'logNorm', 'splogNorm', 'norm', 'spnorm', 'spnormSL', 'spnorm2P', 'norm2PSL', 'spnormB', 'spfix', 'spnormE', 'spEnoise', 'spnormEL'};
% spnormE is same as spnorm, but only excitatory units receive inputs. 
% Note that 'spnorm' takes the input cells out of the exc pool, but doesn't
% renormalized the e-i balance in the rest of the network. 'spnormB'
% rectifies this, by scaling up the excitatory connections. 
network_codes = encodeNetworkMetaparams(nn_g, ng_g, 12 + 0*nn_g, in_g);

network_codes = network_codes(:);

check_pars = cellfun(@(x) generateNetworkParam(x, 1), num2cell(network_codes), ...
    'UniformOutput', false);
check_codes = cellfun(@(x) x.network_tag, check_pars, 'UniformOutput', false);

pulse_width_vals = [ .5 ];
%% This block runs the entire experiment, from input normalization to 
% the multiple frequencies in the same network. 
eff_save = true;
inpCurrent_lowhigh_type = {''};
freq_strs = {}; %{'714', '1020'};        % runs these IN ADDITION to the usual 'f2f'. 
for rng_index = 2:4 %[4 5 8 9 10 13 17] %[ 12:-1:1] %[1 2 3]% 6 7 11 12] %1:15
    for i_nc = 1:length(network_codes) 
        for i_pw = 1:1 %length(pulse_width_vals)

            tic;
            pulse_width = pulse_width_vals(i_pw);
            network_code = network_codes(i_nc);
            try 
                % Find the combinations of input amp/freq that generate the same overall firing rate
                rrn_input_normalization_script(network_code, pulse_width, rng_index, eff_save);
                %%
                %%
                disp('inputs normalized')
                toc
            for i_ft = 1:length(inpCurrent_lowhigh_type)
                input_current_level = inpCurrent_lowhigh_type{i_ft};
                % Simulate experiments
                sim_name = rrn_run_freq2freq_experiment_script(network_code, pulse_width, rng_index, eff_save, input_current_level);
                %%
                input_type = ['f2f' input_current_level];
                run_fitMultiClassifiersModelData(network_code, pulse_width, rng_index, input_type);
                disp('f2f run')
                toc
                close all
                %%
                for i_fs = 1:length(freq_strs)
                    freq_str = freq_strs{i_fs};
%%
                    try sim_name = rrn_run_plus2freqs_experiment_script(...
                        network_code, pulse_width, rng_index, eff_save, freq_str, input_current_level);
                        input_type = ['f' freq_str 'f' input_current_level];
                        run_fitMultiClassifiersModelData(network_code, pulse_width, rng_index, input_type);
                    end
                    disp([freq_str ' run, t = ' num2str(toc, '%1.0f')])
                end
                disp(['done with ' num2str(network_code) ' in ' num2str(toc, '%1.0f') ' seconds'])

%             end
            end
        catch 
                
                disp(['failed: ' num2str(network_code)])
        end
        end
    end
end
beep

%% This block runs a fast version of input normalization 

eff_save = true;
inpCurrent_lowhigh_type = {''};
% freq_strs = {};        % runs these IN ADDITION to the usual 'f2f'. 
for rng_index = 1:1 %[ 12:-1:1] %[1 2 3]% 6 7 11 12] %1:15
    for i_nc = 1:length(network_codes) 
        for i_pw = 1:2 %length(pulse_width_vals)

            tic;
            pulse_width = pulse_width_vals(i_pw);
            network_code = network_codes(i_nc);
%             try 
                % Find the combinations of input amp/freq that generate the same overall firing rate
                rrn_input_normalization_script(network_code, pulse_width, rng_index, eff_save);
                %%
                disp('inputs normalized')
        end
    end
end

%% Simulate experiments - separate from input normalization. 
eff_save = true;
input_current_level = '';
freq_strs = {'714', '1020'};
for i_fs = 1:length(freq_strs)
    freq_str = freq_strs{i_fs};
    for i_nc = 1:length(network_codes)
        for i_pw = 1:1 %length(pulse_width_vals)
            for rng_index = 17 %[1 2 6 7 10]
                pulse_width = pulse_width_vals(i_pw);
                network_code = network_codes(i_nc);

                sim_name = rrn_run_freq2freq_experiment_script(network_code, pulse_width, rng_index, eff_save);
                input_type = ['f2f' input_current_level];
                run_fitMultiClassifiersModelData(network_code, pulse_width, rng_index, input_type);

                try sim_name = rrn_run_plus2freqs_experiment_script(...
                        network_code, pulse_width, rng_index, eff_save, freq_str, input_current_level);
                    input_type = ['f' freq_str 'f' input_current_level];
                    run_fitMultiClassifiersModelData(network_code, pulse_width, rng_index, input_type);

                end
            end
        end
    end
end


%% Fit decoders - separate from input normalization and network simulation.
eff_save = true;
f2f_types = {'f2f'}; %{'f1020f'}; %{'f714f', 'f1020f'}; %      % 'f2f' : original. 'f2flow' : takes lower firing rate set point
                         % 'f714f' : compares 7 Hz to 14 Hz
input_levels = {''};
pulse_width_vals = [0.5]; % 1.5 5];
for i_nc = 1:length(network_codes)
    for i_pw = 1:length(pulse_width_vals)
        for i_ft = 1:length(f2f_types)
            f2f_type = f2f_types{i_ft};
            for i_icl = 1:length(input_levels)
                input_current_level = input_levels{i_icl};

                for rng_index = 2:7 % 8:8 % [3:10 11 12 13 17]  % 1 %[ 1 2 7 10] %1:10 [4 5 8 9]
                    pulse_width = pulse_width_vals(i_pw);
                    network_code = network_codes(i_nc);
                    try
                        disp(['Running for network code ' ...
                            num2str(network_code) ' pulsewidth ' num2str(pulse_width*100) ... 
                            ' rng index ' num2str(rng_index) input_current_level])
                        input_type = [f2f_type input_current_level];
                        run_fitMultiClassifiersModelData(network_code, pulse_width, rng_index, input_type);
                    catch
                        disp(['results missing, probably not run yet. Network code ' ...
                            num2str(network_code) ' pulsewidth ' num2str(pulse_width*100) ... 
                            ' rng index ' num2str(rng_index) input_current_level])
                    end
                end
            end
        end
    end
end

%% Fit decoders for multi-thresholds. Requires that {'f2f', 'f714f', 'f1020f'}
%  have all been previously run. 
eff_save = true;
                         % 'f714f' : compares 7 Hz to 14 Hz
input_levels = {''};
pulse_width_vals = [0.5];% 1.5 5];
for i_nc = 1:length(network_codes)
    for i_pw = 1:length(pulse_width_vals)
            for i_icl = 1:length(input_levels)
                input_current_level = input_levels{i_icl};

                for rng_index =  1:13 % 1 %[ 1 2 7 10] %1:10 [4 5 8 9]
                    pulse_width = pulse_width_vals(i_pw);
                    network_code = network_codes(i_nc);
                    try
                        run_fitMultiFreqThrClassifiersModelData(network_code, pulse_width, rng_index, input_current_level);
                    catch
                        disp(['results missing, probably not run yet. Network code ' ...
                            num2str(network_code) ' pulsewidth ' num2str(pulse_width*100) ... 
                            ' rng index ' num2str(rng_index) input_current_level])
                    end
                end
            end
        
    end
end
   
%% Adaptive : re-run the classifier fitting, for all results generated after
% some point in time (specified). This works by looking in the directory of
% results for particular strings. 

% Pseudo code: determine what type of sims need to be re-analyzed
% (spnorm_g400_n500 , for example) and get all networks with that network
% code.

% if the results file already exists but was created before some time
% point, re-run the analysis. 
%% Analyze decoder weights