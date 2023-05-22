% Script to manage all simulations. 

%% Define network meta-parameters of interest. 
% num_neurons = 100;
% inh_netwt_frac = 0.55;
% network_gain = 1.1;
% network_index = 3; 

num_neurons = [50];
inh_netwt_frac = [0.5];
network_gain = [4 ];

[nn_g, in_g, ng_g] = meshgrid(num_neurons, inh_netwt_frac, network_gain);
%  network_gain, network_index, inh_netwt_frac
network_codes = encodeNetworkMetaparams(nn_g, ng_g, 4 + 0*nn_g, in_g);

network_codes = network_codes(:);

check_pars = cellfun(@(x) generateNetworkParam(x, 1), num2cell(network_codes), ...
    'UniformOutput', false);
check_codes = cellfun(@(x) x.network_tag, check_pars, 'UniformOutput', false);

pulse_width_vals = [  1.5 ];
%% This block runs the entire experiment, from input normalization to 
% the multiple frequencies in the same network. 
eff_save = true;
inpCurrent_lowhigh_type = {'', 'med','low'};
freq_strs = []; %{'714', '1020'};        % runs these IN ADDITION to the usual 'f2f'. 
for rng_index = 1:1
    for i_nc = 1:length(network_codes) 
        for i_pw = 1:length(pulse_width_vals)

            tic;
            pulse_width = pulse_width_vals(i_pw);
            network_code = network_codes(i_nc);
%             try 
                % Find the combinations of input amp/freq that generate the same overall firing rate
%                 rrn_input_normalization_script(network_code, pulse_width, rng_index, eff_save);
                rrn_input_adaptation_script(network_code, pulse_width, rng_index, eff_save);
                %%
            for i_ft = 1:length(inpCurrent_lowhigh_type)
                input_current_level = inpCurrent_lowhigh_type{i_ft};
                % Simulate experiments
                sim_name = rrn_run_freq2freq_experiment_script(network_code, pulse_width, rng_index, eff_save, input_current_level);
                %%
                input_type = ['f2f' input_current_level];
                run_fitMultiClassifiersModelData(network_code, pulse_width, rng_index, input_type);
                close all
                for i_fs = 1:length(freq_strs)
                    freq_str = freq_strs{i_fs};
%%
                    try sim_name = rrn_run_plus2freqs_experiment_script(...
                        network_code, pulse_width, rng_index, eff_save, freq_str, input_current_level);
                        input_type = ['f' freq_str 'f' input_current_level];
                        run_fitMultiClassifiersModelData(network_code, pulse_width, rng_index, input_type);
                    end
                end
                disp(['done with ' num2str(network_code) ' in ' num2str(toc, '%1.0f') ' seconds'])
%             catch 
                
%                 disp(['failed: ' num2str(network_code)])
%             end
            end
        end
    end
end
beep
%% Simulate experiments - separate from input normalization. 
eff_save = true;
input_current_level = '';
freq_strs = {'714', '1020'};
for i_fs = 1:2
    freq_str = freq_strs{i_fs};
    for i_nc = 1:length(network_codes)
        for i_pw = 1:length(pulse_width_vals)
            for rng_index = [1 2 6 7 10]
                pulse_width = pulse_width_vals(i_pw);
                network_code = network_codes(i_nc);

    %             sim_name = rrn_run_freq2freq_experiment_script(network_code, pulse_width, rng_index, eff_save);
                try sim_name = rrn_run_plus2freqs_experiment_script(...
                        network_code, pulse_width, rng_index, eff_save, freq_str, input_current_level);
                end
            end
        end
    end
end


%% Fit decoders - separate from input normalization. 
eff_save = true;
f2f_types ={'f2f'}; % {'f714f', 'f2f', 'f1020f'};      % 'f2f' : original. 'f2flow' : takes lower firing rate set point
                         % 'f714f' : compares 7 Hz to 14 Hz
input_levels = {''};
pulse_width_vals = [0.5 1.5 5];
for i_nc = 1:length(network_codes)
    for i_pw = 3:length(pulse_width_vals)
        for i_ft = 1:1 %length(f2f_types)
            f2f_type = f2f_types{i_ft};
            for i_icl = 1:length(input_levels)
                input_current_level = input_levels{i_icl};

                for rng_index =  3 %[ 1 2 7 10] %1:10 [4 5 8 9]
                    pulse_width = pulse_width_vals(i_pw);
                    network_code = network_codes(i_nc);
                    try
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
         
%% Analyze decoder weights