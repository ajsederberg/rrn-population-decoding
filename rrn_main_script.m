% Script to run all simulations and analyses. 
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

%% Define network meta-parameters of interest. 
% num_neurons = 100;
% inh_netwt_frac = 0.55;
% network_gain = 1.1;
% network_index = 3; 

% these can be vectors. this script will run over all combinations
num_neurons = [500];
inh_netwt_frac = [0.5];
network_gain = [4];

[nn_g, in_g, ng_g] = meshgrid(num_neurons, inh_netwt_frac, network_gain);
%  network_gain, network_index, inh_netwt_frac
%   connect codes    1          2           3       4           5           6           7           8         9        10         11           12 
% connect_tags = {'logNorm', 'splogNorm', 'norm', 'spnorm', 'spnormSL', 'spnorm2P', 'norm2PSL', 'spnormB', 'spfix', 'spnormE', 'spEnoise', 'spnormEL'};
% spnormE is same as spnorm, but only excitatory units receive inputs. 
% Note that 'spnorm' takes the input cells out of the exc pool, but doesn't
% renormalized the e-i balance in the rest of the network. 'spnormB'
% rectifies this, by scaling up the excitatory connections. 
network_codes = encodeNetworkMetaparams(nn_g, ng_g, 4 + 0*nn_g, in_g);

network_codes = network_codes(:);

check_pars = cellfun(@(x) generateNetworkParam(x, 1), num2cell(network_codes), ...
    'UniformOutput', false);
check_codes = cellfun(@(x) x.network_tag, check_pars, 'UniformOutput', false);

pulse_width_vals = [ .5 ];
%% This block runs the entire experiment, from input normalization to 
% the multiple frequencies in the same network. 
eff_save = true;
inpCurrent_lowhigh_type = {''};
freq_strs = {'714', '1020'};        % runs these IN ADDITION to the usual 'f2f'. 
for rng_index = [1 2 3 4 5 6 7 8 9 10 13 17] 
    for i_nc = 1:length(network_codes) 
        for i_pw = 1:length(pulse_width_vals)

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
                % Simulate experiments : this returns the simulated network
                % activity. 
                sim_name = rrn_run_freq2freq_experiment_script(network_code, pulse_width, rng_index, eff_save, input_current_level);
                %%
                % These analyze simulated experiments (selectivity,
                % classifier performance, etc. )
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
