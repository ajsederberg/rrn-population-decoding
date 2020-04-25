function run_fitMultiFreqThrClassifiersModelData(network_code, pulse_width, rng_index, input_current_level)
% Decoder fitting script - take results of 'f2f', 'f714', and 'f1020' and
% fits multiple classifiers with different thresholds. In 'f2f', the
% frequencies are 8 and 16. So we have responses to 7, 8, 10, 14, 16, and
% 20 Hz stimuli; the decision threshold can be placed at [9, 12, 15]
% Fits four classifiers: using all cells, using all excitatory
% cells, using a subset of excitatory cells, and using all inhibitory
% cells. 

network_param = lookUpNetworkParam(network_code, rng_index);

% matfile_dir = network_param.matsdir;
% results_dir = '/Volumes/home/asederberg6/projects/random_network_readouts/results/';

% num_neur = network_param.numNeur;
file_tag = network_param.network_tag; %


% Plan: generate activity struct for each of 'f2f', 'f714', 'f1020'. In
% this function, mix up the activity structs for each frequency decision
% threshold. Generate a threshold-specific file name ("output_filename")
% for each of these thresholds. 
% Runs single-cell AUC analysis and saves reduced stim_resp matrix
f2f_types = {'f2f', 'f714f', 'f1020f'};
all_act_struct = cell(size(f2f_types));
for i_ft = 1:length(f2f_types)
    input_type = [f2f_types{i_ft} input_current_level];
    
    activity_struct = preClassificationAnalysisF2F(network_code, pulse_width, rng_index, input_type);
    all_act_struct{i_ft} = activity_struct;
end
keyboard

%% output file name
output_filename = ['classifierFTh3_' file_tag '_' input_type 'discExp_pw' num2str(100*pulse_width)];
eff_save_dir = [network_param.save_dir 'matfiles/'];
full_save_dir = network_param.matsdir;

if exist(output_filename, 'file')
    disp([output_filename ' already exists'])
   return;
else
    % run the usual classifier, but feed it different activity structs
    % representing different choices of thresholds NOT DONE YET - 02/18/19
    fitMultiClassifiersModelData(activity_struct, output_filename, eff_save_dir, full_save_dir)
end