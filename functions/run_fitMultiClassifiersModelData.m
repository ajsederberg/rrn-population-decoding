function run_fitMultiClassifiersModelData(network_code, pulse_width, rng_index, input_type)
% Decoder fitting script - freq vs 2freq experiments (only two stims)
% Fits four classifiers: using all cells, using all excitatory
% cells, using a subset of excitatory cells, and using all inhibitory
% cells. 

network_param = lookUpNetworkParam(network_code, rng_index);

% matfile_dir = network_param.matsdir;
% results_dir = '/Volumes/home/asederberg6/projects/random_network_readouts/results/';

% num_neur = network_param.numNeur;
file_tag = network_param.network_tag; %

% Runs single-cell AUC analysis and saves reduced stim_resp matrix
activity_struct = preClassificationAnalysisF2F(network_code, pulse_width, rng_index, input_type);


%% output file name
output_filename = ['classifierAEIS_' file_tag '_' input_type 'discExp_pw' num2str(100*pulse_width)];
eff_save_dir = [network_param.save_dir 'matfiles/'];
full_save_dir = network_param.matsdir;

% if exist(output_filename, 'file')
%     disp([output_filename ' already exists'])
%    return;
% else
    fitMultiClassifiersModelData(activity_struct, output_filename, eff_save_dir, full_save_dir)
% end