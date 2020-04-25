function activity_struct = preClassificationAnalysisF2F(network_code, pulse_width, rng_index, input_type)
% Analysis script for a single freq vs 2freq experiments (only two stims)
% Performs single-cell AUC analysis and prepares data for use with the
% classifier. 

network_param = lookUpNetworkParam(network_code, rng_index);

matfile_dir = network_param.matsdir;
% results_dir = '/Volumes/home/asederberg6/projects/random_network_readouts/results/';

% num_neur = network_param.numNeur;
file_tag = network_param.network_tag; %


% load('/Users/audreysederberg/Dropbox/postdoctoral projects/stanley lab work/code/theory code/results/norm_g100_n100_ei80_rng1/matfiles/sim_data_norm_g100_n100_ei80_rng1_f2fdiscExp_pw150.mat')
x_res = load([matfile_dir '/sim_data_' file_tag '_' input_type 'discExp_pw' num2str(100*pulse_width) '.mat']);

proc_simdata = x_res.proc_simdata;
input_params = x_res.input_params;

network_param = x_res.network_param;

%% can compute and save the AUC and network stats
computeSelXCJ_stats(proc_simdata, network_param, input_type, pulse_width, false);

t_vec = proc_simdata.peri_stim_time;

%%
if strcmp(input_params.ampfreq_combo, 'paired')
    num_stims = size(proc_simdata.stimulus_info_array, 2);   
    
    % for the paired stim, the stim infor array has ones on the diagonal
    % entries corresponding to what the stim ID is. 
    trial_inds = cellfun(@(x) find(proc_simdata.stimulus_info_array(:, x, x) == 1), ...
        num2cell(1:num_stims)', 'UniformOutput', false);
%     
%     % find lower-rate and higher-rate trials. 
%     trial_inds = {find(cellfun(@length, proc_simdata.precise_impulse_times) <= 6); find(cellfun(@length, proc_simdata.precise_impulse_times) > 6)};
%     
    % this puts all stim 1's first followed by all stim 2's
    stim_y = cell2mat(cellfun(@(x) x*ones(sum(proc_simdata.stimulus_info_array(:, x, x) == 1), 1), ...
        num2cell(1:num_stims)', 'UniformOutput', false));
end



stim_resps_fullt = cell2mat(cellfun(@(x) proc_simdata.dr_out_peristim(x, :, :), trial_inds, ...
    'UniformOutput', false));


% Downsample: fit a new classifer every 1 tau 
c_timestep = 1; % in units of tau, how often to fit a new classifier
skip_dt_units = round(c_timestep/diff(t_vec(1:2)));
num_classifiers = floor(input_params.trial_length/c_timestep);


stim_resps = zeros(size(stim_resps_fullt, 1), size(stim_resps_fullt, 2), num_classifiers);

for i_tp = 1:num_classifiers
    
    % 
    stim_resps(:, :, i_tp) = squeeze(mean(stim_resps_fullt(:, :, (i_tp-1)*skip_dt_units + (1:skip_dt_units)), 3));


end


if ~(contains(network_param.network_tag, 'spnormE') || contains(network_param.network_tag, 'spEnoise'))
    % use responses from cells that are not directly stimulated
    postL1_cells = network_param.input_pattern == 0;
else
    % if spnormE is the network tag, then take all neurons as input. 
    postL1_cells = true(size((network_param.input_pattern)));
end
%% inhibitory cells logical vector
inh_cells = all(network_param.J <= 0, 1)';


% Separate the excitatory and inhibitory responses from the full stim resp
% mat. 
inh_stim_resps = stim_resps(:, inh_cells & postL1_cells, :);
exc_stim_resps = stim_resps(:, ~inh_cells & postL1_cells, :);
% Generic classifier: all post-L1 cells
stim_resps = stim_resps(:, postL1_cells, :);

% select a random subset of excitatory cells, matching the comparison
% between exc-cell and inh-cell classifiers. 
% % Modified 2/17/19: cannot randomize this step, because the
% multi-threshold analysis needs consistent exc_ISM pools. Modification:
% cells are not ordered, so take equally spaced cells over the full range.

% % old version
% exc_subset = sort(randperm(size(exc_stim_resps, 2), size(inh_stim_resps, 2)));
%  new version
exc_subset = ceil(linspace(1, size(exc_stim_resps, 2), size(inh_stim_resps, 2)));
exc_stim_resps_ISM = exc_stim_resps(:, exc_subset, :);


% Output and save
activity_struct.stim_y = stim_y;
activity_struct.stim_resps = stim_resps;
activity_struct.exc_stim_resps = exc_stim_resps;
activity_struct.exc_stim_resps_ISM = exc_stim_resps_ISM;
activity_struct.inh_stim_resps = inh_stim_resps;
activity_struct.postL1_cells = postL1_cells;
activity_struct.inh_cells = inh_cells;
activity_struct.t_vec = t_vec;
activity_struct.c_timestep = c_timestep;

stimulus_experiment_tag = [input_type 'discExp_pw' num2str(100*pulse_width)];
save(['results/' file_tag '/matfiles/stim_resps_' stimulus_experiment_tag], 'activity_struct')