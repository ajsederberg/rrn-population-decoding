function fitMultiClassifiersModelData(activity_struct, output_filename, eff_save_dir, full_save_dir)
% % Decoder fitting script - freq vs 2freq experiments (only two stims)
% % Fits four classifiers: using all cells, using all excitatory
% % cells, using a subset of excitatory cells, and using all inhibitory
% % cells. 
% 
% network_param = lookUpNetworkParam(network_code, rng_index);
% 
% matfile_dir = network_param.matsdir;
% % results_dir = '/Volumes/home/asederberg6/projects/random_network_readouts/results/';
% 
% % num_neur = network_param.numNeur;
% file_tag = network_param.network_tag; %
% 
% %%%% or set values individually, if running as a script
% % num_neur = 200;
% % gain_neur = 110;
% % rng_val = 2;
% % ei_ratio = 100; % 92 or 108;  80 or 125
% % pulse_width = 1.50; % 50 or 150 or 500
% % connect_tag = 'norm';
% % f2f_type = '';
% 
% % file_tag = [connect_tag '_g' num2str(gain_neur) '_n' num2str(num_neur) '_ei' num2str(ei_ratio) '_rng' num2str(rng_val)];
% %%%%
% 
% % output file name
% output_filename = [network_param.save_dir 'matfiles/classifierAEIS_' file_tag '_' f2f_type 'discExp_pw' num2str(100*pulse_width)];
% 
% 
% % load('/Users/audreysederberg/Dropbox/postdoctoral projects/stanley lab work/code/theory code/results/norm_g100_n100_ei80_rng1/matfiles/sim_data_norm_g100_n100_ei80_rng1_f2fdiscExp_pw150.mat')
% x_res = load([matfile_dir '/sim_data_' file_tag '_' f2f_type 'discExp_pw' num2str(100*pulse_width) '.mat']);
% 
% proc_simdata = x_res.proc_simdata;
% input_params = x_res.input_params;
% 
% network_param = x_res.network_param;
% 
% %% can compute and save the AUC and network stats
% computeSelXCJ_stats(proc_simdata, network_param, f2f_type, false);
% 
% if exist(output_filename, 'file')
%     disp([output_filename ' already exists'])
%    return;
% end
% %%
% if strcmp(input_params.ampfreq_combo, 'paired')
%     num_stims = size(proc_simdata.stimulus_info_array, 2);   
%     
%     % for the paired stim, the stim infor array has ones on the diagonal
%     % entries corresponding to what the stim ID is. 
%     trial_inds = cellfun(@(x) find(proc_simdata.stimulus_info_array(:, x, x) == 1), ...
%         num2cell(1:num_stims)', 'UniformOutput', false);
% %     
% %     % find lower-rate and higher-rate trials. 
% %     trial_inds = {find(cellfun(@length, proc_simdata.precise_impulse_times) <= 6); find(cellfun(@length, proc_simdata.precise_impulse_times) > 6)};
% %     
%     % this puts all stim 1's first followed by all stim 2's
%     stim_y = cell2mat(cellfun(@(x) x*ones(sum(proc_simdata.stimulus_info_array(:, x, x) == 1), 1), ...
%         num2cell(1:num_stims)', 'UniformOutput', false));
% end
% 
% stim_resps = cell2mat(cellfun(@(x) proc_simdata.dr_out_peristim(x, :, :), trial_inds, ...
%     'UniformOutput', false));
% % pick out a time to look at the classifier. Probably between 1 and 50
% % I_TC = 43;      % this is capitalized so that it is not used as a four-loop index
% 
% % use responses from cells that are not directly stimulated
% postL1_cells = network_param.input_pattern == 0;
% 
% % inhibitory cells logical vector
% inh_cells = all(network_param.J <= 0, 1)';
% 
% 
% % Separate the excitatory and inhibitory responses from the full stim resp
% % mat. 
% inh_stim_resps = stim_resps(:, inh_cells & postL1_cells, :);
% exc_stim_resps = stim_resps(:, ~inh_cells & postL1_cells, :);
% % Generic classifier: all post-L1 cells
% stim_resps = stim_resps(:, postL1_cells, :);
% 
% % select a random subset of excitatory cells, matching the comparison
% % between exc-cell and inh-cell classifiers. 
% exc_subset = sort(randperm(size(exc_stim_resps, 2), size(inh_stim_resps, 2)));
% exc_stim_resps_ISM = exc_stim_resps(:, exc_subset, :);
%% LOad from activity struct

stim_y = activity_struct.stim_y;
stim_resps = activity_struct.stim_resps;
exc_stim_resps = activity_struct.exc_stim_resps;
exc_stim_resps_ISM = activity_struct.exc_stim_resps_ISM;
inh_stim_resps = activity_struct.inh_stim_resps;
postL1_cells = activity_struct.postL1_cells;
inh_cells = activity_struct.inh_cells;
t_vec = activity_struct.t_vec;
c_timestep = activity_struct.c_timestep;
%% Split trials into test and train sets
train_inds = (1:2:size(stim_y, 1))';
test_inds = setdiff((1:size(stim_y, 1))', train_inds);

% c_timestep = 1; % in units of tau, how often to fit a new classifier
% skip_dt_units = round(c_timestep/diff(t_vec(1:2)));


% compute the number of classifiers to fit 
% num_classifiers = floor(input_params.trial_length/c_timestep);
num_classifiers = size(stim_resps, 3);

% Initialize classifier results arrays
% % all cells pooled together (but no L1 cells)
normClassifier = cell(num_classifiers, 1);
normValidationAccuracy = cell(num_classifiers, 1);

% % all inhibitory cells
inhClassifier = cell(num_classifiers, 1);
inhValidationAccuracy = cell(num_classifiers, 1);

% % all excitatory cells
excClassifier = cell(num_classifiers, 1);
excValidationAccuracy = cell(num_classifiers, 1);

% % subset of excitatory cells (same # as inhibitory cells)
excISMClassifier = cell(num_classifiers, 1);
excISMValidationAccuracy = cell(num_classifiers, 1);

% Set classifier options: 'linearSVM', don't apply PCA, and full_save is
% true. 
classifier_options.classifier_type = 'linearSVM'; %'pLDA';
classifier_options.applyPCA = false;
classifier_options.full_save = true;

% Time vec for classifiers 
c_time_vec = (1:num_classifiers)*c_timestep;
%% stim labels are always the same
table_y = stim_y(train_inds);
for i_tp = 1:num_classifiers
    
    % train for the post-L1 network
    table_x = stim_resps(train_inds, :, i_tp);
    data_table = table(table_x, table_y);    
    [normClassifier{i_tp}, normValidationAccuracy{i_tp}] = ...
        trainGenericClassifier(data_table, classifier_options);
    
    
    % train for the inhibitory cells only classifier
    table_x = inh_stim_resps(train_inds, :, i_tp); 
    data_table = table(table_x, table_y);    
    [inhClassifier{i_tp}, inhValidationAccuracy{i_tp}] = ...
        trainGenericClassifier(data_table, classifier_options);
    
    % train for the excitatory cells model
    table_x = exc_stim_resps(train_inds, :, i_tp);
    data_table = table(table_x, table_y);
    [excClassifier{i_tp}, excValidationAccuracy{i_tp}] = ...
        trainGenericClassifier(data_table, classifier_options);
    
    % train for the excitatory cells subset only
    table_x = exc_stim_resps_ISM(train_inds, :, i_tp);
    data_table = table(table_x, table_y);
    [excISMClassifier{i_tp}, excISMValidationAccuracy{i_tp}] = ...
        trainGenericClassifier(data_table, classifier_options);
    
   
end

%%
if strcmp(classifier_options.classifier_type, 'pLDA') % basically never use this option
    pinv_eps = 0;
    classifier_weights = cell2mat(cellfun(@(x) (pinv(x.ClassificationEnsemble.Sigma, pinv_eps)*diff(x.ClassificationEnsemble.Mu, 1, 1)')', normClassifier, 'UniformOutput', false));
else
    % classifier weights 

    % for original classifiers
    [classifier_weights, classifier_bs, classifier_betas, classifier_biases] = ...
        cellfun(@(x) linearSVMWeightsUnZscored(x.ClassificationEnsemble.BinaryLearners), ...
        normClassifier, 'uniformoutput', false);

    classifier_betas = cell2mat(classifier_betas);
    classifier_weights = cell2mat(classifier_weights);
    classifier_bs = cell2mat(classifier_bs);

    % for inh-only 
    [inh_classifier_weights, inh_classifier_bs, inh_classifier_betas, inh_classifier_biases] = ...
        cellfun(@(x) linearSVMWeightsUnZscored(x.ClassificationEnsemble.BinaryLearners), ...
        inhClassifier, 'uniformoutput', false);

    inh_classifier_betas = cell2mat(inh_classifier_betas);
    inh_classifier_weights = cell2mat(inh_classifier_weights);
    inh_classifier_bs = cell2mat(inh_classifier_bs);

    % for exc-only
    [exc_classifier_weights, exc_classifier_bs, exc_classifier_betas, exc_classifier_biases] = ...
        cellfun(@(x) linearSVMWeightsUnZscored(x.ClassificationEnsemble.BinaryLearners), ...
        excClassifier, 'uniformoutput', false);

    exc_classifier_betas = cell2mat(exc_classifier_betas);
    exc_classifier_weights = cell2mat(exc_classifier_weights);
    exc_classifier_bs = cell2mat(exc_classifier_bs);

    % for subset of exc only
    [excISM_classifier_weights, excISM_classifier_bs, excISM_classifier_betas, excISM_classifier_biases] = ...
        cellfun(@(x) linearSVMWeightsUnZscored(x.ClassificationEnsemble.BinaryLearners), ...
        excISMClassifier, 'uniformoutput', false);

    excISM_classifier_betas = cell2mat(excISM_classifier_betas);
    excISM_classifier_weights = cell2mat(excISM_classifier_weights);
    excISM_classifier_bs = cell2mat(excISM_classifier_bs);


end

%% classifier prediction accuracy using mis-matched classifiers
classIdataJ_valAcc = zeros(num_classifiers, num_classifiers);
exc_classIdataJ_valAcc = zeros(num_classifiers, num_classifiers);
inh_classIdataJ_valAcc = zeros(num_classifiers, num_classifiers);
excISM_classIdataJ_valAcc = zeros(num_classifiers, num_classifiers);

% Use test-set data now
table_y = stim_y(test_inds);

for i_tc = 1:num_classifiers
    for i_td = 1:num_classifiers
        %%
        table_x = stim_resps(test_inds, :, i_td);
        data_table = table(table_x, table_y);
        pred_y = normClassifier{i_tc}.predictFcn(data_table);
        classIdataJ_valAcc(i_tc, i_td) = mean(pred_y == stim_y(test_inds)); 
        
        
        % first for excitation (all) 
        table_x = exc_stim_resps(test_inds, :, i_td);
        data_table = table(table_x, table_y);
        pred_y = excClassifier{i_tc}.predictFcn(data_table);
        exc_classIdataJ_valAcc(i_tc, i_td) = mean(pred_y == stim_y(test_inds)); 

        % next for inhibition (all)
        table_x = inh_stim_resps(test_inds, :, i_td);
        data_table = table(table_x, table_y);
        pred_y = inhClassifier{i_tc}.predictFcn(data_table);
        inh_classIdataJ_valAcc(i_tc, i_td) = mean(pred_y == stim_y(test_inds)); 
    
        % finally for the subset of excitation 
        table_x = exc_stim_resps_ISM(test_inds, :, i_td);
        data_table = table(table_x, table_y);
        pred_y = excISMClassifier{i_tc}.predictFcn(data_table);
        excISM_classIdataJ_valAcc(i_tc, i_td) = mean(pred_y == stim_y(test_inds)); 
    end
end
%% Estimate errorbars at tau = 50. 

% V1 (as in the paper) : error on test set for a single classifier,
% repeated 50 times (fitting new classifiers)

% V2 (alt) : error on test set for a collection of classifiers, by
% drawing different subsets of input cells. 
%%
tp_sel = 50;    % only repeate CV draws for the last time point (50)

num_cv_draws = 50;

rep_inds_arr = cell(num_cv_draws, 1);

rep_classifier_options = classifier_options;
rep_classifier_options.full_save = false;

% initialize variables: save classifiers, validation accuracy, and testset
% accuracy
rep_normClassifier = cell(num_cv_draws, 1);
rep_normValidationAccuracy = cell(num_cv_draws, 1);
rep_normTestAccuracy = zeros(num_cv_draws, 1);

rep_excClassifier = cell(num_cv_draws, 1);
rep_excValidationAccuracy = cell(num_cv_draws, 1);
rep_excTestAccuracy = zeros(num_cv_draws, 1);

rep_excISMClassifier = cell(num_cv_draws, 1);
rep_excISMValidationAccuracy = cell(num_cv_draws, 1);
rep_excISMTestAccuracy = zeros(num_cv_draws, 1);

rep_inhClassifier = cell(num_cv_draws, 1);
rep_inhValidationAccuracy = cell(num_cv_draws, 1);
rep_inhTestAccuracy = zeros(num_cv_draws, 1);


for i_draw = 1:num_cv_draws
    % new train/test split
    rep_train_inds = sort(randperm(length(stim_y), length(train_inds)));
    rep_test_inds = setdiff(1:length(stim_y), rep_train_inds);
    
    rep_inds.train = rep_train_inds;
    rep_inds.test = rep_test_inds;
    rep_inds_arr{i_draw} = rep_inds;
    
    % fit new classifier for each cv draw
    table_y = stim_y(rep_train_inds);
    
    % train for the post-L1 network
    table_x = stim_resps(rep_train_inds, :, tp_sel);
    data_table = table(table_x, table_y);    
    [rep_normClassifier{i_draw}, rep_normValidationAccuracy{i_draw}] = ...
        trainGenericClassifier(data_table, rep_classifier_options);
    % test set accuracy
    table_x = stim_resps(rep_test_inds, :, tp_sel);
    data_table = table(table_x, table_y);
    pred_y = rep_normClassifier{i_draw}.predictFcn(data_table);
    rep_normTestAccuracy(i_draw) = mean(pred_y == stim_y(rep_test_inds)); 

    % train for the inhibitory cells only classifier
    table_x = inh_stim_resps(rep_train_inds, :, tp_sel); 
    data_table = table(table_x, table_y);    
    [rep_inhClassifier{i_draw}, rep_inhValidationAccuracy{i_draw}] = ...
        trainGenericClassifier(data_table, rep_classifier_options);
    % test accuracy for inhibition (all)
    table_x = inh_stim_resps(rep_test_inds, :, tp_sel);
    data_table = table(table_x, table_y);
    pred_y = rep_inhClassifier{i_draw}.predictFcn(data_table);
    rep_inhTestAccuracy(i_draw) = mean(pred_y == stim_y(rep_test_inds)); 

    % train for the excitatory cells model
    table_x = exc_stim_resps(rep_train_inds, :, tp_sel);
    data_table = table(table_x, table_y);
    [rep_excClassifier{i_draw}, rep_excValidationAccuracy{i_draw}] = ...
        trainGenericClassifier(data_table, rep_classifier_options);
    % test accuracy for excitation (all) 
    table_x = exc_stim_resps(rep_test_inds, :, tp_sel);
    data_table = table(table_x, table_y);
    pred_y = rep_excClassifier{i_draw}.predictFcn(data_table);
    rep_excTestAccuracy(i_draw) = mean(pred_y == stim_y(rep_test_inds));

    % train for the excitatory cells subset only
    table_x = exc_stim_resps_ISM(rep_train_inds, :, tp_sel);
    data_table = table(table_x, table_y);
    [rep_excISMClassifier{i_draw}, rep_excISMValidationAccuracy{i_draw}] = ...
        trainGenericClassifier(data_table, rep_classifier_options);
    % test accuracy for excitation (subset) 
    table_x = exc_stim_resps_ISM(rep_test_inds, :, tp_sel);
    data_table = table(table_x, table_y);
    pred_y = rep_excISMClassifier{i_draw}.predictFcn(data_table);
    rep_excISMTestAccuracy(i_draw) = mean(pred_y == stim_y(rep_test_inds));
    
   
end

%% fit classifier over different draws of cells
cell_set_sizes = [10 50 100 200 300 400];
popSizeRep_normClassifiers = cell(num_cv_draws, length(cell_set_sizes));
popSizeRep_normValidationAccuracy = cell(num_cv_draws, length(cell_set_sizes));
popSizeRep_normTestAccuracy = zeros(num_cv_draws, length(cell_set_sizes));

cell_set_draws = cell(num_cv_draws, length(cell_set_sizes));

% fit new classifier for each cv draw
table_y = stim_y(train_inds);
for i_draw = 1:num_cv_draws
    for i_setsize = 1:length(cell_set_sizes)
        
    cell_inds = sort(randperm(size(stim_resps, 2), cell_set_sizes(i_setsize)));
    cell_set_draws{i_draw, i_setsize} = cell_inds;
    %% V2: one test set, many draws of cells from full population
        % train for the excitatory cells subset only
        table_x = stim_resps(train_inds, cell_inds, tp_sel);
        data_table = table(table_x, table_y);
        [popSizeRep_normClassifiers{i_draw, i_setsize}, popSizeRep_normValidationAccuracy{i_draw, i_setsize}] = ...
            trainGenericClassifier(data_table, rep_classifier_options);
        % test accuracy for excitation (subset) 
        table_x = stim_resps(test_inds, cell_inds, tp_sel);
        data_table = table(table_x, table_y);
        pred_y = popSizeRep_normClassifiers{i_draw, i_setsize}.predictFcn(data_table);
        popSizeRep_normTestAccuracy(i_draw, i_setsize) = mean(pred_y == stim_y(test_inds));


    end
end

%% 

% Group results for easier look-up
cellGroup_string = {'all', 'exc', 'excISM', 'inh'};

cls_res(length(cellGroup_string)) = struct();

% % determined (2/17/19) that 'full_cls' will never be used and it will
% always be easier to just re-run the analysis. 
% full_cls(length(cellGroup_string)) = struct();
% all cell results
cls_res(1).cellGroup = cellGroup_string{1};
% full_cls(1).classifier = normClassifier;    % put these separate - large variable
% full_cls(1).validationAccuracy = normValidationAccuracy;
cls_res(1).classIdataJ_valAcc = classIdataJ_valAcc;
cls_res(1).classifier_weights = classifier_weights;
cls_res(1).classifier_bs = classifier_bs;
cls_res(1).classifier_betas = classifier_betas;
cls_res(1).classifier_biases = classifier_biases;
% exc cell results
cls_res(2).cellGroup = cellGroup_string{2};
% full_cls(2).classifier = excClassifier; % put these separate - large variable
% full_cls(2).validationAccuracy = excValidationAccuracy;
cls_res(2).classIdataJ_valAcc = exc_classIdataJ_valAcc;
cls_res(2).classifier_weights = exc_classifier_weights;
cls_res(2).classifier_bs = exc_classifier_bs;
cls_res(2).classifier_betas = exc_classifier_betas;
cls_res(2).classifier_biases = exc_classifier_biases;
% exc-ISM cell results
cls_res(3).cellGroup = cellGroup_string{3};
% full_cls(3).classifier = excISMClassifier;  % put these separate - large variable
% full_cls(3).validationAccuracy = excISMValidationAccuracy;
cls_res(3).classIdataJ_valAcc = excISM_classIdataJ_valAcc;
cls_res(3).classifier_weights = excISM_classifier_weights;
cls_res(3).classifier_bs = excISM_classifier_bs;
cls_res(3).classifier_betas = excISM_classifier_betas;
cls_res(3).classifier_biases = excISM_classifier_biases;
% inh cell results
cls_res(4).cellGroup = cellGroup_string{4};
% full_cls(4).classifier = inhClassifier;     % put these separate - large variables
% full_cls(4).validationAccuracy = inhValidationAccuracy;
cls_res(4).classIdataJ_valAcc = inh_classIdataJ_valAcc;
cls_res(4).classifier_weights = inh_classifier_weights;
cls_res(4).classifier_bs = inh_classifier_bs;
cls_res(4).classifier_betas = inh_classifier_betas;
cls_res(4).classifier_biases = inh_classifier_biases;

saveInParfor([eff_save_dir output_filename], t_vec, ...
    train_inds, test_inds, classifier_options, c_time_vec, ...
    cls_res)
saveInParfor([full_save_dir output_filename], t_vec, ...
    train_inds, test_inds, activity_struct, classifier_options, c_time_vec, ...
    cls_res)