function cls_res_struct = fitMultiClassifiers_helper(activity_struct, cellGroup_string, classifier_options)
% Fits classifiers and tests on reserved set for the specified cell set. 
% % Decoder fitting script - freq vs 2freq experiments (only two stims)
% % Fits four classifiers: using all cells, using all excitatory
% % cells, using a subset of excitatory cells, and using all inhibitory
% % cells. 

%% Load from activity struct

stim_y = activity_struct.stim_y;

if strcmp(cellGroup_string, 'all')
    stim_resps = activity_struct.stim_resps;
    cellInds = 1:size(stim_resps, 2);
elseif strcmp(cellGroup_string, 'exc')
    stim_resps = activity_struct.exc_stim_resps;
    cellInds = 1:size(stim_resps, 2);
elseif strcmp(cellGroup_string, 'excISM')
	stim_resps = activity_struct.exc_stim_resps_ISM;
    cellInds = [];
elseif strcmp(cellGroup_string, 'inh')
	stim_resps = activity_struct.inh_stim_resps;
    cellInds = [];
elseif strcmp(cellGroup_string, 'exc50')
    % choose 50 of the excitatory cells (connectivity random, should not
    % matter which 50)
    cellInds = 5*(1:50);
    stim_resps = activity_struct.exc_stim_resps(:, cellInds, :);
elseif strcmp(cellGroup_string, 'inh50')
    % choose 50 of the inhibitory cells (connectivity random, should not
    % matter which 50)
    cellInds = 2*(1 : 50);
    stim_resps = activity_struct.inh_stim_resps(:, cellInds, :);    
elseif strcmp(cellGroup_string, 'exc50b')
    % choose 50 of the excitatory cells (connectivity random, should not
    % matter which 50)
    cellInds = -2 + 5*(1:50);
    stim_resps = activity_struct.exc_stim_resps(:, cellInds, :);
elseif strcmp(cellGroup_string, 'inh50b')
    % choose 50 of the inhibitory cells (connectivity random, should not
    % matter which 50)
    cellInds = -1 + 2*(1 : 50);
    stim_resps = activity_struct.inh_stim_resps(:, cellInds, :);    
end
% postL1_cells = activity_struct.postL1_cells;
% inh_cells = activity_struct.inh_cells;
% t_vec = activity_struct.t_vec;
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

% Time vec for classifiers 
c_time_vec = (1:num_classifiers)*c_timestep;
%% stim labels are always the same
table_y = stim_y(train_inds);

% option to run only on certain time points.
if isfield(classifier_options, 'timepoint_subset')
    timepoint_subset = classifier_options.timepoint_subset;
else
    timepoint_subset = 1:num_classifiers;
end
for i_tp = timepoint_subset
    
    % train for the post-L1 network
    table_x = stim_resps(train_inds, :, i_tp);
    data_table = table(table_x, table_y);    
    [normClassifier{i_tp}, normValidationAccuracy{i_tp}] = ...
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

end

%% classifier prediction accuracy using mis-matched classifiers
classIdataJ_valAcc = zeros(num_classifiers, num_classifiers);

% Use test-set data now
table_y = stim_y(test_inds);

% recall that timepoint_subset selects a subset of classifiers time points
% (not all 50). 
for i_tc = timepoint_subset
    for i_td = timepoint_subset
        %%
        table_x = stim_resps(test_inds, :, i_td);
        data_table = table(table_x, table_y);
        pred_y = normClassifier{i_tc}.predictFcn(data_table);
        classIdataJ_valAcc(i_tc, i_td) = mean(pred_y == stim_y(test_inds)); 
        
       
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
% rep_classifier_options.full_save = false;

% initialize variables: save classifiers, validation accuracy, and testset
% accuracy
rep_normClassifier = cell(num_cv_draws, 1);
rep_normValidationAccuracy = cell(num_cv_draws, 1);
rep_normTestAccuracy = zeros(num_cv_draws, 1);

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
    [rep_normClassifier{i_draw}, full_valAcc] = ...
        trainGenericClassifier(data_table, rep_classifier_options);
    % partitionedModel is huge and never used. 
    rep_normValidationAccuracy{i_draw} = rmfield(full_valAcc, 'partitionedModel');
    
    % test set accuracy
    table_x = stim_resps(rep_test_inds, :, tp_sel);
    data_table = table(table_x, table_y);
    pred_y = rep_normClassifier{i_draw}.predictFcn(data_table);
    rep_normTestAccuracy(i_draw) = mean(pred_y == stim_y(rep_test_inds)); 



end

% extract weights
[rep_classifier_weights, rep_classifier_bs, rep_classifier_betas, rep_classifier_biases] = ...
    cellfun(@(x) linearSVMWeightsUnZscored(x.ClassificationEnsemble.BinaryLearners), ...
    rep_normClassifier, 'uniformoutput', false);

rep_classifier_betas = cell2mat(cellfun(@(x) shiftdim(x, -1), rep_classifier_betas, 'uniformoutput', false));
rep_classifier_weights = cell2mat(cellfun(@(x) shiftdim(x, -1), rep_classifier_weights, 'uniformoutput', false));
rep_classifier_bs = cell2mat(rep_classifier_bs);
%% fit classifier over different draws of cells
cell_set_sizes = ceil((0.1:0.2:1)*size(stim_resps, 2));
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
        % train for the new subset
        table_x = stim_resps(train_inds, cell_inds, tp_sel);
        data_table = table(table_x, table_y);
        [popSizeRep_normClassifiers{i_draw, i_setsize}, full_valAcc] = ...
            trainGenericClassifier(data_table, rep_classifier_options);
         popSizeRep_normValidationAccuracy{i_draw, i_setsize} = rmfield(full_valAcc, 'partitionedModel');
        
        % test accuracy for subset
        table_x = stim_resps(test_inds, cell_inds, tp_sel);
        data_table = table(table_x, table_y);
        pred_y = popSizeRep_normClassifiers{i_draw, i_setsize}.predictFcn(data_table);
        popSizeRep_normTestAccuracy(i_draw, i_setsize) = mean(pred_y == stim_y(test_inds));


    end
end

%% 

% % Group results for easier look-up
% cellGroup_string = {'all', 'exc', 'excISM', 'inh'};

% all cell results
cls_res_struct.cellGroup = cellGroup_string;
cls_res_struct.cellInds = cellInds;
cls_res_struct.classIdataJ_valAcc = classIdataJ_valAcc;
cls_res_struct.classifier_weights = classifier_weights;
cls_res_struct.classifier_bs = classifier_bs;
cls_res_struct.classifier_betas = classifier_betas;
cls_res_struct.classifier_biases = classifier_biases;
cls_res_struct.c_time_vec = c_time_vec;
cls_res_struct.test_inds = test_inds;
cls_res_struct.train_inds = train_inds;

cls_res_struct.cvDraws_valAcc = rep_normValidationAccuracy;
cls_res_struct.cvDraws_testAcc = rep_normTestAccuracy;
cls_res_struct.cvDraws_testtrain = rep_inds_arr;
cls_res_struct.cvclassifier_weights = rep_classifier_weights;
cls_res_struct.cvclassifier_bs = rep_classifier_bs;
cls_res_struct.cvclassifier_betas = rep_classifier_betas;
cls_res_struct.cvclassifier_biases = rep_classifier_biases;

cls_res_struct.cvSetDraws_valAcc = popSizeRep_normValidationAccuracy;
cls_res_struct.cvSetDraws_testAcc = popSizeRep_normTestAccuracy;
cls_res_struct.cvSetDraws_inds = cell_set_draws;

