% LAST SCRIPT : use many frequencies, pick a threshold. 
network_gain = 4;
num_neurons = 500;
gain_neur = network_gain*100;
% rng_val = 3;
ei_ratio = 100; % 92 or 108;  80 or 125
pulse_width = .50; % 50 or 150 or 500
connect_tag = 'spnorm'; % 'spnorm2P' for two-pool
file_tag = [connect_tag '_g' num2str(gain_neur) '_n' ...
    num2str(num_neurons) '_ei' num2str(ei_ratio) '_rng'];

%%
rng_val = 1;
f2f_types = {'f714f', 'f2f', 'f1020f'};
as_all = cell(length(f2f_types), 1);

for i_ft = 1:length(f2f_types)
    input_type = [f2f_types{i_ft} 'discExp_pw' num2str(100*pulse_width)];
    x = load(['results/' file_tag num2str(rng_val) '/matfiles/stim_resps_' ...
        input_type]);
    as_all{i_ft} = x.activity_struct;
    
end
%% frequencies for {'f714f', 'f2f', 'f1020f'}; are {'7, 14', '8, 16', ; '10, 20'}
v1_task = {[7 8], [10 14]};
v2_task = {[8 10], [14 16]};
v3_task = {[10 14], [16 20]};

stim_freqs = {[7 14], [8 16], [10 20]};
% act_str_1
% idea: label stims by frequency
for i_ft = 1:length(as_all)
    stim_f = zeros(size(as_all{i_ft}.stim_y));
    stim_f(as_all{i_ft}.stim_y == 1) = stim_freqs{i_ft}(1);
    stim_f(as_all{i_ft}.stim_y == 2) = stim_freqs{i_ft}(2);
    as_all{i_ft}.stim_f = stim_f;
end

%% what follows is CRAPPY code with magic numbers that probably has mistakes in it
% and needs to be generated more automatically. 
field_names = {'stim_resps', 'exc_stim_resps', 'exc_stim_resps_ISM', ...
    'inh_stim_resps'};
v1_as = as_all{i_ft};
for i_fn = 1:length(field_names)
    v1_stim1_resps = [as_all{1}.(field_names{i_fn})(as_all{1}.stim_y == 1, :, :); ... 
        as_all{2}.(field_names{i_fn})(as_all{2}.stim_y == 1, :, :)];
    v1_stim2_resps = [as_all{1}.(field_names{i_fn})(as_all{1}.stim_y == 2, :, :); ... 
        as_all{3}.(field_names{i_fn})(as_all{3}.stim_y == 1, :, :)];

    v1_as.(field_names{i_fn}) = [v1_stim1_resps; v1_stim2_resps];
end
v1_as.stim_y = [ones(size(v1_stim1_resps, 1), 1); ...
    2*ones(size(v1_stim2_resps, 1), 1)];

v2_as = as_all{i_ft};
for i_fn = 1:length(field_names)
    v2_stim1_resps = [as_all{2}.(field_names{i_fn})(as_all{2}.stim_y == 1, :, :); ... 
        as_all{3}.(field_names{i_fn})(as_all{3}.stim_y == 1, :, :)];
    v2_stim2_resps = [as_all{1}.(field_names{i_fn})(as_all{1}.stim_y == 2, :, :); ... 
        as_all{2}.(field_names{i_fn})(as_all{2}.stim_y == 2, :, :)];

    v2_as.(field_names{i_fn}) = [v2_stim1_resps; v2_stim2_resps];
end
v2_as.stim_y = [ones(size(v2_stim1_resps, 1), 1); ...
    2*ones(size(v2_stim2_resps, 1), 1)];
%% fit v1 anvd v2 classifiers. 
classifier_options.classifier_type = 'linearSVM'; %'pLDA';
classifier_options.applyPCA = false;
classifier_options.full_save = true;
% cycle through cell groups
cellGroup_string = {'all', 'exc', 'excISM', 'inh'};
%% run from here
tic
v1_cls_res_struct = fitMultiClassifiers_helper(v1_as, 'all', classifier_options);
toc
v2_cls_res_struct = fitMultiClassifiers_helper(v2_as, 'all', classifier_options);
toc
%%
i_TC = 48;
wx_v1 = v1_cls_res_struct.classifier_weights(i_TC, :)*squeeze(v1_as.stim_resps(:, :, i_TC))';

f_uniq = unique(v1_as.stim_f);
for i_f = 1:length(f_uniq)
wx_vs_f(i_f) = mean(wx_v1(v1_as.stim_f == f_uniq(i_f)));
end