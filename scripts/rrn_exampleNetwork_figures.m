% Make figure set from a single run. 


% num_neur = network_param.numNeur;
% file_tag = network_param.network_tag; %

%%%% or set values individually, if running as a script
num_neur = 200;
gain_neur = 200;
rng_val = 1;
ei_ratio = 100; % 92 or 108;  80 or 125
pulse_width = 1.50; % 50 or 150 or 500
connect_tag = 'spnorm2P';
f2f_type = '';

file_tag = [connect_tag '_g' num2str(gain_neur) '_n' num2str(num_neur) '_ei' num2str(ei_ratio) '_rng' num2str(rng_val)];


% results_dir = '/Volumes/home/asederberg6/projects/random_network_readouts/results/';
results_dir = 'temp/';


classifier_file_tag = ['classifier_' file_tag '_f2f' f2f_type 'discExp_pw' num2str(100*pulse_width)];

x_res = load([results_dir file_tag '/matfiles/sim_data_' file_tag '_f2f' f2f_type 'discExp_pw' num2str(100*pulse_width) '.mat']);

proc_simdata = x_res.proc_simdata;
input_params = x_res.input_params;

network_param = x_res.network_param;

x_class = load([network_param.matsdir classifier_file_tag]);
%% load necessary variables from classifier results

t_vec = proc_simdata.peri_stim_time;


stim_resps = x_class.stim_resps;
stim_y = x_class.stim_y;
classifier_options = x_class.classifier_options; 
c_timestep = x_class.c_timestep;
train_inds = x_class.train_inds;
test_inds = x_class.test_inds;
skip_dt_units = x_class.skip_dt_units;

allClassifier = x_class.normClassifier;
allValidationAccuracy = x_class.normValidationAccuracy;
classIdataJ_valAcc = x_class.classIdataJ_valAcc;
%%
% postL1 cells are not direct recipinets of the input current
postL1_cells = network_param.input_pattern == 0;
l1_inds = find(~postL1_cells);
l2_inds = find(postL1_cells);

% inhibitory cells logical vector
inh_cells = all(network_param.J <= 0, 1);

% all inhibitory cells are post-L1; this gets the logical vector of where
% inhibitory cells are once the post-L1 cells have been separated from the
% full network. 
postL1_inh_cells = inh_cells(postL1_cells);

% pick out a time to look at the classifier. Probably between 1 and 50
I_TC = 43;      % this is capitalized so that it is not used as a four-loop index

% Separate the excitatory and inhibitory responses from the full stim resp
% mat. 
inh_stim_resps = stim_resps(:, postL1_inh_cells, :);
exc_stim_resps = stim_resps(:, ~postL1_inh_cells, :);

% select a random subset of excitatory cells, matching the comparison
% between exc-cell and inh-cell classifiers. 
exc_subset = ceil(linspace(1, sum(~postL1_inh_cells), sum(postL1_inh_cells)));
exc_stim_resps_ISM = exc_stim_resps(:, exc_subset, :);

%% input current analysis and plots

%% Methods plots
offset_time = -10;
exp_length = input_params.isi*input_params.num_reps*2;
% plot examples of the input functions. sample_freq = 4 is about as high as
% this can go without crashing the computer; try serially computting
% in_vals instead of on the entire t_vals vector 
sample_freq = 4;
t_vals = input_params.isi + input_params.burn_in_time + offset_time + (1/sample_freq:1/sample_freq:exp_length);
% in_vals = input_params.input_fun(t_vals');
t_mat = reshape(t_vals, sample_freq*input_params.isi, 2*input_params.num_reps);
in_mat = zeros(sample_freq*input_params.isi, 2*input_params.num_reps);

for i_rep = 1:size(t_mat, 2)
    in_mat(:, i_rep) = input_params.input_fun(t_mat(:, i_rep))';
end
%%

netwk_hz = 50;      % 1/network_tau = 50 Hz

% these come directly from proc_simdata: stim_IDs is the correct stimulus
% label
stim_IDs = input_params.input_ampfreq(:, 2) == input_params.input_frequency{1};
in1_freq = input_params.input_frequency{1};
in2_freq = input_params.input_frequency{2};
in1_mat = imfilter(in_mat(:,  stim_IDs), fspecial('average', [10 1]));
in2_mat = imfilter(in_mat(:, ~stim_IDs), fspecial('average', [10 1]));

% these come from 'stim_resps', which has been sorted. use 'stim_y', which
% is sorted the same way. 
% average network-wide firing rates
aveS1_netfr = mean(mean(stim_resps(stim_y==1, :, :), 3), 2);
aveS2_netfr = mean(mean(stim_resps(stim_y==2, :, :), 3), 2);


% fourier transform analysis of input current 
pulse_fun = @(t, a) (t.^2/a^2).*exp((t>0).*(1-t/a)).*(t > 0);

t_trial = 1/sample_freq:1/sample_freq:input_params.isi;
stim_on_time = t_trial > -offset_time & t_trial <= input_params.trial_length-offset_time;
% pulse_fft = fft(pulse_fun(t_trial(stim_on_time), pulse_width));
% 
% in1_fft = fft(in1_mat(stim_on_time, :));
% in1_asd = abs(bsxfun(@(x,y) x./y, in1_fft(1:(end/2+1), :), pulse_fft(1:(end/2+1))'));
% in1_asd(2:end-1, :) = in1_asd(2:end-1, :)*2;
% 
% in2_fft = fft(in2_mat(stim_on_time, 2:end));
% in2_asd = abs(bsxfun(@(x,y) x./y, in2_fft(1:(end/2+1), :), pulse_fft(1:(end/2+1))'));
% in2_asd(2:end-1, :) = in2_asd(2:end-1, :)*2;
% nu_vals = (0:size(in1_asd, 1)-1)*sample_freq*netwk_hz/2/size(in1_asd, 1);

% binary input mat: +1 is an impulse time on a low-freq trial, -1 is an
% impulse time on a high-freq trial
round_fac = 5;
in1_times = round(round_fac*(input_params.input_times{1} - input_params.burn_in_time - input_params.isi + 10));
in2_times = round(round_fac*(input_params.input_times{2} - input_params.burn_in_time - input_params.isi + 10));

binary_inputs = sparse(1, in1_times, 1, 1, round_fac*input_params.isi*input_params.num_reps*2) - ...
    sparse(1, in2_times, 1, 1, round_fac*input_params.isi*input_params.num_reps*2);
binary_mat = reshape(binary_inputs, input_params.isi*round_fac, 2*input_params.num_reps);

bin1_mat = binary_mat(:, stim_IDs);
bin2_mat = binary_mat(:, ~stim_IDs);

makeMyFigure(30, 20);

sh_in1 = subplot(4, 3, 7);
sh_in2 = subplot(4, 3, 10);
sh_inPulse = subplot(3, 3, 3);
sh_Jrep = subplot(2, 3, 1);
sh_Jdist = subplot(3, 6, 4);

% fft_in12 = subplot(2, 4, 5);
netfr_in12 = subplot(2, 3, 6);
% sh_wtrd = subplot(2, 2, 2);
% sh_valAc = subplot(2, 2, 4);

sh_fr1 = subplot(4, 3, 8);
sh_fr2 = subplot(4, 3, 11);

stim_colors = lines(4);
stim_colors = stim_colors(3:4, :);

% Plot input stuff
set(gcf, 'CurrentAxes', sh_inPulse);
ex_pulse_widths = [0.5 1.5];
line_styles = {'-', '-.', '--'};
line_widths = 0.5 + 0*ex_pulse_widths;
line_widths(ex_pulse_widths == pulse_width) = 1.5;
hold on
for i_pw = 1:length(ex_pulse_widths)
    plot(t_vec, pulse_fun(t_vec, ex_pulse_widths(i_pw)), 'k',...
        'linestyle', line_styles{i_pw}, 'linewidth', line_widths(i_pw))
end
xlim([0 10])
set(gca, 'color', 'none', 'ytick', [])
lh = legend(num2str(ex_pulse_widths', 'a = %1.1f'), 'Position', [0.8400 0.8695 0.0812 0.0467]);
axis square
xlabel('time')
ylabel('input pulse waveform')
% title(func2str(input_params.pulse_fun))
title('(t/a)^2 exp(1 - t/a) for t > 0; else 0')

set(gcf, 'CurrentAxes', sh_in1);
colormap(sh_in1, flipud(gray(2)));
% hold on
% plot(t_trial, in1_mat(:, ex_trs), 'k')
imagesc(t_trial, [], bin1_mat', [0 1])
title(['input frequency: ' num2str(netwk_hz*in1_freq)])
xlim([0 20 + input_params.trial_length])
ylabel('trials')
xlabel('time')
yyaxis right
mPMsd = mean(in1_mat, 2) + 0*std(in1_mat, [], 2)*[-1 0 1];
plot(t_trial, mPMsd, '-', 'linewidth', 2, 'color', stim_colors(1, :))
ylim([0 max(mPMsd(:))])
set(gca, 'YColor', stim_colors(1, :));

v = axis();

set(gcf, 'CurrentAxes', sh_in2);
colormap(sh_in2, flipud(gray(2)));
% hold on
% plot(t_trial, in2_mat(:, ex_trs), 'k')
imagesc(t_trial, [], -bin2_mat', [0 1])
title(['input frequency: ' num2str(netwk_hz*in2_freq)])
xlim([0 20 + input_params.trial_length])
ylabel('trials')
xlabel('time')

yyaxis right
mPMsd = mean(in2_mat, 2) + 0*std(in2_mat, [], 2)*[-1 0 1];

plot(t_trial, mPMsd, '-', 'linewidth', 2, 'color', stim_colors(2, :))
ylim([0 max(v(4), max(mPMsd(:)))])
ylabel('ave input current over 10 \tau +/- sd across trials')
set(gca, 'YColor', stim_colors(2, :));
set(get(gca, 'YLabel'), 'Position', [80 4 0], 'Rotation', -90, 'color', [0 0 0]);




set(gcf, 'CurrentAxes', netfr_in12);
hold on
histogram(aveS1_netfr, 'FaceColor', stim_colors(1, :))
histogram(aveS2_netfr, 'FaceColor', stim_colors(2, :))
xlabel('network average activation (fr) over time and units')
ylabel('trial counts')
axis square
set(gca, 'color', 'none')
legend({num2str(netwk_hz*in1_freq, 'freq : %1.0f'),...
    num2str(netwk_hz*in2_freq, 'freq : %1.0f')}, 'Location', 'northwest')


set(gcf, 'CurrentAxes', sh_Jrep);
neuron_colors = [(1 + lines(1))/2; lines(2)];

if strcmp(network_param.network_tag(1:8), 'spnorm2P')
%     theta_all = sort(2*pi*rand(num_neur, 1)) - pi/3;
    
    nE = network_param.numNeurExc;
    nI = network_param.numNeurInh;
    % pool order is E1, I1, E2, I2
%     pool_locs = [1:nE/2 nE + (1:nI/2) nE/2 + (1:nE/2) nE + nI/2 + (1:nI/2)];
%     [~, pool_order] = sort(pool_locs);
%     theta_all = theta_all(pool_order);
    % center these at 0, pi/2, pi and 3pi/2
    theta_all = zeros(nE + nI, 1);
    dth = pi/9;
    theta_all(1:nE/2) = dth*randn(nE/2, 1);
    theta_all( nE + (1:nI/2) ) = pi/2 + dth*randn(nI/2, 1);
    theta_all(nE/2 + (1:nE/2)) = pi + dth*randn(nE/2, 1);
    theta_all(nE + nI/2 + (1:nI/2)) = 3*pi/2 + dth*randn(nI/2, 1);

    r_all = sqrt(0.5 + 0.5*rand(num_neur, 1));

else
    theta_all = 2*pi*rand(num_neur, 1);
    r_all = sqrt(1.1*rand(num_neur, 1));

end
x_pos = r_all.*cos(theta_all);
y_pos = r_all.*sin(theta_all);
theta_left = pi + 0.3*linspace(-1, 1, length(l1_inds));
% x_pos(l1_inds) = 1.5*cos(theta_left);
% y_pos(l1_inds) = 1.5*sin(theta_left);


hold on

p_draw = 0.3;
for i_neur = 1:num_neur
    connected_cells = find(network_param.J(:, i_neur) ~= 0 & rand(num_neur, 1) < p_draw);
    
        x_vals = [x_pos(i_neur) + 0*connected_cells x_pos(connected_cells)]';
        y_vals = [y_pos(i_neur) + 0*connected_cells y_pos(connected_cells)]';
        if inh_cells(i_neur)
            plot(x_vals, y_vals, 'color', neuron_colors(3, :))
        else
            plot(x_vals, y_vals, 'color', neuron_colors(2, :))
        end
end
plot(x_pos(l1_inds), y_pos(l1_inds), 'ko', 'MarkerFaceColor', neuron_colors(1, :));
plot(x_pos(~inh_cells' & postL1_cells), y_pos(~inh_cells' & postL1_cells), 'ko', 'MarkerFaceColor', neuron_colors(2, :));
plot(x_pos(inh_cells), y_pos(inh_cells), 'ko', 'MarkerFaceColor', neuron_colors(3, :));
axis(1.1*[-1 1 -1 1])
axis square
set(gca, 'color', 'none', 'visible', 'off')
th = text(-2, 1.3,  ['network connectivity (' num2str(p_draw*100, '%1.0f') '% of connections drawn)'], ...
    'fontweight', 'bold'); 

set(gcf, 'CurrentAxes', sh_Jdist);
hold on
histogram(network_param.J(network_param.J > 0), 'EdgeColor', 'none')
histogram(network_param.J(network_param.J < 0), 'EdgeColor', 'none')
set(gca, 'xtick', [-network_param.gII; 0; network_param.gEE ]/sqrt(num_neur), 'XTickLabel', {'g_I','0','g_E'})
set(gca, 'color', 'none')
axis tight
axis square
xlabel('non-zero weights')
title(['sparsity of connections: ' num2str(network_param.sparsity)])


ex_ind1 = 2;
ex_freq1 = length(proc_simdata.precise_impulse_times{ex_ind1});
ex_ind2 = 800;
ex_freq2 = length(proc_simdata.precise_impulse_times{ex_ind2});

disp(['freq 1 : ' num2str(ex_freq1) ' and freq 2 : ' num2str(ex_freq2)])

ave_inh_L2_fr1 = squeeze(mean(proc_simdata.dr_out_peristim(ex_ind1, inh_cells, :), 2));
ave_exc_L1_fr1 = squeeze(mean(proc_simdata.dr_out_peristim(ex_ind1, ~postL1_cells, :), 2));
ave_exc_L2_fr1 = squeeze(mean(proc_simdata.dr_out_peristim(ex_ind1, postL1_cells & ~inh_cells', :), 2));
impulsetimes_ex1 = cell2mat(proc_simdata.precise_impuse_peritimes{ex_ind1});

ave_inh_L2_fr2 = squeeze(mean(proc_simdata.dr_out_peristim(ex_ind2, inh_cells, :), 2));
ave_exc_L1_fr2 = squeeze(mean(proc_simdata.dr_out_peristim(ex_ind2, ~postL1_cells, :), 2));
ave_exc_L2_fr2 = squeeze(mean(proc_simdata.dr_out_peristim(ex_ind2, postL1_cells & ~inh_cells', :), 2));
impulsetimes_ex2 = cell2mat(proc_simdata.precise_impuse_peritimes{ex_ind2});


% plot single-trial examples of exc/inh firing rate
set(gcf, 'CurrentAxes', sh_fr1);
ph = plot(t_vec,[ave_exc_L1_fr1  ave_exc_L2_fr1  ave_inh_L2_fr1], 'linewidth', 1.5); 
assignColorsToLines(ph, neuron_colors);
title(['stimulus freq ' num2str(ex_freq1) ' single trial, pop. ave.'], 'Color', stim_colors(stim_y(ex_ind1), :))
set(gca, 'color', 'none')

set(gcf, 'CurrentAxes', sh_fr2);
ph = plot(t_vec,[ave_exc_L1_fr2  ave_exc_L2_fr2  ave_inh_L2_fr2], 'linewidth', 1.5); 
assignColorsToLines(ph, neuron_colors);
title(['stimulus freq ' num2str(ex_freq2) ' single trial, pop. ave.'], 'Color', stim_colors(stim_y(ex_ind2), :))
set(gca, 'color', 'none')



print(gcf, '-dpdf', [network_param.plotdir classifier_options.classifier_type 'MethodsClear_' network_param.network_tag 'pw' num2str(100*pulse_width)]);


% % % %% extra analysis: optimization of classifiers : 
% % % i_tp = 24;
% % % table_x = squeeze(mean(exc_stim_resps(train_inds, :, (i_tp-1)*skip_dt_units + (1:skip_dt_units)), 3));
% % % table_y = stim_y(train_inds);
% % %    
% % % % mdl_opt = fitcsvm(table_x, table_y, 'OptimizeHyperparameters', 'auto');
% % % %     % "regularize" table_x with a small bit of noise
% % % % %     table_x = table_x + 0.001*randn(size(table_x));
% % % % 
% % % % %     data_table = table(table_x, table_y);
% % % % %%
% % % c = cvpartition(size(table_x, 1), 'kfold', 10);
% % % box = optimizableVariable('box',[1e-5,1e5],'Transform','log');
% % % 
% % % minfn = @(z)kfoldLoss(fitcsvm(table_x,table_y,'CVPartition',c,...
% % %     'KernelFunction','linear','BoxConstraint',z.box));
% % % 
% % % results = bayesopt(minfn,[box],'IsObjectiveDeterministic',true,...
% % %     'AcquisitionFunctionName','expected-improvement-plus');
% % % 
% % % % opts2 = struct('Optimizer', 'bayesopt', 'ShowPlots', true, ... 
% % % %     'CVPartition', c);
% % % % svmmod = fitcsvm(table_x, table_y, 'kernelfunction', 'linear', ... 
% % % %     'optimizehyperparameters', 'auto', 'hyperparameteroptimizationoptions', opts2);
%% extra analysis: readout of excitatory or inhibitory neurons only

num_classifiers = floor(input_params.trial_length/c_timestep);
inhClassifier = cell(num_classifiers, 1);
inhValidationAccuracy = cell(num_classifiers, 1);
excClassifier = cell(num_classifiers, 1);

excValidationAccuracy = cell(num_classifiers, 1);
excISMClassifier = cell(num_classifiers, 1);
excISMValidationAccuracy = cell(num_classifiers, 1);

% classifier_options.classifier_type = 'linearSVM'; %'pLDA';
% classifier_options.boxconstraint = 100;     % not critically important from 1 to 10000 (higher is better (opt ~4000) but not by much and it takes longer
% classifier_options.applyPCA = false;
% classifier_options.full_save = true;

c_time_vec = (1:num_classifiers)*c_timestep;

table_y = stim_y(train_inds);
for i_tp = 1:num_classifiers
    table_x = squeeze(mean(inh_stim_resps(train_inds, :, (i_tp-1)*skip_dt_units + (1:skip_dt_units)), 3));
    
    % "regularize" table_x with a small bit of noise
%     table_x = table_x + 0.001*randn(size(table_x));
    

    data_table = table(table_x, table_y);
    
    [inhClassifier{i_tp}, inhValidationAccuracy{i_tp}] = ...
        trainGenericClassifier(data_table, classifier_options);
    
    table_x = squeeze(mean(exc_stim_resps(train_inds, :, (i_tp-1)*skip_dt_units + (1:skip_dt_units)), 3));
    
    % "regularize" table_x with a small bit of noise
%     table_x = table_x + 0.001*randn(size(table_x));

    data_table = table(table_x, table_y);
    
    [excClassifier{i_tp}, excValidationAccuracy{i_tp}] = ...
        trainGenericClassifier(data_table, classifier_options);
    
      
    table_x = squeeze(mean(exc_stim_resps_ISM(train_inds, :, (i_tp-1)*skip_dt_units + (1:skip_dt_units)), 3));
    
    % "regularize" table_x with a small bit of noise
%     table_x = table_x + 0.001*randn(size(table_x));

    data_table = table(table_x, table_y);
    
    [excISMClassifier{i_tp}, excISMValidationAccuracy{i_tp}] = ...
        trainGenericClassifier(data_table, classifier_options);
    
    disp(['i_tp = ' num2str(i_tp) ' done'])
end

%% compute the weights, accounting for the z-score. 
 
% for original classifiers
[classifier_weights, classifier_bs, classifier_betas, classifier_biases] = ...
    cellfun(@(x) linearSVMWeightsUnZscored(x.ClassificationEnsemble.BinaryLearners), ...
    allClassifier, 'uniformoutput', false);

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

%%
% inh_classifier_weights = cell2mat(cellfun(@(x) x.ClassificationEnsemble.BinaryLearners{1}.Beta', inhClassifier, 'UniformOutput', false));
% inh_classifier_biases = cellfun(@(x) x.ClassificationEnsemble.BinaryLearners{1}.Bias, inhClassifier);

% exc_classifier_weights = cell2mat(cellfun(@(x) x.ClassificationEnsemble.BinaryLearners{1}.Beta', excClassifier, 'UniformOutput', false));
% exc_classifier_biases = cellfun(@(x) x.ClassificationEnsemble.BinaryLearners{1}.Bias, excClassifier);


% excISM_classifier_weights = cell2mat(cellfun(@(x) x.ClassificationEnsemble.BinaryLearners{1}.Beta', excISMClassifier, 'UniformOutput', false));
% excISM_classifier_biases = cellfun(@(x) x.ClassificationEnsemble.BinaryLearners{1}.Bias, excISMClassifier);


%% 

inh_valAcc_vs_time = cellfun(@(x) x.validationAccuracy, inhValidationAccuracy);
exc_valAcc_vs_time = cellfun(@(x) x.validationAccuracy, excValidationAccuracy);
excISM_valAcc_vs_time = cellfun(@(x) x.validationAccuracy, excISMValidationAccuracy);

valAcc_vs_time = cellfun(@(x) x.validationAccuracy, allValidationAccuracy);

%% classifier prediction accuracy using mis-matched classifiers
exc_classIdataJ_valAcc = zeros(num_classifiers, num_classifiers);
inh_classIdataJ_valAcc = zeros(num_classifiers, num_classifiers);
excISM_classIdataJ_valAcc = zeros(num_classifiers, num_classifiers);
%%
table_y = stim_y(test_inds);
% for i_tc = 1:num_classifiers
    for i_td = 1:num_classifiers
        %%
        % first for excitation (all) 
        table_x = squeeze(nanmean(exc_stim_resps(test_inds, :, (i_td-1)*skip_dt_units + (1:skip_dt_units)), 3));
        data_table = table(table_x, table_y);
        pred_y = excClassifier{I_TC}.predictFcn(data_table);
        exc_classIdataJ_valAcc(I_TC, i_td) = mean(pred_y == stim_y(test_inds)); 

        % next for inhibition (all)
        table_x = squeeze(mean(inh_stim_resps(test_inds, :, (i_td-1)*skip_dt_units + (1:skip_dt_units)), 3));
        data_table = table(table_x, table_y);
        pred_y = inhClassifier{I_TC}.predictFcn(data_table);
        inh_classIdataJ_valAcc(I_TC, i_td) = mean(pred_y == stim_y(test_inds)); 
    
        % finally for the subset of excitation 
        table_x = squeeze(mean(exc_stim_resps_ISM(test_inds, :, (i_td-1)*skip_dt_units + (1:skip_dt_units)), 3));
        data_table = table(table_x, table_y);
        pred_y = excISMClassifier{I_TC}.predictFcn(data_table);
        excISM_classIdataJ_valAcc(I_TC, i_td) = mean(pred_y == stim_y(test_inds)); 
        
        disp(['i_td = ' num2str(i_td) ' done'])
        
    end
% end


%% readout accuracy and weight distributions

w_example = classifier_weights(I_TC, :);
% w_example = normClassifier{I_TC}.ClassificationEnsemble.BinaryLearners{1}.Beta;
% b_ex = normClassifier{I_TC}.ClassificationEnsemble.BinaryLearners{1}.Bias;
% compute activity projections into readout space
use_inds = test_inds;
table_y = stim_y(use_inds);
table_x = squeeze(nanmean(stim_resps(use_inds, :, (I_TC-1)*skip_dt_units + (1:skip_dt_units)), 3));
data_table = table(table_x, 0*table_y);
pred_y = allClassifier{I_TC}.predictFcn(data_table);
        
% wx_resps = imfilter(squeeze(sum(bsxfun(@times, stim_resps(use_inds, :, :), w_example), 2)), ones(1, skip_dt_units)/skip_dt_units);
low_pass_resps = imfilter(stim_resps, ones(1, 1, skip_dt_units)/skip_dt_units);
wx_resps = squeeze(sum(bsxfun(@times, low_pass_resps(use_inds, :, :), w_example), 2));

%%
% 
% mu_vals = normClassifier{I_TC}.ClassificationEnsemble.BinaryLearners{1}.Mu;
% sigma_vals = normClassifier{I_TC}.ClassificationEnsemble.BinaryLearners{1}.Sigma;
% lp_resp_z = bsxfun(@times, bsxfun(@minus, low_pass_resps, mu_vals), 1./sigma_vals);
% 
% wx_resps = squeeze(sum(bsxfun(@times, lp_resp_z(use_inds, :, :), w_example), 2));

correct_trials = {pred_y.*table_y == 1, pred_y.*table_y == 4};
incorrect_trials = {pred_y == 2 & table_y == 1, pred_y == 1 & table_y == 2};

%% activity in pca-space (collapsing cell dimension)
cell_resp_ave = squeeze(mean(stim_resps(:, :, 100:500), 3));
[pc_cell, cf_cell, ~, ~, vE_cell] = pca(cell_resp_ave);
use_resps = stim_resps(use_inds, :, :);

parsedTR_ave_resp = cellfun(@(x) squeeze(mean(use_resps(x, :, :), 1)), ... 
    [correct_trials incorrect_trials], 'UniformOutput', false);
%%
parsedTR_pca_x = cellfun(@(x) pc_cell(:, 1)'*x, parsedTR_ave_resp, ... 
    'UniformOutput', false);
parsedTR_pca_y = cellfun(@(x) pc_cell(:, 2)'*x, parsedTR_ave_resp, ... 
    'UniformOutput', false);
% stim1_ave_resp = squeeze(mean(stim_resps(stim_y == 1, :, :), 1));
% stim2_ave_resp = squeeze(mean(stim_resps(stim_y == 2, :, :), 1));


%%
stim1_pc_traj = pc_cell(:, 1:2)'*stim1_ave_resp;
stim2_pc_traj = pc_cell(:, 1:2)'*stim2_ave_resp;
%%


[w_sorted, wt_ord] = sort(w_example);
l2_ord = l2_inds(wt_ord);
%%
% % this uses the cross-validation set for accuracy 
% all_valAccs = [exc_valAcc_vs_time inh_valAcc_vs_time excISM_valAcc_vs_time valAcc_vs_time];
% % this uses the reserved test set for accuracy (cross val worked)
all_valAccs = [exc_classIdataJ_valAcc inh_classIdataJ_valAcc excISM_classIdataJ_valAcc classIdataJ_valAcc];
all_valAccs = reshape(all_valAccs(I_TC, :), length(c_time_vec), 4);

all_weights_raw = {exc_classifier_weights, inh_classifier_weights, excISM_classifier_weights, classifier_weights};
        
all_weights = cellfun(@(x) diag(1./max(abs(x), [], 2))*x, all_weights_raw, 'UniformOutput', false);

num_units = cellfun(@(x) size(x, 2), all_weights, 'uniformoutput', false);
all_labels = {'exc units', 'inh units', 'exc units', 'exc+inh units'};
legent = cellfun(@(x,y) [num2str(x) ' ' y], num_units, all_labels, 'UniformOutput', false);

cmap = lines(size(all_valAccs, 2));

makeMyFigure(24, 16);
sh_valAcc = subplot(5, 3, [2 5]);
sh_wts = cell(4, 1);
for i_sw = 1:length(sh_wts)
    sh_wts{i_sw} = subplot(4, 3, 3*i_sw); 
end
sh_wtDist = subplot(5, 3, [11 14]);
sh_netAct = subplot(5, 3, [10 13]);
sh_pcaAct = subplot(5, 3, [1 4]);

% plot activity in pca space 
line_styles = {'-', '-.', '-', '-.'};
line_colors = [0 0 0; 0 0 0; 1 1 1; 1 1 1]*.5;
line_widths = [1.5 1.5 .5 .5 ];

[~, t_vec_ind] = min(abs(t_vec - c_time_vec(I_TC)));

set(gcf, 'currentaxes', sh_pcaAct);
hold on
for i_ptr = 1:2 %length(parsedTR_pca_x)
    plot(parsedTR_pca_x{i_ptr}, parsedTR_pca_y{i_ptr}, ...
        'linestyle', line_styles{i_ptr}, ... 
        'color', line_colors(i_ptr, :), ... 
        'LineWidth', line_widths(i_ptr));
    
    plot(parsedTR_pca_x{i_ptr}(t_vec_ind), parsedTR_pca_y{i_ptr}(t_vec_ind), ... 
        'o', 'color', cmap(4, :), 'linewidth', 3)
    
end
set(gca, 'color', 'none');
xlabel(['neuron pca dim 1 (' num2str(vE_cell(1), '%1.0f') '% of var)'])
ylabel(['neuron pca dim 2 (' num2str(vE_cell(2), '%1.0f') '% of var)'])

% plot the weighted/summed activity for correct trials
set(gcf, 'CurrentAxes', sh_netAct);
wx_mean = cellfun(@(x) mean(wx_resps(x, :), 1)', [correct_trials incorrect_trials], ...
    'uniformoutput', false);
hold on
for i_ptr = 1:length(parsedTR_pca_x)
    
    plot(t_vec, wx_mean{i_ptr}, ...
        'linestyle', line_styles{i_ptr}, ... 
        'color', line_colors(i_ptr, :), ... 
        'LineWidth', line_widths(i_ptr));
end
% mean_wx_correct = cell2mat(cellfun(@(x) mean(wx_resps(x, :), 1)', correct_trials, ...
%     'uniformoutput', false));
% sd_wx_correct = cell2mat(cellfun(@(x) mean(wx_resps(x, :), 1)', correct_trials, ...
%     'uniformoutput', false));
% mean_wx_incorrect = cell2mat(cellfun(@(x) mean(wx_resps(x, :), 1)', incorrect_trials, ...
%     'uniformoutput', false));
% sd_wx_incorrect = cell2mat(cellfun(@(x) mean(wx_resps(x, :), 1)', incorrect_trials, ...
%     'uniformoutput', false));
% hold on
% ph = plot(t_vec, mean_wx_correct, 'k', 'linewidth', 1);
% ph(1).LineStyle = '-';
% ph(2).LineStyle = '-.';
% ph2 = plot(t_vec, mean_wx_incorrect, 'color', 0.5*[1 1 1], 'linewidth', 1);
% ph2(1).LineStyle = '-';
% ph2(2).LineStyle = '-.';

xlim([0 50])
set(gca, 'color', 'none')

v = axis();
plot(c_time_vec(I_TC)*[1 1], v(3:4), 'color', cmap(4, :), 'linewidth', 1.5)
lh = legend({'low-freq (correct)', 'high-freq (correct)', 'low-freq (classifer wrong)', ... 
    'high-freq (classifier wrong)', 'classifier fit time'}, 'position', [0.15 0.42 0.1878 0.1100], ... 
    'box', 'off'); 
ylabel('weighted sum of firing rates')
xlabel('time')


set(gcf, 'CurrentAxes', sh_valAcc);
hold on
plot(c_time_vec, all_valAccs, 'LineWidth', 1.5)
hold off
set(gca, 'color', 'none')
axis square
% legend(legent, 'location', 'best')

ylabel('fraction correct (test set)')
xlabel('time relative to stim onset')

for i_wt = 1:4
    [~, ord] = sort(abs(all_weights{i_wt}(I_TC, :)), 'descend');
    these_weights = all_weights{i_wt}(I_TC, ord);
    set(gcf, 'CurrentAxes', sh_wts{i_wt});
    bar(these_weights, 'FaceColor', cmap(i_wt, :))
    title(['weights for ' num2str(num_units{i_wt}) ' ' all_labels{i_wt}])
    axis tight
    ylim([-1 1])
    ylabel('readout weight') 
    set(gca, 'color', 'none')

end

% subplot(3, 2, 2)
% bar([zeros(1, length(l1_inds)) w_sorted], 'FaceColor', cmap(i_wt, :))
% xlim([0.5 num_neur+.5])
% ylabel('readout weight')
% xlabel('neuron (sorted by weight; directly stimulated cells not included)')
% axis square
% ch = colorbar;
% ch.Visible = 'off';
% set(gca, 'color', 'none')
% title(['weights for ' num2str(num_units{i_wt}) ' ' all_labels{i_wt}])
% 

set(gcf, 'currentaxes', sh_wtDist);
hold on
br_val = 3*std(all_weights{4}(I_TC, :));
wt_bins = br_val*linspace(-1, 1, 61);
p_norm = normcdf(wt_bins', mean(all_weights{4}(I_TC, :)), std(all_weights{4}(I_TC, :)));

for i_wt = 1:4
    histogram(all_weights{i_wt}(I_TC, :), wt_bins,'Normalization', 'cdf', ...
        'DisplayStyle', 'stairs', 'LineWidth', 1.5)
end
xlabel('readout weights')
ylabel('cdf')
plot(wt_bins, p_norm, 'k', 'LineWidth', 1.5)
ylim([0 1])
xlim(wt_bins([1 end]))
set(gca, 'color', 'none')

axis square
wts_legent = cellfun(@(x,y) [num2str(x) ' ' y], num_units, all_labels, 'UniformOutput', false);
lh = legend([wts_legent, 'normal cdf'], 'Position', [0.4194 0.4146 0.1877 0.11], 'box', 'off');

print(gcf, '-dpdf', [network_param.plotdir 'weightDists' classifier_options.classifier_type '_' network_param.network_tag 'pw' num2str(100*pulse_width)]);




%% connectivity analysis and plots

% simple matrix analysis
full_weights = zeros(num_neur, 1);
full_weights(postL1_cells) = w_example;


num_l1_exc_inputs = sum(network_param.J(:, ~postL1_cells) > 0, 2);
strength_direct_inputs = sum(network_param.J(:, ~postL1_cells), 2);

num_l2_exc_inputs = sum(network_param.J(:, postL1_cells) > 0, 2);



connections_from_inh = sum(network_param.J(:, :) < 0, 2);
connections_to_inh = sum(sign(network_param.J(inh_cells, :)), 1)';

factor_mat = [full_weights num_l1_exc_inputs strength_direct_inputs num_l2_exc_inputs connections_from_inh connections_to_inh];

[xc_mat, p_xc_mat] = corrcoef(factor_mat(postL1_cells, :));
%% 
figure()
if strcmp(network_param.network_tag(1:8), 'spnorm2P')

dr_out = proc_simdata.dr_out_peristim;
e1_resp = squeeze(mean(dr_out(:, 1:nE/2, :), 2));
e2_resp = squeeze(mean(dr_out(:, nE/2 + (1:nE/2), :), 2));
i1_resp = squeeze(mean(dr_out(:, nE + (1:nI/2), :), 2));
i2_resp = squeeze(mean(dr_out(:, nE + nI/2 + (1:nI/2), :), 2));

pop_resps = {e1_resp, e2_resp, i1_resp, i2_resp};

pop_stim1_ave_resps = cellfun(@(x) squeeze(mean(x(stim_IDs, :), 1)), pop_resps, ...
    'UniformOutput', false);
pop_stim2_ave_resps = cellfun(@(x) squeeze(mean(x(~stim_IDs, :), 1)), pop_resps, ...
    'UniformOutput', false);

pop_line_types = {'-', '--', ':', '-.'};
legents = cell(8, 1);
pop_name = {'E1', 'E2', 'I1', 'I2'};
subplot(1, 2, 1)
hold on
for i_pop = 1:4
    plot(t_vec, pop_stim1_ave_resps{i_pop}, 'LineStyle', pop_line_types{i_pop}, ...
        'color', stim_colors(1, :), 'linewidth', 1.5)
    plot(t_vec, pop_stim2_ave_resps{i_pop}, 'LineStyle', pop_line_types{i_pop}, ...
        'color', stim_colors(2, :), 'linewidth', 1.5)
   
    legents{i_pop*2 - 1} = ['stim 1, ' pop_name{i_pop}];
    legents{i_pop*2} = ['stim 2, ' pop_name{i_pop}];
end

axis square
set(gca, 'color', 'none')
legend(legents, 'location', 'best')
 
end
subplot(3, 2, 2)
bar([zeros(1, length(l1_inds)) w_sorted], 'FaceColor', cmap(4, :))
xlim([0.5 num_neur+.5])
ylabel('readout weight')
xlabel('neuron (sorted by weight; directly stimulated cells not included)')
axis square
ch = colorbar;
ch.Visible = 'off';
set(gca, 'color', 'none')
title({'correlation between weight and # of direct inputs:'; ...
    [num2str(xc_mat(1, 2), '%1.2f') num2str(p_xc_mat(1, 2), '(p %1.2f)') ] } )


subplot(3, 2, 4)
neuron_order = [l1_inds; l2_ord];
imagesc(1:num_neur, 1:num_neur, network_param.J(neuron_order, neuron_order));
axis square
hold on
plot(length(l1_inds)*[0 1 1 0 0]+.5, length(l1_inds)*[1 1 0 0 1]+.5, 'w', 'linewidth', 2)
% plot([1 1], [1 length(l1_inds)], 'w', 'linewidth', 2)
xlabel('neuron j (sorted by classifier weight)')
ylabel('neuron i (sorted by classifier weight)')
colorbar
title({'connection strength from j to i'; 'no obvious structure related to weights'})



% weights and activity analysis. 
% downsample/average stim resps

% 
% [ntr, nnr, ~] = size(stim_resps);
% ds_rs_stim_resps = reshape(stim_resps(:, :, (1:skip_dt_units*num_classifiers)), [ntr nnr skip_dt_units num_classifiers]);
% ds_stim_resps = squeeze(nanmean(ds_rs_stim_resps, 3));
% 
% use_inds = test_inds;
% resp_arr = ds_stim_resps(use_inds, :, :);
% stim_labels = stim_y(use_inds);
% w_example = classifier_weights(i_tc, :);

% table_x = squeeze(resp_arr(:, :, i_tc));
% resp_tc = table(table_x);
% classifier_prediction = normClassifier{i_tc}.predictFcn(resp_tc);
% classifier_decision = normValidationAccuracy{i_tc}.validationPredictions;
% weighted_out = squeeze(sum(bsxfun(@times, resp_arr, w_example), 2));
% validation_accuracy_tc = classIdataJ_valAcc(i_tc, :);
% 
% % figure()
% % subplot(2, 1, 1)
% set(gcf, 'CurrentAxes', sh_wtrd);
% hc_wx_s1 = [];
% hold on
% ph1 = plot(c_time_vec, mean(weighted_out(stim_labels == 1 & classifier_prediction == 1, :)), 'k', 'markersize', 4);
% ph2 = plot(c_time_vec+.05, mean(weighted_out(stim_labels == 2 & classifier_prediction == 1, :)), 'k-.', 'markersize', 4);
% ph3 = plot(c_time_vec+.1, mean(weighted_out(stim_labels == 1 & classifier_prediction == 2, :)), 'r', 'markersize', 4);
% ph4 = plot(c_time_vec+.15, mean(weighted_out(stim_labels == 2 & classifier_prediction == 2, :)), 'r-.', 'markersize', 4);
% % % xlim(i_tc + [-1.5 1.5])
% axis tight
% v = axis();
% plot(c_time_vec(i_tc)*[1 1], v(3:4))
% ylabel('weighted readout')
% set(gca, 'color', 'none')
% 
% legend([ph1(1);ph2(1);ph3(1);ph4(1)], ...
%     {'stim 1 (pred: 1)', 'stim 2 (pred: 1)', 'stim 1 (pred: 2)', 'stim 2 (pred: 2)'}, ...
%     'location', 'south')



