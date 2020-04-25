% function fitClassifierModelData(network_code, pulse_width, rng_index, f2f_type)
% Decoder fitting script - freq vs 2freq experiments (only two stims)


% network_param = generateNetworkParam(network_code, rng_index);

results_dir = '/Volumes/home/asederberg6/projects/random_network_readouts/results/';

% num_neur = network_param.numNeur;
% file_tag = network_param.network_tag; %

%%%% or set values individually, if running as a script
num_neur = 100;
gain_neur = 400;
rng_val = 2;
ei_ratio = 100; % 92 or 108;  80 or 125
pulse_width = .50; % 50 or 150 or 500
connect_tag = 'spnorm';
f2f_type = 'low';

file_tag = [connect_tag '_g' num2str(gain_neur) '_n' num2str(num_neur) '_ei' num2str(ei_ratio) '_rng' num2str(rng_val)];
%%%%

% load('/Users/audreysederberg/Dropbox/postdoctoral projects/stanley lab work/code/theory code/results/norm_g100_n100_ei80_rng1/matfiles/sim_data_norm_g100_n100_ei80_rng1_f2fdiscExp_pw150.mat')
x_res = load([results_dir file_tag '/matfiles/sim_data_' file_tag '_f2f' f2f_type 'discExp_pw' num2str(100*pulse_width) '.mat']);

proc_simdata = x_res.proc_simdata;
input_params = x_res.input_params;

network_param = x_res.network_param;
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

stim_resps = cell2mat(cellfun(@(x) proc_simdata.dr_out_peristim(x, :, :), trial_inds, ...
    'UniformOutput', false));

% option to only use responses from cells that are not directly stimulated
postL1_cells = network_param.input_pattern == 0;
stim_resps = stim_resps(:, postL1_cells, :);

train_inds = (1:2:size(stim_y, 1))';
test_inds = setdiff((1:size(stim_y, 1))', train_inds);

c_timestep = 1; % in units of tau, how often to fit a new classifier
skip_dt_units = round(c_timestep/diff(proc_simdata.peri_stim_time(1:2)));
%
num_classifiers = floor(input_params.trial_length/c_timestep);
normClassifier = cell(num_classifiers, 1);
normValidationAccuracy = cell(num_classifiers, 1);

classifier_options.classifier_type = 'linearSVM'; %'pLDA';
classifier_options.applyPCA = false;
classifier_options.full_save = true;

c_time_vec = (1:num_classifiers)*c_timestep;
for i_tp = 1:num_classifiers
    table_x = squeeze(mean(stim_resps(train_inds, :, (i_tp-1)*skip_dt_units + (1:skip_dt_units)), 3));
    
    % "regularize" table_x with a small bit of noise
%     table_x = table_x + 0.01*randn(size(table_x));
    
    table_y = stim_y(train_inds);
    data_table = table(table_x, table_y);
    
    [normClassifier{i_tp}, normValidationAccuracy{i_tp}] = ...
        trainGenericClassifier(data_table, classifier_options);
    
   
end

% can use bayesopt function to set box param 
% % % getting the parameters for the separation hyperplane from the linear
% svm struct. 
% % % parameters from svmStruct
% % w1 = dot(svmStruct.Alpha, svmStruct.SupportVectors(:,1));
% % w2 = dot(svmStruct.Alpha, svmStruct.SupportVectors(:,2));
% % bias = svmStruct.Bias;
% % 
% % % y = a*x + b
% % a = -w1/w2;
% % b = -svmStruct.Bias/w2;
%%
if strcmp(classifier_options.classifier_type, 'pLDA')
    pinv_eps = 0;
    classifier_weights = cell2mat(cellfun(@(x) (pinv(x.ClassificationEnsemble.Sigma, pinv_eps)*diff(x.ClassificationEnsemble.Mu, 1, 1)')', normClassifier, 'UniformOutput', false));
else
    classifier_weights = cell2mat(cellfun(@(x) x.ClassificationEnsemble.BinaryLearners{1}.Beta', normClassifier, 'UniformOutput', false));
    classifier_biases = cellfun(@(x) x.ClassificationEnsemble.BinaryLearners{1}.Bias, normClassifier);
end
%% classifier prediction accuracy using mis-matched classifiers
classIdataJ_valAcc = zeros(num_classifiers, num_classifiers);
for i_tc = 1:num_classifiers
    for i_td = 1:num_classifiers
        %%
        table_x = squeeze(nanmean(stim_resps(test_inds, :, (i_td-1)*skip_dt_units + (1:skip_dt_units)), 3));
        
%         table_x = table_x + randn(size(table_x));
        
        table_y = stim_y(test_inds);
        data_table = table(table_x, table_y);
        pred_y = normClassifier{i_tc}.predictFcn(data_table);
        classIdataJ_valAcc(i_tc, i_td) = mean(pred_y == stim_y(test_inds)); 
    end
end

%%
valAcc_vs_time = cellfun(@(x) x.validationAccuracy, normValidationAccuracy);

%% 
t_vec = proc_simdata.peri_stim_time;
t_order = 150;
[~, full_ord] = sort( squeeze(mean(stim_resps(stim_y == 1, :, t_order), 1)));
ord = full_ord();
is_inh = all(network_param.J <= 0, 1);
is_inh_ordered = is_inh(ord);
stim1_ave_resp = bsxfun(@plus, proc_simdata.baseline_activity(ord)', squeeze(mean(stim_resps(stim_y == 1, ord, :), 1)));
stim2_ave_resp = bsxfun(@plus, proc_simdata.baseline_activity(ord)', squeeze(mean(stim_resps(stim_y == 2, ord, :), 1)));
stim1_std_resp = squeeze(std(stim_resps(stim_y == 1, ord, :),[], 1));
stim2_std_resp = squeeze(std(stim_resps(stim_y == 2, ord, :),[], 1));

makeMyFigure(30, 40);
subplot(531)
imagesc(t_vec, [], stim1_ave_resp, [-.3 1])
hold on
plot(t_vec(t_order)*[1 1], [1 num_neur], 'r')
plot(t_vec(t_order), find(is_inh_ordered), 'wo')
title({'stim 1 (white o -> inh)';'ordered by stim1 activation order'})

subplot(532)
imagesc(t_vec, [], stim2_ave_resp, [-.3 1])
title({'stim 2';'ordered by stim1 activation order'})

subplot(533)

c_lim = max(max(abs(stim1_ave_resp - stim2_ave_resp)))*[-1 1]*.8;
imagesc(t_vec, [1 size(stim1_ave_resp, 1)], stim1_ave_resp - stim2_ave_resp, c_lim)
hold on
plot(t_vec(t_order), find(is_inh_ordered), 'w.')
% axis image
title({'stim 1 minus stim 2';'ordered by stim1 activation order'})

try
% calculate average peri-event triggered response
peri_event_resps = cellfun(@(x) shiftdim(squeeze(mean(cell2mat(x), 1)), -1), ...
    proc_simdata.peri_event_resps_cell, 'uniformoutput', false);
stim_aves = cellfun(@(x) (peri_event_resps(x)), trial_inds, 'UniformOutput', false);
stim_aves = cellfun(@(x) squeeze(mean(cell2mat(x(cellfun(@(y) size(y, 3) > 0, x))), 1)),...
    stim_aves, 'UniformOutput', false);

subplot(5, 3, 4)
ph1 = plot(stim_aves{1}(ord, :)', 'linewidth', 0.25);
assignColorsToLines(ph1, parula(length(ph1)));
title('input-triggered average: stim 1 (low-freq)')
axis tight



subplot(5, 3, 5)
ph1 = plot(stim_aves{2}(ord, :)', 'linewidth', 0.25);
assignColorsToLines(ph1, parula(length(ph1)));
title('input-triggered average: stim 2 (high-freq)')
axis tight


stimave1_z = zscore(stim_aves{1}, [], 2);
subplot(5, 3, 7)
ph1 = plot(stimave1_z(ord, :)', 'linewidth', 0.25);
assignColorsToLines(ph1, parula(length(ph1)));
title('z-scored input-triggered average: stim 1 (low-freq)')
axis tight

stimave2_z = zscore(stim_aves{2}, [], 2);
subplot(5, 3, 8)
ph1 = plot(stimave2_z(ord, :)', 'linewidth', 0.25);
assignColorsToLines(ph1, parula(length(ph1)));
title('z-scored input-triggered average: stim 2 (high-freq)')
axis tight
end

subplot(5, 3, 6)
hold on
xx = stim1_ave_resp(:, 1:10:300);
yy = stim2_ave_resp(:, 1:10:300);
dxx = stim1_std_resp(:, 1:10:300);
dyy = stim2_std_resp(:, 1:10:300);
diff_s1minus2 = xx(~is_inh_ordered, :) - yy(~is_inh_ordered, :);
d_bins = 0.05*linspace(-1, 1, 21) + .025;
cmap = jet(size(diff_s1minus2, 2));
hold on
for i_t = 1:size(diff_s1minus2, 2)
    histogram(diff_s1minus2(:, i_t), d_bins, 'FaceColor', cmap(i_t, :));
end
% ph2 = plot(xx(is_inh_ordered, :), yy(is_inh_ordered, :), '.');
% assignColorsToLines(ph2, jet(length(ph2)));
% eqline
xlabel('activity (stim1 minus stim2) : excitatory cells')
% ylabel('stim 2 activity : excitatory cells')
axis square
title('ave activity excitatory cells')


subplot(5, 3, 9)
diff_s1minus2 = xx(is_inh_ordered, :) - yy(is_inh_ordered, :);
d_bins = 0.025*linspace(-1, 1, 21) + .0125;
hold on
for i_t = 1:size(diff_s1minus2, 2)
    histogram(diff_s1minus2(:, i_t), d_bins, 'FaceColor', cmap(i_t, :));
end
% ph2 = plot(xx(is_inh_ordered, :), yy(is_inh_ordered, :), '.');
% assignColorsToLines(ph2, jet(length(ph2)));
% eqline
xlabel('activity (stim1 minus stim2) : inhibitory cells')
% ylabel('stim 2 activity : inhibitory cells')
axis square
title('ave diff activity inhibitory cells')



% plot classifier weights and classifer accuracy
[~, best_vA_ind] = max(valAcc_vs_time);
[~, last_vA_ind] = max(valAcc_vs_time(round(length(valAcc_vs_time)/2) : end));
last_vA_ind = last_vA_ind + round(length(valAcc_vs_time)/2) - 1;
subplot(5, 3, 10)
plot(c_time_vec, valAcc_vs_time)
hold on
plot(c_time_vec([best_vA_ind last_vA_ind]), valAcc_vs_time([best_vA_ind last_vA_ind]), 'o')
xlabel('time (tau units)')
ylabel('accuracy of optimal classifier at each time')
title(network_param.network_tag, 'Interpreter', 'none');

c_wt_lims = max(max(abs(classifier_weights(10:end, :)), [], 2));
subplot(5, 3, 13)
plot(c_time_vec, classifier_weights)
xlabel('time (tau units)')
ylabel('classifier weights')
ylim(c_wt_lims*[-1 1])
title('moment-by-moment classifier weights')

subplot(5, 3, 11)
plot(c_time_vec, classIdataJ_valAcc(best_vA_ind, :))
xlabel('time (tau units)')
ylabel('accuracy of best late classifier')
title(['classifier fit at time ' num2str(c_time_vec(best_vA_ind))])

subplot(5, 3, 14)
[~, c_weight_inds] = sort(abs(classifier_weights(best_vA_ind, :)), 'descend');
bar(classifier_weights(best_vA_ind, c_weight_inds))
title(['weights fit at time ' num2str(c_time_vec(best_vA_ind))])
ylim(c_wt_lims*[-1 1])


subplot(5, 3, 12)
plot(c_time_vec, classIdataJ_valAcc(last_vA_ind, :))
xlabel('time (tau units)')
ylabel('accuracy of best classifier')
title(['classifier fit at time ' num2str(c_time_vec(last_vA_ind))])

subplot(5, 3, 15)
[~, c_weight_inds] = sort(abs(classifier_weights(last_vA_ind, :)), 'descend');
bar(classifier_weights(last_vA_ind, c_weight_inds))
title(['weights fit at time ' num2str(c_time_vec(last_vA_ind))])
ylim(c_wt_lims*[-1 1])

print(gcf, '-dpdf', [network_param.plotdir 'classifier' classifier_options.classifier_type 'Overview_' network_param.network_tag  f2f_type 'pw' num2str(100*pulse_width)]);

% %% stim input variability - a check
% 
% makeMyFigure(20, 10);
% 
% time_bins = linspace(0, input_params.trial_length, 17);
% 
% for i_sp = 1:2
%     input_isp_times = cellfun(@(x,y) x - y, proc_simdata.precise_impulse_times(trial_inds{i_sp}), ...
%         num2cell(input_params.stim_start_times{i_sp}'), 'uniformoutput', false);
%     
%     subplot(3, 2, i_sp)
%     histogram(cell2mat(input_isp_times'), time_bins)
%     
%     subplot(3, 2, i_sp + 2)
%     histogram(cellfun(@length, input_isp_times))
%     
%     iii_dist = cellfun(@(x) diff(x), input_isp_times, 'uniformoutput', false);
%     subplot(3, 2, i_sp + 4)
%     histogram(cell2mat(iii_dist'))
%     
% end
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


stim_IDs = input_params.input_ampfreq(:, 2) == input_params.input_frequency{1};
in1_freq = input_params.input_frequency{1};
in2_freq = input_params.input_frequency{2};
in1_mat = in_mat(:,  stim_IDs);
in2_mat = in_mat(:, ~stim_IDs);

% average network-wide firing rates
aveS1_netfr = mean(mean(stim_resps(stim_IDs, :, :), 3), 2);
aveS2_netfr = mean(mean(stim_resps(~stim_IDs, :, :), 3), 2);


% fourier transform analysis of input current 
pulse_fun = @(t, a) (t.^2/a^2).*exp((t>0).*(1-t/a)).*(t > 0);

t_trial = 1/sample_freq:1/sample_freq:input_params.isi;
stim_on_time = t_trial > -offset_time & t_trial <= input_params.trial_length-offset_time;
pulse_fft = fft(pulse_fun(t_trial(stim_on_time), pulse_width));

in1_fft = fft(in1_mat(stim_on_time, :));
in1_asd = abs(bsxfun(@(x,y) x./y, in1_fft(1:(end/2+1), :), pulse_fft(1:(end/2+1))'));
in1_asd(2:end-1, :) = in1_asd(2:end-1, :)*2;

in2_fft = fft(in2_mat(stim_on_time, 2:end));
in2_asd = abs(bsxfun(@(x,y) x./y, in2_fft(1:(end/2+1), :), pulse_fft(1:(end/2+1))'));
in2_asd(2:end-1, :) = in2_asd(2:end-1, :)*2;
nu_vals = (0:size(in1_asd, 1)-1)*sample_freq*netwk_hz/2/size(in1_asd, 1);

% binary input mat: +1 is an impulse time on a low-freq trial, -1 is an
% impulse time on a high-freq trial
round_fac = 2;
in1_times = round(round_fac*(input_params.input_times{1} - input_params.burn_in_time - input_params.isi + 10));
in2_times = round(round_fac*(input_params.input_times{2} - input_params.burn_in_time - input_params.isi + 10));

binary_inputs = sparse(1, in1_times, 1, 1, round_fac*input_params.isi*input_params.num_reps*2) - ...
    sparse(1, in2_times, 1, 1, round_fac*input_params.isi*input_params.num_reps*2);
binary_mat = reshape(binary_inputs, input_params.isi*round_fac, 2*input_params.num_reps);

bin1_mat = binary_mat(:, stim_IDs);
bin2_mat = binary_mat(:, ~stim_IDs);

makeMyFigure(30, 20);

sh_in1 = subplot(4, 2, 1);
sh_in2 = subplot(4, 2, 3);
fft_in12 = subplot(2, 4, 5);
netfr_in12 = subplot(2, 4, 6);
sh_wtrd = subplot(2, 2, 2);
sh_valAc = subplot(2, 2, 4);


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
mPMsd = mean(in1_mat, 2) + std(in1_mat, [], 2)*[-1 0 1];
plot(t_trial, mPMsd, '-', 'linewidth', 2)
ylim([0 4])

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
mPMsd = mean(in2_mat, 2) + std(in2_mat, [], 2)*[-1 0 1];

plot(t_trial, mPMsd, '-', 'linewidth', 2)
ylim([0 4])
ylabel('ave input current +/- sd across trials')
set(get(gca, 'YLabel'), 'Position', [76 5.7460 0], 'Rotation', -90);

set(gcf, 'CurrentAxes', fft_in12);
hold on
plot(nu_vals, log10(mean(in1_asd.^2, 2)/sum(mean(in1_asd.^2, 2))))
plot(nu_vals, log10(mean(in2_asd.^2, 2)/sum(mean(in2_asd.^2, 2))))
xlabel('frequency (Hz)')
ylabel('ave pow (normalized; logscale)')
axis tight
xlim([0 50])
axis square
set(gca, 'color', 'none')
legend({num2str(netwk_hz*in1_freq, 'freq : %1.0f'), num2str(netwk_hz*in2_freq, 'freq : %1.0f')})

set(gcf, 'CurrentAxes', netfr_in12);
hold on
histogram(aveS1_netfr)
histogram(aveS2_netfr)
xlabel('network average firing rate (au)')
ylabel('trial counts')
axis square
set(gca, 'color', 'none')
legend({num2str(netwk_hz*in1_freq, 'freq : %1.0f'), num2str(netwk_hz*in2_freq, 'freq : %1.0f')})

% weights and activity analysis. 
% downsample/average stim resps


[ntr, nnr, ~] = size(stim_resps);
ds_rs_stim_resps = reshape(stim_resps(:, :, 1:skip_dt_units*num_classifiers), [ntr nnr skip_dt_units num_classifiers]);
ds_stim_resps = squeeze(nanmean(ds_rs_stim_resps, 3));

resp_arr = ds_stim_resps(test_inds, :, :);
stim_labels = stim_y(test_inds);
i_tc = 50;
w_example = classifier_weights(i_tc, :);

table_x = squeeze(resp_arr(:, :, i_tc));
resp_tc = table(table_x);
classifier_prediction = normClassifier{i_tc}.predictFcn(resp_tc);
classifier_decision = normValidationAccuracy{i_tc}.validationPredictions;
weighted_out = squeeze(sum(bsxfun(@times, resp_arr, w_example), 2));
validation_accuracy_tc = classIdataJ_valAcc(i_tc, :);

% figure()
% subplot(2, 1, 1)
set(gcf, 'CurrentAxes', sh_wtrd);
hc_wx_s1 = [];
hold on
ph1 = plot(c_time_vec, mean(weighted_out(stim_labels == 1 & classifier_prediction == 1, :)), 'k', 'markersize', 4);
ph2 = plot(c_time_vec+.05, mean(weighted_out(stim_labels == 2 & classifier_prediction == 1, :)), 'k-.', 'markersize', 4);
ph3 = plot(c_time_vec+.1, mean(weighted_out(stim_labels == 1 & classifier_prediction == 2, :)), 'r', 'markersize', 4);
ph4 = plot(c_time_vec+.15, mean(weighted_out(stim_labels == 2 & classifier_prediction == 2, :)), 'r-.', 'markersize', 4);
% % xlim(i_tc + [-1.5 1.5])
axis tight
v = axis();
plot(c_time_vec(i_tc)*[1 1], v(3:4))
ylabel('weighted readout')
set(gca, 'color', 'none')

legend([ph1(1);ph2(1);ph3(1);ph4(1)], ...
    {'stim 1 (pred: 1)', 'stim 2 (pred: 1)', 'stim 1 (pred: 2)', 'stim 2 (pred: 2)'}, ...
    'location', 'south')

% subplot(2, 1, 2)
set(gcf, 'CurrentAxes', sh_valAc);
hold on
plot(c_time_vec, validation_accuracy_tc, 'k', 'linewidth', 1)
plot(c_time_vec(i_tc)*[1 1], [.45 1])
axis tight
ylabel('classifier accuracy')
xlabel('time (tau units)')
set(gca, 'color', 'none')


print(gcf, '-dpdf', [network_param.plotdir 'classifier' classifier_options.classifier_type 'MethodsOverview_' network_param.network_tag 'pw' num2str(100*pulse_width)]);


%% connectivity, sorted by weights
l1_inds = find(~postL1_cells);

l2_inds = find(postL1_cells);
[w_sorted, wt_ord] = sort(w_example);
l2_ord = l2_inds(wt_ord);

makeMyFigure(20, 30);

subplot(3, 1, 1)
bar([zeros(1, length(l1_inds)) w_sorted])
xlim([0.5 num_neur+.5])
ylabel('readout weight')
xlabel('neuron (sorted by weight; directly stimulated cells not included)')
ch = colorbar;
ch.Visible = 'off';
subplot(3, 1, [2 3])
neuron_order = [l1_inds; l2_ord];
imagesc(1:num_neur, 1:num_neur, network_param.J(neuron_order, neuron_order));
hold on
plot(length(l1_inds)*[0 1 1 0 0]+.5, length(l1_inds)*[1 1 0 0 1]+.5, 'w', 'linewidth', 2)
% plot([1 1], [1 length(l1_inds)], 'w', 'linewidth', 2)
xlabel('neuron j (sorted by classifier weight)')
ylabel('neuron i (sorted by classifier weight)')
colorbar
title('connection strength from j to i')
%% keyboard

clear x_res
clear proc_simdata
clear peri_event_resps

close all

save([network_param.matsdir 'classifier_' file_tag '_f2f' f2f_type 'discExp_pw' num2str(100*pulse_width)], ...
    'stim_resps', 'stim_y', 'train_inds', 'test_inds', 'c_timestep', 'skip_dt_units', ...
    'classifier_options', 'normClassifier', 'normValidationAccuracy', ...
    'classifier_weights', 'classIdataJ_valAcc')