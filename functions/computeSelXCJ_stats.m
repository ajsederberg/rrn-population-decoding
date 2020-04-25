function selXCJ_res = computeSelXCJ_stats(proc_simdata, network_param, input_type, pulse_width, keep_full_xc_vs_t)

% taking this at the very end means that hte high-rate stimulus tends to
% have more inputs in the window. So, we take it toward the end, but before
% the last 1/f interval. (50 tau, 8 stims, so before the last ~7 tau. 
preChoice_bin = proc_simdata.peri_stim_time > 38 & proc_simdata.peri_stim_time <= 43; 
all_endPoint = mean(proc_simdata.dr_out_peristim(:, :, preChoice_bin), 3);


stim_is1 = squeeze(proc_simdata.stimulus_info_array(:, 1, 1) == 1);

stim1_endPoint = all_endPoint(stim_is1 == 1, :);
stim2_endPoint = all_endPoint(stim_is1 ~= 1, :);

num_neur = size(network_param.J, 1);

inh_cells = all(network_param.J <= 0, 1)';
if  ~(contains(network_param.network_tag, 'spnormE') || contains(network_param.network_tag, 'spEnoise'))
    % if both e and i cells receive inptus, don't keep the 'input' units
    postL1_cells = network_param.input_pattern == 0;
else
    postL1_cells = true(size(network_param.input_pattern));
end
%% START WORKING HERE (11/17) 
cell_auroc_val = zeros(num_neur, 1);

for i_cell = 1:num_neur
    cell_auroc_val(i_cell) = auroc(stim1_endPoint(:, i_cell), stim2_endPoint(:, i_cell), 0);

end
%% old version
% num_shufs = 100;
% cell_auroc_shuf0 = zeros(num_neur, num_shufs);
% 
% for i_shuf = 1:num_shufs
%     for i_cell = 1:num_neur
%      
%         shuf1 = randperm(length(stim_is1), sum(stim_is1 == 1));
%         shuf2 = setdiff(1:length(stim_is1), shuf1);
%         cell_auroc_shuf0(i_cell, i_shuf) = auroc(all_endPoint(shuf1, i_cell), all_endPoint(shuf2, i_cell), 0);       
%     end
% end
%% new version

% take only 200 trials per stimulus 
num_shufs = 100;
num_trials_draw = 100;
cell_auroc_shuf = zeros(num_neur, num_shufs);

for i_shuf = 1:num_shufs
    for i_cell = 1:num_neur
     
        shufTrials = ceil(length(stim_is1)*rand(2*num_trials_draw, 1));
        shuf1 = shufTrials(1:num_trials_draw); % randperm(length(stim_is1), sum(stim_is1 == 1));
        shuf2 = shufTrials(1+num_trials_draw: end); %setdiff(1:length(stim_is1), shuf1);
        cell_auroc_shuf(i_cell, i_shuf) = auroc(all_endPoint(shuf1, i_cell), all_endPoint(shuf2, i_cell), 0);       
    end
end
%%

% First subplot: distribution of AUROC values across cells
shuf_inh_aurocs = cell_auroc_shuf(inh_cells, :);
shuf_exc_aurocs = cell_auroc_shuf(~inh_cells, :);
prctl_2tail = [2.5 97.5]; %[.5 99.5] ; %[2.5 97.5];

inh_2tail = prctile(shuf_inh_aurocs(:), prctl_2tail);
exc_2tail = prctile(shuf_exc_aurocs(:), prctl_2tail);
% inh_2tail = [0.35 0.65]; %prctile(shuf_inh_aurocs(:), prctl_2tail);
% exc_2tail = [0.38 0.62]; %prctile(shuf_exc_aurocs(:), prctl_2tail);
%%
% determine cell selectivity
all_cell_selectivity = zeros(num_neur, 1);
num_excL2 = sum(~inh_cells & postL1_cells);
num_inh = sum(inh_cells);
% coded: -1 means significant AUROC < 0.5 and 
% +1 means significant AUROC > 0.5
all_cell_selectivity(~inh_cells) = -1*(cell_auroc_val(~inh_cells) < exc_2tail(1)) + ... 
    1*(cell_auroc_val(~inh_cells) > exc_2tail(2));
all_cell_selectivity(inh_cells) = -1*(cell_auroc_val(inh_cells) < inh_2tail(1)) + ... 
    1*(cell_auroc_val(inh_cells) > inh_2tail(2));
% compute fraction of cells that are selective
frac_sel_EI = [mean(all_cell_selectivity(~inh_cells & postL1_cells) ~= 0) ... 
    mean(all_cell_selectivity(inh_cells & postL1_cells) ~= 0)];
sd_num_sel_EI =  sqrt([sum(all_cell_selectivity(~inh_cells & postL1_cells) ~= 0) ... 
    sum(all_cell_selectivity(inh_cells & postL1_cells) ~= 0)]);
sd_frac_sel_EI = sd_num_sel_EI./[num_excL2 num_inh];
% compute average choice selectivity (|AUC - 0.5|) ** note in Najafi paper,
% this is 2*|AUC-0.5|, so multiply by 2 for comparison. 
ave_choice_EI = [mean(abs(cell_auroc_val(~inh_cells & postL1_cells)-0.5)) ...
    mean(abs(cell_auroc_val(inh_cells)-0.5))];
se_choice_EI = [std(abs(cell_auroc_val(~inh_cells & postL1_cells)-0.5)) ...
    std(abs(cell_auroc_val(inh_cells)-0.5))]./sqrt([num_excL2 num_inh]-1);
%%
% Compute the PDF of single-cell auroc values
auroc_bin_edges = linspace(0.2, 0.8, 13);
dbin = mean(diff(auroc_bin_edges));
auroc_bins = edges2bins(auroc_bin_edges)';
inh_auroc_pdf = histcounts(cell_auroc_val(inh_cells), auroc_bin_edges, 'Normalization', 'pdf');
exc_auroc_pdf = histcounts(cell_auroc_val(~inh_cells & postL1_cells), auroc_bin_edges, 'Normalization', 'pdf');

% Compute NOISE CORRELATIONS
[xc_vs_time, xc_vs_time_hist]  = simdataNoiseCorr(proc_simdata);
%% take a slice, 
xc_time_window = proc_simdata.peri_stim_time > 10 & proc_simdata.peri_stim_time <= 40; %100:400; % 100:400 is chosen for this network b/c cross corrs stabilize there


mean_xc_vals = squeeze(mean(xc_vs_time(:, :, xc_time_window), 3));

if ~keep_full_xc_vs_t
    xc_time = 1:10:size(xc_vs_time, 3);
    xc_vs_time = xc_vs_time(:, :, 1:10:end);
end

% separate by: E-I pairs with same and different selectivity. 
etoiPairs_logical = inh_cells & (~inh_cells & postL1_cells)';
itoePairs_logical = etoiPairs_logical';
% also: e-e and i-i pairs with same and different selectivity
etoePairs_logical = (~inh_cells & postL1_cells) & (~inh_cells & postL1_cells)';
itoiPairs_logical = inh_cells & inh_cells';

% id with same selectivity
pairs_sameSel = all_cell_selectivity*all_cell_selectivity' == 1;
pairs_oppSel = all_cell_selectivity*all_cell_selectivity' == -1;

% compute same/oppSel and cross-cor
etoi_pairs_sameSel_xc = mean_xc_vals(etoiPairs_logical & pairs_sameSel);
etoi_pairs_oppSel_xc = mean_xc_vals(etoiPairs_logical & pairs_oppSel);
itoe_pairs_sameSel_xc = mean_xc_vals(itoePairs_logical & pairs_sameSel);
itoe_pairs_oppSel_xc = mean_xc_vals(itoePairs_logical & pairs_oppSel);
etoe_pairs_sameSel_xc = mean_xc_vals(etoePairs_logical & pairs_sameSel);
etoe_pairs_oppSel_xc = mean_xc_vals(etoePairs_logical & pairs_oppSel);
itoi_pairs_sameSel_xc = mean_xc_vals(itoiPairs_logical & pairs_sameSel);
itoi_pairs_oppSel_xc = mean_xc_vals(itoiPairs_logical & pairs_oppSel);


mean_xc_SaOp_etoi = [mean(etoi_pairs_sameSel_xc) mean(etoi_pairs_oppSel_xc)];
sem_xc_SaOp_etoi  = [std(etoi_pairs_sameSel_xc)  std(etoi_pairs_oppSel_xc)]/sqrt(num_neur-1);
mean_xc_SaOp_itoe = [mean(itoe_pairs_sameSel_xc) mean(itoe_pairs_oppSel_xc)];
sem_xc_SaOp_itoe  = [std(itoe_pairs_sameSel_xc)  std(itoe_pairs_oppSel_xc)]/sqrt(num_neur-1);
mean_xc_SaOp_etoe = [mean(etoe_pairs_sameSel_xc) mean(etoe_pairs_oppSel_xc)];
sem_xc_SaOp_etoe  = [std(etoe_pairs_sameSel_xc)  std(etoe_pairs_oppSel_xc)]/sqrt(num_neur-1);
mean_xc_SaOp_itoi = [mean(itoi_pairs_sameSel_xc) mean(itoi_pairs_oppSel_xc)];
sem_xc_SaOp_itoi  = [std(itoi_pairs_sameSel_xc)  std(itoi_pairs_oppSel_xc)]/sqrt(num_neur-1);

% connectivity between EI pairs with same vs. opposite selectivity
sameSel_etoi_J_vals = network_param.J(etoiPairs_logical & pairs_sameSel);
opp_Sel_etoi_J_vals = network_param.J(etoiPairs_logical & pairs_oppSel);
sameSel_itoe_J_vals = network_param.J(itoePairs_logical & pairs_sameSel);
opp_Sel_itoe_J_vals = network_param.J(itoePairs_logical & pairs_oppSel);

% connectivity between EE and II pairs with same vs. opposite selectivity
sameSel_etoe_J_vals = network_param.J(etoePairs_logical & pairs_sameSel);
opp_Sel_etoe_J_vals = network_param.J(etoePairs_logical & pairs_oppSel);
sameSel_itoi_J_vals = network_param.J(itoiPairs_logical & pairs_sameSel);
opp_Sel_itoi_J_vals = network_param.J(itoiPairs_logical & pairs_oppSel);


probConnect_SaOp_labels = {'E->I, same', 'E->I, diff', 'I->E, same', 'I->E, diff', ...
    'E->E, same', 'E->E, diff', 'I->I, same', 'I->I, diff'};
probConnect_SaOp = [mean(sameSel_etoi_J_vals ~= 0) mean(opp_Sel_etoi_J_vals ~= 0) ...
    mean(sameSel_itoe_J_vals ~= 0) mean(opp_Sel_itoe_J_vals ~= 0)...
    mean(sameSel_etoe_J_vals ~= 0) mean(opp_Sel_etoe_J_vals ~= 0) ...
    mean(sameSel_itoi_J_vals ~= 0) mean(opp_Sel_itoi_J_vals ~= 0)];

numConnect_SaOp = [sum(sameSel_etoi_J_vals ~= 0) sum(opp_Sel_etoi_J_vals ~= 0) ...
    sum(sameSel_itoe_J_vals ~= 0) sum(opp_Sel_itoe_J_vals ~= 0)...
    sum(sameSel_etoe_J_vals ~= 0) sum(opp_Sel_etoe_J_vals ~= 0) ...
    sum(sameSel_itoi_J_vals ~= 0) sum(opp_Sel_itoi_J_vals ~= 0)];


wt_dists{1} = nonzeros(sameSel_etoi_J_vals);
wt_dists{2} = nonzeros(opp_Sel_etoi_J_vals);
weightConnect_SaOp = [mean(nonzeros(sameSel_etoi_J_vals)) mean(nonzeros(opp_Sel_etoi_J_vals))];

%% number of shared pre-synaptic cells
[shared_preSyn_count, corr_preSyn] = countSharedPreSynapticCells(network_param.J);
%%
sameSel_sharedPre = shared_preSyn_count(etoiPairs_logical & pairs_sameSel);
oppSel_sharedPre = shared_preSyn_count(etoiPairs_logical & pairs_oppSel);
sameSel_corrPre = corr_preSyn(etoiPairs_logical & pairs_sameSel);
oppSel_corrPre = corr_preSyn(etoiPairs_logical & pairs_oppSel);

clear proc_simdata
stimulus_experiment_tag = [input_type 'discExp_pw' num2str(100*pulse_width)];

file_name = [network_param.save_dir 'matfiles/selXCJ_stats_' network_param.network_tag stimulus_experiment_tag];
save(file_name);
selXCJ_res = load(file_name);
