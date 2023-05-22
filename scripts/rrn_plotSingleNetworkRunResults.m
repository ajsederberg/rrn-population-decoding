% Script to plot results from single experiment
sim_name = 'spnormE_g300_n500_ei100_rng3';
load(['results/' sim_name '/matfiles/selXCJ_stats_' ...
    sim_name 'f2fdiscExp_pw50.mat']);
load(['results/' sim_name '/matfiles/newclassifierAEIS_' ...
    sim_name '_f2fdiscExp_pw50.mat'])

%%
makeMyFigure(30, 20);
auc_bins = linspace(0, 1, 21);

subplot(221)
hold on
histogram(cell_auroc_val(~inh_cells), auc_bins, 'Normalization', 'pdf')
histogram(cell_auroc_val(inh_cells), auc_bins, 'Normalization', 'pdf')
y_lim = ylim;
plot([exc_2tail; exc_2tail], [y_lim' y_lim'], 'k', 'linewidth', 1.5)
xlim([0.2 0.8])
legend({'excitatory', 'inhibitory'})
xlabel('AUC')
ylabel('pdf')

subplot(4, 2, 2)
hold on
for ii = 1:2
    errorbar(ii, ave_choice_EI(ii), se_choice_EI(ii), 'o')
end
ylim([0 0.2])
xlim([0.5 2.5])
ylabel('ave |AUC - 0.5|')
set(gca, 'xtick', [1 2], 'xticklabel', {'exc', 'inh'})

subplot(4, 2, 4)
saOp_fields = {'sameSel_etoi_J_vals', 'opp_Sel_etoi_J_vals', ...
    'sameSel_itoe_J_vals', 'opp_Sel_itoe_J_vals'};
dy_val = [0 0 0 0];
for ii = 1:4
    dy_val(ii) = eval(['sqrt(sum(' saOp_fields{ii} '~=0))/length(' saOp_fields{ii} ')']);
end
p_conn = mean(network_param.J(:) ~= 0);
errorbar(1:4, probConnect_SaOp, dy_val, 'o')
hold on
plot([1 4], p_conn*[ 1 1], 'k')
ylabel('p(connection)')
set(gca, 'xtick', 1:4, 'XTickLabel', probConnect_SaOp_labels)

subplot(2, 2, 3)
hold on
plot(diag(cls_res{1}.classIdataJ_valAcc), 'k')
plot(diag(cls_res{3}.classIdataJ_valAcc), 'b')
plot(diag(cls_res{4}.classIdataJ_valAcc), 'r')
ylabel('classifier accuracy')
xlabel('time (stim starts after 0, ends by 50)')

subplot(4, 2, 6)
i_tc = 10;
w = cls_res{1}.classifier_weights(i_tc, :);
[~, wt_ord] = sort(abs(w), 'descend');
inh_ord = inh_cells(wt_ord);
hold on
stem(find(~inh_ord), w(wt_ord(~inh_ord)), '.')
stem(find(inh_ord), w(wt_ord(inh_ord)), '.')
ylim([-1 1])
ylabel('weights')

subplot(4, 2, 8)
hold on
for ii = [1 2 3 4]
    errorbar(ii, mean(cls_res{ii}.cvDraws_testAcc), std(cls_res{ii}.cvDraws_testAcc), 'o')
end
labels = cellfun(@(x) x.cellGroup, cls_res, 'UniformOutput', false);
set(gca, 'xtick', 1:4, 'XTickLabel', labels)
ylim([0.5 1])
% hold on
% histogram(cls_res{1}.classifier_weights(i_tc, ~inh_cells), 'Normalization', 'cdf', 'DisplayStyle', 'stairs')
% histogram(cls_res{1}.classifier_weights(i_tc, inh_cells), 'Normalization', 'cdf', 'DisplayStyle', 'stairs')
% xlim(0.3*[-1 1])