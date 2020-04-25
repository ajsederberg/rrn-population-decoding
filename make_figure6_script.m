% Script to make Figure 5: changes in selectivity under changes in task
% parameters. 
num_neur = 500;
gain_neur = 400;
ei_ratio = 100; % 92 or 108;  80 or 125
connect_tag = 'spnorm'; % 'spnorm2P' for two-pool

i_set = 1;  % choose 1,2, 3, 4 for all/exc/excISM/inh

rng_list = [1:13 17];

all_cls_res = cell(1, length(rng_list));
all_cls_names = cell(1, length(rng_list));
all_selXCJ_res = cell(1, length(rng_list));
all_selXCJ_names = cell(1, length(rng_list));
all_selXCJ = cell(1, length(rng_list));
%%
for i_rv = 1:length(rng_list)
    rng_val = rng_list(i_rv);
file_tag = [connect_tag '_g' num2str(gain_neur) '_n' num2str(num_neur) '_ei' num2str(ei_ratio) '_rng' num2str(rng_val)];

local_res_dir = '/Volumes/Samsung_T5/projects/random_rec_network_project/results/'; % classifier results are here

local_mats_dir = [local_res_dir file_tag '/matfiles/'];
%% 
file_list = dir(local_mats_dir);
is_cAEIS = arrayfun(@(x) strcmp(x.name(1), 'n') && contains(x.name, 'newclassifierAEIS'), file_list);
classifier_res = file_list(is_cAEIS);

all_cls_res{i_rv} = cell(length(classifier_res), 1);
all_cls_names{i_rv} = arrayfun(@(x) x.name, classifier_res, 'UniformOutput', false);
for i_acr = 1:length(classifier_res)
    
    x_in = load([local_mats_dir classifier_res(i_acr).name]);
    all_cls_res{i_rv}{i_acr} = x_in.cls_res{i_set}; % this only loads the decoding results from the full population
end
%%


%%
is_selXCJ = arrayfun(@(x) strcmp(x.name(1), 's') && contains(x.name, 'selXCJ'), file_list);
selXCJ_res = file_list(is_selXCJ);

all_selXCJ_res{i_rv} = cell(length(selXCJ_res), 1);
all_selXCJ_names{i_rv} = arrayfun(@(x) x.name, selXCJ_res, 'UniformOutput', false);
for i_acr = 1:length(selXCJ_res)
    all_selXCJ_res{i_rv}{i_acr} = load([local_mats_dir selXCJ_res(i_acr).name], 'cell_auroc_val', ...
        'inh_2tail', 'exc_2tail', 'auroc_bins', 'auroc_bin_edges', 'exc_auroc_pdf', 'inh_auroc_pdf', 'dbin', ...
        'probConnect_SaOp_labels', 'probConnect_SaOp', 'all_cell_selectivity', 'postL1_cells');
end
%%

all_selXCJ{i_rv} = cell2mat(cellfun(@(x) x.cell_auroc_val', all_selXCJ_res{i_rv}, 'UniformOutput', false));
end


all_cls_name_tags = cellfun(@(x) cellfun(@(y) regexp(y, '_f(\w*)fdiscExp', 'tokens'), ...
    x, 'uniformoutput', false), all_cls_names, 'UniformOutput', false);
all_cls_name_tags = cellfun(@(x) cellfun(@(y) y{1}{1}, x, 'uniformoutput', false), ...
    all_cls_name_tags, 'uniformoutput', false);


%% compute fraction selective (-1, 0, 1)

frac_sel = cellfun(@(x) cell2mat(cellfun(@(y) histcounts(y.all_cell_selectivity(y.postL1_cells), -1.5:1:1.5, 'Normalization', 'probability'), ... 
    x, 'UniformOutput', false)), all_selXCJ_res, 'uniformoutput', false);

isPostL1 = cellfun(@(y) cell2mat(cellfun(@(x) x.postL1_cells, y', 'uniformoutput', false)), all_selXCJ_res, 'uniformoutput', false);

has_three = cellfun(@(x) size(x, 1) >= 3, all_selXCJ_res);
has_two = cellfun(@(x) size(x, 1) >= 2, all_selXCJ_res);

use_sims = has_two;
exp_pair = [1 2];
sel_crossTab = cellfun(@(x, inds_l2) ...
    crosstab(x{exp_pair(1)}.all_cell_selectivity(inds_l2(:,exp_pair(1)) ), ...
    x{exp_pair(2)}.all_cell_selectivity(inds_l2(:,exp_pair(1))))/sum(inds_l2(:, 1)), ... 
    all_selXCJ_res(use_sims), isPostL1(use_sims), 'uniformoutput', false);

sel1_name = cellfun(@(x) x{exp_pair(1)},all_selXCJ_names(use_sims)', 'UniformOutput', false);
sel2_name = cellfun(@(x) x{exp_pair(2)},all_selXCJ_names(use_sims)', 'UniformOutput', false);

% check that sel1 and sel2 names are all the same by viewing them

ave_crossTab = reshape(mean(cell2mat(cellfun(@(x) x(:), sel_crossTab, 'UniformOutput', false)), 2), [3 3]);

ex_inds = find(use_sims);

%% compute over decoding accuracy

dec_acc = cellfun(@(x) cellfun(@(y) mean(y.cvDraws_testAcc), x), all_cls_res, 'UniformOutput', false);
se_dec_acc = cellfun(@(x) cellfun(@(y) std(y.cvDraws_testAcc), x), all_cls_res, 'UniformOutput', false);

dec_acc3 = cell2mat(dec_acc(has_three));
dec_acc2 = cell2mat(cellfun(@(x) x(1:2, :), dec_acc, 'UniformOutput', false));

se_dec_acc3 = cell2mat(se_dec_acc(has_three));
se_dec_acc2 = cell2mat(cellfun(@(x) x(1:2, :), se_dec_acc, 'UniformOutput', false));

%% compare weights after dividing out the value of the bias (put all on the same scale. 

% this divides out the bias b (* warning * if biases flip sign, this will
% not average over the same classifier rule, so check when function is
% called) and then multiplies by the average bias. 
% in practice, most of the biases are very close tot he same thing, and
% this has little effect (whether you divide out the bias or not). 
wt_arr_fun = @(w, b) mat2cell(diag(1./(b))*squeeze(w)*mean(b), size(w, 1), ones(1, size(w, 3)));
wt_arr_fun0 = @(w, b) mat2cell(diag(1./(1+0*b))*squeeze(w), size(w, 1), ones(1, size(w, 3)));

b_vals = cellfun(@(x) cellfun(@(y) mean(y.cvclassifier_bs), x), ...
    all_cls_res, 'UniformOutput', false);
b_vals_se = cellfun(@(x) cellfun(@(y) std(y.cvclassifier_bs), x), ...
    all_cls_res, 'UniformOutput', false);
%%

% structure of this is: wt_arrs{i}{j}{k} is the i-th simulation (network
% topology) with the j-th frequency test (7 vs 14, etc.) and the k-th cell.
% 
wt_arrs = cellfun(@(x) cellfun(@(y) wt_arr_fun(...
    y.cvclassifier_weights(y.cvclassifier_bs > 0, :, :), y.cvclassifier_bs(y.cvclassifier_bs > 0)), x, 'UniformOutput', false), ...
    all_cls_res, 'UniformOutput', false);
b_arrs = cellfun(@(x) cellfun(@(y) y.cvclassifier_bs(y.cvclassifier_bs > 0), x, 'UniformOutput', false), ...
    all_cls_res, 'UniformOutput', false);
%% ranksum of weights for each neuron in the second and third (or first and second if only two)
wt_ranksum23 = cellfun(@(y) cellfun(@(z1, z2) ranksum(z1, z2), y{end-1}, y{end}),  wt_arrs, 'UniformOutput', false);
% ranksum of weights for each neuron in the first and third experiments (or first and second if only two)
wt_ranksum13 = cellfun(@(y) cellfun(@(z1, z2) ranksum(z1, z2), y{1}, y{end}),  wt_arrs, 'UniformOutput', false);
% ranksum of weights for each neuron in the first and second experiments
wt_ranksum12 = cellfun(@(y) cellfun(@(z1, z2) ranksum(z1, z2), y{1}, y{2}),  wt_arrs, 'UniformOutput', false);
%%
wt_ranksum11 = cellfun(@(y) cellfun(@(z1, z2) ranksum(z1(1:round(end/2)), z2(round(end/2)+1:end)), y{1}, y{1}),  wt_arrs, 'UniformOutput', false);

[~, wt_ords] = cellfun(@(x) sort(x), wt_ranksum12, 'UniformOutput', false);
wts_rs_12 = cell2mat(cellfun(@(x, ord) x(ord), wt_ranksum12(use_sims), wt_ords(use_sims), 'UniformOutput', false)');
wts_rs_23 = cell2mat(cellfun(@(x, ord) x(ord), wt_ranksum23(use_sims), wt_ords(use_sims), 'UniformOutput', false)');
wts_rs_13 = cell2mat(cellfun(@(x, ord) x(ord), wt_ranksum13(use_sims), wt_ords(use_sims), 'UniformOutput', false)');


% mean/se weights in each set
mean_cls_wts = cellfun(@(x) cellfun(@(y) ...
    cellfun(@(z) mean(z), y), x, 'UniformOutput', false), wt_arrs, 'UniformOutput', false);
se_cls_wts = cellfun(@(x) cellfun(@(y) ...
    cellfun(@(z) std(z), y), x, 'UniformOutput', false), wt_arrs, 'UniformOutput', false);
wts_nz = cellfun(@(x1, x2) cellfun(@(y1, y2) abs(y1) - 3*y2 > 0, x1, x2, 'UniformOutput', false), ...
    mean_cls_wts, se_cls_wts, 'UniformOutput', false);
any_exp_wt_nz = cellfun(@(x) any(cell2mat(x), 1), wts_nz, 'UniformOutput', false);

% fraction of cells whose weights changed
p_rs_thr = 0.001;
frac_nz_chg_12 = cellfun(@(x, is_nz) mean(x(is_nz) < p_rs_thr), wt_ranksum12, any_exp_wt_nz);
frac_nz_chg_13 = cellfun(@(x, is_nz) mean(x(is_nz) < p_rs_thr), wt_ranksum13, any_exp_wt_nz);
frac_nz_chg_23 = cellfun(@(x, is_nz) mean(x(is_nz) < p_rs_thr), wt_ranksum23, any_exp_wt_nz);

frac_all_chg_12 = cellfun(@(x, is_nz) mean(x < p_rs_thr), wt_ranksum12, any_exp_wt_nz);
frac_all_chg_13 = cellfun(@(x, is_nz) mean(x < p_rs_thr), wt_ranksum13, any_exp_wt_nz);
frac_all_chg_23 = cellfun(@(x, is_nz) mean(x < p_rs_thr), wt_ranksum23, any_exp_wt_nz);

%%
ex_ind = 1;
x = mean_cls_wts{ex_ind}{2};
dx = se_cls_wts{ex_ind}{2};
y = mean_cls_wts{ex_ind}{1};
dy = se_cls_wts{ex_ind}{1};
% these labels are correct if x is from exp 2 ('f2f') and y is from exp 1
% ('f1020f')
exp_labels = {'8 Hz vs. 16 Hz (original)', '10 Hz vs. 20 Hz (new)'};

num_sigs = 2;

dh = sqrt(dx.^2 + dy.^2);
[~, ord] = sort(x);

plot_x = x(ord);
plot_y = y(ord);
plot_dh = dh(ord);

nonzero_x = abs(plot_x)./dx(ord) > num_sigs & abs(plot_x) > 0.05;
nonzero_y = abs(plot_y)./dy(ord) > num_sigs & abs(plot_y) > 0.05;
sig_diff = abs(plot_y - plot_x)./plot_dh;
plot_bold = sig_diff > num_sigs;

lines_cmap = lines(7);
exp_colors = [0.5*[1 1 1]; lines_cmap(4, :) ];
%%

plot_DE = false;

if plot_DE
    hFig = makeMyFigure(8.5*2.54*(11/7), 2*2.25*2.54*11/7);
else
    
    hFig = makeMyFigure(8.5*2.54*(11/7), 2.25*2.54*11/7);
end
% set color standards
stim_colors = lines(4);
stim_colors = stim_colors(3:4, :);
neuron_colors = [(1 + lines(1))/2; lines(2)];
fsz_labels = 16;
fsz_axis = 12;


if ~plot_DE
    a_w = 0.18;
bcd_w = 0.18;
y_pos = 0.2;
h_val = 0.6;
gap_bcd = 0.12;
x_b = 2*gap_bcd + a_w;
x_c = x_b + gap_bcd + bcd_w;



hA = axes(hFig, 'Position', [gap_bcd y_pos a_w h_val ]);
hB = axes(hFig, 'Position', [x_b y_pos bcd_w h_val]);
hC = axes(hFig, 'Position', [x_c y_pos bcd_w h_val]);

else 
    
    a_w = 0.18;
bcd_w = 0.18;
y_pos = 0.6;
h_val = 0.35;
gap_bcd = 0.12;
x_b = 2*gap_bcd + a_w;
x_c = x_b + gap_bcd + bcd_w;



hA = axes(hFig, 'Position', [gap_bcd y_pos a_w h_val ]);
hB = axes(hFig, 'Position', [x_b y_pos bcd_w h_val]);
hC = axes(hFig, 'Position', [x_c y_pos bcd_w h_val]);

    x_d = x_c + gap_bcd + bcd_w;
x_e = x_d + gap_bcd + bcd_w;
h_e = h_val/4;
dh_e = 0.8*h_e;
d_w = 0.3;
y_pos = 0.1;
h_val = 0.35;
gap_def = gap_bcd;
x_e = 2*gap_def + d_w;
x_f = x_e + gap_def + d_w;
% x_d = x_c + gap_def + def_w;
% x_e = x_d + gap_def + def_w;
h_e = h_val/4;
dh_e = 0.8*h_e;

hD = axes(hFig, 'Position', [gap_def y_pos d_w h_val ]);
hE = axes(hFig, 'Position', [x_e y_pos d_w h_val]);
% hF = axes(hFig, 'Position', [x_f y_pos def_w h_val]);

% subplot(1, 3, 1)
set(hFig, 'CurrentAxes', hD);
hold on
eh_dec = errorbar(dec_acc2([2 1], :)', se_dec_acc2([2 1], :)', '.');
assignColorsToLines(eh_dec, exp_colors);
ylim([0.5 1])
set(gca, 'color', 'none', 'FontSize', fsz_axis)

% axis square
lh1 = legend(eh_dec, exp_labels);
lh1.Position(1:2) = [gap_def*2 0.46];
xlabel('network index')
ylabel('population decoding accuracy')
text(-0.1, 1.1, 'D', 'units', 'normalized', 'FontSize', fsz_labels, 'FontWeight', 'bold')
text(-0.2, 0.2, 'population decoding', 'units', 'normalized', 'FontSize', fsz_labels, ... 
    'Rotation', 90)


set(hFig, 'CurrentAxes', hE);
set(gca, 'color', 'none', 'FontSize', fsz_axis)

hold on
eh1 = errorbar(find(~plot_bold), plot_y(~plot_bold), plot_dh(~plot_bold), '.', ...
    'color', (1 + exp_colors(2, :))/2);
eh2 = errorbar(find(plot_bold), plot_y(plot_bold), plot_dh(plot_bold), '.', ...
    'color', exp_colors(2, :), 'markersize', 8, 'linewidth', 1.5);
eh1.CapSize = 0;
eh2.CapSize = 0;
plot(plot_x, 'k.', 'linewidth', 1.5)
% axis square
axis tight
xlabel('cell (ordered by weight, original experiment)')
ylabel('readout weight')
lh = legend({['new weight (< ' num2str(num_sigs, '%1.0f') '\sigma change)'],...
    ['new weight (> ' num2str(num_sigs, '%1.0f') '\sigma change)'], 'original weight'});
lh.Position(1:2) = [gap_def+x_e 0.46];
text(-0.1, 1.1, 'E', 'units', 'normalized', 'FontSize', fsz_labels, 'FontWeight', 'bold')
end


%%% compute # of changed weights (fraction) over all networks

all_nets_x = cellfun(@(x) x{2}, mean_cls_wts, 'UniformOutput', false);
all_nets_dx = cellfun(@(x) x{2}, se_cls_wts, 'UniformOutput', false);
all_nets_y = cellfun(@(x) x{1}, mean_cls_wts, 'UniformOutput', false);
all_nets_dy = cellfun(@(x) x{1}, se_cls_wts, 'UniformOutput', false);

has_changed = cellfun(@(x,y, dx, dy) (abs(x-y)./sqrt(dx.^2 + dy.^2) > num_sigs), ...
    all_nets_x, all_nets_y, all_nets_dx, all_nets_dy, 'UniformOutput', false);

is_nonzero =  cellfun(@(x,y, dx, dy) (abs(x)./dx > 2) | (abs(y)./dy > num_sigs), ...
    all_nets_x, all_nets_y, all_nets_dx, all_nets_dy, 'UniformOutput', false);

num_changed = cellfun(@(is_changed) sum(is_changed), has_changed);
num_changed_nz = cellfun(@(is_changed, nz) sum(is_changed(nz)), has_changed, is_nonzero);
num_wts = cellfun(@(is_changed) length(is_changed), has_changed);
num_nz_wts = cellfun(@(nz) sum(nz), is_nonzero);

dnum_changed = sqrt(num_changed);
frac_changed = num_changed./num_wts;
dfrac_changed = dnum_changed./num_wts;

i_exnet = 1;

net_ID = ex_inds(i_exnet);

for i_e = 1:length(exp_pair)
    i_exp = exp_pair(i_e);
    x_aucstats = all_selXCJ_res{net_ID}{i_exp};
    freqTitleS = regexp(all_selXCJ_names{i_rv}{i_exp}, 'f(\w*)f' , 'tokens');
    
    freqTitle{i_e} = freqTitleS{1}{1};
    auroc_bins{i_e} = x_aucstats.auroc_bins;
    auroc_bin_edges{i_e} = x_aucstats.auroc_bin_edges;
    auc_vals{i_e} = x_aucstats.cell_auroc_val;
    exc_2tail{i_e} = x_aucstats.exc_2tail;
    inh_2tail{i_e} = x_aucstats.inh_2tail;
    dbin{i_e} = x_aucstats.dbin;
    probConnect_SaOp_labels{i_e} = x_aucstats.probConnect_SaOp_labels;
    probConnect_SaOp{i_e} = x_aucstats.probConnect_SaOp;

end  

freqTitle_check = freqTitle;
freqTitle = {'10 Hz vs. 20 Hz', '8 Hz vs. 16 Hz'};

sel_ticks = {'low f', 'non', 'high f'};

% compute the 2d hist
bin_edges = auroc_bin_edges{1}(1):.01:auroc_bin_edges{1}(end);
auc_hist2_cts = histcounts2(auc_vals{1}, auc_vals{2}, bin_edges, ...
    bin_edges);
minmax = @(x) [min(x) max(x)];
auc_lims = minmax(bin_edges(sum(auc_hist2_cts, 1) > 0));


    

% hD = axes(hFig, 'Position', [x_d y_pos bcd_w h_val]);



set(hFig, 'CurrentAxes', hA);

% colormap(hA, flipud(gray));
imagesc(bin_edges, bin_edges, auc_hist2_cts, [0 10])
ch = colorbar;
ch.Label.String = 'counts';
ch.Label.Rotation = 270;
ch.Ticks = [0 10];
ch.TickLabels = {'0', '>10'};
xlabel([freqTitle{2} ' AUC values'])
ylabel([freqTitle{1} ' AUC values'])
xlim(auc_lims)
ylim(auc_lims)
set(gca, 'ydir', 'normal', 'FontSize', fsz_axis)
hold on
plot(auc_lims, auc_lims, 'w', 'linewidth', 1.5)
axis square
title('example network')
text(-0.2, 1.2, 'A', 'units', 'normalized', 'FontSize', fsz_labels, 'FontWeight', 'bold')
text(-0.5, 0, 'single-cell selectivity', 'units', 'normalized', 'FontSize', fsz_labels, ... 
    'Rotation', 90)

% plot the cross-tabs
set(hFig, 'CurrentAxes', hB)
% colormap(hB, flipud(gray))
imagesc(1:3, 1:3, sel_crossTab{i_exnet}, [0 0.3])
ch = colorbar;
ch.Label.String = 'fraction';
ch.Label.Rotation = 270;
ch.Ticks = [0 0.1 0.2 0.3];
ch.TickLabels = {'0','10%','20%', '>30%'};
xlabel([freqTitle{2} ' selectivity'])
ylabel([freqTitle{1} ' selectivity'])
set(gca, 'ydir', 'normal', 'Xtick', 1:3, 'YTick', 1:3, ...
    'XTickLabel', sel_ticks, 'YTickLabel', sel_ticks, 'FontSize', fsz_axis)
axis square
title('example network')
text(-0.2, 1.2, 'B', 'units', 'normalized', 'FontSize', fsz_labels, 'FontWeight', 'bold')

% plot the cross-tabs
set(hFig, 'CurrentAxes', hC)
% colormap(hC, flipud(gray))
imagesc(1:3, 1:3, ave_crossTab, [0 0.3])
ch = colorbar;
ch.Label.String = 'fraction';
ch.Label.Rotation = 270;
ch.Ticks = [0 0.1 0.2 0.3];
ch.TickLabels = {'0','10%','20%', '>30%'};
xlabel([freqTitle{2} ' selectivity'])
ylabel([freqTitle{1} ' selectivity'])
set(gca, 'ydir', 'normal', 'Xtick', 1:3, 'YTick', 1:3, ...
     'XTickLabel', sel_ticks, 'YTickLabel', sel_ticks, 'FontSize', fsz_axis)
axis square
title('average over networks')
text(-0.2, 1.2, 'C', 'units', 'normalized', 'FontSize', fsz_labels, 'FontWeight', 'bold')




print(gcf, '-dpdf', 'writeup/figures/fig6_final')


fraction_change_sel = cellfun(@(x) sum(diag(x, -1))+sum(diag(x, 1)), sel_crossTab);
fraction_same_sel = cellfun(@(x) sum(diag(x)), sel_crossTab);

mean(fraction_change_sel)
mean(fraction_same_sel)
%%
disp(['ave fraction selectivity change : ' num2str(mean(fraction_change_sel), '%1.2f')])
disp(['se fraction selectivity change : ' num2str(std(fraction_change_sel), '%1.2f')])
disp(['ave fraction weights change : ' num2str(mean(frac_changed), '%1.2f')])

