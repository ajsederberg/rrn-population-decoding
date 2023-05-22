% Selectivity in g400, g600, and spfix networks
% Creates SFig 1. 
file_tags = { 'spnorm_g400_n500_ei100_', 'spnorm_g600_n500_ei100_', 'spfix_g400_n500_ei100_',};
file_labels = {'main results', 'g_{E,I}->1.5*g_{E,I}', '\sigma_{JE,JI} = 0'};
par_colors = [0 0 0; 0.5 0.5 0.5; 0.2 0.5 0.2];

rng_tag = 'rng3';

aeis_filename = @(f_tag) ['results/' f_tag rng_tag ...
    '/matfiles/newclassifierAEIS_' f_tag rng_tag '_f2fdiscExp_pw50.mat'];
par_filename = @(f_tag) ['results/' f_tag rng_tag ...
    '/matfiles/network_paramfile.mat'];
selxcj_filename = @(f_tag) ['results/' f_tag rng_tag ... 
    '/matfiles/selXCJ_stats_' f_tag rng_tag 'f2fdiscExp_pw50.mat'];

sel_res = cellfun(@(f_tag) load(aeis_filename(f_tag)), file_tags, 'UniformOutput', false);
par_res = cellfun(@(f_tag) load(par_filename(f_tag)), file_tags, 'UniformOutput', false);
xcj_res = cellfun(@(f_tag) load(selxcj_filename(f_tag)), file_tags, 'UniformOutput', false);
ii_rows = sort(randperm(500, 25));

%
J_mats = cellfun(@(x) x.network_param.J, par_res, 'UniformOutput', false);
J_lim = [-1 0.5]*max(cellfun(@(x) max(abs(x(:))), J_mats));

sel_cellgps = [1 3 4];
% compute the mean accuracy on the testset trials across draws of the
% test/train datasets, also compute the std dev
mean_testAcc = cell2mat(cellfun(@(x) cellfun(@(y) mean(y.cvDraws_testAcc), x.cls_res(sel_cellgps)), ... 
    sel_res, 'UniformOutput', false)');
sd_testAcc = cellfun(@(x) cellfun(@(y) std(y.cvDraws_testAcc), x.cls_res(sel_cellgps)), ... 
    sel_res, 'UniformOutput', false);

% compute the mean accuracy on a single test set trials drawing over
% different sets of cells (of different sizes)
mean_setsTestAcc = cellfun(@(x) cell2mat(cellfun(@(y) mean(y.cvSetDraws_testAcc,1)', ...
    x.cls_res(sel_cellgps), 'uniformoutput', false)), ... 
    sel_res, 'UniformOutput', false);
sd_setsTestAcc = cellfun(@(x) cell2mat(cellfun(@(y) std(y.cvSetDraws_testAcc,[],1)', ...
    x.cls_res(sel_cellgps), 'uniformoutput', false)), ... 
    sel_res, 'UniformOutput', false);
% this estimates the overall uncertainty by adding in quadrature
sd_sets_est = cellfun(@(x, y) bsxfun(@(x1,y1) sqrt(x1.^2 + y1.^2), x, y), ... 
    sd_setsTestAcc, sd_testAcc, 'UniformOutput', false);

% later, need SE of mean accuracy as a matrix for plotting
sd_testAcc = cell2mat(sd_testAcc');
% extract number of cells in each of the repetitions of classifer fitting
% for different size cell sets
numCells = cellfun(@(x) cell2mat(cellfun(@(y) cellfun(@(z) length(z), y.cvSetDraws_inds(1, :))', ...
    x.cls_res(sel_cellgps), 'uniformoutput', false)), sel_res, 'UniformOutput', false);

% compute stim1 and stim2 firing rates
stim1_fr = cellfun(@(x) mean(x.stim1_endPoint, 2), xcj_res, 'UniformOutput', false);
stim2_fr = cellfun(@(x) mean(x.stim2_endPoint, 2), xcj_res, 'UniformOutput', false);


cellGps = cellfun(@(x) x.cellGroup, sel_res{1}.cls_res(sel_cellgps), 'uniformoutput', false);
%% all nets

rng_tags = {'rng1'; 'rng2'; 'rng3'; 'rng4'; 'rng5'; 'rng6'; 'rng7'; 'rng8'; 'rng9'; 'rng10'; 'rng11'; 'rng12'; 'rng13'; 'rng17'};
n_f = length(file_tags); 
n_r = length(rng_tags);

aeis_i_filename = @(f_tag, r_tag) ['results/' f_tag r_tag ...
    '/matfiles/newclassifierAEIS_' f_tag r_tag '_f2fdiscExp_pw50.mat'];
par_i_filename = @(f_tag, r_tag) ['results/' f_tag r_tag ...
    '/matfiles/network_paramfile.mat'];

all_sel_res = cellfun(@(f_tag, r_tag) load(aeis_i_filename(f_tag, r_tag)), ...
    repmat(file_tags, n_r, 1), repmat(rng_tags, 1, n_f), 'UniformOutput', false);
all_par_res = cellfun(@(f_tag, r_tag) load(par_i_filename(f_tag, r_tag)), ...
    repmat(file_tags, n_r, 1), repmat(rng_tags, 1, n_f), 'UniformOutput', false);

mean_Acc = cellfun(@(x) mean(x.cls_res{1}.cvDraws_testAcc), all_sel_res);
se_Acc = cellfun(@(x) std(x.cls_res{1}.cvDraws_testAcc), all_sel_res);


%% generate plots
lines_cmap = lines(4);
stim_colors = lines_cmap(3:4, :);
par_ind = {'-i', '-ii', '-iii'};
makeMyFigure(20, 30);

for i_par = 1:length(file_tags)
    subplot(6, 3, i_par) 
    imagesc(J_mats{i_par}(ii_rows,ii_rows), 0.4*[-1 1])
    axis square
    title({'J (partial)'; file_labels{i_par} }, ... 
        'color', par_colors(i_par, :))
    ch = colorbar();
    ch.Ticks = [-0.4 0 0.4];
    xlabel('pre')
    ylabel('post')
    set(gca, 'xtick', [1 length(ii_rows)], 'ytick', [1 length(ii_rows)])
    text(-0.5, 1.3, ['A' par_ind{i_par}], 'FontSize', 12, 'units', 'normalized');
%     ch.Position %
    ch.Position = [ch.Position(1) + 0.02 0.8386 0.0101 0.0783]; 

    subplot(6, 3, i_par + 3)
    hold on
    histogram(nonzeros(J_mats{i_par}), 'Normalization', 'pdf')
    xlabel('synaptic weights')
    ylabel('pdf')
    xlim(J_lim)
    set(gca, 'color', 'none')
    axis square
    y_lim = ylim;
    text(-0.5, 1.2, ['B' par_ind{i_par}], 'FontSize', 12, 'units', 'normalized');
    title('nonzero J')
    
    subplot(6, 3, i_par + 6)
    hold on
    histogram(stim1_fr{i_par}, 'Normalization', 'pdf', 'FaceColor', stim_colors(1, :))
    histogram(stim2_fr{i_par}, 'Normalization', 'pdf', 'FaceColor', stim_colors(2, :))
%     hold on
%     histogram(nonzeros(J_mats{i_par}), 'Normalization', 'pdf')
%     xlabel('synaptic weight')
%     ylabel('pdf')
%     xlim(J_lim)
    xlabel('pop ave rate')
    ylabel('pdf')
    title('population FR')
    set(gca, 'color', 'none')
    axis square
    text(-0.5, 1.2, ['C' par_ind{i_par}], 'FontSize', 12, 'Units', 'normalized');

end


sh1 = subplot(10, 2, [13 15]);
errorbar(1:3, mean_testAcc(:, 1), sd_testAcc(:, 1), 'o')
set(gca, 'xtick', 1:3, 'xticklabel', file_labels)
ylim([0.5 0.9])
xlim([0.5 3.5])
ylabel('decoding accuracy')
title([ cellGps{1} ' cells'])
set(gca, 'color', 'none')
text(-0.1, 1.1, 'D', 'FontSize', 12, 'Units', 'normalized');

sh2 = subplot(10, 2, [14 16]);
errorbar(bsxfun(@plus, (1:3)', [0 0.2]), mean_testAcc(:, [2 3]), sd_testAcc(:, [2 3]), 'o')
set(gca, 'xtick', 1:3, 'xticklabel', file_labels)
ylim([0.5 0.9])
xlim([0.5 3.5])
ylabel('decoding accuracy')
title(['E-cells and I-cells'])
set(gca, 'color', 'none')
legend({'E-cells', 'I-cells'}, 'location', 'southeast')
text(-0.1, 1.1, 'E', 'FontSize', 12, 'Units', 'normalized');

sh1.Position(2) = 0.3;
sh2.Position(2) = 0.3;


subplot(5, 1, 5)
eh = errorbar(mean_Acc, se_Acc, 'o');
assignColorsToLines(eh, par_colors);
ylim([0.5 0.9])
xlim([0.5 size(mean_Acc, 1) + 0.5])
ylabel('decoding accuracy')
xlabel('network #')
set(gca, 'xtick', 1:size(mean_Acc, 1))
title('all-cell decoding for all network topologies in main results')
set(gca, 'color', 'none')
text(-0.1, 1.1, 'F', 'FontSize', 12, 'Units', 'normalized');


% cg_label = {'G', 'H', 'I'};
% for i_cg = 1:3
%     subplot(6, 3, 15+i_cg)
%     hold on
%     for i_par = 1:3
%         eh = errorbar(numCells{i_par}(:, i_cg), mean_setsTestAcc{i_par}(:, i_cg), ... 
%             sd_setsTestAcc{i_par}(:, i_cg), 'o-');
%         set(eh, 'color', par_colors(i_par, :), ... 
%             'linewidth', 1)
%     end
%     title([ cellGps{i_cg} ' cells'])
%     xlabel('# of cells drawn')
%     ylabel('decoding accuracy')
%     ylim([0.6 0.9])
%     set(gca, 'color', 'none')
%     text(-0.1, 1.2, cg_label{i_cg}, 'FontSize', 12, 'Units', 'normalized');
% 
% 
% end
print(gcf, '-dpdf', 'writeup/figures/sfig1_params_v1')

%%
figure()
errorbar(mean_Acc, se_Acc, 'o-')
