% plot some single-trial results; generates figures for talks
file_tag = 'spnormE_g200_n500_ei100_rng2';
load(['results/' file_tag '/matfiles/sim_data_' ...
    file_tag '_f2fdiscExp_pw50.mat']);

sel_info = load(['results/' file_tag '/matfiles/selXCJ_stats_' ...
    file_tag 'f2fdiscExp_pw50']);
%%
line_cmap = lines(4);
input_colors = line_cmap(3:4, :);


%%

% tt1 = input_params.input_times{1}(1) + proc_simdata.peri_stim_time' - 1;
tt1 = input_params.stim_start_times{1}(1) + proc_simdata.peri_stim_time' - 1;

in1 = input_params.input_fun(tt1);

% tt2 = input_params.input_times{2}(1) + proc_simdata.peri_stim_time' - 1;
tt2 = input_params.stim_start_times{2}(1) + proc_simdata.peri_stim_time' - 1;
in2 = input_params.input_fun(tt2);

figure()
hold on
histogram(in1, 'FaceColor', input_colors(1, :))
histogram(in2, 'FaceColor', input_colors(2, :))

ave1 = squeeze(mean(proc_simdata.stim_ave_dr_out(1, 1, :, :), 4));
ave2 = squeeze(mean(proc_simdata.stim_ave_dr_out(2, 2, :, :), 4));

% [~, ord] = sort(ave1 - ave2);
[~, ord] = sort(sel_info.cell_auroc_val);

cmap = lines(4);
input_colors = cmap(3:4, :);
ei_colors = cmap(1:2, :);

figure()

subplot(3, 1, 1)
ph = plot(tt1, in1, tt1, in2);
assignColorsToLines(ph, input_colors);

subplot(3, 1, 2)
hold on
ph = plot(tt1, squeeze( proc_simdata.dr_out_peristim(1, ord(1), :)), ...
    tt1, squeeze( proc_simdata.dr_out_peristim(2, ord(1), :)))
ph2 = plot(tt1, squeeze( proc_simdata.stim_ave_dr_out(1,1, ord(1), :)), ...
    tt1, squeeze( proc_simdata.stim_ave_dr_out(2,2, ord(1), :)), 'linewidth', 1)

assignColorsToLines(ph, input_colors);
assignColorsToLines(ph2, input_colors);
title('single trial and average, selective for 1')

subplot(3, 1, 3)
hold on
ph = plot(tt1, squeeze( proc_simdata.dr_out_peristim(1, ord(end), :)), ...
    tt1, squeeze( proc_simdata.dr_out_peristim(2, ord(end), :)))
ph2 = plot(tt1, squeeze( proc_simdata.stim_ave_dr_out(1,1, ord(end), :)), ...
    tt1, squeeze( proc_simdata.stim_ave_dr_out(2,2, ord(end), :)), 'linewidth', 1)

assignColorsToLines(ph, input_colors);
assignColorsToLines(ph2, input_colors);

title('single trial, selective for 2')

%% Plot single-cell selectivity 
ex_cell_set = ord([ 20 298 439 500]);
nR = length(ex_cell_set);
makeMyFigure(7,30);

for i_ex = 1:length(ex_cell_set)

    s1_resps = sel_info.stim1_endPoint(:, ex_cell_set(i_ex));
    s2_resps = sel_info.stim2_endPoint(:, ex_cell_set(i_ex));
    
    min_resp = min(min(s1_resps(:)), min(s2_resps(:)));
    max_resp = max(max(s1_resps(:)), max(s2_resps(:)));
    resp_bins = linspace(min_resp, max_resp, 21);

    subplot(nR, 1, i_ex)
    hold on
    histogram(s1_resps, resp_bins, 'FaceColor', input_colors(1, :), ... 
        'normalization', 'pdf')
    histogram(s2_resps, resp_bins, 'FaceColor', input_colors(2, :), ... 
        'normalization', 'pdf')
    
    axis square
    set(gca, 'color', 'none', 'fontsize', 14)
    ylabel('pdf')
    xlabel('activity')
    title(['AUC = ' num2str(sel_info.cell_auroc_val(ex_cell_set(i_ex)), '%1.2f')])
end

print(gcf, '-dpdf', ['plots/talk_plots/single_cell_histogram' file_tag])

%%
makeMyFigure(10,20);
stim_info = proc_simdata.stimulus_info_array(:, 1, 1);

% use_trials = 183:282;
use_trials = 220 + [-19:20];

ind_cell =  ex_cell_set(end);

trial_resps = 0.7*squeeze( proc_simdata.dr_out_peristim(use_trials, ind_cell, :));
[~, trial_ord] = sort(stim_info(use_trials), 'descend'); % this puts stim 1 at beginninngn, stim 2 at ennd
trial_resps = trial_resps(trial_ord, :);

% rescale trial_resps
trial_resps = trial_resps/max(trial_resps(:));

num_s1 = sum(stim_info(use_trials));

plot(tt1 - tt1(1), bsxfun(@plus, (1:size(trial_resps, 1))', trial_resps), 'k')
hold on
choice_time = tt1(sel_info.preChoice_bin) - tt1(1);

plot([choice_time(1); choice_time ; flipud(choice_time)], ...
    0.5 + [num_s1; 0*choice_time ; num_s1 + 0*choice_time], ...
    'color', input_colors(1, :), 'linewidth', 1)
plot([choice_time(1); choice_time ; flipud(choice_time)], ...
    0.5 + [length(use_trials); 0*choice_time + num_s1; length(use_trials) + 0*choice_time], ...
    'color', input_colors(2, :), 'linewidth', 1)

set(gca, 'color', 'none', 'box', 'off', 'fontsize', 14)
axis tight
ylabel('trials')
xlabel('time')
title(['AUC = ' num2str(sel_info.cell_auroc_val(ind_cell), '%1.2f')])

print(gcf, '-dpdf', ['plots/talk_plots/single_cell_raster_ex' num2str(ind_cell)])
print(gcf, '-dsvg', ['plots/talk_plots/single_cell_raster_ex' num2str(ind_cell)])

%% selectivity distributions
is_exc = any(network_param.J > 0, 1);

auc_bins = linspace(0.2, 0.8, 21);
makeMyFigure(10, 10);
hold on
histogram(sel_info.cell_auroc_val(is_exc), auc_bins, 'Normalization', 'pdf')
histogram(sel_info.cell_auroc_val(~is_exc), auc_bins, 'Normalization', 'pdf')




%% Plot population activity projected on eigenvectors of connectivity matrix. 


%% how does connectivity eigenvector correspond to activity covariance eigenvector?
[i_trial,j_stim]=find(proc_simdata.stimulus_info_array==1);

stim1_trials= i_trial(j_stim == 1);
stim2_trials= i_trial(j_stim ~= 1);
%%
[V,D] =eig(network_param.J);
aa = V*D/(V);

figure()
subplot(221)
% have to take real b/c numerical inverse leaves small imaginary parts
imagesc(real(aa))
title('real(VLV^{-1})')
colorbar 

subplot(223)
imagesc(imag(aa), [-1 1])
title('imag(V*L*V^{-1})')
colorbar

subplot(222)
imagesc(network_param.J)
title('J')
colorbar

subplot(224)
imagesc(network_param.J - real(aa), [-1 1])
title('J - V*L*V^{-1}')
colorbar

%% Compute (V^{-1}g): 

inp_patt_transformed = V\network_param.input_pattern;


%%

% % sort by the size of the real part
% [~, ord] = sort(real(diag(D)), 'descend');
% 
% D = D(ord, ord);
% V = V(:, ord);
%%
[nStim, nCells, nTimes] = size(proc_simdata.dr_out_peristim);
%%
% also get the eigenvectors of the covariance matrix of activity
dr_out_all = reshape(permute(proc_simdata.dr_out_peristim,[3 1 2]), ... 
    nTimes*nStim, nCells);
%%
[pc_vecs, pc_wts] = pca(dr_out_all, 'NumComponents', 5);
%% Analyze activity projections

pc_wts_trials = cellfun(@(ind) reshape(pc_wts(:, ind), nTimes, nStim), ... 
    num2cell(1:size(pc_wts, 2)), 'UniformOutput', false);

pc_wts_stim1 = cellfun(@(x) x(:, stim1_trials), ... 
    pc_wts_trials, 'UniformOutput', false);
pc_wts_stim2 = cellfun(@(x) x(:, stim2_trials), ... 
    pc_wts_trials, 'UniformOutput', false);

%% project activity onto the associated e-vects of real-valued eigenvalues
is_real = abs(imag(diag(D))) < 1e-16;
real_evects = V(:, is_real);

j_realVects_stim1 = cell(1, 5);
j_realVects_stim2 = cell(1, 5);

for i_ctr = 1:5

    dr_out_vec_proj = reshape(dr_out_all*real_evects(:, i_ctr), nTimes, nStim);
    
    j_realVects_stim1{i_ctr} = dr_out_vec_proj(:, stim1_trials);
    j_realVects_stim2{i_ctr} = dr_out_vec_proj(:, stim2_trials);
    
    
end

j_topEVects_stim1 = cell(1, 10);
j_topEVects_stim2 = cell(1, 10);

i_ctr = 0;
for i_vect = 1:10
i_ctr = i_ctr + 1;
    if D(i_ctr, i_ctr) == conj(D(i_ctr+1, i_ctr+1))
        % take pair of eigenvectors
        dr_v1 = reshape(dr_out_all*V(:, i_ctr), nTimes, nStim);

        j_topEVects_stim1{i_ctr} = real(dr_v1(:, stim1_trials));
        j_topEVects_stim2{i_ctr} = real(dr_v1(:, stim2_trials));
        j_topEVects_stim1{i_ctr+1} = imag(dr_v1(:, stim1_trials));
        j_topEVects_stim2{i_ctr+1} = imag(dr_v1(:, stim2_trials));
        i_ctr = i_ctr + 1;
    else 
        % jsut take one
    
        dr_out_vec_proj = reshape(dr_out_all*V(:, i_ctr), nTimes, nStim);

        j_topEVects_stim1{i_ctr} = dr_out_vec_proj(:, stim1_trials);
        j_topEVects_stim2{i_ctr} = dr_out_vec_proj(:, stim2_trials);
    end
    i_ctr
end

%%
figure()
% hold on
% plot(pc_wts_stim1{1}, pc_wts_stim1{2},'color', input_colors(1, :))
% plot(pc_wts_stim2{1}, pc_wts_stim2{2},'color', input_colors(2, :))

for i_sp = 1:length(pc_wts_stim1)
    subplot(length(pc_wts_stim1), 3, 3*i_sp - 2)
    hold on
    plot(mean(pc_wts_stim1{i_sp}, 2),'color', input_colors(1, :)) 
    plot(mean(pc_wts_stim2{i_sp}, 2),'color', input_colors(2, :)) 
    
    
    
    subplot(length(pc_wts_stim1), 3, 3*i_sp-1)
    hold on
    plot(mean(j_topEVects_stim1{i_sp}, 2),'color', input_colors(1, :)) 
    plot(mean(j_topEVects_stim2{i_sp}, 2),'color', input_colors(2, :)) 
    title(['\lambda = ' num2str(D(i_sp, i_sp))])
    
    
    subplot(length(pc_wts_stim1), 3, 3*i_sp)
    if ~isreal(V(:, i_sp))
        if mod(i_sp, 2) == 1
            plot(sel_info.cell_auroc_val, real(V(:, i_sp)), '.')
            xc_val = corr(sel_info.cell_auroc_val, real(V(:, i_sp)));
        else
            plot(sel_info.cell_auroc_val, imag(V(:, i_sp)), '.')
            xc_val = corr(sel_info.cell_auroc_val, imag(V(:, i_sp)));

        end
        

    else
        plot(sel_info.cell_auroc_val, (V(:, i_sp)), '.')
        xc_val = corr(sel_info.cell_auroc_val, (V(:, i_sp)));

    end
    xlabel('auroc by cell')
    ylabel('eigenvec entry (real/imag for complex \lambda)')
    title(['xc : ' num2str(xc_val, '%1.2f')])

end
    
%% relationship between single-cell selectivity and eigenvectors of the connectivity matrix? 

% compute the linear correlation between eigenvector entries and
% single-cell selectivity for each eigenvector (for complex, take in pairs)
% plot this values against the real part of the eigenvalue
sel_evec_linear_corr = 0*sel_info.cell_auroc_val;
sel_evec_linear_pcorr = 0*sel_info.cell_auroc_val;

for i_v = 1:size(V, 2)
    if ~isreal(D(i_v, i_v))
        if mod(i_v, 2) == 1
            [xc_val, p] = corr(sel_info.cell_auroc_val, real(V(:, i_v)));
            xc_val = abs(xc_val);   % don't know if we took +/- imag part

        else
            [xc_val, p] = corr(sel_info.cell_auroc_val, imag(V(:, i_v)));
            xc_val = abs(xc_val);   % don't know if we took +/- imag part

        end
        

    else
        [xc_val, p] = corr(sel_info.cell_auroc_val, (V(:, i_v)));

    end
    sel_evec_linear_corr(i_v) = xc_val;
    sel_evec_linear_pcorr(i_v) = p;
end
xc_thresh = min(abs(sel_evec_linear_corr(sel_evec_linear_pcorr < 1/500)));
%%
figure(), 

subplot(221)
L_bins = linspace(-8, 8, 100);
c_bins = linspace(-.2, .4, 50);
dens_plot = histcounts2(sel_evec_linear_corr, real(diag(D)), c_bins, L_bins);
dens_plot_smear = imfilter(dens_plot, fspecial('gauss', [3 3], 1));
imagesc(L_bins, c_bins, dens_plot)
hold on
plot(L_bins, xc_thresh + 0*L_bins, 'w')
plot(L_bins, -xc_thresh + 0*L_bins, 'w')
xlabel('Re(\lambda)')
ylabel('linear correlation between selectivity index vector and eigenvector')
% plot(real(diag(L)), sel_evec_linear_corr, '.')
ch = colorbar;
ch.Label.String = 'counts';
set(gca, 'ydir', 'normal');

subplot(222)
histogram(sel_info.cell_auroc_val)

% plot the original eigenvalues colored by how correlated they are with
% selectivity. 
subplot(223)
num_bins = 8;
evec_groups = max(abs(sel_evec_linear_corr))*linspace(0, 1.1, num_bins + 1);
[e_cts, ~, evec_bin_locs] = histcounts(abs(sel_evec_linear_corr), evec_groups);
e_vals = diag(D);
hold on
cmap = flipud(hot(num_bins));
for i_group = 1:(length(evec_groups)-1)
    plot(e_vals(evec_bin_locs == i_group), '.', 'color', cmap(i_group, :), ...
        'markersize', 10);
end
set(gca, 'color', 'none')
axis square
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')
title('color is the linear correlationn btw selectivity and eigenvector pattern')
% plot the input pattern phase factors (this is a number associated with
% each eigenvector, computed as the (V^{-1})*input_pattern

subplot(2,2,4)
num_bins = 8;
inp_pat_var = real(inp_patt_transformed);
inp_groups = max(abs(inp_pat_var))*linspace(0, 1.1, num_bins + 1);
[e_cts, ~, inp_bin_locs] = histcounts(abs(inp_pat_var), inp_groups);
hold on
cmap = flipud(hot(num_bins));
for i_group = 1:(length(inp_groups)-1)
    plot(e_vals(inp_bin_locs == i_group), '.', 'color', cmap(i_group, :), ...
        'markersize', 10);
end
set(gca, 'color', 'none')
axis square
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')
title('color is the real part of the input pattern projected on inv(V))')