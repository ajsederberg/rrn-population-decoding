% Multi-file e-vec density plot, numbers for paper. 

% results_dir = [dropboxDirectory 'Dropbox/postdoctoral_projects/stanley_lab_work/code/github_repositories/theory code/random_rec_network_project/results/'];
% rng_list = [1 2 3 4 5 6 7 8 9 10 11 12];

results_dir = '/Volumes/Samsung_T5/projects/random_rec_network_project/results/';
rng_list = [2 3 4 5 6 7];     % spnormE
% rng_list = [2 3 4]; % spnormEL

x_evec_info = cell(length(rng_list), 1);

for i_rng = 1:length(rng_list)
    
    % onedrive files:
    file_tag = ['spnormE_g200_n500_ei100_rng' num2str(rng_list(i_rng))];
%     file_tag = ['spnormEL_g100_n500_ei100_rng' num2str(rng_list(i_rng))];

    % dropboxDir results files:
%     file_tag = ['spnorm_g400_n500_ei100_rng' num2str(rng_list(i_rng))];

    % experiment specification tag
    expsim_tag = 'f2fdiscExp_pw50';
    
    x_evec_info{i_rng} = load([results_dir file_tag '/matfiles/eigenvectoranalysis_' ...
    file_tag '_' expsim_tag]);

end

%% plot a density of eigenvalues

num_bins = 500;
imag_bins = linspace(-1.5, 1.5, 500);
real_bins = linspace(-1.5, 1.5, 500);
smth_fac = 0.1*num_bins/range(imag_bins);
title_fsz = 7;
label_fsz = 8;

norm_evallims = 1.05*[-1 1 -1 1];
x_incirc = real_bins(abs(real_bins) <= 1);
x_circ = [x_incirc fliplr(x_incirc(1:end-1))];
y_circ = [sqrt(1 - x_incirc.^2) -sqrt(1 - fliplr(x_incirc(1:end-1)).^2)];

circ_mask = bsxfun(@(x, y) x.^2 + y.^2 < 1, edges2bins(real_bins), edges2bins(imag_bins)');
smth_fun = @(mat) imfilter(mat, fspecial('gaussian', ... 
    ceil(6*smth_fac)*[1 1], smth_fac));
bin_size = (range(real_bins)/num_bins)^2;

e_vals = cellfun(@(x) diag(x.linD), x_evec_info, 'UniformOutput', false);
% normalize e-vals 
eval_norm_size = cellfun(@(x) max(abs(x(:))), e_vals);

disp(['mean eigenvalue normalization: ' ...
    num2str(mean(eval_norm_size)) ', +/- ' num2str(std(eval_norm_size))])

e_vals = cellfun(@(x) x/max(abs(x(:))), e_vals, 'UniformOutput', false);

cts_eval = cellfun(@(x) histcounts2(imag(x), real(x), imag_bins, real_bins), ... 
    e_vals, 'UniformOutput', false);

% density plot of all eigenvalues
mean_cts_eval = squeeze(mean(cell2mat(cellfun(@(x) shiftdim(x, -1), cts_eval, 'UniformOutput', false))));
% smooth_eval_pmf = imfilter(mean_cts_eval, fspecial('gaussian', ... 
%     ceil(4*smth_fac)*[1 1], smth_fac));
smooth_eval_pmf = smth_fun(mean_cts_eval);
sum_smooth_pmf = sum(smooth_eval_pmf(:));
smooth_eval_pdf = smooth_eval_pmf/sum_smooth_pmf/(bin_size);

% pull out the eigenvalues with eigenvectors correlated singificantly with
% selectivity
p_cutoff = 0.05;
evecs_w_xcSel = cellfun(@(x) x.sel_evec_linear_pcorr < p_cutoff/500, x_evec_info, ...
    'UniformOutput', false);
% density plot of correlated eigenvalues
cts_xceval = cellfun(@(x, inds) histcounts2(imag(x(inds)), real(x(inds)), imag_bins, real_bins), ... 
    e_vals, evecs_w_xcSel, 'UniformOutput', false);
mean_cts_xceval = squeeze(mean(cell2mat(cellfun(@(x) shiftdim(x, -1), cts_xceval, 'UniformOutput', false))));
% smooth_xceval_pmf = imfilter(mean_cts_xceval, fspecial('gaussian', ... 
%     ceil(4*smth_fac)*[1 1], smth_fac));
smooth_xceval_pmf = smth_fun(mean_cts_xceval);
smooth_xceval_pmf = 0.5*[smooth_xceval_pmf + flipud(smooth_xceval_pmf)];
smooth_xceval_pdf = smooth_xceval_pmf/sum_smooth_pmf/bin_size;

makeMyFigure(7.5*2.54, 3*2.54);
subplot(131)
imagesc(real_bins, imag_bins, log10(circ_mask.*smooth_eval_pdf), [-4 0])
hold on
plot(0*real_bins, imag_bins, 'k', 'linewidth', 1.5)
plot(real_bins, 0*imag_bins, 'k', 'linewidth', 1.5)
% plot(x_circ, y_circ, 'k', 'linewidth', 1)
axis tight
axis square
xlabel('Re(\lambda)/a')
ylabel('Im(\lambda)/a')
ch = colorbar;
ch.Label.String = 'log_{10} density';
axis(norm_evallims)
title(['all \lambda, N = ' num2str(length(x_evec_info)) ' networks'], 'FontSize', title_fsz)
set(gca, 'ydir', 'normal', 'xtick', [-1 0 1], 'xticklabel', {'-1','0', '1'}, ...
    'ytick', [-1 0 1], 'yticklabel', {'-1','0', '1'}, 'fontsize', label_fsz)
text(-0.25, 1.2, 'A', 'Units', 'normalized', 'FontSize', 14)

subplot(132)
imagesc(real_bins, imag_bins, log10(circ_mask.*smooth_xceval_pdf), [-4 0])
hold on
plot(0*real_bins, imag_bins, 'k', 'linewidth', 1.5)
plot(real_bins, 0*imag_bins, 'k', 'linewidth', 1.5)
% plot(x_circ, y_circ, 'k', 'linewidth', 1)

axis tight
axis square
xlabel('Re(\lambda)/a')
ylabel('Im(\lambda)/a')
ch = colorbar;
ch.Label.String = 'log_{10} density';
axis(norm_evallims)

title('selectivity-correlated \lambda', 'fontsize', title_fsz)
set(gca, 'ydir', 'normal', 'xtick', [-1 0 1], 'xticklabel', {'-1','0', '1'}, ...
    'ytick', [-1 0 1], 'yticklabel', {'-1','0', '1'}, 'fontsize', label_fsz)
text(-0.25, 1.2, 'B', 'Units', 'normalized', 'FontSize', 14)

subplot(133)
frac_xceval_of_eval = smooth_xceval_pmf./(eps + smooth_eval_pmf);

imagesc(real_bins, imag_bins, circ_mask.*frac_xceval_of_eval)
hold on
plot(0*real_bins, imag_bins, 'k', 'linewidth', 1.5)
plot(real_bins, 0*imag_bins, 'k', 'linewidth', 1.5)
axis tight
axis square
xlabel('Re(\lambda)/a')
ylabel('Im(\lambda)/a')
ch = colorbar;
ch.Label.String = 'fraction (by bin)';
axis(norm_evallims)

% title('')
set(gca, 'ydir', 'normal', 'xtick', [-1 0 1], 'xticklabel', {'-1','0', '1'}, ...
    'ytick', [-1 0 1], 'yticklabel', {'-1','0', '1'}, 'fontsize', label_fsz)
text(-0.25, 1.2, 'C', 'Units', 'normalized', 'FontSize', 14)


% print(gcf, '-dpdf', ['plots/eigenvalues_' file_tag])
% print(gcf, '-dsvg', ['plots/eigenvalues_' file_tag])
print(gcf, '-dpdf', 'writeup/figures/fig5_final')
print(gcf, '-dsvg', 'writeup/figures/fig5_final')

%%  are there really eigenvectors that are correlated with selectivity but 
% have negative real parts of the eigenvalue? yes. 
figure(),
eval_bins = linspace(-1, 1, 21);
eval_cts = zeros(length(e_vals), length(eval_bins)-1);

corr_vals = cell(length(e_vals), 1);
eval_vals = cell(length(e_vals), 1);

for i_sim = 1:length(e_vals)

    eval_cts(i_sim, :) = histcounts(...
        real(e_vals{i_sim}(...
        x_evec_info{i_sim}.sel_evec_linear_pcorr < .0001)), ...
        eval_bins);
    
    eval_vals{i_sim} = (e_vals{i_sim}(...
        x_evec_info{i_sim}.sel_evec_linear_pcorr < .0001));
    corr_vals{i_sim} = real(x_evec_info{i_sim}.sel_evec_linear_corr(...
        x_evec_info{i_sim}.sel_evec_linear_pcorr < .0001));

end
subplot(121)
bar(eval_bins(2:end), eval_cts, 'stacked')

num_real_evals = zeros(size(eval_vals));
num_imag_evals = zeros(size(eval_vals));

subplot(122)
hold on
for i_sim = 1:length(eval_vals)
    
    real_eigenvalues = abs(imag(eval_vals{i_sim})) < 1e-16;
    
    num_real_evals(i_sim) = sum(real_eigenvalues);
    num_imag_evals(i_sim) = length(unique(real(eval_vals{i_sim}(~real_eigenvalues))));
    sum(real_eigenvalues)
    plot(real(eval_vals{i_sim}(real_eigenvalues)), corr_vals{i_sim}(real_eigenvalues), 'o')
end
ylabel('linear correlation btw eigenvector and selectivity')
xlabel('Re [\lambda] for corresponding eigenvalue')

mean_num_pairs_imag = mean(num_imag_evals)
mean_num_real = mean(num_real_evals)

num_imag_evals
num_real_evals

max_pos_corr = cellfun(@(x) max(x( x > 0)), corr_vals)
mean_max_pos_corr = mean(max_pos_corr)
std_max_pos_corr = std(max_pos_corr)

num_neg_corr = cellfun(@(x) sum( x < 0), corr_vals)

ave_pos_corr = cellfun(@(x) mean(x( x > 0)), corr_vals)
mean_ave_pos_corr = mean(ave_pos_corr)
std_mean_ave_pos_corr = std(ave_pos_corr)


num_neg_corr = cellfun(@(x) sum(x < 0 ), corr_vals)
num_neg_corr_real_eval = cellfun(@(c, l) sum(c < 0 & abs(imag(l)) < 1e-16), ... 
    corr_vals, eval_vals)
num_pos_corr_real_eval = cellfun(@(c, l) sum(c > 0 & abs(imag(l)) < 1e-16), ... 
    corr_vals, eval_vals)
num_neg_corr_imag_eval = cellfun(@(c, l) sum(c < 0 & abs(imag(l)) > 1e-16), ... 
    corr_vals, eval_vals)
%%
sym_xceval = cellfun(@(x) x(1:249, :) | flipud(x(251:end, :)), cts_xceval, 'UniformOutput', false)
num_xceval = cellfun(@(x) sum(x(:)), sym_xceval)
num_real_xceval = cellfun(@(x) sum(x(250, :)), cts_xceval)
num_imag_pairs = (num_xceval - num_real_xceval)/2

%%
b = 1.2;
drdx = @(x) 0.5./(cosh(x - b).^2);
a_b = proc_simdata.baseline_activity';

J = network_param.J;

x = J*a_b;

r_x = 0.5*(1+tanh(x - b));

lin_J = J*diag(drdx(x));

figure(), 
plot(eig(lin_J), '.')
hold on
% plot(eig(J), '.')


%% questions:
% do the linearized eigenvectors predict selectivity? can you find
% them/track them? 
% about the fixed point, do the linearized eigenvectors predict
% selectivity? what is the period? if range is omega ~ 0.1, that is a time
% window of T = 60 tau, which is very long. What if we aimed for T = 10
% tau? 
%%
i_trial = 3;
figure();
%%


sel_cells = [ 128 144 220 244 260 429 485];
for i_step = 1:281
    r1 = a_b' + squeeze(proc_simdata.dr_out_peristim(i_trial, :, i_step));
 
    
    subplot(2, 1, 2)
    hold on
    plot(i_step, r1(sel_cells), '.')

    
    subplot(2, 1, 1)
    cla;
    plot(eig(J*diag(drdx(r1))), '.')
    xlim([-1.1 1.1])
    ylim([-1.1 1.1])
    
    pause
    
end

%%

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



%% Scaling with n plot
p = 0.2;
pL = 0.1;
p_connect_vals = logspace(-3, -1, 7);
n = logspace(2, 5);
% this assumes ~normal distributions (pn is large) with binomial
% distribution sigma values 
num_cells_above_cutoff = @(pp,nn,frac) nn.*erfc(frac*pp*nn./sqrt(2*pp*(1-pp)*nn))/2;


% compute how many inhibitory cells would have an excess of 15% of connections 
frac_cutoff = 0.15;
% y = n.*erfc(frac_cutoff*p*n./sqrt(2*p*(1-p)*n))/2

figure(),
hold on
for i_pcv = 1:length(p_connect_vals)
    plot((n), num_cells_above_cutoff(p_connect_vals(i_pcv), n, frac_cutoff), 'linewidth', 1)
end
y_lims = ylim;
eqline
ylim(y_lims)
% plot(log10(n), num_cells_above_cutoff(pL, n, frac_cutoff), 'linewidth', 1)
ylabel('# of potentially selective i-cells')
xlabel('# of i-cells')
set(gca, 'fontsize', 14)
legend(num2str(p_connect_vals', 'p = %1.3f'), 'Location', 'northwest')