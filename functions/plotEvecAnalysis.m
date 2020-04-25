function plotEvecAnalysis(file_tag)

load(['results/' file_tag '/matfiles/eigenvectoranalysis_' ...
    file_tag '_f2fdiscExp_pw50'])

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
histogram(cell_auroc_val)

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
