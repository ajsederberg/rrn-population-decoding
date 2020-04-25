function [varCE, cell_aveVarCE, cell_stimAveVarCE] = computeVarCE_rateModel(stim_resps, stim_y)

stim_vals = num2cell(unique(stim_y));
varCE = cellfun(@(x) var(stim_resps(stim_y == x, :, :), [], 1), stim_vals, ...
    'UniformOutput', false);

ave_varCE = squeeze(mean(cell2mat(varCE), 1));

cell_aveVarCE = mean(ave_varCE, 1);
cell_stimAveVarCE = squeeze(mean(cell2mat(varCE), 2));