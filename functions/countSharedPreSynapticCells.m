function [shared_pre_J, pre_J_corr] = countSharedPreSynapticCells(J)
% Count the number of shared pre-synaptic cells and the correlation between
% the pre-synaptic weights for all cell pairs. 

shared_pre_J = nan(size(J));
pre_J_corr = nan(size(J));

for ii = 2:size(J, 1) 
    for jj = 1:ii-1
        pre_ii = J(ii, :);
        pre_jj = J(jj, :);
        shared_pre_J(ii, jj) = sum(pre_ii ~= 0 & pre_jj ~= 0);
        if sum(pre_ii ~= 0 & pre_jj ~= 0) ~= 0
            pre_J_corr(ii, jj) = offDiag(corrcoef([pre_ii' pre_jj']));
        end
    end
end

shared_pre_J = shared_pre_J + shared_pre_J';
pre_J_corr = pre_J_corr + pre_J_corr';
