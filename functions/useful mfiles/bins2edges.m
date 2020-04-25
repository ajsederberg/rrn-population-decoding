function edges = bins2edges(bins)
% gets bin centers from edges vector


edges = zeros(1, length(bins)+1);
edges(2:end-1) = (bins(1:end-1) + bins(2:end))/2;
edges(1) = bins(1) - 0.5*(bins(2)-bins(1));
edges(end) = bins(end) + 0.5*(bins(end)-bins(end-1));