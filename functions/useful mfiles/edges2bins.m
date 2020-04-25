function bins = edges2bins(edges)
% gets bin centers from edges vector

bins = (edges(1:end-1) + edges(2:end))/2;