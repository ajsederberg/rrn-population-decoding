function plot_handles = assignColorsToLines(plot_handles, color_map_matrix)

if size(color_map_matrix, 1) < length(plot_handles)
    error('color map matrix should be a Nx3 matrix of rgb colors, with N >= length(plot_handles)')
end

for ii = 1:length(plot_handles)
    try
        plot_handles(ii).Color = color_map_matrix(ii, :);
    catch
        set(plot_handles(ii), 'color', color_map_matrix(ii, :));
    end
end