function [wt_dist_tag, g0, numNeur, ei_ratio] = parseNetworkTag(network_tag)
% % network tags:
% %   (%w)_g(%w)_n(%w)_ei(%w) : first (%w) is the weight distribution tag; second
% %   (%w) is the gainx100; third (%w) is the number of neuron; fourth (%w)
% is the ei ratio x 100. 

parsed_tag_cell = regexp(network_tag, '(\w*)_g(\w*)_n(\w*)_ei(\w*)', 'tokens');
wt_dist_tag = parsed_tag_cell{1}{1};
g0 = str2double(parsed_tag_cell{1}{2})/100;
numNeur = str2double(parsed_tag_cell{1}{3});
ei_ratio = str2double(parsed_tag_cell{1}{4})/100;