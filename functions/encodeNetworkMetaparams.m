function network_code = encodeNetworkMetaparams(num_neurons, network_gain, network_index, inh_netwt_frac)
% network_code = encodeNetworkMetaparams(num_neurons, network_gain, network_index, inh_netwt_frac)
%   This function controls the encoding of the entwork parameters such that
%   "generateNetworkParam" will work correctly. 
% num_neurons = 100;        % number of neurons
% inh_netwt_frac = 0.55;    % average fraction of input current from
                                % inhibition
% network_gain = 1.1;       % scales the strength of Jij entries
% network_index = 3;        % toggles between weight distributions
                            % {
network_code = 1e6*num_neurons + 1e3*network_gain + network_index + inh_netwt_frac; 