function network_param = initializeNetworkParams(network_code, rng_index)
% Initializes network and saves network parameters. 


network_param = generateNetworkParam(network_code, rng_index);

% this will create new directories if they do not yet exist for this
% network type
if ~exist(network_param.save_dir, 'dir')
    mkdir(network_param.save_dir);
    mkdir([network_param.save_dir 'matfiles/']);
end

if ~exist(network_param.plotdir, 'dir')
    mkdir(network_param.plotdir);
end
if ~exist(network_param.matsdir, 'dir')
    mkdir(network_param.matsdir);
end

save([network_param.matsdir 'network_paramfile'], 'network_param')
save([network_param.save_dir 'matfiles/network_paramfile'], 'network_param')