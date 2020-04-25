function network_param = lookUpNetworkParam(network_code, rng_index)
% Looks up the network_params from the saved file
% This controls the network better for multiple runs. 

network_param = generateNetworkParam(network_code, rng_index);

try load([network_param.save_dir 'matfiles/network_paramfile'], 'network_param');
    
catch
    try
        % this looks up the network_param from previously run files, where
        % 'network_param' was not automatically saved. 
        file_list = dir(network_param.matsdir);
        is_sim_data = arrayfun(@(x) contains(x.name, 'sim'), file_list);
        sim_file = file_list(is_sim_data);
        x_res = load([network_param.matsdir sim_file(1).name], 'network_param');

        network_param = x_res.network_param;

        save([network_param.save_dir 'matfiles/network_paramfile'], 'network_param');
    catch
        %%
        % this looks up the network_param from previously run files, where
        % 'network_param' was not automatically saved. 
        server_dir = ['/Volumes/home/asederberg6/projects/random_network_readouts/'...
            'results/' network_param.network_tag '/matfiles/'];
        file_list = dir(server_dir);
        is_sim_data = arrayfun(@(x) contains(x.name, 'sim'), file_list);
        sim_file = file_list(is_sim_data);
        x_res = load([server_dir sim_file(1).name], 'network_param');

        network_param = x_res.network_param;
        
        network_param.old_matsdir = network_param.matsdir;
        network_param.matsdir = [network_param.save_dir 'matfiles/'];
%%
        save([network_param.save_dir 'matfiles/network_paramfile'], 'network_param');        
    end
end