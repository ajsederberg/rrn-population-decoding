function input_params = lookUpInputParam(network_code, rng_index, sim_name)
% Looks up the network_params from the saved file
% This controls the network better for multiple runs. 

network_param = lookUpNetworkParam(network_code, rng_index);

try load([network_param.matsdir 'sim_data_' sim_name], 'input_params');
    

catch
        %%
        % this looks up the network_param from previously run files, where
        % 'network_param' was not automatically saved. 
        server_dir = ['/Volumes/home/asederberg6/projects/random_network_readouts/'...
            'results/' network_param.network_tag '/matfiles/'];
        x_res = load([server_dir 'sim_data_' sim_name], 'input_params');

        input_params = x_res.input_params;
        
%%
end