function network_param = generateNetworkParam(network_code, rng_index)
% Sets network parameters according to named conventions. RNG_INDEX is used
% to generate instantiations of each random network. 

% Set the rng seed to the rng_index. 
rng(rng_index);

% encode the fraction I/(E+I) in the decimal of network_code
i_frac = rem(network_code, 1);
ei_val = 1/i_frac - 1;
ei_tag = ['_ei' num2str(round(100*ei_val))];

% encode the network size in the network_code (*1e6)
net_size = floor((network_code + eps)/1e6);
size_tag = ['_n' num2str(net_size)];

% subtract off the net_size*1e6
network_code = floor(network_code - 1e6*net_size);

% encode the network gain in the network code (*1e3)
network_gain = floor(network_code/1e2)/10;
gain_tag = ['_g' num2str(100*network_gain)];

% subtract off the net_size to give the network_index
network_index = floor(network_code - 1000*network_gain);

% Allowed network_tags: 'logNorm', 'splogNorm', 'norm', 'spnorm',
% 'spnormSL', 'spnorm2P', 'norm2PSL', 'spnormB', 'spfix', 'spnormE', ...
% 'spEnoise', 'spnormEL', 'veritas'}
connect_tags = {'logNorm', 'splogNorm', 'norm', 'spnorm', 'spnormSL', ...
    'spnorm2P', 'norm2PSL', 'spnormB', 'spfix', 'spnormE', 'spEnoise', ...
    'spnormEL', 'veritas'};
connect_tag = connect_tags{network_index};

% concatenate. This is compatible with 'parseNetworkTag'
network_tag = [connect_tag gain_tag size_tag ei_tag];
network_param = setDefaultRNNModelParam(network_tag);
network_param.network_tag = [network_tag '_rng' num2str(rng_index)];
network_param.save_dir = ['results/' network_param.network_tag '/'];
network_param.plotdir = [network_param.save_dir 'plots/'];

% matsdir is where the large simulation results will go (on the server, to
% prevent filling up local storage). However, if running this on the
% laptop without a direct network connection, files can be saved in a temp
% folder instead and later transferred. This isn't controlled ... manually
% uncomment the lines below. 
% network_param.matsdir = ['temp/' network_param.network_tag '/matfiles/'];

% network_param.matsdir = ['/Volumes/home/asederberg6/projects/random_network_readouts/'...
%     'results/' network_param.network_tag '/matfiles/'];

% save locally only
network_param.matsdir = ['results/' network_param.network_tag '/matfiles/'];