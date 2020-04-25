function param = setDefaultRNNModelParam(network_tag, varargin)
% If this function is behaving unexpectedly with old code, try the legacy
% version: 'setDefaultRNNModelParam_legacy'
% This works with 'generateNetworkParam," which is only called once (when
% the amplitude-normalization is calculated). The 'lookUpNetworkParam'
% function uses this to find the correct directories, but then recovers
% parameters from a saved file. 
% % network tags:
% %   (%w)_g(%w)_n(%w)_ei(%w) : first (%w) is the weight distribution tag; second
% %   (%w) is the gainx100; third (%w) is the number of neuron; fourth (%w)
% is the ei ratio x 100. 
% Default values
fracE = 0.8;

[wt_dist_tag, g0, numNeur, ei_ratio] = parseNetworkTag(network_tag);

% network equations: fr network, with dr/dt = -r + Jx + inp and x = f(r). 
param.numNeur = numNeur;
param.numNeurExc = round(fracE*numNeur);
param.numNeurInh = numNeur - param.numNeurExc;
param.network_gain = g0;

% numerical EI ratio is the # exc neuron / # inh neurons
numEI_ratio = param.numNeurExc/param.numNeurInh;

% set the i-o function. [alternative here: supra-linear? saturation
% nonlinearity ensures stability for now. ]
% This includes a constant bias, which sets the 0-input firing rate of the
% network to 0.05. Without the bias, the 0-input firing rate is 0.5, and
% activity gets saturated. (Change made: 4/19/2018)
% Change 5/4/18: bias from 1.5 to 2. 
% Change 5/8/18: allow this to be changed in the 'wt_dist_tag'. This is the
% default function, but can use a supralinear function instead. 

% this lowers the bias in the i-o function. 
if strcmp(wt_dist_tag, 'spnormEL')
    b = 1.2;
else
    b = 2;
end
param.b = b;
param.fun_io = @(x) (1 + tanh(x - b))/2;
param.iofun_name = 'tanh';

% Give an input pattern to *20%* of excitatory neurons. (! used to be 50%!)
% Because connectivity is random, these are just the first 20% of
% excitatory cells. 
if ~(strcmpi(wt_dist_tag, 'spnormE') || strcmpi(wt_dist_tag, 'spnormEL') || strcmpi(wt_dist_tag, 'spEnoise') || strcmpi(wt_dist_tag, 'veritas'))
    param.input_pattern = false(numNeur, 1);
    num_inputs = 0.2*param.numNeurExc;
    % param.input_pattern(1:param.numNeurExc) = rand(param.numNeurExc, 1) < 0.2;
    param.input_pattern(1:num_inputs) = true;
else
    % for 'spEnoise' and 'spnormE', 'spnormEL', all E-cells get inputs. 
    param.input_pattern = false(numNeur, 1);
    num_inputs = param.numNeurExc;
    % param.input_pattern(1:param.numNeurExc) = rand(param.numNeurExc, 1) < 0.2;
    param.input_pattern(1:num_inputs) = true;        
end
if strcmpi(wt_dist_tag, 'logNorm')
    % define connectivity matrix, with log-normal weight distribution
    % g0 = 1;
    gEE = g0;
    gEI = (numEI_ratio/ei_ratio)*g0;
    gIE = g0;        % in the log-norm weight network, setting gIE = gEI pushes the outline e-val
                        % onto the real axis. The imaginary part grows with
                        % gIE/gEI. 
    gII = (numEI_ratio/ei_ratio)*g0;         % large gII in log-norm weights: large negative real e-vals
    % log-normal
    numNeurExc = param.numNeurExc;
    numNeurInh = param.numNeurInh;
    param.J = 1/sqrt(numNeur)*[ gEE*exp(randn(numNeurExc)) -gEI*exp(randn(numNeurExc, numNeurInh)); ...
         gIE*exp(randn(numNeurInh, numNeurExc)) -gII*exp(randn(numNeurInh))];
    % define nonlinear i-o function . 
elseif strcmpi(wt_dist_tag, 'splogNorm')
    % define connectivity matrix, with log-normal weight distribution
    % g0 = 1;
    gEE = g0;
    gEI = (numEI_ratio/ei_ratio)*g0;
    gIE = g0;        % in the log-norm weight network, setting gIE = gEI pushes the outline e-val
                        % onto the real axis. The imaginary part grows with
                        % gIE/gEI. 
    gII = (numEI_ratio/ei_ratio)*g0;         % large gII in log-norm weights: large negative real e-vals
    % log-normal
    numNeurExc = param.numNeurExc;
    numNeurInh = param.numNeurInh;
    param.J = 1/sqrt(numNeur)*[ gEE*exp(randn(numNeurExc)) -gEI*exp(randn(numNeurExc, numNeurInh)); ...
         gIE*exp(randn(numNeurInh, numNeurExc)) -gII*exp(randn(numNeurInh))];
     
    % 20% sparsity
    param.sparsity = 0.2;
    sparse_mask = double(rand(param.numNeur) < 0.2);
    param.J = param.J.*sparse_mask/param.sparsity;
elseif strcmpi(wt_dist_tag, 'norm')
    % define connectivity matrix, with normal weight distribution
    % g0 = 1;
    param.sparsity = 1; % 100% of connections made
%     frac_inputs = sum(param.input_pattern)/param.numNeurExc;

    gEE = g0;
    gEI = (numEI_ratio/ei_ratio)*g0; % sets 'balance' condition if ei_ratio = 1
    gIE = g0;        
    gII = (numEI_ratio/ei_ratio)*g0; % (1-frac_inputs) accounts for including the 
                                                     % cells with direct
                                                     % inputs in the
                                                     % excitatory neuron
                                                     % count. 
    
    numNeurExc = param.numNeurExc;
    numNeurInh = param.numNeurInh;
    % rectify weights so that e is positive and i is negative. All-to-all
    % connections. 
    param.J = 1/sqrt(numNeur)*[ abs(gEE + randn(numNeurExc)) -abs(gEI + randn(numNeurExc, numNeurInh)); ...
         abs(gIE + randn(numNeurInh, numNeurExc)) -abs(gII + randn(numNeurInh))];
elseif  strcmpi(wt_dist_tag, 'spnormE') || strcmpi(wt_dist_tag, 'spnormEL') || strcmpi(wt_dist_tag, 'spEnoise')
    %%
    numNeurExc = param.numNeurExc;
    numNeurInh = param.numNeurInh;

% 20% sparsity
    p_connect = 0.2;
    
    % define connectivity matrix : all excitatory neurons receive inputs
    % g0 = 1;
    gEE = g0;
    % set gEI to approximately balance: (Jee-Jei)*r = -1 [the -1 comes from
    % "input" size, measured in the tanh scale (r = tanh(x - b))
    c0 = 4; % was set to 1;
    gEI = (c0 + p_connect*numNeurExc*gEE)/(p_connect*numNeurInh);

    % set gII and gEI to balance approximately
    gII = (numNeurExc/numNeurInh)*gEE; % was *.4;
    gIE = (numNeurInh/numNeurExc)*gII;
    
%     %%  latham par
%     gEE = 0.25;
%     gEI = 0.87;
%     gIE = 0.87;
%     gII = 2;
%     
    rtvarI = 2;
    rtvarE = 2;

%     gEI = iToe_factor*g0;
%     gIE = g0;        % in the norm weight network, setting gIE = gEI pushes the outline e-val
%                         % onto the real axis. The imaginary part grows with
%                         % gIE/gEI. 
%     gII = iToe_factor*g0;         % large gII in log-norm weights: large negative real e-vals
    % log-normal

param.J = 1/sqrt(p_connect*numNeur)*[ rect(gEE + rtvarE*(randn(numNeurExc))) -rect(gEI + rtvarI*(randn(numNeurExc, numNeurInh))); ...
     rect(gIE + rtvarE*(randn(numNeurInh, numNeurExc))) -rect(gII + rtvarI*(randn(numNeurInh)))];

% sparseness enforcement
param.sparsity = p_connect;
sparse_mask = double(rand(param.numNeur) < p_connect);
param.J = param.J.*sparse_mask; %/param.sparsity;

param.gEE = gEE;
param.gEI = gEI;
param.gIE = gIE;
param.gII = gII;

    
elseif strcmpi(wt_dist_tag, 'spnorm') || strcmp(wt_dist_tag, 'spnormSL') || strcmp(wt_dist_tag, 'spnormB') || strcmpi(wt_dist_tag, 'spfix')
    % define connectivity matrix, with normal weight distribution
    % g0 = 1;
    gEE = g0;
    
    if strcmp(wt_dist_tag, 'spnormB')
        frac_inputs = sum(param.input_pattern)/param.numNeurExc;
        iToe_factor = (numEI_ratio/ei_ratio)*(1-frac_inputs);
        
        rtvarE = 1/iToe_factor;
        rtvarI = 1;
    elseif strcmp(wt_dist_tag, 'spfix')
        % no spread in weights
        iToe_factor = (numEI_ratio/ei_ratio);
        rtvarI = 0;
        rtvarE = 0;        
    elseif strcmpi(wt_dist_tag, 'spnormE') || strcmpi(wt_dist_tag, 'spnormEL') 
        iToe_factor = (numEI_ratio/ei_ratio);   % this needs to be 
                                                % multiplied by an extra 
                                                % factor to balance thhe 
                                                % E/I cells when only E 
                                                % cells get direct inputs
        rtvarI = 1;
        rtvarE = 1;
    else
        iToe_factor = (numEI_ratio/ei_ratio);
        rtvarI = 1;
        rtvarE = 1;
    end
    gEI = iToe_factor*g0;
    gIE = g0;        % in the norm weight network, setting gIE = gEI pushes the outline e-val
                        % onto the real axis. The imaginary part grows with
                        % gIE/gEI. 
    gII = iToe_factor*g0;         % large gII in log-norm weights: large negative real e-vals
    % log-normal
    numNeurExc = param.numNeurExc;
    numNeurInh = param.numNeurInh;
%     param.J = 1/sqrt(numNeur)*[ gEE*abs(randn(numNeurExc)) -gEI*abs(randn(numNeurExc, numNeurInh)); ...
%          gIE*abs(randn(numNeurInh, numNeurExc)) -gII*abs(randn(numNeurInh))];
%      
    param.J = 1/sqrt(numNeur)*[ abs(gEE + rtvarE*(randn(numNeurExc))) -(gEI + rtvarI*(randn(numNeurExc, numNeurInh))); ...
         abs(gIE + rtvarE*(randn(numNeurInh, numNeurExc))) -1*(gII + rtvarI*(randn(numNeurInh)))];
     %%
    % 20% sparsity
    param.sparsity = 0.2;
    sparse_mask = double(rand(param.numNeur) < 0.2);
    param.J = param.J.*sparse_mask; %/param.sparsity;
    
    param.gEE = gEE;
    param.gEI = gEI;
    param.gIE = gIE;
    param.gII = gII;
    
    if strcmp(wt_dist_tag, 'spnormSL')
        param.fun_io = @(x) 0.04*rect(x).^2; % SSN with k = 0.01 and n = 2
        param.iofun_name = 'SSN';
    elseif strcmp(wt_dist_tag, 'spnormB')
        % remove recurrent connections back onto the inputs
        param.J(param.input_pattern, :) = 0;

    end
    
elseif strcmp(wt_dist_tag, 'spnorm2P') || strcmp(wt_dist_tag, 'spnorm2PSL')
    gEE = g0;
    gEI = (numEI_ratio/ei_ratio)*g0;
    gIE = g0;        % in the norm weight network, setting gIE = gEI pushes the outline e-val
                        % onto the real axis. The imaginary part grows with
                        % gIE/gEI. 
    gII = (numEI_ratio/ei_ratio)*g0;         % large gII in log-norm weights: large negative real e-vals
    % log-normal
    numNeurExc = param.numNeurExc;
    numNeurInh = param.numNeurInh;
%     param.J = 1/sqrt(numNeur)*[ gEE*abs(randn(numNeurExc)) -gEI*abs(randn(numNeurExc, numNeurInh)); ...
%          gIE*abs(randn(numNeurInh, numNeurExc)) -gII*abs(randn(numNeurInh))];
%      
    param.J = 1/sqrt(numNeur)*[ abs(gEE + (randn(numNeurExc))) -(gEI + (randn(numNeurExc, numNeurInh))); ...
         abs(gIE + (randn(numNeurInh, numNeurExc))) -(gII + (randn(numNeurInh)))];
     %%
    % force E1 -> I1 -> E2 and E2 -> I2 -> E1 model. inputs? 
    nE1 = floor(numNeurExc/2);
    nE2 = numNeurExc - nE1;
    nI1 = floor(numNeurInh/2);
    nI2 = numNeurInh - nI1;
    
    pool_mask = [ ones(nE1, nE1) zeros(nE1, nE2)  ones(nE1, nI1) zeros(nE1, nI2); ... 
                 zeros(nE2, nE1)  ones(nE2, nE2) zeros(nE2, nI1)  ones(nE2, nI2); ...
                 zeros(nI1, nE1)  ones(nI1, nE2) zeros(nI1, nI1) zeros(nI1, nI2); ... 
                  ones(nI2, nE1) zeros(nI2, nE2) zeros(nI2, nI1) zeros(nI2, nI2)];
    % 40% sparsity
    param.sparsity = 0.4;
    sparse_mask = double(rand(param.numNeur) < param.sparsity );    
              
              
    param.J = param.J.*pool_mask.*sparse_mask; %/param.sparsity;
    
    param.gEE = gEE;
    param.gEI = gEI;
    param.gIE = gIE;
    param.gII = gII;   
    
    if strcmp(wt_dist_tag(end-1:end), 'SL')
        param.fun_io = @(x) 0.04*rect(x).^2; % SSN with k = 0.01 and n = 2
        param.iofun_name = 'SSN';
    end
    

    % finally, balance the input pattern so that one pool does not get more
    % inputs than the other
    % Give an input pattern to *20%* of excitatory neurons. (! used to be 50%!)
    param.input_pattern = false(numNeur, 1);
    num_inputs_per_pool = round(nE1*0.2);
    param.input_patterns(1:num_inputs_per_pool) = true;
    param.input_patterns(nE1 + (1:num_inputs_per_pool)) = true;
    
%     param.input_pattern(randperm(nE1, num_inputs_per_pool)) = 1;
%     param.input_pattern(nE1 + randperm(nE2, num_inputs_per_pool)) = 1;

elseif strcmpi(wt_dist_tag, 'veritas')
    % this is a set of parameters chosen to closely resemble probabbilities
    % of connection within S1 and V1 : mainly, pEE is low (0.1) while pEI
    % and pIE are high (80%) and the synaptic weights are drawn from a
    % heavy-tailed distribution
    % check Avermann et al for typical mean/median/ranges for PSPs between
    % cell types: idea here is to have high probability of connection, but
    % most connections are weak, and it is the heavy tails that generate
    % selectivity [ are these distribbutions consistent - if we only had
    % ~50 observations, would we get similar relationships between
    % mean/median and range? 
    %%
    gEE = g0;
    gEI = (numEI_ratio/ei_ratio)*g0;
    gIE = g0;        % in the log-norm weight network, setting gIE = gEI pushes the outline e-val
                        % onto the real axis. The imaginary part grows with
                        % gIE/gEI. 
    gII = (numEI_ratio/ei_ratio)*g0;         % large gII in log-norm weights: large negative real e-vals
    % log-normal
    syn_sig_val = 1;
    syn_mu_val = -0.5;
    syn_randn = @(nR,nC) syn_mu_val + syn_sig_val*randn(nR, nC);
    numNeurExc = param.numNeurExc;
    numNeurInh = param.numNeurInh;
%     param.J = 1/sqrt(numNeur)*[ gEE*exp(syn_sig_val*randn(numNeurExc)) -gEI*exp(syn_sig_val*randn(numNeurExc, numNeurInh)); ...
%          gIE*exp(syn_sig_val*randn(numNeurInh, numNeurExc)) -gII*exp(syn_sig_val*randn(numNeurInh))];
%      
    param.J = 1/sqrt(numNeur)*[gEE*exp(syn_randn(numNeurExc, numNeurExc)) -gEI*exp(syn_randn(numNeurExc, numNeurInh)); ...
         gIE*exp(syn_randn(numNeurInh, numNeurExc)) -gII*exp(syn_randn(numNeurInh, numNeurInh))];
    param.gEE = gEE;
%     param.gEI = gEI;
    param.gIE = gIE;
%     param.gII = gII; 
 %%
    % set connection probabilities, by Avermann
    pEE = 0.17;
    pEI = 0.6;
    pIE = 0.57;
    pII = 0.55;
    % save connection probbabilities
    param.pEE = pEE;
    param.pEI = pEI;
    param.pIE = pIE;
    param.pII = pII; 
    
    p_mat = [pEE*ones(numNeurExc) pEI*ones(numNeurExc, numNeurInh); ... 
        pIE*ones(numNeurInh, numNeurExc) pII*ones(numNeurInh)];
    sparse_mask = rand(numNeur) < p_mat;
    
    param.J = param.J.*sparse_mask;
    
    ei_rescale_factor = ei_ratio*sum(mean(param.J(:, 1:numNeurExc)))/...
        sum(mean(param.J(:, numNeurExc+1:end)));
    param.J(:, numNeurExc+1:end) = abs(ei_rescale_factor)*param.J(:, numNeurExc+1:end);

%     
end

% Overwrite if necessary
num_pairs = (nargin-1)/2;
for i_va = 1:num_pairs
    param.(varargin{2*i_va - 1}) = varargin{2*i_va};
end