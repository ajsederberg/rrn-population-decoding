function param = setDefaultRNNModelParam(numNeur, g0, fracE, gIEratio, varargin)

% Default values
% network equations: fr network, with dr/dt = -r + Jx + inp and x = f(r). 
param.numNeur = numNeur;

% have to make some neurons E and some I [argument is that having info in
% the I population implies selective inhibition or selective pools of
% inhibition]. Initially made 2 E populations ; now not sure why... but
% could be put back in 
% fracE = 0.8;
param.numNeurExc = round(fracE*numNeur/2);
param.numNeurInh = numNeur - 2*param.numNeurExc;
numEI_ratio = 2*param.numNeurExc/param.numNeurInh;
% define connectivity matrix. Could be a lot of things..... 
% g0 = 1;
gEE = g0;
gEI = (gIEratio*numEI_ratio)*g0;
gIE = g0;        % in the log-norm weight network, setting gIE = gEI pushes the outline e-val
                    % onto the real axis. The imaginary part grows with
                    % gIE/gEI. 
gII = (gIEratio*numEI_ratio)*g0;         % large gII in log-norm weights: large negative real e-vals
% log-normal? hmm. 
% J = 1/sqrt(numNeur)*[ gEE*abs(randn(numNeurExc)) gIE*abs(randn(numNeurExc, numNeurInh)); ...
%     -gEI*abs(randn(numNeurInh, numNeurExc)) -gII*abs(randn(numNeurInh))];

numNeurExc = param.numNeurExc;
numNeurInh = param.numNeurInh;
param.J = 1/sqrt(numNeur)*[ gEE*exp(randn(numNeurExc)) gEE*exp(randn(numNeurExc)) -gEI*exp(randn(numNeurExc, numNeurInh)); ...
    gEE*exp(randn(numNeurExc)) gEE*exp(randn(numNeurExc)) -gEI*exp(randn(numNeurExc, numNeurInh)); ...
    gIE*exp(randn(numNeurInh, numNeurExc)) gIE*exp(randn(numNeurInh, numNeurExc)) -gII*exp(randn(numNeurInh))];


% define nonlinear i-o function . 
param.fun_io = @(x) (1 + tanh(3*x))/2;


% Overwrite if necessary
num_pairs = (nargin-4)/2;
for i_va = 1:num_pairs
    param.(varargin{2*i_va - 1}) = varargin{2*i_va};
end