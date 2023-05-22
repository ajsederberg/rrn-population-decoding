% random matrix eigenvectors and params
param.numNeurExc = 400;
param.numNeurInh = 100;
param.numNeur = param.numNeurExc + param.numNeurInh;
g0 = 3;
%%

numNeurExc = param.numNeurExc;
numNeurInh = param.numNeurInh;
numNeur = param.numNeur;

    % 20% sparsity
    p_connect = 0.2;
    
    % define connectivity matrix : all excitatory neurons receive inputs
    % g0 = 1;
    gEE = g0;
    % set gEI to approximately balance: (Jee-Jei)*r = -1 [the -1 comes from
    % "input" size, measured in the tanh scale (r = tanh(x - b))
    c0 = 4;
    gEI = (c0 + p_connect*numNeurExc*gEE)/(p_connect*numNeurInh);

    % set gII and gEI to balance approximately
    gII = (numNeurExc/numNeurInh)*gEE;
    gIE = (numNeurInh/numNeurExc)*gII;
    
%     %%  latham par
%     gEE = 0.25;
%     gEI = 0.87;
%     gIE = 0.87;
%     gII = 2;
%     
    rtvarI = 0;
    rtvarE = 0;

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
%%

figure()
[v,l] = eig(param.J);
[~, ord] = sort(real(diag(l)), 'descend');
subplot(211)
plot(diag(l), 'o')
title(['c0: ' num2str(c0) ', gEE: ' num2str(gEE) ', gEI: ' num2str(gEI) ... 
    ', gIE: ' num2str(gIE) ', gII: ' num2str(gII)])
axis equal
subplot(212)
plot(v(:, ord(1)))
hold on
vE = mean(v(1:numNeurExc, ord(1)));
vI = mean(v(1+numNeurExc:end, ord(1)));
plot([1 numNeurExc], vE*[1 1])
plot([1+numNeurExc numNeur], vI*[1 1])
title(['first eigenvector, vE = ' num2str(vE, '%1.2f') ...
    ', vI = ' num2str(vI, '%1.2f') ', ratio: ' num2str(vI/vE, '%1.2f')])

%% find input-responses
input_levels = linspace(0.5, 3, 11);

r_fun = @(x, b) (1 + tanh(x - b))/2;
% r_fun = @(x,b) x;
i0_vec = [ones(numNeurExc, 1); zeros(numNeurInh, 1)];
b0 = 2;
r_min_fun = @(x, i0) sum((x + i0*i0_vec - param.J*r_fun(x, b0)).^2);

x_sol = zeros(param.numNeur, length(input_levels));
min_val = 0*input_levels;
for ii = 1:length(input_levels)
    r0 = i0_vec + 0.1*rand(size(i0_vec));
    [x_sol(:, ii), min_val(ii)] = fminunc(@(x) r_min_fun(x, input_levels(ii)), r0);
end
%%
r_sol = r_fun(x_sol, b0);
min_worked = min_val < 1e6;
r_sol_perI = r_sol*diag(1./input_levels);

is_exc = [true(numNeurExc, 1); false(numNeurInh, 1)];
figure()
subplot(221)
hold on
plot(input_levels(min_worked), r_sol(is_exc, min_worked)', 'b-')
ylabel('exc. responses by cell')

yyaxis right
plot(input_levels(min_worked), r_sol(~is_exc, min_worked)', 'r-')
xlabel(['input (b = ' num2str(b0) ', r = ' func2str(r_fun) ])
ylabel('inh. responses by cell')

subplot(222)
hold on
plot(input_levels(min_worked), zscore(r_sol(is_exc, min_worked), [], 2)', 'b')
plot(input_levels(min_worked), zscore(r_sol(~is_exc, min_worked), [], 2)', 'r')
xlabel(['input (b = ' num2str(b0) ', r = ' func2str(r_fun) ])
ylabel('z-scored (each cell) responses by cell')

subplot(223)
hold on
plot(input_levels(min_worked), r_sol_perI(is_exc, min_worked)', 'b')
ylabel(' responses by cell/input current (exc)')

yyaxis right
plot(input_levels(min_worked), r_sol_perI(~is_exc, min_worked)', 'r-')
xlabel(['input (b = ' num2str(b0) ', r = ' func2str(r_fun) ])
ylabel(' responses by cell/input current (inh)')

subplot(4, 2, 6)
hold on
ph = plot(input_levels(min_worked), zscore(r_sol_perI(is_exc, min_worked), [], 2)')
assignColorsToLines(ph, cool(length(ph)));
title('excitatory')

subplot(4, 2, 8)
ph = plot(input_levels(min_worked), zscore(r_sol_perI(~is_exc, min_worked), [], 2)')
assignColorsToLines(ph, spring(length(ph)))
title('inhibitory')
xlabel(['input (b = ' num2str(b0) ', r = ' func2str(r_fun) ])
ylabel('z-scored (each cell) responses by cell/input current')