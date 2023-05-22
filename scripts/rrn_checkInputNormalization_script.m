sim_name = 'spEnoise_g300_n500_ei100_rng1';

load(['results/' sim_name '/matfiles/sim_data_' sim_name '_AFsweep_pw50.mat'])
%% plot a few inputs trains
figure()
hold on

for i_tr = 1:10
    tt1 = input_params.trial_start_times(i_tr);
    t_len = (0:0.1:input_params.trial_length)';


    ph = plot(t_len, input_params.input_fun(tt1 + t_len));
    plot(t_len, mean(input_params.input_fun(tt1 + t_len)) + 0*t_len, 'color', ph.Color)
end

%% check over trials with same input - how variable??
trial_mean_resp = squeeze(mean(proc_simdata.stim_ave_dr_out, 4));
i_a1 = 4;
j_f1 = 2;
resp_ave1 = squeeze(trial_mean_resp(i_a1, j_f1, :));
f1_lab = ['f = ' num2str(input_params.freq_range(j_f1)) ', a = ' num2str(input_params.amp_range(i_a1))];

i_a2 = 2;
j_f2 = 4;
resp_ave2 = squeeze(trial_mean_resp(i_a2, j_f2, :));
f2_lab = ['f = ' num2str(input_params.freq_range(j_f2)) ', a = ' num2str(input_params.amp_range(i_a2))];
is_exc = network_param.input_pattern;
figure()
subplot(121)
hold on
plot(resp_ave1(is_exc), resp_ave2(is_exc), '.')
plot(resp_ave1(~is_exc), resp_ave2(~is_exc), '.')

eqline
xlabel(f1_lab)
ylabel(f2_lab)
subplot(122)
hold on
histogram(resp_ave1)
histogram(resp_ave2)
legend({f1_lab, f2_lab})

%%

is_s_fun = @(a_ind, f_ind) input_params.input_ampfreq(:, 1) == input_params.amp_range(a_ind) & ...
    input_params.input_ampfreq(:, 2) == input_params.freq_range(f_ind);

is_s1 = is_s_fun(i_a1, j_f1);
is_s2 = is_s_fun(i_a2, j_f2);

std_resp_s1 = squeeze(std(mean(proc_simdata.dr_out_peristim(is_s1,:, :), 3), [], 1));
std_resp_s2 = squeeze(std(mean(proc_simdata.dr_out_peristim(is_s2,:, :), 3), [], 1));

mean_resp_s1 = squeeze(mean(mean(proc_simdata.dr_out_peristim(is_s1,:, :), 3), 1));
mean_resp_s2 = squeeze(mean(mean(proc_simdata.dr_out_peristim(is_s2,:, :), 3), 1));

figure()
hold on
histogram((mean_resp_s1(is_exc) - mean_resp_s2(is_exc))./sqrt(std_resp_s1(is_exc).^2 +std_resp_s2(is_exc).^2), 'Normalization', 'pdf')
histogram((mean_resp_s1(~is_exc) - mean_resp_s2(~is_exc))./sqrt(std_resp_s1(~is_exc).^2 +std_resp_s2(~is_exc).^2), 'Normalization', 'pdf')

%% check eigenvalues of J
param = network_param;
figure()
[v,l] = eig(param.J);
subplot(211)
plot(diag(l), 'o')
title(['c0: ' num2str(1) ', gEE: ' num2str(param.gEE) ', gEI: ' num2str(param.gEI) ... 
    ', gIE: ' num2str(param.gIE) ', gII: ' num2str(param.gII)])
axis equal
subplot(212)
plot(v(:, 1))
hold on
vE = mean(v(1:param.numNeurExc, 1));
vI = mean(v(1+param.numNeurExc:end, 1));
plot([1 param.numNeurExc], vE*[1 1])
plot([1+param.numNeurExc param.numNeur], vI*[1 1])
title(['first eigenvector, vE = ' num2str(vE, '%1.2f') ...
    ', vI = ' num2str(vI, '%1.2f') ', ratio: ' num2str(vI/vE, '%1.2f')])

%% estimate input-responses and compare to theh actual sweeps
param = network_param;
input_levels = linspace(0.5, 3, 11);

r_fun = @(x, b) (1 + tanh(x - b))/2;
% r_fun = @(x,b) x;
i0_vec = [ones(param.numNeurExc, 1); zeros(param.numNeurInh, 1)];
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

is_exc = [true(param.numNeurExc, 1); false(param.numNeurInh, 1)];
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