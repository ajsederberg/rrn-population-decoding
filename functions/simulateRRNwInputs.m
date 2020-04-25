function [sim_data_table, stim_param] = simulateRRNwInputs(param, stim_param)
% Wrapper function to simulate the recurrent firing rate network defined by
% the 'recurrent_fr_network_drdt'. This version uses the new stimulus
% generation methods. For the old version, see
% simuluateRRNwInputs_legacyStims
% Params extracted below. 


% load parameters
J = param.J;
fun_io = param.fun_io;


% x0 is in the initial condition of neurons; typically have a burn-in
% period so not highly relevant. 
x0 = 0*fun_io(randn(size(J, 1), 1));
% add a tiny bit of noise to IC's
x0_sigma = 0.1;


% numNeurExc = param.numNeurExc;
% numNeurInh = param.numNeurInh;
% numNeur = size(J, 1);

% these are the start times of each stimulus (stim 1, stim 2). These are
% used for calculating STAs and peri-stimulus responses.
% input_t = param.input_t;
% inp1_times = param.inp1_times;
% inp2_times = param.inp2_times;

% ode_fun = @(t, y) recurrent_fr_network_drdt(t, y, J, fun_io, input_times, input_mat);


% %  input_mat (if all inputs are spatially identical) can be appended at
% the ned; otherwise, it should enter in the middle of this. That'd be a
% large matrix multiplication at each step, which is better avoided if
% possible... 
% input_fun = @(t) sum(cell2mat(cellfun(@(in_time, in_pat) sum(exp(-(in_time - t).^2/(2*pulse_width^2)), 2).*in_pat, ...
%     input_times, input_mat, 'uniformoutput', false)), 2);

% input fun is a function that returns, as a function of time (t), the
% input to each of the N neurons at that point in time. 

%%%%%%%%%% uncomment this line to go back to old version
% stim_param needs the number of neurons
stim_param.numNeur = length(x0);
[input_fun, input_times, input_amplitudes, t0tf, stim_param] = generateRRNInputFunction(stim_param);

% %%%%%%%%%%%%%%% testing adaptation
% [input_fun, input_times, input_amplitudes, t0tf, stim_param] = generateRRNAdaptingInputFunction(stim_param);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stim_param.input_times = input_times;
stim_param.input_amplitudes = input_amplitudes;
stim_param.input_fun = input_fun;
% t_span is [start_time stop_time] for the simulation. 
t_span = t0tf;

% input_fun = @(t) sum(exp(-(pulse_times - t).^2/(2*pulse_width^2)))*input_mat(:, 1);
input_pattern = param.input_pattern;
ode_fun = @(t, y) driven_recurrent_fr_network_drdt(t, y, J, fun_io, input_fun, input_pattern);

options = odeset('RelTol', 1e-3, 'AbsTol', 1e-3);

% First integrate to get initial conditions
[~,x_outIC] = ode45(ode_fun, [0 stim_param.burn_in_time], x0 + x0_sigma*randn(size(x0)), options);
x0 = x_outIC(end, :);

% % noise check
% dt = 0.01;
% sig2 = stim_param.sig2;
% [~, noisy_x_outIC2] = myRK2SODESolver(ode_fun, [0 stim_param.burn_in_time], ...
%     x0 + x0_sigma*randn(size(x0)), dt, sig2);
% 

% r_out = fun_io(x_out);
% sim_data_table = table(t_out, x_out, r_out);

% Now integrate around each stimulus time, starting from the same IC. 
num_intervals = 1 + length(stim_param.trial_start_times);
t_out = cell(num_intervals, 1);
x_out = cell(num_intervals, 1);

if ~contains(param.network_tag, 'spEnoise')
    % integrate over the stimulus periods
    t_start = 0;

    for i_int = 1:num_intervals
        t_interval = t_start + stim_param.trial_length*[-0.1 1.2];
        [t_out{i_int}, x_out{i_int}] = ode45(ode_fun, t_interval, x0, options);

    %     dt = 0.01;
    %     [t_out{i_int}, x_out{i_int}] = myODESolver(ode_fun, t_interval, x0, dt);


        % reset the initial conditions
    %     x0 = x_out{i_int}(end, :)' + x0_sigma*randn(size(x0));

        if i_int < num_intervals
            % because the first interval is the burn-in period, i_int = 2
            % corresponds to trial_start_times(1). No need to reset t_start on
            % the last trial (Also, trial_start_times is only
            % length(num_intervals-1), by construction)
            t_start = stim_param.trial_start_times(i_int);
        end

    end

else
    % integrate over the stimulus periods
    t_start = 0;

    for i_int = 1:num_intervals
        t_interval = t_start + stim_param.trial_length*[-0.1 1.2];
    %     [t_out{i_int}, x_out{i_int}] = ode45(ode_fun, t_interval, x0, options);
%%
        dt = 0.02;  % was 0.002, but trial is 50 and that is excessive.
        sig2 = stim_param.sig2;
%         [t_out{i_int}, x_out{i_int}] = myRK2SODESolver(ode_fun, t_interval, x0, dt, sig2);
        [tt, xx] = myRK2SODESolver(ode_fun, t_interval, x0, dt, sig2);
%%
        % return this sampled at dt = 0.1; otherwise have out-of-memory
        % errors
        dt_sample = 0.1;
        t_out{i_int} = (tt(1):dt_sample:tt(end))';
        x_out{i_int} = interp1(tt, xx, t_out{i_int});

        % reset the initial conditions
    %     x0 = x_out{i_int}(end, :)' + x0_sigma*randn(size(x0));

        if i_int < num_intervals
            % because the first interval is the burn-in period, i_int = 2
            % corresponds to trial_start_times(1). No need to reset t_start on
            % the last trial (Also, trial_start_times is only
            % length(num_intervals-1), by construction)
            t_start = stim_param.trial_start_times(i_int);
        end
        
        if i_int == 1 % save the full data for the first 'segment' (pre-stimulus)
            save([param.matsdir 'baselineNoise_' param.network_tag], 'tt', 'xx', 'stim_param', 'param');
        end
    end

    
end

% r_out = cellfun(@(x) fun_io(x), x_out, 'uniformoutput', false);

t_out = cell2mat(t_out);
x_out = cell2mat(x_out);
% r_out = cell2mat(r_out);
sim_data_table = table(t_out, x_out);

% % euler method integration. 
% int_delta_t = (input_times(2,1) - input_times(1, 1))/10;
% [TOUT_M, XOUT_M] = myODESolver(ode_fun, t_span, y0, int_delta_t);
%%

