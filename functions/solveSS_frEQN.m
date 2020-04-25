function r_SS = solveSS_frEQN(network_param, input_params)
% finds the steady state solution in the population model (in which the E,
% I, and X populations are summed into one element)
% % % Last updated version: 11/23/18. Not giving quantitative agreement
% with simulated firing rates, and may not be finding all the solutions to
% the FR SS equations. 

%%
% from the simulated connectivity matrix, get the locations of the
% inhibitory and directly stimulated cells
is_inh = all(network_param.J <= 0, 1)';
is_inp = network_param.input_pattern ~= 0;

% these get assigned to cell groups: X (inputs), E (exc, non-input), and I
cell_groups = {is_inp, ~is_inp & ~is_inh, is_inh};

% build up the equivalent J matrix by taking the average connection
% strength between groups of cells, and then multiplying by the number of
% cells in the SOURCE cell group
J_equiv = zeros(3);
for ii_cg = 1:3
    for jj_cg = 1:3
        J_equiv(ii_cg, jj_cg) = sum(cell_groups{jj_cg})* ...
            mean(mean(network_param.J(cell_groups{ii_cg}, cell_groups{jj_cg})));
    end
end

% compute the area under the impulse current: this is a proxy for the
% charge delivered per pulse
dt = input_params.pulse_dt;
t_vec = dt:dt:input_params.pulse_width*10;
pulse_v_t = input_params.pulse_fun(t_vec, input_params.pulse_width);
integrated_pulse = dt*sum(pulse_v_t);

% find the equivalent rate: charge*amplitude*frequency for input
% current/tau. 
input_amp_equiv = integrated_pulse*input_params.input_amplitudes{1}*input_params.input_frequency{1};

% Explore around the "equivalent amplitude" by solving the SS equations.
% This solver is not very intelligent - unresolved whether there are
% actually multiple solutions. 
input_amp_range = linspace(0, 2*input_amp_equiv, 51);
r_SS_vals = zeros(3, length(input_amp_range));
stop_vals = 0*input_amp_range;
exit_flag = 0*input_amp_range;
for i_in = 1:length(input_amp_range)

    s_input = [input_amp_range(i_in); 0 ; 0];

    fr_ss_eqn = @(r_XEI) sum((r_XEI - network_param.fun_io(J_equiv*r_XEI + s_input)).^2);

    [r_SS, stop_vals(i_in), exit_flag(i_in)] = fminunc(fr_ss_eqn, 0.12*[1; .4; 1] + 0.004*randn(3,1));
    
    r_SS_vals(:, i_in) = r_SS;
end

% plot population firing rates as input current is ramped up. 
figure()
plot(input_amp_range(exit_flag == 1), network_param.fun_io(r_SS_vals(:, exit_flag == 1)))
legend({'X', 'E', 'I'})