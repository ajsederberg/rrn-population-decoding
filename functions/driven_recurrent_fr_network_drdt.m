function dx_dt = driven_recurrent_fr_network_drdt(t, x, J_mat, nl_fun, input_fun, input_unit_pattern)
% simulates a recurrent network, dx/dt = -x + Jr + inp and r = f(x) for
% some nonlinear function provided in nl_fun. (This function was named
% with 'drdt' mistakenly: it should be dxdt.)

I_in = J_mat*nl_fun(x) + input_fun(t)*input_unit_pattern;

dx_dt = -x + I_in;