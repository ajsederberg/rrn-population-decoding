function dx_dt = recurrent_fr_network_drdt(t, x, J_mat, nl_fun, input_t12, input_mat)
% simulates a recurrent network, dr/dt = -r + Jx + inp and x = f(r) for
% some nonlinear function provided in nl_fun. 
t_ind = find(t > input_t12(1, :) & t < input_t12(2, :));

if ~isempty(t_ind)
    inp = input_mat(:, t_ind);
else
    inp = 0*x;
end


I_in = J_mat*nl_fun(x) + inp;

dx_dt = -x + I_in;