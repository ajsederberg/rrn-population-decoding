function dr_dt = rect_linear_network_drdt(t, r, W_mat, thr_theta, input_t, input_mat)

t_ind = find(t > input_t(:, 1) & t < input_t(:, 2));
if ~isempty(t_ind)
    inp = input_mat(t_ind, :)';
else
    inp = 0*r;
end


I_in = W_mat*r + thr_theta + inp;
I_in(I_in < 0) = 0;

dr_dt = -r + I_in;