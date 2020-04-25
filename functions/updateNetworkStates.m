function [x_activity, s_state] = updateNetworkStates(t, x_activity, s_state, inp_current,  J_ij, J_state, h_state, T_state, num_flops)
% Inputs:
% t : simulation time step 
% x_activity: current activity pattern 
% s_state: current state pattern (+/- 1)
% inp_current: handle to a function that generates a vector (of size 
% size(x_activity)) that provides input current to the network
% J_ij: recurrent connectivity between neurons 
% J_state: interaction matrix between states 
% h_state: handle to a function that generates a vector (of size 
% size(x_activity)) that provides the value of the external field 
% T_state: units in which to measure energy differences of state vectors
% num_flops: number of sites in the state vector to flip on a dynamical
% update step
% State changes are simulated with Metropolis dynamics, meaning that a
% state change is accepted if it moves to a more probable state and
% accepted probabilistically if it moves to a less probable state (with the
% probability set by the factor P_old/P_new
% Activity is McCulloch-Pitts neurons updated on each step with a
% time-dependent input current and recurrent connectivity matrix. 


 
% First update states
E_oldstate = s_state'*J_state*s_state + h_state(t)'*s_state;
flop_inds = randperm(length(s_state), num_flops);
old_state = s_state; 
s_state(flop_inds) = -1*s_state(flop_inds);
E_newstate = s_state'*J_state*s_state + h_state(t)'*s_state;

% reject the change probabilistically if the energy change is
% positive
if E_oldstate - E_newstate < 0 && rand < exp(((E_oldstate-E_newstate)/num_flops)/T_state)
    s_state = old_state;
end 

% Compute activity
x_activity = (s_state+1).*(inp_current(t) + J_ij*x_activity);