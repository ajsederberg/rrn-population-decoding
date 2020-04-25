function y = stateNetworkStimulationFunction(t, stim_times, stim_funs)
% t is the current simulation time (an integer)
% cell arrays stim_times and stim_funs represent stimulus presentations
% each element of the cell array stim_times are the times corresponding to
% that presentation ("trial")
% each element of the cell array stim_funs the function giving the stimulus
% input at stim_funs(t_index)
% For example: if stim_times{4} = 41:50, and stim_funs{4} = @(x) 10-x, then
% stateNetworkStimulationFunction(44, stim_times, stim_funs) = 6; 
% because 44 is the 4th element of (41:50) and 10-4 = 6. 

[which_stim, t_index] = cellfun(@(x) ismember(t, x), stim_times);

y = stim_funs{which_stim}(t_index(which_stim));