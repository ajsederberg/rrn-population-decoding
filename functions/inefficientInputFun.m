function in_vals = inefficientInputFun(t, input_times, input_amplitudes, pulse_fun, pulse_width)

relevant_times = cellfun(@(x) t > x & (t - x) < 30*pulse_width, input_times, ...
    'uniformoutput', false);
relevant_stim = cellfun(@(x) any(x), relevant_times);

use_input_times = cellfun(@(x,ind) x(ind), input_times(relevant_stim), relevant_times(relevant_stim), ...
    'uniformoutput', false);
use_amplitudes = input_amplitudes(relevant_stim);

if sum(relevant_stim) > 0
    in_vals = sum(cell2mat(cellfun(@(in_time, in_pat) sum(pulse_fun(t - in_time, pulse_width), 2).*in_pat, ...
        use_input_times, use_amplitudes, 'uniformoutput', false)), 2);
else
    in_vals = 0;
end