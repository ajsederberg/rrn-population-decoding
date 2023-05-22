function ampfreq_lookup = splineParamSweepSim(proc_simdata, input_params, network_param, sim_name)

% % old way: works with 1 repetition
% evoke_fr = mean(sum(proc_simdata.dr_out_peristim, 3), 2);
% ave_abs_evoke = mean(abs(evoke_fr), 2);

amp_vals = input_params.input_ampfreq(:, 1);
freq_vals = input_params.input_ampfreq(:, 2);
% 
% [~, a_ord] = sort(amp_vals);
% 
% amp_y = amp_vals(a_ord);
% freq_x = freq_vals(a_ord);
% aae_z = ave_abs_evoke(a_ord);
% 
% [~, f_ord] = sort(freq_x);
% amp_y = amp_y(f_ord);
% freq_x = freq_x(f_ord);
% aae_z = aae_z(f_ord);
% 
% 
freq_orig = input_params.freq_range;
amp_orig = input_params.amp_range;

% amp_mat = reshape(amp_y, length(freq_orig), length(amp_orig));
% freq_mat = reshape(freq_x, length(freq_orig), length(amp_orig));
% aae_mat = reshape(aae_z, length(freq_orig), length(amp_orig));

%% uses stim-averages 
% evoke_fr_cell_ave = mean(proc_simdata.stim_ave_dr_out(:, :, :, proc_simdata.peri_stim_time < input_params.trial_length), 4);

evoke_fr = mean(mean(proc_simdata.stim_ave_dr_out(:, :, :, proc_simdata.peri_stim_time < input_params.trial_length), 4), 3);
% first index is amplitude, second index is frequency 
amp_mat = repmat(amp_orig(:), 1, length(freq_orig));
freq_mat = permute(repmat(freq_orig(:), 1, length(amp_orig)), [2 1]);
aae_mat = evoke_fr;

evokefr_fun = @(ff, aa) interp2(freq_mat, amp_mat, aae_mat, ff, aa, 'spline');
%%

[freq_x, amp_y] = meshgrid(freq_orig(1):.001:freq_orig(end), amp_orig(1):.01:amp_orig(end));

evoke_fr_mesh = evokefr_fun(freq_x, amp_y);

makeMyFigure(30,12);
subplot(121)
% imagesc(freq_x(1, :), amp_y(:, 1), evoke_fr_mesh);
contourf(freq_x, amp_y, evoke_fr_mesh, 10)
xlabel('frequency')
ylabel('amplitude')
ch = colorbar;
ch.Label.String = 'firing rate';
set(gca, 'ydir', 'normal')
axis square
title('contour plot of average evoked network firing rate')
hold on

v_fr = linspace(.2*max(evoke_fr_mesh(:)),.8*max(evoke_fr_mesh(:)), 7);
%%
contour_lines = cell(length(v_fr),1 );
keep_vfr = 0*v_fr;
ctr = 0;
for i_c = 1:length(v_fr)

    C_mat = contourc(freq_x(1, :), amp_y(:, 1), evoke_fr_mesh, v_fr(i_c)*[1 1]);
    % if any of C_mat is out of bounds of the original amplitude range,
    % drop it
    out_of_amp_range = (C_mat(2, :) > max(amp_vals) | C_mat(2, :) < min(amp_vals) ... 
        | C_mat(1, :) > max(freq_vals) | C_mat(1, :) < min(freq_vals));
    C_mat = C_mat(:, ~out_of_amp_range);
    
    if ~isempty(C_mat)
        ctr = ctr + 1;
        contour_lines{ctr} = C_mat;
        keep_vfr(ctr) = v_fr(i_c);
        plot(C_mat(1, :), C_mat(2, :), '.', 'MarkerSize', 10)
    end
end
contour_lines = contour_lines(1:ctr);
v_fr = keep_vfr(1:ctr);

xlim(freq_x([1 end]))
ylim(amp_y([1 end]))
% % as a function of frequency, find the amplitude at which 2*freq has the
% % same firing rate
%%
target_freqs = 0.1:0.01:0.25;
amp2_vals = nan(length(contour_lines), length(target_freqs));
amp1_vals = nan(length(contour_lines), length(target_freqs));

% target_freq1 = 0.16;


for i_freq = 1:length(target_freqs)
    target_freq1 = target_freqs(i_freq);
    for i_c = 1:length(contour_lines)

        % check if there is a value on the contour satisfying freq =
        % 2*target_freq1
        freq_c = contour_lines{i_c}(1, :);
        amp_c = contour_lines{i_c}(2, :);
        
        %% first, some contours are disconnected. find connected segments of
        % the contour that have the full frequency range
            % find monotonically increasing segments of freq_c
        df_c = sign(diff(freq_c));
        df_c_switches =[0 ( find(diff(df_c) ~= 0) + 1) length(freq_c)];
        keep_seg = false(length(df_c_switches) - 1, 1);
        seg_L = 0*keep_seg;
        for i_seg = 1:length(df_c_switches)-1
            freq_seg = freq_c(df_c_switches(i_seg)+1 : df_c_switches(i_seg+1));
            keep_seg(i_seg) = target_freq1 >= min(freq_seg) ...
                && 2*target_freq1 <= max(freq_seg);
            seg_L(i_seg) = length(freq_seg);
        end
        
        if sum(keep_seg) > 1
            % take the longest segment
            [~, ind] = max(seg_L(keep_seg));
            keep_seg = find(keep_seg);

            keep_seg_ind = keep_seg(ind);
        else
            keep_seg_ind = find(keep_seg);
        end
        
        seg_inds = df_c_switches(keep_seg_ind)+1 : df_c_switches(keep_seg_ind+1);
        freq_c = freq_c(seg_inds);
        amp_c = amp_c(seg_inds);
        %%
        
        try amp_2freq1 = interp1(freq_c, amp_c, 2*target_freq1, 'linear', 'extrap');
            amp_freq1 = interp1(freq_c, amp_c, target_freq1, 'linear', 'extrap');
        catch
            amp_2freq1 = nan;
        end

        if ~isnan(amp_2freq1)
            amp2_vals(i_c, i_freq) = amp_2freq1;
            amp1_vals(i_c, i_freq) = amp_freq1;
        end
    end
end
%%
subplot(2, 2, 2);
hold on
ph = plot(target_freqs, amp2_vals, 'linewidth', 1.5);
ph = assignColorsToLines(ph, lines(length(ph)));

legend(ph, num2str(v_fr', 'fr = %1.2f'), 'Location', 'eastoutside')
xlabel('freq1')
ylabel('amplitude at 2*freq1')
ylim(amp_y([1 end]))
title('look-up values of amplitude of input w/ frequency 2*freq1')
subplot(2, 2, 4);
hold on
ph = plot(target_freqs, amp1_vals, 'linewidth', 1.5);
ph = assignColorsToLines(ph, lines(length(ph)));

legend(ph, num2str(v_fr', 'fr = %1.2f'), 'Location', 'eastoutside')
xlabel('freq1')
ylabel('amplitude at freq1')
ylim(amp_y([1 end]))
title('look-up values of amplitude of input w/ frequency freq1')

% print(gcf, '-dpdf', [network_param.plotdir 'amplitude_lookup_' sim_name])
exportgraphics(gcf, [network_param.plotdir 'amplitude_lookup_' sim_name '.pdf'])


% save in the matsfile dir 
saveInParfor([network_param.matsdir 'amplitude_lookup_' sim_name], ...
    amp_vals, freq_vals, evokefr_fun, evoke_fr, amp1_vals, ...
    amp2_vals, target_freqs, v_fr, contour_lines, ...
    freq_x, amp_y, evoke_fr_mesh, input_params);
% save in the save_dir if this points somewhere other than the matfiles dir
if ~strcmp([network_param.save_dir 'matfiles/'], network_param.matsdir)
    if ~exist([network_param.save_dir 'matfiles/'], 'dir')
        mkdir([network_param.save_dir 'matfiles/']);
    end
    saveInParfor([network_param.save_dir 'matfiles/amplitude_lookup_' sim_name], ...
        amp_vals, freq_vals, evokefr_fun, evoke_fr, amp1_vals, ...
        amp2_vals, target_freqs, v_fr, contour_lines, ...
        freq_x, amp_y, evoke_fr_mesh, input_params);
end

ampfreq_lookup.amp1_vals = amp1_vals;
ampfreq_lookup.amp2_vals = amp2_vals;
ampfreq_lookup.freq1_vals = target_freqs;
ampfreq_lookup.freq2_vals = 2*target_freqs;
ampfreq_lookup.fr_vals = v_fr;
ampfreq_lookup.input_param = input_params;