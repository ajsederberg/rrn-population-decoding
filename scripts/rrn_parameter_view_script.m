% Script to look at the parameters used in the simulations in the paper. 
num_neurons = [500];
inh_netwt_frac = [0.5];
network_gain = [4 ]; % 4 or 6 OK

[nn_g, in_g, ng_g] = meshgrid(num_neurons, inh_netwt_frac, network_gain);
%  network_gain, network_index, inh_netwt_frac
%   connect codes    1          2           3       4           5           6           7           8   
% connect_tags = {'logNorm', 'splogNorm', 'norm', 'spnorm', 'spnormSL', 'spnorm2P', 'norm2PSL', 'spnormB'};
% Note that 'spnorm' takes the input cells out of the exc pool, but doesn't
% renormalized the e-i balance in the rest of the network. 'spnormB'
% rectifies this, by scaling up the excitatory connections. 
network_codes = encodeNetworkMetaparams(nn_g, ng_g, 4 + 0*nn_g, in_g);

network_codes = network_codes(:);

check_pars = cellfun(@(x) generateNetworkParam(x, 1), num2cell(network_codes), ...
    'UniformOutput', false);
check_codes = cellfun(@(x) x.network_tag, check_pars, 'UniformOutput', false);

pulse_width_vals = [ .5 1.5 5];
network_code = network_codes(1);

% rng_vals = [1:12 17];
rng_vals = [1 2 3 6 7 11 12];
all_network_param = cellfun(@(x) lookUpNetworkParam(network_code, x), ... 
    num2cell(rng_vals), 'UniformOutput', false);
%% load all the AUC stats
pulse_width = pulse_width_vals(1);
all_aucstats = cell(length(rng_vals), 1);
stimulus_experiment_tag = ['f2fdiscExp_pw' num2str(100*pulse_width)];
for i_rv = 1:length(rng_vals)
    new_rng_str = ['rng' num2str(rng_vals(i_rv))];
    network_param = all_network_param{i_rv};
    auc_file = [ network_param.save_dir 'matfiles/selXCJ_stats_' ...
        network_param.network_tag stimulus_experiment_tag];

    
    rng_locs = regexp(auc_file, 'rng\w*', 'split');
    new_net_auc_file = [rng_locs{1} new_rng_str rng_locs{2} new_rng_str stimulus_experiment_tag];
    all_aucstats{i_rv} = load(new_net_auc_file);
end
sel_bias = cellfun(@(x) sum(x.all_cell_selectivity)/sum(x.all_cell_selectivity ~= 0), all_aucstats);
%%
gain_neur = network_gain*100;
% rng_val = 3;
ei_ratio = 100; % 92 or 108;  80 or 125
pulse_width = .50; % 50 or 150 or 500
connect_tag = 'spnorm'; % 'spnorm2P' for two-pool
f2f_type = 'f2f';
file_tag = [connect_tag '_g' num2str(gain_neur) '_n' ...
    num2str(num_neurons) '_ei' num2str(ei_ratio) '_rng'];
input_type = [f2f_type 'discExp_pw' num2str(100*pulse_width)];

all_input_param = cellfun(@(x) lookUpInputParam(network_code, x, ...
    [file_tag num2str(x) '_' input_type]), num2cell(rng_vals), 'UniformOutput', false);
%% load classifier results

load(['multiRNG_' file_tag '3_20190821.mat'], 'all_cls_res')
%%
amp1_vals = cellfun(@(x) x.amp_range(1), all_input_param);
amp2_vals = cellfun(@(x) x.amp_range(2), all_input_param);


%% 
ampfreq_lookup = cell(length(all_network_param), 1);
for ii = 1:length(all_network_param)
    ampfreq_simname = [all_network_param{ii}.network_tag ...
        '_AFsweep_pw' num2str(round(100*pulse_width))];

    ampfreq_lookup{ii} = load([all_network_param{ii}.matsdir '/amplitude_lookup_' ampfreq_simname], ...
        'amp1_vals', 'amp2_vals', 'target_freqs', 'v_fr', 'amp_vals');
end
%%
targetfreq_lookup = cellfun(@(x) x.amp1_vals(:, x.target_freqs == 0.16), ampfreq_lookup, ... 
    'UniformOutput', false);

fr_sim = cellfun(@(af_str, all_amp1, amp1) af_str.v_fr(all_amp1 == amp1), ... 
    ampfreq_lookup, targetfreq_lookup, num2cell(amp1_vals'));
%%
figure()
cmap = cool(2);
hold on
for ii = 1:length(all_network_param)
    histogram(nonzeros(all_network_param{ii}.J), 'DisplayStyle', 'stairs', 'Normalization', 'cdf')
end
%% 
popAcc = cellfun(@(x) max(offDiag(x{1}.classIdataJ_valAcc)), all_cls_res); 

%%
figure()
subplot(121)
hold on
plot(sel_bias, amp1_vals, 'o')
plot(sel_bias, amp2_vals, 'o')
% plot(sel_bias, 2*amp2_vals, 'o')

subplot(122)
plot(sel_bias, fr_sim, 'o')

%% Display parameters for methods

mean_Jexc = mean(cellfun(@(x) mean(x.J(x.J > 0)), all_network_param))
std_Jexc = mean(cellfun(@(x) std(x.J(x.J > 0)), all_network_param))

mean_Jinh = mean(cellfun(@(x) mean(x.J(x.J < 0)), all_network_param))
std_Jinh = mean(cellfun(@(x) std(x.J(x.J < 0)), all_network_param))


%% full width half max of input pulse and input amplitudes
t = -2:.1:10; 
y = all_input_param{1}.pulse_fun(t', 0.5)

first_nonzeros = find(y > 0, 1)
[ymax, ind_max] = max(y)

t1 = interp1(y(first_nonzeros:ind_max), t(first_nonzeros:ind_max), 0.5*y(ind_max));


t2 = interp1(y(ind_max:end), t(ind_max:end), 0.5*y(ind_max));

pf_fun2 = @(t,a)(t.^2/a^2).*exp((t>0).*(-t/a)).*(t>0);

makeMyFigure(30, 15);
subplot(121)
hold on
plot(t, all_input_param{1}.pulse_fun(t', 0.5))
plot(t, pf_fun2(t', 0.5))

plot([t1 t2], all_input_param{1}.pulse_fun([t1 t2]', 0.5), 'o')
title(['full width half max : ' num2str(t2 - t1, '%1.2f')])

subplot(122)
set(gca, 'visible', 'off')
hold on

for i_net = 1:length(all_input_param)
    inp_amps = exp(1)*cell2mat(all_input_param{i_net}.input_amplitudes);
    
    text(0, i_net/length(all_input_param), ['net ' num2str(i_net) ' inputs: ' ... 
        num2str(inp_amps, '%1.1f  ')])
end
%%
all_inp_amps = exp(1)*cell2mat(cellfun(@(x) cell2mat(x.input_amplitudes), all_input_param', 'UniformOutput', false));

mean_inputs = mean(all_inp_amps)
std_inputs = std(all_inp_amps)

mean_ratio = mean(all_inp_amps(:, 1)./all_inp_amps(:, 2))
std_ratio = std(all_inp_amps(:, 1)./all_inp_amps(:, 2))
%% plot histogram of nonzero synaptic weights, distribution of firing rates, 
% and classifier accuracy for each network
stim_colors = lines(4);
stim_colors = stim_colors(3:4, :);
neuron_colors = [(1 + lines(1))/2; lines(2)];
samediff_colors = [0 0 0; 1 1 1]*0.5;
for ii = 1:length(all_network_param)
makeMyFigure(30, 10);
net_J = all_network_param{ii}.J;
is_inp = all_network_param{ii}.input_pattern;
is_inh = any(net_J < 0, 1)';

auc_J = all_aucstats{ii}.network_param.J;
disp([num2str(mean(auc_J(:) == net_J(:))*100, '%1.0f') '% of J vals match across auc and par files'])
% netInps = sum(all_network_param{ii}.J(:, is_inp), 2);
% netExc = sum(all_network_param{ii}.J(:, ~is_inp & ~is_inh), 2);
% netInh = sum(all_network_param{ii}.J(:, is_inh), 2);
% 
subplot(1, 4, 1)
hold on
histogram(nonzeros(net_J(:, ~is_inh)), 'Normalization', 'pdf', 'EdgeColor', 'none')
histogram(nonzeros(net_J(:, is_inh)), 'Normalization', 'pdf', 'EdgeColor', 'none')
% histogram(netExc)
% histogram(netInh)
fr_bins = linspace(-0.1, 1, 101);

axis square
set(gca, 'color', 'none')
legend({'e-weight hist', 'i-weight hist'}, 'location', 'northoutside')
xlabel(['synaptic weights (' num2str(100*mean(net_J(:)~=0), '%1.1f') '% nonzero)'])
ylabel('pdf')

subplot(1, 4, 2)
hold on
histogram(mean(all_aucstats{ii}.stim1_endPoint, 2),fr_bins,'normalization', 'cdf', 'linewidth', 1.5, 'DisplayStyle', 'stairs', 'EdgeColor', stim_colors(1, :))
histogram(mean(all_aucstats{ii}.stim2_endPoint, 2),fr_bins,'normalization', 'cdf', 'linewidth', 1.5, 'DisplayStyle', 'stairs', 'EdgeColor', stim_colors(2, :))
ylabel('cdf')
xlabel('firing rates')
axis square
set(gca, 'color', 'none')
legend({'ave fr (stim 1)', 'ave fr (stim 2)'}, 'location', 'northoutside')

subplot(1, 4, 3)
hold on
histogram(mean(all_aucstats{ii}.all_endPoint(:, ~is_inh & ~is_inp), 1),fr_bins,'normalization', 'cdf', 'linewidth', 1.5, 'DisplayStyle', 'stairs', 'EdgeColor', neuron_colors(2, :))
histogram(mean(all_aucstats{ii}.all_endPoint(:, is_inp), 1),fr_bins,'normalization', 'cdf', 'linewidth', 1.5,'DisplayStyle', 'stairs', 'EdgeColor', neuron_colors(1, :))
histogram(mean(all_aucstats{ii}.all_endPoint(:, is_inh), 1),fr_bins,'normalization', 'cdf', 'linewidth', 1.5,'DisplayStyle', 'stairs', 'EdgeColor', neuron_colors(3, :))
ylabel('cdf')
xlabel('firing rates')
axis square
set(gca, 'color', 'none')
legend({'ave fr (e-cells)', 'ave fr (input e-cell)', 'ave fr (i-cells)'}, 'location', 'northoutside')

subplot(1, 4, 4)
% plot mean +/- SE of test accuracy
x_groups = [1 3 4];
x_tic_labels = cellfun(@(x) x.cellGroup, all_cls_res{ii}(x_groups), 'uniformoutput', false);
x_mean = cellfun(@(x) mean(x.cvDraws_testAcc), all_cls_res{ii}(x_groups));
x_se = cellfun(@(x) std(x.cvDraws_testAcc), all_cls_res{ii}(x_groups));
errorbar(1:length(x_groups), x_mean, x_se, 'o', 'linewidth', 1.5)
ylim([0.5 1])
xlim([0 4])
set(gca, 'xtick', 1:length(x_groups), 'XTickLabel', x_tic_labels)
axis square
set(gca, 'color', 'none')
% pop_acc_bins = linspace(0.5, 1, 21);
% hold on
% histogram(all_cls_res{ii}{1}.cvDraws_testAcc, pop_acc_bins)
% histogram(all_cls_res{ii}{3}.cvDraws_testAcc, pop_acc_bins)
% histogram(all_cls_res{ii}{4}.cvDraws_testAcc, pop_acc_bins)
ylabel('population decode accuracy')
xlabel('cell groups')
legend([], 'location', 'northoutside', 'visible', 'off')

hsup = suptitle(all_network_param{ii}.network_tag);
set(hsup, 'interpreter', 'none')

print(gcf, '-dpdf', ['plots/' all_network_param{ii}.network_tag '_param_view'])

end
%% net synaptic inputs from X pop

is_inp = all_network_param{ii}.input_pattern;
is_inh = any(all_network_param{ii}.J < 0, 1)';
netInps = sum(all_network_param{ii}.J(:, is_inp), 2);
netExc = sum(all_network_param{ii}.J(:, ~is_inp & ~is_inh), 2);
netInh = sum(all_network_param{ii}.J(:, is_inh), 2);
%%
figure()
hold on
histogram(netInps)
histogram(netExc)
histogram(netInh)
legend({'inputs (e)', 'non-input E', 'inh'})
