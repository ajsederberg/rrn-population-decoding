% Make Figure 4: results from simulations with inputs to excitatory units
% only. 
% num_neur = network_param.numNeur;
% file_tag = network_param.network_tag; %

%%%% or set values individually, if running as a script
num_neur = 500;
gain_neur = 200; % usually 400
rng_val = 2; % usually 2
ei_ratio = 100; % 92 or 108;  80 or 125
pulse_width = .50; % 50 or 150 or 500
connect_tag = 'spnormE'; % 
f2f_type = 'f2f';

file_tag = [connect_tag '_g' num2str(gain_neur) '_n' num2str(num_neur) '_ei' num2str(ei_ratio) '_rng' num2str(rng_val)];

% replace with path to the simulation results
local_res_dir = '/Volumes/Samsung_T5/projects/random_rec_network_project/results/'; % classifier results are here


input_type = [f2f_type 'discExp_pw' num2str(100*pulse_width)];
%% required for Fig 1
stimulus_experiment_tag = [f2f_type 'discExp_pw' num2str(100*pulse_width)];
x_res = load([local_res_dir file_tag '/matfiles/sim_data_' file_tag '_' stimulus_experiment_tag '.mat']);

%%
proc_simdata = x_res.proc_simdata;

input_params = x_res.input_params;

network_param = x_res.network_param;


%% load classifier processed results. 
classifier_file_tag = ['newclassifierAEIS_' file_tag '_' input_type];

x_class = load([local_res_dir file_tag '/matfiles/' classifier_file_tag]);
x_act = load([local_res_dir file_tag '/matfiles/' ... 
            'stim_resps_' input_type], 'activity_struct');
        
        
activity_struct = x_act.activity_struct;



%% load all available classifier results
% if true

if ~exist(['multiRNG_' file_tag '_20191220.mat'], 'file')
rng_vals = [2 3 4 5 6 7] ; %[1 2 3 4 5 6 7 8 9 10 11 12 13 17] ; %[1 2 3 4 5 6 7 8 9 10];
all_cls_res = cell(length(rng_vals), 1);
clear all_cls_ei_inds
all_cls_ei_inds(length(rng_vals)) = struct;
all_varCE = cell(length(rng_vals), 1);

for i_rv = 1:length(rng_vals)
    
%%
        new_rng_str = ['rng' num2str(rng_vals(i_rv))];
        rng_locs = regexp(file_tag, 'rng\w*', 'split');
        new_file_tag = [rng_locs{1} new_rng_str rng_locs{2}];
        new_classifier_file_tag = [local_res_dir new_file_tag '/matfiles/newclassifierAEIS_' ...
            new_file_tag '_' input_type '.mat'];
        new_stim_resps_tag = [local_res_dir new_file_tag '/matfiles/' ... 
            'stim_resps_' input_type '.mat'];
        xxcr = load(new_classifier_file_tag, 'cls_res', 't_vec');
        x_act = load(new_stim_resps_tag, 'activity_struct');
        x_act = x_act.activity_struct;
        all_cls_res{i_rv} = xxcr.cls_res;
      %%  
       
        all_varCE{i_rv} = computeVarCE_rateModel(x_act.stim_resps, x_act.stim_y);

        % keep record of which cells are excitatory and inhibitory, and
        % which are directly stimulated or not. 
        postl1_inds = find(x_act.postL1_cells);
        inh_inds = find(x_act.inh_cells);
        [~, inh_inds] = intersect(postl1_inds, inh_inds);
        exc_inds = setdiff(1:length(postl1_inds), inh_inds);
        all_cls_ei_inds(i_rv).exc_inds = exc_inds;
        all_cls_ei_inds(i_rv).inh_inds = inh_inds;
        
         % compute average firing rates
        all_cls_ei_inds(i_rv).pop_ave_fr = mean(x_act.stim_resps(:, :, 10:50), 3);
        all_cls_ei_inds(i_rv).stim_id = x_act.stim_y;
 
end

% save
save(['multiRNG_' file_tag '_20191220'], 'all_cls_res', 'all_cls_ei_inds', 'all_varCE', 'rng_vals')
%%
else
%     load(['multiRNG_' file_tag], 'all_cls_res_noVA');
    load(['multiRNG_' file_tag '_20191220'], 'all_cls_res', 'all_cls_ei_inds', 'all_varCE', 'rng_vals')

end
%% Load AUC stats

% load these from the saved data files
auc_file = [local_res_dir network_param.network_tag filesep ...
    'matfiles' filesep 'selXCJ_stats_' ...
    network_param.network_tag stimulus_experiment_tag];

x_aucstats = load(auc_file);

% Load other dataset info
all_aucstats = cell(length(rng_vals), 1);

for i_rv = 1:length(rng_vals)
    new_rng_str = ['rng' num2str(rng_vals(i_rv))];
    rng_locs = regexp(auc_file, 'rng\w*', 'split');
    new_net_auc_file = [rng_locs{1} new_rng_str rng_locs{2} new_rng_str stimulus_experiment_tag];
    all_aucstats{i_rv} = load(new_net_auc_file);
end

%% load necessary variables from classifier results

t_vec = proc_simdata.peri_stim_time;
%%

stim_resps = activity_struct.stim_resps;
stim_y = activity_struct.stim_y;
classifier_options = x_class.classifier_options; 
c_timestep = activity_struct.c_timestep;
train_inds = x_class.cls_res{1}.train_inds;
test_inds = x_class.cls_res{1}.test_inds;
% skip_dt_units = x_class.skip_dt_units;

% Classifier info
all_units_index = 1;
all_units_label = x_class.cls_res{all_units_index}.cellGroup;
classIdataJ_valAcc = x_class.cls_res{all_units_index}.classIdataJ_valAcc;
%%
% input info
netwk_hz = 50;  % this is interpretation only - comes from setting overall network tau = 1 [ 20 ms ]
in1_freq = input_params.input_frequency{1};
in2_freq = input_params.input_frequency{2};


postL1_cells = activity_struct.postL1_cells;
inh_cells = activity_struct.inh_cells;

%% set color standards

stim_colors = lines(4);
stim_colors = stim_colors(3:4, :);
neuron_colors = [(1 + lines(1))/2; lines(2)];
samediff_colors = [0 0 0; 1 1 1]*0.5;

fsz_axes = 14;
fsz_labels = 16;


%% Figure 2: single-cell AUC, etc. 

% which networks have pairs with opposite selectivity? 
% if not, a sign that the input scaling didn't work - flag for follow up

has_opp_sel_cells = cellfun(@(x) mean(mean(x.pairs_oppSel(x.postL1_cells, x.postL1_cells))), ... 
    all_aucstats);

has_same_sel_cells = cellfun(@(x) mean(mean(x.pairs_sameSel(x.postL1_cells, x.postL1_cells))), ... 
    all_aucstats);

stim1_ave = arrayfun(@(x) mean(x.pop_ave_fr(x.stim_id == 1, :), 1), all_cls_ei_inds, 'UniformOutput', false);
stim2_ave = arrayfun(@(x) mean(x.pop_ave_fr(x.stim_id == 2, :), 1), all_cls_ei_inds, 'UniformOutput', false);

same_ave_fr = cellfun(@(x,y) ranksum(x,y), stim1_ave', stim2_ave');


keep_nets = has_opp_sel_cells > 0 & same_ave_fr > 1e-4;

x_nets_aucstats = all_aucstats(keep_nets);

%%
hFig = makeMyFigure(6.5*2.54*(11/7), 4.5*2.54*11/7);

x_aucstats = x_nets_aucstats{2};

x_lims = [1 length(x_nets_aucstats)];
x_ticks= 1:1:length(x_nets_aucstats);
auroc_bins = x_aucstats.auroc_bins;
exc_2tail = x_aucstats.exc_2tail;
% exc_auroc_pdf = x_aucstats.exc_auroc_pdf;
inh_2tail = x_aucstats.inh_2tail;
% inh_auroc_pdf = x_aucstats.inh_auroc_pdf;

% smooth auroc pdf with gaussian kernel with sigma = SE of auroc vals. 
auroc_smooth_sigma = mean(std(x_aucstats.cell_auroc_shuf, 0, 2));
exc_aurocs = x_aucstats.cell_auroc_val(~x_aucstats.inh_cells & x_aucstats.postL1_cells);
inh_aurocs = x_aucstats.cell_auroc_val(x_aucstats.inh_cells);
auc_kernel = @(auc, auc_bins) exp(-(auc_bins - auc).^2/(2*auroc_smooth_sigma^2));
auroc_bins = auroc_bins(1):0.01:auroc_bins(end);
exc_auroc_pdf2 = bsxfun(@(x,y) auc_kernel(x,y), exc_aurocs, auroc_bins);
exc_auroc_pdf = sum(exc_auroc_pdf2)/sqrt(2*pi*auroc_smooth_sigma^2)/length(exc_aurocs);

inh_auroc_pdf2 = bsxfun(@(x,y) auc_kernel(x,y), inh_aurocs, auroc_bins);
inh_auroc_pdf = sum(inh_auroc_pdf2)/sqrt(2*pi*auroc_smooth_sigma^2)/length(inh_aurocs);





dbin = x_aucstats.dbin;
probConnect_SaOp_labels = x_nets_aucstats{1}.probConnect_SaOp_labels;
% this traces out the portion of the pdf that extends beyond the [2.5,
% 97.5] range of the shuffled distribution for excitatory cells
exc_low_rng_SIGaur = unique([auroc_bins(1):dbin:exc_2tail(1)  exc_2tail(1)]);
exc_high_rng_SIGaur = unique([auroc_bins(end):-dbin:exc_2tail(2)  exc_2tail(2)]);
exc_sig_aurocLO_pdf = interp1(auroc_bins, exc_auroc_pdf, exc_low_rng_SIGaur);
exc_sig_aurocHI_pdf = interp1(auroc_bins, exc_auroc_pdf, exc_high_rng_SIGaur);
% this traces out the portion of the pdf that extends beyond the [2.5,
% 97.5] range of the shuffled distribution for inhibitory cells
inh_low_rng_SIGaur = unique([auroc_bins(1):dbin:inh_2tail(1) inh_2tail(1)]);
inh_high_rng_SIGaur = unique([auroc_bins(end):-dbin:inh_2tail(2) inh_2tail(2)]);
inh_sig_aurocLO_pdf = interp1(auroc_bins, inh_auroc_pdf, inh_low_rng_SIGaur);
inh_sig_aurocHI_pdf = interp1(auroc_bins, inh_auroc_pdf, inh_high_rng_SIGaur);

abc_w = 0.2;
de_w = 0.35;
y_pos2 = 0.1;
h_val = 0.38;
h_gap = 0.1;
y_pos1 = y_pos2 + h_val + h_gap;

gap_bcd = 0.08;
x_a = gap_bcd;
x_b = 2*gap_bcd + abc_w;
x_c = x_b + gap_bcd + abc_w;
x_d = x_a;
x_e = x_d + de_w + gap_bcd; %x_b + (abc_w + gap_bcd)*[0 0 1 1] ;
h_e = h_val/2*[1 2 1 2];
dh_e = 0.8*h_e(1);
hB = axes(hFig, 'Position', [x_b y_pos1 abc_w h_val]);
hC = axes(hFig, 'Position', [x_c y_pos1 abc_w h_val]);
hD = axes(hFig, 'Position', [x_d y_pos2 de_w h_val]);
hE = axes(hFig, 'Position', [x_e y_pos2 de_w h_val]);


hA = axes(hFig, 'Position', [gap_bcd y_pos1 abc_w h_val ]);



hold on
ph(2) = plot(auroc_bins, exc_auroc_pdf, 'Color', neuron_colors(2, :), 'linewidth', 1.5);

patch([exc_low_rng_SIGaur exc_low_rng_SIGaur([end 1])], [exc_sig_aurocLO_pdf 0 0], ...
    neuron_colors(2, :), 'facealpha', 0.4)
patch([exc_high_rng_SIGaur exc_high_rng_SIGaur([end 1])], [exc_sig_aurocHI_pdf 0 0], ...
    neuron_colors(2, :), 'facealpha', 0.4)



ph(1) = plot(auroc_bins, inh_auroc_pdf, 'Color', neuron_colors(3, :), 'linewidth', 1.5);

xlabel('AUC')
ylabel('probability density')
set(gca, 'color', 'none', 'Fontsize', fsz_axes);

text(0.225, 12, 'Excitatory', 'Color', neuron_colors(2, :), ...
    'FontSize', fsz_axes, 'HorizontalAlignment', 'left', ...
    'FontWeight', 'bold')
text(0.225, 11, 'Inhibitory', 'color', neuron_colors(3, :), ...
    'FontSize', fsz_axes, 'HorizontalAlignment', 'left', ...
    'FontWeight', 'bold')
% % % % Now add results from multiple network runs by cycling through
% x_nets_aucstats
saOp_fields = {'sameSel_etoi_J_vals', 'opp_Sel_etoi_J_vals', ...
    'sameSel_itoe_J_vals', 'opp_Sel_itoe_J_vals'};

% keep the "equivalent connectivity" as probbaiblity of connection between
% pools (ordered: [E1 E2 I1 I2]. Treat E1/E2 symmetrically.)
equiv_poolConnect = cellfun(@(x) zeros(4), x_nets_aucstats, 'UniformOutput', false);
ei_pool_order = {{[1 3], [2 4]}, {[1 4], [2 3]}, ... 
    {[3 1], [4 2]}, {[4 1], [3 2]}};


for i_net = 1:length(x_nets_aucstats)
    aucstat = x_nets_aucstats{i_net};
    
    set(hFig, 'currentAxes', hB);
    hold on
    errorbar(i_net-0.1, aucstat.frac_sel_EI(1), aucstat.sd_frac_sel_EI(1), 'o', ...
        'color', neuron_colors(2, :), 'linewidth', 1)
    errorbar(i_net+0.1, aucstat.frac_sel_EI(2), aucstat.sd_frac_sel_EI(2), 'o', ...
        'color', neuron_colors(3, :), 'linewidth', 1)
    xlabel('network instance #')
    ylabel('fraction of cells selective')
    xlim(x_lims + [-0.5 0.5])
    ylim([0 1])
    set(gca, 'color', 'none', 'xtick', x_ticks, 'Fontsize', fsz_axes);

    set(hFig, 'currentAxes', hC);
    hold on
    errorbar(i_net-0.1, aucstat.ave_choice_EI(1), aucstat.se_choice_EI(1),...
        'o', 'color', neuron_colors(2, :), 'linewidth', 1)
    errorbar(i_net+0.1, aucstat.ave_choice_EI(2), aucstat.se_choice_EI(2),...
        'o', 'color', neuron_colors(3, :), 'linewidth', 1)
    xlabel('network instance #')
    ylabel('choice selectivity')
    xlim(x_lims + [-0.5 0.5])
    ylim([0 0.2])
    set(gca, 'color', 'none', 'xtick',x_ticks, 'Fontsize', fsz_axes);

% 
    for ii = 1:4
%         set(hFig, 'currentAxes', hE(ii));
        
        j_ii = aucstat.(saOp_fields{ii});
        y_val = aucstat.probConnect_SaOp(ii);
        dy_val = sqrt(sum(j_ii~=0))/length(j_ii);
        % save the equiv_connectivity values
        for i_sym = 1:2
            eqC_i = ei_pool_order{ii}{i_sym}(2);
            eqC_j = ei_pool_order{ii}{i_sym}(1);
            equiv_poolConnect{i_net}(eqC_i, eqC_j) = y_val;
        end
    end
        
end

% make panel labels
lABC_y = 6.83;
lABC_x = -2;
set(hFig, 'currentAxes', hA);
text(lABC_x, lABC_y, 'A', 'units', 'centimeters', 'FontSize', fsz_labels, 'FontWeight', 'bold')

set(hFig, 'currentAxes', hB);
text(lABC_x, lABC_y, 'B', 'units', 'centimeters', 'FontSize', fsz_labels, 'FontWeight', 'bold')

set(hFig, 'currentAxes', hC);
text(lABC_x, lABC_y, 'C', 'units', 'centimeters', 'FontSize', fsz_labels, 'FontWeight', 'bold')


keep_cls_res = all_cls_res(keep_nets);
keep_cls_einds = all_cls_ei_inds(keep_nets);
num_nets = length(keep_cls_res);
ex_net = 1;
i_TC = 50;
ex_cIdJ_cell = keep_cls_res{ex_net};
t_vec = 1:size(ex_cIdJ_cell{1}.classIdataJ_valAcc, 2);
plot_groups = [1 7 8];
group_colors = [0 0 0; neuron_colors([2 3], :)];

wt_lims = 0.15*[-1 1];


lh = 0*plot_groups;
for ii = 1:length(plot_groups)
    set(gcf, 'CurrentAxes', hD)
    hold on
    ex_cIdJ_mat = ex_cIdJ_cell{plot_groups(ii)}.classIdataJ_valAcc;
    lh(ii) = plot(t_vec, ex_cIdJ_mat(i_TC, :), 'color', group_colors(ii, :), ...
        'linewidth', 1);
    plot(t_vec(i_TC), ex_cIdJ_mat(i_TC, i_TC), 'o', 'color', group_colors(ii, :));
    
    leg_ent{ii} = ex_cIdJ_cell{plot_groups(ii)}.cellGroup;
    

    disp(['Ave decoding, ' leg_ent{ii} ' last half: ' num2str(mean(ex_cIdJ_mat(i_TC, end/2:end)), '%1.3f')])
end


xlabel('time at which classifier was tested (\tau)')
ylabel('classifier accuracy')
ylim([0.5 1])

set(gcf, 'CurrentAxes',  hD)
hold on

text(-0.2, 1, 'D', 'units', 'normalized', 'FontSize', fsz_labels, 'FontWeight', 'bold')
set(hD, 'color', 'none', 'fontsize', fsz_axes, ...
    'box', 'off', 'ytick', 0.6:.2:1, 'xtick', 10:10:50)
pushTicksOut(hD);
% legend(lh, leg_ent, 'location', 'southeast')
text(2, 1.04, 'All', 'Color', [0 0 0], ...
    'FontSize', fsz_axes, 'HorizontalAlignment', 'left', ...
    'FontWeight', 'bold')
text(2, 0.97, 'Excitatory', 'Color', neuron_colors(2, :), ...
    'FontSize', fsz_axes, 'HorizontalAlignment', 'left', ...
    'FontWeight', 'bold')
text(2, 0.92, 'Inhibitory', 'color', neuron_colors(3, :), ...
    'FontSize', fsz_axes, 'HorizontalAlignment', 'left', ...
    'FontWeight', 'bold')

set(gcf, 'CurrentAxes', hE)
text(-0.2, 1, 'E', 'units', 'normalized', 'FontSize', fsz_labels, 'FontWeight', 'bold')

x_net = 1:length(keep_cls_res);
hold on
for ii = 1:length(plot_groups)
    valAcc_ais = cellfun(@(x) x{plot_groups(ii)}.classIdataJ_valAcc(i_TC, i_TC), ...
        keep_cls_res);
    sd_valAcc = cellfun(@(x) std(x{plot_groups(ii)}.cvDraws_testAcc), ...
        keep_cls_res);
    mean_valAcc = cellfun(@(x) mean(x{plot_groups(ii)}.cvDraws_testAcc), ...
        keep_cls_res);
    ph = errorbar(x_net, valAcc_ais, sd_valAcc, 'o', ...
        'linewidth', 1, 'color', group_colors(ii, :));
    disp([num2str(plot_groups(ii)) ' val Acc: ' num2str(mean_valAcc')])
end
ylim([0.5 1.0])
xlim(x_net([1 end]) + [-0.5 0.5])
xlabel('network sim #')
ylabel('classifier accuracy')
set(gca, 'color', 'none', 'fontsize', fsz_axes, 'box', 'off')
pushTicksOut(hE);

print(gcf, '-dpdf', 'writeup/figures/fig4_final')


%%

%% New plots (12/18/19) : show the differences in connectivity between E1-E2
% and I1-I2 

hFig2_S = makeMyFigure(6.5*2.54*(11/7), 4.5*2.54*11/7);

gap_bcd = 0.08;
x_a = gap_bcd;
x_b = 2*gap_bcd + abc_w;
x_c = x_b + gap_bcd + abc_w;
x_d = x_a;
x_e = x_a + (abc_w + gap_bcd)*[0 0 1 1] ;
h_e = h_val/2*[1 2 1 2];
dh_e = 0.8*h_e(1);

for ii = 1:4
    hE(ii) = axes(hFig2_S, 'Position', [x_e(ii) 1 - dh_e - h_e(ii) abc_w dh_e]);
%     title(num2str(ii))
end


labels = {'sameSel_etoe', 'oppSel_etoe', 'sameSel_itoi', 'oppSel_itoi'};

% keep the "equivalent connectivity" as probability of connection between
% pools (ordered: [E1 E2 I1 I2]. Treat E1/E2 symmetrically.) It is *very
% important* that these match the labels
eeii_pool_order = {{[1 1], [2 2]}, {[1 2], [2 1]}, ... 
    {[3 3], [4 4]}, {[4 3], [3 4]}};


for i_net = 1:length(x_nets_aucstats)
    aucstat = x_nets_aucstats{i_net};
    
    is_exc2 = aucstat.postL1_cells & ~aucstat.inh_cells;
    is_inh = aucstat.inh_cells;
    net_J = aucstat.network_param.J;
    % this has 1's in the E-E connections and I-I connections
    poolJ_logical = {is_exc2*is_exc2', is_exc2*is_exc2', ... 
        is_inh*is_inh', is_inh*is_inh'};
    poolJ_logical = cellfun(@(x) logical(x), poolJ_logical, 'UniformOutput', false);
    %%
   for ii = 1:4
        set(hFig2_S, 'currentAxes', hE(ii));
         
        if ii == 1 || ii == 3
                % take pairs with same selectivity
            sameOpp_logical = aucstat.pairs_sameSel;
        else
            sameOpp_logical = aucstat.pairs_oppSel;
        end
        % compute probability of connection and the SE
        poolJ_bySel = net_J(poolJ_logical{ii} & sameOpp_logical);
        y_val = mean(poolJ_bySel ~= 0);
        dy_val = sqrt(sum(poolJ_bySel~=0))/length(poolJ_bySel);
        
        % save the equiv_connectivity values
        for i_sym = 1:2
            eqC_i = eeii_pool_order{ii}{i_sym}(2);
            eqC_j = eeii_pool_order{ii}{i_sym}(1);
            equiv_poolConnect{i_net}(eqC_i, eqC_j) = y_val;
        end
        
        hold on 
        errorbar(i_net, y_val,dy_val, 'o', ...
            'Color', samediff_colors(mod(ii-1, 2)+1, :), 'linewidth', 1)
        if i_net == 1
            plot(x_lims + [-0.5 0.5], network_param.sparsity*[1 1], 'k-.')
            if ii == 4 || ii == 2
                xlabel('network instance #')
                if ii == 2
                    h_y = ylabel('Probability of connection', 'Position', [-2.15 0.4 -1]);
                end
                set(gca, 'xtick', x_ticks)
            else
                set(gca, 'xtick', []);
            end
            set(gca, 'color', 'none', 'Fontsize', fsz_axes)
            xlim(x_lims + [-0.5 0.5])
        end
        if i_net == length(x_nets_aucstats)
            % comment this line out, but uncomment to check labels are
            % correct
%             title(probConnect_SaOp_labels{ii})
            if ii == 1
                text(2, 0.35, 'E', 'Color', neuron_colors(2, :), ...
                    'FontSize', fsz_axes, 'HorizontalAlignment', 'center', ...
                    'FontWeight', 'bold')
                text(5, 0.35, 'E', 'color', neuron_colors(2, :), ...
                    'FontSize', fsz_axes, 'HorizontalAlignment', 'center', ...
                    'FontWeight', 'bold')
                h_arrow1 = annotation( 'arrow', 0.1553 + [0 0.04], 0.8097*[1 1]);
            elseif ii == 3
                text(5, 0.35, 'I', 'Color', neuron_colors(3, :), ...
                    'FontSize', fsz_axes, 'HorizontalAlignment', 'center', ...
                    'FontWeight', 'bold')
                text(2, 0.35, 'I', 'color', neuron_colors(3, :), ...
                    'FontSize', fsz_axes, 'HorizontalAlignment', 'center', ...
                    'FontWeight', 'bold')
                h_arrow2 = annotation( 'arrow', 0.4333+[0 0.04], 0.8097*[1 1]);
                
                
                text(length(x_nets_aucstats)*1.2, 0.2, {'same';'selectivity'}, ...
                    'Color', samediff_colors(1, :), 'HorizontalAlignment', 'center', ...
                    'FontSize', fsz_axes, 'Rotation', -90)
            elseif ii == 4
                h_text = text(length(x_nets_aucstats)*1.2, 0.2, {'opposite';'selectivity'}, ...
                    'Color', samediff_colors(2, :), 'HorizontalAlignment', 'center', ...
                    'FontSize', fsz_axes, 'Rotation', -90);
            end
            ylim([0.05 0.35])
        end
   end
end

hAL = text(-9.1, 7.1, 'A', 'units', 'centimeters', 'FontSize', fsz_labels, 'FontWeight', 'bold');

hBL =  text(7.3, 7.1, 'B', 'units', 'centimeters', 'FontSize', fsz_labels, 'FontWeight', 'bold');


print(gcf, '-dpdf', ['writeup/figures/fig2_connectivity_sfig_' connect_tag]);
%% print data for Table S1
net_sparsity = cellfun(@(x) mean(x.network_param.J(:) ~= 0), all_aucstats)
no_nan = cellfun(@(x) ~any(isnan(x(:))), equiv_poolConnect);

ave_equiv_connect = reshape(mean(cell2mat(cellfun(@(x) x(:)', ...
    equiv_poolConnect(no_nan), 'UniformOutput', false)), 1), [4 4]);
se_equiv_connect = reshape(std(cell2mat(cellfun(@(x) x(:)', ...
    equiv_poolConnect(no_nan), 'UniformOutput', false)), 1), [4 4]);

net_sparsityE = cellfun(@(x) mean(mean(x.network_param.J(:, 1:400) ~= 0)), all_aucstats)
net_sparsityI = cellfun(@(x) mean(mean(x.network_param.J(:, 1+400:end) ~= 0)), all_aucstats)

norm_equiv_connect = [ave_equiv_connect(:, 1:2)/mean(net_sparsityE) ave_equiv_connect(:, 3:4)/mean(net_sparsityI) ]
se_norm_equiv_connect = [se_equiv_connect(:, 1:2)/mean(net_sparsityE) se_equiv_connect(:, 3:4)/mean(net_sparsityI) ]

tex_str_XY = @(s1,s2) ['$' s2 '\rightarrow ' s1 ' $ & & $'];

disp([tex_str_XY('E_1','E_1')  num2str(norm_equiv_connect(1, 1), '%1.2f') ' \pm ' num2str(se_norm_equiv_connect(1, 1), '%1.2f') '$ \\'])
disp([tex_str_XY('E_1','E_2')  num2str(norm_equiv_connect(1, 2), '%1.2f') ' \pm ' num2str(se_norm_equiv_connect(1, 2), '%1.2f') '$ \\'])
disp([tex_str_XY('E_1','I_1')  num2str(norm_equiv_connect(1, 3), '%1.2f') ' \pm ' num2str(se_norm_equiv_connect(1, 3), '%1.2f') '$ \\'])
disp([tex_str_XY('E_1','I_2')  num2str(norm_equiv_connect(1, 4), '%1.2f') ' \pm ' num2str(se_norm_equiv_connect(1, 4), '%1.2f') '$ \\'])
disp([tex_str_XY('I_1','E_1')  num2str(norm_equiv_connect(3, 1), '%1.2f') ' \pm ' num2str(se_norm_equiv_connect(3, 1), '%1.2f') '$ \\'])
disp([tex_str_XY('I_1','E_2')  num2str(norm_equiv_connect(3, 2), '%1.2f') ' \pm ' num2str(se_norm_equiv_connect(3, 2), '%1.2f') '$ \\'])
disp([tex_str_XY('I_1','I_1')  num2str(norm_equiv_connect(3, 3), '%1.2f') ' \pm ' num2str(se_norm_equiv_connect(3, 3), '%1.2f') '$ \\'])
disp([tex_str_XY('I_1','I_2')  num2str(norm_equiv_connect(3, 4), '%1.2f') ' \pm ' num2str(se_norm_equiv_connect(3, 4), '%1.2f') '$ \\'])
%%

%% connectivity strengths 
disp(['E->E: ' num2str(mean(nonzeros(x_aucstats.network_param.J(1:400, 1:400))))])
disp(['E->I: ' num2str(mean(nonzeros(x_aucstats.network_param.J(1+400:end, 1:400))))])
disp(['I->E: ' num2str(mean(nonzeros(x_aucstats.network_param.J(1:400, 401:end))))])
disp(['I->I: ' num2str(mean(nonzeros(x_aucstats.network_param.J(401:end, 401:end))))])

disp(['SE E->E: ' num2str(std(nonzeros(x_aucstats.network_param.J(1:400, 1:400))))])
disp(['SE E->I: ' num2str(std(nonzeros(x_aucstats.network_param.J(1+400:end, 1:400))))])
disp(['SE I->E: ' num2str(std(nonzeros(x_aucstats.network_param.J(1:400, 401:end))))])
disp(['SE I->I: ' num2str(std(nonzeros(x_aucstats.network_param.J(401:end, 401:end))))])

%% prob of connections

prob_connect_all = cell2mat(cellfun(@(x) x.probConnect_SaOp, all_aucstats, 'uniformoutput', false));

mean_p_connect = mean(prob_connect_all, 1);
std_p_connect = std(prob_connect_all,[], 1);

labels = x_aucstats.probConnect_SaOp_labels;

disp(labels)
disp(num2str(mean_p_connect, '%1.2f '))
disp(num2str(std_p_connect, '%1.2f '))

%% fraction selective
all_frac_sel = cell2mat(cellfun(@(x) x.frac_sel_EI, all_aucstats, 'Uniformoutput', false));

all_abs_sel = cell2mat(cellfun(@(x) x.ave_choice_EI, all_aucstats, 'uniformoutput', false));

disp(['Min frac sel: ' num2str(min(all_frac_sel, [], 1), '%1.2f ')])
disp(['Max frac sel: ' num2str(max(all_frac_sel, [], 1), '%1.2f ')])
disp(['Mean frac sel: ' num2str(mean(all_frac_sel, 1), '%1.2f ')])
disp(['StdDev frac sel: ' num2str(std(all_frac_sel,[], 1), '%1.2f ')])

disp(['Min abs sel: ' num2str(min(all_abs_sel, [], 1), '%1.2f ')])
disp(['Max abs sel: ' num2str(max(all_abs_sel, [], 1), '%1.2f ')])
disp(['Mean abs sel: ' num2str(mean(all_abs_sel, 1), '%1.2f ')])
disp(['StdDev abs sel: ' num2str(std(all_abs_sel, [], 1), '%1.2f ')])

