% Make Figures 1 through 3 

%%%% set parameters [ must run simulations first]
num_neur = 500;
gain_neur = 400; % usually 400; for spnormE, set to 200. 
rng_val = 1; % usualy 1
ei_ratio = 100; % 92 or 108;  80 or 125
pulse_width = .50; % 50 or 150 or 500
connect_tag = 'spnorm'; %  
f2f_type = 'f2f';

file_tag = [connect_tag '_g' num2str(gain_neur) '_n' num2str(num_neur) '_ei' num2str(ei_ratio) '_rng' num2str(rng_val)];

local_res_dir = '/Volumes/Samsung_T5/projects/random_rec_network_project/results/'; % classifier results are here

input_type = [f2f_type 'discExp_pw' num2str(100*pulse_width)];
%% required for Fig 1
stimulus_experiment_tag = [f2f_type 'discExp_pw' num2str(100*pulse_width)];
x_res = load([local_res_dir file_tag '/matfiles/sim_data_' file_tag '_' stimulus_experiment_tag '.mat']);

%%
proc_simdata = x_res.proc_simdata;
% %% short load
% x_res = load([local_res_dir file_tag '/matfiles/sim_data_' file_tag '_' ...
%     input_type '.mat'], 'input_params', ...
%     'network_param');
%%
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

if ~exist(['multiRNG_' file_tag '_20190903.mat'], 'file')
rng_vals = [1 2 3 4 5 6 7 8 9 10 11 12 13 17] ; %[1 2 3 4 5 6 7 8 9 10];
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
save(['multiRNG_' file_tag '_20190903'], 'all_cls_res', 'all_cls_ei_inds', 'all_varCE', 'rng_vals')
%%
else
%     load(['multiRNG_' file_tag], 'all_cls_res_noVA');
    load(['multiRNG_' file_tag '_20190903'], 'all_cls_res', 'all_cls_ei_inds', 'all_varCE', 'rng_vals')

end
%% Load AUC stats

% load these from the saved data files : REPLACE WITH PATH TO THE RESULTS
% ON THE DRIVE. 
auc_file = [local_res_dir network_param.network_tag filesep ...
    'matfiles' filesep 'selXCJ_stats_' ...
    network_param.network_tag stimulus_experiment_tag];

x_aucstats = load(auc_file);

% Load other dataset info
% rng_vals = [1:10];    % rng_vals should be set in the first block of this
% script when the classifier results are loaded. 
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
% allClassifier = x_class.full_cls(all_units_index).classifier;
% allValidationAccuracy = x_class.cls_res(all_units_index).validationAccuracy;
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
%%
%%
% For proc_simdata, use stimulus_info_array to get stimulus label (will
% have 1 in the (:, 1, 1) spot)
stim_is1 = squeeze(proc_simdata.stimulus_info_array(:, 1, 1) == 1);
integrationBin = proc_simdata.peri_stim_time > 5 & proc_simdata.peri_stim_time <= 50;
popave_L2_fr_S1 = squeeze(mean(mean(proc_simdata.dr_out_peristim(stim_is1 == 1, :, integrationBin), 3), 2)); 
popave_L2_fr_S2 = squeeze(mean(mean(proc_simdata.dr_out_peristim(stim_is1 ~= 1, :, integrationBin), 3), 2)); 

%%

% Most of Figure 1 is a cartoon; Panel E shows network activity
% distribution (averaged over units). 
hFig1E = makeMyFigure(6, 6);
hold on
histogram(popave_L2_fr_S1, 'FaceColor', stim_colors(1, :))
histogram(popave_L2_fr_S2, 'FaceColor', stim_colors(2, :))
xlabel('average network activity')
ylabel('trial counts')
axis square
set(gca, 'color', 'none')
legend({num2str(netwk_hz*in1_freq, 'freq : %1.0f'),...
    num2str(netwk_hz*in2_freq, 'freq : %1.0f')}, 'Location', 'northwest')


print(gcf, '-dpdf', 'writeup/figures/fig1_final');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

x_aucstats = x_nets_aucstats{1};

x_lims = [1 length(x_nets_aucstats)];
x_ticks= 1:2:length(x_nets_aucstats);
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

abcd_w = 0.2;
y_pos2 = 0.1;
h_val = 0.38;
h_gap = 0.1;
y_pos1 = y_pos2 + h_val + h_gap;

gap_bcd = 0.08;
x_a = gap_bcd;
x_b = 2*gap_bcd + abcd_w;
x_c = x_b + gap_bcd + abcd_w;
x_d = x_a;
x_e = x_b + (abcd_w + gap_bcd)*[0 0 1 1] ;
h_e = h_val/2*[1 2 1 2];
dh_e = 0.8*h_e(1);
hB = axes(hFig, 'Position', [x_b y_pos1 abcd_w h_val]);
hC = axes(hFig, 'Position', [x_c y_pos1 abcd_w h_val]);
hD = axes(hFig, 'Position', [x_d y_pos2 abcd_w h_val]);

for ii = 1:4
    hE(ii) = axes(hFig, 'Position', [x_e(ii) 1 - y_pos2 - h_val - h_e(ii) abcd_w dh_e]);
%     title(num2str(ii))
end

hA = axes(hFig, 'Position', [gap_bcd y_pos1 abcd_w h_val ]);



hold on
ph(2) = plot(auroc_bins, exc_auroc_pdf, 'Color', neuron_colors(2, :), 'linewidth', 1.5);
% plot([exc_low_rng_SIGaur exc_low_rng_SIGaur([end 1])], [exc_sig_aurocLO_pdf 0 0], 'Color', neuron_colors(2, :))
% plot([exc_high_rng_SIGaur exc_high_rng_SIGaur([end 1])], [exc_sig_aurocHI_pdf 0 0], 'Color', neuron_colors(2, :))
patch([exc_low_rng_SIGaur exc_low_rng_SIGaur([end 1])], [exc_sig_aurocLO_pdf 0 0], ...
    neuron_colors(2, :), 'facealpha', 0.4)
patch([exc_high_rng_SIGaur exc_high_rng_SIGaur([end 1])], [exc_sig_aurocHI_pdf 0 0], ...
    neuron_colors(2, :), 'facealpha', 0.4)


% exc_high_vertices = [exc_high_rng_SIGaur exc_high_rng_SIGaur([end 1]); ...
%     exc_sig_aurocHI_pdf 0 0]';
% exc_low_vertices = [exc_low_rng_SIGaur exc_low_rng_SIGaur([end 1]); ...
%     exc_sig_aurocLO_pdf 0 0]';
% h_excH = drawpolygon('position', exc_high_vertices, 'Color', neuron_colors(2, :), ...
%     'HandleVisibility', 'off');
% h_excL = drawpolygon('position', exc_low_vertices, 'Color', neuron_colors(2, :));


ph(1) = plot(auroc_bins, inh_auroc_pdf, 'Color', neuron_colors(3, :), 'linewidth', 1.5);
% plot([inh_low_rng_SIGaur inh_low_rng_SIGaur([end 1])], [inh_sig_aurocLO_pdf 0 0], 'Color', neuron_colors(3, :))
% plot([inh_high_rng_SIGaur inh_high_rng_SIGaur([end 1])], [inh_sig_aurocHI_pdf 0 0], 'Color', neuron_colors(3, :))
xlabel('AUC')
ylabel('probability density')
set(gca, 'color', 'none', 'Fontsize', fsz_axes);
% legend(ph, {'inhibitory', 'excitatory'}, 'Position', ...
%     [0.20 0.90 0.11 0.052], 'box', 'off', 'color','none')
text(0.225, 12, 'Excitatory', 'Color', neuron_colors(2, :), ...
    'FontSize', fsz_axes, 'HorizontalAlignment', 'left', ...
    'FontWeight', 'bold')
text(0.225, 11, 'Inhibitory', 'color', neuron_colors(3, :), ...
    'FontSize', fsz_axes, 'HorizontalAlignment', 'left', ...
    'FontWeight', 'bold')
% % % % Now add results from multiple network runs by cycling through
% x_netw_aucstats
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
    % NOTE: 'ave_choic_EI = |auc - 0.5| from AJS code; in Najafi paper, it
    % is 2*|auc - 0.5|, so that has been corrected here to match. 
    errorbar(i_net-0.1, 2*aucstat.ave_choice_EI(1), 2*aucstat.se_choice_EI(1),...
        'o', 'color', neuron_colors(2, :), 'linewidth', 1)
    errorbar(i_net+0.1, 2*aucstat.ave_choice_EI(2), 2*aucstat.se_choice_EI(2),...
        'o', 'color', neuron_colors(3, :), 'linewidth', 1)
    xlabel('network instance #')
    ylabel('choice selectivity')
    xlim(x_lims + [-0.5 0.5])
    ylim([0 0.2])
    set(gca, 'color', 'none', 'xtick',x_ticks, 'Fontsize', fsz_axes);

    set(hFig, 'currentAxes', hD);
    hold on
    ph(1) = errorbar(i_net-0.1, aucstat.mean_xc_SaOp_etoi(1), aucstat.sem_xc_SaOp_etoi(1),...
        'ko', 'linewidth', 1);
    ph(2) = errorbar(i_net+0.1, aucstat.mean_xc_SaOp_etoi(2), aucstat.sem_xc_SaOp_etoi(2),...
        'o', 'color', 0.5*[1 1 1], 'linewidth', 1);
    xlabel('network instance #')
    ylabel('corr. coef')

    xlim(x_lims + [-0.5 0.5])
    ylim([-1 1])
    set(gca, 'color', 'none', 'xtick', x_ticks, 'Fontsize', fsz_axes);
%     if i_net == length(x_nets_aucstats)
%         legend(ph, {'same Sel', 'opp. Sel'}, 'location', 'best')
%     end

    for ii = 1:4
        set(hFig, 'currentAxes', hE(ii));
        
        j_ii = aucstat.(saOp_fields{ii});
        y_val = aucstat.probConnect_SaOp(ii);
        % save the equiv_connectivity values
        for i_sym = 1:2
            eqC_i = ei_pool_order{ii}{i_sym}(2);
            eqC_j = ei_pool_order{ii}{i_sym}(1);
            equiv_poolConnect{i_net}(eqC_i, eqC_j) = y_val;
        end
        dy_val = sqrt(sum(j_ii~=0))/length(j_ii);
        hold on 
        errorbar(i_net, y_val,dy_val, 'o', ...
            'Color', samediff_colors(mod(ii-1, 2)+1, :), 'linewidth', 1)
        if i_net == 1
            plot(x_lims + [-0.5 0.5], network_param.sparsity*[1 1], 'k-.')
            if ii == 4 || ii == 2
                xlabel('network instance #')
                if ii == 2
                    h_y = ylabel('probability of connection', 'Position', [-2.15 0.4 -1]);
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
                text(5, 0.35, 'E', 'Color', neuron_colors(2, :), ...
                    'FontSize', fsz_axes, 'HorizontalAlignment', 'center', ...
                    'FontWeight', 'bold')
                text(9, 0.35, 'I', 'color', neuron_colors(3, :), ...
                    'FontSize', fsz_axes, 'HorizontalAlignment', 'center', ...
                    'FontWeight', 'bold')
                h_arrow1 = annotation( 'arrow', 0.434 + [0 0.04], 0.4816*[1 1]);
            elseif ii == 3
                text(9, 0.35, 'E', 'Color', neuron_colors(2, :), ...
                    'FontSize', fsz_axes, 'HorizontalAlignment', 'center', ...
                    'FontWeight', 'bold')
                text(5, 0.35, 'I', 'color', neuron_colors(3, :), ...
                    'FontSize', fsz_axes, 'HorizontalAlignment', 'center', ...
                    'FontWeight', 'bold')
                h_arrow2 = annotation( 'arrow', 0.712+[0 0.04], 0.4816*[1 1]);
                
                
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

% make panel labels
lABC_y = 6.83;
lABC_x = -2;
set(hFig, 'currentAxes', hA);
text(lABC_x, lABC_y, 'A', 'units', 'centimeters', 'FontSize', fsz_labels, 'FontWeight', 'bold')

set(hFig, 'currentAxes', hB);
text(lABC_x, lABC_y, 'B', 'units', 'centimeters', 'FontSize', fsz_labels, 'FontWeight', 'bold')

set(hFig, 'currentAxes', hC);
text(lABC_x, lABC_y, 'C', 'units', 'centimeters', 'FontSize', fsz_labels, 'FontWeight', 'bold')


set(hFig, 'currentAxes', hD);
text(lABC_x, lABC_y, 'D', 'units', 'centimeters', 'FontSize', fsz_labels, 'FontWeight', 'bold')

text(1, -0.5, 'Same selectivity', 'Color', samediff_colors(1, :), ...
    'FontSize', fsz_axes, 'HorizontalAlignment', 'left', ...
    'FontWeight', 'bold')
text(1, -0.7, 'Opposite selectivity', 'color', samediff_colors(2, :), ...
    'FontSize', fsz_axes, 'HorizontalAlignment', 'left', ...
    'FontWeight', 'bold')


set(hFig, 'currentAxes', hE(1));
text(lABC_x, 2.66, 'E', 'units', 'centimeters', 'FontSize', fsz_labels, 'FontWeight', 'bold')

print(gcf, '-dpdf', 'writeup/figures/fig2_final')


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

all_abs_sel = cell2mat(cellfun(@(x) 2*x.ave_choice_EI, all_aucstats, 'uniformoutput', false));

disp(['Min frac sel: ' num2str(min(all_frac_sel, [], 1), '%1.2f ')])
disp(['Max frac sel: ' num2str(max(all_frac_sel, [], 1), '%1.2f ')])
disp(['Mean frac sel: ' num2str(mean(all_frac_sel, 1), '%1.2f ')])
disp(['StdDev frac sel: ' num2str(std(all_frac_sel,[], 1), '%1.2f ')])

disp(['Min abs sel: ' num2str(min(all_abs_sel, [], 1), '%1.2f ')])
disp(['Max abs sel: ' num2str(max(all_abs_sel, [], 1), '%1.2f ')])
disp(['Mean abs sel: ' num2str(mean(all_abs_sel, 1), '%1.2f ')])
disp(['StdDev abs sel: ' num2str(std(all_abs_sel, [], 1), '%1.2f ')])
%% Fig. S2: show the differences in connectivity between E1-E2
% and I1-I2, also compute the 'effective connectivity' in Table S1. 

hFig2_S = makeMyFigure(6.5*2.54*(11/7), 4.5*2.54*11/7);

gap_bcd = 0.08;
x_a = gap_bcd;
x_b = 2*gap_bcd + abcd_w;
x_c = x_b + gap_bcd + abcd_w;
x_d = x_a;
x_e = x_a + (abcd_w + gap_bcd)*[0 0 1 1] ;
h_e = h_val/2*[1 2 1 2];
dh_e = 0.8*h_e(1);
% hB = axes(hFig2_S, 'Position', [x_b y_pos1 abcd_w h_val]);
% hC = axes(hFig2_S, 'Position', [x_c y_pos1 abcd_w h_val]);
% hD = axes(hFig2_S, 'Position', [x_d y_pos2 abcd_w h_val]);

for ii = 1:4
    hE(ii) = axes(hFig2_S, 'Position', [x_e(ii) 1 - dh_e - h_e(ii) abcd_w dh_e]);
%     title(num2str(ii))
end


labels = {'sameSel_etoe', 'oppSel_etoe', 'sameSel_itoi', 'oppSel_itoi'};

% keep the "equivalent connectivity" as probbaiblity of connection between
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
        
%         j_ii = aucstat.(saOp_fields{ii});
%         y_val = aucstat.probConnect_SaOp(ii);
%         dy_val = sqrt(sum(j_ii~=0))/length(j_ii);
%         
%         
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
                text(5, 0.35, 'E', 'Color', neuron_colors(2, :), ...
                    'FontSize', fsz_axes, 'HorizontalAlignment', 'center', ...
                    'FontWeight', 'bold')
                text(9, 0.35, 'E', 'color', neuron_colors(2, :), ...
                    'FontSize', fsz_axes, 'HorizontalAlignment', 'center', ...
                    'FontWeight', 'bold')
                h_arrow1 = annotation( 'arrow', 0.1553 + [0 0.04], 0.8097*[1 1]);
            elseif ii == 3
                text(9, 0.35, 'I', 'Color', neuron_colors(3, :), ...
                    'FontSize', fsz_axes, 'HorizontalAlignment', 'center', ...
                    'FontWeight', 'bold')
                text(5, 0.35, 'I', 'color', neuron_colors(3, :), ...
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


print(gcf, '-dpdf', 'writeup/figures/sfig2_connectivity');
%%
hB = subplot(1, 3, 3);
cla;
no_nan = cellfun(@(x) ~any(isnan(x(:))), equiv_poolConnect);

ave_equiv_connect = reshape(mean(cell2mat(cellfun(@(x) x(:)', ...
    equiv_poolConnect(no_nan), 'UniformOutput', false)), 1), [4 4]);
se_equiv_connect = reshape(std(cell2mat(cellfun(@(x) x(:)', ...
    equiv_poolConnect(no_nan), 'UniformOutput', false)), 1), [4 4]);

net_sparsityE = cellfun(@(x) mean(mean(x.network_param.J(:, 1:400) ~= 0)), all_aucstats)
net_sparsityI = cellfun(@(x) mean(mean(x.network_param.J(:, 1+400:end) ~= 0)), all_aucstats)

norm_equiv_connect = [ave_equiv_connect(:, 1:2)/mean(net_sparsityE) ave_equiv_connect(:, 3:4)/mean(net_sparsityI) ]
se_norm_equiv_connect = [se_equiv_connect(:, 1:2)/mean(net_sparsityE) se_equiv_connect(:, 3:4)/mean(net_sparsityI) ]

tex_str_XY = @(s1,s2) ['$' s2 '\rightarrow ' s1 ' $ & $'];

% print out strings for making table in tex
disp([tex_str_XY('E_1','E_1')  num2str(norm_equiv_connect(1, 1), '%1.2f') ' \pm ' num2str(se_norm_equiv_connect(1, 1), '%1.2f') '$ '])
disp([tex_str_XY('E_1','E_2')  num2str(norm_equiv_connect(1, 2), '%1.2f') ' \pm ' num2str(se_norm_equiv_connect(1, 2), '%1.2f') '$ '])
disp([tex_str_XY('E_1','I_1')  num2str(norm_equiv_connect(1, 3), '%1.2f') ' \pm ' num2str(se_norm_equiv_connect(1, 3), '%1.2f') '$ '])
disp([tex_str_XY('E_1','I_2')  num2str(norm_equiv_connect(1, 4), '%1.2f') ' \pm ' num2str(se_norm_equiv_connect(1, 4), '%1.2f') '$ '])
disp([tex_str_XY('I_1','E_1')  num2str(norm_equiv_connect(3, 1), '%1.2f') ' \pm ' num2str(se_norm_equiv_connect(3, 1), '%1.2f') '$ '])
disp([tex_str_XY('I_1','E_2')  num2str(norm_equiv_connect(3, 2), '%1.2f') ' \pm ' num2str(se_norm_equiv_connect(3, 2), '%1.2f') '$ '])
disp([tex_str_XY('I_1','I_1')  num2str(norm_equiv_connect(3, 3), '%1.2f') ' \pm ' num2str(se_norm_equiv_connect(3, 3), '%1.2f') '$ '])
disp([tex_str_XY('I_1','I_2')  num2str(norm_equiv_connect(3, 4), '%1.2f') ' \pm ' num2str(se_norm_equiv_connect(3, 4), '%1.2f') '$ '])
% connectivity strengths
disp(['E->E: ' num2str(mean(nonzeros(x_aucstats.network_param.J(1:400, 1:400))))])
disp(['E->I: ' num2str(mean(nonzeros(x_aucstats.network_param.J(1+400:end, 1:400))))])
disp(['I->E: ' num2str(mean(nonzeros(x_aucstats.network_param.J(1:400, 401:end))))])
disp(['I->I: ' num2str(mean(nonzeros(x_aucstats.network_param.J(401:end, 401:end))))])


% %%
% line_weights = 1 + 3*log2(ave_equiv_connect);
% node_locs = [0.2 0.8; 0.8 0.8; 0.1 0.2; 0.9 0.2];
% plot(node_locs(:, 1), node_locs(:, 2), 'o', 'markersize', 20)
% hold on
% dxy = 0.1*[1 1];
% cmap = jet(16);
% for i_edge = 1:4
%     for j_edge = 1:4
%         if i_edge ~= j_edge
%             if i_edge < j_edge
%                 mdxy = -1;
%             else
%                 mdxy = 1;
%             end
%             x1y1 = node_locs(j_edge, :);
%             x2y2 = node_locs(i_edge, :);
%             xmym = (x1y1 + x2y2)/2 + mdxy*dxy;
%             x = [x1y1(1) ; xmym(1); x2y2(1) ]
%             y = [x1y1(2) ; xmym(2); x2y2(2) ]
%             xq = x(1):.001:x(3);
%             s = spline(x,y,xq)
%             plot(xq, s, 'color', cmap(i_edge*4 + j_edge - 4, :))
%             plot(x, y, 'o', 'color', cmap(i_edge*4 + j_edge - 4, :))
%         end
%     end
% end
% 
% 
% 
% axis([0 1 0 1])
% axis square
% set(gca, 'visible', 'off')

% %% plot equiv circuits
% bg = biograph(equiv_poolConnect{1}', {'E1', 'E2', 'I1', 'I2'})
% % draw edges according to weight
% for i_bg = 1:length(bg.Edges)
%     bg.Edges(i_bg).LineWidth = 1 + 3*log2(bg.Edges(i_bg).Weight/0.2);
%     bg.Edges(i_bg).Description = ['p = ' num2str(bg.Edges(i_bg).Weight, '%1.2f')];
% end
% 
% view(bg)
% dolayout(bg, 'layouttype', 'equilibrium')
% %%
% node_positions = [19 120; 100 120; 19 20; 100 20];c
% for i_bg = 1:length(bg.Nodes)
%     bg.Nodes(i_bg).Position = node_positions(i_bg, :);
% end
% % dolayout(bg, 'PathsOnly', true)
% view(bg)
% %% Compare 7 vs. 14 and 8 vs. 16 selectivity
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% exp_strings = {'f714f', 'f2f', 'f1020f'};
% pw_tag = ['discExp_pw' num2str(pulse_width*100) '.mat'];
% 
% selXCJ = cell(length(exp_strings), 1);
% stimResp = cell(length(exp_strings), 1);
% for i_str = 1:length(exp_strings)
%     selXCJ{i_str} = load([local_res_dir file_tag '/matfiles/selXCJ_stats_' file_tag ... 
%         exp_strings{i_str} pw_tag]);
%     
%     xxa = load([local_res_dir file_tag '/matfiles/stim_resps_' exp_strings{i_str} ... 
%         pw_tag]);
%     stimResp{i_str} = xxa.activity_struct;
% end
% %% combine stim resps
% %% 
% stim_order = [7 14 8 16 10 20];
% all_stim_endPoints = cell2mat(cellfun(@(x) [x.stim1_endPoint; x.stim2_endPoint], selXCJ([1 3]), ... 
%     'uniformoutput', false));
% 
% cell_aucs = cell2mat(cellfun(@(x) x.cell_auroc_val, selXCJ', 'UniformOutput', false));
% %% 
% all_stim1_mean_resps = cellfun(@(x) squeeze(mean(x.stim_resps(x.stim_y == 1, :, :), 1)), ... 
%     stimResp([1 3]), 'UniformOutput', false);
% all_stim2_mean_resps = cellfun(@(x) squeeze(mean(x.stim_resps(x.stim_y == 2, :, :), 1)), ... 
%     stimResp([1 3]), 'UniformOutput', false);
% 
% all_stim_mean_resps = [all_stim1_mean_resps; all_stim2_mean_resps];
% 
% [~, ord] = sort(mean(all_stim1_mean_resps{1}, 2));
% 
% all_stim_sort_resps = cellfun(@(x) x(ord, :) - all_stim1_mean_resps{1}(ord, :), all_stim_mean_resps, 'UniformOutput', false);
% %% 
% figure(),imagesc(cell2mat(all_stim_sort_resps), 0.1*[-1 1])
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 3: readout accuracy and weight distributions


keep_cls_res = all_cls_res(keep_nets);
keep_cls_einds = all_cls_ei_inds(keep_nets);
hFig3 = makeMyFigure(6.5*2.54*11/7, 4*2.54*11/7);
num_nets = length(keep_cls_res);
ex_net = 1;
i_TC = 50;
ex_cIdJ_cell = keep_cls_res{ex_net};
t_vec = 1:size(ex_cIdJ_cell{1}.classIdataJ_valAcc, 2);
plot_groups = [1 3 4];
group_colors = [0 0 0; neuron_colors([2 3], :)];

wt_lims = 0.15*[-1 1];

hA = subplot(2, 2, 1);
hold on
hB = subplot(2, 2, 2);
hold on
lh = 0*plot_groups;
for ii = 1:length(plot_groups)
    set(gcf, 'CurrentAxes', hA)
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

set(gcf, 'CurrentAxes',  hA)
text(-0.2, 1, 'A', 'units', 'normalized', 'FontSize', fsz_labels, 'FontWeight', 'bold')
set(gca, 'color', 'none', 'fontsize', fsz_axes, ...
    'box', 'off', 'ytick', 0.6:.2:1, 'xtick', 10:10:50)
pushTicksOut(hA);
% legend(lh, leg_ent, 'location', 'southeast')
text(2, 1.04, 'All', 'Color', [0 0 0], ...
    'FontSize', fsz_axes, 'HorizontalAlignment', 'left', ...
    'FontWeight', 'bold')
text(2, 0.97, 'Excitatory', 'Color', neuron_colors(2, :), ...
    'FontSize', fsz_axes, 'HorizontalAlignment', 'left', ...
    'FontWeight', 'bold')
text(2, 0.9, 'Inhibitory', 'color', neuron_colors(3, :), ...
    'FontSize', fsz_axes, 'HorizontalAlignment', 'left', ...
    'FontWeight', 'bold')

set(gcf, 'CurrentAxes', hB)
text(-0.2, 1, 'B', 'units', 'normalized', 'FontSize', fsz_labels, 'FontWeight', 'bold')

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
    
end
ylim([0.5 1.0])
xlim(x_net([1 end]) + [-0.5 0.5])
xlabel('network sim #')
ylabel('classifier accuracy')
set(gca, 'color', 'none', 'fontsize', fsz_axes, 'box', 'off')
pushTicksOut(hB);


% plot weights
exc_inds = keep_cls_einds(ex_net).exc_inds;
inh_inds = keep_cls_einds(ex_net).inh_inds;

hC = subplot(2, 2, 3);
text(-0.2, 1, 'C', 'units', 'normalized', 'FontSize', fsz_labels, 'FontWeight', 'bold')

hold on
all_exc_wts = ex_cIdJ_cell{1}.classifier_weights(:, exc_inds);
all_inh_wts = ex_cIdJ_cell{1}.classifier_weights(:, inh_inds);
plot(t_vec', mean(abs(all_exc_wts), 2), 'color', neuron_colors(2, :), 'linewidth', 1)
plot(t_vec', mean(abs(all_inh_wts), 2), 'color', neuron_colors(3, :), 'linewidth', 1)
xlabel('time at which classifier was fit')
ylabel('ave abs weight')
set(gca, 'color', 'none', 'fontsize', fsz_axes, 'box', 'off')
y_lims = ylim;
ylim([0 y_lims(2)])
% legend({'exc', 'inh'}, 'location', 'southeast')
pushTicksOut(hC);
% axis square 
% title('')

hD = subplot(2, 4, 7); cla;
text(-0.2, 1.1, 'D', 'units', 'normalized', 'FontSize', fsz_labels, 'FontWeight', 'bold')

exc_inds = arrayfun(@(x) x.exc_inds, keep_cls_einds', 'UniformOutput', false);
inh_inds = arrayfun(@(x) x.inh_inds, keep_cls_einds', 'UniformOutput', false);
exc_wts = cellfun(@(x, inds) x{1}.classifier_weights(i_TC, inds), keep_cls_res, ...
    exc_inds, 'uniformoutput', false);
inh_wts = cellfun(@(x, inds) x{1}.classifier_weights(i_TC, inds), keep_cls_res, ...
    inh_inds, 'uniformoutput', false);

sd_wtsE = cellfun(@(x) std(x), exc_wts);
mean_wtsE = cellfun(@(x) mean(x), exc_wts);
sd_wtsI = cellfun(@(x) std(x), inh_wts);
mean_wtsI = cellfun(@(x) mean(x), inh_wts);


% cdfs
wt_range = linspace(-1, 1, 201);
exc_cdfs = cell2mat(cellfun(@(x) histcounts(x, wt_range, 'Normalization', 'cdf'), exc_wts(:), ...
    'UniformOutput', false));
inh_cdfs = cell2mat(cellfun(@(x) histcounts(x, wt_range, 'Normalization', 'cdf'), inh_wts(:), ...
    'UniformOutput', false));

[~, exc_inh_ks2_test] = cellfun(@(x,y) kstest2(x,y), exc_wts, inh_wts);
[~, exc_lillie_test] = cellfun(@(x) lillietest(x), exc_wts);
[~, inh_lillie_test] = cellfun(@(x) lillietest(x), inh_wts);
%
exc_kurt = cellfun(@(x) (bootstrp(200, @(y) ...
    kurtosis(y(randperm(length(y), 100)),0), x)), exc_wts, 'UniformOutput', false);
inh_kurt = cellfun(@(x) (bootstrp(200, @(y) ...
    kurtosis(y,0), x)), inh_wts, 'Uniformoutput' ,false);

mean_exc_kurt = cellfun(@mean, exc_kurt);
mean_inh_kurt = cellfun(@mean, inh_kurt);
std_exc_kurt = cellfun(@std, exc_kurt);
std_inh_kurt = cellfun(@std, inh_kurt);


% display stats for the average 
ex_all_testacc = (ex_cIdJ_cell{1}.cvDraws_testAcc);
ex_exc_testacc = (ex_cIdJ_cell{3}.cvDraws_testAcc);
ex_inh_testacc = (ex_cIdJ_cell{4}.cvDraws_testAcc);
show_res = [1 3 4];
for i_res = show_res
    ex_acc = ex_cIdJ_cell{i_res};
    display(['Example network, ' ex_acc.cellGroup ', mean: ' ...
        num2str(mean(ex_acc.cvDraws_testAcc), '%1.2f') ', se: ' num2str(std(ex_acc.cvDraws_testAcc))])
%     display(['Example network, ' ex_acc.cellGroup ', mean: ' ...
%         num2str(mean(ex_acc), '%1.2f') ', se: ' num2str(std(ex_acc))])

end

disp(['KS2 test between E and I weights: ' num2str((exc_inh_ks2_test)', '%1.2f ')])
disp(['Lilliefors test E weights: ' num2str((exc_lillie_test)', '%1.2d ')])
disp(['Lilliefors test I weights: ' num2str((inh_lillie_test)', '%1.2d ')])


disp(['Ave kurtosis E weights: ' num2str(mean(mean_exc_kurt)', '%1.2f ')])
disp(['Ave kurtosis I weights: ' num2str(mean(mean_inh_kurt)', '%1.2f ')])
disp(['SE kurtosis E weights: ' num2str(std(mean_exc_kurt)', '%1.2f ')])
disp(['SE kurtosis I weights: ' num2str(std(mean_inh_kurt)', '%1.2f ')])

wt_bins = edges2bins(wt_range);
norm_distE_cdf = normcdf(wt_bins, mean(mean_wtsE), mean(sd_wtsE));
norm_distI_cdf = normcdf(wt_bins, mean(mean_wtsI), mean(sd_wtsI));

hold on
plot(wt_bins, exc_cdfs, 'color', neuron_colors(2, :))
plot(wt_bins, inh_cdfs, 'color', neuron_colors(3, :))
xlabel('classifier weights')
ylabel('cdf')
axis tight
xlim(wt_lims)
ylim([0 1])
set(gca, 'color', 'none', 'fontsize', fsz_axes);
title('individual networks', 'fontweight', 'normal')
pushTicksOut(hD);
axis square 


hE = subplot(2, 4, 8);
text(-0.2, 1.1, 'E', 'units', 'normalized', 'FontSize', fsz_labels, 'FontWeight', 'bold')

hold on
plot(wt_bins, mean(exc_cdfs, 1), 'color', neuron_colors(2, :), 'linewidth', 1.5)
plot(wt_bins, mean(inh_cdfs, 1), 'color', neuron_colors(3, :), 'linewidth', 1.5)
plot(wt_bins, norm_distE_cdf, 'color', 0.5*neuron_colors(2, :), 'linewidth', 1)
plot(wt_bins, norm_distI_cdf, 'color', 0.5*neuron_colors(3, :), 'linewidth', 1)
xlabel('classifier weights')
ylabel('cdf')
title('average', 'fontweight', 'normal', 'fontsize', fsz_axes-2)
axis tight
xlim(wt_lims)
ylim([0 1])
lh = legend({'e-cell wts', 'i-cell wts', 'normal cdf'}, 'color', 'none',...
    'box', 'off') ;
lh.Position = [0.83 0.2 0.17 0.11];
set(gca, 'color', 'none', 'fontsize', fsz_axes);
pushTicksOut(hE);
axis square 



print(gcf, '-dpdf', 'writeup/figures/fig3_final')


%%

hold on



subplot(3, 4, 12);
exc_inds = arrayfun(@(x) x.exc_inds, all_cls_ei_inds, 'UniformOutput', false);
inh_inds = arrayfun(@(x) x.inh_inds, all_cls_ei_inds, 'UniformOutput', false);
exc_wts = cellfun(@(x, inds) x(1).classifier_weights(i_TC, inds), all_cls_res, ...
    exc_inds, 'uniformoutput', false);
inh_wts = cellfun(@(x, inds) x(1).classifier_weights(i_TC, inds), all_cls_res, ...
    inh_inds, 'uniformoutput', false);



hold on
histogram(cell2mat(exc_wts), 'Normalization', 'pdf')
histogram(cell2mat(inh_wts), 'Normalization', 'pdf')



%% Figure 4: translated training times

makeMyFigure(6.5*2.54*1.7, 3*2.54*1.7);

num_nets = length(all_cls_res);
plot_groups = [1 3 4];

for i_net = 1:num_nets
    for i_group = 1:3
        
        valAcc_ij = all_cls_res{i_net}(i_group).classIdataJ_valAcc;
        subplot(num_nets, 3, i_group + (i_net - 1)*3)
        imagesc(valAcc_ij, [0.5 1])
        colorbar
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 3: old version. 
classifier_weights = x_class.cls_res(1).classifier_weights;
I_TC = 42;

w_example = classifier_weights(I_TC, :);
% w_example = normClassifier{I_TC}.ClassificationEnsemble.BinaryLearners{1}.Beta;
% b_ex = normClassifier{I_TC}.ClassificationEnsemble.BinaryLearners{1}.Bias;
% compute activity projections into readout space
use_inds = test_inds;
table_y = stim_y(use_inds);
table_x = squeeze(nanmean(stim_resps(use_inds, :, (I_TC-1)*skip_dt_units + (1:skip_dt_units)), 3));
data_table = table(table_x, 0*table_y);
pred_y = allClassifier{I_TC}.predictFcn(data_table);
        
% wx_resps = imfilter(squeeze(sum(bsxfun(@times, stim_resps(use_inds, :, :), w_example), 2)), ones(1, skip_dt_units)/skip_dt_units);
low_pass_resps = imfilter(stim_resps, ones(1, 1, skip_dt_units)/skip_dt_units);
wx_resps = squeeze(sum(bsxfun(@times, low_pass_resps(use_inds, :, :), w_example), 2));

%%
% 
% mu_vals = normClassifier{I_TC}.ClassificationEnsemble.BinaryLearners{1}.Mu;
% sigma_vals = normClassifier{I_TC}.ClassificationEnsemble.BinaryLearners{1}.Sigma;
% lp_resp_z = bsxfun(@times, bsxfun(@minus, low_pass_resps, mu_vals), 1./sigma_vals);
% 
% wx_resps = squeeze(sum(bsxfun(@times, lp_resp_z(use_inds, :, :), w_example), 2));

correct_trials = {pred_y.*table_y == 1, pred_y.*table_y == 4};
incorrect_trials = {pred_y == 2 & table_y == 1, pred_y == 1 & table_y == 2};

%% activity in pca-space (collapsing cell dimension)
cell_resp_ave = squeeze(mean(stim_resps(:, :, 100:500), 3));
[pc_cell, cf_cell, ~, ~, vE_cell] = pca(cell_resp_ave);
use_resps = stim_resps(use_inds, :, :);

parsedTR_ave_resp = cellfun(@(x) squeeze(mean(use_resps(x, :, :), 1)), ... 
    [correct_trials incorrect_trials], 'UniformOutput', false);
%%
parsedTR_pca_x = cellfun(@(x) pc_cell(:, 1)'*x, parsedTR_ave_resp, ... 
    'UniformOutput', false);
parsedTR_pca_y = cellfun(@(x) pc_cell(:, 2)'*x, parsedTR_ave_resp, ... 
    'UniformOutput', false);
% stim1_ave_resp = squeeze(mean(stim_resps(stim_y == 1, :, :), 1));
% stim2_ave_resp = squeeze(mean(stim_resps(stim_y == 2, :, :), 1));


%%
% stim1_pc_traj = pc_cell(:, 1:2)'*stim1_ave_resp;
% stim2_pc_traj = pc_cell(:, 1:2)'*stim2_ave_resp;
%%


[w_sorted, wt_ord] = sort(w_example);
l2_ord = l2_inds(wt_ord);
%%
% % this uses the cross-validation set for accuracy 
all_valAccs = cell2mat(arrayfun(@(x) x.classIdataJ_valAcc(I_TC, :)', x_class.cls_res, 'UniformOutput', false));

% classifier_weights has been "un-z-scored"
all_weights_raw = arrayfun(@(x) x.classifier_weights, x_class.cls_res, 'UniformOutput', false);        
all_weights = cellfun(@(x) diag(1./max(abs(x), [], 2))*x, all_weights_raw, 'UniformOutput', false);

% classifier_betas are the weights from the z-scored activity
% all_weights = arrayfun(@(x) x.classifier_betas, x_class.cls_res, 'UniformOutput', false);
%%
num_units = cellfun(@(x) size(x, 2), all_weights, 'uniformoutput', false);

%%
all_labels = arrayfun(@(x) x.cellGroup, x_class.cls_res, 'uniformoutput', false); % {'exc units', 'inh units', 'exc units', 'exc+inh units'};
%% This is OK, but should create an easier-to-read version, with ALL/Exc/Inh/SubExc 
%  separately plotted, so it isn't so many things on one panel. 
% Also, need to start showing results more networks..... 
legent = cellfun(@(x,y) [num2str(x) ' ' y], num_units, all_labels, 'UniformOutput', false);
c_time_vec = x_class.c_time_vec;
cmap = lines(size(all_valAccs, 2));

makeMyFigure(24, 16);
sh_valAcc = subplot(5, 3, [2 5]);
sh_wts = cell(4, 1);
for i_sw = 1:length(sh_wts)
    sh_wts{i_sw} = subplot(4, 3, 3*i_sw); 
end
sh_wtDist = subplot(5, 3, [11 14]);
sh_netAct = subplot(5, 3, [10 13]);
sh_pcaAct = subplot(5, 3, [1 4]);

% plot activity in pca space 
line_styles = {'-', '-.', '-', '-.'};
line_colors = [0 0 0; 0 0 0; 1 1 1; 1 1 1]*.5;
line_widths = [1.5 1.5 .5 .5 ];

[~, t_vec_ind] = min(abs(t_vec - c_time_vec(I_TC)));

set(gcf, 'currentaxes', sh_pcaAct);
hold on
for i_ptr = 1:2 %length(parsedTR_pca_x)
    plot(parsedTR_pca_x{i_ptr}, parsedTR_pca_y{i_ptr}, ...
        'linestyle', line_styles{i_ptr}, ... 
        'color', line_colors(i_ptr, :), ... 
        'LineWidth', line_widths(i_ptr));
    
    plot(parsedTR_pca_x{i_ptr}(t_vec_ind), parsedTR_pca_y{i_ptr}(t_vec_ind), ... 
        'o', 'color', cmap(4, :), 'linewidth', 3)
    
end
set(gca, 'color', 'none');
xlabel(['neuron pca dim 1 (' num2str(vE_cell(1), '%1.0f') '% of var)'])
ylabel(['neuron pca dim 2 (' num2str(vE_cell(2), '%1.0f') '% of var)'])

% plot the weighted/summed activity for correct trials
set(gcf, 'CurrentAxes', sh_netAct);
wx_mean = cellfun(@(x) mean(wx_resps(x, :), 1)', [correct_trials incorrect_trials], ...
    'uniformoutput', false);
hold on
for i_ptr = 1:length(parsedTR_pca_x)
    
    plot(t_vec, wx_mean{i_ptr}, ...
        'linestyle', line_styles{i_ptr}, ... 
        'color', line_colors(i_ptr, :), ... 
        'LineWidth', line_widths(i_ptr));
end
% mean_wx_correct = cell2mat(cellfun(@(x) mean(wx_resps(x, :), 1)', correct_trials, ...
%     'uniformoutput', false));
% sd_wx_correct = cell2mat(cellfun(@(x) mean(wx_resps(x, :), 1)', correct_trials, ...
%     'uniformoutput', false));
% mean_wx_incorrect = cell2mat(cellfun(@(x) mean(wx_resps(x, :), 1)', incorrect_trials, ...
%     'uniformoutput', false));
% sd_wx_incorrect = cell2mat(cellfun(@(x) mean(wx_resps(x, :), 1)', incorrect_trials, ...
%     'uniformoutput', false));
% hold on
% ph = plot(t_vec, mean_wx_correct, 'k', 'linewidth', 1);
% ph(1).LineStyle = '-';
% ph(2).LineStyle = '-.';
% ph2 = plot(t_vec, mean_wx_incorrect, 'color', 0.5*[1 1 1], 'linewidth', 1);
% ph2(1).LineStyle = '-';
% ph2(2).LineStyle = '-.';

xlim([0 50])
set(gca, 'color', 'none')

v = axis();
plot(c_time_vec(I_TC)*[1 1], v(3:4), 'color', cmap(4, :), 'linewidth', 1.5)
lh = legend({'low-freq (correct)', 'high-freq (correct)', 'low-freq (classifer wrong)', ... 
    'high-freq (classifier wrong)', 'classifier fit time'}, 'position', [0.15 0.42 0.1878 0.1100], ... 
    'box', 'off'); 
ylabel('weighted sum of firing rates')
xlabel('time')


set(gcf, 'CurrentAxes', sh_valAcc);
hold on
plot(c_time_vec, all_valAccs, 'LineWidth', 1.5)
hold off
set(gca, 'color', 'none')
axis square
% legend(legent, 'location', 'best')

ylabel('fraction correct (test set)')
xlabel('time relative to stim onset')

for i_wt = 1:4
    [~, ord] = sort(abs(all_weights{i_wt}(I_TC, :)), 'descend');
    these_weights = all_weights{i_wt}(I_TC, ord);
    set(gcf, 'CurrentAxes', sh_wts{i_wt});
    bar(these_weights, 'FaceColor', cmap(i_wt, :))
    title(['weights for ' num2str(num_units{i_wt}) ' ' all_labels{i_wt}])
    axis tight
    ylim([-1 1])
    ylabel('readout weight') 
    set(gca, 'color', 'none')

end

% subplot(3, 2, 2)
% bar([zeros(1, length(l1_inds)) w_sorted], 'FaceColor', cmap(i_wt, :))
% xlim([0.5 num_neur+.5])
% ylabel('readout weight')
% xlabel('neuron (sorted by weight; directly stimulated cells not included)')
% axis square
% ch = colorbar;
% ch.Visible = 'off';
% set(gca, 'color', 'none')
% title(['weights for ' num2str(num_units{i_wt}) ' ' all_labels{i_wt}])
% 

set(gcf, 'currentaxes', sh_wtDist);
hold on
br_val = 3*std(all_weights{4}(I_TC, :));
wt_bins = br_val*linspace(-1, 1, 61);
p_norm = normcdf(wt_bins', mean(all_weights{4}(I_TC, :)), std(all_weights{4}(I_TC, :)));

for i_wt = 1:4
    histogram(all_weights{i_wt}(I_TC, :), wt_bins,'Normalization', 'cdf', ...
        'DisplayStyle', 'stairs', 'LineWidth', 1.5)
end
xlabel('readout weights')
ylabel('cdf')
plot(wt_bins, p_norm, 'k', 'LineWidth', 1.5)
ylim([0 1])
xlim(wt_bins([1 end]))
set(gca, 'color', 'none')

axis square
wts_legent = cellfun(@(x,y) [num2str(x) ' ' y], num_units, all_labels, 'UniformOutput', false);
lh = legend([wts_legent, 'normal cdf'], 'Position', [0.4194 0.4146 0.1877 0.11], 'box', 'off');

% print(gcf, '-dpdf', [network_param.plotdir 'weightDists' classifier_options.classifier_type '_' network_param.network_tag 'pw' num2str(100*pulse_width)]);
print(gcf, '-dpdf', 'writeup/figures/fig3_final')



