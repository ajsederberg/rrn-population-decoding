function readoutOutputs = analyzeSimulationResults(TOUT, ROUT, input_params, sim_params)
% THis function takes teh output of the network simulationg and the
% stimulus information and analyzes the results. 

numNeurExc = sim_params.numNeurExc;
numNeurInh = sim_params.numNeurInh;
numNeur = size(sim_params.J, 1);


input_t = sim_params.input_t;

input_times = input_params{1};

% these are the start times of each stimulus (stim 1, stim 2). These are
% used for calculating STAs and peri-stimulus responses. You may wish to
% generalize this to handle more than 2 stims. 
inp1_times = input_times{1};
inp2_times = input_times{2};

first_stim_time = min(cellfun(@min, input_times));
burn_in_T = max(1, sum(TOUT < (first_stim_time - 20)));

% z_ROUT = zscore(ROUT(burn_in_T:end, :), [], 1);

%%
% hFig = figure();
% hFig.OuterPosition = [20 20 1200 900];
% hFig.InnerPosition = [20 20 1200 880];
makeMyFigure(40, 30);
%%
subplot(4, 2, 2)
cla;
% imagesc(TOUT(burn_in_T:end), [], ROUT(burn_in_T:end, :)', [0 1])
% xlabel('time')
% ylabel('z-scored output')
% colorbar

wE1 = [ones(numNeurExc, 1)/numNeurExc; zeros(numNeurExc + numNeurInh, 1)];
wE2 = [zeros(numNeurExc, 1); ones(numNeurExc, 1)/numNeurExc; zeros(numNeurInh, 1)];
wI = [zeros(2*numNeurExc, 1); ones(numNeurInh, 1)/numNeurInh];

yEEI =  ROUT*[wE1 + wE2 wI];

% subplot(4, 2, 4);
hold on
plot(TOUT(burn_in_T:end), yEEI(burn_in_T:end, :))
hold on
plot(input_t(inp2_times), 0*inp2_times, '*')
plot(input_t(inp1_times), 0*inp1_times, 'o')
axis tight
xlim(TOUT([burn_in_T end]))
ylabel('summed activity (E blue, I red)')
% compute STA of input 1 and input 2 generated activity 
peristim_time = -1:.01:3;
%%
% define a function to pull out the responses near the stim times
peristim_resp_fun =  @(inp_times, peristim_time) cell2mat( ...
    cellfun(@(x) shiftdim( interp1(TOUT, ROUT, x + (peristim_time)), -1), ...
    num2cell(input_t(inp_times)'), ...
    'uniformoutput', false));
% compute the average stim-triggered responses
sta_fun = @(inp_times, peristim_time) ...
    squeeze(mean(peristim_resp_fun(inp_times, peristim_time) , 1) );

sd_sta_fun = @(inp_times, peristim_time) ...
    squeeze(std( peristim_resp_fun(inp_times, peristim_time), [], 1) );

sta_resp1 = sta_fun(inp1_times, peristim_time);
sta_resp2 = sta_fun(inp2_times, peristim_time);

sta_resp1_ev = bsxfun(@minus, sta_resp1, sta_resp1(1, :));
sta_resp2_ev = bsxfun(@minus, sta_resp2, sta_resp1(1, :));

sd_sta_resp1 = sd_sta_fun(inp1_times, peristim_time);
sd_sta_resp2 = sd_sta_fun(inp2_times, peristim_time);

resp1_neur = ( max(sta_resp1, [], 1) - sta_resp1(1, :) > sd_sta_resp1(1, :));
resp2_neur = ( max(sta_resp2, [], 1) - sta_resp2(1, :) > sd_sta_resp2(1, :));


%%
[~, ind_peak] = max(sum(sta_resp1 + sta_resp2, 2));
postT = peristim_time(ind_peak);
resp1_mat = interp1(TOUT, ROUT, input_t(inp1_times) + postT) ;%- interp1(TOUT, ROUT, input_t(inp1_times));

resp2_mat = interp1(TOUT, ROUT, input_t(inp2_times) + postT) ;%- interp1(TOUT, ROUT, input_t(inp2_times));


x_data = [resp1_mat; resp2_mat];
y_data = [ones(size(resp1_mat, 1), 1); 2*ones(size(resp2_mat, 1), 1)];


% data table from entire population
table_x = x_data;
table_y = y_data;
sim_data_table = table(table_x, table_y);
% try decoding from I population only
table_x = x_data(:, 2*numNeurExc + 1 : end);
table_y = y_data;
inh_data_table = table(table_x, table_y);


J = param.J;
eigJ = eig(J);
subplot(4, 4, 1)
plot(eigJ, 'o')
hold on
plot( [1 1], [-1 1])
axis equal
xlabel('re \lambda')
ylabel('im \lambda')

subplot(4, 4, 2)
plot(1:numNeur, mean(J, 2), 'o')
xlabel('neuron #')
ylabel('sum input weights')
axis square

subplot(4, 4, 5)
imagesc(J, [-1 1])
colorbar
axis square
xlabel('input i')
ylabel('output j')

subplot(4, 4, 6)
set(gca, 'visible', 'off')
emp_eiratio = mean(-sum(param.J(:, 1:2*numNeurExc), 2)./sum(param.J(:, 2*numNeurExc+1:end), 2));
text(0, 1, ['E/I ratio : ' num2str(emp_eiratio, '%1.2f')])
text(0, .8, ['i/o function: ' func2str(fun_io)])
text(0, .6, ['network \tau: ' num2str(param.net_tau)])
text(0, .4, ['stimulus integration window: ' num2str(param.trial_length)])


subplot(4, 4, 9)
plot(peristim_time,  sta_resp1_ev(:, 1:2*numNeurExc)')
xlabel('time')
ylabel('exc neuron STA responses')
title('stim 1')
subplot(4, 4, 13)
plot(peristim_time,  sta_resp1_ev(:, 2*numNeurExc+1:end)')
xlabel('time')
ylabel('inh neuron STA responses')
title('stim 1')

subplot(4, 4, 10)
plot(peristim_time, sta_resp2_ev(:, 1:2*numNeurExc)')
xlabel('time')
ylabel('exc neuron STA responses')
title('stim 2') 
subplot(4, 4, 14)
plot(peristim_time, sta_resp2_ev(:, 2*numNeurExc+1:end)')
xlabel('time')
ylabel('inh neuron STA responses')
title('stim 2') 
% subplot(4, 4, 11)
% cla;
% hold on
% plot(yEEI(burn_in_T:end, 1)+yEEI(burn_in_T:end, 2), yEEI(burn_in_T:end, 3))
% plot(resp1_mat*[wE1 + wE2], resp1_mat*wI, 'o')
% plot(resp2_mat*[wE1 + wE2], resp2_mat*wI, 'o')
% xlabel('sum E population')
% ylabel('sum I population')


%% fit classifiers at one time point;  test over the full 
mean_resp1 = mean(resp1_mat, 1);
mean_resp2 = mean(resp2_mat, 1);

w_diffmean = mean_resp1 - mean_resp2;

inh_inds = (2*numNeurExc + 1) : numNeur;
inh_w = w_diffmean(inh_inds);
exc_inds = randperm(numNeurExc*2, numNeurInh);
exc_w = w_diffmean(exc_inds);

resp_vt_stim1 = peristim_resp_fun(inp1_times, peristim_time);
resp_vt_stim2 = peristim_resp_fun(inp2_times, peristim_time);

w_dot_inhresp1 = sum(bsxfun(@times, resp_vt_stim1(:, :, inh_inds), reshape(inh_w, [1 1 numNeurInh])), 3);
w_dot_excresp1 =  sum(bsxfun(@times, resp_vt_stim1(:, :, exc_inds), reshape(exc_w, [1 1 numNeurInh])), 3);
w_dot_inhresp2 = sum(bsxfun(@times, resp_vt_stim2(:, :, inh_inds), reshape(inh_w, [1 1 numNeurInh])), 3);
w_dot_excresp2 =  sum(bsxfun(@times, resp_vt_stim2(:, :, exc_inds), reshape(exc_w, [1 1 numNeurInh])), 3);

%%
auroc_vt_E = zeros(1, size(w_dot_excresp1, 2));
auroc_vt_I = zeros(1, size(w_dot_inhresp1, 2));
for i_t = 1:length(auroc_vt_E)
    auroc_vt_E(i_t) = auroc(w_dot_excresp2(:, i_t), w_dot_excresp1(:, i_t), 0);
    auroc_vt_I(i_t) = auroc(w_dot_inhresp2(:, i_t), w_dot_inhresp1(:, i_t), 0);
    
end

%%
subplot(4, 4, 15)
hold on
plot(peristim_time, w_dot_inhresp1', 'k')
plot(peristim_time, w_dot_inhresp2', 'r')
xlabel('time')
ylabel('stim response ( w*r_I )')
subplot(4, 4, 16)
bar(sort(inh_w, 'descend'))
ylabel('weights (inh)')

subplot(4, 4, 11)
hold on
plot(peristim_time, w_dot_excresp1', 'k')
plot(peristim_time,w_dot_excresp2', 'r')
xlabel('time')
ylabel('stim response ( w*r_E )')
subplot(4, 4, 12)
bar(sort(exc_w, 'descend'))
ylabel('weights (exc)')

subplot(4, 4, 7)
cla;
plot(peristim_time, auroc_vt_E)
hold on
plot(peristim_time, auroc_vt_I)
ylabel('auroc vs t under LDA classification')
legend({'excitatory', 'inhibitory'})
%% fit classifiers
%  to use the generic fitting function, table columns have names "table_x"
%  and "table_y"
table_x = sim_data_table.table_x(:, randperm(2*param.numNeurExc, param.numNeurInh));
table_y = sim_data_table.table_y;
exc_data_table = table(table_x, table_y);


classifier_options.classifier_type = 'pLDA';
classifier_options.applyPCA = false;
classifier_options.full_save = true;

%  ok, with these settings: 89% accuracy. Also true if the inh neurons receive no direct
%  inputs. Should replace this function with my more generalized version. 
[inhClassifier, inhValidationAccuracy] = trainGenericClassifier(inh_data_table, classifier_options);
[excClassifier, excValidationAccuracy] = trainGenericClassifier(exc_data_table, classifier_options);

readoutOutputs.inhClassifier = inhClassifier;
readoutOutputs.inhValidationAccuracy = inhValidationAccuracy;
readoutOutputs.excClassifier = excClassifier;
readoutOutputs.excValidationAccuracy = excValidationAccuracy;

readoutOutputs.trialBytrial_stim1 = resp_vt_stim1;
readoutOutputs.trialBytrial_stim2 = resp_vt_stim2;


print(gcf, '-dpdf', ['pdf printouts/rrn_driven_discrim_' param.sim_name])