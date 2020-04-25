function [w, b, beta_val, bias_val] = linearSVMWeightsUnZscored(binary_learner)
% a funciton to return the weights (w.x > b) that perform the
% classification specified by (beta.z > bias) where z = (x - mu)/s is the
% element-wise z-scored activity. 

if iscell(binary_learner)
    binary_learner = binary_learner{1};
end

beta_val = binary_learner.Beta';
bias_val = binary_learner.Bias;

mu_val = binary_learner.Mu;
sigma_val = binary_learner.Sigma;

w = beta_val./sigma_val;
b = bias_val + sum(beta_val.*mu_val./sigma_val);

% normalize so max(abs(w)) = 1
b = b/max(abs(w));
w = w/max(abs(w));