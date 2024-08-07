%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script comparise the results of optimizations with tolerance of 1e-3
% and 1e-4. 
%
% By: Huawei Wang
% Date: August 8, 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear

subj = 6;
phase = 2;

W2_list = [0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, ...
           0.225, 0.25, 0.275, 0.3, 0.35, 0.4, 0.45 0.5, 0.6, 0.8, 1, 1.5, 2];

% data path of 1e-4 tolerance
homeSavePath1 = sprintf('D:/HuaweiWang/ReflexIdentResults/Walk/midd/Subj%02d/%dPhases', ...
    subj, phase);

% data path of 1e-3 tolerance
homeSavePath2 = sprintf('D:/HuaweiWang/ReflexIdentResults2/Walk/midd/Subj%02d/%dPhases', ...
    subj, phase);

% load results
for iW2 = 1:length(W2_list)
    
    W2 = W2_list(iW2);
    res1 = importdata(sprintf('%s/W2_%d_W3_%d/best_res.mat', homeSavePath1, W2, 0));
    res2 = importdata(sprintf('%s/W2_%d_W3_%d/best_res.mat', homeSavePath2, W2, 0));

    states_mean1 = mean(abs(res1.states));
    
    rel_states_err(iW2, :) = mean(abs(res1.states - res2.states))./states_mean1;
    
    rel_param_err(iW2, :) = abs(res1.parameters - res2.parameters); %./(abs(res1.parameters)  + 1e-8);
    
end

    h1 = figure();
    boxplot(rel_states_err');
    hold on;
    
    h2 = figure();
    boxplot(rel_param_err');
    ylim([0, 10])
    hold on


