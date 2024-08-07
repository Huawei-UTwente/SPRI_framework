%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is to plot the 3D-reflexes gains averaging between N subjects,
% M walking speeds.
%
% By: Huawei Wang
% Date: 04-19-2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load optimized results
summarizing_path = 'ResultsSummarize/';
subj1trial1 = load(sprintf('%s1subj1trial.mat', summarizing_path));
subj1trialAll = load(sprintf('%s1subjAlltrial.mat', summarizing_path));
subjAlltrialAll = load(sprintf('%sAllsubjAlltrial.mat', summarizing_path));

% calculate average and std of the identified reflex gains
ave_subj1trial1_prf = mean(subj1trial1.Prf)
std_subj1trial1_prf = std(subj1trial1.Prf)

ave_subj1trialAll_prf = mean(subj1trialAll.Prf)
std_subj1trialAll_prf = std(subj1trialAll.Prf)

ave_subjAlltrialAll_prf = mean(subjAlltrialAll.Prf)
std_subjAlltrialAll_prf = std(subjAlltrialAll.Prf)

% plot the 3d reflex gains (mean + std)
phase_occupied = [1, 2, 1, 1.5, 1.5];  % the phase lengths of each period
cphase = 5;

hf = reflex_3d_map(ave_subj1trial1_prf, M, cphase, phase_occupied);

hf = reflex_3d_map(ave_subj1trialAll_prf, M, cphase, phase_occupied);

hf = reflex_3d_map(ave_subjAlltrialAll_prf, M, cphase, phase_occupied);