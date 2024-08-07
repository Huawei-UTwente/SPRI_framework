%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is to check the derivative of the reflex-synergy control
%
% By: Huawei Wang
% Date: August 20, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

M = 4;  % number of muscles
totalPhase = 100;  % total number of the gait phases
Prf = M*M*totalPhase*2+3*M;  % number of reflex control parameters

lce = 0.05 + 0.05*rand(1, M);  % randomize muscle fiber lengths
fse = 800 + 400*rand(1, M);    % randomize muscle forces

lce_opt = zeros(1, M) + 0.08;  % optimal fiber lengths
Fmax = zeros(1, M) + 1200;     % maximum muscle forces

lce_rel = lce./lce_opt/2;  % muscle relative fiber lengths, devide 2 to keep in the range in [0 - 1];
fse_rel = fse./Fmax/1.5;  % muscle relative forces, devide 1.5 to keep in the range in [0 - 1];

% initialize the reflex control parameters: gains [-5, 5]; ref [0- 1]
% initialize the synergy control parameters: h:[0-1], b:[0-1], c:[0-1],
% weights[0-1]
par_rf_fse = -10 + 5*rand(1, M*M);
par_rf_lce = -10 + 5*rand(1, M*M);
par_rf_res = rand(1, 3*M);

% calcuate reflex control derivatives given the current controllers;
[dur_dlce, dur_dfse, dur_dpar_rf_fse, dur_dpar_rf_lce, dur_dpar_rf_res] ...
    = fullReflexController_diff(lce_rel, fse_rel, par_rf_fse, par_rf_lce, par_rf_res, M);
                                                    
% finite differenances of the reflex and synergy controls
delta = 1e-4;

% reflex control
dur_dlce_pd = zeros(M, M);
dur_dfse_pd = zeros(M, M);
dur_dpar_rf_fse_pd = zeros(M, M);
dur_dpar_rf_lce_pd = zeros(M, M);
dur_dpar_rf_res_pd = zeros(M, 3*M);

for m = 1:M
    
    delta_x = lce_rel(m)*delta;
    
    lce_rel(m) = lce_rel(m) + delta_x;
    u_p = fullReflexController(lce_rel, fse_rel, par_rf_fse, par_rf_lce, par_rf_res, M);
    
    lce_rel(m) = lce_rel(m) - 2*delta_x;
    u_m = fullReflexController(lce_rel, fse_rel, par_rf_fse, par_rf_lce, par_rf_res, M);
    
    dur_dlce_pd(:, m) = (u_p - u_m)/(2*delta_x);
    
    lce_rel(m) = lce_rel(m) + delta_x;
                                                
end   

for m = 1:M
    
    delta_x = fse_rel(m)*delta;
    
    fse_rel(m) = fse_rel(m) + delta_x;
    u_p = fullReflexController(lce_rel, fse_rel, par_rf_fse, par_rf_lce, par_rf_res, M);
    
    fse_rel(m) = fse_rel(m) - 2*delta_x;
    u_m = fullReflexController(lce_rel, fse_rel, par_rf_fse, par_rf_lce, par_rf_res, M);
    
    dur_dfse_pd(:, m) = (u_p - u_m)/(2*delta_x);
    
    fse_rel(m) = fse_rel(m) + delta_x;
                                                
end

for p = 1:M*M

    delta_x = par_rf_fse(p)*delta;
    
    par_rf_fse(p) = par_rf_fse(p) + delta_x;
    u_p = fullReflexController(lce_rel, fse_rel, par_rf_fse, par_rf_lce, par_rf_res, M);
    
    par_rf_fse(p) = par_rf_fse(p) - 2*delta_x;
    u_m = fullReflexController(lce_rel, fse_rel, par_rf_fse, par_rf_lce, par_rf_res, M);
    
    dur_dpar_rf_fse_pd(:, p) = (u_p - u_m)/(2*delta_x);
    
    par_rf_fse(p) = par_rf_fse(p) + delta_x;
                                             
end

for p = 1:M*M

    delta_x = par_rf_lce(p)*delta;
    
    par_rf_lce(p) = par_rf_lce(p) + delta_x;
    u_p = fullReflexController(lce_rel, fse_rel, par_rf_fse, par_rf_lce, par_rf_res, M);
    
    par_rf_lce(p) = par_rf_lce(p) - 2*delta_x;
    u_m = fullReflexController(lce_rel, fse_rel, par_rf_fse, par_rf_lce, par_rf_res, M);
    
    dur_dpar_rf_lce_pd(:, p) = (u_p - u_m)/(2*delta_x);
    
    par_rf_lce(p) = par_rf_lce(p) + delta_x;
                                       
end

for p = 1:3*M
    
    delta_x = par_rf_res(p)*delta;
    
    par_rf_res(p) = par_rf_res(p) + delta_x;
    u_p = fullReflexController(lce_rel, fse_rel, par_rf_fse, par_rf_lce, par_rf_res, M);
    
    par_rf_res(p) = par_rf_res(p) - 2*delta_x;
    u_m = fullReflexController(lce_rel, fse_rel, par_rf_fse, par_rf_lce, par_rf_res, M);
    
    dur_dpar_rf_res_pd(:, p) = (u_p - u_m)/(2*delta_x);
    
    par_rf_res(p) = par_rf_res(p) + delta_x;
                                                
end

tolerance = 1e-4;

% difference check of dur_dlce
[errorid_dur_dlce, diff_dur_dlce, rdiff_dur_dlce] = diffEvaluate(dur_dlce, dur_dlce_pd, tolerance);
% difference check of dur_dfse
[errorid_dur_dfse, diff_dur_dfse, rdiff_dur_dfse] = diffEvaluate(dur_dfse, dur_dfse_pd, tolerance);
% difference check of dur_dpar_rf_res
[errorid_dur_dpar_rf_fse, diff_dur_dfpar_rf_fse, rdiff_dur_dfpar_rf_fse] = diffEvaluate(dur_dpar_rf_fse, dur_dpar_rf_fse_pd, tolerance);
% difference check of dur_dpar_rf_res
[errorid_dur_dpar_rf_lce, diff_dur_dfpar_rf_lce, rdiff_dur_dfpar_rf_lce] = diffEvaluate(dur_dpar_rf_res, dur_dpar_rf_res_pd, tolerance);
% difference check of dur_dpar_rf_res
[errorid_dur_dpar_rf_res, diff_dur_dfpar_rf_res, rdiff_dur_dfpar_rf_res] = diffEvaluate(dur_dpar_rf_res, dur_dpar_rf_res_pd, tolerance);

if ~isempty(errorid_dur_dlce)
   fprintf('Differentiations in dur_dlce beyond thresholds\n')
end

if ~isempty(errorid_dur_dfse)
   fprintf('Differentiations in dur_dfse beyond thresholds\n')
end

if ~isempty(errorid_dur_dpar_rf_fse)
   fprintf('Differentiations in dpar_rf_fse beyond thresholds\n')
end

if ~isempty(errorid_dur_dpar_rf_lce)
   fprintf('Differentiations in dpar_rf_lce beyond thresholds\n')
end

if ~isempty(errorid_dur_dpar_rf_res)
   fprintf('Differentiations in dpar_rf_res beyond thresholds\n')
end