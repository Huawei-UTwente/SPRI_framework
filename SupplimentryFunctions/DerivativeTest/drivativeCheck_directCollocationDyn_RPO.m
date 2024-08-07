%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivative test of the direct collocation dynamics
%
% By: Huawei Wang
% Date: 12/06/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

% initlize parameters
T = 2;  % number of data trials
M = 4;  % number of muscles
S = 5;  % number of states in each muscle: a, da, lce, dlce, s
Nphase = 100;  % total phases of reflex control


%% generate initial guesses
% initialize the optimizing parameters, current frame states
a = rand(1, M);
da = 1 - 2*rand(1, M);

lce = 0.03 + 0.02*rand(1, M);
dlce = -0.05 + 0.1*rand(1, M);

s = rand(1, M);

x = [a, da, lce, dlce, s];

% derivative of current frame states
da = 1 - 2*rand(1, M);
dda = 2 - 4*rand(1, M);

dlce = -0.05 + 0.1*rand(1, M);
ddlce = -0.1 + 0.2*rand(1, M);

ds = 1 - 2*rand(1, M);

dx = [da, dda, dlce, ddlce, ds];

% states at the _em frame
a_em = rand(1, M);
da_em = 1 - 2*rand(1, M);

lce_em = 0.03 + 0.02*rand(1, M);
dlce_em = -0.05 + 0.1*rand(1, M);

s_em = rand(1, M);

x_em = [a_em, da_em, lce_em, dlce_em, s_em];

% states at the _rf frame
a_rf = rand(1, M);
da_rf = 1 - 2*rand(1, M);

lce_rf = 0.03 + 0.02*rand(1, M);
dlce_rf = -0.05 + 0.1*rand(1, M);

s_rf = rand(1, M);

x_rf = [a_rf, da_rf, lce_rf, dlce_rf, s_rf];

% muscle parameters
lce_opt0 = 0.03 + 0.02*rand(1, M);
lt_slack0 = 0.3 + 0.2*rand(1, M);
theta0 = pi/18 + pi/18*rand(1, M);
Fmax0 = 500 + 500*rand(1, M);

par_mus = [lce_opt0, lt_slack0, theta0, Fmax0];

% muscle tendon length at current and _rf frames
lmt = 0.45 + 0.15*rand(1, M);
lmt_rf  = 0.45 + 0.15*rand(1, M);

% reflex control gains: length/force reflex gains, length offset, force
% offset, and pre-stimulation.
par_rf_fse = 1 - 2*rand(M, M);
par_rf_lce = 1 - 2*rand(M, M);
par_rf_res = [rand(1, M), 0.5 + rand(1, M), 0.5*rand(1, M)];

% the derivative functions
[df_dx, df_ddx, df_dx_em, df_dx_rf, dfs_dpar_rf_fse, dfs_dpar_rf_lce, dfs_dpar_rf_res] = ...
         directCollocationDyn_diff_RPO(x, x_em, x_rf, lmt, lmt_rf,...
             par_rf_fse, par_rf_lce, par_rf_res, par_mus, M, S);

% finite differentiation
delta = 1e-4;
% differentiation of input parameters
for ia = 1:length(x)
    
   % get values with upper change
   delta_x = x(ia)*delta;
   
   x(ia) = x(ia) + delta_x;
   f_up = directCollocationDyn_RPO(x, dx, x_em, x_rf, lmt, lmt_rf,...
             par_rf_fse, par_rf_lce, par_rf_res, par_mus, M);
   
   % get values with upper change
   x(ia) = x(ia) - 2*delta_x;
   f_dw = directCollocationDyn_RPO(x, dx, x_em, x_rf, lmt, lmt_rf,...
             par_rf_fse, par_rf_lce, par_rf_res, par_mus, M);

   % change back to the original value
   x(ia) = x(ia) + delta_x;

   % calculate the finite differentiation
   df_dx_fd(:, ia) = (f_up - f_dw)/(2*delta_x);

end

% differentiation of input parameters
for ia = 1:length(dx)
    
   % get values with upper change
   delta_x = dx(ia)*delta;
   
   dx(ia) = dx(ia) + delta_x;
   f_up = directCollocationDyn_RPO(x, dx, x_em, x_rf, lmt, lmt_rf,...
             par_rf_fse, par_rf_lce, par_rf_res, par_mus, M);
   
   % get values with upper change
   dx(ia) = dx(ia) - 2*delta_x;
   f_dw = directCollocationDyn_RPO(x, dx, x_em, x_rf, lmt, lmt_rf,...
             par_rf_fse, par_rf_lce, par_rf_res, par_mus, M);

   % change back to the original value
   dx(ia) = dx(ia) + delta_x;

   % calculate the finite differentiation
   df_ddx_fd(:, ia) = (f_up - f_dw)/(2*delta_x);

end

delta = 1e-2;
% differentiation of input parameters
for ia = 1:length(x_em)
    
   % get values with upper change
   delta_x = delta;
   
   x_em(ia) = x_em(ia) + delta_x;
   f_up = directCollocationDyn_RPO(x, dx, x_em, x_rf, lmt, lmt_rf,...
             par_rf_fse, par_rf_lce, par_rf_res, par_mus, M);
   
   % get values with upper change
   x_em(ia) = x_em(ia) - 2*delta_x;
   f_dw = directCollocationDyn_RPO(x, dx, x_em, x_rf, lmt, lmt_rf,...
             par_rf_fse, par_rf_lce, par_rf_res, par_mus, M);

   % change back to the original value
   x_em(ia) = x_em(ia) + delta_x;

   % calculate the finite differentiation
   df_dx_em_fd(:, ia) = (f_up - f_dw)/(2*delta_x);
end

delta = 1e-3;
% differentiation of input parameters
for ia = 1:length(x_rf)
    
   % get values with upper change
   delta_x = x_rf(ia)*delta;
   
   x_rf(ia) = x_rf(ia) + delta_x;
   f_up = directCollocationDyn_RPO(x, dx, x_em, x_rf, lmt, lmt_rf,...
             par_rf_fse, par_rf_lce, par_rf_res, par_mus, M);
   
   % get values with upper change
   x_rf(ia) = x_rf(ia) - 2*delta_x;
   f_dw = directCollocationDyn_RPO(x, dx, x_em, x_rf, lmt, lmt_rf,...
             par_rf_fse, par_rf_lce, par_rf_res, par_mus, M);

   % change back to the original value
   x_rf(ia) = x_rf(ia) + delta_x;

   % calculate the finite differentiation
   df_dx_rf_fd(:, ia) = (f_up - f_dw)/(2*delta_x);

end

delta = 1e-3;
% differentiation of input parameters
for ia = 1:size(par_rf_fse, 1)
    for ib = 1:size(par_rf_fse, 2)
    
       % get values with upper change
       delta_x = par_rf_fse(ia, ib)*delta;

       par_rf_fse(ia, ib) = par_rf_fse(ia, ib) + delta_x;
       f_up = directCollocationDyn_RPO(x, dx, x_em, x_rf, lmt, lmt_rf,...
                 par_rf_fse, par_rf_lce, par_rf_res, par_mus, M);

       % get values with upper change
       par_rf_fse(ia, ib) = par_rf_fse(ia, ib) - 2*delta_x;
       f_dw = directCollocationDyn_RPO(x, dx, x_em, x_rf, lmt, lmt_rf,...
                 par_rf_fse, par_rf_lce, par_rf_res, par_mus, M);

       % change back to the original value
       par_rf_fse(ia, ib) = par_rf_fse(ia, ib) + delta_x;

       % calculate the finite differentiation
       dfs_dpar_rf_fse_fd(:, ia, ib) = (f_up(end-M+1:end) - f_dw(end-M+1:end))/(2*delta_x);
   
    end
end


delta = 1e-3;
% differentiation of input parameters
for ia = 1:size(par_rf_lce, 1)
    for ib = 1:size(par_rf_lce, 2)
    
       % get values with upper change
       delta_x = par_rf_lce(ia, ib)*delta;

       par_rf_lce(ia, ib) = par_rf_lce(ia, ib) + delta_x;
       f_up = directCollocationDyn_RPO(x, dx, x_em, x_rf, lmt, lmt_rf,...
                 par_rf_fse, par_rf_lce, par_rf_res, par_mus, M);

       % get values with upper change
       par_rf_lce(ia, ib) = par_rf_lce(ia, ib) - 2*delta_x;
       f_dw = directCollocationDyn_RPO(x, dx, x_em, x_rf, lmt, lmt_rf,...
                 par_rf_fse, par_rf_lce, par_rf_res, par_mus, M);

       % change back to the original value
       par_rf_lce(ia, ib) = par_rf_lce(ia, ib) + delta_x;

       % calculate the finite differentiation
       dfs_dpar_rf_lce_fd(:, ia, ib) = (f_up(end-M+1:end) - f_dw(end-M+1:end))/(2*delta_x);
   
    end
end


delta = 1e-3;
% differentiation of input parameters
for ia = 1:length(par_rf_res)
    
       % get values with upper change
       delta_x = par_rf_res(ia)*delta;

       par_rf_res(ia) = par_rf_res(ia) + delta_x;
       f_up = directCollocationDyn_RPO(x, dx, x_em, x_rf, lmt, lmt_rf,...
                 par_rf_fse, par_rf_lce, par_rf_res, par_mus, M);

       % get values with upper change
       par_rf_res(ia) = par_rf_res(ia) - 2*delta_x;
       f_dw = directCollocationDyn_RPO(x, dx, x_em, x_rf, lmt, lmt_rf,...
                 par_rf_fse, par_rf_lce, par_rf_res, par_mus, M);

       % change back to the original value
       par_rf_res(ia) = par_rf_res(ia) + delta_x;

       % calculate the finite differentiation
       dfs_dpar_rf_res_fd(:, ia) = (f_up(end-M+1:end) - f_dw(end-M+1:end))/(2*delta_x);
   
end


% check the drivative differences
tolerance = 5e-4;

% difference check of df_dx
[errorid_df_dx, diff_df_dx, rdiff_df_dx] = diffEvaluate(df_dx, df_dx_fd, tolerance);
% difference check of df_ddx
[errorid_df_ddx , diff_df_ddx, rdiff_df_ddx] = diffEvaluate(df_ddx, df_ddx_fd, tolerance);
% difference check of df_dx_em
[errorid_df_dx_em, diff_df_dx_em, rdiff_df_dx_em] = diffEvaluate(df_dx_em, df_dx_em_fd, tolerance);
% difference check of df_dx_rf
[errorid_df_dx_rf, diff_df_dx_rf, rdiff_df_dx_rf] = diffEvaluate(df_dx_rf, df_dx_rf_fd, tolerance);
% difference check of df_dpar_rf_fse
[errorid_df_dpar_rf_fse, diff_df_dpar_rf_fse, rdiff_df_dpar_rf_fse] = ...
    diffEvaluate(dfs_dpar_rf_fse, dfs_dpar_rf_fse_fd, tolerance);
% difference check of df_dpar_rf_lce
[errorid_df_dpar_rf_lce, diff_df_dpar_rf_lce, rdiff_df_dpar_rf_lce] = ...
    diffEvaluate(dfs_dpar_rf_lce, dfs_dpar_rf_lce_fd, tolerance);
% difference check of df_dpar_rf_res
[errorid_df_dpar_rf_res, diff_df_dpar_rf_res, rdiff_df_dpar_rf_res] = ...
    diffEvaluate(dfs_dpar_rf_res, dfs_dpar_rf_res_fd, tolerance);

if ~isempty(errorid_df_dx)
   fprintf('Differentiations in df_dx beyond thresholds\n')
end

if ~isempty(errorid_df_ddx)
   fprintf('Differentiations in df_ddx beyond thresholds\n')
end

if ~isempty(errorid_df_dx_em)
   fprintf('Differentiations in df_dx_em beyond thresholds\n')
end

if ~isempty(errorid_df_dx_rf)
   fprintf('Differentiations in df_dx_rf beyond thresholds\n')
end

if ~isempty(errorid_df_dpar_rf_fse)
   fprintf('Differentiations in df_dpar_rf_fse beyond thresholds\n')
end

if ~isempty(errorid_df_dpar_rf_lce)
   fprintf('Differentiations in df_dpar_rf_lce beyond thresholds\n')
end

if ~isempty(errorid_df_dpar_rf_res)
   fprintf('Differentiations in df_dpar_rf_res beyond thresholds\n')
end