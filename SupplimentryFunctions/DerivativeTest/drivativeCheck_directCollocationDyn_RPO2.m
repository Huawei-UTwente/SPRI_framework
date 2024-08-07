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
C = 5;  % constraints of each muscle
Nphase = 5;  % total phases of reflex control

%% generate initial guesses
% initialize the optimizing parameters, current frame states
a = 0.8*rand(1, M) + 0.01;
da = 1 - 2*rand(1, M);

lce = 0.03 + 0.02*rand(1, M);
dlce = -0.1 + 0.2*rand(1, M);

s = 0.8*rand(1, M) + 0.01;

x = [a, da, lce, dlce, s];

% derivative of current frame states
da = 1 - 2*rand(1, M);
dda = 2 - 4*rand(1, M);

dlce = -0.5 + 0.2*rand(1, M);
ddlce = -0.2 + 0.4*rand(1, M);

ds = 1 - 2*rand(1, M);

dx = [da, dda, dlce, ddlce, ds];

% states at the _em frame
a_em = 0.8*rand(1, M) + 0.01;
da_em = 1 - 2*rand(1, M);

lce_em = 0.03 + 0.02*rand(1, M);
dlce_em = -0.1 + 0.2*rand(1, M);

s_em = 0.8*rand(1, M) + 0.01;

x_em = [a_em, da_em, lce_em, dlce_em, s_em];

% muscle parameters
lce_opt0 = 0.035 + 0.01*rand(1, M);
lt_slack0 = 0.4 + 0.05*rand(1, M);
theta0 = pi/18 + pi/18*rand(1, M);
Fmax0 = 500 + 500*rand(1, M);

par_mus = [lce_opt0, lt_slack0, theta0, Fmax0];

% muscle tendon length at current and _rf frames
lmt = 0.45 + 0.05*rand(1, M);
lmt_rf  = 0.45 + 0.05*rand(1, M);

% the derivative functions
[df_dx, df_ddx, df_dx_em] = directCollocationDyn_diff_RPO2(x, dx, x_em, ...
                                                        lmt, par_mus, M, S);

% finite differentiation
delta = 1e-4;
% differentiation of input parameters
for ia = 1:length(x)
    
   % get values with upper change
   delta_x = x(ia)*delta;
   
   x(ia) = x(ia) + delta_x;
   f_up = directCollocationDyn_RPO2(x, dx, x_em, lmt, par_mus, M);
   
   % get values with upper change
   x(ia) = x(ia) - 2*delta_x;
   f_dw = directCollocationDyn_RPO2(x, dx, x_em, lmt, par_mus, M);

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
   f_up = directCollocationDyn_RPO2(x, dx, x_em, lmt, par_mus, M);
   
   % get values with upper change
   dx(ia) = dx(ia) - 2*delta_x;
   f_dw = directCollocationDyn_RPO2(x, dx, x_em, lmt, par_mus, M);

   % change back to the original value
   dx(ia) = dx(ia) + delta_x;

   % calculate the finite differentiation
   df_ddx_fd(:, ia) = (f_up - f_dw)/(2*delta_x);

end

delta = 1e-4;
% differentiation of input parameters
for ia = 1:length(x_em)
    
   % get values with upper change
   delta_x = delta;
   
   x_em(ia) = x_em(ia) + delta_x;
   f_up = directCollocationDyn_RPO2(x, dx, x_em, lmt, par_mus, M);
   
   % get values with upper change
   x_em(ia) = x_em(ia) - 2*delta_x;
   f_dw = directCollocationDyn_RPO2(x, dx, x_em, lmt, par_mus, M);

   % change back to the original value
   x_em(ia) = x_em(ia) + delta_x;

   % calculate the finite differentiation
   df_dx_em_fd(:, ia) = (f_up - f_dw)/(2*delta_x);
end

% check the drivative differences
tolerance = 1e-4;

% difference check of df_dx
[errorid_df_dx, diff_df_dx, rdiff_df_dx] = diffEvaluate(df_dx, df_dx_fd, tolerance);
% difference check of df_ddx
[errorid_df_ddx , diff_df_ddx, rdiff_df_ddx] = diffEvaluate(df_ddx, df_ddx_fd, tolerance);
% difference check of df_dx_em
[errorid_df_dx_em, diff_df_dx_em, rdiff_df_dx_em] = diffEvaluate(df_dx_em, df_dx_em_fd, tolerance);

if ~isempty(errorid_df_dx)
   fprintf('Differentiations in df_dx beyond thresholds\n')
end

if ~isempty(errorid_df_ddx)
   fprintf('Differentiations in df_ddx beyond thresholds\n')
end

if ~isempty(errorid_df_dx_em)
   fprintf('Differentiations in df_dx_em beyond thresholds\n')
end