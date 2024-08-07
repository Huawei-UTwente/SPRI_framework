%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivative test of the tendon force
%
% By: Huawei Wang
% Date: 12/06/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

% initlize parameters
M = 2;  % number of muscles
lmt = [0.45, 0.35].*(0.9 + 0.2*rand(1, 2));  % muscle lengths
lce = [0.08, 0.06].*(0.9 + 0.2*rand(1, 2));  % drivative of muscle activation
lce_opt = [0.06, 0.07];  % muscle excitation
theta0 = [pi/8, pi/6];  % activation time constant
lt_slack = [0.4, 0.3];  % deactivation time constant

% calculate tendon force derivatives
[~, dFse_dlce]...
             = tendenForce_Groote_diff_RPO(lmt, lce, lce_opt, lt_slack, theta0);

dFse_dlce_equ = diag(dFse_dlce);

% finite differentiation
delta = 1e-4;

% differentiation of muscle fiber lengths lce
for ia = 1:length(lce)
    
   % get values with upper change
   delta_lce = lce(ia)*delta;
   
   lce(ia) = lce(ia) + delta_lce;
   Fse_up = tendenForce_Groote_RPO(lmt, lce, lce_opt, lt_slack, theta0);
   
   % get values with upper change
   lce(ia) = lce(ia) - 2*delta_lce;
   Fse_dw = tendenForce_Groote_RPO(lmt, lce, lce_opt, lt_slack, theta0);

   % change back to the original value
   lce(ia) = lce(ia) + delta_lce;

   % calculate the finite differentiation
   dFse_dlce_fd(:, ia) = (Fse_up - Fse_dw)/(2*delta_lce);

end

% check the drivative differences
tolerance = 1e-4;

% difference check of dFse_dlce
errorid_dFse_dlce = diffEvaluate(dFse_dlce_equ, dFse_dlce_fd, tolerance);

if ~isempty(errorid_dFse_dlce)
   fprintf('Differentiations in dFse_dlce beyond thresholds\n')
end