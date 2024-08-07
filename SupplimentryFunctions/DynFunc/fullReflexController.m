%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code define the fully connected reflex control loop.
% 
% The reflex control gains vary along with the gait phases, from 1 to 100
% 
% By: Huawei Wang
% Date: August 15, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = fullReflexController(lce, fse, par_rf_fse, par_rf_lce, par_rf_res, M)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Inputs:
    %	lce: (1 x M), muscle fiber length, unit: lce_opt
    %	fse: (1 x M), muscle force, unit: Fmax
    %	par_rf_fse: (1 x M*M), reflex control gains, 1/Fmax
    %	par_rf_lce: (1 x M*M), reflex control gains, 1/lce_opt
    %   par_rf_res: (3 x M), reflex thoresholds and base stimulations
    %   M: integer, number of muscles
    %
    % Outputs:
    %	u: (1 x M), estimated muscle activation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    smooth_delta = 1e-6; % smooth parameters

    tfse = par_rf_res(1:M);  % throshold of fse
    tlce = par_rf_res(M+1:2*M); % throshold of lce
    u0 = par_rf_res(2*M+1:3*M);  % the baseline activation (stimulation)

    % reflex control
    % smooth by fse and lce throshold
    sfse = substract_positive_smooth(fse, tfse, smooth_delta);
    slce = substract_positive_smooth(lce, tlce, smooth_delta);

    % calculate the activation from the reflex control loop
    u = u0 + sfse*reshape(par_rf_fse, M, M) + slce*reshape(par_rf_lce, M, M);

end