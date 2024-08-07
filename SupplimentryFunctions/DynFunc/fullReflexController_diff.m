%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code define reflex control derivatives. 
%
% The reflex control gains vary along with the gait phases, from 1 to 100
% 
% By: Huawei Wang
% Date: August 15, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [du_dlce, du_dfse, du_dpar_rf_fse, du_dpar_rf_lce, du_dpar_rf_res] ...
    = fullReflexController_diff(lce, fse, par_rf_fse, par_rf_lce, par_rf_res, M)
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
    %	du_dlce: (M, M), derivative of estimated muscle activation respects to muscle length
    %	du_dfse: (M, M), derivative of estimated muscle activation respects to muscle force
    %	du_dpar_rf_fse: (M, M*M), derivative of estimated muscle activation respects to reflex control gains    
    %	du_dpar_rf_lce: (M, M*M), derivative of estimated muscle activation respects to reflex control gains
    %	du_dpar_rf_res: (M, 3*M), derivative of estimated muscle activation
    %	respects to reflex control thresholds and pre-stimulation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    smooth_delta = 1e-6; % smooth parameters

    tfse = par_rf_res(1:M);  % throshold of fse
    tlce = par_rf_res(M+1:2*M); % throshold of lce
    % u0 = par_rf_res(2*M+1:3*M);  % the baseline activation (stimulation)

    % define reflex output, derivative matrixes
    du_dpar_rf_fse = zeros(M, M*M);
    du_dpar_rf_lce = zeros(M, M*M);
    du_dpar_rf_res = zeros(M, 3*M);
    
    % smooth by fse and lce throshold
    [sfse, dsfse_dfse, dsfse_dtfse] = ...
                    substract_positive_smooth_diff(fse, tfse, smooth_delta);
    [slce, dslce_dlce, dslce_dtlce] = ...
                    substract_positive_smooth_diff(lce, tlce, smooth_delta);
    
    % derivative of the fse and lce from the output u
    du_dfse = dsfse_dfse.*reshape(par_rf_fse, M, M)';
    du_dlce = dslce_dlce.*reshape(par_rf_lce, M, M)';
	
    for i = 1:M
        du_dpar_rf_fse(i, (i-1)*M+1:i*M) = sfse;
        du_dpar_rf_lce(i, (i-1)*M+1:i*M) = slce;
    end

    du_dpar_rf_res(:, 1:M) = dsfse_dtfse.*reshape(par_rf_fse, M, M)';
    du_dpar_rf_res(:, M+1:2*M) = dslce_dtlce.*reshape(par_rf_lce, M, M)';
    du_dpar_rf_res(:, 2*M+1:3*M) = eye(M);

end