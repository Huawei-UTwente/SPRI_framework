%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Direct collocation format of the trajectory optimization in muscle
% dynamics and reflex controllers. Neural transfer delay and 
% electricalmechanical delay are considered
%
% By: Huawei Wang
% Date: August 2, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [df_dx, df_ddx] = ...
                directCollocationDyn_diff_RPO1(x, dx, lmt, par_mus, M, S)

% function [df_dx, df_ddx, df_dx_em, dfs_dx_rf, dfs_rf_dpar_rf_fse, ...
%     dfs_rf_dpar_rf_lce, dfs_rf_dpar_rf_res] = ...
%     directCollocationDyn_diff_RPO1(x, dx, lmt, par_mus, M, S)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
    % x: the state variables at the current frame.
    % dx: the state derivatives at the current frame.
    % x_em: the state for muscle activation with the electricalmechanical delay.
    % x_rf: the state for reflex feeback with the neural transfer delay.
    % S: number of state in each muscle model, they should be the same
    %    among muscles.
    % M: number of muscles in each data trials, they should maintain the same
    %    number among trials at current optimization framwork.
    % lmt: experimental muscle tendon length at current frame.
    % lmt_rf: experimental muscle tendon length at the reflex control frame.
    % par_rf_fse: the reflex control gains for the force feedback loop
    % par_rf_lce: the reflex control gains for the fiber length feedback loop.
    % par_rf_res: the throesholds and base stimulations for the reflex
    %             controller
    % par_mus: the optimized muscle parameters, not optimize here anymore.
    % M: number of muscles
    
% OUTPUTS:
    % df_dx: derivative of constraints with respect to state at current
    % frame.
    % df_ddx: derivative of constraints with respect to state derivative
    % at current frame.
    % df_dx_em: derivative of constraints with respect to state at the
    % electricalmechanical delay frame.
    % dfs_dx_rf: derivative of constraints with respect to state at the
    % reflex feeback delay frame.
    % dfs_dpar_rf_fse: (M, M, M), derivative of estimated muscle activation respects to reflex control gains    
    % dfs_dpar_rf_lce: (M, M, M), derivative of estimated muscle activation respects to reflex control gains
    % dfs_dpar_rf_res: (M, 3*M), derivative of estimated muscle activation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % muscle model parameters
    lce_opt = par_mus(1 : M);
    lt_slack = par_mus(M+1 : 2*M);
    theta0 = par_mus(2*M+1 : 3*M);
    
%     dfs_dx_em = zeros(M, M*S);
%     dfs_dx_rf = zeros(M, M*S);
%     dfs_rf_dpar_rf_fse = zeros(M, M*M);
%     dfs_rf_dpar_rf_lce = zeros(M, M*M);
%     dfs_rf_dpar_rf_res = zeros(M, 3*M);    

    a = x(1 : M);  % muscle activation at current frame
    da = x(M+1 : 2*M);  % muscle activation differentiation at current frame
    lce = x(2*M+1 : 3*M);  % length of contract element at current frame
    dlce = x(3*M + 1 : 4*M);  % length of contract element differentiation at current frame

%     d_a = dx(1 : M);  % muscle activation differentiation at current frame
%     d_lce = dx(2*M + 1 : 3*M);  % length of contract element differentiation at current frame
    
    % constraints of finite differetiations, and their derivatives
    % fd = [da - d_a, dlce - d_lce];
    dfd_dx = zeros(2*M, S*M);
    dfd_ddx = zeros(2*M, S*M);
    
    dfd_dx(1 : M, M+1 : 2*M) = eye(M);
    dfd_ddx(1 : M, 1 : M) = -eye(M);
    
    dfd_dx(M+1 : 2*M, 3*M+1 : 4*M) = eye(M);
    dfd_ddx(M+1 : 2*M, 2*M+1 : 3*M) = -eye(M);
    
    % fa constraints are zeros in this situation
%     dfa_dx = zeros(M, M*S);
%     dfa_dx_em = zeros(M, M*S);

    % activiation nonlinearity
    % a_non = activationNonlinearty(a);
    % da_non_da = activationNonlinearty_diff(a);
   
    % muscle force dynamics
%     [~, dfm_da, dfm_dlce, dfm_ddlce] = ...
%         contractionDyn_Groote_diff_RPO(lmt, a, lce, dlce, lce_opt, lt_slack, theta0);
    
    [~, dfm_da, dfm_dlce, dfm_ddlce] = ...
        contractionDyn_Groote_diff_RPO(lmt, a, lce, dlce, lce_opt, lt_slack, theta0);
    
    dfm_dx = zeros(M, S*M);
    dfm_dx(:, 1 : M) = diag(dfm_da);
    dfm_dx(:, 2*M+1 : 3*M) = diag(dfm_dlce);
    dfm_dx(:, 3*M+1 : 4*M) = diag(dfm_ddlce);    

    % generate the combined constraints
    % f = [fd, fa, fm, fs];
    
%     df_dx = [dfd_dx; dfa_dx; dfm_dx; zeros(M, S*M)];
%     df_ddx = [dfd_ddx; zeros(M, S*M); zeros(M, S*M); zeros(M, S*M)];
%     df_dx_em = [zeros(2*M, S*M); dfa_dx_em; zeros(M, S*M); dfs_dx_em];

    df_dx = [dfd_dx; dfm_dx];
    df_ddx = [dfd_ddx; zeros(M, S*M);];
    % df_dx_em = [zeros(2*M, S*M); zeros(M, S*M); dfs_dx_em];
    
end