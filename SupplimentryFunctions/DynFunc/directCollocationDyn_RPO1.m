%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Direct collocation format of the trajectory optimization in muscle
% dynamics and reflex controllers. The most basic muscle dynamics.
%
% By: Huawei Wang
% Date: August 2, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function f = directCollocationDyn_RPO1(x, dx, lmt, par_mus, M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
    % x: the state variables at the current frame.
    % dx: the state derivatives at the current frame.
    % lmt: experimental muscle tendon length at current frame.
    % par_mus: the optimized muscle parameters, not optimize here anymore.
    % M: number of muscles
% OUTPUTS:
    % f: constraints of reflex control and muscle dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % muscle model parameters
    lce_opt = par_mus(1 : M);
    lt_slack = par_mus(M+1 : 2*M);
    theta0 = par_mus(2*M+1 : 3*M);
    
    % fs = zeros(1, M);
    
    a = x(1 : M);  % muscle activation at current frame
    da = x(M+1 : 2*M);  % muscle activation differentiation at current frame
    lce = x(2*M+1 : 3*M);  % length of contract element at current frame
    dlce = x(3*M + 1 : 4*M);  % length of contract element differentiation at current frame

    d_a = dx(1 : M);  % muscle activation differentiation at current frame
    d_lce = dx(2*M + 1 : 3*M);  % length of contract element differentiation at current frame

    % constraints of finite differetiations
    fd = [da - d_a, dlce - d_lce];

    % muscle activation dynamics, not included, since no neural stimulation
    % inputs
    % fa = zeros(1, M); % activationDyn_Groote(a, da, s_em);

    % activiation nonlinearity, [not needed for the leg muscles]
    % a_non = activationNonlinearty(a);
    
    % muscle force dynamics
%     fm = contractionDyn_Groote(lmt, a, lce, dlce, lce_opt, ...
%                                lt_slack, theta0);
                           
    fm = contractionDyn_Groote_RPO(lmt, a, lce, dlce, lce_opt, ...
                               lt_slack, theta0);

    % generate the combined constraints
    % f = [fd, fa, fm, fs];
    f = [fd, fm];

end