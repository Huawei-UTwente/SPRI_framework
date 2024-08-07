function obj =  objective_RPO_L0(x, M, S, N, nPar_rf, lmt, par_mus,...
     tPhase, W1, W2, W3, torque_m, act_m, d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The callback for calculating the objective
% Input: 
%        x: optimizing states [lce_opt, lt_slack, theta0]
%        num_nodes: direct collocation nodes (data nodes) in trajectory
%                   optimization
%        d:
%
% By: Huawei Wang
% Date: August 2, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    % muscle parameters
    lce_opt0 = par_mus(1:M);
    lt_slack0 = par_mus(M+1:2*M);
    theta0 = par_mus(2*M+1:3*M);
    Fmax0 = par_mus(3*M+1:4*M);  % nromalized by mass

    % reflex control gains
    par_rf_fse = x(end - nPar_rf + 1:end - nPar_rf + tPhase*M*M);
    par_rf_lce = x(end - nPar_rf + tPhase*M*M + 1:end - nPar_rf + 2*tPhase*M*M);
    par_rf = x(end - nPar_rf + 1:end - nPar_rf + 2*tPhase*M*M);
    % par_rf_res = x(end-3*M+1:end);
    
    
    % fit of muscle activations and joint torques
    sum_act = 0;
    sum_tor = 0;
    % sum_act_smo = 0;
    
    for t = 1:length(N)
        sta_st = sum(N(1:t-1))*M*S;
        mea_st = sum(N(1:t-1));
        for n = 1:N(t)
            
            sta_st_n = (n-1)*M*S;
            mea_st_n = mea_st + n;
            
            % extract optimized states and measured data
            act_opt = x(sta_st + sta_st_n + 4*M + 1:sta_st + sta_st_n + 5*M);
            lce = x(sta_st + sta_st_n + 2*M + 1: sta_st + sta_st_n + 3*M);
            lmt_mea = lmt(mea_st_n, :);
            act_mea = act_m(mea_st_n, :);
            d_mea = d(mea_st_n, :);
            tor_mea = torque_m(mea_st_n, :);
            
            % calculate muscle force and joint torques
            Fse = tendenForce_Groote_PRO(lmt_mea, lce, lce_opt0, lt_slack0, theta0);
            tor_opt = Fse.*Fmax0*d_mea';
            
            % calculate the sum of square fit
            sum_tor = sum_tor + sum((tor_opt - tor_mea).^2);
            sum_act = sum_act + sum((act_opt - act_mea).^2);
            
            % if n < N(t)
            %     act_opt_nt = x(sta_st + sta_st_n + M*S + 4*M + 1: ...
            %                   sta_st + sta_st_n + M*S + 5*M);
            %    sum_act_smo = sum_act_smo + sum((act_opt_nt - act_opt).^2);
            % end
        end
    end
    sum_tor = sum_tor/(sum(N)*M);
    sum_act = sum_act/(sum(N)*M);
    % sum_act_smo = sum_act_smo/(sum(N)*M);
    
    %% sparsity of the reflex control gains
    % since the fse and lce inputs has different amplitude and change
    % ranges, optimization will only prefer on type of reflexes (lce) when
    % using sum(sqrt(|gains|)).
    
    % One way to solve this is to make sure each gain that larger than a
    % specific number (for example 1e-3) have similar weights on the
    % objective function.  considering the gains' range is roughly [-1, 1].
    % Then we can use  g_nor = (exp(mul*g) - exp(-mul*g))/(exp(mul*g) + exp(mul*g))
    % where mul is the multiplier to change the threshold of the changes.
    % eventually, the g_nor should be a function of g in the following form
    %                                  ^  g_nor
    %                                1 |  ****************
    %                                  | *                 
    %                                  |*                      
    %                  ----------------*-------------->  g
    %                                 *|
    %                                * |
    %                               *  |
    %                  *************   |  -1
    
    
    mul = 10;
    g_nor = atan(mul*par_rf);
    % g_nor = (exp(mul*par_rf) - exp(-mul*par_rf))./(exp(mul*par_rf) + exp(-mul*par_rf));
    
    smooth_delta = 1e-4;
    % the sparsity term of reflex controllers
    sum_spar = sum((g_nor.^2 + smooth_delta).^0.25);
    
    
%     % the sparsity term of reflex controllers
%     smooth_delta = 1e-4;
%     sum_spar = sum((par_rf_fse.^2 + smooth_delta).^0.25) + ...
%                sum((par_rf_lce.^2 + smooth_delta).^0.25);
    
    
    % the smoothness term of reflex controllers
    sum_smo = 0;
    for p = 1:tPhase-1
        
        par_rf_fse_p1 = par_rf_fse((p - 1)*M*M + 1:p*M*M);
        par_rf_fse_p2 = par_rf_fse(p*M*M + 1:(p + 1)*M*M);
        
        par_rf_lce_p1 = par_rf_lce((p - 1)*M*M + 1:p*M*M);
        par_rf_lce_p2 = par_rf_lce(p*M*M + 1:(p + 1)*M*M);
        
        sum_smo = sum_smo + sum((par_rf_fse_p2 - par_rf_fse_p1).^2 + (par_rf_lce_p2 - par_rf_lce_p1).^2);
        
    end
    
    % sum all the objective terms together
    obj = (W1*(sum_tor + 5*sum_act) + W2*sum_spar/(2*tPhase*M*M) + W3*sum_smo/(2*tPhase*M*M));
    
end