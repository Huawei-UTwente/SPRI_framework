function obj =  objective_RPO_Fit(x, M, S, N, lmt, par_mus, W1, ...
                                  w11, w12, w13, w14, torque_m, act_m, d)
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

    % fit of muscle activations and joint torques
    sum_act = 0;
    sum_tor = 0;
    sum_act_smo = 0;
    sum_fse_smo = 0;
    
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
            Fse = tendenForce_Groote_RPO(lmt_mea, lce, lce_opt0, lt_slack0, theta0);
            tor_opt = Fse.*Fmax0*d_mea';
            
            % calculate the sum of square fit
            sum_tor = sum_tor + sum(sum((tor_opt - tor_mea).^2));
            sum_act = sum_act + sum(sum((act_opt - act_mea).^2));
            
            if n < N(t)
                act_opt_nt = x(sta_st + sta_st_n + M*S + 4*M + 1: ...
                              sta_st + sta_st_n + M*S + 5*M);
               sum_act_smo = sum_act_smo + sum((act_opt_nt - act_opt).^2);
               
               lce_nt = x(sta_st + sta_st_n + M*S + 2*M + 1: sta_st + sta_st_n + M*S + 3*M);
               lmt_mea_nt = lmt(mea_st_n + 1, :);
               
               Fse_nt = tendenForce_Groote_RPO(lmt_mea_nt, lce_nt, lce_opt0, lt_slack0, theta0);
               sum_fse_smo = sum_fse_smo + sum((Fse_nt - Fse).^2);
               
            end
        end
    end

    % sum all the objective terms together
    obj = W1/(sum(N)*M)*(w11*sum_tor + w12*sum_act + w13*sum_act_smo + w14*sum_fse_smo);
    
end