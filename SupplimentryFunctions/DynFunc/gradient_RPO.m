function grad = gradient_RPO(x, M, S, N, nPar_rf, lmt, par_mus,...
     tPhase, W1, W2, W3, w11, w12, w13, w14, torque_m, act_m, d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The callback for calculating the gradient
%
% By: Huawei Wang
% Date: August 2, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    grad = zeros(size(x));
    
    % muscle parameters
    lce_opt0 = par_mus(1:M);
    lt_slack0 = par_mus(M+1:2*M);
    theta0 = par_mus(2*M+1:3*M);
    Fmax0 = par_mus(3*M+1:4*M);  % nromalized by mass

    % reflex control gains
    par_rf_fse = x(end - nPar_rf + 1:end - nPar_rf + tPhase*M*M);
    par_rf_lce = x(end - nPar_rf + tPhase*M*M + 1:end - nPar_rf + 2*tPhase*M*M);
    % par_rf_res = x(end-3*M+1:end);
    
    % gradient of the fit objective term
    % fit of muscle activations and joint torques
    w1 = W1/(sum(N)*M);
    w1_1 = w1*w11;
    w1_2 = w1*w12;
    w1_3 = w1*w13;
    w1_4 = w1*w14;
    
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
            [Fse, dFse_dlce] = tendenForce_Groote_diff_RPO(lmt_mea, ...
                lce, lce_opt0, lt_slack0, theta0);
            
            tor_opt = Fse.*Fmax0*d_mea';
            
            grad(sta_st + sta_st_n + 2*M + 1: sta_st + sta_st_n + 3*M) = ...
                grad(sta_st + sta_st_n + 2*M + 1: sta_st + sta_st_n + 3*M) + ...
                2*w1_1*(tor_opt - tor_mea)*Fmax0.*d_mea.*dFse_dlce;
            
            % calculate gradients
            grad(sta_st + sta_st_n + 4*M + 1:sta_st + sta_st_n + 5*M) = ...
                grad(sta_st + sta_st_n + 4*M + 1:sta_st + sta_st_n + 5*M) ... 
                + 2*w1_2*(act_opt - act_mea);
            
%             activation smoothness gradients
            if n < N(t)

                act_opt_nt = x(sta_st + sta_st_n + M*S + 4*M + 1: ...
                            sta_st + sta_st_n + M*S + 5*M);

                grad(sta_st + sta_st_n + M*S + 4*M + 1:sta_st + sta_st_n + M*S + 5*M) = ...
                grad(sta_st + sta_st_n + M*S + 4*M + 1:sta_st + sta_st_n + M*S + 5*M) ...
                + 2*w1_3*(act_opt_nt - act_opt);
            
                grad(sta_st + sta_st_n + 4*M + 1:sta_st + sta_st_n + 5*M) = ...
                grad(sta_st + sta_st_n + 4*M + 1:sta_st + sta_st_n + 5*M) ...
                - 2*w1_3*(act_opt_nt - act_opt);
            
                lce_nt = x(sta_st + sta_st_n + M*S + 2*M + 1: sta_st + sta_st_n + M*S + 3*M);
                lmt_mea_nt = lmt(mea_st_n + 1, :);
                [Fse_nt, dFse_nt_dlce] = ...
                    tendenForce_Groote_diff_RPO(lmt_mea_nt, lce_nt, lce_opt0, lt_slack0, theta0);
                
                dsum_fse_dFse = -2*w1_4*(Fse_nt - Fse);
                dsum_fse_dFse_nt = 2*w1_4*(Fse_nt - Fse);
                
                % gradient of the muscle force changes
                grad(sta_st + sta_st_n + 2*M + 1: sta_st + sta_st_n + 3*M) = ...
                grad(sta_st + sta_st_n + 2*M + 1: sta_st + sta_st_n + 3*M) + ... 
                        dsum_fse_dFse.*dFse_dlce;

                % gradient of the muscle force changes
                grad(sta_st + sta_st_n + M*S + 2*M + 1: sta_st + sta_st_n + M*S + 3*M) = ...
                grad(sta_st + sta_st_n + M*S + 2*M + 1: sta_st + sta_st_n + M*S + 3*M) + ...
                        dsum_fse_dFse_nt.*dFse_nt_dlce;
                
            end
            
        end
    end
    
    
    % the sparsity term of reflex controllers
    w2 = W2/(2*tPhase*M*M);
    
    smooth_delta = 1e-4;
    
    grad(end - nPar_rf + 1:end - nPar_rf + tPhase*M*M) = ...
        w2*par_rf_fse./(2*(par_rf_fse.^2 + smooth_delta).^0.75);
    
    grad(end - nPar_rf + tPhase*M*M + 1:end - nPar_rf + 2*tPhase*M*M) = ...
        w2*par_rf_lce./(2*(par_rf_lce.^2 + smooth_delta).^0.75);
    
    
    % the smoothness term of reflex controllers
    w3 = W3/(2*tPhase*M*M);
    for p = 1:tPhase-1
        
        par_rf_fse_p1 = par_rf_fse((p - 1)*M*M + 1:p*M*M);
        par_rf_fse_p2 = par_rf_fse(p*M*M + 1:(p + 1)*M*M);
        
        par_rf_lce_p1 = par_rf_lce((p - 1)*M*M + 1:p*M*M);
        par_rf_lce_p2 = par_rf_lce(p*M*M + 1:(p + 1)*M*M);
        
        grad(end - nPar_rf + (p - 1)*M*M + 1:end - nPar_rf + p*M*M) = ...
        grad(end - nPar_rf + (p - 1)*M*M + 1:end - nPar_rf + p*M*M) - ...
        2*w3.*(par_rf_fse_p2 - par_rf_fse_p1);
    
        grad(end - nPar_rf + p*M*M + 1:end - nPar_rf + (p + 1)*M*M) = ...
        grad(end - nPar_rf + p*M*M + 1:end - nPar_rf + (p + 1)*M*M) + ...
        2*w3.*(par_rf_fse_p2 - par_rf_fse_p1);
    
        grad(end - nPar_rf + tPhase*M*M + (p - 1)*M*M + 1:end - nPar_rf + tPhase*M*M + p*M*M) = ...
        grad(end - nPar_rf + tPhase*M*M +  (p - 1)*M*M + 1:end - nPar_rf + tPhase*M*M + p*M*M) - ...
        2*w3.*(par_rf_lce_p2 - par_rf_lce_p1);
    
        grad(end - nPar_rf + tPhase*M*M + p*M*M + 1:end - nPar_rf + tPhase*M*M + (p + 1)*M*M) = ...
        grad(end - nPar_rf + tPhase*M*M + p*M*M + 1:end - nPar_rf + tPhase*M*M + (p + 1)*M*M) + ...
        2*w3.*(par_rf_lce_p2 - par_rf_lce_p1);
 
    end
    
end