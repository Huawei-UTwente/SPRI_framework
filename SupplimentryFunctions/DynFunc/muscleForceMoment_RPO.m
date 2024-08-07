%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is to generate muscle forces and joint moments based on the
% optimized muscle parameters.
%
% By: Huawei Wang
% Date: August 8, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mom, force] = muscleForceMoment_RPO(x, lmt, d, M, N, S, J, par_mus)

    lce_opt0 = par_mus(1:M);
    lt_slack0 = par_mus(M + 1:2*M);
    theta0 = par_mus(2*M + 1:3*M);
    Fmax0 = par_mus(3*M + 1:4*M);
    
    mom = zeros(sum(N), J);
    force = zeros(sum(N), M);
    
    for t = 1:length(N)
        sta_st = sum(N(1:t-1))*M*S;
        mea_st = sum(N(1:t-1));
        for n = 1:N(t)
            sta_st_n = (n-1)*M*S;
            
            lce_opt = x(sta_st + sta_st_n + 2*M + 1: sta_st + sta_st_n + 3*M);
            lmt_mea = lmt(mea_st+n, :);
            d_mea = d(mea_st+n, :);
            
            % calculate muscle force and joint torques
            force(mea_st+n, :) = tendenForce_Groote_RPO(lmt_mea, lce_opt, lce_opt0, lt_slack0, theta0);
            mom(mea_st+n, :) = force(mea_st+n, :).*Fmax0*d_mea';

        end
    end
end