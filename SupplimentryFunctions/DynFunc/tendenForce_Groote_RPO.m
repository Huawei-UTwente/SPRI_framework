function Fse_S = tendenForce_Groote_RPO(lmt, lce, lce_opt, lt_slack, theta0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% This is the tenden force calculation
%
% By: Huawei Wang
% Date: August 1, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % some default parameters, from Groote 2016 paper
    c1 = 0.2;
    c2 = 1;    %0.995;
    c3 = 0.2;  % 0.25
    % kT = 35;
    
    kT = 17.5; %10% strain
    
    delta = 1e-4; % smooth parameters 
    
    % pennation angle calculation and its derivatives
    cos_theta = pennationAngSmooth(lce, lce_opt, theta0);
        
    % calcualte tenden lenght
    lt = lmt - lce.*cos_theta;

    % calculate normalized tendon length 
    lt_nor = lt./lt_slack;
    
    % calculate tendon force
    Fse = c1*exp(kT.*(lt_nor - c2)) - c3;
    
    % smooth the Fse & removing negative Fse values
    Fse_S = (sqrt(Fse.^2 + delta) + Fse)/2;
    
    % else
    % Fse_S = Fse;
    
end