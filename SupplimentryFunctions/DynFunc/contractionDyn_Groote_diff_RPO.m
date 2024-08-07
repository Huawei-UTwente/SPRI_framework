
function [f, df_da, df_dlce, df_ddlce]...
         = contractionDyn_Groote_diff_RPO(lmt, a, lce, dlce, lce_opt, lt_slack, theta0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the contraction dynamics of Hill's muscle model
%
% By: Huawei Wang
% Date: August 1, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % some default parameters, from Groote 2016 paper
    b11 = 0.815;
    b21 = 1.055;
    b31 = 0.162;
    b41 = 0.063;
    b12 = 0.433;
    b22 = 0.717;
    b32 = -0.030;
    b42 = 0.200;
    b13 = 0.100;
    b23 = 1.000;
    b33 = 0.354;
    b43 = 0.000;
    kpe = 4.0;
    e0 = 0.6;
    d1 = -0.318;
    d2 = -8.149;
    d3 = -0.374;
    d4 = 0.886;
    
    delta = 1e-4;
    
    % normalize muscle fiber lengths
    lce_nor = lce./lce_opt;
    
    dlce_nor_dlce = 1./lce_opt;
    
    % active force length relationship and its differentiatons
    fce = gaussianFunctionAct(lce_nor, b11, b21, b31, b41) + ...
          gaussianFunctionAct(lce_nor, b12, b22, b32, b42) + ...
          gaussianFunctionAct(lce_nor, b13, b23, b33, b43);
      
    dfce_dlce = gaussianFunctionAct_diff(lce_nor, b11, b21, b31, b41).*dlce_nor_dlce + ...
        gaussianFunctionAct_diff(lce_nor, b12, b22, b32, b42).*dlce_nor_dlce + ...
        gaussianFunctionAct_diff(lce_nor, b13, b23, b33, b43).*dlce_nor_dlce;
        
   % passive force length relationship and its differentiations
   fpee1 = exp(kpe*(lce_nor - 1)/e0);
   dfpee1_dlce = fpee1*(kpe/e0).*dlce_nor_dlce;
   
   fpee2 = (fpee1 - 1)./(exp(kpe) - 1);
   
   dfpee2_dlce = dfpee1_dlce/(exp(kpe) - 1);
   
   %% smooth fpee values & larger than 0
   fpee = (sqrt(fpee2.^2 + delta) + fpee2)/2;
   
   dfpee_dfpee2 = 0.5*fpee2./sqrt(fpee2.^2 + delta) + 0.5;

   dfpee_dlce = dfpee_dfpee2.*dfpee2_dlce;
   
%    % else
%    fpee = fpee2;
%    dfpee_dlce = dfpee2_dlce;
   
   
   %% force velocity relationship and its differentiations.
   dlceMax = 10*lce_opt;
   
   dlce_nor = dlce./dlceMax;
   
   ddlce_nor_ddlce = 1./dlceMax;
   
   fv_logfun = (d2*dlce_nor + d3) + sqrt(((d2*dlce_nor + d3).^2) + 1);
   fv_log = log(fv_logfun);
   dfv_log_ddlce = 1./fv_logfun.*(d2.*ddlce_nor_ddlce +...
       1./sqrt(((d2*dlce_nor + d3).^2) + 1)*d2.*(d2*dlce_nor + d3).*ddlce_nor_ddlce);
   
   fv = d1*fv_log + d4;
   dfv_ddlce = d1*dfv_log_ddlce;
   
   [cos_theta, dcos_theta_dlce] ...
        = pennationAngSmooth_diff_RPO(lce, lce_opt, theta0);
    
    % calcualte the force of the contraction element and PEE together
    Fce = (a.*fce.*fv + fpee).*cos_theta;
    
    % calculate overall derivatives
    dFce_da = fce.*fv.*cos_theta;
    dFce_dlce = (a.*fv.*dfce_dlce + dfpee_dlce).*cos_theta... 
                 + (a.*fce.*fv + fpee).*dcos_theta_dlce;  
    dFce_ddlce = (a.*fce.*dfv_ddlce).*cos_theta;

    [Fse, dFse_dlce] = ...
        tendenForce_Groote_diff_RPO(lmt, lce, lce_opt, lt_slack, theta0);
    
    % force of the contraction element should equal to the force of the tense unit
    f = (Fce - Fse);
    
    % derivative of the muslce force constraints
    df_da = dFce_da;
    df_dlce = (dFce_dlce - dFse_dlce);    
    df_ddlce = dFce_ddlce;

end 