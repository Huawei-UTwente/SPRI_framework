function [cos_theta, dcos_theta_dlce] ...
                = pennationAngSmooth_diff_RPO(lce, lce_opt, theta0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the pennation angle cos_theta and corresponding derivative,
% assuming the muscle parameters are determined already.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % smooth function parameter that determine the smoothness of the
    % changes
    smooth_delta = 1e-8;

    % angle of the contraction element
    h = lce_opt.*sin(theta0);
    
    % calcualte the ce length in the parallel direction
    y = lce.^2 - h.^2; 
    
    % calcualte it's derivatives
    dy_dlce = 2*lce;
    
    dy = (sqrt(y.^2 + smooth_delta) + y)/2;
    
    ddy_dy = y./sqrt(y.^2 + smooth_delta)/2 + 1/2; 
    
    ddy_dlce = ddy_dy.*dy_dlce;
    
    cos_theta = sqrt(dy)./lce;  % calcualte cos(theta) at current lce length
    
    % calculate its derivatives    
    dcos_theta_ddy = 1./(2*sqrt(dy).*lce);
    
    dcos_theta_dlce = dcos_theta_ddy.*ddy_dlce - sqrt(dy)./lce.^2;
   
end