% smooth function of the positive curve
    function [px, dpx_dx, dpx_dthroshold] = substract_positive_smooth_diff(x, throshold, smooth_delta)

        dx = x - throshold;
        
        px = (sqrt(dx.^2 + smooth_delta) + dx)/2;

        dpx_dx = 1./(2*sqrt(dx.^2 + smooth_delta)).*dx + 1/2;
        dpx_dthroshold = -1./(2*sqrt(dx.^2 + smooth_delta)).*dx - 1/2;

    end