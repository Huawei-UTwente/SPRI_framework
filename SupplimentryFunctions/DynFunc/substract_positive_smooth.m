    % smooth function of the positive curve
    function px = substract_positive_smooth(x, throshold, smooth_delta)

        dx = x - throshold;
        px = (sqrt(dx.^2 + smooth_delta) + dx)/2;
        
    end
