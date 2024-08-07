%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get walking/running phase according to the GRF and ankle joint motions
%
% By: Huawei Wang
% Date: April 05, 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [phase, fig, l1, l2] = getPhase(ankle_ik, fy)

    % first plot them, then get the mannual detection of the phase separation
    % points, finally generate the phase.

    goodFlag = 0;

    fprintf('Click the point to separate phases ...\n')
    fprintf('The number of phases are depends on the number of clicking points + 1 ...\n')
    fprintf('Click the green dot to finish the phase separation ...\n')

    mean_fy = mean(fy);
    phase_id = [];
    
    fig = figure();
    l1 = plot(ankle_ik, '-o', 'linewidth', 2);
    hold on
    l2 = plot(fy, '-', 'linewidth', 2);
    hold on
    plot(0, mean(fy), 'go');
    title('PhaseSeparation')
    hold on
    
    while ~goodFlag

        [x,y]=ginput(1);

        FyError = 1;

        if sqrt(x^2 + (y-mean_fy)^2) > FyError
            plot(x, y, 'b*')
            hold on
            plot([x, x], [min(min(ankle_ik), min(fy)), max(max(ankle_ik), max(fy))], ...
                'b-', 'linewidth', 2)
            hold on

            phase_id = [phase_id, x];

        else
            break
        end
    end
    phase = zeros(length(phase_id), 1);
    for phase_i = 1:length(phase_id)
        if phase_i == 1
            phase(1:round(phase_id(phase_i))) = phase_i;
        else
            phase(round(phase_id(phase_i-1)) + 1:round(phase_id(phase_i))) = phase_i;
        end
        
        phase(round(phase_id(end)) + 1:length(fy)) = length(phase_id) + 1;

    end

end