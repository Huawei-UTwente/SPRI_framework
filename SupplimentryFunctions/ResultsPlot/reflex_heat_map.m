%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot reflex gains as heat map
%
% By: Huawei Wang
% Date: 04-11-2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function hf = reflex_heat_map(ref_gains, M, cphase, phase_occupied)

    hf = figure();

    prf_fse_matrix = reshape(ref_gains(1:M*M*cphase), M, M, cphase);  % fse reflex gains
    prf_lce_matrix = reshape(ref_gains(M*M*cphase+1:M*M*cphase*2), M, M, cphase);  % lce reflex gains

    ratio_fill_blank = 2.5;  % ratio between the filled reflex blocks and the blank area
    vblock = 1/(M*ratio_fill_blank + M + 1);  % vertical height of the blank area
    vmblock = vblock*ratio_fill_blank/M;  % vertical height of each muscle feedback row

    hphase = zeros(1, length(phase_occupied)+1);
    for l = 1:length(phase_occupied)
        hphase(l+1) = sum(phase_occupied(1:l))/sum(phase_occupied);
    end
    hphase(1) = 1e-3;
    hphase(end) = 1 - 1e-3;

    MaxColor = 10;  % the maximum gain for the color

    muscle_name = ["LGas", "MGas", "Sol", "TA"];
    hmuscle = zeros(1, M);
    phase_name = ["E Stance", "M Stance", "L Stance", "E Swing", "L Swing"];

    ncol = 24;
    
    % plot the muscle force relfex gains
    fse_g = subplot(4, ncol, [0*ncol+1:1*ncol-3, 1*ncol+1:2*ncol-3]);
    box on
    for m = 1:M
        hmuscle(m) = m*vblock + (m-1)*ratio_fill_blank*vblock + M/2*vmblock;
        for mm = 1:M
            for p = 1:cphase
                if prf_fse_matrix(mm, m, p) > 0
                    block_color = [1, 1- prf_fse_matrix(mm, m, p)/MaxColor,...
                                      1- prf_fse_matrix(mm, m, p)/MaxColor];
                else
                    block_color = [1 + prf_fse_matrix(mm, m, p)/MaxColor,...
                                   1 + prf_fse_matrix(mm, m, p)/MaxColor, 1];
                end

                rectangle('Position',[hphase(p),...
                    m*vblock + (m-1)*ratio_fill_blank*vblock + (mm-1)*vmblock,...
                    hphase(p+1)-hphase(p), vmblock], 'FaceColor', block_color,...
                    'EdgeColor', 'k', 'LineWidth', 0.5)
                hold on

            end
        end
    end
    set(fse_g, 'box','on','XTick',[],'YTick',[])
    %set(gca,'xtick', (hphase(1:end-1) + hphase(2:end))/2 ,'xticklabel',phase_name)
    set(gca,'ytick', hmuscle, 'yticklabel',muscle_name, 'fontweight','normal')
    title('Muscle Force Reflex Gains', 'fontweight','bold')
    

    % plot the muscle fibre length relfex gains
    lce_g = subplot(4, ncol, [2*ncol+1:3*ncol-3, 3*ncol+1:4*ncol-3]);
    box on
    for m = 1:M
        hmuscle(m) = m*vblock + (m-1)*ratio_fill_blank*vblock + M/2*vmblock;
        for mm = 1:M
            for p = 1:cphase
                if prf_lce_matrix(mm, m, p) > 0
                    block_color = [1, 1- prf_lce_matrix(mm, m, p)/MaxColor,...
                                      1- prf_lce_matrix(mm, m, p)/MaxColor];
                else
                    block_color = [1 + prf_lce_matrix(mm, m, p)/MaxColor,...
                                   1 + prf_lce_matrix(mm, m, p)/MaxColor, 1];
                end

                rectangle('Position',[hphase(p),...
                    m*vblock + (m-1)*ratio_fill_blank*vblock + (mm-1)*vmblock,...
                    hphase(p+1)-hphase(p), vmblock], 'FaceColor', block_color,...
                    'EdgeColor', 'k', 'LineWidth', 0.5)
                hold on

            end
        end
    end
    set(lce_g, 'box','on','XTick',[],'YTick',[])
    set(gca,'ytick', hmuscle, 'yticklabel',muscle_name, 'fontweight','bold')
    set(gca,'xtick', (hphase(1:end-1) + hphase(2:end))/2 ,'xticklabel',...
        phase_name, 'fontweight','normal')
    title('Muscle Length Reflex Gains', 'fontweight','bold')

    tfse = subplot(4, ncol, [ncol-2:ncol-1, 2*ncol-2:2*ncol-1]);
    for m = 1:M

        rectangle('Position',[0.02,...
                    m*vblock + (m-1)*ratio_fill_blank*vblock,...
                    0.04, vmblock*M*ref_gains(M*M*5*2+m)],...
                    'FaceColor', [0, 0, 0], 'LineStyle', 'none')
        hold on
        plot([0.08, 0.08], [m*vblock + (m-1)*ratio_fill_blank*vblock,...
            m*vblock + m*ratio_fill_blank*vblock], 'color', [0.5, 0.5, 0.5],...
            'linewidth', 0.5)
        hold on
        text(0.09, m*vblock + (m-1)*ratio_fill_blank*vblock + vmblock/2,...
            num2str(0), 'FontSize', 7)
        text(0.09, m*vblock + m*ratio_fill_blank*vblock - vmblock/2,...
            num2str(1), 'FontSize', 7)
        
    end
    xlim([-0.0, 0.12])
    ylim([0, 1])
    set(tfse, 'box','off','XTickLabel',[],'XTick',[],...
        'YTickLabel',[],'YTick',[], 'xcolor',[1 1 1],'ycolor',[1 1 1])

    tlce = subplot(4, ncol, [3*ncol-2:3*ncol-1, 4*ncol-2:4*ncol-1]);
    for m = 1:M

        rectangle('Position',[0.02,...
                    m*vblock + (m-1)*ratio_fill_blank*vblock,...
                    0.04, vmblock*M*ref_gains(M*M*5*2+M+m)],...
                    'FaceColor', [0.25, 0.25, 0.25],...
                    'LineStyle', 'none')
        hold on
        plot([0.08, 0.08], [m*vblock + (m-1)*ratio_fill_blank*vblock,...
            m*vblock + m*ratio_fill_blank*vblock], 'color', [0.5, 0.5, 0.5],...
            'linewidth', 0.5)
        hold on
        
        text(0.09, m*vblock + (m-1)*ratio_fill_blank*vblock + vmblock/2,...
            num2str(0), 'FontSize', 7)
        text(0.09, m*vblock + m*ratio_fill_blank*vblock - vmblock/2,...
            num2str(1), 'FontSize', 7)
        
    end
    xlim([-0.0, 0.12])
    ylim([0, 1])
    text(-0.03, -0.08, 'Threshold', 'FontSize', 8)
    set(tlce, 'box','off', 'XTickLabel', [], 'XTick', [],...
        'YTickLabel',[],'YTick',[], 'xcolor',[0 0 0],'ycolor',[1 1 1])
    

    cbar = subplot(4, ncol, [ncol, 2*ncol, 3*ncol, 4*ncol]);
    hc = -1:0.001:1;
    for h = 1:length(hc)

        if hc(h) > 0
            block_color = [1, 1- hc(h),...
                              1- hc(h)];
        else
            block_color = [1 + hc(h),...
                           1 + hc(h), 1];
        end

        rectangle('Position',[0,...
                    hc(h),...
                    0.1, 0.001],...
                    'FaceColor', block_color, 'LineStyle', 'none')
        hold on
    end
    xlim([0, 0.1])
    ylim([-1, 1])
    text(0.15, -1, num2str(-MaxColor), 'FontSize', 8)
    text(0.15, 0, num2str(0), 'FontSize', 8)
    text(0.15, 1, num2str(MaxColor), 'FontSize', 8)
    
    set(cbar, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
    set(gca,'Visible','off')
    
end
