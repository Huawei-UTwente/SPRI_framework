%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot reflex gains as heat map
%
% By: Huawei Wang
% Date: 04-11-2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function hf = reflex_heat_map_std_dots(ref_gains_mean, ref_gains_std, M, mus_names, cphase, phase_occupied, phase_names)

    hf = figure();
    
    % set(gcf,'position',[10 10 600 450])

    prf_fse_matrix = reshape(ref_gains_mean(1:M*M*cphase), M, M, cphase);  % fse reflex gains
    prf_lce_matrix = reshape(ref_gains_mean(M*M*cphase+1:M*M*cphase*2), M, M, cphase);  % lce reflex gains
    
    prf_fse_matrix_std = reshape(ref_gains_std(1:M*M*cphase), M, M, cphase);  % fse reflex gains
    prf_lce_matrix_std = reshape(ref_gains_std(M*M*cphase+1:M*M*cphase*2), M, M, cphase);  % lce reflex gains

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

    hmuscle = zeros(1, M);
    % phase_name = ["E Stance", "M Stance", "L Stance", "E Swing", "L Swing"];

    ncol = 24;
    
    % plot the muscle force relfex gains
    fse_g = subplot(M, ncol, [0*ncol+1:1*ncol-M, 1*ncol+1:2*ncol-M]);
    box on
    for m = 1:M
        hmuscle(m) =  1 - (m*vblock + (m-1)*ratio_fill_blank*vblock + M/2*vmblock);
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
                    1 - (m*vblock + (m-1)*ratio_fill_blank*vblock + mm*vmblock),...
                    hphase(p+1)-hphase(p), vmblock], 'FaceColor', block_color,...
                    'EdgeColor', 'k', 'LineWidth', 0.5)
                hold on
                
                if prf_fse_matrix_std(m, mm, p) > 5
                    prf_fse_matrix_std(m, mm, p) = 5;
                end
                
                block_color_std = [1-prf_fse_matrix_std(m, mm, p)/MaxColor*2,...
                    1-prf_fse_matrix_std(m, mm, p)/MaxColor*2,...
                    1-prf_fse_matrix_std(m, mm, p)/MaxColor*2];
                
                rectangle('Position',[(hphase(p) + hphase(p+1))/2-vmblock/2, ...
                    1 - (m*vblock + (m-1)*ratio_fill_blank*vblock + (mm-0.1)*vmblock),...
                    vmblock/3, vmblock*0.8], 'FaceColor', block_color_std, 'LineStyle', 'none',...
                    'Curvature',[1 1])
                hold on
            end
        end
    end
    xlim([0, 1])
    set(fse_g, 'box','on','XTick',[],'YTick',[])
    %set(gca,'xtick', (hphase(1:end-1) + hphase(2:end))/2 ,'xticklabel',phase_name)
    set(gca,'ytick', flip(hmuscle), 'yticklabel', flip(mus_names), 'fontweight','normal')
    title('Muscle Force Reflex Gains', 'fontweight','bold')
    

    % plot the muscle fibre length relfex gains
    lce_g = subplot(M, ncol, [2*ncol+1:3*ncol-M, 3*ncol+1:4*ncol-M]);
    box on
    for m = 1:M
        hmuscle(m) = 1 - (m*vblock + (m-1)*ratio_fill_blank*vblock + M/2*vmblock);
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
                    1 - (m*vblock + (m-1)*ratio_fill_blank*vblock + mm*vmblock),...
                    hphase(p+1)-hphase(p), vmblock], 'FaceColor', block_color,...
                    'EdgeColor', 'k', 'LineWidth', 0.5)
                hold on
                
                if prf_lce_matrix_std(mm, m, p) > 5
                    prf_lce_matrix_std(mm, m, p) = 5;
                end
                
                block_color_std = [1-prf_lce_matrix_std(mm, m, p)/MaxColor*2,...
                    1-prf_lce_matrix_std(mm, m, p)/MaxColor*2,...
                    1-prf_lce_matrix_std(mm, m, p)/MaxColor*2];
                
                rectangle('Position',[(hphase(p) + hphase(p+1))/2-vmblock/2, ...
                    1 - (m*vblock + (m-1)*ratio_fill_blank*vblock + (mm-0.1)*vmblock),...
                    vmblock/3, vmblock*0.8], 'FaceColor', block_color_std, 'LineStyle', 'none',...
                    'Curvature',[1, 1])
                hold on
            end
        end
    end
    xlim([0, 1])
    % axis([0, 1, 0, 1])
    set(lce_g, 'box','on','XTick',[],'YTick',[])
    set(gca,'ytick', flip(hmuscle), 'yticklabel',flip(mus_names), 'fontweight','bold')
    set(gca,'xtick', (hphase(1:end-1) + hphase(2:end))/2 ,'xticklabel',...
        phase_names, 'fontweight','normal')
    title('Muscle Length Reflex Gains', 'fontweight','bold')

    % plot thresholds of force and lengths
    tfse = subplot(M, ncol, [ncol-3:ncol-2, 2*ncol-3:2*ncol-2]);
    for m = 1:M
        
        box_mean = 1 - (m*vblock + m*ratio_fill_blank*vblock) + ...
            vmblock*M*ref_gains_mean(M*M*cphase*2+m);
        box_std = vmblock*M*ref_gains_std(M*M*cphase*2+m);
        
        rectangle('Position',[0.02,...
                    1- (m*vblock + m*ratio_fill_blank*vblock),...
                    0.04, vmblock*M*ref_gains_mean(M*M*cphase*2+m)],...
                    'FaceColor', [0, 0, 0], 'LineStyle', 'none')
        hold on
        plot([0.08, 0.08], [1 - (m*vblock + (m-1)*ratio_fill_blank*vblock),...
            1- (m*vblock + m*ratio_fill_blank*vblock)], 'color', [0.5, 0.5, 0.5],...
            'linewidth', 0.5)
        hold on
        text(0.09, 1- (m*vblock + m*ratio_fill_blank*vblock - vmblock/2),...
            num2str(0), 'FontSize', 7)
        text(0.09, 1 - (m*vblock + (m-1)*ratio_fill_blank*vblock + vmblock/2),...
            num2str(1), 'FontSize', 7)
        
        errorbar(0.04, box_mean, box_std, 'r-')
        
%         block_color_std = [1-ref_gains_std(M*M*cphase*2+m),...
%                     1-ref_gains_std(M*M*cphase*2+m),...
%                     1-ref_gains_std(M*M*cphase*2+m)];
%                 
%         rectangle('Position',[0.03, ...
%             m*vblock + (m-1)*ratio_fill_blank*vblock + vmblock,...
%             0.02, vmblock], 'FaceColor', block_color_std, 'LineStyle', 'none',...
%             'Curvature',[1, 1])
        hold on
        
    end
    xlim([-0.0, 0.12])
    ylim([0, 1])
    set(tfse, 'box','off','XTickLabel',[],'XTick',[],...
        'YTickLabel',[],'YTick',[], 'xcolor',[1 1 1],'ycolor',[1 1 1])

    tlce = subplot(M, ncol, [3*ncol-3:3*ncol-2, 4*ncol-3:4*ncol-2]);
    for m = 1:M
        
        box_mean = 1 - (m*vblock + (m)*ratio_fill_blank*vblock) + ...
            vmblock*M*(ref_gains_mean(M*M*cphase*2+M+m) - 0.0);
        box_std = vmblock*M*ref_gains_std(M*M*cphase*2+M+m);

        rectangle('Position',[0.02,...
                    1 - (m*vblock + (m)*ratio_fill_blank*vblock),...
                    0.04, vmblock*M*(ref_gains_mean(M*M*cphase*2+M+m) - 0.0)],...
                    'FaceColor', [0.25, 0.25, 0.25],...
                    'LineStyle', 'none')
        hold on
        plot([0.08, 0.08], [1- (m*vblock + (m-1)*ratio_fill_blank*vblock),...
            1 - (m*vblock + m*ratio_fill_blank*vblock)], 'color', [0.5, 0.5, 0.5],...
            'linewidth', 0.5)
        hold on
        
        text(0.09, 1- (m*vblock + m*ratio_fill_blank*vblock - vmblock/2),...
            num2str(0.0), 'FontSize', 7)
        text(0.09, 1- (m*vblock + (m-1)*ratio_fill_blank*vblock + vmblock/2),...
            num2str(1.5), 'FontSize', 7)
        
        errorbar(0.04, box_mean, box_std, 'r-')
        
%         block_color_std = [1-ref_gains_std(M*M*cphase*2+M+m)*2,...
%                     1-ref_gains_std(M*M*cphase*2+M+m)*2,...
%                     1-ref_gains_std(M*M*cphase*2+M+m)*2];
%                 
%         rectangle('Position',[0.03, ...
%             m*vblock + (m-1)*ratio_fill_blank*vblock + vmblock,...
%             0.02, vmblock], 'FaceColor', block_color_std, 'LineStyle', 'none',...
%             'Curvature',[1, 1])
        hold on
        
    end
    xlim([-0.0, 0.12])
    ylim([0, 1.5])
    text(-0.03, -0.08, 'Threshold', 'FontSize', 8)
    set(tlce, 'box','off', 'XTickLabel', [], 'XTick', [],...
        'YTickLabel',[],'YTick',[], 'xcolor',[0 0 0],'ycolor',[1 1 1])
    

    cbar = subplot(M, ncol, [ncol-1, 2*ncol-1, 3*ncol-1, 4*ncol-1]);
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
    text(0.01, -0.95, num2str(-MaxColor), 'FontSize', 7, 'color', 'r')
    text(0.03, 0, num2str(0), 'FontSize', 7)
    text(0.02, 0.95, num2str(MaxColor), 'FontSize', 7, 'color', 'b')
    text(-0.05, 1.05, 'Mean', 'FontSize', 8)
    set(cbar, 'box','off','XTickLabel','mean','XTick',0.05,'YTickLabel',[],'YTick',[])
    set(gca,'Visible','off')
    
    cbar = subplot(M, ncol, [ncol, 2*ncol, 3*ncol, 4*ncol]);
    hc = 0:0.001:1;
    for h = 1:length(hc)

    block_color = [1 - hc(h),...
                   1 - hc(h), 1 - hc(h)];

    rectangle('Position',[0,...
                hc(h),...
                0.1, 0.001],...
                'FaceColor', block_color, 'LineStyle', 'none')
    hold on
    end
    xlim([0, 0.1])
    ylim([0, 1])
    text(0.12, 0.025, '0', 'FontSize', 7)
    %text(0.15, 0, num2str(0), 'FontSize', 8)
    text(0.12, 0.975, '5(1)', 'FontSize', 7)
    text(-0.01, -0.035, 'Std', 'FontSize', 8)
    set(cbar, 'box', 'off', 'XTickLabel', 'std', 'XTick', 0.05, 'YTickLabel', [], 'YTick', [])
    set(gca,'Visible','off')
    
end
