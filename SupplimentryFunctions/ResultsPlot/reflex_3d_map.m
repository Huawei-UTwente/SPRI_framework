%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot reflex gains as heat map
%
% By: Huawei Wang
% Date: 04-11-2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function hf = reflex_3d_map(ref_gains, M, cphase, phase_occupied)

    hf = figure();
    box on
    grid on
    set(gcf,'renderer','zbuffer')

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
    set(gcf,'renderer','zbuffer')
    for m = 1:M
        hmuscle(m) = m*vblock + (m-1)*ratio_fill_blank*vblock + M/2*vmblock;
        for mm = 1:M
            for p = 1:cphase
                if prf_fse_matrix(m, mm, p) > 0
                    block_color = [1, 1- prf_fse_matrix(m, mm, p)/MaxColor,...
                                      1- prf_fse_matrix(m, mm, p)/MaxColor];
                else
                    block_color = [1 + prf_fse_matrix(m, mm, p)/MaxColor,...
                                   1 + prf_fse_matrix(m, mm, p)/MaxColor, 1];
                end
                
%                 p1 = [hphase(p), m*vblock + (m-1)*ratio_fill_blank*vblock + (mm-1)*vmblock, prf_fse_matrix(m, mm, p)];
%                 p2 = [hphase(p), m*vblock + (m-1)*ratio_fill_blank*vblock + mm*vmblock, prf_fse_matrix(m, mm, p)];
%                 p3 = [hphase(p+1), m*vblock + (m-1)*ratio_fill_blank*vblock + (mm-1)*vmblock, prf_fse_matrix(m, mm, p)];
%                 p4 = [hphase(p+1), m*vblock + (m-1)*ratio_fill_blank*vblock + mm*vmblock, prf_fse_matrix(m, mm, p)];
                
                x = (hphase(p) + hphase(p+1))/2;
                y = m*vblock + (m-1)*ratio_fill_blank*vblock + (mm-0.5)*vmblock;
                z = prf_fse_matrix(m, mm, p);
                
                width_x = hphase(p+1)-hphase(p);
                width_y = vmblock;
                
                edge_alpha = 0;
                
                h(1)=patch([-width_x -width_x width_x width_x]+x,[-width_y width_y width_y -width_y]+y, [0 0 0 0], 'b', 'EdgeAlpha', 1);
                h(2)=patch(width_x.*[-1 -1 1 1]+x, width_y.*[-1 -1 -1 -1]+y, z.*[0 1 1 0], 'b', 'EdgeAlpha', edge_alpha);
                h(3)=patch(width_x.*[-1 -1 -1 -1]+x, width_y.*[-1 -1 1 1]+y, z.*[0 1 1 0], 'b', 'EdgeAlpha', edge_alpha);
                h(4)=patch([-width_x -width_x width_x width_x]+x, [-width_y width_y width_y -width_y]+y, [z z z z], 'b', 'EdgeAlpha', edge_alpha);
                h(5)=patch(width_x.*[-1 -1 1 1]+x, width_y.*[1 1 1 1]+y, z.*[0 1 1 0], 'b', 'EdgeAlpha', edge_alpha);
                h(6)=patch(width_x.*[1 1 1 1]+x, width_y.*[-1 -1 1 1]+y, z.*[0 1 1 0], 'b', 'EdgeAlpha', edge_alpha);
                set(h, 'facecolor', block_color)
                
                % scatterbar3(x, y, z, hphase(p+1)-hphase(p), vmblock, block_color)

                %fill3(x, y, z, 'grouped', block_color)

%                 rectangle('Position',[hphase(p),...
%                     m*vblock + (m-1)*ratio_fill_blank*vblock + (mm-1)*vmblock,...
%                     hphase(p+1)-hphase(p), vmblock], 'FaceColor', block_color,...
%                     'EdgeColor', 'k', 'LineWidth', 0.5)
                hold on

            end
        end
    end
    axis([0 1 0 1 -10, 10])
    caxis([-10, 10])
    view(3)
    
    
%     set(fse_g, 'box','on','XTick',[],'YTick',[])
%     %set(gca,'xtick', (hphase(1:end-1) + hphase(2:end))/2 ,'xticklabel',phase_name)
%     set(gca,'ytick', hmuscle, 'yticklabel',muscle_name, 'fontweight','normal')
    title('Muscle Force Reflex Gains', 'fontweight','bold')
    

    % plot the muscle fibre length relfex gains
    lce_g = subplot(4, ncol, [2*ncol+1:3*ncol-3, 3*ncol+1:4*ncol-3]);
    set(gcf,'renderer','zbuffer')
    box on
    for m = 1:M
        hmuscle(m) = m*vblock + (m-1)*ratio_fill_blank*vblock + M/2*vmblock;
        for mm = 1:M
            for p = 1:cphase
                if prf_lce_matrix(m, mm, p) > 0
                    block_color = [1, 1- prf_lce_matrix(m, mm, p)/MaxColor,...
                                      1- prf_lce_matrix(m, mm, p)/MaxColor];
                else
                    block_color = [1 + prf_lce_matrix(m, mm, p)/MaxColor,...
                                   1 + prf_lce_matrix(m, mm, p)/MaxColor, 1];
                end
                
%                 p1 = [hphase(p), m*vblock + (m-1)*ratio_fill_blank*vblock + (mm-1)*vmblock, prf_lce_matrix(m, mm, p)];
%                 p2 = [hphase(p), m*vblock + (m-1)*ratio_fill_blank*vblock + mm*vmblock, prf_lce_matrix(m, mm, p)];
%                 p3 = [hphase(p+1), m*vblock + (m-1)*ratio_fill_blank*vblock + (mm-1)*vmblock, prf_lce_matrix(m, mm, p)];
%                 p4 = [hphase(p+1), m*vblock + (m-1)*ratio_fill_blank*vblock + mm*vmblock, prf_lce_matrix(m, mm, p)];
%                 
%                 x = [p1(1) p2(1) p3(1) p4(1)];
%                 y = [p1(2) p2(2) p3(2) p4(2)];
%                 z = [p1(3) p2(3) p3(3) p4(3)];
                
                x = (hphase(p) + hphase(p+1))/2;
                y = m*vblock + (m-1)*ratio_fill_blank*vblock + (mm-0.5)*vmblock;
                z = prf_lce_matrix(m, mm, p);
                
                                width_x = hphase(p+1)-hphase(p);
                width_y = vmblock;
                
                edge_alpha = 0;
                
                h(1)=patch([-width_x -width_x width_x width_x]+x,[-width_y width_y width_y -width_y]+y, [0 0 0 0], 'b', 'EdgeAlpha', 1);
                h(2)=patch(width_x.*[-1 -1 1 1]+x, width_y.*[-1 -1 -1 -1]+y, z.*[0 1 1 0], 'b', 'EdgeAlpha', edge_alpha);
                h(3)=patch(width_x.*[-1 -1 -1 -1]+x, width_y.*[-1 -1 1 1]+y, z.*[0 1 1 0], 'b', 'EdgeAlpha', edge_alpha);
                h(4)=patch([-width_x -width_x width_x width_x]+x, [-width_y width_y width_y -width_y]+y, [z z z z], 'b', 'EdgeAlpha', edge_alpha);
                h(5)=patch(width_x.*[-1 -1 1 1]+x, width_y.*[1 1 1 1]+y, z.*[0 1 1 0], 'b', 'EdgeAlpha', edge_alpha);
                h(6)=patch(width_x.*[1 1 1 1]+x, width_y.*[-1 -1 1 1]+y, z.*[0 1 1 0], 'b', 'EdgeAlpha', edge_alpha);
                set(h, 'facecolor', block_color)
                
                % scatterbar3(x, y, z, hphase(p+1)-hphase(p), vmblock, block_color)
                
                % fill3(x, y, z, 'grouped', block_color)

%                 rectangle('Position',[hphase(p),...
%                     m*vblock + (m-1)*ratio_fill_blank*vblock + (mm-1)*vmblock,...
%                     hphase(p+1)-hphase(p), vmblock], 'FaceColor', block_color,...
%                     'EdgeColor', 'k', 'LineWidth', 0.5)
                hold on

            end
        end
    end
    axis([0 1 0 1 -10, 10])
    caxis([-10, 10])
    view(3)
    
%     set(lce_g, 'box','on','XTick',[],'YTick',[])
%     set(gca,'ytick', hmuscle, 'yticklabel',muscle_name, 'fontweight','bold')
%     set(gca,'xtick', (hphase(1:end-1) + hphase(2:end))/2 ,'xticklabel',...
%         phase_name, 'fontweight','normal')
    title('Muscle Length Reflex Gains', 'fontweight','bold')

%     tfse = subplot(4, ncol, [ncol-2:ncol-1, 2*ncol-2:2*ncol-1]);
%     for m = 1:M
% 
%         rectangle('Position',[0.02,...
%                     m*vblock + (m-1)*ratio_fill_blank*vblock,...
%                     0.04, vmblock*M*ref_gains(M*M*5*2+m)],...
%                     'FaceColor', [0, 0, 0], 'LineStyle', 'none')
%         hold on
%         plot([0.08, 0.08], [m*vblock + (m-1)*ratio_fill_blank*vblock,...
%             m*vblock + m*ratio_fill_blank*vblock], 'color', [0.5, 0.5, 0.5],...
%             'linewidth', 0.5)
%         hold on
%         text(0.09, m*vblock + (m-1)*ratio_fill_blank*vblock + vmblock/2,...
%             num2str(0), 'FontSize', 7)
%         text(0.09, m*vblock + m*ratio_fill_blank*vblock - vmblock/2,...
%             num2str(1), 'FontSize', 7)
%         
%     end
%     xlim([-0.0, 0.12])
%     ylim([0, 1])
%     set(tfse, 'box','off','XTickLabel',[],'XTick',[],...
%         'YTickLabel',[],'YTick',[], 'xcolor',[1 1 1],'ycolor',[1 1 1])
% 
%     tlce = subplot(4, ncol, [3*ncol-2:3*ncol-1, 4*ncol-2:4*ncol-1]);
%     for m = 1:M
% 
%         rectangle('Position',[0.02,...
%                     m*vblock + (m-1)*ratio_fill_blank*vblock,...
%                     0.04, vmblock*M*ref_gains(M*M*5*2+M+m)],...
%                     'FaceColor', [0.25, 0.25, 0.25],...
%                     'LineStyle', 'none')
%         hold on
%         plot([0.08, 0.08], [m*vblock + (m-1)*ratio_fill_blank*vblock,...
%             m*vblock + m*ratio_fill_blank*vblock], 'color', [0.5, 0.5, 0.5],...
%             'linewidth', 0.5)
%         hold on
%         
%         text(0.09, m*vblock + (m-1)*ratio_fill_blank*vblock + vmblock/2,...
%             num2str(0), 'FontSize', 7)
%         text(0.09, m*vblock + m*ratio_fill_blank*vblock - vmblock/2,...
%             num2str(1), 'FontSize', 7)
%         
%     end
%     xlim([-0.0, 0.12])
%     ylim([0, 1])
%     set(tlce, 'box','off', 'XTickLabel', 'Threshold', 'XTick', 0.06,...
%         'YTickLabel',[],'YTick',[], 'xcolor',[0 0 0],'ycolor',[1 1 1])
%     
% 
%     cbar = subplot(4, ncol, [ncol, 2*ncol, 3*ncol, 4*ncol]);
%     hc = -1:0.001:1;
%     for h = 1:length(hc)
% 
%         if hc(h) > 0
%             block_color = [1, 1- hc(h),...
%                               1- hc(h)];
%         else
%             block_color = [1 + hc(h),...
%                            1 + hc(h), 1];
%         end
% 
%         rectangle('Position',[0,...
%                     hc(h),...
%                     0.1, 0.001],...
%                     'FaceColor', block_color, 'LineStyle', 'none')
%         hold on
%     end
%     xlim([0, 0.1])
%     ylim([-1, 1])
%     text(0.15, -1, num2str(-MaxColor), 'FontSize', 8)
%     text(0.15, 0, num2str(0), 'FontSize', 8)
%     text(0.15, 1, num2str(MaxColor), 'FontSize', 8)
%     set(cbar, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
%     set(gca,'Visible','off')
    
end
