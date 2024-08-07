%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Validation of the identified reflex controllers and muscle parameters
% 
% By: Huawei Wang
% Date: April 13th 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function obj_val = fitnessHeatMap(homeDataPath, musOptPath, homeOrgPath, ...
    homeValPath, subj, oriTrialNames, valTrialNames, Phase_list, W2_list, homeSavePath)

    % initlize parameters
    T = length(oriTrialNames);                      % number of data trials
    vT = length(valTrialNames);
    N = 100 + zeros(1, T);                  % number of data node at each data trial
    t_em = 0.1 + zeros(1, T);   % electromechanical delay
    t_rf = 0.03 + zeros(1, T);   % reflex control delay
    J = 1;                      % number of joints
    M = 4;                      % number of muscles
    S = 5;                      % number of muscle states
    C = 5;                      % number of constraints of the model dynamics
    % tPhase = 5;                 % total number of phases
    phase = {};                 % phase of the reflex controller
    
    % weights of the terms of the objective function
    W1 = 25;
    W3_list = 0;
    w11 = 10;
    w12 = 5;
    w13 = 10;
    w14 = 10;
    
    obj_ide = zeros(length(Phase_list), length(W2_list));
    obj_val = zeros(length(Phase_list), length(W2_list));

    %% load the experimental data
    for iphase = 1:length(Phase_list)
        
        tPhase = Phase_list(iphase);

        % optimization saving path
        save_path_ori = sprintf('%s/Subj%02d/%dPhases', ...
                            homeOrgPath, subj,  tPhase);
        save_path_val = sprintf('%s/Subj%02d/%dPhases', ...
                            homeValPath, subj, tPhase);

        % initialize the matrix of optimization inputs
        mus_act = zeros(sum(N), M);
        torque = zeros(sum(N), J);
        lmt = zeros(sum(N), M);
        ma = zeros(sum(N), M);
        hs = zeros(1, T);

        for t = 1:T  % load trial data and get muscle parameters from the averaged gait data

            trial = oriTrialNames(t);

            musPar = load(sprintf('%s/Subj%02d/Subj%02d_%s.mat', ...
                homeDataPath, subj, subj, trial));

            % get muscle activation
            mus_act(sum(N(1:t-1))+1:sum(N(1:t)), :) = musPar.idParaData.mus_act(1:N(t), 1:M);

            % get joint torques
            torque(sum(N(1:t-1))+1:sum(N(1:t)), :) = musPar.idParaData.torque(1:N(t), 1:J);

            % get muscle length and moment arms
            lmt(sum(N(1:t-1))+1:sum(N(1:t)), :) = musPar.idParaData.lmt(1:N(t), 1:M);
            ma(sum(N(1:t-1))+1:sum(N(1:t)), :) = musPar.idParaData.ma(1:N(t), 1:M);

            % get the time interval
            hs(t) = musPar.idParaData.hs;

            phase{t} = musPar.idParaData.(sprintf('phase%d', tPhase));

        end

        % phase_names = ["Early Stance", "Mid Stance", "Late Stance", "Early Swing", "Late Swing"];
          musOptPar = load(sprintf('%s/Subj%02d/mus_par.mat', musOptPath, subj));
         par_mus = musOptPar.mus_par; 

        % number of totoal phase of the reflex controller. Right now, only stance
        % and swing are separated.

        nPar_rf = 2*tPhase*M*M + 3*M;  % number of reflex control parameters
        
        for iW2 = 1:length(W2_list)
            for iW3 = 1:length(W3_list)
                
                W2 = W2_list(iW2);
                W3 = W3_list(iW3);
                
                folder_ide = sprintf('%s/W2_%d', save_path_ori, W2);
                folder_val = sprintf('%s/W2_%d', save_path_val, W2);

                ide_res = load(sprintf('%s/best_res.mat', folder_ide));
                x = reshape(ide_res.res.states', sum(N)*M*S, 1)';

                obj_ide(iphase, iW2) ...
                    = objective_RPO_Fit(x, M, S, N, lmt, par_mus, W1, ...
                                  w11, w12, w13, w14, torque, mus_act, ma);
                
                val_res_file = sprintf('%s/best_res.mat', folder_val);
                if exist(val_res_file, 'file')
                    val_res = load(sprintf('%s/best_res.mat', folder_val));
                    obj_val(iphase, iW2) = val_res.res.obj;
                else
                    obj_val(iphase, iW2) = NaN;
                end
                
            end
        end
    end
    
    fig1 = objPlotHeatMap(obj_ide, obj_val, Phase_list, W2_list);
    save_fig_file = sprintf('%s/resSummary/Subj%02d/WalkComb%d/Obj_ide_val_HeatMap.fig', ...
                    homeSavePath, subj, T);
    mkdir(sprintf('%s/resSummary/Subj%02d/WalkComb%d/', ...
                    homeSavePath, subj, T));
    savefig(fig1, save_fig_file);
    
end

function fig = objPlotHeatMap(obj_ide, obj_val, phase_list, W2_list)

    fig = figure();
    
    subplot(1, 2, 1)
    heatmap(phase_list, W2_list, obj_ide')
    xlabel('Gait phase division')
    ylabel('Sparsity weight')
    title('Identification')
    
    subplot(1, 2, 2)
    heatmap(phase_list, W2_list, obj_val')
    xlabel('Gait phase division')
    ylabel('Sparsity weight')
    title('Validation')
    
end

% function fig = objPlotSummarize(obj_ide, obj_val, phase_list, W2_list)
% 
%     % plot the objective function when W3 = 0
%     fig = figure();
% %     color_list = ["r", "b", "m", "k", "g", "c"];
% 
%     [W2_list_mesh, phase_list_mesh] = meshgrid(W2_list, phase_list);
%     
%     subplot(1, 2, 1)
%     s1 = surface(phase_list_mesh, W2_list_mesh, obj_ide, 'FaceAlpha', 0.5);
%     legend('Ide')
%     xlabel('Phase Number')
%     ylabel('W2 Weight')
%     zlabel('Fittness')
%     subplot(1, 2, 2)
%     s2 = surface(phase_list_mesh, W2_list_mesh, obj_val, 'FaceAlpha', 0.5);
%     xlabel('Phase Number')
%     ylabel('W2 Weight')
%     zlabel('Fittness')
%     legend('Val')
%     sgtitle('Identification & Validation Objectives')
% 
% end


