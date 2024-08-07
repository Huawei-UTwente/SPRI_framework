%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Validation of the identified reflex controllers and muscle parameters
% 
% By: Huawei Wang
% Date: April 13th 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;

% define global parameters
global lmt torque mus_act ma par_mus N M S C T W1 W2 W3 tPhase nPar_rf ...
       t_em t_rf hs phase;
   
subj = 6;

% specify the data trials
trialNames = ["walk_09", "walk_18", "walk_27", "walk_36", "walk_45", "walk_54", ...
              "run_63", "run_81", "run_99"];
% trialNames = ["run_63", "run_81", "run_99"]; 

for iT = 1:8

% initlize parameters
T = 1;                      % number of data trials
N = 100 + zeros(1, T);                  % number of data node at each data trial
t_em = 0.0 + zeros(1, T);   % electromechanical delay
t_rf = 0.0 + zeros(1, T);   % reflex control delay
J = 1;                      % number of joints
M = 4;                      % number of muscles
S = 5;                      % number of muscle states
C = 5;                      % number of constraints of the model dynamics
% tPhase = 5;                 % total number of phases
phase = {};                 % phase of the reflex controller

    phase_list = [2, 5, 10, 25, 50, 100];

    %% load the experimental data
    for tPhase = phase_list

        % optimization saving path
        save_path = sprintf('D:/HuaweiWang/ReflexOptResults2/LocomComb%d/%dPhases/Subj%02d/trail%d', ...
                            T, tPhase, subj, iT);

        % initialize the matrix of optimization inputs
        mus_act = zeros(sum(N), M);
        torque = zeros(sum(N), J);
        lmt = zeros(sum(N), M);
        ma = zeros(sum(N), M);
        hs = zeros(1, T);

        for t = 1:T  % load trial data and get muscle parameters from the averaged gait data

            trial = trialNames(iT);

            load(sprintf('D:/HuaweiWang/ReflexIdParaData/Subj%02d/Subj%02d_%s.mat', ...
                subj, subj, trial));

            % get muscle activation
            mus_act(sum(N(1:t-1))+1:sum(N(1:t)), :) = idParaData.mus_act(1:N(t), 1:M);

            % get joint torques
            torque(sum(N(1:t-1))+1:sum(N(1:t)), :) = idParaData.torque(1:N(t), 1:J);

            % get muscle length and moment arms
            lmt(sum(N(1:t-1))+1:sum(N(1:t)), :) = idParaData.lmt(1:N(t), 1:M);
            ma(sum(N(1:t-1))+1:sum(N(1:t)), :) = idParaData.ma(1:N(t), 1:M);

            % get the time interval
            hs(t) = idParaData.hs;

            phase{t} = idParaData.(sprintf('phase%d', tPhase));

        end

        phase_names = ["Early Stance", "Mid Stance", "Late Stance", "Early Swing", "Late Swing"];

         par_mus = idParaData.mus_par;
         mus_names = idParaData.mus_names;

        % number of totoal phase of the reflex controller. Right now, only stance
        % and swing are separated.

        nPar_rf = 2*tPhase*M*M + 3*M;  % number of reflex control parameters

        % weights of the terms of the objective function
        W1 = 25;
        W2_list = [0.0, 0.5, 2, 8, 25];
        W3_list = 0;
        
        obj.locom.(sprintf('trial%d', iT)).(sprintf('tPhase%d', tPhase)) = zeros(length(W2_list), length(W3_list));

        for iW2 = 1:length(W2_list)
            for iW3 = 1:length(W3_list)
                
                W2 = W2_list(iW2);
                W3 = W3_list(iW3);
                
                w11 = 1;
                w12 = 5;
                w13 = 10;
                w14 = 10;

                folder = sprintf('%s/W2_%d_W3_%d', save_path, W2, W3);

                load(sprintf('%s/best_res.mat', folder));
                
                x = reshape(res.states', sum(N)*M*S, 1)';

                obj.locom.(sprintf('trial%d', iT)).(sprintf('tPhase%d', tPhase))(iW2, iW3) ...
                    = objective_RPO_Fit(x, M, S, N, lmt, par_mus, W1, ...
                                  w11, w12, w13, w14, torque, mus_act, ma);
            end
        end
    end
    
    fig1 = objPlotSummarize(obj.locom.(sprintf('trial%d', iT)), phase_list, W2_list, W3_list);
    sgtitle(trialNames(iT))
    save_fig_file = sprintf('D:/HuaweiWang/ReflexOptResults2/LocomComb%d/Subj%02d/Obj_trail%d.fig', ...
                    T, subj, iT);
    mkdir(sprintf('D:/HuaweiWang/ReflexOptResults2/LocomComb%d/Subj%02d', ...
                    T, subj));
    savefig(fig1, save_fig_file);
    
end



function fig1 = objPlotSummarize(obj_tr, phase_list, W2_list, W3_list)

    % plot the objective function when W3 = 0
    fig1 = figure();
%     color_list = ["r", "b", "m", "k", "g", "c"];
    for iw3 = 1:length(W3_list)
        if length(W3_list) == 1
            subplot(1, 1, 1)
        else
            subplot(ceil(length(W3_list)/2), 2, iw3)
        end
        obj_mat = zeros(length(phase_list), length(W2_list));
        for ip = 1:length(phase_list)
            tPhase = phase_list(ip);
            for iw2 = 1:length(W2_list)
                obj_mat(ip, iw2) = ...
                    obj_tr.(sprintf('tPhase%d', tPhase))(iw2, iw3);
            end
        end
        
        [W2_list_mesh, phase_list_mesh] = meshgrid(W2_list, phase_list);
        s = surface(phase_list_mesh, W2_list_mesh, obj_mat, 'FaceAlpha', 0.5);
        xlabel('Phase Number')
        ylabel('W2 Weight')
        zlabel('Fittness')
        title(sprintf('W3-%d', W3_list(iw3)))
    end

end


