%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is to optimize the ankle joint reflex control with the
% experimental data.
%
% The experimental data consists with joint angles, joint torques, and
% muscle activations (EMG).
%
% The ankle joint reflex will general muscle activations for the ankle
% muscles based on the reflex feedback. 
%
% The ankle joint muscles including: Soleus; Tibiais Anterior;
% M/L Gastrocneius. 
%
% Optimization is based on the averaged half gait cycle's data. gradient
% based optimizer is used here to find the best reflex control gains that
% can explain the experimental data.
%
% By: Huawei Wang
% Date: Augst 1, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all

% define global parameters
global lmt torque mus_act ma par_mus N M S C T W1 W2 W3 tPhase nPar_rf ...
       t_em t_rf hs phase;
   
% specify the data trials
trialNames = ["walk_09", "walk_18", "walk_27", "walk_36", "walk_45", "walk_54", ...
              "run_63", "run_81", "run_99"];
          
musOptPath = 'D:\HuaweiWang\MuscleParamOptResults';  % path for saving results

% trialNames = ["run_63", "run_81", "run_99"];

subj = 5;

for iT = 5                      % number of data trials

    % initlize parameters
    T = 1;                      % number of data trials
    N = 100 + zeros(1, T);                  % number of data node at each data trial
    t_em = 0.0 + zeros(1, T);   % electromechanical delay
    t_rf = 0.0 + zeros(1, T);   % reflex control delay
    J = 1;                      % number of joints
    M = 4;                      % number of muscles
    S = 5;                      % number of muscle states
    C = 5;                      % number of constraints of the model dynamics
%     tPhase = 50;                 % total number of phases
%     phase = {};                 % phase of the reflex controller

    for tPhase = [2]  % , 50, 100
        % tPhase = 100;                 % total number of phases
        phase = {};                 % phase of the reflex controller


        %% load the experimental data

        % optimization saving path
        save_path = sprintf('D:/HuaweiWang/ReflexIdentResults/IndividualTest/Subj%02d/%dPhases/trial_%s', ...
                            subj, tPhase, trialNames(iT));

        % initialize the matrix of optimization inputs
        mus_act = zeros(sum(N), M);
        torque = zeros(sum(N), J);
        lmt = zeros(sum(N), M);
        ma = zeros(sum(N), M);
        hs = zeros(1, T);

        for t = 1:T  % load trial data and get muscle parameters from the averaged gait data

            trial = trialNames(1 + t);

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

        phase_names = string(linspace(1, tPhase, tPhase));

         musOptPar = load(sprintf('%s/Subj%02d/mus_par.mat', musOptPath, subj));
         par_mus = musOptPar.mus_par; 
         mus_names = idParaData.mus_names;

        % number of totoal phase of the reflex controller. Right now, only stance
        % and swing are separated.

        nPar_rf = 2*tPhase*M*M + 3*M;  % number of reflex control parameters

        % weights of the terms of the objective function
        W1 = 25;
        for W2 = [0] % [0.0, 0.5, 2, 8, 25] % [0.25, 0.75, 1, 1.5, 4, 6, 15] %  % [0.0, 0.25, 0.5, 0.75, 1]
            for W3 = 0
                w11 = 1;
                w12 = 5;
                for w13 = 10
                    for w14 = 10

                folder = sprintf('%s/W2_%d_W3_%d', save_path, W2, W3);

                    %% load optimized results
                    
                    clear par_rf

                    for opt = 1:100

                        saving_names = sprintf('%s/optimization_res%02d.mat', folder, opt);
                        res = load(saving_names);

                        par_rf(opt, :) = res.parameters;

                        Obj_res(opt) = res.obj;
                        Status_res(opt) = res.status;
                        Time_res(opt) = res.time;
                    end

                    %% plot optimized parameters
                    plot_sign = -2;
                    succ_index = find(Status_res == plot_sign);
                    
                    if isempty(succ_index)
                        continue
                    end
                    
                    [Obj_sort, sort_index] = sort(Obj_res(succ_index));

                    par_rf_sort = par_rf(succ_index(sort_index), :);

                    Time_res_sort = Time_res(succ_index(sort_index))';

                    fig4 = figure();

                    angle = -45;

                    plot(Obj_sort, 'ro')
                    title('Objective Values')
                    xlabel('Optimization trails')

                    savefig(fig4, sprintf('%s/objectiveFunction_OptimizedParameters.fig', folder))
                    close(fig4);
                    % reshape the reflex control gains for the plots
                    nPlot = min(1, length(Obj_sort));
                    par_rf_fse_m = zeros(nPlot, tPhase, M);
                    par_rf_lce_m = zeros(nPlot, tPhase, M);
                    par_rf_res_m = zeros(nPlot, 3);


                    par_rf_mean = mean(par_rf_sort(1:nPlot, :), 1);
                    if nPlot == 1
                        par_rf_std = zeros(size(par_rf_mean));
                    else
                        par_rf_std = std(par_rf_sort(1:nPlot, :), 1);
                    end
                    
                    phase_occupied = zeros(tPhase, 1);
                    for p = 1:tPhase
                        phase_occupied(p) = length(find(phase{1} == p));
                    end

                    hf_fig = reflex_heat_map_std_dots(par_rf_mean, par_rf_std, M, mus_names, ...
                        tPhase, phase_occupied, phase_names);
                    savefig(hf_fig, sprintf('%s/reflexGains.fig', folder))
                    close(hf_fig);
                    % plot joint torque comparison
                    plot_ind = succ_index(sort_index);

                    color_vec = ['r', 'b', 'g', 'm', 'c', 'y'];

                    color_ind = 1;

                    h1 = figure();
                    joint_names = ["ANKLE MOMENT Nm"];

                    st = 1;
                    ed = 1;

                    for opt = plot_ind(st:ed)

                        saving_names = sprintf('%s/optimization_res%02d.mat', folder, opt);
                        res = load(saving_names);

                        save(sprintf('%s/best_res.mat', folder), 'res');

                        subj_h1 = ceil(sqrt(J));

                        for j = 1:J
                            subplot(subj_h1, subj_h1, j)
                            for t = 1:T
                                plot(1:N(t), torque(sum(N(1:t-1))+1:sum(N(1:t)), :),...
                                    'k-', 'linewidth', 2.5)
                                hold on
                                plot(1:N(t), res.mom_res(sum(N(1:t-1))+1:sum(N(1:t)), 1),...
                                    '--', 'linewidth', 2.5, 'color', color_vec(color_ind))
                                hold on

                                if t == 1
                                    title(joint_names(j))
                                end

                                if j == 1
                                    trail_name = trialNames(t);
                                    ylabel(trail_name)
                                end            
                                if t == T
                                    legend("experimental Joint Torques", ...
                                           "optimized Joint Torques")
                                end
                            end
                        end
                        color_ind = color_ind + 1;
                    end

                    savefig(h1, sprintf('%s/JointTorqueFits.fig', folder))
                    close(h1);

                    h2 = figure();
                    color_ind = 1;

                    for opt = plot_ind(st:ed)

                        saving_names = sprintf('%s/optimization_res%02d.mat', folder, opt);
                        res = load(saving_names);

                        for t = 1:T
                            for m = 1:M
                                subplot(T, M, (t-1)*M + m)
                                plot(1:N(t), res.force_res(sum(N(1:t-1)) + 1:sum(N(1:t)), m),...
                                    '--', 'linewidth', 2.5, 'color', color_vec(color_ind))
                                hold on

                                if t == 1
                                    title(mus_names(m))
                                end

                                if m == 1
                                    trail_name = trialNames(t);
                                    ylabel(trail_name)
                                end

                            end
                            sgtitle('Muscle Forces')
                        end
                        color_ind = color_ind + 1;
                    end
                    savefig(h2, sprintf('%s/MuscleForces.fig', folder))
                    close(h2);
                    %%    
                    h3 = figure();
                    color_ind = 1;

                    for opt = plot_ind(st:ed)

                        saving_names = sprintf('%s/optimization_res%02d.mat', folder, opt);
                        res = load(saving_names);

                        for t = 1:T
                            activations = zeros(N(t), M);
                            for n = 1:N(t)
                                activations(n, :) = res.states(sum(N(1:t - 1)) + n, (S - 1)*M + 1:S*M);
                            end

                            for m = 1:M

                                subplot(T, M, (t-1)*M + m)
                                plot(1:N(t), mus_act(sum(N(1:t - 1)) + 1:sum(N(1:t)), m),...
                                    'k-', 'linewidth', 2.5)
                                hold on
                                plot(1:N(t), res.states(sum(N(1:t - 1)) + 1:sum(N(1:t)), (S - 1)*M + m),...
                                    '--', 'linewidth', 2.5, 'color', color_vec(color_ind))

                                ylim([0, 1])

                                if t == 1
                                    title(mus_names(m))
                                end

                                if m == 1
                                    trail_name = trialNames(t);
                                    ylabel(trail_name)
                                end

                            end
                            sgtitle('Muscle Activation')
                        end

                        color_ind = color_ind + 1;
                    end
                        savefig(h3, sprintf('%s/MuscleActivation.fig', folder))
                        close(h3);
                    end
                end
            end
        end
    end
end
