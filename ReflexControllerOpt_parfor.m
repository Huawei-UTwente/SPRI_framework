%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimization of muscle reflex controllers, with electromechanical +
% neural transition delays.
%
% By: Huawei Wang
% Date: 12/06/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

% global lmt torque mus_act ma par_mus N M S C T W1 W2 W3 tPhase nPar_rf ...
%        t_em t_rf hs phase rCons cCons;

subj = 6;
% specify the data trials
trialNames = ["walk_09", "walk_18", "walk_27", "walk_36", "walk_45", "walk_54", ...
               "run_63", "run_81", "run_99"]; 

% trialNames = ["run_63", "run_81", "run_99"]; 

% initlize parameters
for iT = 1            % number of data trials
    T = 4;
    N = 100 + zeros(1, T);      % number of data node at each data trial
    t_rf = 0.03 + zeros(1, T);   % reflex control delay
    J = 1;                      % number of joints
    M = 4;                      % number of muscles
    S = 5;                      % number of muscle states
    C = 5;                      % number of constraints of the model dynamics
    
    for tPhase = [2, 5, 10, 25]
        % tPhase = 100;             % total number of phases
        phase = {};                 % phase of the reflex controller

        %% load the experimental data
        
        % initialize the matrix of optimization inputs
        mus_act = zeros(sum(N), M);
        torque = zeros(sum(N), J);
        lmt = zeros(sum(N), M);
        ma = zeros(sum(N), M);
        hs = zeros(1, T);
        
        % optimization saving path
        save_path = sprintf('D:/HuaweiWang/ReflexOptResults2/WalkComb%d/%dPhases/Subj%02d/trail%d', ...
                            T, tPhase, subj, iT);
        mkdir(save_path);

        for t = 1:T  % load trial data and get muscle parameters from the averaged gait data

            trial = trialNames(1+t);

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

            t_em(t) = hs(t)*6;   % electromechanical delay, 6 data points delay

            phase{t} = idParaData.(sprintf('phase%d', tPhase));

        end

         par_mus = idParaData.mus_par; 

        % number of totoal phase of the reflex controller. Right now, only stance
        % and swing are separated.

        nPar_rf = 2*tPhase*M*M + 3*M;  % number of reflex control parameters

        %% optimization boundaries and setups
        x_lb1 = [zeros(sum(N), M) + 0.001,...          % activation
                zeros(sum(N), M) - 30,...                  % d_activation
                zeros(sum(N), M) + par_mus(1:M)*0.5,...    % lce
                zeros(sum(N), M) - par_mus(1:M)*15, ....   % dlce
                zeros(sum(N), M) + 0.001];                 % nerual stimulation

        x_lb2 = reshape(x_lb1', [1, sum(N)*M*S]);

        par_rf_lb = [-10*ones(1, 2*tPhase*M*M), zeros(1, M), 0.5 + zeros(1, M), zeros(1, M)];

        x_lb3 = [x_lb2, par_rf_lb];

        % set up optimizing parameter upper bounds.  
        % lce can be maximum 2 time the lce_opt
        x_ub1 = [zeros(sum(N), M) + 1,...                  % activation
                 zeros(sum(N), M) + 30,...                 % d_activation
                 zeros(sum(N), M) + par_mus(1:M)*1.5,...   % lce
                 zeros(sum(N), M) + par_mus(1:M)*15, ...   % dlce
                 zeros(sum(N), M) + 1];                    % nerual stimulation

        x_ub2 = reshape(x_ub1', [1, sum(N)*M*S]);

        % baseline activation shouldn't larger than measured minimual
        % activations
        par_rf_ub = [10*ones(1, 2*tPhase*M*M), zeros(1, M), ...  % no force thresholds
                                        1.5 + zeros(1, M), max([min(mus_act); zeros(1, M)])];  
                                    

        x_ub3 = [x_ub2, par_rf_ub];

        %% run optimization
        W1 = 25;
        for W2 = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.5, 0.75, 1, 1.5, 2, 4, 6, 8, 15, 25]
            for W3 = 0
                w11 = 1;
                w12 = 5;
                for w13 = 10
                    for w14 = 10

                        auxdata.M = M;
                        auxdata.S = S;
                        auxdata.C = C;
                        auxdata.N = N;
                        auxdata.T = T;
                        auxdata.J = J;
                        auxdata.nPar_rf = nPar_rf;
                        auxdata.lmt = lmt;
                        auxdata.par_mus = par_mus;
                        auxdata.t_em = t_em;
                        auxdata.t_rf = t_rf;
                        auxdata.hs = hs;
                        auxdata.phase = phase;
                        auxdata.tPhase = tPhase;
                        auxdata.torque = torque;
                        auxdata.mus_act = mus_act;
                        auxdata.ma = ma;
                        auxdata.W1 = W1;
                        auxdata.W2 = W2;
                        auxdata.W3 = W3;
                        auxdata.w11 = w11;
                        auxdata.w12 = w12;
                        auxdata.w13 = w13;
                        auxdata.w14 = w14;

                        [row, col] = jacobianstructure_ipopt_RPO_rc(auxdata);
                        auxdata.row = row;
                        auxdata.col = col;

                        options.lb = x_lb3;
                        options.ub = x_ub3;

                        rCons = M*C*sum(N - ceil((t_em + t_rf)./hs) - 1) ...
                                + M*(C-1)*sum(ceil((t_em + t_rf)./hs) - ceil(t_em./hs)) ...
                                + M*(C-2)*sum(ceil(t_em./hs));

                        cCons = sum(N)*M*S + nPar_rf;

                        options.cl = zeros(1, rCons);
                        options.cu = zeros(1, rCons);

                        % Set the IPOPT options.
                        options.ipopt.hessian_approximation = 'limited-memory';
                        options.ipopt.mu_strategy           = 'adaptive';
                        options.ipopt.tol                   = 1e-4;
                        options.ipopt.max_iter              = 5000;
                        options.ipopt.linear_solver         = 'mumps';

                        % The callback functions.
                        funcs.objective         = @(x) objective_ipopt_RPO(x, auxdata);
                        funcs.constraints       = @(x) constraints_ipopt_RPO(x, auxdata);
                        funcs.gradient          = @(x) gradient_ipopt_RPO(x, auxdata);
                        funcs.jacobian          = @(x) jacobian_ipopt_RPO(x, auxdata);
                        funcs.jacobianstructure = @() jacobianstructure_ipopt_RPO(auxdata);

                        auxdata.options = options;
                        auxdata.funcs = funcs;

                        folder = sprintf('%s/W2_%d_W3_%d', save_path, W2, W3);
                        mkdir(folder);

                        auxdata.folder = folder;

                        parfor opt = 1:100
                            do_optimization_mus(opt, auxdata)
                        end
                        
                    end
                end
            end
        end
    end
end
