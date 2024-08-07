%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Validation of the identified reflex controllers and muscle parameters
% 
% By: Huawei Wang
% Date: April 13th 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function obj_fit = fitnessValueCal(homeDataPath, oriTrialNames, musOptPath, ...
                                   res, subj, T)

    % T: number of data trials

    % initlize parameters
    
    N = 100 + zeros(1, T);      % number of data node at each data trial
    J = 1;                      % number of joints
    M = 4;                      % number of muscles
    S = 5;                      % number of muscle states
    
    % weights of the terms of the objective function
    W1 = 25;
    w11 = 10;
    w12 = 5;
    w13 = 10;
    w14 = 10;

    %% load the experimental data

    % initialize the matrix of optimization inputs
    mus_act = zeros(sum(N), M);
    torque = zeros(sum(N), J);
    lmt = zeros(sum(N), M);
    ma = zeros(sum(N), M);
    
    % load trial data and get muscle parameters from the averaged gait data
    for t = 1:T  

        trial = oriTrialNames(t);  % 

        musPar = load(sprintf('%s/Subj%02d/Subj%02d_%s.mat', ...
            homeDataPath, subj, subj, trial));

        % get muscle activation
        mus_act(sum(N(1:t-1))+1:sum(N(1:t)), :) = musPar.idParaData.mus_act(1:N(t), 1:M);

        % get joint torques
        torque(sum(N(1:t-1))+1:sum(N(1:t)), :) = musPar.idParaData.torque(1:N(t), 1:J);

        % get muscle length and moment arms
        lmt(sum(N(1:t-1))+1:sum(N(1:t)), :) = musPar.idParaData.lmt(1:N(t), 1:M);
        ma(sum(N(1:t-1))+1:sum(N(1:t)), :) = musPar.idParaData.ma(1:N(t), 1:M);
        
    end

        % phase_names = ["Early Stance", "Mid Stance", "Late Stance", "Early Swing", "Late Swing"];
         musOptPar = load(sprintf('%s/Subj%02d/mus_par.mat', musOptPath, subj));
         par_mus = musOptPar.mus_par;
        
         x = reshape(res.states', sum(N)*M*S, 1)';

         obj_fit = objective_RPO_Fit(x, M, S, N, lmt, par_mus, W1, ...
                                     w11, w12, w13, w14, torque, mus_act, ma);

end
