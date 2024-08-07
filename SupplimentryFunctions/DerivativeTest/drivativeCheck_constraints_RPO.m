%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivative test of the contaction dynamics
%
% By: Huawei Wang
% Date: 12/06/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

subj = 8;

% specify the data trials
trialNames = ["walk_36"];

% initlize parameters
T = 1;                      % number of data trials
N = 100 + zeros(1, T);       % number of data node at each data trial
t_rf = 0.03 + zeros(1, T);   % reflex control delay
t_em = 0.1 + zeros(1, T);   % electromechanical delay
J = 1;                      % number of joints
M = 4;                      % number of muscles
S = 5;                      % number of muscle states
C = 5;                      % number of constraints of the model dynamics
phase = {};                 % phase of the reflex controller
tPhase = 5;                 % total number of phases

%% load the experimental data
% initialize the matrix of optimization inputs
mus_act = zeros(sum(N), M);
torque = zeros(sum(N), J);
lmt = zeros(sum(N), M);
ma = zeros(sum(N), M);
hs = zeros(1, T);

for t = 1:T  % load trial data and get muscle parameters from the averaged gait data

    trial = trialNames(t);

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

par_mus0 = idParaData.mus_par0;
musOptPar = load(sprintf('D:/HuaweiWang/MuscleParamOptResults8/Subj%02d/mus_par.mat', subj));
par_mus = musOptPar.mus_par; 

% number of totoal phase of the reflex controller. Right now, only stance
% and swing are separated.

nPar_rf = 2*tPhase*M*M + 3*M;  % number of reflex control parameters

%% generate initial guesses
% initialize the optimizing parameters
mus_a = mus_act(1:sum(N), :);
mus_da = zeros(size(mus_a));
for t = 1:T
    if t== 1
        mus_da(2:N(t), :) = (mus_act(2:N(t), :) - mus_act(1:N(t)-1, :))./hs(t);
    else
        mus_da(sum(N(1:t-1))+2:sum(N(1:t)), :) = ...
            (mus_act(sum(N(1:t-1))+2:sum(N(1:t)), :) ...
            - mus_act(sum(N(1:t-1))+1:sum(N(1:t))-1, :))./hs(t);
    end
end

lce = 0.05 + 0.02*rand(size(lmt(1:sum(N), :)));
dlce = -0.05 + 0.1*rand(size(lmt(1:sum(N), :)));

x0_1 = [mus_a(1:sum(N), :), mus_da(1:sum(N), :),...
        lce(1:sum(N), :), dlce(1:sum(N), :), mus_act(1:sum(N), :)];

x0_2 = reshape(x0_1', [1, sum(N)*M*S]);

par_rf_fse = rand(1, tPhase*M*M) - 2;
par_rf_lce = rand(1, tPhase*M*M) - 2;
par_rf_res = rand(1, 3*M);

x = [x0_2, par_rf_fse, par_rf_lce, par_rf_res];


Jac = jacobian_RPO(x, M, S, C, N, nPar_rf, lmt, par_mus, t_em, t_rf, hs, phase, tPhase);
[row, col] = jacobianStructure_RPO(M, S, C, N, nPar_rf, lmt, par_mus, t_em, t_rf, hs, phase, tPhase);

rCons = M*C*sum(N - ceil((t_em + t_rf)./hs) - 1) ...
            + M*(C-1)*sum(ceil((t_em + t_rf)./hs) - ceil(t_em./hs)) ...
            + M*(C-2)*sum(ceil(t_em./hs));

Jac_equ = zeros(rCons, sum(N)*M*S + nPar_rf);

for i = 1:length(Jac)
    Jac_equ(row(i), col(i)) = Jac(i);
end

% finite differentiation
delta = 1e-6;
% differentiation of muscle fiber lengths lce
for ia = 1:length(x)
    
   % get values with upper change
   delta_x = x(ia)*delta;
   
   if delta_x < delta
       delta_x = delta;
   end
   
   x(ia) = x(ia) + delta_x;
   cons_up = constraints_RPO(x, M, S, C, N, nPar_rf, lmt, par_mus, t_em, t_rf, hs, phase, tPhase);
   
   % get values with upper change
   x(ia) = x(ia) - 2*delta_x;
   cons_dw = constraints_RPO(x, M, S, C, N, nPar_rf, lmt, par_mus, t_em, t_rf, hs, phase, tPhase);

   % change back to the original value
   x(ia) = x(ia) + delta_x;

   % calculate the finite differentiation
   Jac_fd(:, ia) = (cons_up - cons_dw)/(2*delta_x);

end

% check the drivative differences
tolerance = 1e-4;

% difference check of df_da
errorid_Jac = diffEvaluate(Jac_equ, Jac_fd, tolerance);

if ~isempty(errorid_Jac)
   fprintf('Differentiations in Jac beyond thresholds\n')
end