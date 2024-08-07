%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivative test of the objective function
%
% By: Huawei Wang
% Date: 12/06/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

% initlize parameters
% initlize parameters
T = 2;                      % number of data trials
N = [50, 60];                  % number of data node at each data trial
t_em = 0.1 + zeros(1, T);   % electromechanical delay
t_rf = 0.03 + zeros(1, T);   % reflex control delay
J = 1;                      % number of joints
M = 4;                      % number of muscles
S = 5;                      % number of muscle states
C = 5;                      % number of constraints of the model dynamics
tPhase = 5;                 % total number of phases
phase = {};                 % phase of the reflex controller

% specify the data trials
trialNames = ["walk_36", "run_81"];

%% load the experimental data

% initialize the matrix of optimization inputs
mus_act = zeros(sum(N), M);
torque = zeros(sum(N), J);
lmt = zeros(sum(N), M);
ma = zeros(sum(N), M);
hs = zeros(1, T);

for t = 1:T  % load trial data and get muscle parameters from the averaged gait data

    trial = trialNames(t);

    load(sprintf('idParaData/%s.mat', trial));

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
    
    par_mus = idParaData.mus_par; 

end

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

W1 = 1;
W2 = 1;
W3 = 1;

grad_equ = gradient_RPO_L0(x, M, S, N, nPar_rf, lmt, par_mus,...
                            tPhase, W1, W2, W3, torque, mus_act, ma);


% finite differentiation
delta = 1e-4;
% differentiation of muscle fiber lengths lce
for ia = 1:length(x)
    
   % get values with upper change
   delta_x = x(ia)*delta;
   
   if abs(delta_x) < 1e-6
       delta_x = 1e-6;
   end
   
   x(ia) = x(ia) + delta_x;
   obj_up = objective_RPO_L0(x, M, S, N, nPar_rf, lmt, par_mus,...
     tPhase, W1, W2, W3, torque, mus_act, ma);
   
   % get values with upper change
   x(ia) = x(ia) - 2*delta_x;
   obj_dw = objective_RPO_L0(x, M, S, N, nPar_rf, lmt, par_mus,...
     tPhase, W1, W2, W3, torque, mus_act, ma);

   % change back to the original value
   x(ia) = x(ia) + delta_x;

   % calculate the finite differentiation
   grad_fd(ia) = (obj_up - obj_dw)/(2*delta_x);

end

diff_grad = [grad_equ', grad_fd', grad_equ' - grad_fd', (grad_equ' - grad_fd')./grad_equ'];

% check the drivative differences
tolerance = 1e-4;

% difference check of df_da
errorid_grad = diffEvaluate(grad_equ, grad_fd, tolerance);

if ~isempty(errorid_grad)
   fprintf('Differentiations in grad beyond thresholds\n')
end