function saveOptResults(x, info, i, save_path)

global lmt ma moment mus_act muscle_par0 L J N M P S T W1 W2 W3 W4 W5;

%% save results

% generate the joint moments and muscle force and activations

[mom_res, force_res] = muscleForceMoment(x, lmt, ma, L, P, T, J, M, N, S);

obj = objective_MPO(x, lmt, moment, mus_act, ma, T, L, J, N, M, S, P,...
                muscle_par0, W1, W2, W3, W4, W5);
time = info.cpu;
status = info.status;

states = reshape(x(1:S*M*N*L*T), S*M*L*T, N)';

activations = reshape(x(S*M*N*L*T + 1:S*M*N*L*T + M*N*L*T), M*L*T, N)';

parameters = x(end-P+1:end);

saving_names = sprintf('%s/optimization_res%02d.mat', save_path, i);

save(saving_names, 'states', 'activations', 'parameters', 'mom_res',...
'force_res', 'obj', 'time', 'status');

end