function save_res_RPO(saving_names, x, info, auxdata)
% save results

    M = auxdata.M;
    S = auxdata.S;
    N = auxdata.N;
    J = auxdata.J;
    nPar_rf = auxdata.nPar_rf;
    lmt = auxdata.lmt;
    par_mus = auxdata.par_mus;
    tPhase = auxdata.tPhase;
    W1 = auxdata.W1;
    W2 = auxdata.W2;
    W3 = auxdata.W3;
    w11 = auxdata.w11;
    w12 = auxdata.w12;
    w13 = auxdata.w13;
    w14 = auxdata.w14;
    torque = auxdata.torque;
    mus_act = auxdata.mus_act;
    ma = auxdata.ma;

    % generate the joint moments and muscle force and activations

    [mom_res, force_res] = muscleForceMoment_RPO(x, lmt, ma, M, N, S, J, par_mus);

    obj = objective_RPO(x, M, S, N, nPar_rf, lmt, par_mus,...
                            tPhase, W1, W2, W3, w11, w12, w13, w14, ...
                            torque, mus_act, ma);
    time = info.cpu;
    status = info.status;

    states = reshape(x(1:sum(N)*M*S), M*S, sum(N))';

    parameters = x(end-nPar_rf+1:end);

    save(saving_names, 'states', 'parameters', 'mom_res',...
        'force_res', 'obj', 'time', 'status');
            
end