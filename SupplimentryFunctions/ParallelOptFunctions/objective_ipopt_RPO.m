function f = objective_ipopt_RPO(x, auxdata)
% objective function of the id problem
    M = auxdata.M;
    S = auxdata.S;
    N = auxdata.N;
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
    
    % calculate objective function
    f =  objective_RPO(x, M, S, N, nPar_rf, lmt, par_mus,...
                            tPhase, W1, W2, W3, w11, w12, w13, w14, ...
                            torque, mus_act, ma);

end