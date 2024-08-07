function J = jacobian_ipopt_RPO(x, auxdata)
% sparse values of the jacobian matrix for the optimization

    M = auxdata.M;
    S = auxdata.S;
    C = auxdata.C;
    N = auxdata.N;
    nPar_rf = auxdata.nPar_rf;
    lmt = auxdata.lmt;
    par_mus = auxdata.par_mus;
    t_em = auxdata.t_em;
    t_rf = auxdata.t_rf;
    hs = auxdata.hs;
    phase = auxdata.phase;
    tPhase = auxdata.tPhase;
    row = auxdata.row;
    col = auxdata.col;

    jac_nz = jacobian_RPO(x, M, S, C, N, nPar_rf, lmt, par_mus, t_em, ...
                             t_rf, hs, phase, tPhase);
                         
    if length(row) ~= length(jac_nz)
        jac_nz = [jac_nz, zeros(1, length(row) - length(jac_nz))];
    end

    J = sparse(row, col, jac_nz);

end