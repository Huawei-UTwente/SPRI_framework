function [row, col] = jacobianstructure_ipopt_RPO_rc(auxdata)
% row and col of sparse jacobian structure for the identification problem
    
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
    
    [row, col] = jacobianStructure_RPO(M, S, C, N, nPar_rf, lmt, par_mus, t_em, t_rf, hs, phase, tPhase);

end