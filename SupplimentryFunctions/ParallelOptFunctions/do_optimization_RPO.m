function do_optimization_RPO(opt, auxdata)

    folder = auxdata.folder;
    options = auxdata.options;
    funcs = auxdata.funcs;
    
    M = auxdata.M;
    S = auxdata.S;
    N = auxdata.N;
    T = auxdata.T;
    hs = auxdata.hs;
%     nPar_rf = auxdata.nPar_rf;
%     lmt = auxdata.lmt;
    par_mus = auxdata.par_mus;
    tPhase = auxdata.tPhase;
%     W1 = auxdata.W1;
%     W2 = auxdata.W2;
%     W3 = auxdata.W3;
%     torque = auxdata.torque;
    mus_act = auxdata.mus_act;
%     ma = auxdata.ma;
    lb = options.lb;
    ub = options.ub;
    nPar_rf = auxdata.nPar_rf;

    %% generate initial guesses
    rng(opt)
    % initialize the optimizing parameters
    mus_a = (0.95 + 0.1*rand(sum(N), M)).*mus_act(1:sum(N), :);
    mus_da = zeros(size(mus_a));
    for t = 1:length(N)
        if t== 1
            mus_da(2:N(t), :) = (mus_act(2:N(t), :) - mus_act(1:N(t)-1, :))./hs(t);
        else
            mus_da(sum(N(1:t-1))+2:sum(N(1:t)), :) = ...
                (mus_act(sum(N(1:t-1))+2:sum(N(1:t)), :) ...
                - mus_act(sum(N(1:t-1))+1:sum(N(1:t))-1, :))./hs(t);
        end
    end

    mus_da(1, :) = (0.95 + 0.1*rand(1, M)).*mus_da(2, :);
    mus_s = (0.9 + 0.2*rand(sum(N), M)).*mus_act(1:sum(N), :);

    lce = (0.95 + 0.1*rand(sum(N), M)).*par_mus(1:M);
    dlce = (0.5 - 1*rand(sum(N), M)).*par_mus(1:M);

    x0_1 = [mus_a(1:sum(N), :), mus_da(1:sum(N), :),...
        lce(1:sum(N), :), dlce(1:sum(N), :), mus_s];

    x0_2 = reshape(x0_1', [1, sum(N)*M*S]);
    
%     par_rf_fse = -20*rand(1, tPhase*M*M) + 10;
%     par_rf_lce = -20*rand(1, tPhase*M*M) + 10;
%     par_rf_res = [zeros(1, M), 0.5 + rand(1, M), 0.5*rand(1, M)];
    
    par_rf0 = lb(end - nPar_rf + 1:end) ....
        + (ub(end - nPar_rf + 1:end) - lb(end - nPar_rf + 1:end)).*rand(1, nPar_rf);

    x0 = [x0_2, par_rf0];

    %% do optimizations with ipopt

    [x, info] = ipopt(x0, funcs, options);

    saving_names = sprintf('%s/optimization_res%02d.mat', folder, opt);
    
    save_res_RPO(saving_names, x, info, auxdata)
    
end