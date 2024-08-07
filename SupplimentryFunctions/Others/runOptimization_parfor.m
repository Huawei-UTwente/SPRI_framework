function runOptimization_parfor(i, x0, options, save_path)

    %% generate initial guesses

    % The callback functions.
    funcs.objective         = @objective_ipopt;
    funcs.constraints       = @constraints_ipopt;
    funcs.gradient          = @gradient_ipopt;
    funcs.jacobian          = @jacobian_ipopt;
    funcs.jacobianstructure = @jacobianstructure_ipopt;

    [x, info] = ipopt(x0, funcs, options);

    saveOptResults(x, info, i, save_path)  % save optimized results
        
end