%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function defined the effects of electromechancial delay on the jacobian
% structure
%
% By: Huawei Wang
% Date: April 14th, 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function df_dx_cb = jacDelay_rf(n_em, n_rf, df_dx_em1, df_dx_em2, ...
                                df_dx_em3, df_dx_rf1, df_dx_rf2, ...
                                df_dx_rf3, M, S)

    % only df_dx_em and df_dx_rf have nonzero element 
    % between the row 4*M to 5*M.

    % depends on how large is the _rf delay, x_rf and x_em
    % maybe constructed using some common data nodes.
    % They should be combined together when extracting
    % the nonzero Jac elements.

    % if n_em and (n_em + n_rf) is an integer multiple the time interval (hs)
    if ceil(n_em) == n_em && ceil(n_em + n_rf) == n_em + n_rf  

        df_dx_em = [df_dx_em2, df_dx_em3];  % df_dx_em1 does not exist
        df_dx_rf = [df_dx_rf2, df_dx_rf3];  % df_dx_rf1 does not exist
        
        % minimal size of the combine jac block is 2*M*S, then increase
        % depends on n_rf
        rbolck = 2 + ceil(n_rf);
        df_dx_cb = zeros(M, rbolck*M*S);

        % df_dx_em are always at the end
        df_dx_cb(:, (rbolck - 2)*M*S + 1:end) = df_dx_em;

        % df_dx_rf are always at the beginning
        df_dx_cb(:, 1:2*M*S) = df_dx_cb(:, 1:2*M*S) + df_dx_rf;

    % if n_em + n_rf is not integer times of data node anymore
    elseif ceil(n_em) == n_em && ceil(n_em + n_rf) ~= n_em + n_rf

        df_dx_em = [df_dx_em2, df_dx_em3];  % df_dx_em1 does not exist
        df_dx_rf = [df_dx_rf1, df_dx_rf2, df_dx_rf3];  % df_dx_rf1 does exist
        
        % minimal size of the combine jac block is 3*M*S, then increase
        % depends on n_rf
        rbolck = 3 + floor(n_rf);
        df_dx_cb = zeros(M, rbolck*M*S);

        % df_dx_em are always at the end
        df_dx_cb(:, (rbolck - 2)*M*S + 1:end) = df_dx_em;

        % df_dx_rf are always at the beginning
        df_dx_cb(:, 1:3*M*S) = df_dx_cb(:, 1:3*M*S) + df_dx_rf;
        
    elseif ceil(n_em) ~= n_em && ceil(n_em + n_rf) == n_em + n_rf
        
        df_dx_em = [df_dx_em1, df_dx_em2, df_dx_em3];  % df_dx_em1 does exist
        df_dx_rf = [df_dx_rf2, df_dx_rf3];  % df_dx_rf1 does not exist
        
        % minimal size of the combine jac block is 3*M*S, then increase
        % depends on n_rf
        rbolck = 3 + ceil(n_em + n_rf) - ceil(n_em);
        df_dx_cb = zeros(M, rbolck*M*S);

        % df_dx_em are always at the end
        df_dx_cb(:, (rbolck - 3)*M*S + 1:end) = df_dx_em;

        % df_dx_rf are always at the beginning
        df_dx_cb(:, 1:2*M*S) = df_dx_cb(:, 1:2*M*S) + df_dx_rf;
        
    else
        
        df_dx_em = [df_dx_em1, df_dx_em2, df_dx_em3];  % df_dx_em1 does exist
        df_dx_rf = [df_dx_rf1, df_dx_rf2, df_dx_rf3];  % df_dx_rf1 does exist
        
        % minimal size of the combine jac block is 3*M*S, then increase
        % depends on n_rf
        rbolck = 3 + ceil(n_em + n_rf) - ceil(n_em);
        df_dx_cb = zeros(M, rbolck*M*S);

        % df_dx_em are always at the end
        df_dx_cb(:, (rbolck - 3)*M*S + 1:end) = df_dx_em;

        % df_dx_rf are always at the beginning
        df_dx_cb(:, 1:3*M*S) = df_dx_cb(:, 1:3*M*S) + df_dx_rf;

    end
end
