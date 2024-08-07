function [row, col] = jacobianStructure_RPO(M, S, C, N, nPar_rf, lmt, ...
                                    par_mus, t_em, t_rf, hs, phase, tPhase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structure of the Jacobian matrix. Specifically, row and col indexes of
% nonzero elements. Jacobian should be stored in sparse matrix structure,
% to avoid store large matrixes. These row and col indexes will be used to
% formulate the sparse Jacobian matrix for optimizers, like Ipopt.
%
% By: Huawei Wang
% Date: March 13, 2022
%
% INPUTS:
%   M, integer: number of muscles that are included in the musculoskeletal
%               mdoel
%   S, integer: number of states that included in each muscle model
%   C, integer: number of constraints that included in each muscle model
%   N, vector: number of data nodes in each data trial
%   nPar_rf, integer: number of control parameters in the reflex control model
%   lmt, matrix: the muscle-tendon unit lengths of optimizing muscles in
%                different data trials, predetermined by the joint angles
%   par_mus, vector: personalized muscle parameters
%   t_em, vector: electricalmechanical delay of each data trial
%   t_rf, vector: reflex control delay of each data trial
%   hs, vector: time intervals of each data trials
%   phase, cell: the phase of the gait of each data trial {1:N1; 1:N2, ...}
%   tPhase, integer: total number of phases of the identified reflex
%                    controller

% OUTPUT:
%   row, vector: the row indexes of nonzero vector of the jacobian matrix of the
%   optimization problem.
%   col, vector: the col indexes of nonzero vector of the jacobian matrix of the
%   optimization problem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % number of repetitions in generate the jac matrix. Normally larger
    % than 1, in case some location get false zero due to the random
    % selected states x.
    num_test = 2;  

    % Initilize the empty row and col vector for nonzero Jacobian vector
    row = [];
    col = [];

    % Get the controller parameters, reshape the fse, lce reflex control 
    % gains from vector to matrix
    
    for t = 1:length(N)       % Run through all data trials
        
        hs_t = hs(t);         % Time interval between the data nodes
        n_em = t_em(t)/hs_t;  % Number of data nodes that caused by electromechanical delay (_em)
        n_rf = t_rf(t)/hs_t;  % Number of data nodes that caused by reflex delay (_rf)
        
        x_st = sum(N(1:t-1))*M*S; % The index of the first state parameter of current data trial 
        lmt_st = sum(N(1:t-1));   % The index of the first lmt data of current data trial
        
        % Looking back to the data nodes that caused by _em and _rf delays 
        stn_em1 = ceil(n_em)*M*S;
        % stn_em2 = floor(n_em)*M*S;
        
        stn_rf1 = ceil(n_em + n_rf)*M*S;
        % stn_rf2 = floor(n_em + n_rf)*M*S;
        
        % Index of the first constraint vector of the current data trial
        cons_st = M*C*sum(N(1:t-1) - 1 - ceil((t_em(1:t-1) + t_rf(1:t-1))./hs(1:t-1))) ...
            + M*(C-1)*sum(ceil((t_em(1:t-1) + t_rf(1:t-1))./hs(1:t-1)) - ceil(t_em(1:t-1)./hs(1:t-1))) ...
            + M*(C-2)*sum(ceil(t_em(1:t-1)./hs(1:t-1)));
        
        % Number of constraint vectors that do not have either _em and _rf
        % be aware that they have different number of constraints in each 
        % secinarios. With neither _em and _rf, number of constraints is
        % C - 2; with only _em, number of constraints is C - 1.
        cons_st_em = ceil(n_em)*M*(C-2);
        cons_st_rf = (ceil(n_em + n_rf) - ceil(n_em))*M*(C-1);
                
        for n = 1 : N(t)-1  % Run through all data nodes (using middle point method)
            
            x_stn = x_st + (n-1)*M*S; % The index right before the first state of current data node
            % x_stn1 = x_st + n*M*S;      % The index of the last state of current data node
            % x_stn2 = x_st + (n+1)*M*S;  % The index of last state of next data node
                     
            % If the number of current data node is smaller than the data node
            % required by the _em delay. Then the activation
            % dynamic and reflex control constraints are not included,
            % since there are no state can be used to formulate these
            % feedbacks.
            if n <= ceil(n_em)
                
                df_dx1 = zeros(M*(C-2), M*S);
                df_dx2 = zeros(M*(C-2), M*S);
                
                for test = 1:num_test   % run multiple repetitions
                    
                    % generate randomized x1 and x2, instead of extracting
                    % from the optimizing variables x. Randomization can
                    % largly avoid false zero elements.
                    x1 = rand(1, M*S);
                    x2 = rand(1, M*S);

                    x1(2*M+1:3*M) = x1(2*M+1:3*M)*0.01+0.05;
                    x1(3*M+1:4*M) = -x1(3*M+1:4*M)*0.2+ 0.1;

                    x2(2*M+1:3*M) = x2(2*M+1:3*M)*0.01+0.06;
                    x2(3*M+1:4*M) = -x2(3*M+1:4*M)*0.3 + 0.2;

                    % Calculate the x, dx, and lmt based on the middle point method
                    x_tn = (x1 + x2)/2;
                    dx_tn = (x1 - x2)./hs_t;

                    lmt_tn = (lmt(lmt_st + n, :) + lmt(lmt_st + n + 1, :))/2; 
                    
                    % directCollocationDyn_diff_RPO1 is the derivative function
                    % for directCollocationDyn_RPO1. Only df_dx, and fx_ddx are
                    % calculate, since x_em and x_rf are not available
                    [df_dx_tn, df_ddx_tn] = ...
                    directCollocationDyn_diff_RPO1(x_tn, dx_tn, ...
                        lmt_tn, par_mus, M, S);

                    % (middle point method) assign back derivatives to data nodes
                    df_dx1 = df_dx1 + df_dx_tn./2 - df_ddx_tn./hs_t;
                    df_dx2 = df_dx2 + df_dx_tn./2 + df_ddx_tn./hs_t;
                    
                end
                    
                % extract the indexes of nonzero elements in the jacobian 
                % matrix, row by row
                row_st = cons_st + (n - 1)*M*(C-2);

                for r = 1:M*(C-2)
                    [row_i, col_i] = find([df_dx1(r, :), df_dx2(r, :)]);

                    row_i = row_i + row_st + r - 1;
                    col_i = col_i + x_stn;

                    row = [row, row_i];
                    col = [col, col_i];
                end        

            % Else if the number of current data node is larger than the
            % _em, but smaller than _em + _rf. Then the activation dynamics
            % is included, the reflex control constraint will not be 
            % included.
            elseif n <= ceil(n_em + n_rf)
                
                df_dx1 = zeros(M*(C-1), M*S);
                df_dx2 = zeros(M*(C-1), M*S);

                df_dx_em1 = zeros(M*(C-1), M*S);
                df_dx_em2 = zeros(M*(C-1), M*S);
                df_dx_em3 = zeros(M*(C-1), M*S);
                
                % If n_em is not integer, the _em delay data nodes are
                % calculated using linear interpolation between the near by
                % data nodes. If n_em is integer, then directly find the
                % delayed data nodes.
                if n_em == floor(n_em)
                    w1_em = 0;
                    w2_em = 1;
                else
                    w1_em = n_em - floor(n_em);
                    w2_em = ceil(n_em) - n_em;
                end
                   
                for test = 1:num_test
                    
                    % generate randomized x1, x2, and x_em, instead of extracting
                    % from the optimizing variables x. Randomization can
                    % largly avoid false zero elements.
                    x1 = rand(1, M*S);
                    x2 = rand(1, M*S);

                    x1(2*M+1:3*M) = x1(2*M+1:3*M)*0.01+0.05;
                    x1(3*M+1:4*M) = -x1(3*M+1:4*M)*0.2+ 0.1;

                    x2(2*M+1:3*M) = x2(2*M+1:3*M)*0.01+0.06;
                    x2(3*M+1:4*M) = -x2(3*M+1:4*M)*0.3 + 0.2;

                    % calculate the x and dx based on the middle point method
                    x_tn = (x1 + x2)/2;
                    dx_tn = (x1 - x2)./hs_t;

                    lmt_tn = (lmt(lmt_st + n, :) + lmt(lmt_st + n + 1, :))/2; 
                    
                    x_em_tn = rand(1, M*S);

                    % directCollocationDyn_diff_RPO2 is the derivative function
                    % for directCollocationDyn_RPO2. df_dx_em is included here,
                    % however, df_dx_rf is still not included, since x_rf is
                    % not available yet.
                    [df_dx_tn, df_ddx_tn, df_dx_em_tn] = ...
                    directCollocationDyn_diff_RPO2(x_tn, dx_tn, x_em_tn, ...
                        lmt_tn, par_mus, M, S);
                    
                    % (middle point method) assign back derivatives to data nodes
                    df_dx1 = df_dx1 + df_dx_tn./2 - df_ddx_tn./hs_t;
                    df_dx2 = df_dx2 + df_dx_tn./2 + df_ddx_tn./hs_t;

                    df_dx_em1 = df_dx_em1 + df_dx_em_tn.*w1_em/2;
                    df_dx_em2 = df_dx_em2 + df_dx_em_tn/2;
                    df_dx_em3 = df_dx_em3 + df_dx_em_tn*w2_em/2;
                    
                end
                    
                % extract nonzero elements in the jacobian matrix, row by
                % row
                row_st = cons_st + cons_st_em + (n - ceil(n_em) - 1)*M*(C-1);

                for r = 1:M*(C-1)
                    if r <= 2*M  % derivative constraints
                        
                        [row_i, col_i] = find([df_dx1(r, :), df_dx2(r, :)]);

                        row_i = row_i + row_st + r - 1;
                        col_i = col_i + x_stn;

                    elseif r <= 3*M  % activation constraints
                        
                        df_dx_x_em = jacDelay_em(n_em, df_dx_em1(2*M+1:3*M, :),...
                            df_dx_em2(2*M+1:3*M, :), df_dx_em3(2*M+1:3*M, :),...
                            df_dx1(2*M+1:3*M, :), df_dx2(2*M+1:3*M, :), M, S);
                         
                        [row_i, col_i] = find(df_dx_x_em(r - 2*M, :));
                                
                        row_i = row_i + row_st + r - 1;
                        col_i = col_i + x_stn - stn_em1;
                        
                    elseif r <= 4*M  % muscle dynamic constraints
                        
                        [row_i, col_i] = find([df_dx1(r, :), df_dx2(r, :)]);

                        row_i = row_i + row_st + r - 1;
                        col_i = col_i + x_stn;

                    end
                    row = [row, row_i];
                    col = [col, col_i];
                end  
                
            % Else when there are states can be used to formulate the
            % reflex control loops, include the reflex control constraints.
            else
                
                df_dx1 = zeros(M*C, M*S);
                df_dx2 = zeros(M*C, M*S);

                df_dx_em1 = zeros(M*C, M*S);
                df_dx_em2 = zeros(M*C, M*S);
                df_dx_em3 = zeros(M*C, M*S);

                df_dx_rf1 = zeros(M, M*S);
                df_dx_rf2 = zeros(M, M*S);
                df_dx_rf3 = zeros(M, M*S);

                dfs_rf_dpar_rf_fse = zeros(M, M*M);
                dfs_rf_dpar_rf_lce = zeros(M, M*M);
                dfs_rf_dpar_rf_res0 = zeros(M, 3*M);

                par_rf_tn = (phase{t}(n)-1)*M*M;
                % par_rf_tn1 = (phase{t}(n))*M*M;
                % par_rf_tn2 = (phase{t}(n+1))*M*M;

                % Use the same method used above to reconstruct x_em
                if n_em == floor(n_em)
                    w1_em = 0;
                    w2_em = 1;
                else
                    w1_em = n_em - floor(n_em);
                    w2_em = ceil(n_em) - n_em;
                end

                % If the _em + _rf delay is not an integer mutiple of hs(t),
                % then use interpolation between the nearest two data nodes
                % to reconstruct the x_rf.
                if n_rf + n_em == floor(n_em + n_rf)
                    w1_rf = 0;
                    w2_rf = 1;
                else
                    w1_rf = n_rf + n_em - floor(n_em + n_rf);
                    w2_rf = ceil(n_em + n_rf) - n_rf - n_em;
                end

                for test = 1:num_test
                    
                    % generate randomized x1, x2, x_em, x_rf, and par_rf, 
                    % instead of extracting
                    % from the optimizing variables x. Randomization can
                    % largly avoid false zero elements.
                    x1 = rand(1, M*S);
                    x2 = rand(1, M*S);

                    x1(2*M+1:3*M) = x1(2*M+1:3*M)*0.01+0.05;
                    x1(3*M+1:4*M) = -x1(3*M+1:4*M)*0.2+ 0.1;

                    x2(2*M+1:3*M) = x2(2*M+1:3*M)*0.01+0.06;
                    x2(3*M+1:4*M) = -x2(3*M+1:4*M)*0.3 + 0.2;
                    
                     x_em_tn = rand(1, M*S);

                    x_rf_tn = rand(1, M*S);
                    x_rf_tn(2*M+1:3*M) = x_rf_tn(2*M+1:3*M)*0.01+0.05;
                    x_rf_tn(3*M+1:4*M) = -x_rf_tn(3*M+1:4*M)*0.2+ 0.1;

                    % calculate the x and dx based on the middle point method
                    x_tn = (x1 + x2)/2;
                    dx_tn = (x1 - x2)./hs_t;

                    lmt_tn = (lmt(lmt_st + n, :) + lmt(lmt_st + n+1, :))/2; 

                    lmt_rf1 = lmt(lmt_st + n - ceil(n_em+n_rf), :)*w1_rf ...
                       + lmt(lmt_st + n - floor(n_em+n_rf), :)*w2_rf;
                   
                    lmt_rf2 = lmt(lmt_st + n + 1 - ceil(n_em+n_rf), :)*w1_rf ...
                           + lmt(lmt_st + n + 1 - floor(n_em+n_rf), :)*w2_rf;

                    lmt_rf_tn = (lmt_rf1 + lmt_rf2)/2;

                    par_rf_fse_tn = rand(1, M*M) - 2;
                    par_rf_lce_tn = rand(1, M*M) - 2;
                    par_rf_res = rand(1, 3*M);

                    % directCollocationDyn_diff_RPO3 is the derivative function
                    % for directCollocationDyn_RPO3. All derivatives are
                    % included.
                    [df_dx_tn, df_ddx_tn, df_dx_em_tn, dfs_dx_rf_tn, dfs_rf_dpar_rf_fse_tn, ...
                     dfs_rf_dpar_rf_lce_tn, dfs_rf_dpar_rf_res] = ...
                     directCollocationDyn_diff_RPO3(x_tn, dx_tn, ...
                            x_em_tn, x_rf_tn, lmt_tn, lmt_rf_tn, par_rf_fse_tn, ...
                            par_rf_lce_tn, par_rf_res, par_mus, M, S);

                    % (middle point method) assign back derivatives to data nodes
                    df_dx1 =  df_dx1 + df_dx_tn./2 - df_ddx_tn./hs_t;
                    df_dx2 = df_dx2 + df_dx_tn./2 + df_ddx_tn./hs_t;

                    df_dx_em1 = df_dx_em1 + df_dx_em_tn.*w1_em/2;
                    df_dx_em2 = df_dx_em2 + df_dx_em_tn/2;
                    df_dx_em3 = df_dx_em3 + df_dx_em_tn.*w2_em/2;

                    dfs_dx_rf1 = df_dx_rf1 + dfs_dx_rf_tn.*w1_rf/2;
                    dfs_dx_rf2 = df_dx_rf2 + dfs_dx_rf_tn/2;
                    dfs_dx_rf3 = df_dx_rf3 + dfs_dx_rf_tn.*w2_rf/2;

                    dfs_rf_dpar_rf_fse = dfs_rf_dpar_rf_fse + dfs_rf_dpar_rf_fse_tn;
                    % dfs_rf_dpar_rf_fse2 = dfs_rf_dpar_rf_fse2 + dfs_rf_dpar_rf_fse_tn/2;

                    dfs_rf_dpar_rf_lce = dfs_rf_dpar_rf_lce + dfs_rf_dpar_rf_lce_tn;
                    % dfs_rf_dpar_rf_lce2 = dfs_rf_dpar_rf_lce2 + dfs_rf_dpar_rf_lce_tn/2;

                    dfs_rf_dpar_rf_res0 = dfs_rf_dpar_rf_res0 + dfs_rf_dpar_rf_res;
                end
                
                % extract nonzero elements in the sub jacobian matries,
                % row by row
                row_st = cons_st + cons_st_em + cons_st_rf + (n - ceil(n_em + n_rf) - 1)*M*C;

                for r = 1:M*C
                    
                    if r <= 2*M
                        [row_i, col_i] = find([df_dx1(r, :), df_dx2(r, :)]);

                        row_i = row_i + row_st + r - 1;
                        col_i = col_i + x_stn;

                    elseif r <= 3*M
                        
                        df_dx_x_em = jacDelay_em(n_em, df_dx_em1(2*M+1:3*M, :),...
                            df_dx_em2(2*M+1:3*M, :), df_dx_em3(2*M+1:3*M, :),...
                            df_dx1(2*M+1:3*M, :), df_dx2(2*M+1:3*M, :), M, S);
                         
                        [row_i, col_i] = find(df_dx_x_em(r - 2*M, :));
                                
                        row_i = row_i + row_st + r - 1;
                        col_i = col_i + x_stn - stn_em1;

                    elseif r <= 4*M
                        
                        [row_i, col_i] = find([df_dx1(r, :), df_dx2(r, :)]);

                        row_i = row_i + row_st + r - 1;
                        col_i = col_i + x_stn;

                    else  % 4*M + 1 to 5*M rows
                        
                        % since df_dx_rf only have nonzeros at row 4*M to
                        % 5*M, then only dfs_dx_rf were returned by the
                        % _diff functions. therefore, r_rf is for the
                        % dfs_dx_rf only.
                        
                        % First compute the combined matrix
                        df_dx_em_rf = jacDelay_rf(n_em, n_rf, ...
                            df_dx_em1(4*M + 1:5*M, :), df_dx_em2(4*M + 1:5*M, :), ...
                            df_dx_em3(4*M + 1:5*M, :), dfs_dx_rf1, dfs_dx_rf2, ...
                            dfs_dx_rf3, M, S);
                        
                        % then get the nonzero indexes
                        [row_i1, col_i1] = find(df_dx_em_rf(r - 4*M, :));
                            
                        row_i1 = row_i1 + row_st + r - 1;
                        col_i1 = col_i1 + x_stn - stn_rf1;

                        [row_i3, col_i3] = find(dfs_rf_dpar_rf_fse(r - 4*M, :));
                                             
                        row_i3 = row_i3 + row_st + r - 1;
                        col_i3 = col_i3 + sum(N)*M*S + par_rf_tn;

                        [row_i4, col_i4] = find(dfs_rf_dpar_rf_lce(r - 4*M, :));
                                             
                        row_i4 = row_i4 + + row_st + r - 1;
                        col_i4 = col_i4 + sum(N)*M*S + M*M*tPhase + par_rf_tn;

                        [row_i5, col_i5] = find(dfs_rf_dpar_rf_res0(r - 4*M, :));
                        
                        row_i5 = row_i5 + + row_st + r - 1;
                        col_i5 = col_i5 + sum(N)*M*S + 2*M*M*tPhase;

                        row_i = [row_i1, row_i3, row_i4, row_i5];
                        col_i = [col_i1, col_i3, col_i4, col_i5];

                    end
                    row = [row, row_i];
                    col = [col, col_i];
                end           
            end
        end
    end
end