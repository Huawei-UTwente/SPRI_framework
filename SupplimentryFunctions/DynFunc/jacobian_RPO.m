function Jac = jacobian_RPO(x, M, S, C, N, nPar_rf, lmt, par_mus, t_em, t_rf, hs, phase, tPhase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function to calculate Jacobian of the optimization problem. Only
% nonzero elements will be saved.
%
% By: Huawei Wang
% Date: March 13, 2022

% INPUTS:
%   x, vector: state variables that need to be optimized
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
%   Jac, vector: the nonzero vector of the jacobian matrix of the
%   optimization problem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Extract the fse, lce reflex control gains and other reflex control
    % parameters.
    par_rf = x(end-nPar_rf+1:end);
    par_rf_fse = par_rf(1:M*M*tPhase);
    par_rf_lce = par_rf(M*M*tPhase+1:2*M*M*tPhase);
    par_rf_res = par_rf(2*M*M*tPhase+1:end);
    
    % Initilize the empty nonzero Jacobian vector
    Jac = [];

    for t = 1:length(N)      % Run through all data trials
        
        hs_t = hs(t);        % Time interval between the data nodes
        n_em = t_em(t)/hs_t; % Number of data nodes that caused by electromechanical delay (_em)
        n_rf = t_rf(t)/hs_t; % Number of data nodes that caused by reflex delay (_rf)

        x_st = sum(N(1:t-1))*M*S; % The index of the first state parameter of current data trial 
        lmt_st = sum(N(1:t-1));   % The index of the first lmt data of current data trial
        
        % Looking back to the data nodes that caused by _em and _rf delays 
        stn_em1 = ceil(n_em)*M*S;
        stn_em2 = floor(n_em)*M*S;
        
        stn_rf1 = ceil(n_em + n_rf)*M*S;
        stn_rf2 = floor(n_em + n_rf)*M*S;
        
        for n = 1: N(t)-1   % Run through all data nodes (using middle point method)
            
            x_stn = x_st + (n-1)*M*S; 	% The index right before the first state of current data node
            x_stn1 = x_st + n*M*S;      % The index of the last state of current data node
            x_stn2 = x_st + (n+1)*M*S;  % The index of last state of next data node
            
            % Calculate the x, dx, and lmt based on the middle point method
            x_tn = (x(x_stn1+1: x_stn2) + x(x_stn+1 : x_stn1))/2;
            dx_tn = (x(x_stn1+1: x_stn2) - x(x_stn+1 : x_stn1))./hs_t;
            
            lmt_tn = (lmt(lmt_st + n, :) + lmt(lmt_st + n + 1, :))/2; 
            
            % If the number of current data node is smaller than the data node
            % required by the _em delay. Then the activation
            % dynamic and reflex control constraints are not included,
            % since there are no state can be used to formulate these
            % feedbacks.
            if n <= ceil(n_em)
                
                % directCollocationDyn_diff_RPO1 is the derivative function
                % for directCollocationDyn_RPO1. Only df_dx, and fx_ddx are
                % calculate, since x_em and x_rf are not available
                [df_dx_tn, df_ddx_tn] = ...
                directCollocationDyn_diff_RPO1(x_tn, dx_tn, ...
                    lmt_tn, par_mus, M, S);
                
                % (middle point method) assign back derivatives to data nodes
                df_dx1 = df_dx_tn./2 - df_ddx_tn./hs_t;
                df_dx2 = df_dx_tn./2 + df_ddx_tn./hs_t;

                % extract nonzero elements in the jacobian matrix, row by
                % row
                for r = 1:M*(C-2)
                    Jac = [Jac, nonzeros([df_dx1(r, :), df_dx2(r, :)])'];
                end
                
            % Else if the number of current data node is larger than the
            % _em, but smaller than _em + _rf. Then the activation dynamics
            % is included, the reflex control constraint will not be 
            % included.            
            elseif n <= ceil(n_em + n_rf)
            
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
                
                % _em delayed data node of x1
                x_em1 = x(x_stn-stn_em1+1 : x_stn1-stn_em1).*w1_em  ...
                        + x(x_stn-stn_em2+1 : x_stn1-stn_em2).*w2_em;
                
                % _em delayed data node of x2
                x_em2 = x(x_stn1-stn_em1+1 : x_stn2-stn_em1).*w1_em  ...
                        + x(x_stn1-stn_em2+1 : x_stn2-stn_em2).*w2_em;
                
                % Middle point method to get the _em delayed variable for
                % the activation dynamics
                x_em_tn = (x_em1 + x_em2)/2;
               
                % directCollocationDyn_diff_RPO2 is the derivative function
                % for directCollocationDyn_RPO2. df_dx_em is included here,
                % however, df_dx_rf is still not included, since x_rf is
                % not available yet.
                [df_dx_tn, df_ddx_tn, df_dx_em_tn] = ...
                directCollocationDyn_diff_RPO2(x_tn, dx_tn, x_em_tn, ...
                    lmt_tn, par_mus, M, S);
                
                % (middle point method) assign back derivatives to data nodes
                df_dx1 = df_dx_tn./2 - df_ddx_tn./hs_t;
                df_dx2 = df_dx_tn./2 + df_ddx_tn./hs_t;

                df_dx_em1 = df_dx_em_tn.*w1_em/2;
                df_dx_em2 = df_dx_em_tn.*(w1_em + w2_em)/2;
                df_dx_em3 = df_dx_em_tn.*w2_em/2;

                % extract nonzero elements in the jacobian matrix, row by
                % row
                for r = 1:M*(C-1)
                    if r <= 2*M
                        
                        Jac = [Jac, nonzeros([df_dx1(r, :), df_dx2(r, :)])'];

                    elseif r <= 3*M
                        
                        df_dx_x_em = jacDelay_em(n_em, df_dx_em1(2*M+1:3*M, :),...
                            df_dx_em2(2*M+1:3*M, :), df_dx_em3(2*M+1:3*M, :),...
                            df_dx1(2*M+1:3*M, :), df_dx2(2*M+1:3*M, :), M, S);
                         
                        Jac = [Jac, nonzeros(df_dx_x_em(r - 2*M, :))'];
  
                    else
                        
                        Jac = [Jac, nonzeros([df_dx1(r, :), df_dx2(r, :)])'];

                    end
                end       
              
            % Else when there are states can be used to formulate the
            % reflex control loops, include the reflex control constraints.
            else
                
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
                    w1_rf = n_rf + n_em - floor(n_em + n_rf); % weight of the next data node
                    w2_rf = ceil(n_em + n_rf) - n_rf - n_em;  % weight of the perior data node
                end
                
                % Same as used above
                x_em1 = x(x_stn-stn_em1+1 : x_stn1-stn_em1).*w1_em  ...
                        + x(x_stn-stn_em2+1 : x_stn1-stn_em2).*w2_em;
                    
                x_em2 = x(x_stn1-stn_em1+1 : x_stn2-stn_em1).*w1_em  ...
                        + x(x_stn1-stn_em2+1 : x_stn2-stn_em2).*w2_em;
                
                x_em_tn = (x_em1 + x_em2)/2;
                
                % Similar to x_em reconstruction, calculate x_rf
                x_rf1 = x(x_stn-stn_rf1+1 : x_stn1-stn_rf1).*w1_rf ...
                       + x(x_stn-stn_rf2+1 : x_stn1-stn_rf2).*w2_rf;
                   
                x_rf2 = x(x_stn1-stn_rf1+1 : x_stn2-stn_rf1).*w1_rf ...
                       + x(x_stn1-stn_rf2+1 : x_stn2-stn_rf2).*w2_rf;

                x_rf_tn = (x_rf1 + x_rf2)/2;
                
                % Similarly, calculate the lmt at the corresponding data position
                % of x_rf.
                lmt_rf1 = lmt(lmt_st + n - ceil(n_em+n_rf), :)*w1_rf ...
                       + lmt(lmt_st + n - floor(n_em+n_rf), :)*w2_rf;
                   
                lmt_rf2 = lmt(lmt_st + n + 1 - ceil(n_em+n_rf), :)*w1_rf ...
                       + lmt(lmt_st + n + 1 - floor(n_em+n_rf), :)*w2_rf;

                lmt_rf_tn = (lmt_rf1 + lmt_rf2)/2;

                % Reflex control parameters need to be reconstructed.
                % However, reflex control gains should be extracted based
                % on the phase of current data node.
                % par_rf_tn = (phase{t}(n-1))*M*M;
                par_rf_tn1 = (phase{t}(n)-1)*M*M;
                par_rf_tn2 = (phase{t}(n))*M*M;

                par_rf_fse_tn = par_rf_fse(par_rf_tn1+1 : par_rf_tn2);
                par_rf_lce_tn = par_rf_lce(par_rf_tn1+1 : par_rf_tn2);
                    
                % directCollocationDyn_diff_RPO3 is the derivative function
                % for directCollocationDyn_RPO3. All derivatives are
                % included.
                [df_dx_tn, df_ddx_tn, df_dx_em_tn, dfs_dx_rf_tn, dfs_rf_dpar_rf_fse_tn, ...
                 dfs_rf_dpar_rf_lce_tn, dfs_rf_dpar_rf_res] = ...
                 directCollocationDyn_diff_RPO3(x_tn, dx_tn, ...
                        x_em_tn, x_rf_tn, lmt_tn, lmt_rf_tn, par_rf_fse_tn, ...
                        par_rf_lce_tn, par_rf_res, par_mus, M, S);

                % (middle point method) assign back derivatives to data nodes
                df_dx1 = df_dx_tn./2 - df_ddx_tn./hs_t;
                df_dx2 = df_dx_tn./2 + df_ddx_tn./hs_t;

                df_dx_em1 = df_dx_em_tn.*w1_em/2;
                df_dx_em2 = df_dx_em_tn/2;
                df_dx_em3 = df_dx_em_tn.*w2_em/2;

                dfs_dx_rf1 = dfs_dx_rf_tn.*w1_rf/2;
                dfs_dx_rf2 = dfs_dx_rf_tn/2;
                dfs_dx_rf3 = dfs_dx_rf_tn.*w2_rf/2;

%                 dfs_rf_dpar_rf_fse1 = dfs_rf_dpar_rf_fse_tn/2;
%                 dfs_rf_dpar_rf_fse2 = dfs_rf_dpar_rf_fse_tn/2;
% 
%                 dfs_rf_dpar_rf_lce1 = dfs_rf_dpar_rf_lce_tn/2;
%                 dfs_rf_dpar_rf_lce2 = dfs_rf_dpar_rf_lce_tn/2;
            
                % extract nonzero elements in the sub jacobian matries,
                % row by row
                for r = 1:M*C
                    
                    if r <= 2*M
                        
                        Jac = [Jac, nonzeros([df_dx1(r, :), df_dx2(r, :)])'];

                    elseif r <= 3*M
                        
                        df_dx_x_em = jacDelay_em(n_em, df_dx_em1(2*M+1:3*M, :),...
                            df_dx_em2(2*M+1:3*M, :), df_dx_em3(2*M+1:3*M, :),...
                            df_dx1(2*M+1:3*M, :), df_dx2(2*M+1:3*M, :), M, S);
                         
                        Jac = [Jac, nonzeros(df_dx_x_em(r - 2*M, :))'];

                    elseif r <= 4*M
                        
                        Jac = [Jac, nonzeros([df_dx1(r, :), df_dx2(r, :)])'];

                    else
                        
                        df_dx_em_rf = jacDelay_rf(n_em, n_rf, ...
                            df_dx_em1(4*M + 1:5*M, :), df_dx_em2(4*M + 1:5*M, :), ...
                            df_dx_em3(4*M + 1:5*M, :), dfs_dx_rf1, dfs_dx_rf2, ...
                            dfs_dx_rf3, M, S);
                        
                        Jac = [Jac, nonzeros([df_dx_em_rf(r - 4*M, :), ...
                                        dfs_rf_dpar_rf_fse_tn(r - 4*M, :), ...
                                        dfs_rf_dpar_rf_lce_tn(r - 4*M, :),...
                                        dfs_rf_dpar_rf_res(r - 4*M, :)])'];
                        
                    end
                end
            end
        end
    end
end