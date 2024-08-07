% df_dx_em only have nonzero element between the
                        % row 2*M to row 3*M.
                        
                        % depends on how large is the _em delay, x and x_em
                        % maybe constructed using some common data nodes.
                        % They should be combined together when extracting
                        % the nonzero Jac elements.
                        if n_em == 0  % if no electromechancial delay, then same as normal
                            
                            [row_i, col_i] = find([
                                    df_dx_em2(r, :) + df_dx1(r, :), ...
                                    df_dx_em3(r, :) + df_dx2(r, :)]);
                                
                            row_i = row_i + row_st + r - 1;
                            col_i = col_i + x_stn - stn_em1;
                            
                        elseif n_em < 1  % if the _em delay is no more than one data node
                            
                            [row_i, col_i] = find([df_dx_em1(r, :),...
                                    df_dx_em2(r, :) + df_dx1(r, :), ...
                                    df_dx_em3(r, :) + df_dx2(r, :)]);
                                
                            row_i = row_i + row_st + r - 1;
                            col_i = col_i + x_stn - stn_em1;
                            
                        elseif n_em == 1
                            [row_i, col_i] = find([df_dx_em2(r, :),...
                                    df_dx_em3(r, :) + df_dx1(r, :), ...
                                    df_dx2(r, :)]);
                                
                            row_i = row_i + row_st + r - 1;
                            col_i = col_i + x_stn - stn_em1;

                        elseif n_em < 2  % esle if the _em delay is no more than 2 data nodes
                            [row_i, col_i] = find([df_dx_em1(r, :),...
                                    df_dx_em2(r, :),...
                                    df_dx_em3(r, :) + df_dx1(r, :), ...
                                    df_dx2(r, :)]);
                                
                            row_i = row_i + row_st + r - 1;
                            col_i = col_i + x_stn - stn_em1;
                            
                        elseif n_em == 2
                            
                            [row_i, col_i] = find([df_dx_em2(r, :),...
                                    df_dx_em3(r, :),...
                                    df_dx1(r, :), ...
                                    df_dx2(r, :)]);
                                
                            row_i = row_i + row_st + r - 1;
                            col_i = col_i + x_stn - stn_em1;
                            
                                
                        else  % if larger than 2 data nodes
                            
                            if ceil(n_em) == n_em
                            
                                [row_i1, col_i1] = find([
                                                     df_dx_em2(r, :), ...
                                                     df_dx_em3(r, :),]);
                                                 
                                row_i1 = row_i1 + row_st + r - 1;
                                col_i1 = col_i1 + x_stn - stn_em1;
                            
                            else
                                
                                [row_i1, col_i1] = find([df_dx_em1(r, :), ...
                                                     df_dx_em2(r, :), ...
                                                     df_dx_em3(r, :),]);
                                                 
                                row_i1 = row_i1 + row_st + r - 1;
                                col_i1 = col_i1 + x_stn - stn_em1;
                                
                            end

                            [row_i2, col_i2] = find([df_dx1(r, :), df_dx2(r, :)]);
                            
                            row_i2 = row_i2 + row_st + r - 1;
                            col_i2 = col_i2 + x_stn;

                            row_i = [row_i1, row_i2];
                            col_i = [col_i1, col_i2];
                        end
                        
        
                        
 %%                        
                        
                        
                        % only df_dx_em and df_dx_rf have nonzero element 
                        % between the row 4*M to 5*M.
                        
                        % depends on how large is the _rf delay, x_rf and x_em
                        % maybe constructed using some common data nodes.
                        % They should be combined together when extracting
                        % the nonzero Jac elements.
                        
                        % since df_dx_rf only have nonzeros at row 4*M to
                        % 5*M, then only dfs_dx_rf were returned by the
                        % _diff functions. therefore, r_rf is for the
                        % dfs_dx_rf only.
                        
                        r_rf = r - 4*M;
                        
                        if ceil(n_em + n_rf) - ceil(n_em) == 0
                            
                            if n_em == 0 && n_rf == 0
                                
                                [row_i1, col_i1] = find([dfs_dx_rf2(r_rf, :) +  df_dx_em2(r, :), ...
                                        dfs_dx_rf3(r_rf, :) +  df_dx_em3(r, :)]);
                                    
                            else
                                
                                [row_i1, col_i1] = find([dfs_dx_rf1(r_rf, :) + df_dx_em1(r, :),...
                                        dfs_dx_rf2(r_rf, :) +  df_dx_em2(r, :), ...
                                        dfs_dx_rf3(r_rf, :) +  df_dx_em3(r, :)]);
                                    
                            end
                                    
                            row_i1 = row_i1 + row_st + r - 1;
                            col_i1 = col_i1 + x_stn - stn_rf1;

                            row_i2 = [];
                            col_i2 = [];
                        
                        elseif ceil(n_em + n_rf) - ceil(n_em) == 1  % if _rf and _em delay is one node larger than _em delay
                            
                            if ceil(n_em) == n_em && ceil(n_em + n_rf) == n_em + n_rf
                                
                                [row_i1, col_i1] = find([ ...
                                    dfs_dx_rf2(r_rf, :), ...
                                    dfs_dx_rf3(r_rf, :) +  df_dx_em2(r, :),...
                                    df_dx_em3(r, :)]);
                                
                            elseif ceil(n_em) == n_em && ceil(n_em + n_rf) ~= n_em + n_rf
                                
                                [row_i1, col_i1] = find([ ...
                                    dfs_dx_rf1(r_rf, :), ...
                                    dfs_dx_rf2(r_rf, :) + df_dx_em2(r, :),...
                                    dfs_dx_rf3(r_rf, :) + df_dx_em3(r, :)]);
                                
                            elseif ceil(n_em) ~= n_em && ceil(n_em + n_rf) == n_em + n_rf
                                
                                [row_i1, col_i1] = find([ ...
                                    dfs_dx_rf1(r_rf, :), ...
                                    dfs_dx_rf2(r_rf, :) + df_dx_em1(r, :),...
                                    df_dx_em2(r, :), ...
                                    df_dx_em3(r, :)]);
                                
                            else
                                
                                [row_i1, col_i1] = find([dfs_dx_rf1(r_rf, :), ...
                                    dfs_dx_rf2(r_rf, :) +  df_dx_em1(r, :),...
                                    dfs_dx_rf3(r_rf, :) +  df_dx_em2(r, :), ...
                                    df_dx_em3(r, :)]);
                                
                            end
                                
                            row_i1 = row_i1 + row_st + r - 1;
                            col_i1 = col_i1 + x_stn - stn_rf1;

                            row_i2 = [];
                            col_i2 = [];
                            
                        elseif ceil(n_em + n_rf) - ceil(n_em) == 2
                            if ceil(n_em) == n_em && ceil(n_em + n_rf) == n_em + n_rf
                                
                                [row_i1, col_i1] = find([dfs_dx_rf2(r_rf, :), ...
                                    dfs_dx_rf2(r_rf, :),...
                                    df_dx_em2(r, :), ...
                                    df_dx_em3(r, :)]);
                                
                            elseif ceil(n_em) == n_em && ceil(n_em + n_rf) ~= n_em + n_rf
                                
                                [row_i1, col_i1] = find([dfs_dx_rf1(r_rf, :), ...
                                    dfs_dx_rf2(r_rf, :),...
                                    dfs_dx_rf3(r_rf, :) + df_dx_em2(r, :), ...
                                    df_dx_em3(r, :)]);
                                
                            elseif ceil(n_em) ~= n_em && ceil(n_em + n_rf) == n_em + n_rf
                                
                                [row_i1, col_i1] = find([dfs_dx_rf2(r_rf, :), ...
                                    dfs_dx_rf3(r_rf, :),...
                                    df_dx_em1(r, :), ...
                                    df_dx_em2(r, :), ...
                                    df_dx_em3(r, :)]);
                                
                            else
                                [row_i1, col_i1] = find([dfs_dx_rf1(r_rf, :), ...
                                        dfs_dx_rf2(r_rf, :),...
                                        dfs_dx_rf3(r_rf, :) +  df_dx_em1(r, :),...
                                        df_dx_em2(r, :), ...
                                        df_dx_em3(r, :)]);
                            end
                                
                            row_i1 = row_i1 + row_st + r - 1;
                            col_i1 = col_i1 + x_stn - stn_rf1;

                            row_i2 = [];
                            col_i2 = [];
                            
                        else
                            if ceil(n_em + n_rf) == n_em + n_rf
                                [row_i1, col_i1] = find([dfs_dx_rf2(r_rf, :), ...
                                                     dfs_dx_rf3(r_rf, :)]);
                            else
                                
                                [row_i1, col_i1] = find([dfs_dx_rf1(r_rf, :), ...
                                                     dfs_dx_rf2(r_rf, :), ...
                                                     dfs_dx_rf3(r_rf, :)]);
                            end
                                                 
                            row_i1 = row_i1 + row_st + r - 1;
                            col_i1 = col_i1 + x_stn - stn_rf1;
                            
                            if ceil(n_em) == n_em
                                [row_i2, col_i2] = find([df_dx_em2(r, :), ....
                                                     df_dx_em3(r, :)]);
                            else
                                [row_i2, col_i2] = find([df_dx_em1(r, :), ....
                                                        df_dx_em2(r, :), ....
                                                     df_dx_em3(r, :)]);
                            end
                                                 
                            row_i2 = row_i2 + row_st + r - 1;
                            col_i2 = col_i2 + x_stn - stn_em1;
                        
                        end
                        
                        
                        %%
                        
                                                        
                        row_i = row_i + row_st + r - 1;
                        col_i = col_i + x_stn - stn_em1;

                        
                        % df_dx_em only have nonzero element between the
                        % row 2*M to row 3*M.
                        
                        % depends on how large is the _em delay, x and x_em
                        % maybe constructed using some common data nodes.
                        % They should be combined together when extracting
                        % the nonzero Jac elements.
                        
                        if n_em == 0  % if no electromechancial delay, then same as normal
                            
                            Jac = [Jac, nonzeros([
                                    df_dx_em2(r, :) + df_dx1(r, :), ...
                                    df_dx_em3(r, :) + df_dx2(r, :)])'];
                        
                        
                        elseif n_em < 1

                            Jac = [Jac, nonzeros([df_dx_em1(r, :), ...
                                    df_dx_em2(r, :) + df_dx1(r, :), ...
                                    df_dx_em3(r, :) + df_dx2(r, :)])'];
                                
                        elseif n_em == 1
                            
                            Jac = [Jac, nonzeros([df_dx_em2(r, :), ...
                                    df_dx_em3(r, :) + df_dx1(r, :), ...
                                    df_dx2(r, :)])'];
                                
                        elseif n_em < 2
                            
                            Jac = [Jac, nonzeros([df_dx_em1(r, :), ...
                                    df_dx_em2(r, :), ...
                                    df_dx_em3(r, :) + df_dx1(r, :), ...
                                    df_dx2(r, :)])'];
                                
                        elseif n_em == 2
                            
                            Jac = [Jac, nonzeros([df_dx_em2(r, :), ...
                                    df_dx_em3(r, :), ...
                                    df_dx1(r, :), ...
                                    df_dx2(r, :)])'];

                        else
                            
                            if ceil(n_em) == n_em
                                
                                Jac = [Jac, nonzeros([
                                    df_dx_em2(r, :), ...
                                    df_dx_em3(r, :), ...
                                    df_dx1(r, :), ...
                                    df_dx2(r, :)])'];
                                
                            else
                            
                                Jac = [Jac, nonzeros([df_dx_em1(r, :), ...
                                    df_dx_em2(r, :), ...
                                    df_dx_em3(r, :), ...
                                    df_dx1(r, :), ...
                                    df_dx2(r, :)])'];
                            end
                                
                        end
                        
                        %%
                        
                        % only df_dx_em and df_dx_rf have nonzero element 
                        % between the row 4*M to 5*M.
                        
                        % depends on how large is the _rf delay, x_rf and x_em
                        % maybe constructed using some common data nodes.
                        % They should be combined together when extracting
                        % the nonzero Jac elements.
                        
                        % since df_dx_rf only have nonzeros at row 4*M to
                        % 5*M, then only dfs_dx_rf were returned by the
                        % _diff functions. therefore, r_rf is for the
                        % dfs_dx_rf only.
                        r_rf = r - 4*M;  
                        
                        if ceil(n_em + n_rf) - ceil(n_em) == 0
                            
                            Jac = [Jac, nonzeros([dfs_dx_rf1(r_rf, :) + df_dx_em1(r, :),...
                                    dfs_dx_rf2(r_rf, :) +  df_dx_em2(r, :), ...
                                    dfs_dx_rf3(r_rf, :) +  df_dx_em3(r, :), ...
                                    dfs_rf_dpar_rf_fse_tn(r_rf, :), ...
                                    dfs_rf_dpar_rf_lce_tn(r_rf, :), ...
                                    dfs_rf_dpar_rf_res(r_rf, :)])'];
                      
                        elseif ceil(n_em + n_rf) - ceil(n_em) == 1
                            
                            Jac = [Jac, nonzeros([dfs_dx_rf1(r_rf, :), ...
                                    dfs_dx_rf2(r_rf, :) +  df_dx_em1(r, :),...
                                    dfs_dx_rf3(r_rf, :) +  df_dx_em2(r, :), ...
                                    df_dx_em3(r, :), ...
                                    dfs_rf_dpar_rf_fse_tn(r_rf, :), ...
                                    dfs_rf_dpar_rf_lce_tn(r_rf, :),...
                                    dfs_rf_dpar_rf_res(r_rf, :)])'];
                      
                        elseif ceil(n_em + n_rf) - ceil(n_em) == 2
                            
                            Jac = [Jac, nonzeros([dfs_dx_rf1(r_rf, :), ...
                                    dfs_dx_rf2(r_rf, :),...
                                    dfs_dx_rf3(r_rf, :) +  df_dx_em1(r, :),...
                                    df_dx_em2(r, :), ...
                                    df_dx_em3(r, :), ...
                                    dfs_rf_dpar_rf_fse_tn(r_rf, :), ...
                                    dfs_rf_dpar_rf_lce_tn(r_rf, :),...
                                    dfs_rf_dpar_rf_res(r_rf, :)])'];
                      
                        else
                            
                            Jac = [Jac, nonzeros([dfs_dx_rf1(r_rf, :), ...
                                    dfs_dx_rf2(r_rf, :), ...
                                    dfs_dx_rf3(r_rf, :), ...
                                    df_dx_em1(r, :), ...
                                    df_dx_em2(r, :), ...
                                    df_dx_em3(r, :), ...
                                    dfs_rf_dpar_rf_fse_tn(r_rf, :), ...
                                    dfs_rf_dpar_rf_lce_tn(r_rf, :),...
                                    dfs_rf_dpar_rf_res(r_rf, :)])'];
                      
                        end