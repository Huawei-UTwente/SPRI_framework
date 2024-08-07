%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function plot the detail contribution to the muscle activations from
% the reflex loops
%
% By: Huawei Wang
% Date: April 15th, 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fig1, fig2] = reflex_contribution(par_rf, lce, force, act, act_mean, phases, ...
                                            muscle_list, trial_list, N)

    M = length(muscle_list);  % number of muscles
    tPhase = max(phases{1});  % number of phases

    par_rf_for = par_rf(1 : tPhase*M*M);
    par_rf_lce = par_rf(tPhase*M*M + 1 : 2*tPhase*M*M);
    par_rf_res = par_rf(2*tPhase*M*M + 1 : end);
    par_rf_for_thres = par_rf_res(1 : M);
    par_rf_lce_thres = par_rf_res(M + 1 : 2*M);
    par_rf_u0 = par_rf_res(2*M + 1 : 3*M);

    fig1 = figure();
    for m = 1:M
        subplot(2, M, m)
        for t = 1:length(N)
            full_phase = linspace(1, 100, N(t));
            l(t) = plot(full_phase, force(sum(N(1:t-1)) + 1 : sum(N(1:t)), m), '-', 'linewidth', 2);
            hold on
        end
        kp = plot([1, 100], [par_rf_res(m), par_rf_res(m)], 'k--', 'linewidth', 2);
        if m == 1
            ylabel('Force, 1/Fmax')
        end
        title(muscle_list(t))
        if m == M
            legend(l, trial_list(1:length(N)))
            legend(kp, 'threshold')
        end
        hold off

        subplot(2, M, M + m)
        for t = 1:length(N)
            full_phase = linspace(1, 100, N(t));
            plot(full_phase, lce(sum(N(1:t-1)) + 1 : sum(N(1:t)), m), '-', 'linewidth', 2)
            hold on
        end
        if m == 1
            ylabel('Lce, 1/Lce0')
        end
        xlabel('Phase')
        plot([1, 100], [par_rf_res(M + m), par_rf_res(M + m)], 'k--', 'linewidth', 2)
        
        hold off

    end
    hold off

    
    smooth_delta = 1e-8;
    u_force = zeros(sum(N), M*M);
    u_lce = zeros(sum(N), M*M);

    for t = 1:length(N)
        for n = 1:N(t)
            phase_n = phases{t}(n);
            % force reflex control gains at current phase
            par_rf_force_n = par_rf_for((phase_n-1)*M*M + 1 : phase_n*M*M);
            par_rf_lce_n = par_rf_lce((phase_n-1)*M*M + 1 : phase_n*M*M);

            force_n = force(sum(N(1:t - 1)) + n, :);  % force list at node of n
            lce_n = lce(sum(N(1:t - 1)) + n, :);   % lce list at node of n

            s_froce_n = substract_positive_smooth(force_n, par_rf_for_thres, smooth_delta);
            s_lce_n = substract_positive_smooth(lce_n, par_rf_lce_thres, smooth_delta);
            u_force(sum(N(1:t-1)) + n, :) = reshape(reshape(par_rf_force_n, M, M).*s_froce_n', M*M, 1);
            u_lce(sum(N(1:t-1)) + n, :) = reshape(reshape(par_rf_lce_n, M, M).*s_lce_n', M*M, 1);

        end
    end

    fig2 = figure();
    for t = 1:length(N)
        for m = 1:M
            subplot(length(N), M, (t-1)*M + m)
            full_phase = linspace(1, 100, N(t));
            
            plot([0, 100], [par_rf_u0(m), par_rf_u0(m)], 'k-', 'linewidth', 2)
            hold on
            plot(full_phase, act(sum(N(1:t-1)) + 1 : sum(N(1:t)), m), 'r--', 'linewidth', 2)
            hold on
            plot(full_phase, act_mean(sum(N(1:t-1)) + 1 : sum(N(1:t)), m), 'k--', 'linewidth', 2)

            for mm = 1:M
                l(mm) = plot(full_phase, u_force(sum(N(1:t-1)) + 1 : sum(N(1:t)), (m-1)*M + mm), '-');
                hold on
            end

            for mm = 1:M
                plot(full_phase, u_lce(sum(N(1:t-1)) + 1 : sum(N(1:t)), (m-1)*M + mm), '--')
                hold on
            end

            if m == M && t == 1
                legend(l, muscle_list)
            end

            if t == 1
                title(muscle_list(m))
            end
            
            if m == 1
                ylabel(trial_list(t))
            end
            
            if t == length(N)
                xlabel('Phase')
            end

            hold off
        end
    end
end