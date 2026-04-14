%% 5B3D-Finger-Model-Results: ellipseCalc
%% Author: Alexander B. Ambrose
%% Date: 8-6-2025

%% Description:
%  Calculates the grasp force for specific circular object size
%  Inputs:
%  param -> basic stucture of fixed finger parameters
%  t_1 -> scalar or array of test angles t_1
%  t_2 -> scalar or array of test angles t_2
%  t_3 -> a vector of possible crank/actaotr angles/\theta_3 angles (rad)
%  plotFlag -> toggle for displaying a plot (1 = true)
%
%  Outputs:
%  N_s -> Array of the contact forces: Row -> t_3, Col -> pads (min and max
%  values from test angles t_1 and t_2 are stored in the depth dimension)
%  Nx -> vector of the x-direction contact forces
%  Ny -> vector of the x-direction contact forces
%  F -> vector of the total enveloping grasp force (Fx^2 + Fy^2)^1/2

function [N_s, Fx, Fy, F] = ellipseCalc(param, t_1, t_2, t_3, plotFlag)   
    warning('off', 'all')
    if length(t_1) == 1 && length(t_2) == 1
        [N_s, Fx, Fy, F] = Compute(param, t_1, t_2, t_3)
    else
        for ii = 1:length(t_1)
            for jj = 1:length(t_2)
                if ii == 1 && jj == 1
                    [N_s, Fx, Fy, F] = Compute(param, t_1(ii), t_2(jj), t_3);                  
                    N_min = N_s;
                    N_max = N_s;
                    Fx_min = Fx;
                    Fx_max = Fx;
                    Fy_min = Fy;
                    Fy_max = Fy;
                    F_min = F;
                    F_max = F;
                else
                    [n_s, fx, fy, f] = Compute(param, t_1(ii), t_2(jj), t_3);
                    N_min(:, 1) = min(n_s(:, 1), N_min(:, 1));
                    N_min(:, 2) = min(n_s(:, 2), N_min(:, 2));
                    N_max(:, 1) = max(n_s(:, 1), N_max(:, 1));
                    N_max(:, 2) = max(n_s(:, 2), N_max(:, 2));
                    Fx_min = min(fx, Fx_min);
                    Fx_max = max(fx, Fx_max);
                    Fy_min = min(fy, Fy_min);
                    Fy_max = max(fy, Fy_max);
                    F_min = min(f, F_min);
                    F_max = max(f, F_max);
                end
                fprintf('\nii=%d/%d, jj=%d/%d→ Done\n', ii, length(t_1), jj, length(t_2));
            end
        end
        N_s = cat(3, N_min, N_max);
        Fx = cat(2, Fx_min, Fx_max);
        Fy = cat(2, Fy_min, Fy_max);
        F = cat(2, F_min, F_max);
    end
    if plotFlag        
        if length(t_1) > 1 || length(t_2) > 1
            figure
            hold on
            high = [t_3, Fx_max];
            low = [t_3, Fx_min];
            comb = [high;flipud(low)];
            if any(isnan(comb(:, 2)))
                idx = find(isnan(comb(:, 2)) == 1);
                comb(idx, 2) = 0;
            end
            if any(comb(5:end, 2)) < 0
                idx = find(comb(5:end, 2) < 0, 1, 'first')+4;
                comb = comb(1:idx-1, :); comb = [comb; comb(1, :)];
            end
            fill(comb(:, 1), comb(:, 2), ...
                'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none')

            high = [t_3, Fy_max];
            low = [t_3, Fy_min];
            comb = [high;flipud(low)];
            if any(isnan(comb(:, 2)))
                idx = find(isnan(comb(:, 2)) == 1);
                comb(idx, 2) = 0;
            end
            if any(comb(5:end, 2)) < 0
                idx = find(comb(5:end, 2) < 0, 1, 'first')+4;
                comb = comb(1:idx-1, :); comb = [comb; comb(1, :)];
            end
            fill(comb(:, 1), comb(:, 2), ...
                'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none')

            high = [t_3, F_max];
            low = [t_3, F_min];
            comb = [high;flipud(low)];
            if any(isnan(comb(:, 2)))
                idx = find(isnan(comb(:, 2)) == 1);
                comb(idx, 2) = 0;
            end
            if any(comb(5:end, 2)) < 0
                idx = find(comb(5:end, 2) < 0, 1, 'first')+4;
                comb = comb(1:idx-1, :); comb = [comb; comb(1, :)];
            end
            fill(comb(:, 1), comb(:, 2), ...
                'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
        else
            figure
            hold on
            plot(t_3, mean(N_s(:, 1, :), 3), 'b', 'LineWidth', 2)
            plot(t_3, mean(N_s(:, 2, :), 3), 'r', 'LineWidth', 2)
            plot(t_3, mean(Fx, 2), ':k', 'LineWidth', 2)
            plot(t_3, mean(Fy, 2), '--k', 'LineWidth', 2)
            plot(t_3, mean(F, 2), 'k', 'LineWidth', 2)
            title('Finger Contact Forces by Pad')
            ylabel('Contact Force (N)')
            xlabel('Actuator Angle (deg)')
            legend('Proximal Contact', 'Distal Contact', ...
                'Total X-Direction', 'Total Y-Direction', ...
                'Total Enveloping Force', 'location', 'northwest')
        end
    end
end

function [N_s, Fx, Fy, F] = Compute(param, t_1_, t_2_, t_3)
    t_1 = t_1_;
    t_2 = t_2_;
    n3 = numel(t_3);
    N_s = nan(n3, 2);
    
    % Symbolic contact force magnitudes
    syms N1 N2
    assume(N1, {'real', 'positive'});
    assume(N2, {'real', 'positive'});
    
    % Make a “broadcast” copy of your base parameters
    baseParam = param;
    elapsedTime = 0;
    for kk = 1:n3
        % update parameters for this configuration
        localParam = baseParam;
        localParam = updateParam(localParam, t_1, t_2, t_3(kk));
        % compute kinematics
        Jac = getJacobians(localParam, localParam.t_1_, localParam.t_2_, t_3(kk));
        % compute spring deflection & force
        [cmp, F_cmp, TA, flag] = getSpring(localParam, localParam.t_1_, localParam.t_2_, t_3(kk));
        if flag == 0
            % build torque balance
            Tau_contacts = Jac.J1'*localParam.n1*N1 + Jac.J2'*localParam.n2*N2;
            Tau_h = [0; localParam.K2*(localParam.t_2_ - localParam.l2_r)];
            Tau_C = Jac.JC'*localParam.nC*F_cmp;
            EQ = Tau_contacts + Tau_h + Tau_C == 0;
            if any(isnan(Tau_C)) == 1 || any(isnan(Tau_contacts)) == 1 || any(isnan(Tau_h)) == 1
                N_s(kk, :) = [NaN, NaN];
            else
                [A, B] = equationsToMatrix([EQ(1), EQ(2)], [N1, N2]);
                sol    = vpa(linsolve(A, B))';  %# solve for [N1;N2]
                % discard unstable solutions
                f1 = localParam.n1*sol(1);
                f2 = localParam.n2*sol(2);
                if sol(1) < 0
                    N_s(kk, 1) = 0;
                end
                if sol(2) < 0
                    N_s(kk, 2) = 0;
                end
                if f1(2) + f2(2) < 0
                    N_s(kk, :) = [NaN, NaN];
                else
                    N_s(kk, :) = sol;
                end
            end
        elseif flag == 1
            % spring out of range → stop kk‐loop
            fprintf('\nkk=%d → Done\n', kk);
            break;
        end
        if kk == 801
            F_cmp;
        end
    end
    % Cartesian Forces imparted onto the finger pads by the object
    Fx = N_s(:, 1).*localParam.n1(1) + ...
        N_s(:, 2).*localParam.n2(1);
    Fy = N_s(:, 1).*localParam.n1(2) + ...
        N_s(:, 2).*localParam.n2(2);
    F = sqrt(Fx.^2 + Fy.^2);
end

