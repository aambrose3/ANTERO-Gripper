%% 5B3D-Finger-Model-Results: circleCalc
%% Author: Alexander B. Ambrose
%% Date: 8-6-2025

%% Description:
%  Calculates the grasp force for specific circular object size
%  Inputs:
%  param -> basic stucture of fixed finger parameters
%  R -> the circular object radius in m
%  t_3 -> a vector of possible crank/actaotr angles/\theta_3 angles (rad)
%  plotFlag -> toggle for displaying a plot (1 = true)
%
%  Outputs:
%  N_s -> Array of the contact forces: R -> t_3, C -> pads
%  Nx -> vector of the x-direction contact forces
%  Ny -> vector of the x-direction contact forces
%  F -> vector of the total enveloping grasp force (Fx^2 + Fy^2)^1/2

function [N_s, Fx, Fy, F] = circleCalc(param, R, t_3, plotFlag)
    % Circular Object:
    syms t_1_ p1 t_2_ p2 real
    assume(t_1_ > 0 & t_1_ < pi);
    assume(p1 > 0);
    assume(t_2_ > 0 & t_2_ < pi);
    assume(p2 > 0);
    
    %% first contact point
    ep = t_1_ - pi/2;
    eta = R*[cos(ep);sin(ep)];
    Rv = [0; R];
    l = Rv + eta;
    r_p1 = p1*[cos(t_1_);sin(t_1_)] - param.t*[cos(ep);sin(ep)];
    EQ = r_p1 == param.O + l;
    sol= vpasolve([EQ(1), EQ(2)], [t_1_, p1], [pi/10, 0.005]);
    t_1 = double(sol.t_1_);
    p1 = double(sol.p1);
    
    %% second contact point using first solution
    a = param.a*[cos(t_1); sin(t_1)];
    ep = t_1 + t_2_ - pi/2;
    r_p2 = p2*[cos(t_1 + t_2_); sin(t_1 + t_2_)] - param.t*[cos(ep);sin(ep)];
    eta = R*[cos(ep);sin(ep)];
    EQ = a + r_p2 == param.O + Rv + eta;
    sol= vpasolve([EQ(1), EQ(2)], [t_2_, p2], [pi/10, 0.005]);
    t_2 = double(sol.t_2_);
    p2 = double(sol.p2);
    
    %% Calculate contact forces
    % Initialize contact force array
    n3 = numel(t_3);
    N_s = nan(n3, 2);
    
    % Symbolic contact force magnitudes
    syms N1 N2
    assume(N1, {'real', 'positive'});
    assume(N2, {'real', 'positive'});
    
    % Make a “broadcast” copy of your base parameters
    baseParam = param;
    t_1
    p1
    t_2
    p2
    elapsedTime = 0;
    for kk = 1:n3
        % update parameters for this configuration
        localParam = baseParam;
        localParam = updateParam(localParam, t_1, t_2, t_3(kk));
        localParam.p1 = p1; localParam.p2 = p2;
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
    Fx = N_s(:, 1).*cos(localParam.t_1_ - pi/2) + ...
        N_s(:, 2).*cos(localParam.t_2_ + localParam.t_2_ - pi/2);
    Fy = N_s(:, 1).*sin(localParam.t_1_ - pi/2) + ...
        N_s(:, 2).*sin(localParam.t_2_ + localParam.t_2_ - pi/2);
    F = sqrt(Fx.^2 + Fy.^2);
    if plotFlag
        figure
        hold on
        plot(rad2deg(t_3), N_s(:, 1), 'b', 'LineWidth', 2)
        plot(rad2deg(t_3), N_s(:, 2), 'r', 'LineWidth', 2)
        plot(rad2deg(t_3), Fx, ':k', 'LineWidth', 2)
        plot(rad2deg(t_3), Fy, '--k', 'LineWidth', 2)
        plot(rad2deg(t_3), F, 'k', 'LineWidth', 2)
        title('Finger Contact Forces by Pad')
        ylabel('Contact Force (N)')
        xlabel('Actuator Angle (deg)')
        legend('Proximal Contact', 'Distal Contact', ...
            'Total X-Direction', 'Total Y-Direction', ...
            'Total Enveloping Force', 'location', 'northwest')
    end
end
