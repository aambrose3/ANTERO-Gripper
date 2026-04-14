%% 5B3D-Finger-Model: Main
%% Author: Alexander B. Ambrose
%% Date: 8-6-2025

%% Description:
% This script solves the 5-Bar-3-DOF finger static model for contact 
% force estimation. A parallel pool of 8 workers can solve this model in
% about 2-minutes. The output contact force (N) has two columns, the first
% column is the contact force produced by the proximal and the second 
% column is for the distal phalange.

%% Main
clear all; close all; clc;
warning('off', 'optim:fsolve:NonSquareSystem');

% Load fixed finger paramters
param = getParam();

%% Sweep from hardstop to hardstop of finger states (from CAD)
% Sweep over t_1 every 0.25 degrees
% t_1 = linspace(deg2rad(45), deg2rad(105), (105-45)+1)'; 
t_1 = (1.2335:0.01:1.2535)';
% Sweep over t_2 every 0.25 degrees
% t_2 = linspace(deg2rad(20), deg2rad(105), (105-20)+1)'; 
t_2 = (1.0961:0.01:1.1161)';
% Sweep the actuator angle over its entire range every 0.25 degrees
% t_3 = linspace(deg2rad(-90), deg2rad(40), (90+40)+1)';
t_3 = linspace(-0.6, -0.2, 401);

%% Parfor loop prepeartions
% Preallocate output and timing
n1 = numel(t_1);
n2 = numel(t_2);
n3 = numel(t_3);
N = nan(n1, n2, n3, 2);   % result array
times = zeros(n1, 1);         % per‐ii timing

% Symbolic contact force magnitudes
syms N1 N2
assume(N1, {'real', 'positive'});
assume(N2, {'real', 'positive'});

% Make a “broadcast” copy of your base parameters
baseParam = param;

%% Loop through possible actuator angles and finger states and calculate 
%  the contact forces
parfor ii = 1:n1
    % each worker gets its own copy
    localParam = baseParam;
    elapsedTime = 0;
    for jj = 1:n2
        for kk = 1:n3
            tic;
            % update parameters for this configuration
            localParam = updateParam(localParam, t_1(ii), t_2(jj), t_3(kk));
            % compute kinematics
            Jac = getJacobians(localParam, t_1(ii), t_2(jj), t_3(kk));
            % compute spring deflection & force
            [cmp, F_cmp, TA, flag] = getSpring(localParam, t_1(ii), t_2(jj), t_3(kk));
            if flag == 0
                % build torque balance
                Tau_contacts = Jac.J1'*localParam.n1*N1 + Jac.J2'*localParam.n2*N2;
                Tau_h = [0; localParam.K2*(t_2(jj) - localParam.l2_r)];
                Tau_C = Jac.JC'*localParam.nC*F_cmp;
                EQ = Tau_contacts + Tau_h + Tau_C == 0;
                if any(isnan(Tau_C)) == 1 || any(isnan(Tau_contacts)) == 1 || any(isnan(Tau_h)) == 1
                    N(ii, jj, kk, :) = [NaN, NaN];
                else
                    [A, B] = equationsToMatrix([EQ(1), EQ(2)], [N1, N2]);
                    sol    = vpa(linsolve(A, B))';  %# solve for [N1;N2]
                    % discard unstable solutions
                    f1 = localParam.n1*sol(1);
                    f2 = localParam.n2*sol(2);
                    if sol(1) < 0
                        N(ii, jj, kk, :) = [0, sol(2)];
                    end
                    if sol(2) < 0
                        N(ii, jj, kk, :) = [sol(1), 0];
                    end
                    if f1(2) + f2(2) < 0
                        N(ii, jj, kk, :) = [NaN, NaN];
                    else
                        N(ii, jj, kk, :) = sol;
                    end
                end
            elseif flag == 1
                % spring out of range → stop kk‐loop
                fprintf('\nii=%d, jj=%d → Done\n', ii, jj);
                break;
            end
            elapsedTime = elapsedTime + toc;
            % fprintf('ii=%d, jj=%d → elapsed: %5.2f s\n', ii, jj, elapsedTime);
        end
    end
    times(ii) = elapsedTime;
    % Print out ii-loop timings
    fprintf('\nii=%d done (total %5.3f s)\n', ii, elapsedTime);
end
% optional: report overall time
fprintf('\nTotalotal time: %5.1f s\n', sum(times));

% save the results
save('Results/Results_t1_t2.mat', 'N', 't_1', 't_2', 't_3', 'param'); % will overwirte if necessary