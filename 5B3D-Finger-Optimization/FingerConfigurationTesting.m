%% Run a sweep over t3 at a specific finger configuration

%% T2 K2 = 0.02
% F =
   % -4.6922
   % -0.7556
    % 0.2998
% M =
    % 0.0197
   % -0.0700
    % 0.1304

%% T3 K2 = 0.011
% F =
   % -3.4507
   % -0.2876
    % 0.1063
% M =
    % 0.0051
   % -0.0466
    % 0.1007

%% T4 K2 = 0.0353 F = 5
% F =
   % -3.7545
   % -0.4253
    % 0.2148
% M =
    % 0.0101
   % -0.0536
    % 0.1022

%% T5 K2 = 0.0353 F = 10
% F =
   % -4.0855
   % -0.8384
    % 0.2742
% M =
    % 0.0260
   % -0.0586
    % 0.1076
clear all; close all; clc;
%% T2
% th1 = 1.045; th2 = 1.0472;
%% T3
% th1 = 1.0407; th2 = 1.0363;
%% T4
th1 = 1.0385; th2 = 1.0450;
%% T5
% th1 = 1.0363; th2 = 1.0559;

th3 = linspace(deg2rad(-50), deg2rad(-35), 151); % every 0.1 degrees
th3 =th3(:); % enforce column vectors
[tau, N1, N2, m] = updateFinger(th1, th2, th3);

figure(1)
hold on
plot(th3, N1, 'b', 'LineWidth', 2)
plot(th3, N2, 'r', 'LineWidth', 2)
ylabel('Normal Contact Forces (N)')
xlabel('Actuator Crank Angle (deg)')
legend('N1', 'N2', 'location', 'northwest')


figure(2)
hold on
plot(th3, N1*cos(th1+pi/2) + N2*cos(th1+th2+pi/2), 'b', 'LineWidth', 2)
plot(th3, N1*sin(th1+pi/2) + N2*sin(th1+th2+pi/2), 'r', 'LineWidth', 2)
ylabel('Cartesian Contact Forces (N)')
xlabel('Actuator Crank Angle (deg)')
legend('Nx', 'Ny', 'location', 'northwest')

%% T2
% th3 = -0.7333;
%% T3
% th3 = -0.7446;
%% T4
th3 = -0.7483;
%% T5
% th3 = -0.7139;

[tau, N1, N2, m] = updateFinger(th1, th2, th3);
N1
Fx = N1*cos(th1+pi/2) + N2*cos(th1+th2+pi/2)
N2
Fy = N1*sin(th1+pi/2) + N2*sin(th1+th2+pi/2)
Fx./Fy
fprintf("Should be %.4f\n\n", -4.6922/-.7556)
figure(1)
scatter(th3, N1, 100, '.k')
scatter(th3, N2, 100, '.k')

figure(2)
scatter(th3, N1*cos(th1+pi/2) + N2*cos(th1+th2+pi/2), 100, '.k')
scatter(th3, N1*sin(th1+pi/2) + N2*sin(th1+th2+pi/2), 100, '.k')

%% Finger Model (called from objective function):
function [tau, N1, N2, m] = updateFinger(th1, th2, th3)
    % Project-specific function calls (presumed available in your codebase)
    optFlag = false; % set optimization flag to true
    param = getParam(); % Not vector safe
    % param.K2 = 0.011;
    param = updateParam(param, th1, th2, th3); % Vector safe
    Jac   = getJacobians(param, th1, th2, th3); % Vector safe
    [cmp, F_cmp, TA, ~] = getSpring(param, th1, th2, th3); % Vector safe (cmp in mm, F_cmp in N, TA in Nm, k1 in N/mm)
    [N1v, N2v] = fastSolveN1N2(Jac.J1, Jac.J2, Jac.JC, ... % numerically solve the system of equations
                                  param.n1, param.n2, param.nC, ...
                                  F_cmp, param.K2, th2, param.l2_r, false);
    tau = TA(:); N1 = N1v(:); N2 = N2v(:); m = cmp(:); % output vectors
end