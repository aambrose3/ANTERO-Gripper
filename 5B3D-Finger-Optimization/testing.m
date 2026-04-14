clear;close all;clc;
theta1_range = [deg2rad(50), deg2rad(100)]; % Entire Workspace (D = 50.8:152.4 mm)
theta2_range = [deg2rad(50), deg2rad(100)]; % Entire Workspace
theta3_range = [deg2rad(-70), deg2rad(40)]; % Every possible actuator angle
% Finger discretized angle grid (rad):
N_samples = 8;  % same length as theta3 for vectorization convenience
range3 = rad2deg(theta3_range(2))-rad2deg(theta3_range(1));
d3 = 10; % precision of actuator angle in degrees
theta_1 = linspace(theta1_range(1), theta1_range(2), N_samples+1);  % affects number of objects
theta_2 = linspace(theta2_range(1), theta2_range(2), N_samples+1);  % affects number of objects
theta_3 = linspace(theta3_range(1), theta3_range(2), d3*range3+1);  % dense grid helps robustness

[th1, th2] = meshgrid(theta_1, theta_2);

x = [0.04602; 0.11834; 7.2; 0.0353];
d = x(1);c_r = x(2);k1 = x(3);k2 = x(4);
optFlag = true; % set optimization flag to true
param = getParam(optFlag, x(1), x(2), x(4)); % Not vector safe


% [p1_og, p2_og, a_og, b_og] = fastContactSolver(param, theta_1(1), theta_2(1), [0.04, 0.04, 0.04, 0.04]);
% [p1, p2, a, b] = fastContactSolverTesting(param, theta_1(1), theta_2(1), [0.04, 0.04, 0.04, 0.04]);

th1
th2

[p1_og, p2_og, a_og, b_og] = fastContactSolver(param, th1(:), th2(:), [0.04, 0.04, 0.04, 0.04]);
[p1, p2, a, b] = fastContactSolverTesting(param, th1(:), th2(:), [0.04, 0.04, 0.04, 0.04]);


abs(p1-p1_og)./p1_og
abs(p2-p2_og)./p2_og
abs(a-a_og)./a_og
abs(b-b_og)./b_og