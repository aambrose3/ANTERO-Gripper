%% 5B3D-Finger-Model-Results: Main
%% Author: Alexander B. Ambrose
%% Date: 8-6-2025

%% Description:
% This script helps analyze the results from the 5B3D-Model: Main

clear all; close all; clc;
addpath('../')
% load('Results.mat');


%% Generate a video (.mp4) of the 
% rate = 10;
% name = 'Finger-Stable-Configurations.mp4';
% genVideo(param, N, t_1, t_2, t_3, rate, name);


%% Output Grasp Force Data for specific circular object diameter
% R = 101.6/1000/2; % input object diameter in m (4 in diameter)
% t_3 = linspace(deg2rad(-90), deg2rad(40), (90+40)*10+1)';
% [N_s, Nx, Ny, F] = circleCalc(param, R, t_3, 1);


%% Output Grasp Force Data for specific angles
% t_1 and t_2 are scalars
%% From Test Object in CAD
% t_3 = linspace(-0.8283, -0.6283, 201);
% t_1 = 1.0341;
% t_2 = 1.0363;
%% Another Test Object in CAD
% t_3 = linspace(0.2312, 0.2486, 2);
% t_1 = 1.5664;
% t_2 = 1.5446;
% [N_s, Nx, Ny, F] = ellipseCalc(param, t_1, t_2, t_3, 1);


%% Output the x and y-direction forces given specific enveloping force
% ref = [10 20]; % reference forces in N
% [Fx, Fy] = cartForces(F, Nx, Ny, ref)

load('OptimizationFunctions.mat');
% N1 = matlabFunction(S.N1);
% N2 = matlabFunction(S.N2);

% Minimize this objective function that simultaneously
% - minimizes actuator effort
% - maximizes the workspace
% - maximizes the contact force magnitude
% - averaged over the entire workspace of the gripper (integrate w.r.t t_1
% and t_2)

% Obj = - WS*sqrt(S.N1.^2+S.N2.^2)./TA

% subject to the following contraints:
% - limiting the contact forces to both be positive (N_1 > 0 && N_2 > 0)
% - limiting the actuator effort to a range (TA > 0 && TA < TA_max)
% - limiting spring compression (cmp > 0 && cmp < cmp_max)

% maximum physical ranges for the three angles
t_1 = linspace(deg2rad(30), deg2rad(100), 8); n1 = numel(t_1);
t_2 = linspace(deg2rad(20), deg2rad(100), 9); n2 = numel(t_2);
t_3 = linspace(deg2rad(-90), deg2rad(45), 46); n3 = numel(t_3);

% F = sqrt(S.N1^2+S.N2^2)/TA;
g_1 = S.N1;
g_2 = S.N2;
g_3 = cmp;
g_4 = TA;
clear S cmp TA
syms d l1_r t_3
N1 = num2cell(nan(n1, n2, n3)); N2 = N1; cmp = N1;
for ii = 1:n1
    for jj = 1:n2
        for kk = 1:n3
            N1{ii, jj, kk} = g_1(d, l1_r, t_1(ii), t_2(jj), t_3);
            N2{ii, jj, kk} = g_2(d, l1_r, t_1(ii), t_2(jj), t_3);
            cmp{ii, jj, kk} = g_3(d, l1_r, t_1(ii), t_2(jj), t_3);
        end
    end
end
