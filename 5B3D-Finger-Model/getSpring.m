%% 5B3D-Finger-Model: Main
%% Author: Alexander B. Ambrose
%% Date: 8-6-2025

%% Description:
% The getSpring function calculates the compression spring force (link c). 
% This function outputs the the following
% cmp -> the amount the spring has compressed (mm)
% F_cmp -> the force produced by the spring compression (N)
% TA -> the torque required of the acutator (Nm)
% flag -> (-1, 0, or 1) for constraint handling (maximum spring compression 
% and actuator torque limits)

%% This function takes in the following inputs:
%   param -> the structure of parameters to update
%   t_1_ -> angle \theta_1 (scalar)
%   t_2_ -> angle \theta_2 (scalar)
%   t_3_ -> angle \theta_3 (scalar)

function [cmp, F_cmp, TA, flag] = getSpring(varargin)
    param = varargin{1};
    t_1_ = varargin{2};
    t_2_ = varargin{3};
    t_3_ = varargin{4};
    k1 = 0;
    optFlag = false;
    if nargin > 4
        k1 = varargin{5};
        optFlag = varargin{6};
    end

    if optFlag == false
        % check spring constraints
        if param.c < param.l1_o 
            flag = 1; % stop computation
        elseif param.c > param.l1_r
            flag = -1; % don't need to compute anythin
        else
            flag = 0; % perform computation
        end
    else
        flag = 0;
    end
    cmp = 1000*(param.l1_r-param.c); % in mm
    if optFlag == false
        % nonliear hookes law for spring
        F_cmp = -0.044880367*cmp^2 + 7.73061709*cmp; % Actual stiffness of d = 2.36mm spring
    else
        F_cmp = cmp.*k1; % k [=] N/mm
    end
    phi = param.tC - (t_3_+pi/2);
    TA = param.d.*F_cmp.*cos(phi);
    if optFlag == false && TA > param.TA_max
        flag = 1;
    end
end