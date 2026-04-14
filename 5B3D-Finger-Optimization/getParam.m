%% 5B3D-Finger-Model: Main
%% Author: Alexander B. Ambrose
%% Date: 8-6-2025

%% Description:
% The getParam function declares the fixed parameters of the 5-Bar-3-DOF Finger.

function param = getParam(varargin)
    optFlag = false;
    if nargin >=1
        optFlag = varargin{1};
    end
    if optFlag % for optimization only
        % syms d l1_r k1 k2
        d = varargin{2};
        l1_r = varargin{3};
        K2 = varargin{4};
        param.d = d;
        param.l1_r = l1_r;
        param.K2 = K2;
    else % Actual parameters of physical system
        param.d = 46.018/1000; % crank length (m)
        param.l1_r = .118342; % Resting length of the compression spring (m)
        % param.K2 = 0.0221; % torsion spring stiffness (Nm/rad)
        param.K2 = 0.02;
    end
    %% Fixed link Lengths of the 5-bar-3-dof finger design (m)
    param.a = 81.715/1000;
    param.b = 31.75/1000;
    % param.d set in if statement above
    param.e = 38.1/1000;
    param.t = 12.7/1000; % effective pad width
    
    %% Fixed Angles (rad)
    param.alpha = deg2rad(165); % angle of the link b with relation to theta_2. negated
    param.gamma = deg2rad(75); % global angle of the vector e. negated
    
    %% Fixed Coordinate of the center of the palm w.r.t. J1
    param.O = [-47.3271; 7.9769]/1000;
    
    %% Fixed Constraints
    param.TA_max = 5; % limited max torque in Nm the actuator can muster
    
    param.l1_o = param.l1_r-0.025; % limited minimal length of the compression spring (m)
    %% ideal torsional spring stiffness from Mechinery handbook -> 0.043 Nm/rad
    % wire diameter = 0.031" (0.7874)
    % internal coil diameter = 0.421" (10.69)
    % Material 304 SS cold drawn -> E = 28000-29000 ksi (193-200 GPa)
    % Active coils (with arms) -> N ~= 3.75
    % k = (E*d^4)/(4000*N*D) % all units in standard imperial
    % k = 0.0221-0.0281 Nm/rad
    param.l2_r = 0; % Resting angle of the restoration torsion spring (rad)
    param.l2_o = 11/18*pi; % maximum angle of the torsion spring (~110 deg)
end