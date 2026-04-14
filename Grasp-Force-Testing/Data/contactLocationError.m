clear; clc; close all;

param = getParam();
for ii = 1:5
    try
        clear finger
        clear fieldNames
        clear temp
    end
    if ii == 1
        load('Finger_Data_60x45.mat');
        GT{ii} = [42.282, 40.189, 68.389, 76.409]/1000;
    elseif ii == 2
        load('Finger_Data_60x60.mat');
        GT{ii} = [38.631, 43.083, 63.961, 62.790]/1000;
    elseif ii == 3
        load('Finger_Data_60x75.mat');
        GT{ii} = [35.398, 47.539, 60.490, 52.753]/1000;
    elseif ii == 4
        load('Finger_Data_75x75.mat');
        GT{ii} = [39.317, 43.657, 46.665, 44.124]/1000;
    else
        load('Finger_Data_90x90.mat');
        GT{ii} = [38.496, 47.327, 34.627, 30.519]/1000;
    end
    fieldNames = fieldnames(finger);
    for jj = 1:numel(fieldNames)
        Name = fieldNames{jj};
        idx = 251;
        if jj == 1
            t1{ii} = finger.(Name).t1(idx:end);
            t2{ii} = finger.(Name).t2(idx:end);
        else
            t1{ii} = [t1{ii}; finger.(Name).t1(idx:end)];
            t2{ii} = [t2{ii}; finger.(Name).t2(idx:end)];
        end
        % mean_t1(ii, jj) = mean(finger.(Name).p1(251:5000));
        % mean_t2(ii, jj) = mean(finger.(Name).p2(251:5000));
    end
    temp = table2array(combinations(unique(t1{ii}), unique(t2{ii})));
    t1{ii} = temp(:, 1); t2{ii} = temp(:, 2);
    [p1{ii}, p2{ii}, ea{ii} eb{ii}] = fastContactSolver(param, t1{ii}, t2{ii}, [0.04, 0.04, 0.04, 0.04]);
    p1_mean(ii) = mean(p1{ii});
    p1_std(ii) = std(p1{ii});
    p2_mean(ii) = mean(p2{ii});
    p2_std(ii) = std(p2{ii});
    ea_mean(ii) = mean(ea{ii});
    ea_std(ii) = std(ea{ii});
    eb_mean(ii) = mean(eb{ii});
    eb_std(ii) = std(eb{ii});
    ep1{ii} = (p1{ii} - GT{ii}(1))./GT{ii}(1);
    ep2{ii} = (p2{ii} - GT{ii}(2))./GT{ii}(2);
    eea{ii} = (ea{ii} - GT{ii}(3))./GT{ii}(3);
    eeb{ii} = (eb{ii} - GT{ii}(4))./GT{ii}(4);
end

% overall mean and deviation
% p1_mean = mean(p1)
% p1_std = max([std(p1), std(mean_p1)]) % max of the variance within and between trials
% p2_mean = mean(p2)
% p2_std = max([std(p2), std(mean_p2)]) % max of the variance within and between trials
% 
% 
p1_error = [ep1{1}; ep1{2}; ep1{3}; ep1{4}; ep1{5}];
p1_group = [ones(size(ep1{1})); 2*ones(size(ep1{2})); 3*ones(size(ep1{3})); ...
    4*ones(size(ep1{4})); 5*ones(size(ep1{5}))];
figure(1)
boxplot(p1_error, p1_group, 'Labels', {'60x45', '60x60', '60x75', '75x75', '90x90'})
xlabel('Finger Positon ($\theta_1$x $\theta_2$)', 'Interpreter', 'latex')
ylabel('Contact Location Error (mm)')
title('P1 Contact Location Error')

p2_error = [ep2{1}; ep2{2}; ep2{3}; ep2{4}; ep2{5}];
p2_group = [ones(size(ep2{1})); 2*ones(size(ep2{2})); 3*ones(size(ep2{3})); ...
    4*ones(size(ep2{4})); 5*ones(size(ep2{5}))];
figure(2)
boxplot(p2_error, p2_group, 'Labels', {'60x45', '60x60', '60x75', '75x75', '90x90'})
xlabel('Finger Positon ($\theta_1$x $\theta_2$)', 'Interpreter', 'latex')
ylabel('Contact Location Error (mm)')
title('P2 Contact Location Error')

ea_error = [eea{1}; eea{2}; eea{3}; eea{4}; eea{5}];
ea_group = [ones(size(eea{1})); 2*ones(size(eea{2})); 3*ones(size(eea{3})); ...
    4*ones(size(eea{4})); 5*ones(size(eea{5}))];
figure(3)
boxplot(ea_error, ea_group, 'Labels', {'60x45', '60x60', '60x75', '75x75', '90x90'})
xlabel('Finger Positon ($\theta_1$x $\theta_2$)', 'Interpreter', 'latex')
ylabel('Ellipse Axis Error (mm)')
title('Ellipse Axis Estimate: e_a')

eb_error = [eeb{1}; eeb{2}; eeb{3}; eeb{4}; eeb{5}];
eb_group = [ones(size(eeb{1})); 2*ones(size(eeb{2})); 3*ones(size(eeb{3})); ...
    4*ones(size(eeb{4})); 5*ones(size(eeb{5}))];
figure(4)
boxplot(eb_error, eb_group, 'Labels', {'60x45', '60x60', '60x75', '75x75', '90x90'})
xlabel('Finger Positon ($\theta_1$x $\theta_2$)', 'Interpreter', 'latex')
ylabel('Ellipse Axis Error (mm)')
title('Ellipse Axis Estimate: e_b')

%% Static Functions for kinematics and contact location estimations

function param = getParam(varargin)
    optFlag = false;
    if nargin >=1
        optFlag = varargin{1};
    end
    if optFlag % for optimization only
        % syms l1_r d
        d = varargin{2};
        l1_r = varargin{3};
        K2 = varargin{4};
        param.d = d;
        param.l1_r = l1_r;
        param.K2 = K2;
    else % Actual parameters of physical system
        param.d = 46.018/1000; % crank length (m)
        param.l1_r = .118342; % Resting length of the compression spring (m)
        param.K2 = 0.0221; % torsion spring stiffness (Nm/rad)
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
    param.l2_r = 0; %-deg2rad(30); % Resting angle of the restoration torsion spring (rad)
    param.l2_o = 11/18*pi; % maximum angle of the torsion spring (~110 deg)
end

function [p1, p2, a, b] = fastContactSolver(param, t1, t2, z0)
% Fast numeric solver for p1,p2,a,b using damped Newton.
% Inputs:
%   param : struct with fields a (L), t (thickness offset), O=[xc;yc]
%   t1,t2 : scalars OR column vectors of equal length
%   z0    : optional initial guess [p1;p2;a;b] (used for the first sample)
% Outputs are double vectors (same length as t1).
%
% Requires R2019a+ (uses basic ops only).

    if ~isvector(t1) || ~isvector(t2) || numel(t1)~=numel(t2)
        error('t1,t2 must be vectors of equal length.');
    end
    t1 = t1(:); t2 = t2(:);
    N  = numel(t1);

    p1 = zeros(N,1); p2 = zeros(N,1); a = zeros(N,1); b = zeros(N,1);

    % Initial guess
    if nargin < 4 || isempty(z0)
        % crude but safe: few cm contact, ellipse ~ finger dims
        z_prev = [0.02; 0.02; 0.03; 0.03]; % [p1;p2;a;b] (meters)
    else
        z_prev = z0(:);
    end

    for i = 1:N
        [zi, ok] = newton_contact_one(param, t1(i), t2(i), z_prev);
        if ~ok
            % fallback: try a reset guess once
            z_reset = max([0.01;0.01;0.02;0.02], 0.5*z_prev + 0.5*[0.02;0.02;0.03;0.03]);
            [zi, ok] = newton_contact_one(param, t1(i), t2(i), z_reset);
            if ~ok
                % last resort: accept current iterate even if not converged
                % (OR mark this sample infeasible upstream via soft mask)
            end
        end
        p1(i)=zi(1); p2(i)=zi(2); a(i)=zi(3); b(i)=zi(4);
        z_prev = zi; % warm start next sample
    end
end

function [z, converged] = newton_contact_one(P, t1, t2, z0)
% Damped Newton with log-parameterization (positivity enforced).

    % Unpack params (rename to match your symbols)
    L  = P.a;           % proximal length
    t  = P.t;           % pad thickness offset
    xc = P.O(1);
    yc0= P.O(2);        % yc = yc0 + b

    % unit directions
    c1 = cos(t1); s1 = sin(t1);
    c12= cos(t1+t2); s12= sin(t1+t2);

    % Log-param to enforce positivity
    y  = log(max(z0(:), 1e-6));  % y = log z
    maxIter = 250;
    tolR    = 1e-10;
    tolStep = 1e-10;
    converged = false;

    for k = 1:maxIter
        z  = exp(y);           % z = [p1;p2;a;b] > 0
        p1 = z(1); p2 = z(2); a = z(3); b = z(4);

        % Contact points (from your equations)
        x1 = p1*c1 - t*sin(t1);
        y1 = p1*s1 + t*c1;
        x2 = L*c1 + p2*c12 - t*sin(t1+t2);
        y2 = L*s1 + p2*s12 + t*c12;

        yc = yc0 + b;

        % Residuals: E1..E4 = 0
        Y1 = y1 - yc;  Y2 = y2 - yc;
        E1 = (x1 - xc)^2/a^2 + (Y1)^2/b^2 - 1;
        E2 = (x2 - xc)^2/a^2 + (Y2)^2/b^2 - 1;
        E3 = ((x1 - xc)/a^2)*c1  + (Y1/b^2)*s1;
        E4 = ((x2 - xc)/a^2)*c12 + (Y2/b^2)*s12;

        r  = [E1;E2;E3;E4];

        % Check residual norm
        if norm(r,2) < tolR
            converged = true;  z = [p1;p2;a;b];  return;
        end

        % Analytic Jacobian wrt z = [p1,p2,a,b]
        % Partials for contact 1
        dx1_dp1 = c1;  dy1_dp1 = s1;
        % Partials for contact 2
        dx2_dp2 = c12; dy2_dp2 = s12;

        % dE1/dp1
        dE1_dp1 = 2*(x1 - xc)/a^2 * dx1_dp1 + 2*(Y1)/b^2 * dy1_dp1;
        % dE1/dp2 = 0
        % dE1/da
        dE1_da  = -2*(x1 - xc)^2 / a^3;
        % dE1/db  (note Y1 = y1 - yc0 - b  => dY1/db = -1)
        dE1_db  = (-2*Y1)/b^2 - 2*(Y1^2)/b^3;

        % dE2/dp1 = 0
        % dE2/dp2
        dE2_dp2 = 2*(x2 - xc)/a^2 * dx2_dp2 + 2*(Y2)/b^2 * dy2_dp2;
        % dE2/da
        dE2_da  = -2*(x2 - xc)^2 / a^3;
        % dE2/db
        dE2_db  = (-2*Y2)/b^2 - 2*(Y2^2)/b^3;

        % dE3/dp1
        dE3_dp1 = (dx1_dp1/a^2)*c1 + (dy1_dp1/b^2)*s1;
        % dE3/dp2 = 0
        dE3_da  = -2*(x1 - xc)/a^3 * c1;
        % dE3/db
        dE3_db  = (-1/b^2)*s1 + ( -2*Y1 / b^3 )*s1;

        % dE4/dp1 = 0
        dE4_dp2 = (dx2_dp2/a^2)*c12 + (dy2_dp2/b^2)*s12;
        dE4_da  = -2*(x2 - xc)/a^3 * c12;
        dE4_db  = (-1/b^2)*s12 + ( -2*Y2 / b^3 )*s12;

        Jz = [ dE1_dp1,      0,   dE1_da, dE1_db;
                   0,  dE2_dp2,   dE2_da, dE2_db;
               dE3_dp1,      0,   dE3_da, dE3_db;
                   0,  dE4_dp2,   dE4_da, dE4_db ];

        % Chain rule: y = log z  =>  dz/dy = diag(z)
        % Jacobian wrt y: Jy = Jz * diag(z)
        Jy = Jz .* z.';  % column-wise scale

        % Solve normal eqs with small LM damping (trust-ish)
        % (Jy'Jy + mu I) * dy = -Jy' r
        g   = Jy.' * r;
        H   = Jy.' * Jy;
        mu  = 1e-9 * (trace(H)/4 + 1);  % scale-aware tiny damping
        dy  = -(H + mu*eye(4)) \ g;

        % Backtracking line search on ||r||
        step = 1.0;
        rnorm = norm(r,2);
        for bt = 1:8
            y_try = y + step*dy;
            zt = exp(y_try);
            % evaluate residual at trial (cheap)
            [rt] = residual_only(zt, t1, t2, L, t, xc, yc0);
            if norm(rt,2) < (1 - 1e-4*step)*rnorm
                y = y_try;  break;
            else
                step = step * 0.5;
            end
        end

        if norm(step*dy,2) < tolStep
            % small step; accept
            z = exp(y); converged = true; return;
        end
    end

    z = exp(y);  % return last iterate
end

function r = residual_only(z, t1, t2, L, t, xc, yc0)
    p1=z(1); p2=z(2); a=z(3); b=z(4);
    c1=cos(t1); s1=sin(t1); c12=cos(t1+t2); s12=sin(t1+t2);
    x1 = p1*c1 - t*sin(t1);  y1 = p1*s1 + t*c1;
    x2 = L*c1 + p2*c12 - t*sin(t1+t2);
    y2 = L*s1 + p2*s12 + t*c12;
    yc = yc0 + b;  Y1=y1-yc; Y2=y2-yc;
    E1 = (x1 - xc)^2/a^2 + (Y1)^2/b^2 - 1;
    E2 = (x2 - xc)^2/a^2 + (Y2)^2/b^2 - 1;
    E3 = ((x1 - xc)/a^2)*c1  + (Y1/b^2)*s1;
    E4 = ((x2 - xc)/a^2)*c12 + (Y2/b^2)*s12;
    r  = [E1;E2;E3;E4];
end