%% 5B3D-Finger-Model: contactSolver
%% Author: Alexander B. Ambrose
%% Date: 8-13-2025

%% Description:
% The updateParam function outputs various additinal parameters that are 
% dependent on finger configuration/state:
%   p1 -> contact point 1 along prot_1_imal link aligned with link axis (m)
%   p2 -> contact point 2 along distal link aligned with link axis (m)
%   a -> half of the width of the ellipse, e_a (m)
%   b -> half of the the height of the ellipse, e_b (m)

%% This function takes in the following inputs:
%   param -> the structure of parameters to update
%   t_1_ -> angle \theta_1 (scalar)
%   t_2_ -> angle \theta_2 (scalar)

function [p1, p2, a, b] = contactSolver(varargin) % param, t_1_, t_2_, optFlag
    param = varargin{1};
    t_1_ = varargin{2};
    t_2_ = varargin{3};
    optFlag = false;
    if nargin > 3
        optFlag = varargin{4};
    end
    
    syms p1 p2 a b real
    assume (p1 > 0);
    assume (p2 > 0);
    assume (a > 0);
    assume (b > 0);
    
    t1 = t_1_;
    t2 = t_2_;
    L = param.a;
    t = param.t;
    O = param.O;
    if numel(t1) == 1 && numel(t2) == 1 % if the provided t_1_ and t_2_ are scalars
        x1 = p1*cos(t1) - t*cos(t1-pi/2);
        y1 = p1*sin(t1) - t*sin(t1-pi/2);
        
        x2 = L*cos(t1) + p2*cos(t1+t2) - t*cos(t1+t2-pi/2);
        y2 = L*sin(t1) + p2*sin(t1+t2) - t*sin(t1+t2-pi/2);
        
        xc = O(1);
        yc = O(2) + b;
        
        % Continuity constraints (contact points are on the ellipse)
        EQ1 = (x1 - xc).^2./a^2 + (y1 - yc).^2./b^2 == 1;
        EQ2 = (x2 - xc).^2./a^2 + (y2 - yc).^2./b^2 == 1;
        
        % Tangency contraints
        EQ3 = (x1 - xc)./a^2*cos(t1) + (y1-yc)./b^2*sin(t1) == 0;
        EQ4 = (x2 - xc)./a^2*cos(t1+t2) + (y2-yc)./b^2*sin(t1+t2) == 0;
        
        % Solve system of nonlinear equations
        [p1, p2, a, b] = solve([EQ1, EQ2, EQ3, EQ4], [p1, p2, a, b]);
        %%%%%%%%%%%%%%%%%%%%%% This may result in terrible solutions, only input
        %%%%%%%%%%%%%%%%%%%%%% reasonable values for t_1 and t_2
    elseif numel(t1) > 1 && numel(t2) > 1 && numel(t1) == numel(t2) % if the provided t_1_ and t_2_ are vectors
        x1 = p1*cos(t1) - t*cos(t1-pi/2);
        y1 = p1*sin(t1) - t*sin(t1-pi/2);
        
        x2 = L*cos(t1) + p2*cos(t1+t2) - t*cos(t1+t2-pi/2);
        y2 = L*sin(t1) + p2*sin(t1+t2) - t*sin(t1+t2-pi/2);
        
        xc = O(1);
        yc = O(2) + b;

        % Continuity constraints (contact points are on the ellipse)
        EQ1 = (x1 - xc).^2./a^2 + (y1 - yc).^2./b^2 == 1;
        EQ2 = (x2 - xc).^2./a^2 + (y2 - yc).^2./b^2 == 1;
        
        % Tangency contraints
        EQ3 = (x1 - xc)./a^2.*cos(t1) + (y1-yc)./b^2.*sin(t1) == 0;
        EQ4 = (x2 - xc)./a^2.*cos(t1+t2) + (y2-yc)./b^2.*sin(t1+t2) == 0;
        
        for ii = 1:numel(t1)
            % Solve system of nonlinear equations
            [p1s(ii, 1), p2s(ii, 1), as(ii, 1), bs(ii, 1)] = ...
                solve([EQ1(ii), EQ2(ii), EQ3(ii), EQ4(ii)], [p1, p2, a, b]);
            %%%%%%%%%%%%%%%%%%%%%% This may result in terrible solutions, only input
            %%%%%%%%%%%%%%%%%%%%%% reasonable values for t_1 and t_2
        end
        p1 = p1s; p2 = p2s; a = as; b = bs;
    else % Bad
        fprintf('@contactSolver: input vectors (t1 and t2) do not have the same length.\n')
        while true == true %% Hang here
            pause(1);
        end
    end

    % Convert to double precision
    if optFlag == false
        p1 = vpa(p1);
        p2 = vpa(p2);
        a = vpa(a);
        b = vpa(b);
    end
end


