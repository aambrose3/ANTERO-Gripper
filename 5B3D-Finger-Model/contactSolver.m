%% 5B3D-Finger-Model: Main
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

function [p1, p2, a, b] = contactSolver(param, t_1_, t_2_)

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

% Convert to double precision
p1 = vpa(p1);
p2 = vpa(p2);
a = vpa(a);
b = vpa(b);

end


