clear all; close all; clc;

t = 1001; 
d = 0.001807; % inital guess of in plane thickness
param = genParam(t, d, 0);

fig = figure(1);
jj = 1; kk = 1; 
error = 1; error2 = 1;
epsilon = 1E-6; epsilon2 = 1E-5;
while error > epsilon % optimize maximum thickness -> Stiffness
    clf(fig)
    sol = solveBVP(param);      % sol.x = param.t, sol.y(1, :) = gamma, sol.y(2, :) = gamma'
    param = genParam(length(sol.x), d, 0);
    % param = genParam(length(sol.x), d_max, 0); % only for variable thickness
    L_c = cumsum(sin(param.alpha-sol.y(1, :)))*mean(diff(sol.x)); L_comp(kk) = L_c(end);
    dL(kk) = param.L_n - L_comp(kk);
    adj = 1 + 10*(dL(kk) - param.dL0);
    d = adj*param.d;
    % d_max = adj*param.d_max; % Only for variable thickness
    error = abs(param.dL0-dL(kk));
    kk = kk+1;
    plotSpring(sol, param)
    drawnow
end
thickness = param.d
try
    mass = trapz(param.t, param.A.*param.p)
catch
    mass = param.s0*param.A*param.p
end
param.d_0 = param.d;
% variable width spring version (7.2 N/mm)
param.d = 0.00236;
param.F_s = -100;
sol = solveBVP(param);      % sol.x = param.t, sol.y(1, :) = gamma, sol.y(2, :) = gamma'
plotThickness(sol.x, param);
S = plotStress(sol.x, sol.y, param);
FoS = (param.St*10^-6)/S

save('OptimalSpring.mat', 'param')

clear all;
load('OptimalSpring.mat')
f_s = linspace(-0.005, param.F_s, 100);

figure
for ii = 1:length(f_s)
    dL1(ii) = returnCompression(f_s(ii), param, 0);
    plotStiffness(dL1(ii)*1000, -f_s(ii)) % can only pass scalars
end
dL1 = dL1*1000;
p = polyfitZero(dL1, -f_s, 2)
plot(dL1, p(1)*dL1.^2 + p(2)*dL1, 'b', 'LineWidth', 2)
title('Simulated Beam Spring Force-Displacement', 'FontSize', 14)
ylabel('Compressive Force (N)', 'FontSize', 12)
xlabel('Displacement (mm)', 'FontSize', 12)
figure
k = -f_s./dL1;
plot(dL1, k, 'b', 'LineWidth', 2)
hold on
plot([0 param.dL0*1000], mean(k)*[1 1], 'k', 'LineWidth', 2)
title('Optimal Simulated Beam Spring Stiffness', 'FontSize', 14)
ylabel('Stiffness (N/mm)', 'FontSize', 12)
xlabel('Displacement (mm)', 'FontSize', 12)
param.f_s = -f_s;
param.k = k;
param.dL = dL1;

save('OptimalSpring_Plus.mat', 'param')

function [p,S,mu] = polyfitZero(x,y,degree)
% POLYFITZERO Fit polynomial to data, forcing y-intercept to zero.
%   P = POLYFITZERO(X,Y,N) is similar POLYFIT(X,Y,N) except that the
%   y-intercept is forced to zero, i.e. P(N) = 0. In the same way as
%   POLYFIT, the coefficients, P(1:N-1), fit the data Y best in the least-
%   squares sense. You can also use Y = POLYVAL(PZERO,X) to evaluate the
%   polynomial because the output is the same as POLYFIT.
%
%   [P,S,MU] = POLYFITZERO() Return structure, S, similar to POLYFIT for use
%   with POLYVAL to calculate error estimates of predictions with P.
%
%   [P,S,MU] = POLYFITZERO() Scale X by std(X), returns MU = [0,std(X)].
%
%   See also POLYVAL, POLYFIT
%
%   Also see <a href="http://www.mathworks.com/matlabcentral/fileexchange/34765-polyfitn">POLYFITN by John D'Errico</a>
%
% Copyright (c) 2013 Mark Mikofski
% Version 1-1, 2013-10-15
%   add delta output
%   center and scale
% Version 1-0, 2011-06-29
%% check args
% X & Y should be numbers
assert(isnumeric(x) && isnumeric(y),'polyfitZero:notNumeric', ...
    'X and Y must be numeric.')
dim = numel(x); % number of elements in X
% DEGREE should be scalar positive number between 1 & 10 inclusive
assert(isnumeric(degree) && isscalar(degree) && degree>0 && degree<=10, ...
    'polyfitZero:degreeOutOfRange', ...
    'DEGREE must be an integer between 1 and 10.')
% DEGREE must be less than number of elements in X & Y
assert(degree<dim && degree==round(degree), ...
    'polyfitZero:DegreeGreaterThanDim', 'DEGREE must be less than numel(X)')
% X & Y should be same size vectors
assert(isvector(x) && isvector(y) && dim==numel(y), ...
    'polyfitZero:vectorMismatch', 'X and Y must be vectors of the same length.')
%% solve
% convert X & Y to column vectors
x = x(:); y = y(:);
% Scale X.
% attribution: this is based on code from POLYFIT by The MathWorks Inc.
if nargout > 2
   mu = [0; std(x)];
   x = (x - mu(1))/mu(2);
end
% using pow() is actually as fast or faster than looping, same # of flops!
z = zeros(dim,degree);
for n = 1:degree
    z(:,n) = x.^(degree-n+1);
end
p = z\y; % solve
p = [p;0]; % set y-intercept to zero
%% error estimates
% attribution: this is based on code from POLYFIT by The MathWorks Inc.
if nargout > 1
    V = [z,ones(dim,1)]; % append constant term for Vandermonde matrix
    % Return upper-triangular factor of QR-decomposition for error estimates
    R = triu(qr(V,0));
    r = y - V*p;
    S.R = R(1:size(R,2),:);
    S.df = max(0,length(y) - (degree+1));
    S.normr = norm(r);
end
p = p'; % polynomial output is row vector by convention
end
