function dL = returnCompression(F_s, param, flag)

% load('OptimalSpring.mat')
param.F_s = F_s;
sol = solveBVP(param);
param = genParam(length(sol.x), param.d, flag);
% param = genParam(length(sol.x), param.d_max, flag); % only for variable thickness
% param = genParam(length(sol.x), param.d_max, param.d_min, flag); % only for variable thickness
y_comp = cumsum(sin(param.alpha+sol.y(1, :)))*mean(diff(sol.x));
dL = abs(param.L_n - y_comp(end));

end