
param = genParam();

sol = solveBVP(param);      % sol.x = param.t, sol.y(1, :) = gamma, sol.y(2, :) = gamma'
L = trapz(sol.x, sin(param.alpha + sol.y(1, :)));
x = cumsum(cos(param.alpha))*mean(diff(sol.x));
y = cumsum(sin(param.alpha))*mean(diff(sol.x));
x_comp = cumsum(cos(param.alpha+sol.y(1, :)))*mean(diff(sol.x));
y_comp = cumsum(sin(param.alpha+sol.y(1, :)))*mean(diff(sol.x));
figure
hold on
plot(x, y, '--b', 'LineWidth', 2)
plot(x_comp, y_comp, 'b', 'LineWidth', 2)
axis equal

plotThickness(sol.x, param)
plotStress(sol.x, sol.y, param);