function plotSpring(sol, param)
    x = cumsum(cos(param.alpha))*mean(diff(sol.x));
    y = cumsum(sin(param.alpha))*mean(diff(sol.x));
    x_comp = cumsum(cos(param.alpha-sol.y(1, :)))*mean(diff(sol.x));
    y_comp = cumsum(sin(param.alpha-sol.y(1, :)))*mean(diff(sol.x));
    hold on
    plot(x*1000, y*1000, '--b', 'LineWidth', 2)
    plot(x_comp*1000, y_comp*1000, 'b', 'LineWidth', 2)
    axis equal
    title('Spring Deflection at 120 N')
    ylabel('Position (mm)')
    xlabel('Position (mm)')
    legend('Natural Shape', 'Deflected Shape')
    hold off
end