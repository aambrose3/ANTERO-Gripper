function [sx] = plotStress(X,Y,param)
    E  = param.E;
    rn = param.r_n;
    d = param.d;
    S  = X*0;
    for i = 1:length(X)
        a  = d/(exp(d/rn)-1);
        S(i) = abs(E*(1-rn/a)*rn*Y(2,i));
    end
    sx = (max(S)/1e6);
    
    figure
    hold on
    plot(X,S,'b', 'LineWidth',2);
    plot([X(1),X(end)],[param.St,param.St],'r--', 'LineWidth',2);
    title('Maximum Bending Stress along Beam Spring')
    ylabel('Maximum Bending Stress (Pa)')
    xlabel('Distance along beam spring (m)')
    hold off
end