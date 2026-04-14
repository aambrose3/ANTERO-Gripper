function [sol] = solveBVP(param)
    syms gam(s)
    r_n = param.r_n;
    E = param.E;
    F_s = param.F_s;
    s0 = param.s0;
    d = param.d;
    w_min = param.w_min;
    w_max = param.w_max;
    buff = param.buff;
    w = w_min + (w_max-w_min)/2.*(1+cos(2*pi*(s-buff)./(s0-2*buff) - pi)); % variable beam in-plane thickness
    % ecc = param.ecc;
    %% Only for variable thickness
    % d_min = param.d_min;
    % d_max = param.d_max;                                             % maximum thickness of the beam in the middle
    % d = d_min + (d_max-d_min)/2.*(1+cos(2*pi*s./s0 - pi));     % variable beam in-plane thickness
    % d = d_min + (d_max-d_min).*sin(pi*s./s0); % alternate profile
    A = d.*w;
    alpha = (s./r_n) - param.phi + (pi/2);
    a  = d./(exp(d./r_n)-1);
    b  = a + d;
    r_c = (a+b)/2;
    ecc = r_c - r_n;                                   % how far the neutral surface is from the center line?
    
    %%
    N_hat = -cos(alpha + gam);
    C = simplify(-diff(A)./A); % Variable Thickness
    % C = 0; % Constant Thickness
    D = simplify(-N_hat./A)./E./ecc./r_n;
    %% use for varaible thickness or variable width beams only
    [V] = odeToVectorField(diff(gam, 2) == C*diff(gam) + D*F_s);
    %% % use for constant contact area beams
    % [V] = odeToVectorField(diff(gam, 2) == 0*diff(gam) + D*F_s); 
    
    M = matlabFunction(V, 'vars', {'s', 'Y'});
    
    LMesh = linspace(0,param.L_n,length(param.t));
    yInit = bvpinit(param.t,@(s)initFcn(s,param.s0));

    sol = bvp4c(M,@bcfcn,yInit); 

end