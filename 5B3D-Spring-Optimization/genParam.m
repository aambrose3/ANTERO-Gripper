function param = genParam(t, d, flag)
% function param = genParam(t, d, flag)
    %% Geometrical properties of the Beam Spring (SI)
    param.r_n = 0.06;                                                 % Natrual Radius of Curvature of the beam spring
    param.L_n = 118.342/1000;                                         % Pin-to-pin distance is 0.118342 m
    param.phi = asin(param.L_n/2/param.r_n);                  
    param.s0 = 2*param.phi*param.r_n;                                 % Free arc length of the beam spring
    param.t = linspace(0, param.s0, t);                               % Numerical s
    %% Constant Area beam spring
    param.d = d;                                                    % Constant in plane thickness of the beam spring
    param.w = 0.012446;                                               % Fixed width of the beam spring (out of plane)
    %% Only for variable thickness
    % param.d_min = 0.002;                                              % minimum thickness of the beams near the pins
    % param.d_min = d_min;
    % param.d_max = d;                                                  % maximum thickness of the beam in the middle
    % param.d = param.d_min + ...
        % (param.d_max-param.d_min)/2.*(1+...
        % cos(2*pi*param.t./param.s0 - pi));                            % variable beam in-plane thickness
    % Alternate profile
    % param.d = param.d_min + (param.d_max-param.d_min).*sin(pi*param.t./param.s0);
    %% Only for variable width
    param.w_min = 0.012446;
    param.w_max = 0.0254;
    param.buff = 0.01031875;
    param.w = param.w_min + ...
        (param.w_max-param.w_min)/2.*(1+...
        cos(2*pi*(param.t-param.buff)./(param.s0-2*param.buff) - pi));     % variable beam in-plane thickness
    %%
    param.A = param.d.*param.w;                                       % Cross-sectional Area of the beam spring
    param.alpha = (param.t/param.r_n) - param.phi + (pi/2);
    param.a  = param.d./(exp(param.d./param.r_n)-1);
    param.b  = param.a + param.d;
    param.r_c = (param.a+param.b)/2;
    param.ecc = param.r_c - param.r_n;                                   % how far the neutral surface is from the center line?
    if flag == 1 % flag for changing model
        param.L_n = param.L_n - 2*0.0719*25.4/1000; % shorten pin-to-pin distance
        param.ecc = param.ecc + 0.3041*25.4/1000; % need to increase eccentricity
    end
    
    %% Material Properties of the Beam Spring (SI)
    % % Aluminum 7075
    param.E  = 71.7e9;   % Elasticity,   7075 ~ 72e9,  M250 ~ 186e9,  M300 ~ 190e9,  M350 ~ 200e9,  Ti64 ~ 114e9
    param.St = 503e6;  % Yield,        7075 ~ 503e6, M250 ~ 1710e6, M300 ~ 1972e6, M350 ~ 2317e6, Ti64 ~ 880e6
    param.p  = 2768;   % Density,      7075 ~ 2810,  M250 ~ 8000,   M300 ~ 8000,   M350 ~ 8080,   Ti64 ~ 4430

    % % Ti-3AI-2.5V (Grade 9), alpha annealed
    % param.E = 100e9;
    % param.St = 500e6;
    % param.p = 4480;

    % % Ti-6AI-4V-0.5Ni-0.06Pd (grade 25)
    % param.E = 113.8e9;
    % param.St = 880e6;
    % param.p = 4430;

    % % G-10 
    % param.E = 18.6e9; % length wise
    % param.St = 241e6; % lower failure (cross wise)
    % param.p = 1910;

    % % Nylon 6/6 Extrusion
    % param.E = 2.45e9;
    % param.St = 75e6;
    % param.p = 1150;

    %% Loading Conditions (SI)
    % param.k_ideal = 6655;                       % Optimal Spring Rate in N/m
    % param.dL0 = 15.1352/1000;                      % Ideal compression amount in m
    % param.F_s = -124.5;
    param.k_ideal = 7200;                       % Optimal Spring Rate in N/m
    param.dL0 = 12.7/1000;                      % Ideal compression amount in m
    % param.k_ideal = 1808;                       % Actual measured stiffness 
    % param.dL0 = 12.7/1000;                      % Actual ideal deflection
    param.F_s = -param.k_ideal*param.dL0;        % Spring Compression load
    param.F_s = -102;
end
