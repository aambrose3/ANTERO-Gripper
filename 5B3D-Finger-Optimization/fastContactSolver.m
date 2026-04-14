function [p1, p2, a, b] = fastContactsolver(param, t1, t2, z0)
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
    maxIter = 25;
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
