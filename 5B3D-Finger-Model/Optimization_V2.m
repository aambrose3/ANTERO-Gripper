% main.m
clear all; close all; clc;
rng(42,'twister');                                  % reproducible
p = gcp('nocreate'); if isempty(p), parpool('processes'); end

%% Declare global variables
global theta_1 theta_2 theta_3

% Angle domains (radians)
theta1_range = [deg2rad(35), deg2rad(100)];
theta2_range = [deg2rad(35), deg2rad(100)];
theta3_range = [deg2rad(-70), deg2rad(40)]; % Possible actuator angles
% Example actuator angle grid (rad):
N_samples = 1001;  % same length as theta3 for vectorization convenience
theta_1 = linspace(theta1_range(1), theta1_range(2), N_samples);  % dense grid helps robustness
theta_2 = linspace(theta2_range(1), theta2_range(2), N_samples);  % dense grid helps robustness
theta_3 = linspace(theta3_range(1), theta3_range(2), N_samples);  % dense grid helps robustness

% Design variable bounds: x = [d; c; k1; k2]
lb = [0.045; 0.108; 4; 1e-3];        % [m; m; N/mm; Nm/rad]
ub = [0.060; 0.122; 10; 1e-1];       % [m; m; N/mm; Nm/rad]
x0 = [0.046; 0.118; 7.2; 0.02];      % initial guess

% ---- Options for fmincon ----
opts = optimoptions('fmincon', ...
    'Algorithm','interior-point', ...
    'Display','off', ...
    'MaxIterations', 1500, ...
    'MaxFunctionEvaluations', 2e5, ...
    'SpecifyObjectiveGradient', false, ...
    'SpecifyConstraintGradient', false, ...
    'FiniteDifferenceType','central', ...
    'StepTolerance',1e-10, ...
    'OptimalityTolerance',1e-8);

% ---- Problem definition ----
problem = createOptimProblem('fmincon', ...
    'objective', @(x) objectiveWrapper(x), ...
    'x0', x0, ...
    'lb', lb, 'ub', ub, ...
    'nonlcon', @(x) nonlconWrapper(x), ...
    'options', opts);

% ---- MultiStart with parallel ----
ms = MultiStart('UseParallel', true, 'Display','iter', 'FunctionTolerance', 1e-12);

% Latin Hypercube initial points (good coverage for nonconvex Ω)
nStarts = 200;                                       % bump to 500+ if time allows
Xlhs = lhsdesign(nStarts, numel(lb));
X0   = lb(:)'+ Xlhs.*(ub(:)'-lb(:)');
startpts = arrayfun(@(i) CustomStartPoint(X0(i,:).'), 1:nStarts, 'UniformOutput', false);

% Run
[ticVal, ~] = deal(tic);
[x_best, f_best, exitflag, out, sols] = run(ms, problem, [startpts{:}]); %#ok<ASGLU>
elapsed = toc(ticVal);

fprintf('\n=== Best solution ===\n');
disp(x_best), fprintf('Ω(x*) = %.6g\n', f_best);
fprintf('Elapsed: %.2f s, successful starts: %d/%d\n', elapsed, out.nSucc, nStarts);

save results_multistart.mat x_best f_best out sols lb ub




function val = objectiveWrapper(x)

    %% Archimedean Weights
    Q = [w_1, 0; 0, w_2];
    [tau, N1, N2, m] = updateFinger(x); 
    N = [N1;N2];
    Omega = tau./(N'*Q*N);
    
    % --- Example placeholder (REPLACE) ---
    Omega = toyOmega(x);  % remove this after wiring your function
    
    % If you want soft penalties inside the objective (optional):
    % [c, ceq] = nonlconWrapper(x);
    % pen = 0;
    % if any(c > 0) || any(abs(ceq) > 0)
    %     pen = 1e6 * (sum(max(c,0)) + sum(abs(ceq)));
    % end
    % val = Omega + pen;
    
    val = Omega;
end

function [tau, N1, N2, m] = updateFinger(th1, th2, th3, d, c, k1, k2)
    % Project-specific function calls (presumed available in your codebase)
    N = numel(th1);
    if ~(numel(th2)==N && numel(th3)==N), error('updateFinger: angle vectors must have equal length.'); end

    optFlag = true;
    param   = getParam(optFlag, d, c, k2);
    param   = updateParam(param, th1, th2, th3);
    Jac     = getJacobians(param, th1, th2, th3);
    [cmp, F_cmp, TA, ~] = getSpring(param, th1, th2, th3, k1, true);

    [N1v, N2v] = solveN1N2_linear(Jac.J1, Jac.J2, Jac.JC, ...
                                  param.n1, param.n2, param.nC, ...
                                  F_cmp, param.K2, th2, param.l2_r, true);

    tau = TA(:); N1 = N1v(:); N2 = N2v(:); m = cmp(:);
end

function [N1, N2] = solveN1N2_linear(J1, J2, JC, n1, n2, nC, F_cmp, K2, th2, l2_r, clampNonneg)
% Solve (J1' n1) N1 + (J2' n2) N2 + [0; K2*(th2 - l2_r)] + (JC' nC) F_cmp = 0.
    if nargin < 11 || isempty(clampNonneg), clampNonneg = true; end

    if ndims(J1) == 2
        J1 = reshape(J1,2,2,1); J2 = reshape(J2,2,2,1); JC = reshape(JC,2,2,1);
        if size(n1,2)==1, n1 = reshape(n1,2,1); end
        if size(n2,2)==1, n2 = reshape(n2,2,1); end
        if size(nC,2)==1, nC = reshape(nC,2,1); end
    end

    Npg = size(J1,3);  % pages
    a1 = zeros(2,Npg); a2 = zeros(2,Npg); aC = zeros(2,Npg);

    % broadcast scalars
    F_cmp = reshape(F_cmp, [], 1); if numel(F_cmp)==1, F_cmp = repmat(F_cmp,Npg,1); end
    K2    = reshape(K2,    [], 1); if numel(K2   )==1, K2    = repmat(K2,   Npg,1); end
    th2   = reshape(th2,   [], 1); if numel(th2  )==1, th2   = repmat(th2,  Npg,1); end
    l2_r  = double(l2_r);

    if size(n1,2)==1 && Npg>1, n1 = repmat(n1,1,Npg); end
    if size(n2,2)==1 && Npg>1, n2 = repmat(n2,1,Npg); end
    if size(nC,2)==1 && Npg>1, nC = repmat(nC,1,Npg); end

    % a1 = J1' n1 etc.
    for i=1:Npg
        a1(:,i) = (J1(:,:,i).') * n1(:,i);
        a2(:,i) = (J2(:,:,i).') * n2(:,i);
        aC(:,i) = (JC(:,:,i).') * nC(:,i);
    end

    Tau_h = [zeros(1,Npg); (K2(:).').*(th2(:).' - l2_r)];
    Tau_C = aC .* (F_cmp(:).');
    b     = Tau_h + Tau_C;

    a1x = a1(1,:); a1y = a1(2,:);
    a2x = a2(1,:); a2y = a2(2,:);
    bx  = b(1,:);  by  = b(2,:);

    detA = a1x.*a2y - a2x.*a1y;
    eps_det = 1e-12;
    near = abs(detA) < eps_det;
    detA(near) = sign(detA(near)).*eps_det + (~sign(detA(near)))*eps_det;

    N1v = -( a2y.*bx - a2x.*by ) ./ detA;
    N2v = -( -a1y.*bx + a1x.*by ) ./ detA;

    if clampNonneg
        N1v = max(N1v,0);
        N2v = max(N2v,0);
    end
    N1 = N1v(:); N2 = N2v(:);
end

function [c, ceq] = nonlconWrapper(x)
    global theta_3
    ceq = [];  % no equalities in your spec
    
    % ---- Compute all quantities your constraints need across theta ----
    % Replace the next block with calls to your model, producing arrays over theta:
    %   effort(theta), N1(theta), N2(theta), springComp(theta), ...
    %   plus any flags you need to infer contact and spans for h1, h2.
    [effort, N1, N2, springComp, contactOK] = mockKinematicsOverTheta(x, theta_3); %#ok<ASGLU>
    
    % ---- User thresholds (edit from your design spec) ----
    effortMax      = 1.0;      % g1: actuator effort ≤ limit
    springMaxComp  = 0.010;    % g6: max compression (m or rad equivalent)
    mustPushOnly   = true;     % g2: actuator cannot pull (effort ≥ 0)
    
    % ---- Point-wise g1..g6 -> inequalities (c ≤ 0) ----
    g1 = max(effort - effortMax);                 % ≤0
    g2 = mustPushOnly * max(-effort);             % ≤0 (if must push, no negatives)
    g3 = max(-N1);                                % N1 ≥ 0 -> -N1 ≤ 0
    g4 = max(-N2);                                % N2 ≥ 0
    g5 = max(-springComp);                        % spring compression ≥ 0
    g6 = max(springComp - springMaxComp);         % ≤ 0
    
    % ---- Domain-wide h1, h2 from your spec ----
    % h1: require an angle subset S (contiguous) with span > Δ1 where all pointwise are satisfied
    Delta1 = deg2rad(10);  % edit
    passMask = (effort <= effortMax) & (effort >= 0 | ~mustPushOnly) & ...
               (N1 >= 0) & (N2 >= 0) & (springComp >= 0) & (springComp <= springMaxComp);
    maxSpan = contiguousSpan(theta, passMask);
    h1 = -(maxSpan - Delta1);      % ≤0 when maxSpan ≥ Δ1
    
    % h2: pre-contact margin smaller than Δ2
    % interpret “distance from argmin(S) to θ at first contact” < Δ2
    Delta2 = deg2rad(5);   % edit
    thetaFirstContact = firstContactAngle(theta, springComp); % example proxy
    h2 = (argminS(theta, passMask) - thetaFirstContact) - Delta2; % ≤0 when distance ≤ Δ2
    
    c = [g1; g2; g3; g4; g5; g6; h1; h2];
end