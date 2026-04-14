%% Optimization of the ANTERO Finger link + compliance parameters: d, c, k1, k2
%% Minimize  J(x) = F̄(x) - λ · ρ(x)
% ρ(x) measures 2‑D coverage over (θ1, θ2) bins and is 1 for a bin only if:
%   g7: there exists a consecutive θ3 span ≥ C7 inside [θ3_min, θ3_max], and
%   g8: “no bad predecessor”: for every feasible θ3 in the bin, there is no
%       compressed & infeasible θ3′ < θ3 with (θ3 - θ3′) ≤ C8.
%
% Angles are in radians throughout.
%
% Pointwise constraints g1..g6 (per sample):
%   g1:  τ - C1 < 0          (max actuator torque)
%   g2: -τ < 0               (actuator can only push)
%   g3:  εN - N1 < 0         (contact 1 positive)
%   g4:  εN - N2 < 0         (contact 2 positive)
%   g5: -m < 0               (spring compressed only)
%   g6:  m - C6 < 0          (max spring compression)

clear; close all; clc;
addpath("Results")

% -------------------------------------------------------------------------
% Global contact weights (Archimedean): contact = w_1.*N1 + w_2.*N2
% IMPORTANT: declare global BEFORE assigning.
% Defaults (if not set elsewhere): w_1 = 0.25; w_2 = 0.75
% -------------------------------------------------------------------------
global w_1 w_2
if isempty(w_1) || isempty(w_2)
    w_1 = 0.25;
    w_2 = 1 - w_1;
end

runEverything = false;  % set true to run optimization

if runEverything
    %% ------------------ USER SETTINGS ------------------
    % Angle domains (radians)
    theta1_range = [deg2rad(35), deg2rad(100)];
    theta2_range = [deg2rad(35), deg2rad(100)];
    theta3_range = [deg2rad(-70), deg2rad(40)];  % requested limits

    % Design variable bounds: x = [d; c; k1; k2]
    lb = [0.045; 0.108; 4; 1e-3];        % [m; m; N/mm; Nm/rad]
    ub = [0.060; 0.122; 10; 1e-1];       % [m; m; N/mm; Nm/rad]
    x0 = [0.046; 0.118; 7.2; 0.02];      % initial guess

    % Contact Importance (scalars)
    global w_1 w_2
    w_1 = 0.25;             % in [0,1]
    w_2 = 1 - w_1;          % -> 0.75
    ensure_contact_weights_on_workers(); % propagate to parallel workers if any

    % Pointwise constants
    C1    = 5.0;       % Nm
    C6    = 25;        % max compression (units consistent with getSpring)
    eps_N = 1e-3;      % min positive contact (N or model units)

    % Domain-wise θ3 rules (radians)
    C7 = deg2rad(5);    % required consecutive span for g7
    C8 = deg2rad(1);     % lookback window for g8

    % Coverage shaping on (θ1,θ2)
    nb1 = 13; nb2 = 13;
    w_uniform = 1.0;
    w_worst   = 0.8;
    tau_cov   = 0.05;     % softness for soft-min

    % Target per-bin (pass g7 & g8)
    rho_target     = 1.0;
    mu_cover_start = 0.0;
    mu_cover_end   = 1.0;

    % Objective / SAA
    lambda_base  = 0.10;
    lambda_power = 1.25;
    eps_contact  = 1.0;     % only for F̄ denominator guard (not a feasibility constraint)
    eps_den      = 1e-12;

    % Annealing
    anneal_rounds = 5;
    gamma_start   = 5e-1;
    gamma_end     = 1e-2;
    use_sobol     = true; rng(42);
    N_sched    = [2000, 4000, 6000, 8000, 10000];
    N_validate = 20000;

    % Initial sampling (resampled every round)
    N0 = N_sched(1);
    [Th1, Th2, Th3] = sample_thetas(N0, theta1_range, theta2_range, theta3_range, use_sobol);
    THETA.samples = [Th1, Th2, Th3];
    THETA.range1  = theta1_range;
    THETA.range2  = theta2_range;
    THETA.range3  = theta3_range;

    % θ3 bounds
    PROB.theta3_min = theta3_range(1);
    PROB.theta3_max = theta3_range(2);

    % Optimizer
    opts = optimoptions('fmincon', ...
        'Algorithm','sqp', 'Display','iter-detailed', 'ScaleProblem',true, ...
        'SpecifyObjectiveGradient',false, 'FiniteDifferenceType','forward', ...
        'FiniteDifferenceStepSize',1e-4, 'MaxFunctionEvaluations',5e4, ...
        'MaxIterations',300, 'UseParallel',true);

    % Ensure pool is up and weights reach workers before fmincon starts
    try
        pp = gcp('nocreate');
        if isempty(pp)
            try, parpool('threads'); catch, parpool('local'); end
        end
    catch
        % no parallel available; proceed serially
    end
    ensure_contact_weights_on_workers();

    %% ------------------ ANNEALED OPTIMIZATION ------------------
    x_curr = x0;
    for r = 1:anneal_rounds
        B = log(gamma_end./gamma_start)./(anneal_rounds-1);
        A = gamma_start./exp(B);
        gamma = A*exp(B*r);
        fprintf('\n=== Anneal round %d/%d, gamma = %.4g ===\n', r, anneal_rounds, gamma);

        Nr = N_sched(min(r, numel(N_sched)));
        [Th1, Th2, Th3] = sample_thetas(Nr, theta1_range, theta2_range, theta3_range, use_sobol);
        THETA.samples = [Th1, Th2, Th3];

        lambda_cov_r = lambda_base * (gamma_start / max(gamma,1e-12))^lambda_power;
        mu_cover_r   = mu_cover_start + (mu_cover_end - mu_cover_start) * (r-1)/max(anneal_rounds-1,1);

        PROB.C1 = C1; PROB.C6 = C6; PROB.C7 = C7; PROB.C8 = C8;
        PROB.lambda_cov = lambda_cov_r;
        PROB.gamma = gamma; PROB.eps_den = eps_den;
        PROB.THETA = THETA;
        PROB.nb1 = nb1; PROB.nb2 = nb2;
        PROB.w_uniform = w_uniform; PROB.w_worst = w_worst; PROB.tau_cov = tau_cov;
        PROB.rho_target = rho_target; PROB.mu_cover = mu_cover_r;
        PROB.eps_contact = eps_contact; PROB.eps_N = eps_N;
        PROB.lb = lb; PROB.ub = ub;

        obj = @(x) objective_softmask(x, PROB);
        x_curr = fmincon(obj, x_curr, [],[],[],[], lb, ub, [], opts);
    end

    x_opt = x_curr;
    fprintf('\nOptimized design: d = %.6f, c = %.6f, k1 = %.6f, k2 = %.6f\n', x_opt(1), x_opt(2), x_opt(3), x_opt(4));

    %% ------------------ VALIDATION ------------------
    [Th1v, Th2v, Th3v] = sample_thetas(N_validate, theta1_range, theta2_range, theta3_range, true);
    THETA_V.samples = [Th1v, Th2v, Th3v];
    THETA_V.range1  = theta1_range;
    THETA_V.range2  = theta2_range;
    THETA_V.range3  = theta3_range;

    PROB_V = PROB; PROB_V.THETA = THETA_V; PROB_V.gamma = gamma_end;

    [softJ, softFbar, softRho] = objective_softmask(x_opt, PROB_V);
    [hardJ, hardFbar, hardRho] = objective_hardmask(x_opt, PROB_V);

    fprintf('\nValidation on %d samples:\n', N_validate);
    fprintf(' Soft-mask:  J = %.8g,  Fbar = %.8g,  feasible fraction (2D) = %.4f\n', softJ, softFbar, softRho);
    fprintf(' Hard-mask:  J = %.8g,  Fbar = %.8g,  feasible fraction (2D) = %.4f\n', hardJ, hardFbar, hardRho);

    save('Results/Optimization_Out.mat', 'x_opt', 'PROB_V', 'nb1', 'nb2', 'lb', 'ub');

else
    %% ------------------ VISUALIZATION ONLY ------------------
    % Keep global weights consistent on workers too
    ensure_contact_weights_on_workers();

    load('Results/Optimization_Out.mat');
    if exist('lb','var')~=1 && exist('PROB_V','var') && isfield(PROB_V,'lb'), lb = PROB_V.lb; end
    if exist('ub','var')~=1 && exist('PROB_V','var') && isfield(PROB_V,'ub'), ub = PROB_V.ub; end
    if exist('lb','var')~=1, lb = [0.045; 0.108; 4; 1e-3]; end
    if exist('ub','var')~=1, ub = [0.060; 0.122; 10; 1e-1]; end

    visualizeCoverageMaps(x_opt, PROB_V, nb1, nb2, 'at x^*');
    scatterSoftVsHard(x_opt, PROB_V, nb1, nb2);
    tradeoffGammaAtXstar(x_opt, PROB_V, [0.3 0.2 0.1 0.05 0.02 0.01]);

    % d–c slice
    contourDesignSurface(PROB_V, x_opt, {'d','c'}, ...
        {linspace(lb(1), ub(1), 61), linspace(lb(2), ub(2), 61)}, ...
        'Approx','taylor1', 'UseBinSubsample',true, 'VizMaxPerBin',8, ...
        'Mask','soft','RhoKind','objective', 'UseParallel',true, 'FDStepRel',1e-3);

    % k1–k2 slice
    contourDesignSurface(PROB_V, x_opt, {'k1','k2'}, ...
        {linspace(lb(3), ub(3), 61), linspace(lb(4), ub(4), 61)}, ...
        'Mask','hard','RhoKind','binmean', ...
        'Approx','none', 'UseBinSubsample',true, 'VizMaxPerBin',6, ...
        'UseParallel',true, 'FDStepRel',1e-3);
end

%% =====================================================================
%%                        OPTIMIZATION FUNCTIONS
%% =====================================================================

function [Th1, Th2, Th3] = sample_thetas(N, r1, r2, r3, use_sobol)
    a1=r1(1); b1=r1(2); a2=r2(1); b2=r2(2); a3=r3(1); b3=r3(2);
    try
        if use_sobol
            p = sobolset(3,'Skip',1e3,'Leap',1e2); U = net(p,N);
        else
            U = lhsdesign(N,3,'criterion','maximin','iterations',50);
        end
    catch
        U = rand(N,3);
    end
    Th1 = a1 + (b1-a1)*U(:,1);
    Th2 = a2 + (b2-a2)*U(:,2);
    Th3 = a3 + (b3-a3)*U(:,3);
end

function [J, Fbar, rho2D] = objective_softmask(x, PROB)
    d=x(1); c=x(2); k1=x(3); k2=x(4);
    lambda=PROB.lambda_cov; gamma=PROB.gamma; epsden=PROB.eps_den;

    TH  = PROB.THETA.samples;
    th1 = TH(:,1); th2 = TH(:,2); th3 = TH(:,3);

    [tau,N1,N2,m] = updateFinger(th1, th2, th3, d, c, k1, k2);

    % ---- Contact metric (Archimedean weights) ----
    contact = contact_metric(N1, N2);
    F1 = tau ./ max(contact, PROB.eps_contact);

    g1 =  tau - PROB.C1;
    g2 = -tau;
    g3 =  PROB.eps_N - N1;
    g4 =  PROB.eps_N - N2;
    g5 = -m;
    g6 =  m - PROB.C6;
    G = [g1,g2,g3,g4,g5,g6];

    W = soft_mask(G, gamma);

    feas_base = all(G < 0, 2);
    within3   = (th3 >= PROB.theta3_min) & (th3 <= PROB.theta3_max);

    % binning on (θ1,θ2)
    if isfield(PROB.THETA,'range1'), r1=PROB.THETA.range1; else, r1=[min(th1),max(th1)]; end
    if isfield(PROB.THETA,'range2'), r2=PROB.THETA.range2; else, r2=[min(th2),max(th2)]; end
    e1 = linspace(r1(1), r1(2), PROB.nb1+1); e1(end)=e1(end)+1e-12;
    e2 = linspace(r2(1), r2(2), PROB.nb2+1); e2(end)=e2(end)+1e-12;
    b1 = discretize(th1, e1); b2 = discretize(th2, e2);

    % g7: θ3-span
    rho_bin_g7 = theta3_span_bins( ...
        b1, b2, th3, feas_base & within3 & ~isnan(b1) & ~isnan(b2), ...
        PROB.nb1, PROB.nb2, PROB.theta3_min, PROB.theta3_max, PROB.C7, pi/180);

    % g8: no bad predecessor
    compressed = (m > 0) & (m < PROB.C6);
    ok_g8 = theta3_local_recovery_bins( ...
        b1, b2, th3, feas_base & within3, compressed & within3, ...
        PROB.nb1, PROB.nb2, PROB.theta3_min, PROB.theta3_max, PROB.C8, pi/180);

    % gate samples by (g7 & g8) per-bin and θ3 bounds
    Gate_bin   = rho_bin_g7 .* ok_g8;          % (nb1*nb2)x1
    good12     = ~isnan(b1) & ~isnan(b2);
    Gate_samp  = zeros(size(th1));             % Nx1
    idx_lin    = sub2ind([PROB.nb1, PROB.nb2], b1(good12), b2(good12));
    Gate_samp(good12) = Gate_bin(idx_lin);
    Gate_samp  = Gate_samp .* within3;

    % effective weights for F̄
    W_eff = W .* Gate_samp;
    sumW_eff = sum(W_eff) + epsden;
    Fbar = sum(W_eff .* F1) / sumW_eff;

    % ρ over bins = binary pass/fail per bin
    rho_bin = rho_bin_g7 .* ok_g8;
    rho_binmean = mean(rho_bin);
    rho2D = rho_binmean;

    % coverage reward with soft-min
    worst_softmin = -PROB.tau_cov * log( mean( exp( -rho_bin / PROB.tau_cov ) ) );
    rho_reward = PROB.w_uniform * rho_binmean + PROB.w_worst * worst_softmin;

    % per-bin target penalty
    shortfall = max(0, PROB.rho_target - rho_bin);
    pen_cover = mean(shortfall.^2);

    J = (Fbar - lambda * rho_reward) + PROB.mu_cover * pen_cover;
end

function W = soft_mask(G, gamma)
    if gamma <= 0, gamma = 1e-12; end
    Z = -G ./ gamma;
    S = 1 ./ (1 + exp(-Z));
    W = exp( mean(log(max(S, realmin)), 2) );
end

function [J, Fbar, rho2D] = objective_hardmask(x, PROB)
    d=x(1); c=x(2); k1=x(3); k2=x(4);
    lambda=PROB.lambda_cov;

    TH  = PROB.THETA.samples; th1=TH(:,1); th2=TH(:,2); th3=TH(:,3);
    [tau,N1,N2,m] = updateFinger(th1, th2, th3, d, c, k1, k2);

    contact = contact_metric(N1, N2);
    F1 = tau ./ max(contact, PROB.eps_contact);

    g1 =  tau - PROB.C1;
    g2 = -tau;
    g3 =  PROB.eps_N - N1;
    g4 =  PROB.eps_N - N2;
    g5 = -m;
    g6 =  m - PROB.C6;
    G  = [g1,g2,g3,g4,g5,g6];
    feas = all(G < 0, 2);

    e1 = linspace(PROB.THETA.range1(1), PROB.THETA.range1(2), PROB.nb1+1); e1(end)=e1(end)+1e-12;
    e2 = linspace(PROB.THETA.range2(1), PROB.THETA.range2(2), PROB.nb2+1); e2(end)=e2(end)+1e-12;
    b1 = discretize(th1, e1); b2 = discretize(th2, e2);
    within3 = (th3 >= PROB.theta3_min) & (th3 <= PROB.theta3_max);

    rho_bin_g7 = theta3_span_bins(b1, b2, th3, feas & within3 & ~isnan(b1) & ~isnan(b2), ...
                                  PROB.nb1, PROB.nb2, PROB.theta3_min, PROB.theta3_max, PROB.C7, pi/180);
    compressed = (m > 0) & (m < PROB.C6);
    ok_g8 = theta3_local_recovery_bins(b1, b2, th3, feas & within3, compressed & within3, ...
                                       PROB.nb1, PROB.nb2, PROB.theta3_min, PROB.theta3_max, PROB.C8, pi/180);

    Gate_bin   = rho_bin_g7 .* ok_g8;
    good12     = ~isnan(b1) & ~isnan(b2);
    Gate_samp  = zeros(size(th1));
    idx_lin    = sub2ind([PROB.nb1, PROB.nb2], b1(good12), b2(good12));
    Gate_samp(good12) = Gate_bin(idx_lin);
    Gate_samp  = Gate_samp .* within3;

    pass_all = feas & (Gate_samp > 0);
    if any(pass_all)
        Fbar = mean(F1(pass_all));
    else
        Fbar = inf;
    end

    rho2D = mean(Gate_bin);
    J = Fbar - lambda * rho2D;
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

%% =====================================================================
%%                        VISUALIZATION
%% =====================================================================

function visualizeCoverageMaps(x, PROB, nb1, nb2, title_tag)
    TH  = PROB.THETA.samples; th1=TH(:,1); th2=TH(:,2); th3=TH(:,3);
    [tau,N1,N2,m] = updateFinger(th1, th2, th3, x(1), x(2), x(3), x(4));

    g1 =  tau - PROB.C1;
    g2 = -tau;
    g3 =  PROB.eps_N - N1;
    g4 =  PROB.eps_N - N2;
    g5 = -m;
    g6 =  m - PROB.C6;
    Gsoft = [g1,g2,g3,g4,g5,g6];

    Wsoft = soft_weights_from_G(Gsoft, PROB.gamma);

    feas_base = all(Gsoft<0,2) & (th3>=PROB.theta3_min) & (th3<=PROB.theta3_max);
    BIN0 = make_bin_index(PROB.THETA, PROB.nb1, PROB.nb2, PROB.THETA.range1, PROB.THETA.range2);

    rho_bin_g7 = theta3_span_bins(BIN0.b1, BIN0.b2, th3, feas_base, ...
        PROB.nb1, PROB.nb2, PROB.theta3_min, PROB.theta3_max, PROB.C7, pi/180);
    compressed = (m > 0) & (m < PROB.C6);
    ok_g8 = theta3_local_recovery_bins(BIN0.b1, BIN0.b2, th3, feas_base, compressed, ...
        PROB.nb1, PROB.nb2, PROB.theta3_min, PROB.theta3_max, PROB.C8, pi/180);

    prod_bin = rho_bin_g7 .* ok_g8;
    Whard = zeros(size(Wsoft));
    Whard(BIN0.good) = prod_bin(BIN0.idx);

    if isfield(PROB.THETA,'range1'), r1=PROB.THETA.range1; else, r1=[min(th1),max(th1)]; end
    if isfield(PROB.THETA,'range2'), r2=PROB.THETA.range2; else, r2=[min(th2),max(th2)]; end
    BIN = make_bin_index(PROB.THETA, nb1, nb2, r1, r2);

    rho_bin_soft = perBinFast(BIN, Wsoft);
    rho_bin_hard = perBinFast(BIN, Whard);

    Rsoft = reshape(rho_bin_soft, nb1, nb2);
    Rhard = reshape(rho_bin_hard, nb1, nb2);
    Rdiff = Rsoft - Rhard;

    e1 = BIN.edges{1}; e2 = BIN.edges{2};
    figure('Color','w','Position',[200 200 1200 360]); colormap("gray")
    subplot(1,3,1);
    imagesc(rad2deg(e2), rad2deg(e1), Rsoft); axis xy equal tight;
    colorbar;
    title('Soft Mask Workspace ($\gamma$ = 0.01)', 'Interpreter', 'latex', 'FontSize', 12);
    xlabel('$\theta_2$ (deg)', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$\theta_1$ (deg)', 'Interpreter', 'latex', 'FontSize', 12);

    subplot(1,3,2);
    imagesc(rad2deg(e2), rad2deg(e1), Rhard); axis xy equal tight;
    colorbar;
    title('Hard Mask Workspace ($\gamma \rightarrow$ 0)', 'Interpreter', 'latex', 'FontSize', 12);
    xlabel('$\theta_2$ (deg)', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$\theta_1$ (deg)', 'Interpreter', 'latex', 'FontSize', 12);

    subplot(1,3,3);
    imagesc(rad2deg(e2), rad2deg(e1), Rdiff); axis xy equal tight;
    colorbar;
    title(' Mask Differences', 'Interpreter', 'latex', 'FontSize', 12);
    xlabel('$\theta_2$ (deg)', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$\theta_1$ (deg)', 'Interpreter', 'latex', 'FontSize', 12);
end

function scatterSoftVsHard(x, PROB, nb1, nb2)
    TH  = PROB.THETA.samples; th1=TH(:,1); th2=TH(:,2); th3=TH(:,3);
    [tau,N1,N2,m] = updateFinger(th1, th2, th3, x(1), x(2), x(3), x(4));

    g1 =  tau - PROB.C1;
    g2 = -tau;
    g3 =  PROB.eps_N - N1;
    g4 =  PROB.eps_N - N2;
    g5 = -m;
    g6 =  m - PROB.C6;
    G = [g1,g2,g3,g4,g5,g6];

    Wsoft = soft_weights_from_G(G, PROB.gamma);

    feas_base = all(G<0,2) & (th3>=PROB.theta3_min) & (th3<=PROB.theta3_max);
    BIN0 = make_bin_index(PROB.THETA, PROB.nb1, PROB.nb2, PROB.THETA.range1, PROB.THETA.range2);

    rho_bin_g7 = theta3_span_bins(BIN0.b1, BIN0.b2, th3, feas_base, ...
        PROB.nb1, PROB.nb2, PROB.theta3_min, PROB.theta3_max, PROB.C7, pi/180);
    compressed = (m > 0) & (m < PROB.C6);
    ok_g8 = theta3_local_recovery_bins(BIN0.b1, BIN0.b2, th3, feas_base, compressed, ...
        PROB.nb1, PROB.nb2, PROB.theta3_min, PROB.theta3_max, PROB.C8, pi/180);

    prod_bin = rho_bin_g7 .* ok_g8;
    Whard = zeros(size(Wsoft)); Whard(BIN0.good) = prod_bin(BIN0.idx);

    if isfield(PROB.THETA,'range1'), r1=PROB.THETA.range1; else, r1=[min(th1),max(th1)]; end
    if isfield(PROB.THETA,'range2'), r2=PROB.THETA.range2; else, r2=[min(th2),max(th2)]; end
    BIN = make_bin_index(PROB.THETA, nb1, nb2, r1, r2);

    rs = perBinFast(BIN, Wsoft);
    rh = perBinFast(BIN, Whard);

    figure('Color','w'); hold on; grid on; box on;
    plot([0 1],[0 1],'k--','LineWidth',1);
    scatter(rh, rs, 28, 'filled', 'MarkerFaceAlpha',0.6);
    xlim([0 1]); ylim([0 1]);
    title('Coverage Agreement Between Soft and Hard Mask', 'Interpreter', 'latex', 'FontSize', 12);
    xlabel('Hard Mask Bin Coverage (\gamma \rightarrow 0)', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Soft Mask Bin Coverage (\gamma)', 'Interpreter', 'latex', 'FontSize', 12);
end

function tradeoffGammaAtXstar(xstar, PROB, gamma_list)
    TH  = PROB.THETA.samples; th1=TH(:,1); th2=TH(:,2); th3=TH(:,3);
    [tau,N1,N2,m] = updateFinger(th1, th2, th3, xstar(1), xstar(2), xstar(3), xstar(4));

    contact = contact_metric(N1, N2);
    F1 = tau ./ max(contact, PROB.eps_contact);

    if isfield(PROB.THETA,'range1'), r1=PROB.THETA.range1; else, r1=[min(th1),max(th1)]; end
    if isfield(PROB.THETA,'range2'), r2=PROB.THETA.range2; else, r2=[min(th2),max(th2)]; end
    BIN = make_bin_index(PROB.THETA, PROB.nb1, PROB.nb2, r1, r2);

    F = zeros(numel(gamma_list),1);
    R = zeros(numel(gamma_list),1);
    for i=1:numel(gamma_list)
        gamma = gamma_list(i);
        g1 =  tau - PROB.C1; g2 = -tau;
        g3 =  PROB.eps_N - N1; g4 =  PROB.eps_N - N2; g5 = -m; g6 = m - PROB.C6;
        G = [g1,g2,g3,g4,g5,g6];
        W = soft_weights_from_G(G, gamma);
        F(i) = sum(W .* F1) / (sum(W) + PROB.eps_den);
        R(i) = mean(perBinFast(BIN, W));   % visual metric
    end

    figure('Color','w','Position',[160 180 1200 360]);
    subplot(1,2,1);
    plot(gamma_list, 1./F,'-o'); grid on; set(gca,'XDir','reverse');
    title('Average Mechanical Advantage (1/F)', 'Interpreter', 'latex', 'FontSize', 12);
    xlabel('\gamma', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('1/F', 'Interpreter', 'latex', 'FontSize', 12);

    subplot(1,2,2);
    plot(gamma_list, R,'-o'); grid on; set(gca,'XDir','reverse');
    title('Grasp Workspace Coverage (\rho)', 'Interpreter', 'latex', 'FontSize', 12);
    xlabel('\gamma', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('\rho', 'Interpreter', 'latex', 'FontSize', 12);
end

function contourDesignSurface(PROB, x_star, varnames, grids, varargin)
    % 2‑D slice; grids is a 1x2 cell array {gridForVar1, gridForVar2}
    p = inputParser;
    addParameter(p,'Mask','soft',@(s)ischar(s)||isstring(s));
    addParameter(p,'RhoKind','objective',@(s)ischar(s)||isstring(s));
    addParameter(p,'Levels',20,@(n)isnumeric(n)&&isscalar(n));
    addParameter(p,'UseParallel',true,@(b)islogical(b)||ismember(b,[0 1]));
    addParameter(p,'Approx','taylor1',@(s) any(strcmpi(s,{'none','taylor1'})));
    addParameter(p,'FDStepRel',1e-3,@(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(p,'UseBinSubsample',true,@(b)islogical(b)||ismember(b,[0 1]));
    addParameter(p,'VizMaxPerBin',8,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
    parse(p,varargin{:});
    mask    = lower(string(p.Results.Mask));
    rhoKind = lower(string(p.Results.RhoKind));
    nlev    = p.Results.Levels;
    usePar  = logical(p.Results.UseParallel);
    approx  = lower(string(p.Results.Approx));
    fdrel   = p.Results.FDStepRel;
    doSub   = logical(p.Results.UseBinSubsample);
    kPerBin = p.Results.VizMaxPerBin;

    idx = @(name) switch_idx(name);
    i1 = idx(varnames{1});
    i2 = idx(varnames{2});

    g1v = grids{1}; g1v = g1v(:);
    g2v = grids{2}; g2v = g2v(:);
    n1 = numel(g1v); n2 = numel(g2v);
    JZ = zeros(n1, n2); FZ = JZ; RZ = JZ;

    TH = PROB.THETA.samples;
    th1_all=TH(:,1); th2_all=TH(:,2); th3_all=TH(:,3);
    if isfield(PROB.THETA,'range1'), r1 = PROB.THETA.range1; else, r1=[min(th1_all),max(th1_all)]; end
    if isfield(PROB.THETA,'range2'), r2 = PROB.THETA.range2; else, r2=[min(th2_all),max(th2_all)]; end
    BIN_all = make_bin_index(PROB.THETA, PROB.nb1, PROB.nb2, r1, r2);

    if doSub
        keep_idx = pick_reps_per_bin(BIN_all, kPerBin);
    else
        keep_idx = (1:size(TH,1)).';
    end

    th1=th1_all(keep_idx); th2=th2_all(keep_idx); th3=th3_all(keep_idx);
    THETA_sub.samples = [th1, th2, th3];
    BIN = make_bin_index(THETA_sub, PROB.nb1, PROB.nb2, r1, r2);
    useHard = (mask=="hard");

    if usePar
        try
            pp = gcp('nocreate');
            if isempty(pp)
                try, parpool('processes'); catch, parpool('local'); end
            end
            % propagate contact weights to workers here too (for viz-only use)
            global w_1 w_2
            parfevalOnAll(@set_contact_weights, 0, w_1, w_2);
        catch
            usePar = false;
        end
    end

    if approx == "taylor1"
        hv1 = max(fdrel*max(abs(x_star(i1)),1), 1e-6);
        hv2 = max(fdrel*max(abs(x_star(i2)),1), 1e-6);
        x0=x_star; x1p=x_star; x2p=x_star; x1p(i1)=x1p(i1)+hv1; x2p(i2)=x2p(i2)+hv2;

        [tau0,N10,N20,m0]     = updateFinger(th1, th2, th3, x0(1),  x0(2),  x0(3),  x0(4));
        [tau1p,N11p,N21p,m1p] = updateFinger(th1, th2, th3, x1p(1), x1p(2), x1p(3), x1p(4));
        [tau2p,N12p,N22p,m2p] = updateFinger(th1, th2, th3, x2p(1), x2p(2), x2p(3), x2p(4));

        dTau_dv1 = (tau1p - tau0)/hv1;   dTau_dv2 = (tau2p - tau0)/hv2;
        dN1_dv1  = (N11p - N10)/hv1;     dN1_dv2  = (N12p - N10)/hv2;
        dN2_dv1  = (N21p - N20)/hv1;     dN2_dv2  = (N22p - N20)/hv2;
        dM_dv1   = (m1p  - m0 )/hv1;     dM_dv2   = (m2p  - m0 )/hv2;

        if usePar
            parfor ii = 1:n1
                [JZ(ii,:),FZ(ii,:),RZ(ii,:)] = eval_row_taylor( ...
                    g1v(ii), g2v, x_star, i1, i2, ...
                    tau0,N10,N20,m0, dTau_dv1,dTau_dv2, dN1_dv1,dN1_dv2, dN2_dv1,dN2_dv2, dM_dv1,dM_dv2, ...
                    th1,th2,th3, PROB, BIN, useHard, rhoKind);
            end
        else
            for ii = 1:n1
                [JZ(ii,:),FZ(ii,:),RZ(ii,:)] = eval_row_taylor( ...
                    g1v(ii), g2v, x_star, i1, i2, ...
                    tau0,N10,N20,m0, dTau_dv1,dTau_dv2, dN1_dv1,dN1_dv2, dN2_dv1,dN2_dv2, dM_dv1,dM_dv2, ...
                    th1,th2,th3, PROB, BIN, useHard, rhoKind);
            end
        end

    else
        if usePar
            parfor ii = 1:n1
                Jrow=zeros(1,n2); Frow=Jrow; Rrow=Jrow;
                for jj = 1:n2
                    x = x_star; x(i1)=g1v(ii); x(i2)=g2v(jj);
                    [Jrow(jj),Frow(jj),Rrow(jj)] = local_eval_exact( ...
                        x, th1,th2,th3, PROB, BIN, useHard, rhoKind);
                end
                JZ(ii,:)=Jrow; FZ(ii,:)=Frow; RZ(ii,:)=Rrow;
            end
        else
            for ii = 1:n1
                Jrow=zeros(1,n2); Frow=Jrow; Rrow=Jrow;
                for jj = 1:n2
                    x = x_star; x(i1)=g1v(ii); x(i2)=g2v(jj);
                    [Jrow(jj),Frow(jj),Rrow(jj)] = local_eval_exact( ...
                        x, th1,th2,th3, PROB, BIN, useHard, rhoKind);
                end
                JZ(ii,:)=Jrow; FZ(ii,:)=Frow; RZ(ii,:)=Rrow;
            end
        end
    end

    [X,Y] = meshgrid(g2v, g1v);
    figure('Color','w','Position',[160 180 1200 360]);
    subplot(1,3,1);
    contourf(X,Y,JZ,nlev,'LineColor','none'); colorbar;
    title(sprintf('J(x)  [%s mask, \\rho=%s]', mask, rhoKind), 'Interpreter','latex', 'FontSize', 12);
    xlabel(label_of(i2), 'Interpreter','latex', 'FontSize', 12);
    ylabel(label_of(i1), 'Interpreter','latex', 'FontSize', 12); hold on;
    plot(x_star(i2), x_star(i1), 'kp', 'MarkerFaceColor','w', 'MarkerSize',8);

    subplot(1,3,2);
    contourf(X,Y,FZ,nlev,'LineColor','none'); colorbar;
    title('F̄(x)', 'Interpreter','latex', 'FontSize', 12);
    xlabel(label_of(i2), 'Interpreter','latex', 'FontSize', 12);
    ylabel(label_of(i1), 'Interpreter','latex', 'FontSize', 12); hold on;
    plot(x_star(i2), x_star(i1), 'kp', 'MarkerFaceColor','w', 'MarkerSize',8);

    subplot(1,3,3);
    contourf(X,Y,RZ,nlev,'LineColor','none'); colorbar;
    title(sprintf('\\rho(x)  [%s]', rhoKind), 'Interpreter','latex', 'FontSize', 12);
    xlabel(label_of(i2), 'Interpreter','latex', 'FontSize', 12);
    ylabel(label_of(i1), 'Interpreter','latex', 'FontSize', 12); hold on;
    plot(x_star(i2), x_star(i1), 'kp', 'MarkerFaceColor','w', 'MarkerSize',8);
end

function [Jrow,Frow,Rrow] = eval_row_taylor(v1_fixed, v2_grid, x_star, i1, i2, ...
    tau0,N10,N20,m0, dTau_dv1,dTau_dv2, dN1_dv1,dN1_dv2, dN2_dv1,dN2_dv2, dM_dv1,dM_dv2, ...
    th1,th2,th3, PROB, BIN, useHard, rhoKind)

    dv1 = (v1_fixed - x_star(i1));
    Jrow = zeros(1,numel(v2_grid)); Frow = Jrow; Rrow = Jrow;

    for jj = 1:numel(v2_grid)
        dv2 = (v2_grid(jj) - x_star(i2));

        tau = tau0 + dTau_dv1*dv1 + dTau_dv2*dv2;
        N1  = N10  + dN1_dv1 *dv1 + dN1_dv2 *dv2;
        N2  = N20  + dN2_dv1 *dv1 + dN2_dv2 *dv2;
        m   = m0   + dM_dv1  *dv1 + dM_dv2  *dv2;

        contact = contact_metric(N1, N2);
        F1 = tau ./ max(contact, PROB.eps_contact);

        g1 =  tau - PROB.C1; g2 = -tau;
        g3 =  PROB.eps_N - N1; g4 =  PROB.eps_N - N2; g5 = -m; g6 = m - PROB.C6;
        G  = [g1,g2,g3,g4,g5,g6];

        if ~useHard
            W = soft_weights_from_G(G, PROB.gamma);
            Fbar = sum(W.*F1)/(sum(W)+PROB.eps_den);
            rho_bin = perBinFast(BIN, W);
            if rhoKind=="binmean"
                rho = mean(rho_bin);
            else
                rho_uniform   = mean(rho_bin);
                worst_softmin = -PROB.tau_cov * log( mean( exp( -rho_bin / PROB.tau_cov ) ) );
                rho = PROB.w_uniform*rho_uniform + PROB.w_worst*worst_softmin;
            end
            shortfall = max(0, PROB.rho_target - rho_bin);
            J = (Fbar - PROB.lambda_cov*rho) + PROB.mu_cover*mean(shortfall.^2);
        else
            feas = all(G<0,2);
            if any(feas), Fbar = mean(F1(feas)); else, Fbar = inf; end

            within3 = (th3 >= PROB.theta3_min) & (th3 <= PROB.theta3_max);
            rho_bin_g7 = theta3_span_bins(BIN.b1, BIN.b2, th3, feas & within3, ...
                                          PROB.nb1, PROB.nb2, PROB.theta3_min, PROB.theta3_max, PROB.C7, pi/180);
            compressed = (m > 0) & (m < PROB.C6);
            ok_g8 = theta3_local_recovery_bins(BIN.b1, BIN.b2, th3, feas & within3, compressed & within3, ...
                                               PROB.nb1, PROB.nb2, PROB.theta3_min, PROB.theta3_max, PROB.C8, pi/180);
            if rhoKind=="binmean"
                rho = mean(rho_bin_g7 .* ok_g8);
            else
                rb = rho_bin_g7 .* ok_g8;
                rho_uniform   = mean(rb);
                worst_softmin = -PROB.tau_cov * log( mean( exp( -rb / PROB.tau_cov ) ) );
                rho = PROB.w_uniform*rho_uniform + PROB.w_worst*worst_softmin;
            end
            J = Fbar - PROB.lambda_cov*rho;
        end

        Jrow(jj)=J; Frow(jj)=Fbar; Rrow(jj)=rho;
    end
end

function [J,Fbar,rho] = local_eval_exact(x, th1,th2,th3, PROB, BIN, useHard, rhoKind)
    [tau,N1,N2,m] = updateFinger(th1, th2, th3, x(1), x(2), x(3), x(4));
    contact = contact_metric(N1, N2);
    F1 = tau ./ max(contact, PROB.eps_contact);
    g1 =  tau - PROB.C1; g2 = -tau;
    g3 =  PROB.eps_N - N1; g4 =  PROB.eps_N - N2; g5 = -m; g6 = m - PROB.C6;
    G  = [g1,g2,g3,g4,g5,g6];

    if ~useHard
        W = soft_weights_from_G(G, PROB.gamma);
        Fbar = sum(W.*F1)/(sum(W)+PROB.eps_den);
        rho_bin = perBinFast(BIN, W);
        if rhoKind=="binmean"
            rho = mean(rho_bin);
        else
            rho_uniform   = mean(rho_bin);
            worst_softmin = -PROB.tau_cov * log( mean( exp( -rho_bin / PROB.tau_cov ) ) );
            rho = PROB.w_uniform*rho_uniform + PROB.w_worst*worst_softmin;
        end
        shortfall = max(0, PROB.rho_target - rho_bin);
        J = (Fbar - PROB.lambda_cov * rho) + PROB.mu_cover * mean(shortfall.^2);
    else
        feas = all(G<0,2);
        if any(feas), Fbar = mean(F1(feas)); else, Fbar = inf; end
        within3 = (th3 >= PROB.theta3_min) & (th3 <= PROB.theta3_max);
        rho_bin_g7 = theta3_span_bins(BIN.b1, BIN.b2, th3, feas & within3, ...
                                      PROB.nb1, PROB.nb2, PROB.theta3_min, PROB.theta3_max, PROB.C7, pi/180);
        compressed = (m > 0) & (m < PROB.C6);
        ok_g8 = theta3_local_recovery_bins(BIN.b1, BIN.b2, th3, feas & within3, compressed & within3, ...
                                           PROB.nb1, PROB.nb2, PROB.theta3_min, PROB.theta3_max, PROB.C8, pi/180);
        if rhoKind=="binmean"
            rho = mean(rho_bin_g7 .* ok_g8);
        else
            rb = rho_bin_g7 .* ok_g8;
            rho_uniform   = mean(rb);
            worst_softmin = -PROB.tau_cov * log( mean( exp( -rb / PROB.tau_cov ) ) );
            rho = PROB.w_uniform*rho_uniform + PROB.w_worst*worst_softmin;
        end
        J = Fbar - PROB.lambda_cov * rho;
    end
end

%% =====================================================================
%%                        SMALL HELPERS
%% =====================================================================
function BIN = make_bin_index(THETA, nb1, nb2, r1, r2)
    th1 = THETA.samples(:,1); th2 = THETA.samples(:,2);
    if nargin < 4 || isempty(r1), r1 = [min(th1), max(th1)]; end
    if nargin < 5 || isempty(r2), r2 = [min(th2), max(th2)]; end
    e1 = linspace(r1(1), r1(2), nb1+1); e1(end)=e1(end)+1e-12;
    e2 = linspace(r2(1), r2(2), nb2+1); e2(end)=e2(end)+1e-12;
    b1 = discretize(th1, e1); b2 = discretize(th2, e2);
    good = ~isnan(b1) & ~isnan(b2);
    idx = sub2ind([nb1,nb2], b1(good), b2(good));
    cnt = accumarray(idx, 1, [nb1*nb2,1], @sum, 0);
    BIN.idx = idx; BIN.good = good; BIN.cnt = cnt;
    BIN.nb1 = nb1; BIN.nb2 = nb2; BIN.edges = {e1,e2};
    BIN.r1 = r1; BIN.r2 = r2; BIN.b1 = b1; BIN.b2 = b2;
end

function rho_bin = theta3_span_bins(b1, b2, th3, feas_vec, nb1, nb2, min_rad, max_rad, span_req_rad, rad_step)
    if nargin < 10 || isempty(rad_step), rad_step = pi/180; end
    valid = ~isnan(b1) & ~isnan(b2) & feas_vec;
    rho_bin = zeros(nb1*nb2,1);
    if ~any(valid), return; end

    bin_idx = sub2ind([nb1,nb2], b1(valid), b2(valid));
    edges = min_rad:rad_step:(max_rad + rad_step);
    L = max(1, ceil(span_req_rad / rad_step));

    bins_with_data = unique(bin_idx);
    for kk = 1:numel(bins_with_data)
        bb = bins_with_data(kk);
        sel = (bin_idx == bb);
        thb = th3(valid); thb = thb(sel);
        counts = histcounts(thb, edges);
        pres = counts > 0;                     % presence per θ3 bucket
        if any(movsum(double(pres), [L-1,0]) >= L)
            rho_bin(bb) = 1;
        end
    end
end

function ok_bin = theta3_local_recovery_bins(b1, b2, th3, feas_vec, compressed_vec, nb1, nb2, min_rad, max_rad, C8_rad, rad_step)
% g8: For each (θ1,θ2) bin, 1 iff for every feasible θ3 there is NO compressed & infeasible
% predecessor θ3' within (0, C8_rad] below it.
    if nargin < 11 || isempty(rad_step), rad_step = pi/180; end
    ok_bin = ones(nb1*nb2,1);

    valid = ~isnan(b1) & ~isnan(b2);
    if ~any(valid), return; end

    linbin = sub2ind([nb1,nb2], b1(valid), b2(valid));
    thb    = th3(valid);
    feasb  = feas_vec(valid);
    compb  = compressed_vec(valid);

    badb   = compb & ~feasb;

    edges  = min_rad:rad_step:(max_rad + rad_step);
    [~,~,binIdx] = histcounts(thb, edges);
    Lback = max(1, ceil(C8_rad / rad_step));

    bins = unique(linbin);
    for kk = 1:numel(bins)
        bb = bins(kk);
        sel = (linbin == bb);
        if ~any(sel), continue; end

        feas_idx = binIdx(sel & feasb);
        bad_idx  = binIdx(sel & badb);

        if isempty(feas_idx), continue; end  % vacuously satisfied; g7 will gate coverage

        bad_map = false(1, numel(edges)-1);
        bad_idx = unique(bad_idx(bad_idx>0));
        bad_map(bad_idx) = true;

        ufeas = unique(feas_idx(feas_idx>0)).';
        for k = ufeas
            j0 = max(1, k - Lback);
            if any(bad_map(j0:k-1))
                ok_bin(bb) = 0;
                break;
            end
        end
    end
end

function rho_bin = perBinFast(BIN, W)
    sum_w = accumarray(BIN.idx, W(BIN.good), [BIN.nb1*BIN.nb2, 1], @sum, 0);
    rho_bin = zeros(BIN.nb1*BIN.nb2,1);
    nz = BIN.cnt > 0; rho_bin(nz) = sum_w(nz) ./ BIN.cnt(nz);
end

function W = soft_weights_from_G(G, gamma)
    W = soft_mask(G, gamma);
end

function keep_idx = pick_reps_per_bin(BIN, kPerBin)
    good_idx = find(BIN.good);
    keep_idx = [];
    for b = 1:numel(BIN.cnt)
        if BIN.cnt(b) == 0, continue; end
        pos = find(BIN.idx == b);
        sel = pick_k_spaced(pos, min(kPerBin, numel(pos)));
        keep_idx = [keep_idx; good_idx(sel)]; %#ok<AGROW>
    end
    keep_idx = sort(keep_idx);
end

function sel = pick_k_spaced(vec, k)
    n = numel(vec);
    if k>=n, sel = vec; return; end
    loc = round(linspace(1,n,k)); sel = vec(loc(:));
end


function k = switch_idx(name)
    % Map variable name to index in x = [d;c;k1;k2].
    % Accepts: 'd','c','k1','k2' (char/string), or numeric 1..4.
    if isnumeric(name)
        name = round(name);
        if any(name == [1 2 3 4])
            k = name;
            return;
        else
            error('Unknown numeric index: %g (valid: 1..4)', name);
        end
    end
    name = lower(string(name));
    switch name
        case "d",   k = 1;
        case "c",   k = 2;
        case "k1",  k = 3;
        case "k2",  k = 4;
        otherwise,  error('Unknown variable name: %s', name);
    end
end

function lab = label_of(k)
    switch k
        case 1, lab = 'd (m)';
        case 2, lab = 'c (m)';
        case 3, lab = 'k_1 (N/mm)';
        case 4, lab = 'k_2 (Nm/rad)';
    end
end

function contact = contact_metric(N1, N2)
% Archimedean weighted contact: contact = w_1 .* N1 + w_2 .* N2
% Supports scalar or per-sample weights (column vectors).
    global w_1 w_2
    if isempty(w_1), w_1 = 0.25; end
    if isempty(w_2), w_2 = 0.75; end
    contact = w_1 .* N1 + w_2 .* N2;  % elementwise to avoid dim-mismatch
end

function ensure_contact_weights_on_workers()
% Ensure a pool exists and push global contact weights to all workers.
    global w_1 w_2
    try
        p = gcp('nocreate');
        if isempty(p)
            try, parpool('processes'); catch, parpool('local'); end
            p = gcp('nocreate');
        end
        if ~isempty(p)
            parfevalOnAll(@set_contact_weights, 0, w_1, w_2);
        end
    catch
        % ignore if Parallel Computing Toolbox is unavailable
    end
end

function set_contact_weights(a,b)
% Worker-side setter for global contact weights.
    global w_1 w_2
    w_1 = a; w_2 = b;
end
