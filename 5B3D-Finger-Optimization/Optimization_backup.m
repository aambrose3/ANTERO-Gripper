
%% Optimization of the ANTERO Finger link + compliance parameters: d, c, k1, k2
%% SAA optimization for:  minimize  J(x) = F_bar(x) - lambda * rho(x)
% where F_bar is the weighted average of F1 over the FEASIBLE subset of (theta1,theta2,theta3)
% using a soft mask for differentiability, and rho is a *domain-wide* feasibility
% metric that rewards uniform coverage and lifts weak regions.
%
% Author: Alexander B. Ambrose 
% Date:   10/6/2025

clear; close all; clc;
addpath("Results")

runEverything = false; % flag to run the annealed optimization

if runEverything == true
    %% ------------------ USER SETTINGS ------------------
    % Angle domains (radians) to sample:
    theta1_range = [deg2rad(35), deg2rad(100)];
    theta2_range = [deg2rad(35), deg2rad(100)];
    theta3_range = [deg2rad(-70), deg2rad(40)];
    
    % Design variable bounds: x = [d; c; k1; k2]
    lb = [0.045; 0.0762; 4; 1e-3];        % [m; m; N/mm; Nm/rad]
    ub = [0.0762; 0.150; 10; 1e-1];         % [m; m; N/mm; Nm/rad]
    x0 = [0.046; 0.118; 7.2; 0.02];         % initial w/ guess same units as above
    
    % Constraint constants
    C1 = 5.0;    % for g1:  tau - C1 < 0 (max actuator effort [Nm])
    C6 = 25;     % for g5:  m - C6 < 0  (max spring compression [mm])  (units: keep consistent)
    
    % --- Feasibility pressure controls ---
    lambda_base   = 0.10;    % starting coverage weight
    lambda_power  = 1.25;     % how fast lambda grows as gamma drops (1–2 typical)
    
    % Per-bin coverage shaping (coverage "throughout")
    nb1 = 13; nb2 = 13;        % theta1/theta2 bins
    w_uniform = 1.0;         % weight on uniform (mean) bin coverage
    w_worst   = 0.8;         % weight on worst-bin soft-min (encourage lifting holes)
    tau_cov   = 0.05;        % softness for soft-min (smaller = harsher)
    
    % Coverage target per bin (push toward >= 20% per bin)
    rho_target = 0.50;       % target per-bin coverage (you asked for > 20%)
    mu_cover_start = 0.0;    % penalty weight at early rounds
    mu_cover_end   = 1.0;    % penalty weight by last round (tune 0.5–2)
    
    % Contact magnitude floor (keeps F1 denominator safe)
    eps_contact = 1.0;       % resultant contact threshold (unit consistent with h1,h2)
    
    % SAA / annealing settings
    anneal_rounds  = 5;             % number of annealing stages
    gamma_start    = 5e-1;          % soft-mask smoothness (start)
    gamma_end      = 1e-2;          % soft-mask smoothness (end)
    eps_den        = 1e-12;         % avoid divide-by-zero in averages
    use_sobol      = true;          % fall back to lhsdesign if false or Sobol unavailable
    rng(42);                         % reproducibility
    
    % (Optional) increase samples as the mask sharpens (improves bin stats)
    N_sched = [2000, 4000, 6000, 8000, 10000];  % per anneal round (will clamp to last)
    
    % Validation settings (bigger N to re-evaluate winners)
    N_validate = 20000;
    
    %% ------------------ INITIAL SAMPLING (will be resampled per round) ------------------
    N_samples = N_sched(1);
    [Th1, Th2, Th3] = sample_thetas(N_samples, theta1_range, theta2_range, theta3_range, use_sobol);
    
    % reuse same samples within a round to reduce noise (common random numbers)
    THETA.samples = [Th1, Th2, Th3];
    THETA.range1  = theta1_range;
    THETA.range2  = theta2_range;
    THETA.range3  = theta3_range;
    
    %% ------------------ OPTIMIZATION OPTIONS ------------------
    opts = optimoptions('fmincon',...
        'Algorithm','sqp', ...
        'Display','iter-detailed', ...
        'ScaleProblem', true, ...
        'SpecifyObjectiveGradient', false, ...
        'FiniteDifferenceType','forward', ...      % cheaper early rounds
        'FiniteDifferenceStepSize', 1e-4, ...
        'MaxFunctionEvaluations', 5e4, ...
        'MaxIterations', 300, ...                  % give early rounds more runway
        'UseParallel', true);                      % finite difference in parallel
    
    %% ------------------ ANNEALED OPTIMIZATION ------------------
    x_curr = x0;
    for r = 1:anneal_rounds
        % Anneal gamma via exponential decay
        B = log(gamma_end./gamma_start)./(anneal_rounds-1);
        A = gamma_start./exp(B);
        gamma = A*exp(B*r);
        fprintf('\n=== Anneal round %d/%d, gamma = %.4g ===\n', r, anneal_rounds, gamma);
    
        % (Optional) resample with more points as mask sharpens
        Nr = N_sched(min(r, numel(N_sched)));
        [Th1, Th2, Th3] = sample_thetas(Nr, theta1_range, theta2_range, theta3_range, use_sobol);
        THETA.samples = [Th1, Th2, Th3];
    
        % Coverage weight ramp (more pressure later)
        lambda_cov_r = lambda_base * (gamma_start / max(gamma,1e-12))^lambda_power;
    
        % Per-bin shortfall penalty ramp
        mu_cover_r = mu_cover_start + (mu_cover_end - mu_cover_start) * (r-1)/max(anneal_rounds-1,1);
    
        % Wrap data for the objective
        PROB.C1          = C1;
        PROB.C6          = C6;
        PROB.lambda_cov  = lambda_cov_r;
        PROB.gamma       = gamma;
        PROB.eps_den     = eps_den;
        PROB.THETA       = THETA;
    
        % pass feasibility shaping parameters
        PROB.nb1 = nb1; PROB.nb2 = nb2;
        PROB.w_uniform = w_uniform; PROB.w_worst = w_worst; PROB.tau_cov = tau_cov;
        PROB.rho_target = rho_target; PROB.mu_cover = mu_cover_r;
        PROB.eps_contact = eps_contact;
    
        % fmincon solve (bounds only; feasibility is handled by soft mask inside the objective)
        obj = @(x) objective_softmask(x, PROB);
        x_curr = fmincon(obj, x_curr, [],[],[],[], lb, ub, [], opts);
    end
    
    x_opt = x_curr;
    fprintf('\nOptimized design: d = %.6f, c = %.6f, k1 = %.6f, k2 = %.6f\n', x_opt(1), x_opt(2), x_opt(3), x_opt(4));
    
    %% ------------------ VALIDATION & REPORT ------------------
    % Resample with many points and evaluate both soft and hard objectives
    [Th1v, Th2v, Th3v] = sample_thetas(N_validate, theta1_range, theta2_range, theta3_range, true);
    THETA_V.samples = [Th1v, Th2v, Th3v];
    THETA_V.range1  = theta1_range;
    THETA_V.range2  = theta2_range;
    THETA_V.range3  = theta3_range;
    
    % Evaluate with small gamma for near-hard behavior
    PROB_V = PROB; PROB_V.THETA = THETA_V; PROB_V.gamma = gamma_end;
    
    [softJ, softFbar, softRho] = objective_softmask(x_opt, PROB_V);
    
    % Hard-mask (literal feasible-subset average) using same contact guard
    [hardJ, hardFbar, hardRho] = objective_hardmask(x_opt, PROB_V);
    
    fprintf('\nValidation on %d samples:\n', N_validate);
    fprintf(' Soft-mask:  J = %.8g,  Fbar = %.8g,  feasible fraction = %.4f\n', softJ, softFbar, softRho);
    fprintf(' Hard-mask:  J = %.8g,  Fbar = %.8g,  feasible fraction = %.4f\n', hardJ, hardFbar, hardRho);
    
    save('Results/Optimization_Out.mat', 'x_opt', 'PROB_V', 'nb1', 'nb2');
else
    %% Visualize the optimization metrics and optimal parameters
    load('Results/Optimization_Out.mat');
    
    visualizeCoverageMaps(x_opt, PROB_V, nb1, nb2, 'at x^*');
    scatterSoftVsHard(x_opt, PROB_V, nb1, nb2);
    tradeoffGammaAtXstar(x_opt, PROB_V, [0.3 0.2 0.1 0.05 0.02 0.01]);
    % These design plots can take a long time -> using parallel computing and downsampling to speed them up
    % d–c slice near optimal x
    contourDesignSurface(PROB_V, x_opt, {'d','c'}, ...
        {linspace(0.045,0.060,31), linspace(0.075,0.120,31)}, ...
        'Approx','taylor1', 'UseBinSubsample',true, 'VizMaxPerBin',8, ...
        'UseParallel',true);

    % k1-k2 slice near optimal x
    contourDesignSurface(PROB_V, x_opt, {'k1','k2'}, ...
        {linspace(4,20,33), linspace(1e-3,1,33)}, ...
        'Mask','hard','RhoKind','binmean', ...
        'Approx','none', 'UseBinSubsample',true, 'VizMaxPerBin',6, ...
        'UseParallel',true);
end


%% =====================================================================
%%                        OPTIMIZATION FUNCTIONS
%% =====================================================================

function [Th1, Th2, Th3] = sample_thetas(N, r1, r2, r3, use_sobol)
    a1 = r1(1); b1 = r1(2);
    a2 = r2(1); b2 = r2(2);
    a3 = r3(1); b3 = r3(2);

    try
        if use_sobol
            p = sobolset(3,'Skip',1e3,'Leap',1e2);
            U = net(p,N);
        else
            U = lhsdesign(N,3,'criterion','maximin','iterations',50);
        end
    catch
        % Fallback if sobolset/lhsdesign unavailable
        U = rand(N,3);
    end

    Th1 = a1 + (b1 - a1) * U(:,1);
    Th2 = a2 + (b2 - a2) * U(:,2);
    Th3 = a3 + (b3 - a3) * U(:,3);
end

function [J, Fbar, rho] = objective_softmask(x, PROB)
    % x = [d; c; k1; k2]
    d = x(1); c = x(2); k1 = x(3); k2 = x(4);

    C1     = PROB.C1;
    C6     = PROB.C6;
    lambda = PROB.lambda_cov;
    gamma  = PROB.gamma;
    epsden = PROB.eps_den;

    TH  = PROB.THETA.samples;
    th1 = TH(:,1); th2 = TH(:,2); th3 = TH(:,3);

    % ---- Model evaluations (vectorized) ----
    [tau, h1, h2, m] = updateFinger(th1, th2, th3, d, c, k1, k2);

    % ---- Objective pieces ----
    contact = sqrt(h1.^2 + h2.^2);
    denom   = max(contact, PROB.eps_contact);  % guard against small denom
    F1      = tau ./ denom;

    % ---- Vectorized inequality constraints g_i < 0 ----
    % g1:  tau - C1 < 0
    % g2: -tau      < 0
    % g3:  eps_contact - contact < 0  (require sufficient resultant contact)
    % g4: -m        < 0
    % g5:  m - C6   < 0
    G = [tau - C1, -tau, PROB.eps_contact - contact, -m, m - C6];  % [N x 5]

    % ---- Soft mask (product of logistics) ----
    W   = soft_mask(G, gamma);     % [N x 1], 0..1
    sumW = sum(W) + epsden;
    Fbar = sum(W .* F1) / sumW;

    % --------- UNIFORM COVERAGE OVER (theta1,theta2) BINS ---------
    % robust ranges
    if isfield(PROB.THETA,'range1'), r1=PROB.THETA.range1; else, r1=[min(th1),max(th1)]; end
    if isfield(PROB.THETA,'range2'), r2=PROB.THETA.range2; else, r2=[min(th2),max(th2)]; end

    nb1 = PROB.nb1; nb2 = PROB.nb2;
    e1  = linspace(r1(1), r1(2), nb1+1); e1(end) = e1(end) + 1e-12;
    e2  = linspace(r2(1), r2(2), nb2+1); e2(end) = e2(end) + 1e-12;

    b1 = discretize(th1, e1); b2 = discretize(th2, e2);
    good = ~isnan(b1) & ~isnan(b2);
    idx  = sub2ind([nb1,nb2], b1(good), b2(good));

    cnt   = accumarray(idx, 1, [nb1*nb2, 1], @sum, 0);
    sum_w = accumarray(idx, W(good), [nb1*nb2, 1], @sum, 0);
    rho_bin = zeros(nb1*nb2,1);
    nz = cnt > 0;
    rho_bin(nz) = sum_w(nz)./cnt(nz);  % mean W in each non-empty bin
    % empty bins remain zero (penalized by uniform/worst-bin metrics)

    % Uniform average across bins (each bin equal weight)
    rho_uniform   = mean(rho_bin);

    % Worst-bin soft-min (encourages lifting the weakest regions)
    worst_softmin = -PROB.tau_cov * log( mean( exp( -rho_bin / PROB.tau_cov ) ) );

    % Combine coverage rewards
    rho = PROB.w_uniform * rho_uniform + PROB.w_worst * worst_softmin;

    % Target coverage (per bin) penalty — pushes toward >= rho_target
    shortfall = max(0, PROB.rho_target - rho_bin);
    pen_cover = mean(shortfall.^2);            % smooth-ish (squared hinge)

    % Final objective
    J = (Fbar - lambda * rho) + PROB.mu_cover * pen_cover;
end

function W = soft_mask(G, gamma)
    % Elementwise sigmoid(-g/gamma), product across columns
    Z = -G ./ gamma;                % larger negative g => larger weight
    % use geometric mean of logistics to reduce underflow over many constraints
    S = 1 ./ (1 + exp(-Z));         % logistic
    % original: W = prod(S, 2);
    W = exp( mean(log(max(S, realmin)), 2) );  % geometric mean in (0,1]
end

function [J, Fbar, rho] = objective_hardmask(x, PROB)
    % Validation-only hard mask, using same resultant-contact guard
    d = x(1); c = x(2); k1 = x(3); k2 = x(4);

    C1     = PROB.C1;
    C6     = PROB.C6;
    lambda = PROB.lambda_cov;

    TH  = PROB.THETA.samples;
    th1 = TH(:,1); th2 = TH(:,2); th3 = TH(:,3);

    [tau, h1, h2, m] = updateFinger(th1, th2, th3, d, c, k1, k2);

    contact = sqrt(h1.^2 + h2.^2);
    F1 = tau ./ max(contact, PROB.eps_contact);

    G = [tau - C1, -tau, PROB.eps_contact - contact, -m, m - C6];  % < 0 feasible
    feas = all(G < 0, 2);

    if any(feas)
        Fbar = mean(F1(feas));
        rho  = mean(feas);
    else
        Fbar = inf;
        rho  = 0.0;
    end
    J = Fbar - lambda * rho;
end

function [tau, N1, N2, m] = updateFinger(th1, th2, th3, d, c, k1, k2)
    % FAST numeric solver for the finger kinematics and kinetics — no symbolic solver. 
    % Solves linear 2x2 for N1,N2.
    
    % Returns:
    %   tau = actuator torque vector (TA)
    %   N1  = contact 1 magnitude
    %   N2  = contact 2 magnitude
    %   m   = cmp (spring compression)
    %
    % Assumes th1, th2, th3 are paired column vectors of equal length N.

    N = numel(th1);
    if ~(numel(th2)==N && numel(th3)==N)
        error('updateFinger: angle vectors must have equal length.');
    end

    % ---- Fixed / geometry-dependent parameters
    optFlag = true; % as in your codebase
    param   = getParam(optFlag, d, c, k2);

    % ---- Update kinematics for all samples
    param   = updateParam(param, th1, th2, th3);   % keep vectorized
    Jac     = getJacobians(param, th1, th2, th3);  % expects 2x2xN-ish

    % ---- Spring deflection & actuator torque
    [cmp, F_cmp, TA, ~] = getSpring(param, th1, th2, th3, k1, true);

    % Solve linear contact system for each sample in vectorized form
    [N1v, N2v] = solveN1N2_linear(Jac.J1, Jac.J2, Jac.JC, ...
                                  param.n1, param.n2, param.nC, ...
                                  F_cmp, param.K2, th2, param.l2_r, true);

    % Outputs
    tau = TA(:);
    N1  = N1v(:);
    N2  = N2v(:);
    m   = cmp(:);
end

function [N1, N2] = solveN1N2_linear(J1, J2, JC, n1, n2, nC, F_cmp, K2, th2, l2_r, clampNonneg)
% Solve (J1' n1) N1 + (J2' n2) N2 + [0; K2*(th2 - l2_r)] + (JC' nC) F_cmp = 0
% for N1, N2. Works for scalar or batched inputs (pages).

    if nargin < 11 || isempty(clampNonneg), clampNonneg = true; end

    if ndims(J1) == 2
        % single sample -> treat as one page
        J1 = reshape(J1,2,2,1); J2 = reshape(J2,2,2,1); JC = reshape(JC,2,2,1);
        if size(n1,2)==1, n1 = reshape(n1,2,1); end
        if size(n2,2)==1, n2 = reshape(n2,2,1); end
        if size(nC,2)==1, nC = reshape(nC,2,1); end
    end

    N = size(J1,3);

    % Expand vectors to row form for broadcasting
    F_cmp = reshape(F_cmp, [], 1); if numel(F_cmp)==1, F_cmp = repmat(F_cmp,N,1); end
    K2    = reshape(K2,    [], 1); if numel(K2   )==1, K2    = repmat(K2,   N,1); end
    th2   = reshape(th2,   [], 1); if numel(th2  )==1, th2   = repmat(th2,  N,1); end
    l2_r  = double(l2_r);

    % Make n1,n2,nC as 2xN
    if size(n1,2)==1 && N>1, n1 = repmat(n1,1,N); end
    if size(n2,2)==1 && N>1, n2 = repmat(n2,1,N); end
    if size(nC,2)==1 && N>1, nC = repmat(nC,1,N); end

    % Compute a1 = J1' * n1, a2 = J2' * n2, aC = JC' * nC  (all 2xN)
    J1T = permute(J1,[2 1 3]); J2T = permute(J2,[2 1 3]); JCT = permute(JC,[2 1 3]);
    havePM = exist('pagemtimes','file')==2;

    if havePM
        a1 = squeeze(pagemtimes(J1T, reshape(n1,2,1,[]))); % 2xN
        a2 = squeeze(pagemtimes(J2T, reshape(n2,2,1,[]))); % 2xN
        aC = squeeze(pagemtimes(JCT, reshape(nC,2,1,[]))); % 2xN
    else
        a1 = zeros(2,N); a2 = zeros(2,N); aC = zeros(2,N);
        for i=1:N
            a1(:,i) = J1T(:,:,i) * n1(:,i);
            a2(:,i) = J2T(:,:,i) * n2(:,i);
            aC(:,i) = JCT(:,:,i) * nC(:,i);
        end
    end

    % b = Tau_h + Tau_C
    Tau_h = [zeros(1,N); (K2(:).').*(th2(:).' - l2_r)];     % 2xN
    Tau_C = aC .* (F_cmp(:).');                             % 2xN
    b     = Tau_h + Tau_C;                                  % 2xN

    % Solve A*[N1;N2] = -b, where A = [a1x a2x; a1y a2y]
    a1x = a1(1,:); a1y = a1(2,:);
    a2x = a2(1,:); a2y = a2(2,:);
    bx  = b(1,:);  by  = b(2,:);

    detA = a1x.*a2y - a2x.*a1y;

    % Regularize near-singulars
    eps_det = 1e-12;
    near    = abs(detA) < eps_det;
    detA(near) = sign(detA(near)).*eps_det + (~sign(detA(near)))*eps_det;

    % Closed-form inverse times (-b):
    N1v = -( a2y.*bx - a2x.*by ) ./ detA;
    N2v = -( -a1y.*bx + a1x.*by ) ./ detA;

    if clampNonneg
        N1v = max(N1v, 0);
        N2v = max(N2v, 0);
    end

    N1 = N1v(:);
    N2 = N2v(:);
end

%% =====================================================================
%%                        VISUALIZATION FUNCTIONS
%% =====================================================================

function visualizeCoverageMaps(x, PROB, nb1, nb2, title_tag)
    % Fast: one model eval + precomputed bin indices reused for soft & hard.
    
    TH  = PROB.THETA.samples;
    th1 = TH(:,1); th2 = TH(:,2); th3 = TH(:,3);
    
    % ---- expensive model ONCE ----
    [tau,h1,h2,m] = updateFinger(th1, th2, th3, x(1), x(2), x(3), x(4));
    contact = sqrt(h1.^2 + h2.^2);
    
    % ---- weights ONCE ----
    Gsoft = [tau - PROB.C1, -tau, PROB.eps_contact - contact, -m, m - PROB.C6];
    Wsoft = soft_weights_from_G(Gsoft, PROB.gamma);
    Whard = double(all(Gsoft < 0, 2));
    
    % ---- make bin index ONCE ----
    if isfield(PROB.THETA,'range1'), r1=PROB.THETA.range1; else, r1=[min(th1),max(th1)]; end
    if isfield(PROB.THETA,'range2'), r2=PROB.THETA.range2; else, r2=[min(th2),max(th2)]; end
    BIN = make_bin_index(PROB.THETA, nb1, nb2, r1, r2);
    
    % ---- fast per-bin means (no discretize inside loop) ----
    rho_bin_soft = perBinFast(BIN, Wsoft);
    rho_bin_hard = perBinFast(BIN, Whard);
    
    Rsoft = reshape(rho_bin_soft, nb1, nb2);
    Rhard = reshape(rho_bin_hard, nb1, nb2);
    Rdiff = Rsoft - Rhard;
    
    rho_soft_mean = mean(rho_bin_soft);
    rho_hard_mean = mean(rho_bin_hard);
    
    % ---- plot ----
    e1 = BIN.edges{1}; e2 = BIN.edges{2};
    figure('Color','w','Position',[200 200 1200 360]);
    colormap("gray")

    subplot(1,3,1);
    imagesc(rad2deg(e2), rad2deg(e1), Rsoft); axis xy equal tight;
    cb=colorbar;
    title('Soft Mask Workspace ($\gamma$ = 0.01)', 'Interpreter', 'latex', 'FontSize', 12);
    xlabel('$\theta_2$ (deg)', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$\theta_1$ (deg)', 'Interpreter', 'latex', 'FontSize', 12); 

    
    subplot(1,3,2);
    imagesc(rad2deg(e2), rad2deg(e1), Rhard); axis xy equal tight;
    cb=colorbar;
    title('Hard Mask Workspace ($\gamma \rightarrow$ 0)', 'Interpreter', 'latex', 'FontSize', 12);
    xlabel('$\theta_2$ (deg)', 'Interpreter', 'latex', 'FontSize', 12); 
    ylabel('$\theta_1$ (deg)', 'Interpreter', 'latex', 'FontSize', 12); 

    
    subplot(1,3,3);
    imagesc(rad2deg(e2), rad2deg(e1), Rdiff); axis xy equal tight;
    cb=colorbar;
    title(' Mask Differences', 'Interpreter', 'latex', 'FontSize', 12);
    xlabel('$\theta_2$ (deg)', 'Interpreter', 'latex', 'FontSize', 12); 
    ylabel('$\theta_1$ (deg)', 'Interpreter', 'latex', 'FontSize', 12); 

    
    % if nargin>=4 && ~isempty(title_tag)
    %     sgtitle('Annealed Workspace Maps');
    % end
end

function scatterSoftVsHard(x, PROB, nb1, nb2)
    % Fast: one model eval + shared bin indices.
    
    TH  = PROB.THETA.samples; th1=TH(:,1); th2=TH(:,2); th3=TH(:,3);
    
    % single model eval
    [tau,h1,h2,m] = updateFinger(th1, th2, th3, x(1), x(2), x(3), x(4));
    contact = sqrt(h1.^2 + h2.^2);
    G = [tau - PROB.C1, -tau, PROB.eps_contact - contact, -m, m - PROB.C6];
    Wsoft = soft_weights_from_G(G, PROB.gamma);
    Whard = double(all(G<0,2));
    
    % bin index once
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
    xlabel('Hard Mask Bin Coverage (\gamma = 0.01)', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Soft Mask Bin Coverage (\gamma = 0.2)', 'Interpreter', 'latex', 'FontSize', 12); 

end

function tradeoffGammaAtXstar(xstar, PROB, gamma_list)
    % Fast: evaluate model ONCE at x*, then vary gamma cheaply; reuse bins.
    
    TH  = PROB.THETA.samples; th1=TH(:,1); th2=TH(:,2); th3=TH(:,3);
    
    % one model eval at x*
    [tau,h1,h2,m] = updateFinger(th1, th2, th3, xstar(1), xstar(2), xstar(3), xstar(4));
    contact = sqrt(h1.^2 + h2.^2);
    F1 = tau ./ max(contact, PROB.eps_contact);
    
    % bin index once
    if isfield(PROB.THETA,'range1'), r1=PROB.THETA.range1; else, r1=[min(th1),max(th1)]; end
    if isfield(PROB.THETA,'range2'), r2=PROB.THETA.range2; else, r2=[min(th2),max(th2)]; end
    BIN = make_bin_index(PROB.THETA, PROB.nb1, PROB.nb2, r1, r2);
    
    F = zeros(numel(gamma_list),1);
    R = zeros(numel(gamma_list),1);
    for i=1:numel(gamma_list)
        gamma = gamma_list(i);
        G = [tau - PROB.C1, -tau, PROB.eps_contact - contact, -m, m - PROB.C6];
        W = soft_weights_from_G(G, gamma);
        F(i) = sum(W .* F1) / (sum(W) + PROB.eps_den);
        R(i) = mean(perBinFast(BIN, W));           % bin-mean coverage (visual metric)
    end
    
    figure('Color','w');

    subplot(1,2,1); 
    plot(gamma_list, 1./F,'-o'); grid on; set(gca,'XDir','reverse');
    title('Optimum Average Mechanical Advantage (1/F)', 'Interpreter', 'latex', 'FontSize', 12);
    xlabel('Average Mechanical Advantage (N/Nm)', 'Interpreter', 'latex', 'FontSize', 12); 
    ylabel('\gamma', 'Interpreter', 'latex', 'FontSize', 12); 

    subplot(1,2,2); 
    plot(gamma_list, R,'-o'); grid on; set(gca,'XDir','reverse');
    title('Optimum Grasp Workspoace Coverage (\rho)', 'Interpreter', 'latex', 'FontSize', 12);
    xlabel('Grasp Workspace Coverage (%)', 'Interpreter', 'latex', 'FontSize', 12); 
    ylabel('\gamma', 'Interpreter', 'latex', 'FontSize', 12); 


end

function contourDesignSurface(PROB, x_star, varnames, grids, varargin)
    % Fast 2-D slice: optional θ-subsampling + first-order Taylor surrogate.
    % contourDesignSurface(PROB, x_star, {'d','c'}, {grid1, grid2}, ...
    %   'Mask','soft','RhoKind','objective','Levels',20, ...
    %   'UseParallel',true, 'Approx','taylor1', 'FDStepRel',1e-3, ...
    %   'UseBinSubsample',true, 'VizMaxPerBin',8)
    
    % ---- options ----
    p = inputParser;
    addParameter(p,'Mask','soft',@(s)ischar(s)||isstring(s));
    addParameter(p,'RhoKind','objective',@(s)ischar(s)||isstring(s));
    addParameter(p,'Levels',20,@(n)isnumeric(n)&&isscalar(n));
    addParameter(p,'UseParallel',true,@(b)islogical(b)||ismember(b,[0 1]));
    addParameter(p,'Approx','taylor1',@(s) any(strcmpi(s,{'none','taylor1'})));
    addParameter(p,'FDStepRel',1e-3,@(x)isnumeric(x) && isscalar(x) && x>0);
    addParameter(p,'UseBinSubsample',true,@(b)islogical(b)||ismember(b,[0 1]));
    addParameter(p,'VizMaxPerBin',8,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
    parse(p,varargin{:});
    mask      = lower(string(p.Results.Mask));
    rhoKind   = lower(string(p.Results.RhoKind));
    nlev      = p.Results.Levels;
    usePar    = logical(p.Results.UseParallel);
    approx    = lower(string(p.Results.Approx));
    fdrel     = p.Results.FDStepRel;
    doSub     = logical(p.Results.UseBinSubsample);
    kPerBin   = p.Results.VizMaxPerBin;
    
    % ---- var mapping ----
    idx = @(name) switch_idx(lower(string(name)));
    i1 = idx(varnames{1});
    i2 = idx(varnames{2});
    
    g1 = grids{1}(:); g2 = grids{2}(:);
    n1 = numel(g1); n2 = numel(g2);
    [JZ,FZ,RZ] = deal(zeros(n1, n2));
    
    % ---- precompute THETA, binning, and optional subsample ----
    TH  = PROB.THETA.samples;
    th1_all = TH(:,1); th2_all = TH(:,2); th3_all = TH(:,3);
    if isfield(PROB.THETA,'range1'), r1 = PROB.THETA.range1; else, r1=[min(th1_all),max(th1_all)]; end
    if isfield(PROB.THETA,'range2'), r2 = PROB.THETA.range2; else, r2=[min(th2_all),max(th2_all)]; end
    BIN_all = make_bin_index(PROB.THETA, PROB.nb1, PROB.nb2, r1, r2);
    
    if doSub
        keep_idx = pick_reps_per_bin(BIN_all, kPerBin);        % <= k per bin
    else
        keep_idx = (1:size(TH,1)).';
    end
    
    th1 = th1_all(keep_idx); th2 = th2_all(keep_idx); th3 = th3_all(keep_idx);
    THETA_sub.samples = [th1, th2, th3];
    BIN = make_bin_index(THETA_sub, PROB.nb1, PROB.nb2, r1, r2);  % binning for subset
    
    useHard = (mask=="hard");
    
    % ---- optional parallel pool ----
    if usePar
        try
            pp = gcp('nocreate');
            if isempty(pp)
                try, parpool('threads'); catch, parpool('local'); end
            end
        catch
            usePar = false;
        end
    end
    
    % ---- choose evaluation strategy ----
    if approx == "taylor1"
        % First-order Taylor surrogate around x_star: 3 model calls total.
        hv1 = max(fdrel*max(abs(x_star(i1)),1), 1e-6);
        hv2 = max(fdrel*max(abs(x_star(i2)),1), 1e-6);
    
        x0  = x_star;
        x1p = x_star; x1p(i1) = x1p(i1) + hv1;
        x2p = x_star; x2p(i2) = x2p(i2) + hv2;
    
        [tau0,h10,h20,m0]    = updateFinger(th1, th2, th3, x0(1),  x0(2),  x0(3),  x0(4));
        [tau1p,h11p,h21p,m1p]= updateFinger(th1, th2, th3, x1p(1), x1p(2), x1p(3), x1p(4));
        [tau2p,h12p,h22p,m2p]= updateFinger(th1, th2, th3, x2p(1), x2p(2), x2p(3), x2p(4));
    
        dTau_dv1 = (tau1p - tau0)/hv1;   dTau_dv2 = (tau2p - tau0)/hv2;
        dH1_dv1  = (h11p - h10)/hv1;     dH1_dv2  = (h12p - h10)/hv2;
        dH2_dv1  = (h21p - h20)/hv1;     dH2_dv2  = (h22p - h20)/hv2;
        dM_dv1   = (m1p  - m0 )/hv1;     dM_dv2   = (m2p  - m0 )/hv2;
    
        % row-sliced loops (parfor-friendly)
        if usePar
            parfor ii = 1:n1
                [JZ(ii,:),FZ(ii,:),RZ(ii,:)] = eval_row_taylor( ...
                    g1(ii), g2, x_star, i1, i2, ...
                    tau0,h10,h20,m0, dTau_dv1,dTau_dv2, dH1_dv1,dH1_dv2, dH2_dv1,dH2_dv2, dM_dv1,dM_dv2, ...
                    th1,th2,th3, PROB, BIN, useHard, rhoKind);
            end
        else
            for ii = 1:n1
                [JZ(ii,:),FZ(ii,:),RZ(ii,:)] = eval_row_taylor( ...
                    g1(ii), g2, x_star, i1, i2, ...
                    tau0,h10,h20,m0, dTau_dv1,dTau_dv2, dH1_dv1,dH1_dv2, dH2_dv1,dH2_dv2, dM_dv1,dM_dv2, ...
                    th1,th2,th3, PROB, BIN, useHard, rhoKind);
            end
        end
    
    else % approx == "none"  (exact, but slower)
        if usePar
            parfor ii = 1:n1
                Jrow = zeros(1,n2); Frow = zeros(1,n2); Rrow = zeros(1,n2);
                for jj = 1:n2
                    x = x_star; x(i1)=g1(ii); x(i2)=g2(jj);
                    [Jrow(jj),Frow(jj),Rrow(jj)] = local_eval_exact( ...
                        x, th1,th2,th3, PROB, BIN, useHard, rhoKind);
                end
                JZ(ii,:) = Jrow; FZ(ii,:) = Frow; RZ(ii,:) = Rrow;
            end
        else
            for ii = 1:n1
                Jrow = zeros(1,n2); Frow = zeros(1,n2); Rrow = zeros(1,n2);
                for jj = 1:n2
                    x = x_star; x(i1)=g1(ii); x(i2)=g2(jj);
                    [Jrow(jj),Frow(jj),Rrow(jj)] = local_eval_exact( ...
                        x, th1,th2,th3, PROB, BIN, useHard, rhoKind);
                end
                JZ(ii,:) = Jrow; FZ(ii,:) = Frow; RZ(ii,:) = Rrow;
            end
        end
    end
    
    % ---- plot ----
    [X,Y] = meshgrid(g2, g1);
    figure('Color','w','Position',[160 180 1200 360]);
    
    subplot(1,3,1); contourf(X,Y,JZ,nlev,'LineColor','none'); colorbar;
    title(sprintf('J(x)  [%s mask, \\rho=%s]', mask, rhoKind));
    xlabel(label_of(i2)); ylabel(label_of(i1)); hold on;
    plot(x_star(i2), x_star(i1), 'kp', 'MarkerFaceColor','w', 'MarkerSize',8);
    
    subplot(1,3,2); contourf(X,Y,FZ,nlev,'LineColor','none'); colorbar;
    title('F̄(x)'); xlabel(label_of(i2)); ylabel(label_of(i1)); hold on;
    plot(x_star(i2), x_star(i1), 'kp', 'MarkerFaceColor','w', 'MarkerSize',8);
    
    subplot(1,3,3); contourf(X,Y,RZ,nlev,'LineColor','none'); colorbar;
    title(sprintf('\\rho(x)  [%s]', rhoKind));
    xlabel(label_of(i2)); ylabel(label_of(i1)); hold on;
    plot(x_star(i2), x_star(i1), 'kp', 'MarkerFaceColor','w', 'MarkerSize',8);
    
    sgtitle(sprintf('Design slice: %s vs %s', upper(string(varnames{1})), upper(string(varnames{2}))));

end

%% =====================================================================
%%                        SMALL HELPER FUNCTIONS
%% =====================================================================

function k = switch_idx(name)
    % Map variable name to index in x = [d;c;k1;k2]
    switch name
        case 'd',  k = 1;
        case 'c',  k = 2;
        case 'k1', k = 3;
        case 'k2', k = 4;
        otherwise, error('Unknown variable name: %s', name);
    end
end

function lab = label_of(k)
    % Pretty axis labels with units
    switch k
        case 1, lab = 'd (m)';
        case 2, lab = 'c (m)';
        case 3, lab = 'k_1 (N/mm)';
        case 4, lab = 'k_2 (Nm/rad)';
    end
end

function [J,Fbar,rho] = local_eval_at(v1, v2, x_star, i1, i2, th1, th2, th3, PROB, BIN, useHard, rhoKind)
    % Evaluate J, Fbar, rho at a single grid point (parfor safe).
    
    x = x_star; x(i1) = v1; x(i2) = v2;
    
    [tau,h1,h2,m] = updateFinger(th1, th2, th3, x(1), x(2), x(3), x(4));
    contact = sqrt(h1.^2 + h2.^2);
    F1 = tau ./ max(contact, PROB.eps_contact);
    
    G = [tau - PROB.C1, -tau, PROB.eps_contact - contact, -m, m - PROB.C6];
    
    if ~useHard
        % soft J, Fbar
        W = soft_weights_from_G(G, PROB.gamma);
        sumW = sum(W) + PROB.eps_den;
        Fbar = sum(W .* F1) / sumW;
    
        % rho (objective or binmean)
        rho_bin = perBinFast(BIN, W);
        if rhoKind=="binmean"
            rho = mean(rho_bin);
        else
            rho_uniform   = mean(rho_bin);
            worst_softmin = -PROB.tau_cov * log( mean( exp( -rho_bin / PROB.tau_cov ) ) );
            rho = PROB.w_uniform * rho_uniform + PROB.w_worst * worst_softmin;
        end
    
        shortfall = max(0, PROB.rho_target - rho_bin);
        pen_cover = mean(shortfall.^2);
        J = (Fbar - PROB.lambda_cov * rho) + PROB.mu_cover * pen_cover;
    else
        % hard J, Fbar
        feas = all(G<0,2);
        if any(feas), Fbar = mean(F1(feas)); else, Fbar = inf; end
    
        if rhoKind=="binmean"
            rho = mean(perBinFast(BIN, double(feas)));
        else
            rho = mean(feas); % matches objective_hardmask
        end
    
        J = Fbar - PROB.lambda_cov * rho;
    end
end

function [Jrow,Frow,Rrow] = eval_row_taylor(v1_fixed, v2_grid, x_star, i1, i2, ...
    tau0,h10,h20,m0, dTau_dv1,dTau_dv2, dH1_dv1,dH1_dv2, dH2_dv1,dH2_dv2, dM_dv1,dM_dv2, ...
    th1,th2,th3, PROB, BIN, useHard, rhoKind)

    dv1 = (v1_fixed - x_star(i1));
    Jrow = zeros(1,numel(v2_grid)); Frow = Jrow; Rrow = Jrow;
    
    for jj = 1:numel(v2_grid)
        dv2 = (v2_grid(jj) - x_star(i2));
    
        % First-order approx of model fields at (v1_fixed, v2_grid(jj))
        tau = tau0 + dTau_dv1*dv1 + dTau_dv2*dv2;
        h1  = h10  + dH1_dv1 *dv1 + dH1_dv2 *dv2;
        h2  = h20  + dH2_dv1 *dv1 + dH2_dv2 *dv2;
        m   = m0   + dM_dv1  *dv1 + dM_dv2  *dv2;
    
        contact = sqrt(h1.^2 + h2.^2);
        F1 = tau ./ max(contact, PROB.eps_contact);
        G  = [tau - PROB.C1, -tau, PROB.eps_contact - contact, -m, m - PROB.C6];
    
        if ~useHard
            W = soft_weights_from_G(G, PROB.gamma);
            sumW = sum(W) + PROB.eps_den;
            Fbar = sum(W .* F1) / sumW;
    
            rho_bin = perBinFast(BIN, W);
            if rhoKind=="binmean"
                rho = mean(rho_bin);
            else
                rho_uniform   = mean(rho_bin);
                worst_softmin = -PROB.tau_cov * log( mean( exp( -rho_bin / PROB.tau_cov ) ) );
                rho = PROB.w_uniform * rho_uniform + PROB.w_worst * worst_softmin;
            end
    
            shortfall = max(0, PROB.rho_target - rho_bin);
            pen_cover = mean(shortfall.^2);
            J = (Fbar - PROB.lambda_cov * rho) + PROB.mu_cover * pen_cover;
        else
            feas = all(G<0,2);
            if any(feas), Fbar = mean(F1(feas)); else, Fbar = inf; end
            if rhoKind=="binmean", rho = mean(perBinFast(BIN, double(feas))); else, rho = mean(feas); end
            J = Fbar - PROB.lambda_cov * rho;
        end
    
        Jrow(jj)=J; Frow(jj)=Fbar; Rrow(jj)=rho;
    end
end

function [J,Fbar,rho] = local_eval_exact(x, th1,th2,th3, PROB, BIN, useHard, rhoKind)
    [tau,h1,h2,m] = updateFinger(th1, th2, th3, x(1), x(2), x(3), x(4));
    contact = sqrt(h1.^2 + h2.^2);
    F1 = tau ./ max(contact, PROB.eps_contact);
    G  = [tau - PROB.C1, -tau, PROB.eps_contact - contact, -m, m - PROB.C6];
    
    if ~useHard
        W = soft_weights_from_G(G, PROB.gamma);
        sumW = sum(W) + PROB.eps_den;
        Fbar = sum(W .* F1) / sumW;
    
        rho_bin = perBinFast(BIN, W);
        if rhoKind=="binmean"
            rho = mean(rho_bin);
        else
            rho_uniform   = mean(rho_bin);
            worst_softmin = -PROB.tau_cov * log( mean( exp( -rho_bin / PROB.tau_cov ) ) );
            rho = PROB.w_uniform * rho_uniform + PROB.w_worst * worst_softmin;
        end
    
        shortfall = max(0, PROB.rho_target - rho_bin);
        pen_cover = mean(shortfall.^2);
        J = (Fbar - PROB.lambda_cov * rho) + PROB.mu_cover * pen_cover;
    else
        feas = all(G<0,2);
        if any(feas), Fbar = mean(F1(feas)); else, Fbar = inf; end
        if rhoKind=="binmean", rho = mean(perBinFast(BIN, double(feas))); else, rho = mean(feas); end
        J = Fbar - PROB.lambda_cov * rho;
    end
end

function BIN = make_bin_index(THETA, nb1, nb2, r1, r2)
    % Precompute bin indices & counts so we avoid discretize() every time.
    th1 = THETA.samples(:,1); th2 = THETA.samples(:,2);
    if nargin < 4 || isempty(r1), r1 = [min(th1), max(th1)]; end
    if nargin < 5 || isempty(r2), r2 = [min(th2), max(th2)]; end
    e1 = linspace(r1(1), r1(2), nb1+1); e1(end)=e1(end)+1e-12;
    e2 = linspace(r2(1), r2(2), nb2+1); e2(end)=e2(end)+1e-12;
    b1 = discretize(th1, e1); b2 = discretize(th2, e2);
    good = ~isnan(b1) & ~isnan(b2);
    idx = sub2ind([nb1,nb2], b1(good), b2(good));
    cnt = accumarray(idx, 1, [nb1*nb2,1], @sum, 0);
    
    BIN.idx   = idx;           % N_good x 1 linear bin ids
    BIN.good  = good;          % N x 1 logical
    BIN.cnt   = cnt;           % (nb1*nb2) x 1 counts per bin
    BIN.nb1   = nb1; BIN.nb2 = nb2;
    BIN.edges = {e1, e2};
    BIN.r1    = r1; BIN.r2 = r2;
end

function rho_bin = perBinFast(BIN, W)
    % Mean per bin using precomputed indices/counts.
    sum_w = accumarray(BIN.idx, W(BIN.good), [BIN.nb1*BIN.nb2, 1], @sum, 0);
    rho_bin = zeros(BIN.nb1*BIN.nb2,1);
    nz = BIN.cnt > 0; rho_bin(nz) = sum_w(nz) ./ BIN.cnt(nz);
end

function W = soft_weights_from_G(G, gamma)
    % Geometric-mean of logistics to avoid underflow (matches your objective)
    Z = -G ./ gamma;
    S = 1 ./ (1 + exp(-Z));
    W = exp( mean(log(max(S, realmin)), 2) );
end

function keep_idx = pick_reps_per_bin(BIN, kPerBin)
    % Up to kPerBin representatives per bin (uniformly spaced among members).
    good_idx = find(BIN.good);
    keep_idx = [];
    for b = 1:numel(BIN.cnt)
        if BIN.cnt(b) == 0, continue; end
        pos = find(BIN.idx == b);                   % positions within "good" subset
        sel = pick_k_spaced(pos, min(kPerBin, numel(pos)));
        keep_idx = [keep_idx; good_idx(sel)]; %#ok<AGROW>
    end
    keep_idx = sort(keep_idx);
end

function sel = pick_k_spaced(vec, k)
    % deterministically spread k indices across vec
    n = numel(vec);
    if k>=n, sel = vec; return; end
    loc = round(linspace(1,n,k));
    sel = vec(loc(:));
end