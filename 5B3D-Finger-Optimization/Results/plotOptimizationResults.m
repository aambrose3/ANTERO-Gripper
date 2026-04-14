%% fast_visualize_slices.m
% Super fast Ω slices:
%  - Parallelized across slice points (parfor, threads pool recommended)
%  - Decimates theta grids for speed (configurable)
%  - Streams θ3 in blocks to keep memory low

clear; clc;
%% Run this once if the wrong pool is opening
% close any existing pool
p = gcp('nocreate'); if ~isempty(p), delete(p); end
%% Try custom par pool
c = parcluster('local');
c.NumWorkers = 16;
saveProfile(c);
parpool(c, 16);

pctRunOnAll maxNumCompThreads(1) ;

%% ---- Load .mat ----
[fn, fp] = uigetfile('*.mat','Select saved optimization .mat');
if isequal(fn,0), error('No file selected.'); end
S = load(fullfile(fp,fn));
fprintf('Loaded %s\n', fullfile(fp,fn));

%% ---- Choose angle grids from file (prefer final > ref > base) ----
if all(isfield(S, {'th1_final','th2_final','th3_final'}))
    th1 = S.th1_final(:); th2 = S.th2_final(:); th3 = S.th3_final(:);
    grid_name = 'final';
elseif all(isfield(S, {'th1_ref','th2_ref','th3_ref'}))
    th1 = S.th1_ref(:); th2 = S.th2_ref(:); th3 = S.th3_ref(:);
    grid_name = 'ref';
else
    th1 = S.theta_1(:);  th2 = S.theta_2(:);  th3 = S.theta_3(:);
    grid_name = 'base';
end
fprintf('θ-grid: %s  |θ1|=%d, |θ2|=%d, |θ3|=%d\n', grid_name, numel(th1), numel(th2), numel(th3));

%% ---- Config: SPEED knobs ----
% Decimate the θ grids for visualization only (bigger stride => faster)
stride.t1 = 4;   % 1 = no decimation, 2 = keep every 2nd, etc.
stride.t2 = 4;
stride.t3 = 4;   % θ3 is often the biggest cost; feel free to go 4–10

% Slice grid sizes (coarse for speed, then we can upsample visually)
Nd = round((S.ub(1) - S.lb(1))*10^3)+1; Nc  = round((S.ub(2) - S.lb(2))*10^3)+1;   % d × c
Nk1= round((S.ub(3) - S.lb(3))*4)+1; Nk2 = Nk1;   % k1 × k2

% θ3 streaming block size
B = 128;            % 64–256 is fine

% Use threads pool if available (no-op if already open)
p = gcp('nocreate'); if isempty(p), try parpool("threads"); end, end

%% ---- Apply decimation for viz ----
th1v = th1(1:stride.t1:end);
th2v = th2(1:stride.t2:end);
th3v = th3(1:stride.t3:end);
fprintf('Using decimated θ-grids: |θ1|=%d, |θ2|=%d, |θ3|=%d\n', numel(th1v), numel(th2v), numel(th3v));

%% ---- Bounds and base design ----
lb = S.lb(:); ub = S.ub(:);
if isfield(S,'x_ref'), x_base = S.x_ref(:);
elseif isfield(S,'x_best'), x_base = S.x_best(:);
else, error('No x_ref or x_best found in file.\n'); end
fprintf('\nOptimal Design Parameters:\nc=%.6f m \nd=%.6f m \nk1=%.6f N/mm \nk2=%.6f Nm/rad\n\n', ...
    x_base(2), x_base(1), x_base(3), x_base(4))

%% ---- Precompute θ1×θ2 cross (once), and constants for workers ----
[TH1xy, TH2xy] = ndgrid(th1v, th2v);
th1xy = TH1xy(:); th2xy = TH2xy(:);
th1xyC = parallel.pool.Constant(th1xy);
th2xyC = parallel.pool.Constant(th2xy);
th3C   = parallel.pool.Constant(th3v);


%% ---- Local evaluator (memory-light, streams θ3) ----
function val = evalOmega(x, th1xy, th2xy, th3, B, W)
% evalOmega — objective with domain-wide constraints (post-contact span)
% x      : [d;c;k1;k2]
% th1xy  : (M12×1) vector of θ1 replicated across θ2 (ndgrid(:))
% th2xy  : (M12×1) vector of θ2 replicated across θ1 (ndgrid(:))
% th3    : (K×1)   vector of θ3 samples (monotonic)
% B      : block size for streaming θ3 (e.g., 128–256)

    % ------ Design vars ------
    d = x(1); c = x(2); k1 = x(3); k2 = x(4);

    % ------ Objective weights & pointwise limits (match your optimizer) ------
    w1 = W(1); w2 = W(2);
    eps_denom = 1e-12;
    tau_max = 5.0; m_max = 13; push_only = true;

    % ------ Domain-wide thresholds ------
    Delta1 = deg2rad(3);     % post-contact contiguous span requirement
    Delta2 = deg2rad(3);     % must become feasible within Δ2 after contact
    p_required = 0.50;       % fraction of (θ1,θ2) rows that must pass
    BIG = 1e12;              % penalty if fraction < p_required (for plotting)

    % ------ Shapes ------
    th1xy = th1xy(:); th2xy = th2xy(:); th3 = th3(:);
    M12 = numel(th1xy); K = numel(th3);
    if M12==0 || K==0, val = BIG; return; end

    % ------ Per-row accumulators for Ω across θ3 (only count feasible cells) ------
    sumOmega = zeros(M12,1,'double');
    cntOmega = zeros(M12,1,'uint32');

    % ------ Per-row state for domain-wide constraints ------
    firstContactI  = zeros(M12,1,'uint32');  % k_c: first m>0
    firstFeasibleI = zeros(M12,1,'uint32');  % k_f: first feasible after contact
    h2_failed      = false(M12,1);           % true if gap>Δ2 or never feasible after contact

    % Run tracking ONLY for the span that starts at k_f
    inRun      = false(M12,1);
    runStartI  = zeros(M12,1,'uint32');
    runDone    = false(M12,1);               % stop after first break (ignore later islands)
    span       = zeros(M12,1,'double');

    % Active rows are those that can still change outcome
    active = true(M12,1);

    % ------ Stream θ3 in blocks ------
    for ks = 1:B:K
        ke  = min(ks+B-1, K);
        Kb  = ke-ks+1;
        t3b = th3(ks:ke);

        % Build block vectors
        th1_vec = repmat(th1xy, Kb, 1);
        th2_vec = repmat(th2xy, Kb, 1);
        th3_vec = kron(t3b,   ones(M12,1));

        % Model
        [tau, N1, N2, m] = updateFinger(th1_vec, th2_vec, th3_vec, d, c, k1, k2);
        tau = reshape(tau, M12, Kb);
        N1  = reshape(N1 , M12, Kb);
        N2  = reshape(N2 , M12, Kb);
        m   = reshape(m  , M12, Kb);

        % Pointwise feasibility for this block
        passBlk = (abs(tau) <= tau_max) & ((~push_only) | (tau >= 0)) & ...
                  (N1 >= 0) & (N2 >= 0) & (m >= 0) & (m <= m_max);

        % ---- Accumulate Ω for pointwise-feasible cells (all rows) ----
        denom = max(w1.*(N1.^2) + w2.*(N2.^2), eps_denom);
        Om    = (tau.^2) ./ denom;
        Om(~passBlk) = 0;
        sumOmega = sumOmega + sum(Om, 2);
        cntOmega = cntOmega + uint32(sum(passBlk, 2));

        % ---- Update domain-wide state column-by-column ----
        for j = 1:Kb
            kj = ks + j - 1;
            if ~any(active); break; end

            cur_m    = m(:,j);
            cur_pass = passBlk(:,j);

            % (A) First contact (m>0)
            needContact = active & (firstContactI==0);
            gotContact  = needContact & (cur_m > 0);
            if any(gotContact)
                firstContactI(gotContact) = uint32(kj);
            end

            % (B) First feasible after contact (also require m>0 to avoid pre-contact)
            eligibleFeas = active & (firstContactI>0) & (firstFeasibleI==0) & (cur_m > 0);
            gotFeasible  = eligibleFeas & cur_pass;
            if any(gotFeasible)
                firstFeasibleI(gotFeasible) = uint32(kj);
            end

            % (C) h2: enforce gap ≤ Δ2; early fail if exceeded
            checkGap = active & (firstContactI>0) & (firstFeasibleI==0);
            if any(checkGap)
                gapNow = th3(kj) - th3(double(firstContactI(checkGap)));
                exceeded = gapNow > Delta2;
                if any(exceeded)
                    idx = find(checkGap);
                    h2_failed(idx(exceeded)) = true;
                    active(idx(exceeded))    = false;
                    inRun(idx(exceeded))     = false;
                end
            end

            % (D) Track ONLY the span that begins at k_f
            track = active & (firstFeasibleI>0) & ~runDone;
            if any(track)
                atOrAfterKF = (uint32(kj) >= firstFeasibleI);
                cur = track & cur_pass & (cur_m > 0) & atOrAfterKF;

                % start the run exactly at k_f
                needStart = track & ~inRun & (uint32(kj) == firstFeasibleI);
                if any(needStart)
                    runStartI(needStart) = uint32(kj);
                    inRun(needStart)     = true;
                end

                % end the run on first break -> finalize span and stop tracking this row
                ended = track & inRun & ~cur;
                if any(ended)
                    eIdx = kj - 1;
                    sIdx = double(runStartI(ended));
                    seg  = th3(eIdx) - th3(sIdx);
                    span(ended)    = seg;       % only the first post-contact island
                    inRun(ended)   = false;
                    runDone(ended) = true;
                    active(ended)  = false;
                end
            end
        end
    end

    % Close any runs that lasted to the final θ3
    stillOpen = inRun & ~runDone;
    if any(stillOpen)
        sIdx = double(runStartI(stillOpen));
        eIdx = K;
        seg  = th3(eIdx) - th3(sIdx);
        span(stillOpen)   = seg;
        runDone(stillOpen)= true;
        active(stillOpen) = false;
    end

    % Rows with no contact or no feasible-after-contact => h2 fail
    noContact         = (firstContactI==0);
    hadContactNoFeas  = (firstContactI>0) & (firstFeasibleI==0);
    h2_failed(noContact | hadContactNoFeas) = true;

    % ---- Row pass/fail and final objective ----
    rowPass = (~h2_failed) & (span >= Delta1);

    % Enforce the same global requirement used in optimization
    frac = mean(rowPass);
    if frac < p_required
        val = BIG;      % visualize as large (violated domain-wide constraint)
        return;
    end

    % Compute Ω only over rows that pass both h1 & h2
    goodRows = rowPass & (cntOmega > 0);
    if ~any(goodRows)
        val = BIG;
        return;
    end
    rowAvg = sumOmega(goodRows) ./ double(cntOmega(goodRows));
    val    = mean(rowAvg, 'omitnan');
    if ~isfinite(val), val = BIG; end
end


%% Finger Model (called from objective function):
function [tau, N1, N2, m] = updateFinger(th1, th2, th3, d, c, k1, k2)
    th1 = th1(:); th2 = th2(:); th3 =th3(:); % enforce column vectors
    % Project-specific function calls (presumed available in your codebase)
    N = numel(th1);
    if ~(numel(th2)==N && numel(th3)==N), error('updateFinger: angle vectors must have equal length.'); end
    optFlag = true; % set optimization flag to true
    param = getParam(optFlag, d, c, k2); % Not vector safe
    param = updateParam(param, th1, th2, th3); % Vector safe
    Jac   = getJacobians(param, th1, th2, th3); % Vector safe
    [cmp, F_cmp, TA, ~] = getSpring(param, th1, th2, th3, k1, true); % Vector safe (cmp in mm, F_cmp in N, TA in Nm, k1 in N/mm)
    [N1v, N2v] = fastSolveN1N2(Jac.J1, Jac.J2, Jac.JC, ... % numerically solve the system of equations
                                  param.n1, param.n2, param.nC, ...
                                  F_cmp, param.K2, th2, param.l2_r, true);
    tau = TA(:); N1 = N1v(:); N2 = N2v(:); m = cmp(:); % output vectors
end


%% ---- Build the two slice meshes ----
d_vec  = linspace(lb(1), ub(1), Nd);
c_vec  = linspace(lb(2), ub(2), Nc);
k1_vec = linspace(lb(3), ub(3), Nk1);
k2_vec = linspace(lb(4), ub(4), Nk2);

[DD, CC]   = meshgrid(d_vec,  c_vec);
[KK1, KK2] = meshgrid(k1_vec, k2_vec);

% Pass a few things into local functions for each worker
if isfield(S,'w1') && isfield(S,'w2'), W = [S.w1,S.w2]; 
else, W = [0.25, 0.75]; end


%% ---- Evaluate Ω(d,c) slice in parallel (fix k1,k2 at base) ----
fprintf('Evaluating Ω(d,c) slice in parallel...\n');
Z_dc = nan(numel(CC),1);
parfor q = 1:numel(CC)
    x = [DD(q); CC(q); x_base(3); x_base(4)];
    Z_dc(q) = evalOmega([DD(q); CC(q); x_base(3); x_base(4)], th1xy, th2xy, th3v, B, W);
end
Z_dc = reshape(Z_dc, size(CC));

%% ---- Evaluate Ω(k1,k2) slice in parallel (fix d,c at base) ----
fprintf('Evaluating Ω(k1,k2) slice in parallel...\n');
Z_k = nan(numel(KK1),1);
parfor q = 1:numel(KK1)
    x = [x_base(1); x_base(2); KK1(q); KK2(q)];
    Z_k(q)  = evalOmega([x_base(1); x_base(2); KK1(q); KK2(q)], th1xy, th2xy, th3v, B, W);
end
Z_k = reshape(Z_k, size(KK1));

%% ---- (Optional) visual upsampling for smoothness ----
% Comment out if you prefer the raw coarse grid.
% fprintf('Upscaling for visuals...\n');
% ups = 2;  % 1 = off, 2 = ~2x finer display
% if ups>1
%     d_f = linspace(d_vec(1), d_vec(end), ups*numel(d_vec));
%     c_f = linspace(c_vec(1), c_vec(end), ups*numel(c_vec));
%     [Df, Cf] = meshgrid(d_f, c_f);
%     Z_dc = interp2(DD, CC, Z_dc, Df, Cf, 'linear');
%     DD = Df; CC = Cf;
% 
%     k1_f = linspace(k1_vec(1), k1_vec(end), ups*numel(k1_vec));
%     k2_f = linspace(k2_vec(1), k2_vec(end), ups*numel(k2_vec));
%     [K1f, K2f] = meshgrid(k1_f, k2_f);
%     Z_k = interp2(KK1, KK2, Z_k, K1f, K2f, 'linear');
%     KK1 = K1f; KK2 = K2f;
% end

%% ---- Plots ----
fprintf('Plotting...\n');
figure('Name','Ω slice: d vs c','Color','w');
contourf(DD, CC, Z_dc, 20, 'LineStyle','none'); colorbar
hold on; plot(x_base(1), x_base(2), 'kp', 'MarkerFaceColor','y','MarkerSize',9);
xlabel('d [m]'); ylabel('c [m]'); title(sprintf('Ω(d,c) | θ-grid:%s (decimated)', grid_name));
grid on; box on;

figure('Name','Ω slice: k1 vs k2','Color','w');
contourf(KK1, KK2, Z_k, 20, 'LineStyle','none'); colorbar
hold on; plot(x_base(3), x_base(4), 'kp', 'MarkerFaceColor','y','MarkerSize',9);
xlabel('k_1 [N/mm]'); ylabel('k_2 [N·m/rad]'); title(sprintf('Ω(k_1,k_2) | θ-grid:%s (decimated)', grid_name));
grid on; box on;
