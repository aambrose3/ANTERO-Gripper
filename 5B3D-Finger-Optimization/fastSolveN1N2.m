function [N1, N2] = fastSolveN1N2(J1, J2, JC, n1, n2, nC, F_cmp, K2, th2, l2_r, clampNonneg)
% Solve (J1' n1) N1 + (J2' n2) N2 + [0; K2*(th2 - l2_r)] + (JC' nC) F_cmp = 0
% for N1, N2. Works for scalar or batched inputs (pages).
%
% Inputs (per sample or batched across 3rd dim / 2nd dim for vectors):
%   J1,J2,JC : 2x2  or  2x2xN
%   n1,n2,nC : 2x1  or  2xN
%   F_cmp    : scalar or Nx1
%   K2, th2  : scalar or Nx1
%   l2_r     : scalar
%   clampNonneg (optional, default=true): clamp N1,N2 >= 0
%
% Outputs:
%   N1,N2 : column vectors (N×1)

    if nargin < 11 || isempty(clampNonneg), clampNonneg = false; end

    % Normalize shapes
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
    Tau_h = [zeros(1,N); -(K2(:).').*(th2(:).' - l2_r)];    % 2xN **negative torque when t2 is positive
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
