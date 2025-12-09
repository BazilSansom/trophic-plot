function [X, Y, h] = trophicLayout(W, varargin)
%TROPHICLAYOUT  Soft-constrained trophic layout in 2D.
%
%   [X, Y, h] = TROPHICLAYOUT(W) computes trophic levels
%       h = trophic_levels(W)
%   and returns a 2D layout (X,Y) where:
%       - vertical coordinates are softly constrained to h,
%       - horizontal coordinates are determined by a force-directed
%         algorithm with attractive and repulsive forces,
%       - by default, the final Y are snapped exactly to the (possibly
%         rescaled) trophic levels used for anchoring.
%
%   The algorithm has two phases:
%     (1) a "free" force-directed phase (lambda = 0),
%     (2) a trophic anchoring phase where lambda is gradually increased.
%
%   By default, the initial horizontal coordinates X are computed using a
%   spectral initialisation (Fiedler vector of the Laplacian).
%
%   [X, Y, h] = TROPHICLAYOUT(W, 'Name', Value, ...)
%   supports the following Name/Value options:
%
%   ---- Hierarchy options ----
%   'hProvided'      : custom level vector h0 (N x 1). If empty or omitted,
%                      trophic levels h = trophic_levels(W) are computed
%                      internally.
%
%   'RescaleLevels'  : logical (default = true). If true, the levels h are
%                      centred and scaled to a standardised range for
%                      vertical plotting. The layout remains linear in h.
%
%   'SnapToLevels'   : logical (default = true). If true, after the
%                      optimisation loop finishes, the final Y are set to
%                      the anchoring levels h_plot (the possibly rescaled
%                      version of h). If false, Y remains the softly
%                      constrained result of the dynamics.
%
%   ---- Initialisation options ----
%   'InitMode'   : one of:
%         'spectral'         - deterministic Fiedler-vector init.
%         'random'           - X,Y initialised randomly.
%         'spectral+noise'   - spectral init plus small noise.
%
%   'Seed'       : integer or [] (default = []). If provided, the RNG is
%                  seeded so that random components become reproducible.
%
%   'InitialX'   : override for initial X positions (N x 1).
%   'InitialY'   : override for initial Y positions (N x 1).
%
%   ---- Force-layout parameters ----
%   'NumIters'        : total number of iterations (default = 1000).
%   'AttractStrength' : spring constant for attractive forces (default 0.01).
%   'RepelStrength'   : repulsion strength between all node pairs (default 0.1).
%   'StepSize'        : initial step size for gradient-like updates (default 0.1).
%   'Cooling'         : multiplicative decay of step size per iteration (default 0.99).
%   'Damping'         : velocity damping factor (default 0.9).
%   'LambdaMax'       : maximum strength of trophic anchoring (default 5).
%   'LambdaPower'     : exponent p in lambda(t) = LambdaMax * (t/T)^p (default 2).
%
%   ---- Node weighting ----
%   'NodeWeight'      : node weight vector u (N x 1). If empty (default),
%                       all nodes have weight 1 in the anchoring term.
%
%   Output:
%     X, Y  - final coordinates (N x 1 each).
%            If SnapToLevels = true (default), Y is exactly h_plot.
%     h     - trophic levels (unscaled, as returned by trophic_levels).
%
%   The layout is intended to be used together with TROPHIC_LEVELS and
%   PLOTTFL for visualising direction and feedback in weighted directed
%   networks.
%
% Author:
%   Bazil Sansom (Warwick Business School, University of Warwick)
%   Contact: bazil.sansom@wbs.ac.uk
%
% -------------------------------------------------------------------------

% ---- Input checks ----
if ~ismatrix(W) || size(W,1) ~= size(W,2)
    error('trophicLayout:WNotSquare', 'W must be an N x N adjacency matrix.');
end

% ---- Parse optional parameters ----
params = parseInputsSoft(varargin{:});  % local helper function below

T          = params.NumIters;
k_a        = params.AttractStrength;
k_r        = params.RepelStrength;
eta        = params.StepSize;
cool       = params.Cooling;
alpha      = params.Damping;
u          = params.NodeWeight;
x0         = params.InitialX;
y0         = params.InitialY;
hProvided  = params.hProvided;
doRescale  = params.RescaleLevels;
lam_max    = params.LambdaMax;
lam_pow    = params.LambdaPower;
snapToH    = params.SnapToLevels;

N = size(W,1);

% Default node weights if none given
if isempty(u)
    u = ones(N,1);
else
    u = u(:);
    if numel(u) ~= N
        error('trophicLayout:BadNodeWeight', ...
              'NodeWeight must be empty or of length N.');
    end
end

% ---- Optional local RNG seeding ----
if ~isempty(params.Seed)
    oldRngState = rng;      % save caller's RNG state
    rng(params.Seed);       % use our own local seed
else
    oldRngState = [];       % flag: no local seeding
end

% ----  Trophic levels ----
if ~isempty(hProvided)
    h = hProvided(:);
else
    % Use external trophic_levels function (canonical implementation)
    h = trophic_levels(W);
    h = h(:);
end

% ---- Rescale trophic levels for layout if requested ----
if doRescale
    wSum = sum(u);
    if wSum <= eps
        % fallback: unweighted mean/variance if all weights are zero
        h_centered = h - mean(h);
        s = std(h);
    else
        mu = sum(u .* h) / wSum;
        h_centered = h - mu;
        s2 = sum(u .* (h_centered.^2)) / wSum;
        s  = sqrt(max(s2, eps));
    end
    h_plot = h_centered / s;
else
    h_plot = h;
end

% ---- Initial positions (X, Y) according to InitMode ----

% X initialisation
if ~isempty(x0)
    X = x0(:);
    if numel(X) ~= N
        error('trophicLayout:BadInitialX', 'InitialX must be length N.');
    end
else
    switch params.InitMode
        case 'spectral'
            X = initX_spectral(W);

        case 'spectral+noise'
            X = initX_spectral(W);
            noiseScaleX = 0.3;
            X = X + noiseScaleX * randn(N,1);

        case 'random'
            X = 0.1 * randn(N,1);

        otherwise
            error('trophicLayout:UnknownInitMode', ...
                  'Unknown InitMode: %s', params.InitMode);
    end
end

% Y initialisation
if ~isempty(y0)
    Y = y0(:);
    if numel(Y) ~= N
        error('trophicLayout:BadInitialY', 'InitialY must be length N.');
    end
else
    switch params.InitMode
        case 'spectral'
            % start close to trophic levels (already rescaled h_plot)
            Y = h_plot;

        case 'spectral+noise'
            noiseScaleY = 0.2;
            Y = h_plot + noiseScaleY * randn(N,1);

        case 'random'
            Y = h_plot + 0.2 * randn(N,1);  % anchored around h, but looser

        otherwise
            Y = h_plot;
    end
end

% ---- Prepare adjacency for attractive forces (undirected view) ----
W_und = W + W.';    % symmetrised weights
[ei, ej, wij] = find(W_und);
E = numel(wij);

% ---- Main iteration loop ----
Vx = zeros(N,1);
Vy = zeros(N,1);

for t = 1:T
    % Annealed lambda for vertical anchoring
    if T > 1
        lam = lam_max * ((t-1) / (T-1))^lam_pow;
    else
        lam = lam_max;
    end

    % --- Attractive forces along edges (2D spring) ---
    Fx_attr = zeros(N,1);
    Fy_attr = zeros(N,1);
    for e = 1:E
        i    = ei(e);
        j    = ej(e);
        w_ij = wij(e);
        dx   = X(j) - X(i);
        dy   = Y(j) - Y(i);
        % simple linear spring, no rest length
        fx = k_a * w_ij * dx;
        fy = k_a * w_ij * dy;
        Fx_attr(i) = Fx_attr(i) + fx;
        Fx_attr(j) = Fx_attr(j) - fx;
        Fy_attr(i) = Fy_attr(i) + fy;
        Fy_attr(j) = Fy_attr(j) - fy;
    end

    % --- Repulsive forces between all pairs (Coulomb-like) ---
    DX = X - X.';           % N x N
    DY = Y - Y.';           % N x N
    R2 = DX.^2 + DY.^2;
    R2(1:N+1:end) = Inf;    % no self-force
    R  = sqrt(R2);
    Fmag = k_r ./ R2;       % magnitude ~ 1/r^2
    Fx_rep_mat = Fmag .* (DX ./ (R + eps));
    Fy_rep_mat = Fmag .* (DY ./ (R + eps));
    Fx_rep = sum(Fx_rep_mat, 2);
    Fy_rep = sum(Fy_rep_mat, 2);

    % --- Trophic anchoring force (only vertical) ---
    %   F_y_troph = - d/dY [ lam * sum u_i (Y_i - h_plot_i)^2 ]
    %              = - 2 * lam * u_i * (Y_i - h_plot_i)
    Fy_troph = - 2.0 * lam * u .* (Y - h_plot);

    % --- Total forces ---
    Fx = Fx_attr + Fx_rep;
    Fy = Fy_attr + Fy_rep + Fy_troph;

    % --- Velocity updates (damped) ---
    Vx = alpha * Vx + eta * Fx;
    Vy = alpha * Vy + eta * Fy;

    % Clip max step to avoid blow-ups
    maxStep = 1.0;
    Vx = max(min(Vx, maxStep), -maxStep);
    Vy = max(min(Vy, maxStep), -maxStep);

    % --- Position updates ---
    X = X + Vx;
    Y = Y + Vy;

    % --- Cooling ---
    eta = eta * cool;
end

% ---- Optional hard snap to trophic levels ----
if snapToH
    % Snap Y exactly to the anchoring levels used in the dynamics.
    % This works whether h_plot is rescaled or not.
    Y = h_plot;
end

% ---- Restore RNG state if we changed it ----
if ~isempty(params.Seed)
    rng(oldRngState);
end

end % main function


% =====================================================================
function params = parseInputsSoft(varargin)
%PARSEINPUTSSOFT  Parse name/value pairs for trophicLayout.

% Default parameters
params.NumIters        = 1000;
params.AttractStrength = 0.01;
params.RepelStrength   = 0.1;
params.StepSize        = 0.1;
params.Cooling         = 0.99;
params.Damping         = 0.9;
params.InitialX        = [];
params.InitialY        = [];
params.NodeWeight      = [];
params.RescaleLevels   = true;
params.hProvided       = [];
params.LambdaMax       = 5;
params.LambdaPower     = 2;
params.InitMode        = 'spectral+noise';  % 'spectral' | 'random' | 'spectral+noise'
params.Seed            = [];                % [] => don't touch RNG
params.SnapToLevels    = true;              % snap Y to h_plot at the end

if isempty(varargin)
    return;
end

if mod(numel(varargin), 2) ~= 0
    error('trophicLayout:BadArgs', ...
          'Optional arguments must be name/value pairs.');
end

for k = 1:2:numel(varargin)
    name  = varargin{k};
    value = varargin{k+1};

    if ~ischar(name) && ~isstring(name)
        error('trophicLayout:BadParamName', ...
              'Parameter names must be strings.');
    end

    switch lower(char(name))
        case 'numiters'
            params.NumIters = value;
        case 'attractstrength'
            params.AttractStrength = value;
        case 'repelstrength'
            params.RepelStrength = value;
        case 'stepsize'
            params.StepSize = value;
        case 'cooling'
            params.Cooling = value;
        case 'damping'
            params.Damping = value;
        case 'initialx'
            params.InitialX = value;
        case 'initialy'
            params.InitialY = value;
        case 'nodeweight'
            params.NodeWeight = value;
        case 'rescalelevels'
            params.RescaleLevels = logical(value);
        case 'lambdamax'
            params.LambdaMax = value;
        case 'lambdapower'
            params.LambdaPower = value;
        case 'initmode'
            params.InitMode = lower(char(value));
        case 'seed'
            params.Seed = value;
        case 'hprovided'
            params.hProvided = value;
        case 'snaptolevels'
            params.SnapToLevels = logical(value);
        otherwise
            error('trophicLayout:UnknownParam', ...
                  'Unknown parameter name: %s', name);
    end
end

end  % parseInputsSoft


% =====================================================================
function X0 = initX_spectral(W)
%INITX_SPECTRAL  Deterministic spectral initialisation for X-coordinate.
%
%   X0 = INITX_SPECTRAL(W) returns a 1D spectral embedding based on the
%   Fiedler vector of the symmetrised Laplacian, with deterministic sign
%   fixing. This gives stable left/right separation for the layout.

    % Undirected, unweighted structure
    A = (W + W.') ~= 0;
    A = double(A);

    % Combinatorial Laplacian
    d = sum(A, 2);
    L = diag(d) - A;
    L = 0.5 * (L + L.');   % ensure symmetry

    n = size(L,1);

    % Trivial tiny cases: no need for eigs
    if n <= 2
        % Just spread nodes out slightly on a line
        X0 = linspace(-0.5, 0.5, n).';
        return;
    end

    % --- Try robust eigs with a small positive shift ---
    opts = struct();
    opts.isreal = true;

    gotFiedler = false;
    f = [];

    try
        % Use shift-invert around a small positive sigma so that
        % (L - sigma*I) is nonsingular even though L is singular.
        sigma = 1e-3;
        [V, D] = eigs(L, 2, sigma, opts);   % eigenvalues closest to sigma
        lam    = diag(D);

        % Sort eigenvalues; we expect to recover 0 and lambda2
        [~, idx] = sort(lam, 'ascend');

        if numel(idx) >= 2
            f = V(:, idx(2));  % approximate Fiedler direction
            gotFiedler = true;
        end

    catch
        % Fall through to full eig if eigs fails for any reason
        gotFiedler = false;
    end

    % --- Fallback: full eig if eigs failed ---
    if ~gotFiedler
        Lfull = full(L);
        [Vfull, Dfull] = eig(Lfull);
        lamAll         = diag(Dfull);
        [~, idxAll]    = sort(lamAll, 'ascend');

        if numel(idxAll) >= 2
            f = Vfull(:, idxAll(2));   % exact Fiedler vector
        else
            % Extreme degenerate case: just random
            f = randn(n,1);
        end
    end

    % --- Normalise to roughly [-1,1] range ---
    f = f(:);
    f = f - mean(f);
    maxAbs = max(abs(f));
    if maxAbs > 0
        f = f / maxAbs;
    end

    % ---- Deterministic sign fixing ----
    s = sum(f);

    if abs(s) > 1e-12
        % Case 1: sum determines orientation
        if s < 0
            f = -f;
        end
    else
        % Case 2: symmetric case – use largest-magnitude entry
        [~, idxMax] = max(abs(f));
        if f(idxMax) < 0
            f = -f;
        end
    end

    X0 = f;
end
