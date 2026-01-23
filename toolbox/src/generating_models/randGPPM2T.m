function [W, meta] = randGPPM2T(N, p, T_topo, T_wt, varargin)
%RANDGPPM2T  GPPM-style random directed graph with separate topology and weight coherence.
%
%   W = randGPPM2T(N, p, T_topo, T_wt)
%   [W, meta] = randGPPM2T(N, p, T_topo, T_wt, 'Name', value, ...)
%
% Idea:
%   - Topology (edge existence) is generated using a GPPM-style propensity based on
%     trophic distance, with temperature T_topo.
%   - Weights are then assigned with a separate temperature T_wt, so weights can be
%     more (or less) coherent than the topology.
%
% Propensity:
%   P_ij(T) = exp(-((h_i - h_j - 1)^2) / (2 T^2)), with T=0 as a hard constraint.
%
% INPUTS
%   N      : number of nodes (positive integer)
%   p      : target density in [0,1] (expected fraction of possible directed edges,
%            excluding self-loops by default)
%   T_topo : temperature >= 0 controlling edge existence (smaller => more coherent)
%   T_wt   : temperature >= 0 controlling weight coherence (smaller => more coherent)
%
% NAME–VALUE OPTIONS
%   'NumBasal'          (0)        : number of basal nodes B. If 0, uses B=1.
%   'AllowSelfLoops'    (false)    : allow i->i edges
%   'ForceWeakConn'     (true)     : ensure weak connectivity via skeleton stage
%   'ExactEdges'        ([])       : if provided, use this exact number of edges L instead of p
%   'Seed'              ([])       : RNG seed (scalar) for reproducibility
%
%   'Weighted'          (false)    : if false, returns logical adjacency (topology only)
%   'WeightDist'        ('uniform'): 'uniform'|'lognormal'|'exponential'|'constant'
%   'WeightParams'      ([0 1])    : parameters for WeightDist
%
%   'WeightModel'       (auto)     : 'independent' | 'meanProportionalToP' | 'scaledByP'
%                                   Default when Weighted=true and not provided:
%                                     WeightModel = 'meanProportionalToP'
%   'WeightPExponent'   (1)        : alpha >= 0 used in coupling factor P.^alpha
%
%   'LevelFun'          (@trophic_levels) : function handle: h = LevelFun(A)
%                                          should return an N-vector of levels
%   'LevelFrom'         ('topology'): 'topology' | 'final'
%                                   If 'topology' (default), compute levels on the
%                                   binary skeleton+topology before weighting.
%                                   If 'final', compute levels on the final weighted W
%                                   and use that only for meta output (no iteration).
%
% OUTPUTS
%   W    : logical (if Weighted=false) or double weighted adjacency
%   meta : struct describing parameters and draw
%
% Notes
%   - This generator does NOT iterate levels/weights to self-consistency. That would
%     be a different model (can be added later if desired).
%   - Basal nodes are enforced structurally: no incoming edges into basal nodes.
%
% (GPPM-style skeleton + propensity sampling, with two temperatures.)

% -------------------- parse --------------------
ip = inputParser;
ip.addParameter('NumBasal', 0, @(x)isnumeric(x) && isscalar(x) && x>=0);
ip.addParameter('AllowSelfLoops', false, @(x)islogical(x) && isscalar(x));
ip.addParameter('ForceWeakConn', true, @(x)islogical(x) && isscalar(x));
ip.addParameter('ExactEdges', [], @(x)isempty(x) || (isnumeric(x) && isscalar(x) && x>=0));
ip.addParameter('Seed', [], @(x)isempty(x) || (isnumeric(x) && isscalar(x)));

ip.addParameter('Weighted', false, @(x)islogical(x) && isscalar(x));
ip.addParameter('WeightDist', 'uniform', @(x)ischar(x) || isstring(x));
ip.addParameter('WeightParams', [0 1], @(x)isnumeric(x) && ~isempty(x));
ip.addParameter('WeightModel', 'independent', @(x)ischar(x) || isstring(x));
ip.addParameter('WeightPExponent', 1, @(x)isnumeric(x) && isscalar(x) && x>=0);

ip.addParameter('LevelFun', @trophic_levels, @(f)isa(f,'function_handle'));
ip.addParameter('LevelFrom', 'topology', @(x)ischar(x) || isstring(x));

ip.parse(varargin{:});
opt = ip.Results;

% Default logic: if Weighted and WeightModel not supplied, default to mean-coupled
if opt.Weighted && any(strcmp(ip.UsingDefaults, 'WeightModel'))
    opt.WeightModel = "meanProportionalToP";
end

validateattributes(N, {'numeric'}, {'scalar','integer','>=',1});
validateattributes(p, {'numeric'}, {'scalar','>=',0,'<=',1});
validateattributes(T_topo, {'numeric'}, {'scalar','>=',0});
validateattributes(T_wt, {'numeric'}, {'scalar','>=',0});

if ~isempty(opt.Seed)
    rng(opt.Seed);
end

% Choose B
B = opt.NumBasal;
if B == 0
    B = 1;
end
B = min(B, N);
basal = (1:B).';

% Target number of edges L
if ~isempty(opt.ExactEdges)
    L = round(opt.ExactEdges);
else
    M = N*N;
    if ~opt.AllowSelfLoops
        M = N*(N-1);
    end
    L = round(p * M);
end
L = max(L, 0);

% -------------------- stage 1: skeleton (optional weak connectivity) --------------------
A = zeros(N, N); % binary adjacency during construction (A(to,from)=1)

if opt.ForceWeakConn
    existing = basal.'; % row vector
    for newNode = (B+1):N
        src = existing(randi(numel(existing)));
        A(newNode, src) = 1; % src -> newNode
        existing(end+1) = newNode; %#ok<AGROW>
    end
end

% If L smaller than skeleton, truncate and return
L_skel = nnz(A);
if L <= L_skel
    if L < L_skel
        [ii,jj] = find(A);
        keep = randperm(L_skel, L);
        A(:) = 0;
        A(sub2ind([N N], ii(keep), jj(keep))) = 1;
    end

    if opt.Weighted
        % With too few edges, propensity isn't meaningful; treat as independent unless user insists otherwise
        Pwt = ones(N);
        W = assignWeightsFromProp(A, opt.WeightDist, opt.WeightParams, opt.WeightModel, opt.WeightPExponent, Pwt);
    else
        W = logical(A);
    end

    meta = makeMeta(N,B,L,T_topo,T_wt,opt,W,[]);
    return;
end

% -------------------- levels from current skeleton (for propensity matrices) --------------------
h = opt.LevelFun(A);

% Precompute trophic distance matrix x_hat = h_i - h_j (target i, source j)
[hi, hj] = ndgrid(h, h);
x_hat = hi - hj;

% -------------------- stage 2: topology edge addition using T_topo --------------------
Ptop = propensityFromXhat(x_hat, T_topo);

% Basal constraint: no incoming edges to basal targets
Ptop(basal,:) = 0;

% Self-loops
if ~opt.AllowSelfLoops
    Ptop(1:N+1:end) = 0;
end

% Remove already present edges
Ptop(A ~= 0) = 0;

% Fallback if empty (rare)
if all(Ptop(:) == 0)
    Ptop = ones(N);
    Ptop(basal,:) = 0;
    if ~opt.AllowSelfLoops
        Ptop(1:N+1:end) = 0;
    end
    Ptop(A ~= 0) = 0;
end

need = L - nnz(A);
if need > 0
    A = addEdgesWeightedWithoutReplacement(A, Ptop, need);
end

% -------------------- weights using T_wt --------------------
if opt.Weighted
    % Build a weight propensity matrix using the SAME levels (by design)
    Pwt = propensityFromXhat(x_hat, T_wt);

    % Apply same admissibility constraints (for weight coupling only)
    Pwt(basal,:) = 0;
    if ~opt.AllowSelfLoops
        Pwt(1:N+1:end) = 0;
    end

    % Assign weights on existing edges of A, with coupling to Pwt
    W = assignWeightsFromProp(A, opt.WeightDist, opt.WeightParams, opt.WeightModel, opt.WeightPExponent, Pwt);

else
    W = logical(A);
    Pwt = [];
end

% -------------------- optional: compute final levels for meta only --------------------
levelFrom = lower(string(opt.LevelFrom));
if levelFrom == "final"
    try
        h_final = opt.LevelFun(W);
    catch
        h_final = [];
    end
else
    h_final = [];
end

meta = makeMeta(N,B,L,T_topo,T_wt,opt,W,h_final);

end

% ======================================================================
function P = propensityFromXhat(x_hat, T)
% P = exp(-((x_hat-1)^2)/(2T^2)), with T=0 -> hard constraint.
if T == 0
    tol = 1e-12;
    P = double(abs(x_hat - 1) < tol);
else
    P = exp(-((x_hat - 1).^2) ./ (2*T^2));
end
end

% ======================================================================
function A = addEdgesWeightedWithoutReplacement(A, P, k)
idx = find(P(:) > 0);
wts = P(idx);

m = numel(idx);
k = min(k, m);

% Exponential-keys sampling (robust weighted sampling without replacement)
U = rand(m,1);
keys = -log(U) ./ wts;
[~, ord] = mink(keys, k);
pick = idx(ord);

A(pick) = 1;
end

% ======================================================================
function Ww = assignWeightsFromProp(A, dist, params, weightModel, alpha, P)
% Assign positive weights to existing edges in A.
% weightModel:
%   'independent'         : weights independent of P
%   'meanProportionalToP' : soft coupling via mean ~ P^alpha (stabilised)
%   'scaledByP'           : hard scaling: w = base * P^alpha

[i,j] = find(A);
m = numel(i);

dist = lower(string(dist));
weightModel = lower(string(weightModel));

% ---- base weights ----
switch dist
    case "uniform"
        if numel(params) < 2
            error('randGPPM2T:WeightParams', 'Uniform weights require WeightParams = [a b].');
        end
        a = params(1); b = params(2);
        base = a + (b - a) * rand(m,1);

    case "exponential"
        lambda = params(1);
        if lambda <= 0
            error('randGPPM2T:WeightParams', 'Exponential weights require lambda > 0.');
        end
        base = exprnd(1/lambda, m, 1);

    case "lognormal"
        if numel(params) < 2
            error('randGPPM2T:WeightParams', 'Lognormal weights require WeightParams = [mu sigma].');
        end
        mu = params(1); sigma = params(2);
        if sigma < 0
            error('randGPPM2T:WeightParams', 'Lognormal weights require sigma >= 0.');
        end
        base = lognrnd(mu, sigma, m, 1);

    case "constant"
        base = params(1) * ones(m,1);

    otherwise
        error('randGPPM2T:WeightDist', 'Unknown WeightDist "%s".', dist);
end

% ---- coupling ----
if weightModel == "independent"
    w = base;
else
    pij = P(sub2ind(size(P), i, j));
    pij = max(pij, 0);
    f = pij .^ alpha;

    if weightModel == "scaledbyp"
        w = base .* f;

    elseif weightModel == "meanproportionaltop"
        % Soft coupling: scale by f but stabilise global mean scale so varying T_wt
        % doesn't trivially shrink/expand all weights.
        w = base .* f;

        mw = mean(w);
        mb = mean(base);
        if mw > 0 && mb > 0
            w = w * (mb / mw);
        end

    else
        error('randGPPM2T:WeightModel', 'Unknown WeightModel "%s".', weightModel);
    end
end

Ww = zeros(size(A));
Ww(sub2ind(size(A), i, j)) = w;
end

% ======================================================================
function meta = makeMeta(N,B,L,T_topo,T_wt,opt,W,h_final)
meta = struct();
meta.model = "randGPPM2T";
meta.N = N;
meta.B = B;
meta.L = L;
meta.p_target = [];
meta.T_topo = T_topo;
meta.T_wt = T_wt;

meta.weighted = opt.Weighted;
meta.weightDist = string(opt.WeightDist);
meta.weightModel = string(opt.WeightModel);
meta.weightPExponent = opt.WeightPExponent;

meta.allowSelfLoops = opt.AllowSelfLoops;
meta.forceWeakConn = opt.ForceWeakConn;
meta.seed = opt.Seed;

meta.nnz = nnz(W);
meta.h_final = h_final;
end
