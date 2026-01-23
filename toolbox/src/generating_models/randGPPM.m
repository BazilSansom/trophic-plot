function [W, meta] = randGPPM(N, p, T, varargin)
%randGPPM  Generate directed networks with tunable trophic coherence (GPPM-style).
%
%   W = randGPPM(N, p, T)
%   [W, meta] = randGPPM(N, p, T, 'Name', value, ...)
%
% Convention (matches MATLAB digraph):
%   W(i,j) ~= 0 encodes a directed edge i -> j.
%
% Two-stage construction (Klaise & Johnson, 2016 style):
%   1) Build a weakly connected directed skeleton
%   2) Compute preliminary levels h on the skeleton
%   3) Add remaining edges using propensity P_ij biased towards h(j)-h(i)≈1
%
% Options
%   'NumBasal'       (0)        : if 0 => B=1 basal node; else B basals
%   'AllowSelfLoops' (false)
%   'ForceWeakConn'  (true)
%   'ForceDAG'       (false)    : if true, restrict all added edges to follow
%                                a total order induced by preliminary h
%                                (guarantees acyclic output)
%   'ExactEdges'     ([])
%   'Seed'           ([])
%
%   'Weighted'       (false)
%   'WeightDist'     ('uniform') : 'uniform'|'lognormal'|'exponential'|'constant'
%   'WeightParams'   ([0 1])
%   'WeightModel'    (auto)      : 'independent'|'meanProportionalToP'|'scaledByP'
%   'WeightPExponent'(1)
%
%   'LevelFun'       (@trophic_levels)
%
% Outputs
%   W    : logical adjacency if unweighted; double weights if weighted
%   meta : struct with diagnostics, including dagOrder if ForceDAG=true

% -------------------- parse --------------------
ip = inputParser;
ip.addParameter('NumBasal', 0, @(x)isnumeric(x) && isscalar(x) && x>=0);
ip.addParameter('AllowSelfLoops', false, @(x)islogical(x) && isscalar(x));
ip.addParameter('ForceWeakConn', true, @(x)islogical(x) && isscalar(x));
ip.addParameter('ForceDAG', false, @(x)islogical(x) && isscalar(x));
ip.addParameter('ExactEdges', [], @(x)isempty(x) || (isnumeric(x) && isscalar(x) && x>=0));
ip.addParameter('Seed', [], @(x)isempty(x) || (isnumeric(x) && isscalar(x)));

ip.addParameter('Weighted', false, @(x)islogical(x) && isscalar(x));
ip.addParameter('WeightDist', 'uniform', @(x)ischar(x) || isstring(x));
ip.addParameter('WeightParams', [0 1], @(x)isnumeric(x) && ~isempty(x));
ip.addParameter('WeightModel', 'independent', @(x)ischar(x) || isstring(x));
ip.addParameter('WeightPExponent', 1, @(x)isnumeric(x) && isscalar(x) && x>=0);

ip.addParameter('LevelFun', @trophic_levels, @(f)isa(f,'function_handle'));

ip.parse(varargin{:});
opt = ip.Results;

% If Weighted=true and WeightModel not explicitly passed, default to mean-coupled
if opt.Weighted && any(strcmp(ip.UsingDefaults,'WeightModel'))
    opt.WeightModel = "meanProportionalToP";
end

validateattributes(N, {'numeric'}, {'scalar','integer','>=',1});
validateattributes(p, {'numeric'}, {'scalar','>=',0,'<=',1});
validateattributes(T, {'numeric'}, {'scalar','>=',0});

if ~isempty(opt.Seed)
    rng(opt.Seed);
end

% Basals
B = opt.NumBasal;
if B == 0, B = 1; end
B = min(B, N);
basal = (1:B).';

% Target number of edges
if ~isempty(opt.ExactEdges)
    L = round(opt.ExactEdges);
else
    if opt.AllowSelfLoops
        M = N*N;
    else
        M = N*(N-1);
    end
    L = round(p * M);
end
L = max(L, 0);

% -------------------- stage 1: skeleton --------------------
% Build a weakly connected skeleton with edges from an existing node to the new node
% using *row->col* convention: src -> newNode means W(src,newNode)=1.
W = zeros(N,N);

if opt.ForceWeakConn
    existing = basal(:).';  % row vector
    for newNode = (B+1):N
        src = existing(randi(numel(existing)));
        W(src, newNode) = 1;         % src -> newNode
        existing(end+1) = newNode; %#ok<AGROW>
    end
end

% If L smaller than skeleton, truncate target
L_skel = nnz(W);
if L <= L_skel
    if L < L_skel
        [ii,jj] = find(W);
        keep = randperm(L_skel, L);
        W(:) = 0;
        W(sub2ind([N N], ii(keep), jj(keep))) = 1;
    end

    if opt.Weighted
        Pprop = ones(N); % unused meaningfully here
        W = assignWeights(W, opt.WeightDist, opt.WeightParams, opt.WeightModel, opt.WeightPExponent, Pprop);
    else
        W = logical(W);
    end

    meta = makeMeta(N,B,L,T,opt,W,[]);
    return;
end

% -------------------- stage 2: coherence-biased edge addition --------------------
% Compute preliminary levels on skeleton (must use same convention as W)
h = opt.LevelFun(W);
h = h(:);

% Prefer edges with trophic distance ~ 1: h(j) - h(i) ≈ 1 for i->j
[hi, hj] = ndgrid(h, h);    % hi(i,j)=h(i), hj(i,j)=h(j)
d_hat = hj - hi;            % d_hat(i,j)=h(j)-h(i)

if T == 0
    tol = 1e-12;
    P = double(abs(d_hat - 1) < tol);
else
    P = exp(-((d_hat - 1).^2) ./ (2*T^2));
end

% -------------------- optional: enforce DAG (acyclic) --------------------
% Build a deterministic total order from h, tie-break by node index.
% Allow only edges that go forward in that order: rank(j) > rank(i).
dagOrder = [];
if opt.ForceDAG
    [~, perm] = sortrows([h(:), (1:N)']);  % deterministic total order
    rank = zeros(N,1);
    rank(perm) = 1:N;

    allow = (rank(:).' > rank(:));  % allow(i,j)=true iff rank(j) > rank(i)
    P(~allow) = 0;

    dagOrder = perm(:);
end

% Basal constraint: no incoming edges to basal nodes.
% Incoming to basal b are edges i->b => column b.
P(:, basal) = 0;

% Self-loops
if ~opt.AllowSelfLoops
    P(1:N+1:end) = 0;
end

% Remove already present edges
P(W ~= 0) = 0;

% Fallback if P is empty but we still need edges: uniform over admissible
if all(P(:) == 0)
    P = ones(N);

    % Reapply constraints
    P(:, basal) = 0;
    if ~opt.AllowSelfLoops
        P(1:N+1:end) = 0;
    end
    if opt.ForceDAG
        [~, perm] = sortrows([h(:), (1:N)']);
        rank = zeros(N,1);
        rank(perm) = 1:N;
        allow = (rank(:).' > rank(:));
        P(~allow) = 0;
        dagOrder = perm(:);
    end

    P(W ~= 0) = 0;
end

% Sample additional edges without replacement
need = L - nnz(W);
if need > 0
    W = addEdgesWeightedWithoutReplacement(W, P, need);
end

% -------------------- weights (optional) --------------------
if opt.Weighted
    W = assignWeights(W, opt.WeightDist, opt.WeightParams, opt.WeightModel, opt.WeightPExponent, P);
else
    W = logical(W);
end

meta = makeMeta(N,B,L,T,opt,W,dagOrder);

end

% ======================================================================
function W = addEdgesWeightedWithoutReplacement(W, P, k)
idx = find(P(:) > 0);
wts = P(idx);

m = numel(idx);
k = min(k, m);
if k <= 0
    return;
end

% Weighted sampling without replacement via exponential keys
U = rand(m,1);
keys = -log(U) ./ wts;
[~, ord] = mink(keys, k);
pick = idx(ord);

W(pick) = 1;
end

% ======================================================================
function Ww = assignWeights(A, dist, params, weightModel, alpha, P)
[i,j] = find(A);
m = numel(i);

dist = lower(string(dist));
weightModel = lower(string(weightModel));

% ---- base weights ----
switch dist
    case "uniform"
        if numel(params) < 2
            error('randGPPM:WeightParams','Uniform needs WeightParams=[a b].');
        end
        a = params(1); b = params(2);
        base = a + (b-a)*rand(m,1);

    case "exponential"
        lambda = params(1);
        if lambda <= 0
            error('randGPPM:WeightParams','Exponential needs lambda>0.');
        end
        base = exprnd(1/lambda, m, 1);

    case "lognormal"
        if numel(params) < 2
            error('randGPPM:WeightParams','Lognormal needs [mu sigma].');
        end
        mu = params(1); sigma = params(2);
        if sigma < 0
            error('randGPPM:WeightParams','Lognormal needs sigma>=0.');
        end
        base = lognrnd(mu, sigma, m, 1);

    case "constant"
        base = params(1)*ones(m,1);

    otherwise
        error('randGPPM:WeightDist','Unknown WeightDist "%s".', dist);
end

% ---- couple to propensity ----
if weightModel == "independent"
    w = base;
else
    pij = P(sub2ind(size(P), i, j));
    pij = max(pij, 0);
    f = pij .^ alpha;

    if weightModel == "scaledbyp"
        w = base .* f;

    elseif weightModel == "meanproportionaltop"
        w = base .* f;
        mw = mean(w); mb = mean(base);
        if mw > 0 && mb > 0
            w = w * (mb / mw);
        end
    else
        error('randGPPM:WeightModel','Unknown WeightModel "%s".', weightModel);
    end
end

Ww = zeros(size(A));
Ww(sub2ind(size(A), i, j)) = w;
end

% ======================================================================
function meta = makeMeta(N,B,L,T,opt,W,dagOrder)
meta = struct();
meta.N = N;
meta.B = B;
meta.L = L;
meta.T = T;

meta.weighted = opt.Weighted;
meta.weightDist = string(opt.WeightDist);
meta.weightModel = string(opt.WeightModel);
meta.weightPExponent = opt.WeightPExponent;

meta.allowSelfLoops = opt.AllowSelfLoops;
meta.forceWeakConn  = opt.ForceWeakConn;
meta.forceDAG       = opt.ForceDAG;

meta.seed = opt.Seed;
meta.nnz  = nnz(W);

if opt.ForceDAG
    meta.dagOrder = dagOrder(:);
end
end
