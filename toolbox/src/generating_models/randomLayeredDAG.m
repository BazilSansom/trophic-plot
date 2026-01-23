function [W, layerOfNode, layerNodes, info] = randomLayeredDAG(N, density, numLayers, varargin)
%RANDOMLAYEREDDAG Generate a random layered DAG (perfectly layered: only k -> k+1).
%
%   [W, layerOfNode, layerNodes, info] = randomLayeredDAG(N, density, numLayers, ...)
%
% Creates an NxN adjacency matrix W for a directed acyclic graph with nodes
% arranged into numLayers layers. Edges are only permitted from layer k to
% layer k+1 (no skips, no intra-layer edges). This yields a DAG with a strict
% layering structure (a "perfectly layered" baseline generator).
%
% INPUTS
%   N          : number of nodes (positive integer)
%   density    : density parameter in [0,1]. Interpretation depends on 'Sampling':
%                - 'bernoulli': each allowable edge included with prob=density
%                - 'fixedPerLayerPair': each layer-pair block gets exactly round(density * possibleEdges)
%                - 'fixedTotal': total edges across all layer-pairs is exactly round(density * totalPossible)
%   numLayers  : number of layers. If empty or omitted, uses max(2, round(sqrt(N))).
%
% NAME-VALUE OPTIONS
%   'Sampling'      : 'bernoulli' (default) | 'fixedPerLayerPair' | 'fixedTotal'
%   'EnsureSpine'   : true/false (default: false)
%                     If true, ensures at least one edge from each layer k to k+1
%                     (when both layers are non-empty). Useful for avoiding
%                     completely disconnected successive layers.
%   'LayerSizes'    : [] (default) or vector of length numLayers summing to N
%                     If provided, uses these layer sizes instead of random allocation.
%   'NodeOrder'     : [] (default) or permutation of 1..N
%                     If provided, uses this node order when assigning nodes to layers.
%   'Weighted'      : true/false (default: false)
%                     If true, assigns weights via 'WeightFun' (or all ones if omitted).
%   'WeightFun'     : function handle w = f(m) returning m weights (default: [])
%                     Only used if Weighted==true. If empty, weights are ones.
%   'Sparse'        : true/false (default: true)
%   'Validate'      : true/false (default: true)
%
% OUTPUTS
%   W           : adjacency matrix (sparse by default)
%   layerOfNode : Nx1 integer layer index for each node (1..numLayers)
%   layerNodes  : cell array of nodes in each layer
%   info        : struct with diagnostics:
%                 .N, .numLayers, .layerSizes
%                 .sampling, .density
%                 .nPossibleEdges, .nEdgesRequested, .nEdgesRealized
%                 .ensureSpine, .hadWeights
%
% NOTES
% - This generator produces your "regime 1" baseline (perfect layering, DAG).
% - For "regime 2" (DAG but varying incoherence), you’ll typically generalise
%   to allow forward skips (k -> k+r) with a controllable skip-length distribution.
%
% EXAMPLE
%   [W, layerOfNode] = randomLayeredDAG(50, 0.2, [], 'Sampling','bernoulli');
%
% See also: util.edgelistToAdjacency (if you add it)

% -------------------- parse inputs --------------------
if nargin < 2
    error('Usage: [W, layerOfNode, layerNodes, info] = randomLayeredDAG(N, density, [numLayers], ...)');
end
if nargin < 3 || isempty(numLayers)
    numLayers = max(2, round(sqrt(N)));
end

p = inputParser;
p.FunctionName = 'randomLayeredDAG';

addParameter(p, 'Sampling', 'bernoulli', @(s) ischar(s) || isstring(s));
addParameter(p, 'EnsureSpine', false, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'LayerSizes', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
addParameter(p, 'NodeOrder', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
addParameter(p, 'Weighted', false, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'WeightFun', [], @(f) isempty(f) || isa(f, 'function_handle'));
addParameter(p, 'Sparse', true, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'Validate', true, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));

parse(p, varargin{:});
opts = p.Results;

sampling = lower(string(opts.Sampling));
ensureSpine = logical(opts.EnsureSpine);
useSparse = logical(opts.Sparse);
doValidate = logical(opts.Validate);
weighted = logical(opts.Weighted);

% -------------------- validation --------------------
if doValidate
    if ~(isscalar(N) && isnumeric(N) && isfinite(N) && N == floor(N) && N > 0)
        error('N must be a positive integer.');
    end
    if ~(isscalar(numLayers) && isnumeric(numLayers) && isfinite(numLayers) && numLayers == floor(numLayers) && numLayers >= 2)
        error('numLayers must be an integer >= 2.');
    end
    if N < numLayers
        error('Number of layers cannot exceed number of nodes.');
    end
    if ~(isscalar(density) && isnumeric(density) && isfinite(density) && density >= 0 && density <= 1)
        error('density must be a scalar in [0,1].');
    end
    if ~ismember(sampling, ["bernoulli","fixedperlayerpair","fixedtotal"])
        error('Sampling must be ''bernoulli'', ''fixedPerLayerPair'', or ''fixedTotal''.');
    end
end

% -------------------- assign nodes to layers --------------------
if ~isempty(opts.LayerSizes)
    layerSizes = opts.LayerSizes(:);
    if doValidate
        if numel(layerSizes) ~= numLayers
            error('LayerSizes must have length numLayers.');
        end
        if any(layerSizes < 1) || any(layerSizes ~= floor(layerSizes))
            error('LayerSizes must be positive integers (>=1).');
        end
        if sum(layerSizes) ~= N
            error('LayerSizes must sum to N.');
        end
    end
else
    % random composition of N into numLayers positive parts
    layerSizes = ones(numLayers,1);
    remaining = N - numLayers;
    if remaining > 0
        idx = randi(numLayers, remaining, 1);
        for k = 1:numLayers
            layerSizes(k) = layerSizes(k) + sum(idx == k);
        end
    end
end

if ~isempty(opts.NodeOrder)
    nodeOrder = opts.NodeOrder(:)';
    if doValidate
        if numel(nodeOrder) ~= N || ~isequal(sort(nodeOrder), 1:N)
            error('NodeOrder must be a permutation of 1..N.');
        end
    end
else
    nodeOrder = randperm(N);
end

layerNodes = cell(numLayers,1);
layerOfNode = zeros(N,1);
startIdx = 1;
for k = 1:numLayers
    nK = layerSizes(k);
    nodesK = nodeOrder(startIdx:startIdx+nK-1);
    layerNodes{k} = nodesK(:);
    layerOfNode(nodesK) = k;
    startIdx = startIdx + nK;
end

% -------------------- compute possible edges --------------------
possiblePerPair = zeros(numLayers-1,1);
for k = 1:numLayers-1
    possiblePerPair(k) = numel(layerNodes{k}) * numel(layerNodes{k+1});
end
nPossibleEdges = sum(possiblePerPair);

% -------------------- decide requested edges --------------------
switch sampling
    case "bernoulli"
        nEdgesRequested = round(density * nPossibleEdges); % expected-ish, for info only
    case "fixedperlayerpair"
        nEdgesRequested = sum(round(density * possiblePerPair));
    case "fixedtotal"
        nEdgesRequested = round(density * nPossibleEdges);
end

% -------------------- sample edges block-by-block --------------------
Icells = cell(numLayers-1,1);
Jcells = cell(numLayers-1,1);
Wcells = cell(numLayers-1,1);

% For fixedTotal, allocate per layer-pair using a multinomial-like draw
mPerPair = zeros(numLayers-1,1);
if sampling == "fixedtotal"
    if nPossibleEdges == 0
        mPerPair(:) = 0;
    else
        % Expected allocation proportional to possible edges, then fix rounding to sum exactly
        mFloat = nEdgesRequested * (possiblePerPair / nPossibleEdges);
        mPerPair = floor(mFloat);
        remainder = nEdgesRequested - sum(mPerPair);
        if remainder > 0
            [~, order] = sort(mFloat - mPerPair, 'descend');
            mPerPair(order(1:remainder)) = mPerPair(order(1:remainder)) + 1;
        end
        % Cap by feasibility
        mPerPair = min(mPerPair, possiblePerPair);
        % If capping reduced total, redistribute leftover where possible
        leftover = nEdgesRequested - sum(mPerPair);
        if leftover > 0
            slack = possiblePerPair - mPerPair;
            while leftover > 0 && any(slack > 0)
                k = find(slack > 0, 1, 'first');
                add = min(leftover, slack(k));
                mPerPair(k) = mPerPair(k) + add;
                slack(k) = slack(k) - add;
                leftover = leftover - add;
            end
        end
    end
end

for k = 1:numLayers-1
    src = layerNodes{k};
    dst = layerNodes{k+1};
    nA = numel(src);
    nB = numel(dst);
    nBlock = nA * nB;

    if nBlock == 0
        Icells{k} = zeros(0,1);
        Jcells{k} = zeros(0,1);
        Wcells{k} = zeros(0,1);
        continue;
    end

    switch sampling
        case "bernoulli"
            % Sample each possible edge independently
            M = rand(nA, nB) < density;
            [ia, ib] = find(M);
            i = src(ia);
            j = dst(ib);

        case "fixedperlayerpair"
            m = round(density * nBlock);
            m = min(m, nBlock);
            if m == 0
                i = zeros(0,1); j = zeros(0,1);
            else
                % Sample m unique edges uniformly within block
                idx = randperm(nBlock, m);
                [ia, ib] = ind2sub([nA, nB], idx);
                i = src(ia(:));
                j = dst(ib(:));
            end

        case "fixedtotal"
            m = mPerPair(k);
            if m == 0
                i = zeros(0,1); j = zeros(0,1);
            else
                idx = randperm(nBlock, m);
                [ia, ib] = ind2sub([nA, nB], idx);
                i = src(ia(:));
                j = dst(ib(:));
            end
    end

    Icells{k} = i(:);
    Jcells{k} = j(:);

    if weighted
        if isempty(opts.WeightFun)
            Wcells{k} = ones(numel(i),1);
        else
            Wcells{k} = opts.WeightFun(numel(i));
            if doValidate
                if ~isnumeric(Wcells{k}) || numel(Wcells{k}) ~= numel(i)
                    error('WeightFun must return a numeric vector of length m.');
                end
            end
        end
    else
        Wcells{k} = ones(numel(i),1);
    end
end

I = vertcat(Icells{:});
J = vertcat(Jcells{:});
V = vertcat(Wcells{:});

% Build sparse adjacency
W = sparse(I, J, V, N, N);

% -------------------- ensure spine (optional) --------------------
if ensureSpine
    for k = 1:numLayers-1
        src = layerNodes{k};
        dst = layerNodes{k+1};
        if isempty(src) || isempty(dst)
            continue;
        end
        i = src(randi(numel(src)));
        j = dst(randi(numel(dst)));

        if weighted && ~isempty(opts.WeightFun)
            v = opts.WeightFun(1);
            if ~isnumeric(v) || numel(v) ~= 1
                error('WeightFun must return a scalar when called with m=1.');
            end
        else
            v = 1;
        end

        W(i,j) = max(W(i,j), v); % ensure at least one edge; if present, keep max weight
    end
end

if ~useSparse
    W = full(W);
end

% -------------------- info --------------------
info = struct();
info.N = N;
info.numLayers = numLayers;
info.layerSizes = layerSizes(:);
info.sampling = char(sampling);
info.density = density;
info.nPossibleEdges = nPossibleEdges;
info.nEdgesRequested = nEdgesRequested;
info.nEdgesRealized = nnz(W);
info.ensureSpine = ensureSpine;
info.hadWeights = weighted;

end
