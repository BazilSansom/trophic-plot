function adj = randDG(n, p, connected, varargin)
%RANDDG  Random directed Erdos-Renyi digraph (optionally weighted).
%
%   adj = randDG(n, p)
%   adj = randDG(n, p, connected)
%   adj = randDG(..., 'Name', value, ...)
%
% INPUTS
%   n         : number of nodes (positive integer)
%   p         : edge probability in [0,1]
%   connected : if equal to 'connected' (string/char), force weak connectivity
%
% Name–value options
%   'AllowSelfLoops'  (false)    : allow i->i edges
%   'ForceWeakConn'   (auto)     : true/false; if omitted uses strcmp(connected,'connected')
%   'MaxTries'        (5000)     : rejection-sampling cap for connectivity
%   'Seed'            ([] )      : set RNG seed (scalar) for reproducibility
%
%   'Weighted'        (false)    : assign random weights to present edges
%   'WeightDist'      ('uniform'): 'uniform' | 'exponential' | 'lognormal' | 'constant'
%   'WeightParams'    ([0 1])    : parameters for WeightDist (see below)
%
% WeightDist / WeightParams
%   'uniform'      : [a b]          -> U(a,b)
%   'exponential'  : lambda         -> Exp(lambda)   (mean = 1/lambda)
%   'lognormal'    : [mu sigma]     -> LogNormal(mu,sigma)
%   'constant'     : c              -> constant weight c
%
% OUTPUT
%   adj : n-by-n adjacency matrix
%         - logical if unweighted
%         - double  if weighted
%
% NOTES
%   - Generates a directed ER graph: each ordered pair (i,j), i~=j, included with prob p.
%   - "Weakly connected" means connected after symmetrising (adj|adj').
%

% ---- parse inputs ----
if nargin < 3 || isempty(connected)
    connected = '';
end

ip = inputParser;
ip.addParameter('AllowSelfLoops', false, @(x)islogical(x) && isscalar(x));
ip.addParameter('ForceWeakConn', [], @(x)islogical(x) && isscalar(x));
ip.addParameter('MaxTries', 5000, @(x)isnumeric(x) && isscalar(x) && x>=1);
ip.addParameter('Seed', [], @(x)isempty(x) || (isnumeric(x) && isscalar(x)));

ip.addParameter('Weighted', false, @(x)islogical(x) && isscalar(x));
ip.addParameter('WeightDist', 'uniform', @(x)ischar(x) || isstring(x));
ip.addParameter('WeightParams', [0 1], @(x)isnumeric(x) && ~isempty(x));

ip.parse(varargin{:});
opt = ip.Results;

if isempty(opt.ForceWeakConn)
    opt.ForceWeakConn = (ischar(connected) || isstring(connected)) && strcmpi(string(connected), "connected");
end

if ~isempty(opt.Seed)
    rng(opt.Seed);
end

validateattributes(n, {'numeric'}, {'scalar','integer','>=',1});
validateattributes(p, {'numeric'}, {'scalar','>=',0,'<=',1});

% ---- helper mask for eligible edges ----
mask = true(n);
if ~opt.AllowSelfLoops
    mask(1:n+1:end) = false;
end

% ---- generate (with optional weak-connectivity rejection) ----
if ~opt.ForceWeakConn
    adj = (rand(n) <= p) & mask;
    adj = logical(adj);

    if opt.Weighted
        adj = assignWeights(adj, opt.WeightDist, opt.WeightParams);
    end
    return;
end

for t = 1:opt.MaxTries
    adj = (rand(n) <= p) & mask;
    adj = logical(adj);

    if isWeaklyConnected(adj)
        if opt.Weighted
            adj = assignWeights(adj, opt.WeightDist, opt.WeightParams);
        end
        return;
    end
end

error('randDG:MaxTriesExceeded', ...
    'Failed to generate a weakly connected digraph after %d tries. Try increasing p or MaxTries.', opt.MaxTries);

end

% -------------------------------------------------------------------------
function tf = isWeaklyConnected(adj)
% Weak connectivity of directed graph = connectivity of undirected symmetrisation.
A = adj | adj.';  % undirected adjacency (logical)

% quick reject: isolates
if any(sum(A,2) == 0)
    tf = false;
    return;
end

% Use graph/conncomp if available (faster/cleaner), else fallback BFS.
if exist('graph','file') == 2 && exist('conncomp','file') == 2
    G = graph(A);
    cc = conncomp(G);
    tf = (max(cc) == 1);
else
    % fallback: simple BFS/DFS from node 1
    n = size(A,1);
    seen = false(n,1);
    stack = 1;
    seen(1) = true;

    while ~isempty(stack)
        v = stack(end); stack(end) = [];
        nbr = find(A(v,:)).';
        new = nbr(~seen(nbr));
        seen(new) = true;
        stack = [stack; new]; %#ok<AGROW>
    end
    tf = all(seen);
end
end

% -------------------------------------------------------------------------
function W = assignWeights(A, dist, params)
%ASSIGNWEIGHTS  Assign random positive weights to existing edges.
%
%   W = assignWeights(A, dist, params)
%
%   A : logical adjacency
%   W : weighted adjacency (double)

[i,j] = find(A);
m = numel(i);

dist = lower(string(dist));

switch dist
    case "uniform"
        if numel(params) < 2
            error('randDG:WeightParams', 'Uniform weights require WeightParams = [a b].');
        end
        a = params(1); b = params(2);
        w = a + (b - a) * rand(m,1);

    case "exponential"
        lambda = params(1);
        if lambda <= 0
            error('randDG:WeightParams', 'Exponential weights require lambda > 0.');
        end
        w = exprnd(1/lambda, m, 1); % mean 1/lambda

    case "lognormal"
        if numel(params) < 2
            error('randDG:WeightParams', 'Lognormal weights require WeightParams = [mu sigma].');
        end
        mu = params(1); sigma = params(2);
        if sigma < 0
            error('randDG:WeightParams', 'Lognormal weights require sigma >= 0.');
        end
        w = lognrnd(mu, sigma, m, 1);

    case "constant"
        c = params(1);
        w = c * ones(m,1);

    otherwise
        error('randDG:WeightDist', 'Unknown WeightDist "%s".', dist);
end

W = zeros(size(A));
idx = sub2ind(size(A), i, j);
W(idx) = w;
end
