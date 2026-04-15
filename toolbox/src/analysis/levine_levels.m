function [s, info] = levine_levels(W, varargin)
%LEVINE_LEVELS Levine trophic levels for a directed network.
%
%   s = LEVINE_LEVELS(W) returns the basal-node-based trophic levels
%   introduced by Levine (1980). For non-basal nodes, levels satisfy
%
%       s_i = b + (1 / k_i^{in}) * sum_j W_{ji} s_j,
%
%   where b is the basal level convention (default b = 1), and basal nodes
%   (nodes with zero in-weight) are assigned level b directly.
%
%   This is equivalent to solving
%
%       (diag(v) - W') s = b v,
%
%   where v_i = max(k_i^{in}, 1).
%
%   [s, info] = LEVINE_LEVELS(W) also returns a struct with diagnostic
%   information.
%
%   Name-value options:
%       'BasalLevel'   - Scalar basal-node level convention (default 1).
%       'WarnSymmetric' - If true, warn when W is symmetric
%                         (default true).
%
%   Inputs:
%       W  - N-by-N adjacency matrix. May be weighted or unweighted,
%            dense or sparse. Must be real, finite, and nonnegative.
%            Self-loops are allowed, though interpretation should be made
%            with care.
%
%   Outputs:
%       s     - N-by-1 vector of Levine trophic levels.
%       info  - Struct with fields:
%                   .kIn
%                   .basal
%                   .nBasal
%                   .components
%                   .nComponents
%                   .basalLevel
%
% Notes:
% - Levine levels require at least one basal node in each weakly
%   connected component.
% - For the improved MacKay--Johnson--Sansom levels, which remain defined
%   without basal nodes, use TROPHIC_LEVELS.
%
%   References:
%   Levine, S. H. (1980). Several measures of trophic structure applicable
%   to complex food webs. Journal of Theoretical Biology, 83, 195--207.
%
%   MacKay, R. S., Johnson, S., & Sansom, B. (2020). How directed is a
%   directed network? Royal Society Open Science, 7, 201138.
%
%   See also: trophic_levels, incoherence

    % ------------------------------------------------------------
    % Parse inputs
    % ------------------------------------------------------------
    p = inputParser;
    p.FunctionName = mfilename;

    addRequired(p, 'W', @(x) validateAdjacency_(x));
    addParameter(p, 'BasalLevel', 1, @(x) validateattributes(x, ...
        {'numeric'}, {'real','finite','scalar'}));
    addParameter(p, 'WarnSymmetric', true, @(x) islogical(x) || isnumeric(x));

    parse(p, W, varargin{:});
    opts = p.Results;

    W = double(W);
    n = size(W, 1);

    % ------------------------------------------------------------
    % Optional warning for symmetric network
    % ------------------------------------------------------------
    if opts.WarnSymmetric && isequal(W, W.')
        warning('levine_levels:SymmetricAdjacency', ...
            ['Adjacency matrix is symmetric. ', ...
             'Network may be undirected, in which case trophic levels ', ...
             'may be of limited interpretive value.']);
    end

    % ------------------------------------------------------------
    % Basic quantities
    % ------------------------------------------------------------
    kIn = full(sum(W, 1)).';
    basal = (kIn == 0);

    % Check structural condition for existence:
    % every weakly connected component must contain at least one basal node.
    Aweak = spones(W + W.');
    comp = weakComponents_(Aweak);
    nComp = max(comp);

    badComp = false(nComp, 1);
    for c = 1:nComp
        idx = (comp == c);
        badComp(c) = ~any(basal(idx));
    end

    if any(badComp)
        badList = find(badComp);
        error('levine_levels:NoBasalInComponent', ...
            ['Levine trophic levels are undefined because at least one ', ...
             'weakly connected component has no basal node. ', ...
             'Components without basal nodes: %s. ', ...
             'Use trophic_levels for the improved MacKay2020 levels instead.'], ...
             mat2str(badList.'));
    end

    % ------------------------------------------------------------
    % Solve linear system
    % ------------------------------------------------------------
    v = max(kIn, 1);
    rhs = opts.BasalLevel * v;

    if issparse(W)
        L = spdiags(v, 0, n, n) - W.';
    else
        L = diag(v) - W.';
    end

    s = L \ rhs;
    s = full(s(:));

    % ------------------------------------------------------------
    % Diagnostics
    % ------------------------------------------------------------
    if nargout > 1
        info = struct();
        info.kIn        = kIn;
        info.basal      = basal;
        info.nBasal     = nnz(basal);
        info.components = comp;
        info.nComponents = nComp;
        info.basalLevel = opts.BasalLevel;
    end
end

% =====================================================================
function tf = validateAdjacency_(W)
    tf = true;

    if ~(isnumeric(W) || islogical(W))
        error('levine_levels:InvalidAdjacencyType', ...
            'W must be numeric or logical.');
    end

    if ndims(W) ~= 2 || size(W,1) ~= size(W,2)
        error('levine_levels:AdjacencyNotSquare', ...
            'W must be a square adjacency matrix.');
    end

    if ~isreal(W)
        error('levine_levels:AdjacencyComplex', ...
            'W must be real.');
    end

    if any(~isfinite(W(:)))
        error('levine_levels:AdjacencyNonFinite', ...
            'W must contain only finite values.');
    end

    if any(W(:) < 0)
        error('levine_levels:AdjacencyNegative', ...
            'W must be nonnegative.');
    end
end

% =====================================================================
function comp = weakComponents_(A)
%WEAKCOMPONENTS_ Weakly connected components of an undirected support graph.
%
% Input:
%   A  - N-by-N logical/sparse adjacency for an undirected graph.
%
% Output:
%   comp - component labels in 1,2,...,K

    n = size(A,1);
    comp = zeros(n,1);
    cid = 0;

    for i = 1:n
        if comp(i) ~= 0
            continue;
        end

        cid = cid + 1;
        stack = i;
        comp(i) = cid;

        while ~isempty(stack)
            v = stack(end);
            stack(end) = [];

            nbr = find(A(v,:));
            nbr = nbr(comp(nbr) == 0);

            if ~isempty(nbr)
                comp(nbr) = cid;
                stack = [stack, nbr(:).']; %#ok<AGROW>
            end
        end
    end
end