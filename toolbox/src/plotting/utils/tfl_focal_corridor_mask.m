function corr = tfl_focal_corridor_mask(W, h, focalIdx, opts)
%TFL_FOCAL_CORRIDOR_MASK
% Build a multi-hop trophically aligned corridor around a focal node.
%
% Edge i->j is considered aligned if h(j) > h(i) + MinDeltaH.
%
% Inputs:
%   W        : NxN weighted adjacency
%   h        : Nx1 trophic levels
%   focalIdx : focal node index
%   opts fields:
%       .MaxHops   (default 2)
%       .MinDeltaH (default 1e-6)
%       .MinWeight (default 0)
%
% Outputs:
%   corr.nodeMask
%   corr.edgeMask
%   corr.upMask
%   corr.downMask
%   corr.Aalign
%   corr.distUp
%   corr.distDown

    if nargin < 4 || isempty(opts), opts = struct(); end
    if ~isfield(opts,'MaxHops'),   opts.MaxHops = 2; end
    if ~isfield(opts,'MinDeltaH'), opts.MinDeltaH = 1e-6; end
    if ~isfield(opts,'MinWeight'), opts.MinWeight = 0; end

    W = double(W);
    h = h(:);
    n = size(W,1);

    [ii,jj,vv] = find(W);
    isAligned = (h(jj) > h(ii) + opts.MinDeltaH) & (vv > opts.MinWeight);

    ia = ii(isAligned);
    ja = jj(isAligned);

    Aalign = sparse(ia, ja, true, n, n);

    % Downstream: reachable from focal in aligned graph
    distDown = bfs_reachable_(Aalign, focalIdx, opts.MaxHops);

    % Upstream: nodes that can reach focal in aligned graph
    distUp = bfs_reachable_(Aalign', focalIdx, opts.MaxHops);

    downMask = isfinite(distDown);
    upMask   = isfinite(distUp);

    nodeMask = upMask | downMask;
    nodeMask(focalIdx) = true;

    % Highlight aligned edges internal to upstream/downstream corridor sets
    
    edgeMask = false(size(W));

    for k = 1:numel(ia)
        u = ia(k);
        v = ja(k);

        % Keep any trophically aligned edge internal to the downstream corridor
        inDown = downMask(u) && downMask(v);

        % Keep any trophically aligned edge internal to the upstream corridor
        inUp   = upMask(u) && upMask(v);

        if inDown || inUp
            edgeMask(u,v) = true;
        end
    end
    
%{
    for k = 1:numel(ia)
        u = ia(k);
        v = ja(k);

        inDown = downMask(u) && downMask(v) && (distDown(v) > distDown(u));
        inUp   = upMask(u)   && upMask(v)   && (distUp(u)   > distUp(v));

        if inDown || inUp
            edgeMask(u,v) = true;
        end
    end
%}
    corr = struct();
    corr.nodeMask = nodeMask;
    corr.edgeMask = edgeMask;
    corr.upMask = upMask;
    corr.downMask = downMask;
    corr.Aalign = Aalign;
    corr.distUp = distUp;
    corr.distDown = distDown;
end

function dist = bfs_reachable_(A, source, maxHops)
%BFS_REACHABLE_  Hop distance from source in sparse adjacency A.
%
% Unreachable nodes get Inf.

    n = size(A,1);
    dist = inf(n,1);
    dist(source) = 0;

    frontier = false(n,1);
    frontier(source) = true;

    for d = 1:maxHops
        nbrs = (A' * frontier) > 0;
        nbrs = nbrs & ~isfinite(dist);

        if ~any(nbrs)
            break;
        end

        dist(nbrs) = d;
        frontier = nbrs;
    end
end