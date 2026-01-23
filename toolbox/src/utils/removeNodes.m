function Wsub = removeNodes(W, nodesToRemove)
%REMOVENODES  Remove specified nodes from a directed graph.
%
%   Wsub = REMOVENODES(W, nodesToRemove)
%
%   INPUTS
%     W              : NxN adjacency matrix (directed, weighted or unweighted)
%     nodesToRemove  : vector of node indices to remove
%
%   OUTPUT
%     Wsub           : (N-k)x(N-k) adjacency matrix with the specified nodes removed
%
%   Example:
%     Wsub = removeNodes(W, [3 7 12]);

    % ---- Input checks ----
    if ~ismatrix(W) || size(W,1) ~= size(W,2)
        error('removeNodes:WNotSquare', 'W must be a square adjacency matrix.');
    end

    N = size(W,1);

    nodesToRemove = unique(nodesToRemove(:));

    if any(nodesToRemove < 1 | nodesToRemove > N)
        error('removeNodes:BadNodeIndex', ...
              'nodesToRemove contains invalid node indices.');
    end

    % ---- Nodes to keep ----
    keep = true(N,1);
    keep(nodesToRemove) = false;

    % ---- Subgraph ----
    Wsub = W(keep, keep);
end
