function [edgeList, info] = adjacencyToEdgelist(A, varargin)
%ADJACENCYTOEDGELIST Convert an adjacency matrix to an edge list.
%
%   edgeList = adjacencyToEdgelist(A)
%   edgeList = adjacencyToEdgelist(A, 'Name', Value, ...)
%   [edgeList, info] = adjacencyToEdgelist(A, ...)
%
% INPUT
%   A : n x n adjacency matrix (numeric, sparse or full). Directed by default.
%       For undirected graphs, you may choose how to emit edges via options.
%
% NAME-VALUE OPTIONS
%   'Weighted'     : true/false (default: auto; true if any non-binary weight)
%                   If false, returns Nx2 edge list (i,j) ignoring weights.
%                   If true, returns Nx3 edge list (i,j,w).
%
%   'IncludeZero'  : true/false (default: false). If true, includes zero
%                   entries (dense only) — usually not desired.
%
%   'Directed'     : true/false (default: true). If false, treats A as
%                   undirected and controls duplicates with 'UndirMode'.
%
%   'UndirMode'    : how to emit edges when Directed==false (default: 'upper')
%                   'upper' : emit i<j edges from upper triangle (plus loops)
%                   'lower' : emit i>j edges from lower triangle (plus loops)
%                   'both'  : emit both directions for each undirected edge
%
%   'IncludeSelf'  : true/false (default: true). Include self-loops (i==j)
%                   if present (nonzero) (and/or if IncludeZero==true).
%
%   'NodeIds'      : [] (default) or vector/cellstr/string of length n.
%                   If provided, returns node IDs instead of numeric indices.
%                   - For numeric NodeIds: edgeList columns remain numeric.
%                   - For string/cellstr NodeIds: edgeList is string array
%                     (Nx2) or (Nx3) with weights in a separate numeric vector
%                     is NOT supported (so we instead return a cell array).
%                     See Notes below.
%
%   'SortEdges'    : 'none' (default) | 'lex' | 'stable'
%                   - 'lex' sorts by (src,dst) ascending.
%                   - 'stable' keeps MATLAB find/triu order (already stable).
%
% OUTPUT
%   edgeList : Nx2 or Nx3 edge list.
%              By default uses 1..n node indices.
%
%              If NodeIds are provided:
%                - numeric NodeIds => numeric edgeList, Nx2 or Nx3
%                - string/cellstr NodeIds => edgeList is cell array:
%                     {srcId, dstId, weight} rows if Weighted==true
%                     {srcId, dstId} rows if Weighted==false
%
%   info     : struct with diagnostics:
%               .nNodes, .nEdges, .weighted, .directed, .hadSelfLoops
%
% NOTES
% - For non-numeric NodeIds, mixing labels with numeric weights inside a
%   single MATLAB string array is awkward; we return a cell array to keep
%   weights numeric. This mirrors common MATLAB conventions.
%
% EXAMPLES
%   A = sparse([1 2 3],[2 3 1],[5 1 2],3,3);
%   E = adjacencyToEdgelist(A);                % auto weighted, returns Nx3
%   E = adjacencyToEdgelist(A,'Weighted',false); % returns Nx2
%   E = adjacencyToEdgelist(A,'Directed',false,'UndirMode','upper');
%
% Author:
%   Bazil Sansom (Warwick Business School, University of Warwick)
% -------------------------------------------------------------------------

% -------- Parse options --------
p = inputParser;
p.FunctionName = 'adjacencyToEdgelist';

addParameter(p, 'Weighted', [], @(x) isempty(x) || islogical(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'IncludeZero', false, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'Directed', true, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'UndirMode', 'upper', @(s) ischar(s) || isstring(s));
addParameter(p, 'IncludeSelf', true, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'NodeIds', [], @(x) isempty(x) || isvector(x));
addParameter(p, 'SortEdges', 'none', @(s) ischar(s) || isstring(s));

parse(p, varargin{:});
opts = p.Results;

undirMode = lower(string(opts.UndirMode));
sortMode  = lower(string(opts.SortEdges));

% -------- Validate A --------
if isempty(A)
    error('adjacencyToEdgelist:EmptyInput', 'A is empty.');
end
if ~ismatrix(A) || size(A,1) ~= size(A,2)
    error('adjacencyToEdgelist:BadShape', 'A must be a square matrix.');
end
if ~isnumeric(A)
    error('adjacencyToEdgelist:BadType', 'A must be numeric.');
end

n = size(A,1);
directed = logical(opts.Directed);

% -------- Determine weighting behaviour --------
if isempty(opts.Weighted)
    % Auto: treat as weighted if any entry is not exactly 0/1
    if issparse(A)
        v = nonzeros(A);
    else
        v = A(:);
    end
    v = v(isfinite(v));
    if isempty(v)
        weighted = false;
    else
        weighted = any(v ~= 0 & v ~= 1);
    end
else
    weighted = logical(opts.Weighted);
end

includeZero = logical(opts.IncludeZero);
includeSelf = logical(opts.IncludeSelf);

if includeZero && issparse(A)
    error('adjacencyToEdgelist:IncludeZeroSparse', ...
        'IncludeZero=true is not supported for sparse matrices (would be huge).');
end

% -------- Extract edges --------
if directed
    if includeZero
        % Dense only: include all entries (optionally drop diagonal)
        [ii, jj] = ndgrid(1:n, 1:n);
        ii = ii(:); jj = jj(:);
        if ~includeSelf
            keep = ii ~= jj;
            ii = ii(keep); jj = jj(keep);
        end
        w = A(sub2ind([n n], ii, jj));
        if ~weighted
            edgeList = [ii, jj];
        else
            edgeList = [ii, jj, w(:)];
        end

    else
        % Standard: include nonzeros only
        [ii, jj, w] = find(A);
        if ~includeSelf
            keep = ii ~= jj;
            ii = ii(keep); jj = jj(keep); w = w(keep);
        end
        if ~weighted
            edgeList = [ii, jj];
        else
            edgeList = [ii, jj, w];
        end
    end

else
    % Undirected treatment
    if ~ismember(undirMode, ["upper","lower","both"])
        error('adjacencyToEdgelist:BadUndirMode', ...
            'UndirMode must be ''upper'', ''lower'', or ''both''.');
    end

    if includeZero
        % Dense only
        switch undirMode
            case "upper"
                M = triu(true(n), 0); % incl diagonal
                if ~includeSelf, M = triu(true(n), 1); end
                [ii, jj] = find(M);
                w = A(sub2ind([n n], ii, jj));
            case "lower"
                M = tril(true(n), 0);
                if ~includeSelf, M = tril(true(n), -1); end
                [ii, jj] = find(M);
                w = A(sub2ind([n n], ii, jj));
            case "both"
                % Emit both directions for all off-diagonals; diagonal once
                [iiU, jjU] = find(triu(true(n), 1));
                ii = [iiU; jjU];
                jj = [jjU; iiU];
                w  = [A(sub2ind([n n], iiU, jjU)); A(sub2ind([n n], jjU, iiU))];
                if includeSelf
                    d = (1:n)';
                    ii = [ii; d];
                    jj = [jj; d];
                    w  = [w; diag(A)];
                end
        end

        if ~weighted
            edgeList = [ii, jj];
        else
            edgeList = [ii, jj, w(:)];
        end

    else
        % Nonzeros only
        switch undirMode
            case "upper"
                Au = triu(A, 0);
                if ~includeSelf, Au = triu(A, 1); end
                [ii, jj, w] = find(Au);

            case "lower"
                Al = tril(A, 0);
                if ~includeSelf, Al = tril(A, -1); end
                [ii, jj, w] = find(Al);

            case "both"
                % Emit both directions for each undirected nonzero edge.
                % Approach: take upper-tri nonzeros (off-diagonal) + diagonal.
                Aoff = triu(A, 1);
                [iu, ju, wu] = find(Aoff);

                ii = [iu; ju];
                jj = [ju; iu];
                w  = [wu; wu];   % mirror weight from upper triangle

                if includeSelf
                    [id, jd, wd] = find(spdiags(diag(A), 0, n, n));
                    ii = [ii; id];
                    jj = [jj; jd];
                    w  = [w; wd];
                end
        end

        if ~weighted
            edgeList = [ii, jj];
        else
            edgeList = [ii, jj, w];
        end
    end
end

% -------- Sort edges if requested --------
if ~isempty(edgeList)
    switch sortMode
        case "none"
            % do nothing
        case "stable"
            % find/triu are already stable; keep as-is
        case "lex"
            if size(edgeList,2) >= 2
                [~, ord] = sortrows(edgeList(:,1:2), [1 2]);
                edgeList = edgeList(ord,:);
            end
        otherwise
            error('adjacencyToEdgelist:BadSortEdges', ...
                'SortEdges must be ''none'', ''stable'', or ''lex''.');
    end
end

% -------- Apply NodeIds mapping if provided --------
nodeIds = opts.NodeIds;
if ~isempty(nodeIds)
    if numel(nodeIds) ~= n
        error('adjacencyToEdgelist:BadNodeIds', ...
            'NodeIds must have length n = size(A,1).');
    end

    if isempty(edgeList)
        % Keep empty type consistent
        if isnumeric(nodeIds)
            edgeList = zeros(0, weighted*3 + (~weighted)*2);
        else
            edgeList = cell(0, weighted*3 + (~weighted)*2);
        end
    else
        ii = edgeList(:,1);
        jj = edgeList(:,2);

        if isnumeric(nodeIds)
            src = nodeIds(ii);
            dst = nodeIds(jj);
            if ~weighted
                edgeList = [double(src(:)), double(dst(:))];
            else
                edgeList = [double(src(:)), double(dst(:)), double(edgeList(:,3))];
            end
        else
            % Use cell array output to keep weights numeric and IDs flexible
            src = nodeIds(ii);
            dst = nodeIds(jj);

            if isstring(src), src = cellstr(src); end
            if isstring(dst), dst = cellstr(dst); end
            if iscategorical(src), src = cellstr(string(src)); end
            if iscategorical(dst), dst = cellstr(string(dst)); end

            if ~weighted
                edgeList = [src(:), dst(:)];
            else
                wcol = num2cell(double(edgeList(:,3)));
                edgeList = [src(:), dst(:), wcol];
            end
        end
    end
end

% -------- Info --------
hadSelfLoops = false;
if ~isempty(edgeList)
    hadSelfLoops = any(edgeList(:,1) == edgeList(:,2));
    if ~isnumeric(edgeList) && iscell(edgeList)
        % If NodeIds non-numeric, first two columns are labels; self-loop test
        % is not meaningful without mapping back to indices, so leave false.
        hadSelfLoops = false;
    end
end

info = struct();
info.nNodes       = n;
info.nEdges       = size(edgeList,1);
info.weighted     = weighted;
info.directed     = directed;
info.hadSelfLoops = hadSelfLoops;

end
