function [A, nodeIds, map, info] = edgelistToAdjacency(edgeList, varargin)
%EDGELIST2ADJ Convert an edge list to an adjacency matrix.
%
%   A = edgelistToAdjacency(edgeList)
%   [A, nodeIds] = edgelistToAdjacency(edgeList)
%   [A, nodeIds, map, info] = edgelistToAdjacency(edgeList, Name, Value, ...)
%
% INPUT
%   edgeList : Nx2 or Nx3 edge list.
%              - Columns 1-2: source, target.
%              - Column 3 (optional): weight.
%              Node IDs may be numeric, string, char, categorical, or cellstr.
%
% NAME-VALUE OPTIONS
%   'Weighted'       : true/false (default: auto; true if 3 cols, else false)
%   'DefaultWeight'  : scalar default weight if Nx2 and Weighted==true (default: 1)
%   'Combine'        : how to handle duplicate edges (default: 'sum')
%                      'sum' | 'last' | 'first' | 'max' | 'min'
%   'Directed'       : true/false (default: true). If false, symmetrises.
%   'SymMethod'      : symmetrisation method if Directed==false (default: 'sum')
%                      'sum' | 'max' | 'min' | 'mean'
%   'Sparse'         : true/false (default: true). Use sparse output.
%   'SortNodes'      : 'stable' | 'sorted' (default: 'stable')
%                      Determines the node ordering for non-numeric IDs.
%   'Validate'       : true/false (default: true). Input checks.
%
% OUTPUTS
%   A       : adjacency matrix (sparse by default).
%   nodeIds : node ID list in the order used for A (numeric or string).
%   map     : containers.Map from nodeId -> index (only for non-numeric IDs).
%             Empty if numeric IDs are already 1..n and no remapping needed.
%   info    : struct with basic diagnostics:
%               .nNodes, .nEdgesInput, .nEdgesUnique, .hadWeights,
%               .hadDuplicates, .wasRemapped, .directed
%
% NOTES
% - If edgeList uses arbitrary numeric IDs (e.g., 1001, 2005), this function
%   remaps them to 1..n by default (public-toolbox safe behaviour).
% - If you *know* your numeric IDs are already 1..n, you'll still be fine.
%
% EXAMPLES
%   % Unweighted
%   A = edgelist2adj([1 2; 2 3; 1 3]);
%
%   % Weighted with duplicates summed
%   E = [1 2 5; 1 2 7; 2 1 3];
%   A = edgelist2adj(E, 'Combine','sum');
%
%   % String node IDs
%   E = ["A","B"; "B","C"; "A","C"];
%   [A,nodeIds] = edgelist2adj(E);
%
%   % Undirected
%   A = edgelist2adj([1 2; 2 3], 'Directed', false);
%

% -------- Parse options --------
p = inputParser;
p.FunctionName = 'edgelist2adj';

addParameter(p, 'Weighted', [], @(x) isempty(x) || islogical(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'DefaultWeight', 1, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'Combine', 'sum', @(s) ischar(s) || isstring(s));
addParameter(p, 'Directed', true, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'SymMethod', 'sum', @(s) ischar(s) || isstring(s));
addParameter(p, 'Sparse', true, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'SortNodes', 'stable', @(s) ischar(s) || isstring(s));
addParameter(p, 'Validate', true, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));

parse(p, varargin{:});
opts = p.Results;

combine = lower(string(opts.Combine));
symMethod = lower(string(opts.SymMethod));
sortMode = lower(string(opts.SortNodes));

% -------- Basic validation --------
if opts.Validate
    if isempty(edgeList)
        error('edgelist2adj:EmptyInput', 'edgeList is empty.');
    end
    if ~(ismatrix(edgeList) && (size(edgeList,2) == 2 || size(edgeList,2) == 3))
        error('edgelist2adj:BadShape', 'edgeList must be Nx2 or Nx3.');
    end
end

nIn = size(edgeList,1);
hadWeights = (size(edgeList,2) == 3);

if isempty(opts.Weighted)
    weighted = hadWeights;
else
    weighted = logical(opts.Weighted);
end

% -------- Extract columns --------
src = edgeList(:,1);
dst = edgeList(:,2);

if weighted
    if hadWeights
        w = edgeList(:,3);
    else
        w = repmat(opts.DefaultWeight, nIn, 1);
    end
else
    w = ones(nIn,1);
end

% -------- Normalise types for node IDs --------
% We support numeric IDs or "labels" (strings/categoricals/cellstr/chars).
isNumericIds = isnumeric(src) && isnumeric(dst);

map = [];
wasRemapped = false;

if isNumericIds
    srcNum = double(src);
    dstNum = double(dst);

    if opts.Validate
        if any(~isfinite(srcNum)) || any(~isfinite(dstNum))
            error('edgelist2adj:NonFiniteIds', 'Numeric node IDs must be finite.');
        end
        if any(srcNum < 1 | dstNum < 1 | mod(srcNum,1)~=0 | mod(dstNum,1)~=0)
            error('edgelist2adj:BadNumericIds', 'Numeric node IDs must be positive integers.');
        end
    end

    % Remap arbitrary integer labels to 1..n to avoid huge matrices.
    allIds = [srcNum; dstNum];
    if ~isequal(sort(allIds), (1:max(allIds))') % quick-ish heuristic, still safe
        [nodeIds, ~, idx] = unique(allIds, 'stable');
        srcIdx = idx(1:nIn);
        dstIdx = idx(nIn+1:end);
        wasRemapped = true;
    else
        nodeIds = (1:max(allIds))';
        srcIdx = srcNum;
        dstIdx = dstNum;
    end

else
    % Convert everything to string labels
    srcStr = normalizeToString(src);
    dstStr = normalizeToString(dst);
    allStr = [srcStr; dstStr];

    if sortMode == "sorted"
        nodeIds = unique(allStr);
    else
        nodeIds = unique(allStr, 'stable');
    end

    % Map labels -> indices
    map = containers.Map(cellstr(nodeIds), num2cell(1:numel(nodeIds)));
    srcIdx = cellfun(@(s) map(s), cellstr(srcStr));
    dstIdx = cellfun(@(s) map(s), cellstr(dstStr));
    wasRemapped = true;
end

nNodes = numel(nodeIds);

% -------- Combine duplicates --------
% Build sparse first (even if final dense requested) for combination control.
switch combine
    case "sum"
        A = sparse(srcIdx, dstIdx, double(w), nNodes, nNodes);

    case {"last","first","max","min"}
        % Use grouping on (i,j) then apply reducer
        ij = [srcIdx(:), dstIdx(:)];
        [ijU, ~, g] = unique(ij, 'rows', 'stable');
        w = double(w(:));

        switch combine
            case "last"
                ww = accumarray(g, (1:numel(w))', [], @max); % last index per group
                v = w(ww);
            case "first"
                ww = accumarray(g, (1:numel(w))', [], @min); % first index per group
                v = w(ww);
            case "max"
                v = accumarray(g, w, [], @max);
            case "min"
                v = accumarray(g, w, [], @min);
        end

        A = sparse(ijU(:,1), ijU(:,2), v, nNodes, nNodes);

    otherwise
        error('edgelist2adj:BadCombine', 'Unknown Combine option: %s', combine);
end

hadDuplicates = nnz(A) < nIn; % rough but informative

% -------- Symmetrise if undirected requested --------
directed = logical(opts.Directed);
if ~directed
    switch symMethod
        case "sum"
            A = A + A.';
        case "max"
            A = max(A, A.');
        case "min"
            A = min(A, A.');
        case "mean"
            A = (A + A.')/2;
        otherwise
            error('edgelist2adj:BadSymMethod', 'Unknown SymMethod option: %s', symMethod);
    end
end

% -------- Output dense if requested --------
if ~logical(opts.Sparse)
    A = full(A);
end

% -------- Info --------
info = struct();
info.nNodes       = nNodes;
info.nEdgesInput  = nIn;
info.nEdgesUnique = nnz(A);
info.hadWeights   = hadWeights;
info.hadDuplicates = hadDuplicates;
info.wasRemapped  = wasRemapped;
info.directed     = directed;

end

% ---- helper: robust conversion to string ----
function s = normalizeToString(x)
    if isstring(x)
        s = x;
    elseif iscategorical(x)
        s = string(x);
    elseif iscell(x)
        % cellstr or mixed: convert each element
        s = string(x);
    elseif ischar(x)
        s = string(cellstr(x)); % handle char matrix
    else
        error('edgelist2adj:BadIdType', ...
            'Unsupported node ID type: %s. Use numeric, string, char, categorical, or cellstr.', class(x));
    end
    s = s(:);
end
