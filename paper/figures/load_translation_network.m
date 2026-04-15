function net = load_translation_network(dataDir, varargin)
%LOAD_TRANSLATION_NETWORK Load book-translation network for plotting/TFL.
%
% net = load_translation_network(dataDir)
% net = load_translation_network(dataDir, 'TopN', 80, ...
%                                'KeepLargestWeakComponent', true, ...
%                                'DropSelfLoops', true)
%
% Expected files in dataDir:
%   - books_edges.tsv   (required)
%   - books_nodes.tsv   (optional metadata)
%
% The function tries to detect common column names automatically.
%
% OUTPUT
%   net.W            sparse weighted adjacency matrix
%   net.labels       string array of node labels
%   net.G            MATLAB digraph
%   net.edgeTable    table: Source, Target, Weight
%   net.nodeTable    table with Name + weights (+ optional merged metadata)
%   net.inWeight     column vector
%   net.outWeight    column vector
%   net.totalWeight  in + out
%   net.paths        struct with file paths used
%   net.opts         options used
%
% NOTES
% - Direction is Source -> Target.
% - Duplicate edges are aggregated by summing weights.
% - By default, self-loops are dropped.
% - Node filtering is applied after building the full network.

    p = inputParser;
    addRequired(p, 'dataDir', @(x) ischar(x) || isstring(x));
    addParameter(p, 'EdgeFile', 'books_edges.tsv', @(x) ischar(x) || isstring(x));
    addParameter(p, 'NodeFile', 'books_nodes.tsv', @(x) ischar(x) || isstring(x));
    addParameter(p, 'DropSelfLoops', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'MinEdgeWeight', 0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'MinNodeTotalWeight', 0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'TopN', inf, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'KeepLargestWeakComponent', false, @(x) islogical(x) && isscalar(x));
    parse(p, dataDir, varargin{:});
    opts = p.Results;

    dataDir  = char(opts.dataDir);
    edgePath = fullfile(dataDir, char(opts.EdgeFile));
    nodePath = fullfile(dataDir, char(opts.NodeFile));

    if ~isfile(edgePath)
        error('Edge file not found: %s', edgePath);
    end

    % ------------------ read edges ------------------
    E = read_tsv_robust(edgePath);

    srcVar = pick_var(E.Properties.VariableNames, ...
        {'SourceLanguageName','Source','source','From','from','src','SourceLanguage'});
    tgtVar = pick_var(E.Properties.VariableNames, ...
        {'TargetLanguageName','Target','target','To','to','tgt','TargetLanguage'});
    wtVar  = pick_var(E.Properties.VariableNames, ...
        {'Coocurrences','Cooccurrences','Occurrences','occurrences','Weight','weight','Count','count','N','n'});

    if isempty(srcVar) || isempty(tgtVar) || isempty(wtVar)
        error(['Could not detect source/target/weight columns in %s.\n' ...
               'Found columns:\n  %s'], ...
               edgePath, strjoin(E.Properties.VariableNames, ', '));
    end

    src = normalize_label(E.(srcVar));
    tgt = normalize_label(E.(tgtVar));
    w   = double(E.(wtVar));

    keep = strlength(src) > 0 & strlength(tgt) > 0 & isfinite(w) & (w > opts.MinEdgeWeight);
    src = src(keep);
    tgt = tgt(keep);
    w   = w(keep);

    T = table(src, tgt, w, 'VariableNames', {'Source','Target','Weight'});

    % Aggregate duplicates
    [Gid, keys] = findgroups(T(:, {'Source','Target'}));
    aggW = splitapply(@sum, T.Weight, Gid);
    T = table(keys.Source, keys.Target, aggW, ...
        'VariableNames', {'Source','Target','Weight'});

    % ------------------ build node list and sparse W ------------------
    labels = unique([T.Source; T.Target], 'stable');
    n = numel(labels);

    idxMap = containers.Map(cellstr(labels), num2cell(1:n));

    sIdx = zeros(height(T),1);
    tIdx = zeros(height(T),1);
    for k = 1:height(T)
        sIdx(k) = idxMap(char(T.Source(k)));
        tIdx(k) = idxMap(char(T.Target(k)));
    end

    W = sparse(sIdx, tIdx, T.Weight, n, n);

    if opts.DropSelfLoops
        W = W - spdiags(diag(W), 0, n, n);
        W = sparse(W);
    end

    % ------------------ optional node filtering ------------------
    [inW, outW, u] = node_weights(W);

    keepNodes = true(n,1);

    if opts.MinNodeTotalWeight > 0
        keepNodes = keepNodes & (u >= opts.MinNodeTotalWeight);
    end

    if isfinite(opts.TopN) && opts.TopN < nnz(keepNodes)
        active = find(keepNodes);
        [~, ord] = sort(u(active), 'descend');
        keepTop = false(n,1);
        keepTop(active(ord(1:opts.TopN))) = true;
        keepNodes = keepNodes & keepTop;
    end

    if ~all(keepNodes)
        W = W(keepNodes, keepNodes);
        labels = labels(keepNodes);
        [inW, outW, u] = node_weights(W);
    end

    if opts.KeepLargestWeakComponent
        Gtmp = digraph(spones(W), cellstr(labels));
        bins = conncomp(Gtmp, 'Type', 'weak');
        counts = accumarray(bins(:), 1);
        [~, biggest] = max(counts);
        keepComp = (bins == biggest);

        W = W(keepComp, keepComp);
        labels = labels(keepComp);
        [inW, outW, u] = node_weights(W);
    end

    % ------------------ final edge and node tables ------------------
    [si, ti, wi] = find(W);
    edgeTable = table(labels(si), labels(ti), wi, ...
        'VariableNames', {'Source','Target','Weight'});

    nodeTable = table(labels, inW, outW, u, ...
        'VariableNames', {'Name','InWeight','OutWeight','TotalWeight'});

    % Try to merge optional node metadata
    if isfile(nodePath)
        try
            nodeMeta = read_tsv_robust(nodePath);
            nodeMeta = coerce_node_metadata(nodeMeta, labels);
            if ~isempty(nodeMeta)
                nodeTable = outerjoin(nodeTable, nodeMeta, ...
                    'Keys', 'Name', 'MergeKeys', true, 'Type', 'left');
            end
        catch ME
            warning('Could not merge node metadata from %s: %s', nodePath, ME.message);
        end
    end

    % MATLAB digraph
    G = digraph(W, cellstr(labels));

    % ------------------ output ------------------
    net = struct();
    net.W           = W;
    net.labels      = labels;
    net.G           = G;
    net.edgeTable   = edgeTable;
    net.nodeTable   = nodeTable;
    net.inWeight    = inW;
    net.outWeight   = outW;
    net.totalWeight = u;
    net.paths       = struct('edgeFile', edgePath, 'nodeFile', nodePath);
    net.opts        = opts;
end

% =======================================================================
function T = read_tsv_robust(pathStr)
    opts = detectImportOptions(pathStr, 'FileType', 'text', 'Delimiter', '\t');
    T = readtable(pathStr, opts);

    % Normalise variable names a little
    T.Properties.VariableNames = matlab.lang.makeValidName(T.Properties.VariableNames, ...
        'ReplacementStyle', 'delete');
end

function varName = pick_var(varNames, candidates)
    varName = '';
    lowerVars = lower(string(varNames));
    for i = 1:numel(candidates)
        hit = find(lowerVars == lower(string(candidates{i})), 1);
        if ~isempty(hit)
            varName = varNames{hit};
            return;
        end
    end
end

function x = normalize_label(x)
    if iscellstr(x) || isstring(x) || ischar(x) || iscategorical(x)
        x = string(x);
    elseif isnumeric(x)
        x = string(x);
    else
        x = string(x);
    end
    x = strip(x);
    x = regexprep(x, '\s+', ' ');
end

function [inW, outW, u] = node_weights(W)
    outW = full(sum(W, 2));
    inW  = full(sum(W, 1))';
    u    = inW + outW;
end

function nodeMeta = coerce_node_metadata(T, labels)
    % Attempt to identify a label/name column and keep rows matching labels.
    nodeMeta = table();

    vars = T.Properties.VariableNames;
    nameVar = pick_var(vars, {'Name','Label','Language','LanguageName','lang','node'});
    if isempty(nameVar)
        % fallback: first text-like column
        for k = 1:numel(vars)
            v = T.(vars{k});
            if iscellstr(v) || isstring(v) || iscategorical(v) || ischar(v)
                nameVar = vars{k};
                break;
            end
        end
    end

    if isempty(nameVar)
        return;
    end

    names = normalize_label(T.(nameVar));
    keep = ismember(names, labels);

    if ~any(keep)
        return;
    end

    T = T(keep, :);
    names = names(keep);

    % Remove duplicate names, keeping first occurrence
    [~, ia] = unique(names, 'stable');
    T = T(ia, :);
    names = names(ia);

    T.Name = names;

    % Put Name first; leave all other columns untouched
    otherVars = setdiff(T.Properties.VariableNames, {nameVar, 'Name'}, 'stable');
    nodeMeta = T(:, ['Name', otherVars]);
end