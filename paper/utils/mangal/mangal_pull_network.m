function out = mangal_pull_network(network_id, varargin)
%MANGAL_PULL_NETWORK  Fetch nodes+interactions for one network_id, build W.
%
% out.W        sparse adjacency (n x n)
% out.nodes    cellstr node labels (n x 1)
% out.node_ids numeric Mangal node ids (n x 1)
% out.edges    table(node_from,node_to,value,type,attr_id)
%
% Options:
%   'FlipEdge' (default true)  % if Mangal node_from is predator and node_to prey,
%                              % FlipEdge=true creates prey->predator edges.
%   'DefaultWeight' (default 1) for missing/NaN values
%   'Type' (default '') e.g. 'predation' to filter interactions
%
% Mangal interactions endpoint supports network_id, type, etc. :contentReference[oaicite:5]{index=5}

    p = inputParser;
    p.addParameter('FlipEdge', true, @(x)islogical(x)&&isscalar(x));
    p.addParameter('DefaultWeight', 1, @(x)isnumeric(x)&&isscalar(x));
    p.addParameter('Type', '', @(s)ischar(s)||isstring(s));
    p.parse(varargin{:});
    opt = p.Results;

    % --- nodes
    nodes = mangal_api_get("node", struct('network_id', network_id, 'count', 1000, 'page', 0));
    
    disp(class(nodes));
    whos nodes
    
    node_ids = [nodes.id]';
    n = numel(node_ids);

    % label preference: taxonomy.name if present, else original_name
    labels = cell(n,1);
    for i = 1:n
        if isfield(nodes(i),'taxonomy') && ~isempty(nodes(i).taxonomy) && isfield(nodes(i).taxonomy,'name') && ~isempty(nodes(i).taxonomy.name)
            labels{i} = string(nodes(i).taxonomy.name);
        else
            labels{i} = string(nodes(i).original_name);
        end
    end
    labels = cellstr(labels);

    id2ix = containers.Map('KeyType','double','ValueType','double');
    for i = 1:n
        id2ix(node_ids(i)) = i;
    end

    % --- interactions
    q = struct('network_id', network_id, 'count', 1000, 'page', 0);
    if strlength(string(opt.Type))>0
        q.type = char(opt.Type);
    end
    inter = mangal_api_get("interaction", q);

    % build edge arrays
    m = numel(inter);
    si = zeros(m,1); sj = zeros(m,1); w = zeros(m,1);
    typ = strings(m,1); aid = nan(m,1);

    keep = true(m,1);
    for e = 1:m
        if ~isKey(id2ix, inter(e).node_from) || ~isKey(id2ix, inter(e).node_to)
            keep(e) = false; continue;
        end

        i = id2ix(inter(e).node_from);
        j = id2ix(inter(e).node_to);

        val = inter(e).value;
        if isempty(val) || isnan(val)
            val = opt.DefaultWeight;
        end

        if opt.FlipEdge
            % prey -> predator (if node_from is predator, node_to is prey)
            si(e) = j; sj(e) = i;
        else
            si(e) = i; sj(e) = j;
        end
        w(e) = val;

        if isfield(inter(e),'type') && ~isempty(inter(e).type), typ(e) = string(inter(e).type); end
        if isfield(inter(e),'attr_id') && ~isempty(inter(e).attr_id), aid(e) = inter(e).attr_id; end
    end

    si = si(keep); sj = sj(keep); w = w(keep); typ = typ(keep); aid = aid(keep);

    % drop exactly-zero weights (treat as absent edges)
    nz = (w ~= 0);
    si = si(nz); sj = sj(nz); w = w(nz); typ = typ(nz); aid = aid(nz);

    W = sparse(si, sj, w, n, n);

    out = struct();
    out.network_id = network_id;
    out.W = W;
    out.nodes = labels;
    out.node_ids = node_ids;
    out.edges = table(si, sj, w, typ, aid, 'VariableNames', ...
        {'i','j','value','type','attr_id'});
end