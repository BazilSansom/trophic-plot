function [net, meta, files] = mangal_get_example_network(dataset_id, network_id, varargin)
%MANGAL_GET_EXAMPLE_NETWORK  Pull one chosen Mangal network with metadata.
%
%   [net, meta, files] = mangal_get_example_network(dataset_id, network_id)
%
% Returns:
%   net   : struct with fields
%       .dataset_id, .network_id, .name, .date, .description, .public
%       .nodes, .node_ids, .W, .n, .m, .edges
%
%   meta  : struct with fields
%       .dataset_record
%       .network_record
%       .reference_record
%       .provenance
%
%   files : struct with cache paths used/written
%
% Requires:
%   mangal_api_get.m
%   mangal_pull_network.m

p = inputParser;
p.addParameter('CacheDir', "", @(s)ischar(s)||isstring(s));
p.addParameter('UseCache', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('SaveJSON', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('FlipEdge', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('DefaultWeight', 1, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('Type', "", @(s)ischar(s)||isstring(s));
p.addParameter('Retry', 2, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('Verbose', true, @(x)islogical(x)&&isscalar(x));
p.parse(varargin{:});
opt = p.Results;

files = struct('mat',"",'json',"");

cacheBase = "";
if strlength(string(opt.CacheDir)) > 0
    if ~exist(opt.CacheDir, 'dir')
        mkdir(opt.CacheDir);
    end
    cacheBase = fullfile(char(opt.CacheDir), ...
        sprintf('mangal_dataset_%d_network_%d', dataset_id, network_id));
    files.mat  = string(cacheBase) + ".mat";
    files.json = string(cacheBase) + "_meta.json";
end

% ------------------------ Cache load -------------------------------------
if opt.UseCache && strlength(string(cacheBase)) > 0 && exist(files.mat, 'file') == 2
    S = load(files.mat, 'net', 'meta');
    net  = S.net;
    meta = S.meta;
    if opt.Verbose
        fprintf('Loaded cached example network: %s\n', files.mat);
    end
    return;
end

% ------------------------ Dataset record ---------------------------------
dataset_record = mangal_get_dataset_record_(dataset_id);

% ------------------------ Network record ---------------------------------
network_record = mangal_get_network_record_(dataset_id, network_id);

% ------------------------ Pull adjacency ---------------------------------
ok = false;
lastErr = [];
for attempt = 0:opt.Retry
    try
        out = mangal_pull_network(network_id, ...
            'FlipEdge',      opt.FlipEdge, ...
            'DefaultWeight', opt.DefaultWeight, ...
            'Type',          opt.Type);
        ok = true;
        break;
    catch ME
        lastErr = ME;
        if opt.Verbose
            fprintf('mangal_pull_network failed (attempt %d): %s\n', ...
                attempt+1, ME.message);
        end
        pause(0.25);
    end
end

if ~ok
    error('mangal_get_example_network:PullFailed', ...
        'Failed to pull network_id=%d: %s', network_id, lastErr.message);
end

% ------------------------ Pack net ---------------------------------------
net = struct();
net.dataset_id   = dataset_id;
net.network_id   = network_id;
net.name         = getfield_or_string_(network_record, 'name', "");
net.date         = getfield_or_string_(network_record, 'date', "");
net.description  = getfield_or_string_(network_record, 'description', "");
net.public       = getfield_or_logical_(network_record, 'public', true);

net.nodes        = out.nodes;
net.node_ids     = out.node_ids;
net.W            = out.W;
net.n            = size(out.W,1);
net.m            = nnz(out.W);
net.edges        = out.edges;

% ------------------------ Reference record -------------------------------
ref_id = getfield_or_numeric_(dataset_record, 'ref_id', NaN);
reference_record = struct();
if ~isnan(ref_id)
    reference_record = mangal_get_reference_record_(ref_id);
end

% ------------------------ Metadata ---------------------------------------
meta = struct();
meta.dataset_record   = dataset_record;
meta.network_record   = network_record;
meta.reference_record = reference_record;

meta.provenance = struct();
meta.provenance.generated_utc = char(datetime('now', ...
    'TimeZone','UTC', 'Format','yyyy-MM-dd''T''HH:mm:ss''Z'''));
meta.provenance.dataset_id    = dataset_id;
meta.provenance.network_id    = network_id;
meta.provenance.flip_edge     = opt.FlipEdge;
meta.provenance.default_weight= opt.DefaultWeight;
meta.provenance.type          = char(string(opt.Type));

% ------------------------ Save cache -------------------------------------
if strlength(string(cacheBase)) > 0
    save(files.mat, 'net', 'meta', '-v7.3');

    if opt.SaveJSON
        meta_json = build_jsonable_meta_(net, meta);
        fid = fopen(char(files.json), 'w');
        if fid >= 0
            fprintf(fid, '%s', jsonencode(meta_json));
            fclose(fid);
        end
    end

    if opt.Verbose
        fprintf('Saved example network cache: %s\n', files.mat);
    end
end

end

% ========================================================================
function dataset_record = mangal_get_dataset_record_(dataset_id)

dataset_record = struct();

% Try direct endpoint first
try
    rec = mangal_api_get(sprintf('dataset/%d', dataset_id), struct());
    if ~isempty(rec)
        dataset_record = rec;
        return;
    end
catch
end

% Fallback: query listing and pick exact id
try
    rec = mangal_api_get("dataset", struct('id', dataset_id, 'count', 1000, 'page', 0));
    dataset_record = pick_exact_id_(rec, dataset_id);
catch
    dataset_record = struct();
end

end

% ========================================================================
function network_record = mangal_get_network_record_(dataset_id, network_id)

network_record = struct();
page = 0;
count = 1000;

while true
    batch = mangal_api_get("network", ...
        struct('dataset_id', dataset_id, 'count', count, 'page', page));

    if isempty(batch)
        break;
    end

    ids = extract_ids_(batch);
    hit = find(ids == network_id, 1, 'first');
    if ~isempty(hit)
        network_record = batch(hit);
        return;
    end

    if numel(batch) < count
        break;
    end
    page = page + 1;
end

error('mangal_get_example_network:NetworkNotFound', ...
    'network_id=%d not found in dataset_id=%d.', network_id, dataset_id);

end

% ========================================================================
function reference_record = mangal_get_reference_record_(ref_id)

reference_record = struct();

try
    rec = mangal_api_get(sprintf('reference/%d', ref_id), struct());
    if ~isempty(rec)
        reference_record = rec;
        return;
    end
catch
end

try
    rec = mangal_api_get("reference", struct('id', ref_id, 'count', 1000, 'page', 0));
    reference_record = pick_exact_id_(rec, ref_id);
catch
    reference_record = struct();
end

end

% ========================================================================
function ids = extract_ids_(S)
if isempty(S)
    ids = [];
    return;
end
ids = nan(numel(S),1);
for k = 1:numel(S)
    if isfield(S(k),'id') && ~isempty(S(k).id)
        ids(k) = double(S(k).id);
    end
end
end

% ========================================================================
function rec = pick_exact_id_(S, target_id)
rec = struct();
if isempty(S), return; end

if isstruct(S) && numel(S) == 1 && isfield(S,'id') && isequal(double(S.id), double(target_id))
    rec = S;
    return;
end

if isstruct(S)
    ids = extract_ids_(S);
    hit = find(ids == target_id, 1, 'first');
    if ~isempty(hit)
        rec = S(hit);
    end
end
end

% ========================================================================
function x = getfield_or_string_(s, f, defaultVal)
x = string(defaultVal);
if isstruct(s) && isfield(s, f)
    v = s.(f);
    if ~isempty(v)
        x = string(v);
    end
end
end

function x = getfield_or_numeric_(s, f, defaultVal)
x = defaultVal;
if isstruct(s) && isfield(s, f)
    v = s.(f);
    if ~isempty(v) && isnumeric(v) && isscalar(v)
        x = double(v);
    end
end
end

function x = getfield_or_logical_(s, f, defaultVal)
x = defaultVal;
if isstruct(s) && isfield(s, f)
    v = s.(f);
    if ~isempty(v)
        x = logical(v);
    end
end
end

% ========================================================================
function meta_json = build_jsonable_meta_(net, meta)

meta_json = struct();

meta_json.dataset_id   = net.dataset_id;
meta_json.network_id   = net.network_id;
meta_json.network_name = char(string(net.name));
meta_json.network_date = char(string(net.date));
meta_json.description  = char(string(net.description));
meta_json.n            = net.n;
meta_json.m            = net.m;

meta_json.dataset_record   = jsonable_struct_(meta.dataset_record);
meta_json.network_record   = jsonable_struct_(meta.network_record);
meta_json.reference_record = jsonable_struct_(meta.reference_record);
meta_json.provenance       = jsonable_struct_(meta.provenance);

end

% ========================================================================
function out = jsonable_struct_(s)
out = struct();
if ~isstruct(s), return; end
fns = fieldnames(s);
for i = 1:numel(fns)
    f = fns{i};
    v = s.(f);
    if ischar(v) || isstring(v) || isnumeric(v) || islogical(v)
        out.(f) = v;
    end
end
end