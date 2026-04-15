function [ds, info] = mangal_build_foodweb_ds(dataset_id, varargin)
%MANGAL_BUILD_FOODWEB_DS  Build a local MATLAB dataset of all networks in a Mangal dataset.
%
%   [ds, info] = mangal_build_foodweb_ds(dataset_id, ...)
%
% Returns:
%   ds(k) : struct with fields:
%       .dataset_id, .network_id, .name, .date, .description, .public
%       .nodes, .node_ids, .W, .n, .m, .edges
%   info  : struct with fields:
%       .networks (raw network listing)
%       .failed   (table of failed network_ids + errors)
%
% Options (Name/Value):
%   'Count'         : page size for network listing (default 1000)
%   'PauseSeconds'  : pause between network pulls (default 0.10)
%   'MaxNetworks'   : cap number of networks to pull (default Inf)
%   'OnlyPublic'    : true/false (default true)
%   'FlipEdge'      : passed to mangal_pull_network (default true)
%   'DefaultWeight' : passed to mangal_pull_network (default 1)
%   'Type'          : passed to mangal_pull_network (default '')
%   'Retry'         : retries per failed network pull (default 2)
%   'Verbose'       : true/false (default true)
%   'SaveEvery'     : save every N networks (default 0 = never)
%   'OutFile'       : path to .mat for incremental saving (default '')
%
% Requires:
%   mangal_api_get.m
%   mangal_pull_network.m

p = inputParser;
p.addParameter('Count',        1000, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('PauseSeconds', 0.10, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('MaxNetworks',  inf,  @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('OnlyPublic',   true, @(x)islogical(x)&&isscalar(x));
p.addParameter('FlipEdge',     true, @(x)islogical(x)&&isscalar(x));
p.addParameter('DefaultWeight',1,    @(x)isnumeric(x)&&isscalar(x));
p.addParameter('Type',         "",   @(s)ischar(s)||isstring(s));
p.addParameter('Retry',        2,    @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('Verbose',      true, @(x)islogical(x)&&isscalar(x));
p.addParameter('SaveEvery',    0,    @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('OutFile',      "",   @(s)ischar(s)||isstring(s));
p.parse(varargin{:});
opt = p.Results;

% --- 1) list networks in dataset (paged) ---
nets = mangal_list_networks_(dataset_id, 'Count', opt.Count);

if opt.OnlyPublic && ~isempty(nets) && isfield(nets,'public')
    nets = nets([nets.public] == true);
end

if isempty(nets)
    ds = struct([]);
    info = struct('networks', nets, 'failed', table());
    return;
end

if isfinite(opt.MaxNetworks)
    nets = nets(1:min(numel(nets), opt.MaxNetworks));
end

% Preallocate ds struct array
ds = repmat(struct( ...
    'dataset_id', dataset_id, ...
    'network_id', NaN, ...
    'name', "", ...
    'date', "", ...
    'description', "", ...
    'public', true, ...
    'nodes', {{}}, ...
    'node_ids', [], ...
    'W', sparse(0,0), ...
    'n', 0, ...
    'm', 0, ...
    'edges', table()), numel(nets), 1);

failed_id   = [];
failed_msg  = strings(0,1);

% --- 2) pull each network ---
for k = 1:numel(nets)
    nid = nets(k).id;

    if opt.Verbose
        fprintf('[%d/%d] network_id=%d  name=%s\n', k, numel(nets), nid, string(nets(k).name));
    end

    ok = false;
    lastErr = [];

    for attempt = 0:opt.Retry
        try
            out = mangal_pull_network(nid, ...
                'FlipEdge',      opt.FlipEdge, ...
                'DefaultWeight', opt.DefaultWeight, ...
                'Type',          opt.Type);
            ok = true;
            break;
        catch ME
            lastErr = ME;
            if opt.Verbose
                fprintf('  attempt %d failed: %s\n', attempt+1, ME.message);
            end
            pause(0.25);
        end
    end

    if ~ok
        failed_id(end+1,1)  = nid; %#ok<AGROW>
        failed_msg(end+1,1) = string(lastErr.message); %#ok<AGROW>
        continue;
    end

    % --- pack ds(k) ---
    ds(k).dataset_id   = dataset_id;
    ds(k).network_id   = nid;
    if isfield(nets(k),'name'),        ds(k).name        = string(nets(k).name); end
    if isfield(nets(k),'date') && ~isempty(nets(k).date), ds(k).date = string(nets(k).date); end
    if isfield(nets(k),'description'), ds(k).description = string(nets(k).description); end
    if isfield(nets(k),'public'),      ds(k).public      = logical(nets(k).public); end

    ds(k).nodes    = out.nodes;
    ds(k).node_ids = out.node_ids;
    ds(k).W        = out.W;
    ds(k).n        = size(out.W,1);
    ds(k).m        = nnz(out.W);
    ds(k).edges    = out.edges;

    if opt.PauseSeconds > 0
        pause(opt.PauseSeconds);
    end

    % --- incremental save ---
    if opt.SaveEvery > 0 && strlength(string(opt.OutFile)) > 0
        if mod(k, opt.SaveEvery) == 0
            outFile = char(opt.OutFile);
            try
                save(outFile, 'ds', 'nets', 'dataset_id', '-v7.3');
                if opt.Verbose, fprintf('  saved -> %s\n', outFile); end
            catch ME
                warning('mangal_build_foodweb_ds:SaveFailed', 'Save failed: %s', ME.message);
            end
        end
    end
end

% Remove empty slots where pulls failed (optional but usually nicer)
okMask = ~isnan([ds.network_id]);
ds = ds(okMask);

info = struct();
info.networks = nets;
if isempty(failed_id)
    info.failed = table();
else
    info.failed = table(failed_id, failed_msg, 'VariableNames', {'network_id','message'});
end

end

% =====================================================================
function nets = mangal_list_networks_(dataset_id, varargin)
p = inputParser;
p.addParameter('Count', 1000, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.parse(varargin{:});
count = p.Results.Count;

nets = struct([]);
page = 0;

while true
    batch = mangal_api_get("network", struct('dataset_id', dataset_id, 'count', count, 'page', page));
    if isempty(batch)
        break;
    end
    if isempty(nets)
        nets = batch;
    else
        nets = [nets; batch]; %#ok<AGROW>
    end
    if numel(batch) < count
        break;
    end
    page = page + 1;
end
end