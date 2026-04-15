function [allDs, info] = mangal_build_all_ds(varargin)
%MANGAL_BUILD_ALL_DS  Pull many Mangal datasets (and their networks) into one MATLAB struct array.
%
%   [allDs, info] = mangal_build_all_ds()
%   [allDs, info] = mangal_build_all_ds('MaxDatasets', 5, 'MaxNetworks', 10)
%   [allDs, info] = mangal_build_all_ds('AllDatasets', true)
%   [allDs, info] = mangal_build_all_ds('DatasetIDs', [15 16 17], 'MaxNetworks', 50)
%
% OUTPUT
%   allDs : struct array concatenating ds returned by mangal_build_foodweb_ds for each dataset.
%
%   info  : struct with fields:
%       .datasets        : raw dataset listing (struct array)
%       .datasetsTable   : table summary (id,name,description,public,ref_id,user_id)
%       .datasetSummary  : per-dataset table (id,name,nNetworksListed,nPulled,nFailed)
%       .failedDatasets  : table of dataset-level failures (dataset_id,message)
%       .failedNetworks  : table of network-level failures aggregated across datasets
%
% OPTIONS (Name/Value)
%   Dataset selection:
%     'AllDatasets'     : false (default). If true, ignores MaxDatasets and pulls all.
%     'MaxDatasets'     : 5 (default). Number of datasets to pull (in ascending dataset id order).
%     'DatasetIDs'      : [] (default). If provided, pull exactly these dataset ids.
%     'OnlyPublicDatasets' : true (default). Filters dataset listing by public flag when available.
%
%   Passed through to mangal_build_foodweb_ds (per dataset):
%     'Count'           : default 1000
%     'PauseSeconds'    : default 0.10
%     'MaxNetworks'     : default Inf   (cap networks per dataset)
%     'OnlyPublic'      : default true  (networks)
%     'FlipEdge'        : default true
%     'DefaultWeight'   : default 1
%     'Type'            : default ""
%     'Retry'           : default 2
%     'Verbose'         : default true
%
%   Saving:
%     'SaveEveryDataset': 0 (default). If >0, save after every N datasets.
%     'OutFile'         : "" (default). If set and SaveEveryDataset>0, writes incremental .mat.
%
% NOTES
% - Dataset listing is paged. This function pulls all datasets first, then selects a subset.
% - Some datasets are not predator-prey "food webs" (pollination, seed dispersal, etc.).
%   That’s fine if your goal is “find compelling weighted-vs-topology examples”.
%
% Requires:
%   mangal_api_get.m
%   mangal_build_foodweb_ds.m
%   mangal_pull_network.m

% ------------------------ Parse options ---------------------------------------
p = inputParser;

% Dataset selection
p.addParameter('AllDatasets', false, @(x)islogical(x)&&isscalar(x));
p.addParameter('MaxDatasets', 5, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('DatasetIDs', [], @(x)isnumeric(x) && isvector(x));
p.addParameter('OnlyPublicDatasets', true, @(x)islogical(x)&&isscalar(x));

% Per-dataset (network pull) options (passed through)
p.addParameter('Count',        1000, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('PauseSeconds', 0.10, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('MaxNetworks',  inf,  @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('OnlyPublic',   true, @(x)islogical(x)&&isscalar(x));
p.addParameter('FlipEdge',     true, @(x)islogical(x)&&isscalar(x));
p.addParameter('DefaultWeight',1,    @(x)isnumeric(x)&&isscalar(x));
p.addParameter('Type',         "",   @(s)ischar(s)||isstring(s));
p.addParameter('Retry',        2,    @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('Verbose',      true, @(x)islogical(x)&&isscalar(x));

% Saving
p.addParameter('SaveEveryDataset', 0, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('OutFile', "", @(s)ischar(s)||isstring(s));

p.parse(varargin{:});
opt = p.Results;

% ------------------------ 1) List datasets -----------------------------------
datasets = mangal_list_datasets_('Count', 1000);  % big page size

if opt.OnlyPublicDatasets && ~isempty(datasets) && isfield(datasets,'public')
    datasets = datasets([datasets.public] == true);
end

% Make a table summary for convenience
datasetsTable = struct2table(datasets);
keepCols = intersect(datasetsTable.Properties.VariableNames, ...
    {'id','name','description','public','ref_id','user_id'}, 'stable');
datasetsTable = datasetsTable(:, keepCols);

% ------------------------ 2) Choose dataset IDs -------------------------------
if ~isempty(opt.DatasetIDs)
    dataset_ids = unique(opt.DatasetIDs(:).', 'stable');
else
    % default: smallest ids first
    if isempty(datasets)
        allDs = struct([]);
        info = struct('datasets',datasets,'datasetsTable',datasetsTable, ...
                      'datasetSummary',table(), ...
                      'failedDatasets',table(), ...
                      'failedNetworks',table());
        return;
    end
    dataset_ids = sort([datasets.id]);

    if ~opt.AllDatasets
        dataset_ids = dataset_ids(1:min(numel(dataset_ids), opt.MaxDatasets));
    end
end

% Map dataset_id -> dataset name/desc if available
id2meta = containers.Map('KeyType','double','ValueType','any');
for k = 1:numel(datasets)
    id2meta(datasets(k).id) = datasets(k);
end

% ------------------------ 3) Pull datasets -----------------------------------
allDs = struct([]);
failedDataset_id  = [];
failedDataset_msg = strings(0,1);

% per-dataset summary accumulators
sum_id = zeros(0,1);
sum_name = strings(0,1);
sum_listed = zeros(0,1);
sum_pulled = zeros(0,1);
sum_failed = zeros(0,1);

% aggregate network failures across datasets
failedNetworksAll = table();

for di = 1:numel(dataset_ids)
    dataset_id = dataset_ids(di);

    dname = "";
    if isKey(id2meta, dataset_id)
        meta = id2meta(dataset_id);
        if isfield(meta,'name') && ~isempty(meta.name)
            dname = string(meta.name);
        end
    end

    if opt.Verbose
        fprintf('\n=== DATASET [%d/%d] id=%d  name=%s ===\n', di, numel(dataset_ids), dataset_id, dname);
    end

    try
        [ds, infods] = mangal_build_foodweb_ds(dataset_id, ...
            'Count',         opt.Count, ...
            'PauseSeconds',  opt.PauseSeconds, ...
            'MaxNetworks',   opt.MaxNetworks, ...
            'OnlyPublic',    opt.OnlyPublic, ...
            'FlipEdge',      opt.FlipEdge, ...
            'DefaultWeight', opt.DefaultWeight, ...
            'Type',          opt.Type, ...
            'Retry',         opt.Retry, ...
            'Verbose',       opt.Verbose, ...
            'SaveEvery',     0, ...         % disable inner incremental saves
            'OutFile',       "" );
    catch ME
        failedDataset_id(end+1,1)  = dataset_id; %#ok<AGROW>
        failedDataset_msg(end+1,1) = string(ME.message); %#ok<AGROW>
        if opt.Verbose
            fprintf('!! dataset failed: %s\n', ME.message);
        end
        continue;
    end

    % merge ds into allDs
    if ~isempty(ds)
        if isempty(allDs)
            allDs = ds;
        else
            allDs = [allDs; ds]; %#ok<AGROW>
        end
    end

    % dataset summary numbers
    nListed = 0;
    if isfield(infods,'networks') && ~isempty(infods.networks)
        nListed = numel(infods.networks);
    end
    nPulled = numel(ds);
    nFail   = 0;
    if isfield(infods,'failed') && ~isempty(infods.failed) && istable(infods.failed)
        nFail = height(infods.failed);
        % attach dataset_id for aggregation
        tmp = infods.failed;
        tmp.dataset_id = repmat(dataset_id, height(tmp), 1);
        failedNetworksAll = [failedNetworksAll; tmp]; %#ok<AGROW>
    end

    sum_id(end+1,1)     = dataset_id; %#ok<AGROW>
    sum_name(end+1,1)   = dname;      %#ok<AGROW>
    sum_listed(end+1,1) = nListed;    %#ok<AGROW>
    sum_pulled(end+1,1) = nPulled;    %#ok<AGROW>
    sum_failed(end+1,1) = nFail;      %#ok<AGROW>

    % incremental save (per dataset)
    if opt.SaveEveryDataset > 0 && strlength(string(opt.OutFile)) > 0
        if mod(di, opt.SaveEveryDataset) == 0
            outFile = char(opt.OutFile);
            try
                save(outFile, 'allDs', 'datasets', 'datasetsTable', ...
                    'dataset_ids', 'sum_id','sum_name','sum_listed','sum_pulled','sum_failed', ...
                    'failedDataset_id','failedDataset_msg','failedNetworksAll', '-v7.3');
                if opt.Verbose
                    fprintf('  saved -> %s\n', outFile);
                end
            catch ME
                warning('mangal_build_all_ds:SaveFailed', 'Save failed: %s', ME.message);
            end
        end
    end
end

% ------------------------ 4) Assemble info -----------------------------------
datasetSummary = table(sum_id, sum_name, sum_listed, sum_pulled, sum_failed, ...
    'VariableNames', {'dataset_id','name','nNetworksListed','nPulled','nFailed'});

if isempty(failedDataset_id)
    failedDatasets = table();
else
    failedDatasets = table(failedDataset_id, failedDataset_msg, ...
        'VariableNames', {'dataset_id','message'});
end

info = struct();
info.datasets       = datasets;
info.datasetsTable  = datasetsTable;
info.datasetSummary = datasetSummary;
info.failedDatasets = failedDatasets;
info.failedNetworks = failedNetworksAll;

end

% =====================================================================
function dsAll = mangal_list_datasets_(varargin)
%MANGAL_LIST_DATASETS_  Page through /dataset endpoint.
p = inputParser;
p.addParameter('Count', 1000, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('MaxPages', Inf, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.parse(varargin{:});
count = p.Results.Count;
maxPages = p.Results.MaxPages;

dsAll = struct([]);
page = 0;

while page < maxPages
    batch = mangal_api_get("dataset", struct('count', count, 'page', page));
    if isempty(batch)
        break;
    end
    if isempty(dsAll)
        dsAll = batch;
    else
        dsAll = [dsAll; batch]; %#ok<AGROW>
    end
    if numel(batch) < count
        break;
    end
    page = page + 1;
end
end