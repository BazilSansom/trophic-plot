function out = load_foodweb_FW24(opts)
%LOAD_FOODWEB_FW24  Download and parse FW24 food-web data for plotting.
%
% The source zip contains text files such as FW24.txt with rows like:
%   Macroph            DemHerb            2.576
%
% Output:
%   out.edge_table   : table with source_name, target_name, weight
%   out.nodes        : node names (cellstr)
%   out.node_table   : table(idx, name)
%   out.E            : numeric edge list [i j w]
%   out.W            : weighted adjacency matrix
%   out.files        : paths used
%   out.source       : provenance / citation info
%
% Example:
%   out = load_foodweb_FW24();
%   fig_food_web(out.E, struct('Export', false));

    if nargin < 1 || isempty(opts), opts = struct(); end
    opts = applyDefaults_(opts);

    if ~exist(opts.DataDir, 'dir')
        mkdir(opts.DataDir);
    end

    zipPath  = fullfile(opts.DataDir, char(opts.ZipFilename));
    unzipDir = fullfile(opts.DataDir, char(opts.UnzipFolder));

    % ---------------------------------------------------------------------
    % 1) Download zip if needed
    % ---------------------------------------------------------------------
    if opts.OverwriteDownload || ~exist(zipPath, 'file')
        fprintf('Downloading food-web dataset zip...\n');
        websave(zipPath, opts.ZipURL);
    else
        fprintf('Using existing zip: %s\n', zipPath);
    end

    % ---------------------------------------------------------------------
    % 2) Unzip if needed
    % ---------------------------------------------------------------------
    if opts.OverwriteUnzip && exist(unzipDir, 'dir')
        rmdir(unzipDir, 's');
    end

    txtPath = resolveDatasetFile_(unzipDir, opts.DatasetFilename);

    if isempty(txtPath) || ~exist(txtPath, 'file')
        if ~exist(unzipDir, 'dir')
            mkdir(unzipDir);
        end
        fprintf('Extracting zip to: %s\n', unzipDir);
        unzip(zipPath, unzipDir);

        % resolve again after extraction
        txtPath = resolveDatasetFile_(unzipDir, opts.DatasetFilename);
    else
        fprintf('Using existing extracted file: %s\n', txtPath);
    end

    if isempty(txtPath) || ~exist(txtPath, 'file')
        error('load_foodweb_FW24:DatasetMissing', ...
            'Could not find expected dataset file after unzip under: %s', unzipDir);
    end

    % ---------------------------------------------------------------------
    % 3) Parse FW24.txt
    % ---------------------------------------------------------------------
    edgeTable = parseFoodwebTxt_(txtPath);

    % ---------------------------------------------------------------------
    % 4) Build node list and numeric edge list
    % ---------------------------------------------------------------------
    allNames = unique([edgeTable.source_name; edgeTable.target_name], 'stable');
    n = numel(allNames);

    idxMap = containers.Map('KeyType', 'char', 'ValueType', 'double');
    for k = 1:n
        idxMap(char(allNames{k})) = k;
    end

    I = zeros(height(edgeTable), 1);
    J = zeros(height(edgeTable), 1);
    Wgt = edgeTable.weight;

    for r = 1:height(edgeTable)
        I(r) = idxMap(char(edgeTable.source_name{r}));
        J(r) = idxMap(char(edgeTable.target_name{r}));
    end

    E = [I, J, Wgt];
    W = sparse(I, J, Wgt, n, n);

    nodeTable = table((1:n)', allNames(:), ...
        'VariableNames', {'idx','name'});

    % ---------------------------------------------------------------------
    % 5) Return
    % ---------------------------------------------------------------------
    out = struct();
    out.edge_table = edgeTable;
    out.nodes      = allNames(:);
    out.node_table = nodeTable;
    out.E          = E;
    out.W          = W;

    out.files = struct();
    out.files.zip   = string(zipPath);
    out.files.unzip = string(unzipDir);
    out.files.txt   = string(txtPath);

    out.source = struct();
    out.source.dataset_id = 'FW24';
    out.source.zip_url = string(opts.ZipURL);
    out.source.dataset_filename = string(opts.DatasetFilename);
    out.source.primary_citation = 'Lin, Davis, Jordan, Liu (2024), Food Webs 38:e00336';
    out.source.note = ['Data retrieved from the InteractionFunctionalDiversityInFoodWebs repository ' ...
                       'and parsed from FW24.txt.'];
end

% =========================================================================
function opts = applyDefaults_(opts)
    thisDir  = fileparts(mfilename('fullpath'));
    repoRoot = getRepoRoot_(thisDir);

    d = struct();

    % Default storage location for downloaded + extracted food-web data
    d.DataDir = fullfile(repoRoot, 'paper', 'data', 'foodwebs');

    % Raw download URL, not GitHub blob URL
    d.ZipURL = "https://github.com/WeiChungLiu/InteractionFunctionalDiversityInFoodWebs/raw/refs/heads/main/FoodWebDataTraitData.zip";

    d.ZipFilename = "FoodWebDataSet.zip";
    d.UnzipFolder = "FoodWebDataSet";
    d.DatasetFilename = "FW24.txt";
    d.OverwriteDownload = false;
    d.OverwriteUnzip = false;

    fn = fieldnames(d);
    for k = 1:numel(fn)
        f = fn{k};
        if ~isfield(opts, f) || isempty(opts.(f))
            opts.(f) = d.(f);
        end
    end
end

%{
function opts = applyDefaults_(opts)
    thisDir  = fileparts(mfilename('fullpath'));
    repoRoot = getRepoRoot_(thisDir);

    d = struct();
    %d.DataDir = fullfile(repoRoot, 'data', 'foodwebs');
    d.DataDir = fullfile(repoRoot, 'paper', 'data', 'foodwebs');
    d.ZipURL = "https://github.com/WeiChungLiu/InteractionFunctionalDiversityInFoodWebs/raw/refs/heads/main/FoodWebDataTraitData.zip";
    %d.ZipURL = "https://github.com/WeiChungLiu/InteractionFunctionalDiversityInFoodWebs/blob/main/FoodWebDataTraitData.zip"
    % Also available here: "https://github.com/WeiChungLiu/Species-uniqueness-index-project/raw/refs/heads/main/FoodWebDataSet.zip";
    d.ZipFilename = "FoodWebDataSet.zip";
    d.UnzipFolder = "FoodWebDataSet";
    d.DatasetFilename = "FW24.txt";
    d.OverwriteDownload = false;
    d.OverwriteUnzip = false;

    fn = fieldnames(d);
    for k = 1:numel(fn)
        f = fn{k};
        if ~isfield(opts, f) || isempty(opts.(f))
            opts.(f) = d.(f);
        end
    end
end
%}

% =========================================================================
function T = parseFoodwebTxt_(txtPath)
% Expect whitespace-separated rows:
%   source   target   weight
%
% Handles arbitrary spacing between columns.

    fid = fopen(txtPath, 'r');
    if fid < 0
        error('Could not open file: %s', txtPath);
    end
    cleaner = onCleanup(@() fclose(fid));

    src = {};
    dst = {};
    w   = [];

    lineNo = 0;
    while true
        ln = fgetl(fid);
        if ~ischar(ln), break; end
        lineNo = lineNo + 1;

        ln = strtrim(ln);
        if isempty(ln)
            continue;
        end

        % Split on runs of whitespace
        parts = regexp(ln, '\s+', 'split');

        if numel(parts) < 3
            error('Bad row in %s at line %d: %s', txtPath, lineNo, ln);
        end

        % Assume last token is weight; first two are source/target
        % This matches lines like: Macroph   DemHerb   2.576
        wt = str2double(parts{end});
        if ~isfinite(wt)
            error('Could not parse weight in %s at line %d: %s', txtPath, lineNo, ln);
        end

        src{end+1,1} = parts{1}; %#ok<AGROW>
        dst{end+1,1} = parts{2}; %#ok<AGROW>
        w(end+1,1)   = wt; %#ok<AGROW>
    end

    T = table(src, dst, w, ...
        'VariableNames', {'source_name','target_name','weight'});
end

% =========================================================================
function repoRoot = getRepoRoot_(startDir)
    repoRoot = startDir;
    while true
        if exist(fullfile(repoRoot, '.git'), 'dir') || ...
           (exist(fullfile(repoRoot, 'src'), 'dir') && exist(fullfile(repoRoot, 'data'), 'dir'))
            return;
        end
        parent = fileparts(repoRoot);
        if strcmp(parent, repoRoot)
            error('Could not locate repo root starting from: %s', startDir);
        end
        repoRoot = parent;
    end
end

function txtPath = resolveDatasetFile_(rootDir, datasetFilename)
%RESOLVEDATASETFILE_  Find dataset file under rootDir, allowing nested folders.

    candidate = fullfile(rootDir, datasetFilename);
    if exist(candidate, 'file')
        txtPath = candidate;
        return;
    end

    hits = dir(fullfile(rootDir, '**', datasetFilename));

    if isempty(hits)
        txtPath = "";
        return;
    end

    if numel(hits) > 1
        warning('Multiple matches found for %s. Using first: %s', ...
            datasetFilename, fullfile(hits(1).folder, hits(1).name));
    end

    txtPath = fullfile(hits(1).folder, hits(1).name);
end