function out = fig_early_picture_foodweb_from_mangal(dataset_id, network_id, opts)
%FIG_EARLY_PICTURE_FOODWEB_FROM_MANGAL
% Convenience wrapper: fetch one network, then plot the 1x3 early visual.

if nargin < 3 || isempty(opts), opts = struct(); end

% Cache raw example in a stable place unless caller overrides
if ~isfield(opts, 'CacheDir')
    opts.CacheDir = fullfile('data', 'mangal_examples');
end

[net, meta, sourceFiles] = mangal_get_example_network(dataset_id, network_id, ...
    'CacheDir',      opts.CacheDir, ...
    'UseCache',      true, ...
    'SaveJSON',      true, ...
    'FlipEdge',      true, ...
    'DefaultWeight', 1, ...
    'Verbose',       true);

% Plot using your pure plotting function
out = fig_early_picture_foodweb_net(net, opts);

% Attach source metadata
out.source_meta  = meta;
out.source_files = sourceFiles;

% Ensure out.files exists
if ~isfield(out, 'files') || ~isstruct(out.files)
    out.files = struct();
end

% Expose raw cached source files too
if isfield(sourceFiles, 'mat') && strlength(string(sourceFiles.mat)) > 0
    out.files.source_mat = string(sourceFiles.mat);
end
if isfield(sourceFiles, 'json') && strlength(string(sourceFiles.json)) > 0
    out.files.source_json = string(sourceFiles.json);
end

% Optionally write figure-specific metadata sidecars next to exported figure
if isfield(opts,'Export') && opts.Export && ...
   isfield(opts,'ExportDir') && isfield(opts,'ExportBase')

    [metaMatPath, metaJsonPath] = writeFigureMetaSidecar_( ...
        meta, opts.ExportDir, opts.ExportBase);

    out.files.meta_mat  = string(metaMatPath);
    out.files.meta_json = string(metaJsonPath);
end

end

function [matPath, jsonPath] = writeFigureMetaSidecar_(meta, exportDir, exportBase)

if ~exist(exportDir, 'dir')
    mkdir(exportDir);
end

matPath  = fullfile(exportDir, [exportBase '_meta.mat']);
jsonPath = fullfile(exportDir, [exportBase '_meta.json']);

save(matPath, 'meta');

fid = fopen(jsonPath, 'w');
if fid < 0
    error('Could not open metadata JSON for writing: %s', jsonPath);
end
fprintf(fid, '%s', jsonencode(meta));
fclose(fid);

end