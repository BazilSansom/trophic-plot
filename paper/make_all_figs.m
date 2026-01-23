function out = make_all_figs(opts)
%MAKE_ALL_FIGS  Build all paper figures reproducibly into paper/outputs/figs.
%
% out = make_all_figs()
% out = make_all_figs(struct('ExportPNG',true,'PNGResolution',300))
%
% Assumes figure generators are callable functions (e.g. fig_cyclic_comparison)
% and support opts fields: Export, ExportDir, ExportBase, ExportPDF, ExportPNG,
% PNGResolution, PDFVector (as in fig_cyclic_comparison).
%
% Returns:
%   out.exportDir
%   out.results  : struct of per-figure outputs
%   out.files    : table-like struct array listing produced files

    if nargin < 1 || isempty(opts), opts = struct(); end

    % ------------------ defaults ------------------
    def.ExportPDF     = true;
    def.ExportPNG     = true;
    def.PNGResolution = 300;
    def.PDFVector     = true;

    % Control whether figures pop up while building
    def.CloseFigs     = true;

    opts = applyDefaults_(opts, def);

    % ------------------ paths ------------------
    repoRoot = fileparts(fileparts(mfilename('fullpath'))); % paper/ is one level below root
    srcDir   = fullfile(repoRoot, 'src');

    addpath(genpath(srcDir));  % toolbox
    addpath(genpath(fullfile(repoRoot, 'paper'))); % paper helpers + fig fns (if stored under paper)

    exportDir = fullfile(repoRoot, 'paper', 'outputs', 'figs');
    if ~exist(exportDir, 'dir'), mkdir(exportDir); end

    % ------------------ build all figures ------------------
    results = struct();
    files = [];

    % Helper to call a fig fn with standard export args
    function r = callFig(figFn, baseName, extraOpts)
        if nargin < 3, extraOpts = struct(); end
        o = extraOpts;

        % Standard export options
        o.Export        = true;
        o.ExportDir     = exportDir;
        o.ExportBase    = baseName;
        o.ExportPDF     = opts.ExportPDF;
        o.ExportPNG     = opts.ExportPNG;
        o.PNGResolution = opts.PNGResolution;
        o.PDFVector     = opts.PDFVector;

        r = feval(figFn, o);

        % Collect file paths if present
        if isfield(r,'files')
            if isfield(r.files,'pdf') && strlength(string(r.files.pdf)) > 0
                files = [files; struct('figure',baseName,'type',"pdf",'path',string(r.files.pdf))]; %#ok<AGROW>
            end
            if isfield(r.files,'png') && strlength(string(r.files.png)) > 0
                files = [files; struct('figure',baseName,'type',"png",'path',string(r.files.png))]; %#ok<AGROW>
            end
        end

        if opts.CloseFigs && isfield(r,'fig') && ishghandle(r.fig)
            close(r.fig);
        end
    end

    paperPlotStyle = struct( ...
        'ArrowSize',     0.10, ...      % try 0.10–0.14 if you want larger
        'ArrowSizeMode', 'fixed', ...
        'ArrowMinFrac',  0.03, ...      % keep a sensible floor
        'ArrowMaxFrac',  0 ...          % disable cap for paper stability
    );


    % ---- LIST PAPER FIGURES HERE ----
    % Use canonical basenames that match what Overleaf references.
    results.fig_cyclic_comparison = callFig('fig_cyclic_comparison', 'fig_cyclic_comparison');
    results.fig_general_dag_comparison = callFig('fig_general_dag_comparison','fig_general_dag_comparison');
    results.fig_topology_vs_weighted_flow = callFig('fig_topology_vs_weighted_flow','fig_topology_vs_weighted_flow');
    results.fig_perfectly_layered_comparison = callFig('fig_perfectly_layered_comparison','fig_perfectly_layered_comparison');
    results.fig_ordinal_metric_scatter = callFig( ...
        'fig_scatter_tfl_vs_layered', ...
        'fig_ordinal_metric_scatter', ...
        struct( ...
            'generator','gppm','n',60,'p',0.02,'T',10, ...
            'gppmArgs', {{'ForceDAG',true,'Weighted',false}}, ...
            'showLayouts',true,'seed',13,'visible',false, ...
            'PaperPlotStyle', paperPlotStyle ...
        ) ...
    );
    
    % ------------------ return ------------------
    out = struct();
    out.exportDir = exportDir;
    out.results   = results;
    out.files     = files;
end

% ======================================================================
function s = applyDefaults_(s, def)
    fn = fieldnames(def);
    for k = 1:numel(fn)
        f = fn{k};
        if ~isfield(s,f) || isempty(s.(f))
            s.(f) = def.(f);
        end
    end
end


function s = mergeStructs_(a,b)
s = a;
fn = fieldnames(b);
for k = 1:numel(fn)
    s.(fn{k}) = b.(fn{k});
end
end

function nv = structToNV_(s)
nv = {};
fn = fieldnames(s);
for k = 1:numel(fn)
    nv(end+1:end+2) = {fn{k}, s.(fn{k})}; %#ok<AGROW>
end
end
