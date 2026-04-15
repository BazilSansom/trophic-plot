function out = make_all_figs(opts)
%MAKE_ALL_FIGS  Build all paper figures reproducibly into paper/outputs/figs.
%
% out = make_all_figs()
% out = make_all_figs(struct('ExportPNG',true,'PNGResolution',300))
%
% Assumes figure generators are callable functions and support opts fields:
%   Export, ExportDir, ExportBase, ExportPDF, ExportPNG, PNGResolution,
%   PDFVector
%
% Returns:
%   out.exportDir
%   out.results   : struct of per-figure outputs
%   out.files     : struct array listing produced files

    if nargin < 1 || isempty(opts), opts = struct(); end

    % ------------------ defaults ------------------
    def.ExportPDF     = true;
    def.ExportPNG     = true;
    def.PNGResolution = 300;
    def.PDFVector     = true;
    def.CloseFigs     = true;

    opts = applyDefaults_(opts, def);

    % ------------------ paths ------------------
    thisFile = mfilename('fullpath');
    thisDir  = fileparts(thisFile);

    % Preferred layout:
    %   repoRoot/
    %       toolbox/src/   or src/
    %       paper/
    %           build/
    %               make_all_figs.m
    paperDir = fileparts(thisDir);
    repoRoot = fileparts(paperDir);

    toolboxSrcDir = fullfile(repoRoot, 'toolbox', 'src');
    plainSrcDir   = fullfile(repoRoot, 'src');

    % Fallback if make_all_figs.m is stored directly under paper/
    if ~exist(toolboxSrcDir, 'dir') && ~exist(plainSrcDir, 'dir')
        paperDir = thisDir;
        repoRoot = fileparts(paperDir);

        toolboxSrcDir = fullfile(repoRoot, 'toolbox', 'src');
        plainSrcDir   = fullfile(repoRoot, 'src');
    end

    if exist(toolboxSrcDir, 'dir')
        srcDir = toolboxSrcDir;
    elseif exist(plainSrcDir, 'dir')
        srcDir = plainSrcDir;
    else
        error('make_all_figs:MissingSrcDir', ...
            'Could not find source directory under toolbox/src or src.');
    end

    exportDir = fullfile(paperDir, 'outputs', 'figs');
    cacheDir  = fullfile(repoRoot, 'paper', 'data', 'mangal_examples');

    if ~exist(exportDir, 'dir'), mkdir(exportDir); end
    if ~exist(cacheDir,  'dir'), mkdir(cacheDir);  end

    addProjectTreeSansScratch_(srcDir);
    addProjectTreeSansScratch_(paperDir);

    %{
    % ------------------ paths ------------------
    thisFile = mfilename('fullpath');
    thisDir  = fileparts(thisFile);

    % Preferred layout:
    %   repoRoot/
    %       src/
    %       paper/
    %           build/
    %               make_all_figs.m
    paperDir = fileparts(thisDir);
    repoRoot = fileparts(paperDir);

    % Fallback if make_all_figs.m is stored directly under paper/
    if ~exist(fullfile(repoRoot, 'src'), 'dir')
        paperDir = thisDir;
        repoRoot = fileparts(paperDir);
    end

    exportDir = fullfile(paperDir, 'outputs', 'figs');
    cacheDir  = fullfile(repoRoot, 'paper', 'data', 'mangal_examples');

    toolboxSrcDir = fullfile(repoRoot, 'toolbox', 'src');
    plainSrcDir   = fullfile(repoRoot, 'src');

    if exist(toolboxSrcDir, 'dir')
        srcDir = toolboxSrcDir;
    elseif exist(plainSrcDir, 'dir')
        srcDir = plainSrcDir;
    else
        error('make_all_figs:MissingSrcDir', ...
            'Could not find source directory under toolbox/src or src.');
    end

    %srcDir    = fullfile(repoRoot, 'src');
    %exportDir = fullfile(paperDir, 'outputs', 'figs');
    %cacheDir  = fullfile(repoRoot, 'paper', 'data', 'mangal_examples');

    %srcDir    = fullfile(repoRoot, 'src');
    %exportDir = fullfile(paperDir, 'outputs', 'figs');
    %cacheDir  = fullfile(repoRoot, 'data', 'mangal_examples');

    if ~exist(exportDir, 'dir'), mkdir(exportDir); end
    if ~exist(cacheDir,  'dir'), mkdir(cacheDir);  end

    %addpath(genpath(srcDir));
    %addpath(genpath(paperDir));
    addProjectTreeSansScratch_(srcDir);
    addProjectTreeSansScratch_(paperDir);
    %}


    % ------------------ build all figures ------------------
    results = struct();
    files   = struct('figure', {}, 'type', {}, 'path', {});

    % Helper: collect any nonempty file paths returned by a figure function
        function collectFiles_(r, baseName)
        if ~isfield(r, 'files') || ~isstruct(r.files)
            return;
        end

        fns = fieldnames(r.files);
        for ii = 1:numel(fns)
            f = fns{ii};
            raw = r.files.(f);

            if isempty(raw)
                continue;
            end

            vals = string(raw);
            vals = vals(:);
            vals = vals(strlength(vals) > 0);

            for jj = 1:numel(vals)
                files(end+1,1) = struct( ... %#ok<AGROW>
                    'figure', string(baseName), ...
                    'type',   string(f), ...
                    'path',   vals(jj));
            end
        end
    end
    %{
    function collectFiles_(r, baseName)
        if ~isfield(r, 'files') || ~isstruct(r.files)
            return;
        end
        fns = fieldnames(r.files);
        for ii = 1:numel(fns)
            f = fns{ii};
            p = string(r.files.(f));
            if strlength(p) > 0
                files(end+1,1) = struct( ... %#ok<AGROW>
                    'figure', string(baseName), ...
                    'type',   string(f), ...
                    'path',   p);
            end
        end
    end
    %}

    
    % Helper: standard opts + figure function of form r = figFn(opts)
        function r = callFig(figFn, baseName, extraOpts)
        if nargin < 3, extraOpts = struct(); end
        o = mergeStructs_(standardFigOpts_(baseName), extraOpts);

        figsBefore = findall(groot, 'Type', 'figure');

        r = struct( ...
            'ok', false, ...
            'figure', string(baseName), ...
            'generator', string(figFn), ...
            'files', struct(), ...
            'errorIdentifier', "", ...
            'errorMessage', "");

        if exist(figFn, 'file') ~= 2
            r.errorIdentifier = "make_all_figs:MissingFigureFunction";
            r.errorMessage = "Figure function not found: " + string(figFn);
            warning('%s', r.errorMessage);
            return;
        end

        try
            r0 = feval(figFn, o);
            if isempty(r0)
                r0 = struct();
            end
            if ~isstruct(r0)
                error('make_all_figs:BadFigureReturn', ...
                    'Figure function %s must return a struct.', figFn);
            end

            r = mergeStructs_(r0, r);
            r.ok = true;

            collectFiles_(r, baseName);

            if opts.CloseFigs
                closeReturnedFigs_(r);
                closeNewFigures_(figsBefore);
            end

        catch ME
            r.ok = false;
            r.errorIdentifier = string(ME.identifier);
            r.errorMessage    = string(ME.message);

            warning('Figure build failed for %s (%s): %s', ...
                baseName, figFn, ME.message);

            if opts.CloseFigs
                closeNewFigures_(figsBefore);
            end
        end
    end

    %{
    function r = callFig(figFn, baseName, extraOpts)
        if nargin < 3, extraOpts = struct(); end
        o = mergeStructs_(standardFigOpts_(baseName), extraOpts);

        r = feval(figFn, o);

        collectFiles_(r, baseName);

        if opts.CloseFigs
            closeReturnedFigs_(r);
        end
    end

        function closeReturnedFigs_(r)
        closeFigureHandlesInStruct_(r);
    end

    %}

    function closeFigureHandlesInStruct_(s)
        if ~isstruct(s), return; end

        fn = fieldnames(s);
        for ii = 1:numel(fn)
            v = s.(fn{ii});

            if ishghandle(v)
                try
                    if strcmp(get(v, 'Type'), 'figure')
                        close(v);
                    end
                catch
                end

            elseif isstruct(v)
                closeFigureHandlesInStruct_(v);

            elseif iscell(v)
                for jj = 1:numel(v)
                    tryCloseValue_(v{jj});
                end
            end
        end
    end

    function tryCloseValue_(v)
        if ishghandle(v)
            try
                if strcmp(get(v, 'Type'), 'figure')
                    close(v);
                end
            catch
            end
        elseif isstruct(v)
            closeFigureHandlesInStruct_(v);
        elseif iscell(v)
            for kk = 1:numel(v)
                tryCloseValue_(v{kk});
            end
        end
    end

    %{
    function closeReturnedFigs_(r)
        figFields = {'fig', 'fig_force', 'fig_tfl', 'fig_networks', ...
             'fig_diagnostics', 'fig_layered', 'fig_snapped'};
        %figFields = {'fig', 'fig_networks', 'fig_diagnostics', 'fig_layered', 'fig_snapped'};
        for jj = 1:numel(figFields)
            ff = figFields{jj};
            if isfield(r, ff) && ishghandle(r.(ff))
                close(r.(ff));
            end
        end
    end
    %}
    
    % Helper: standard opts + figure function with positional args:
    %   r = figFn(arg1, arg2, ..., opts)

        function r = callFigWithArgs(figFn, figArgs, baseName, extraOpts)
        if nargin < 4, extraOpts = struct(); end
        o = mergeStructs_(standardFigOpts_(baseName), extraOpts);

        figsBefore = findall(groot, 'Type', 'figure');

        r = struct( ...
            'ok', false, ...
            'figure', string(baseName), ...
            'generator', string(figFn), ...
            'files', struct(), ...
            'errorIdentifier', "", ...
            'errorMessage', "");

        if exist(figFn, 'file') ~= 2
            r.errorIdentifier = "make_all_figs:MissingFigureFunction";
            r.errorMessage = "Figure function not found: " + string(figFn);
            warning('%s', r.errorMessage);
            return;
        end

        try
            r0 = feval(figFn, figArgs{:}, o);
            if isempty(r0)
                r0 = struct();
            end
            if ~isstruct(r0)
                error('make_all_figs:BadFigureReturn', ...
                    'Figure function %s must return a struct.', figFn);
            end

            r = mergeStructs_(r0, r);
            r.ok = true;

            collectFiles_(r, baseName);

            if opts.CloseFigs
                closeReturnedFigs_(r);
                closeNewFigures_(figsBefore);
            end

        catch ME
            r.ok = false;
            r.errorIdentifier = string(ME.identifier);
            r.errorMessage    = string(ME.message);

            warning('Figure build failed for %s (%s): %s', ...
                baseName, figFn, ME.message);

            if opts.CloseFigs
                closeNewFigures_(figsBefore);
            end
        end
    end

    %{
    function r = callFigWithArgs(figFn, figArgs, baseName, extraOpts)
        if nargin < 4, extraOpts = struct(); end
        o = mergeStructs_(standardFigOpts_(baseName), extraOpts);

        r = feval(figFn, figArgs{:}, o);

        collectFiles_(r, baseName);

        if opts.CloseFigs
            closeReturnedFigs_(r);
        end
    end
    %}


    function o = standardFigOpts_(baseName)
        o = struct();
        o.Export        = true;
        o.ExportDir     = exportDir;
        o.ExportBase    = baseName;
        o.ExportPDF     = opts.ExportPDF;
        o.ExportPNG     = opts.ExportPNG;
        o.PNGResolution = opts.PNGResolution;
        o.PDFVector     = opts.PDFVector;
    end

    % ------------------ shared style ------------------
    paperPlotStyle = struct( ...
        'ArrowSize',     0.10, ...
        'ArrowSizeMode', 'fixed', ...
        'ArrowMinFrac',  0.03, ...
        'ArrowMaxFrac',  0 ...
    );

    earlyVisualStyle = struct( ...
        'CacheDir',          cacheDir, ...
        'ForceSeed',         1, ...
        'ShowLabels',        true, ...
        'LabelFontSize',     5, ...
        'LabelColor',        [0.6350, 0.0780, 0.1840], ...
        'LabelHalo',         true, ...
        'LabelHaloPadFrac',  0.006, ...
        'LabelHaloAlpha',    0.30, ...
        'LabelHaloColor',    [1 1 1], ...
        'AdjustLabels',      false, ...
        'ArrowSize',         0.03, ...
        'ArrowSizeMode',     'fixed', ...
        'ArrowMinFrac',      0.03, ...
        'ArrowMaxFrac',      0, ...
        'EdgeWidthMode',     'fixed' ...
    );

    % ------------------ figure list ------------------
    results.fig_cyclic_comparison = callFig( ...
        'fig_cyclic_comparison', ...
        'fig_cyclic_comparison');

    results.fig_general_dag_comparison = callFig( ...
        'fig_general_dag_comparison', ...
        'fig_general_dag_comparison');

    results.fig_topology_vs_weighted_flow = callFig( ...
        'fig_topology_vs_weighted_flow', ...
        'fig_topology_vs_weighted_flow');

    results.fig_perfectly_layered_comparison = callFig( ...
        'fig_perfectly_layered_comparison', ...
        'fig_perfectly_layered_comparison');

    results.fig_ordinal_metric_scatter = callFig( ...
        'fig_scatter_tfl_vs_layered', ...
        'fig_ordinal_metric_scatter', ...
        struct( ...
            'generator',      'gppm', ...
            'n',              60, ...
            'p',              0.02, ...
            'T',              10, ...
            'gppmArgs',       {{'ForceDAG', true, 'Weighted', false}}, ...
            'showLayouts',    true, ...
            'seed',           13, ...
            'visible',        false, ...
            'PaperPlotStyle', paperPlotStyle ...
        ) ...
    );

    % Early visual empirical food-web comparison (Mangal)
    DATASET_ID = 74;
    NETWORK_ID = 1068;

    results.fig_early_visual = callFigWithArgs( ...
        'fig_early_picture_foodweb_from_mangal', ...
        {DATASET_ID, NETWORK_ID}, ...
        'fig_early_visual', ...
        earlyVisualStyle ...
    );

    % Empirical food web (FW24): TFL vs Sugiyama, and TFL vs snapped
    results.fig_foodweb_FW24 = callFig( ...
        'fig_foodweb_FW24', ...
        'fig_foodweb_FW24', ...
        struct( ...
            'ShowLabels', true, ...
            'LabelFontSize', 6, ...
            'LabelHalo', true, ...
            'LabelHaloPadFrac', 0.006, ...
            'LabelHaloAlpha', 0.35, ...
            'LabelHaloColor', [1 1 1], ...
            'LabelColor',       [0.6350, 0.0780, 0.1840],...%[0, 0, 1], ...  
            'AdjustLabels', true, ...
            'ExportBaseLayered', 'fig_foodweb_FW24_tfl_vs_sugiyama', ...
            'ExportBaseSnapped', 'fig_foodweb_FW24_tfl_vs_snapped' ...
        ) ...
    );

    % Bike sharing TfL Santander Cycles data

    % TfL Santander Cycles: morning/evening comparison
    results.fig_tfl_santander_morning_evening_compare = callFig( ...
        'fig_tfl_santander_morning_evening_compare', ...
        'fig_tfl_santander_morning_evening_compare', ...
        struct( ...
            'TargetDate', datetime(2024,6,4), ...
            'MorningWindow', [7 10], ...
            'EveningWindow', [17 20], ...
            'K', 30, ...
            'ExportBaseNetworks', 'fig_tfl_santander_morning_evening_networks', ...
            'ExportBaseDiagnostics', 'fig_tfl_santander_morning_evening_diagnostics' ...
        ) ...
    );

    % TfL Santander Cycles: focal corridor comparison
    results.fig_tfl_santander_focal_corridor_compare = callFig( ...
        'fig_tfl_santander_focal_corridor_compare', ...
        'fig_tfl_santander_focal_corridor_compare', ...
        struct( ...
            'TargetDate', datetime(2024,6,4), ...     % Fixed weekday snapshot for reproducible paper figures
            'MorningWindow', [7 10], ...
            'EveningWindow', [17 20], ...
            'K', 30, ...
            'FocalStationName', "Waterloo Station 3, Waterloo", ...
            'MaxHops', 2 ...
        ) ...
    );

    % Transaltion network: Force-directed vs TFL
    results.fig_translation_force_vs_tfl = callFig( ...
        'fig_translation_force_vs_tfl', ...
        'fig_translation_force_vs_tfl', ...
        struct( ...
            'DataDir', fullfile(repoRoot, 'paper', 'data', 'LanguageNets'), ... % 'DataDir', fullfile(repoRoot, 'data', 'LanguageNets'), ...
            'TopN', 60, ...
            'ForceSeed', 1, ...
            'ShowAllLabels', true, ...
            'ExportBaseForce', 'fig_translation_force_directed', ...
            'ExportBaseTFL',   'fig_translation_tfl' ...
        ) ...
    );

    % UK payments network with CDFD overlay and split
    results.fig_uk_payments_cdfd = callFig( ...
        'fig_uk_payments_cdfd', ...
        'fig_uk_payments_cdfd', ...
        struct( ...
            'ExportBaseOverlay', 'fig_uk_payments_overlay', ...
            'ExportBaseSplit',   'fig_uk_payments_split' ...
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
    if ~isfield(s, f) || isempty(s.(f))
        s.(f) = def.(f);
    end
end
end

function s = mergeStructs_(a, b)
    s = a;
    fn = fieldnames(b);
    for k = 1:numel(fn)
        s.(fn{k}) = b.(fn{k});
    end
end

function closeNewFigures_(figsBefore)
    figsAfter = findall(groot, 'Type', 'figure');
    figsNew = setdiff(figsAfter, figsBefore);
    for ii = 1:numel(figsNew)
        if ishghandle(figsNew(ii))
            try
                close(figsNew(ii));
            catch
            end
        end
    end
end

function addProjectTreeSansScratch_(rootDir)
% Add rootDir recursively, excluding scratch-like folders.

    p = genpath(rootDir);
    parts = strsplit(p, pathsep);

    keep = true(size(parts));
    for i = 1:numel(parts)
        pi = string(parts{i});
        if strlength(pi) == 0
            keep(i) = false;
            continue;
        end

        low = lower(pi);
        if contains(low, filesep + "scratch") || ...
           contains(low, filesep + ".git") || ...
           contains(low, filesep + ".github")
            keep(i) = false;
        end
    end

    parts = parts(keep);
    if ~isempty(parts)
        addpath(strjoin(parts, pathsep));
    end
end