function r = fig_translation_force_vs_tfl(opts)
%FIG_TRANSLATION_FORCE_VS_TFL
% Build two separate figures for the translations-network example:
%   (1) force-directed baseline
%   (2) Trophic Flow Layout (TFL)
%
% Designed to work with make_all_figs / callFig conventions.
%
% Required external functions:
%   load_translation_network
%   tflPlot
%
% Example:
%   r = fig_translation_force_vs_tfl(struct( ...
%       'Export', true, ...
%       'ExportDir', fullfile(pwd,'outputs','figs'), ...
%       'ExportBase', 'fig_translation_network', ...
%       'TopN', 60));

    if nargin < 1 || isempty(opts)
        opts = struct();
    end

    % ------------------------------------------------------------
    % defaults
    % ------------------------------------------------------------
    def = struct();

    % export controls
    def.Export        = false;
    def.ExportDir     = '';
    def.ExportBase    = 'fig_translation_network';
    def.ExportBaseForce = '';
    def.ExportBaseTFL   = '';
    def.ExportPDF     = false;
    def.ExportPNG     = false;
    def.PNGResolution = 300;
    def.PDFVector     = false;

    % data / selection
    def.DataDir       = '';
    def.TopN          = 60;
    def.UseSqrtWeights = true;
    def.ForceSeed     = 1;

    % labels
    def.ShowAllLabels = true;
    def.LabelsToShow  = [];   % optional explicit indices into subset
    def.LabelColor    = [1 0 0];
    def.LabelFontSize = 5;
    def.LabelHalo     = true;
    def.LabelHaloAlpha = 0.25;
    def.LabelHaloPadFrac = 0.004;

    % titles
    def.TitleForce = 'Translations network: force-directed';
    def.TitleTFL   = 'Translations network: Trophic Flow Layout';

    % figure sizes
    def.PositionForce = [100 100 900 900];
    def.PositionTFL   = [150 150 1100 800];

    % TFL styling
    def.NodeSizeRangeTFL = [20 220];
    def.XScaleTFL        = 2.2;

    def.TFLPlotOpts = struct( ...
        'ShowLabels', true, ...
        'AdjustLabels', true, ...
        'NodeSizeScaleMode', 'area', ...
        'ArrowLenBasis', 'ltyp', ...
        'ArrowSizeMode', 'fixed', ...
        'ArrowSize', 0.05, ...
        'ArrowMinFrac', 0, ...
        'ArrowOverNodes', true, ...
        'ArrowNodeGapPts', 1);

    def.TFLLayoutOpts = struct( ...
        'InitMode', 'spectral', ...
        'InitYNoiseScale', 0.2, ...
        'FreePhaseFrac', 0.2, ...
        'Seed', 1, ...
        'NodeSizeData', [], ...
        'SizeAwareDeoverlap', true, ...
        'NodeRadiusRange', [0.04 0.14], ...
        'OverlapBandwidthY', 0.6, ...
        'NodeGap', 0.02, ...
        'DeoverlapPasses', 5);

    opts = mergeStructsLocal_(def, opts);

    % ------------------------------------------------------------
    % locate data directory if not supplied
    % ------------------------------------------------------------
    if strlength(string(opts.DataDir)) == 0
    thisFile = mfilename('fullpath');
    thisDir  = fileparts(thisFile);

    cand = {
        fullfile(thisDir, 'data', 'LanguageNets')
        fullfile(fileparts(thisDir), 'data', 'LanguageNets')
        fullfile(fileparts(fileparts(thisDir)), 'data', 'LanguageNets')
        fullfile(fileparts(fileparts(fileparts(thisDir))), 'data', 'LanguageNets')
        fullfile(pwd, 'data', 'LanguageNets')
        fullfile(fileparts(pwd), 'data', 'LanguageNets')
    };

    found = "";
    for ii = 1:numel(cand)
        if exist(cand{ii}, 'dir')
            found = string(cand{ii});
            break;
        end
    end

    if strlength(found) == 0
        msg = sprintf(['Could not locate data/LanguageNets automatically.\n' ...
            'Tried:\n  %s'], strjoin(cand, newline + "  "));
        error('fig_translation_force_vs_tfl:DataDirMissing', '%s', msg);
    end

    opts.DataDir = char(found);
end

if ~exist(opts.DataDir, 'dir')
    error('fig_translation_force_vs_tfl:DataDirMissing', ...
        'DataDir not found: %s', opts.DataDir);
end

%{
    
    if strlength(string(opts.DataDir)) == 0
        thisFile = mfilename('fullpath');
        thisDir  = fileparts(thisFile);

        % assume function is in repoRoot/src/
        [parentDir, thisLeaf] = fileparts(thisDir);
        if strcmpi(thisLeaf, 'src')
            repoRoot = parentDir;
        else
            repoRoot = thisDir;
        end
        opts.DataDir = fullfile(repoRoot, 'data', 'LanguageNets');
    end

    if ~exist(opts.DataDir, 'dir')
        error('fig_translation_force_vs_tfl:DataDirMissing', ...
            'DataDir not found: %s', opts.DataDir);
    end

    %}

    % ------------------------------------------------------------
    % load data
    % ------------------------------------------------------------
    net = load_translation_network(opts.DataDir);

    [~, ord] = sort(net.totalWeight, 'descend');
    keep = ord(1:min(opts.TopN, numel(ord)));

    Wsub      = net.W(keep, keep);
    labelsSub = string(net.labels(keep));
    uSub      = net.totalWeight(keep);

    if opts.UseSqrtWeights
        WsubPlot = sqrt(Wsub);
    else
        WsubPlot = Wsub;
    end
    WsubPlot = sparse(WsubPlot);

    % ------------------------------------------------------------
    % labels to show
    % ------------------------------------------------------------
    if opts.ShowAllLabels
        labelsForPlot = cellstr(labelsSub);
    else
        labelsForPlot = repmat({''}, numel(labelsSub), 1);
        if isempty(opts.LabelsToShow)
            % default fallback: top 20 by weight
            [~, ordLab] = sort(uSub, 'descend');
            showIdx = ordLab(1:min(20, numel(ordLab)));
        else
            showIdx = opts.LabelsToShow(:);
        end
        labelsForPlot(showIdx) = cellstr(labelsSub(showIdx));
    end

    % ------------------------------------------------------------
    % common node-size scaling
    % ------------------------------------------------------------
    % TFL uses scatter marker area units
    nodeSizeRangeTFL = opts.NodeSizeRangeTFL;

    % digraph plot uses marker sizes in points
    uScaled = sqrt(uSub / max(uSub));
    markerSizeForce = 3 + 7 * uScaled;

    % ------------------------------------------------------------
    % deterministic force-directed coordinates
    % ------------------------------------------------------------
    Gsub = digraph(WsubPlot, cellstr(labelsSub));

    rng(opts.ForceSeed);

    tmpFig = figure('Visible', 'off', 'Color', 'w');
    tmpAx  = axes(tmpFig);
    pTmp   = plot(tmpAx, Gsub, 'Layout', 'force', 'NodeLabel', {});
    xForce = pTmp.XData(:);
    yForce = pTmp.YData(:);
    close(tmpFig);

    % ------------------------------------------------------------
    % TFL options
    % ------------------------------------------------------------
    plotOptsTFL = mergeStructsLocal_(def.TFLPlotOpts, opts.TFLPlotOpts);
    layoutOptsTFL = mergeStructsLocal_(def.TFLLayoutOpts, opts.TFLLayoutOpts);

    plotOptsTFL.Labels = labelsForPlot;
    plotOptsTFL.ShowLabels = true;
    plotOptsTFL.LabelHalo = opts.LabelHalo;
    plotOptsTFL.LabelHaloAlpha = opts.LabelHaloAlpha;
    plotOptsTFL.LabelHaloPadFrac = opts.LabelHaloPadFrac;
    plotOptsTFL.LabelColor = opts.LabelColor;
    plotOptsTFL.LabelFontSize = opts.LabelFontSize;
    plotOptsTFL.XScale = opts.XScaleTFL;
    plotOptsTFL.NodeSizeData = uSub;
    plotOptsTFL.NodeSizeRange = nodeSizeRangeTFL;

    layoutOptsTFL.NodeSizeData = uSub;

    % ------------------------------------------------------------
    % export names
    % ------------------------------------------------------------
    if strlength(string(opts.ExportBaseForce)) == 0
        exportBaseForce = sprintf('%s_force_directed', opts.ExportBase);
    else
        exportBaseForce = char(opts.ExportBaseForce);
    end

    if strlength(string(opts.ExportBaseTFL)) == 0
        exportBaseTFL = sprintf('%s_tfl', opts.ExportBase);
    else
        exportBaseTFL = char(opts.ExportBaseTFL);
    end

    % ------------------------------------------------------------
    % figure 1: force-directed
    % ------------------------------------------------------------
    figForce = figure('Color', 'w', 'Position', opts.PositionForce);
    axForce = axes(figForce);
    set(axForce, 'Color', 'w');

    p = plot(axForce, Gsub, ...
        'XData', xForce, ...
        'YData', yForce, ...
        'NodeLabel', labelsForPlot);

    % style
    p.NodeColor = [0.1 0.3 0.7];
    p.EdgeColor = [0.45 0.45 0.45];
    p.LineWidth = 0.5;
    p.ArrowSize = 8;

    p.NodeLabelColor = opts.LabelColor;
    p.NodeFontSize   = opts.LabelFontSize;

    try
        p.MarkerSize = markerSizeForce;
    catch
        p.MarkerSize = 4;
        for ii = 1:numel(labelsSub)
            highlight(p, ii, 'MarkerSize', markerSizeForce(ii));
        end
    end

    % frame / axes
    xpad = 0.08 * range(xForce);
    ypad = 0.08 * range(yForce);
    if xpad == 0, xpad = 1; end
    if ypad == 0, ypad = 1; end

    xlim(axForce, [min(xForce) - xpad, max(xForce) + xpad]);
    ylim(axForce, [min(yForce) - ypad, max(yForce) + ypad]);

    axis(axForce, 'equal');
    pbaspect(axForce, [1 1 1]);
    set(axForce, 'XTick', [], 'YTick', [], 'XTickLabel', [], 'YTickLabel', []);
    box(axForce, 'on');
    set(axForce, 'LineWidth', 0.75);

    title(axForce, opts.TitleForce, 'FontWeight', 'bold', 'Interpreter', 'none');

    % ------------------------------------------------------------
    % figure 2: TFL
    % ------------------------------------------------------------
    figTFL = figure('Color', 'w', 'Position', opts.PositionTFL);
    axTFL = axes(figTFL);
    set(axTFL, 'Color', 'w');
    axes(axTFL); %#ok<LAXES>

    tflPlot(WsubPlot, ...
        'PlotOpts', plotOptsTFL, ...
        'LayoutOpts', layoutOptsTFL);

    set(axTFL, 'Color', 'w');
    set(axTFL, 'LineWidth', 0.75);
    title(axTFL, opts.TitleTFL, 'FontWeight', 'bold', 'Interpreter', 'none');

    % ------------------------------------------------------------
    % export
    % ------------------------------------------------------------
    files = struct();
    if opts.Export
        if strlength(string(opts.ExportDir)) == 0
            error('fig_translation_force_vs_tfl:MissingExportDir', ...
                'opts.ExportDir must be supplied when Export=true.');
        end
        if ~exist(opts.ExportDir, 'dir')
            mkdir(opts.ExportDir);
        end

        files = struct();

        drawnow;
        files = exportFig_(figForce, opts, exportBaseForce, files, 'force');
        files = exportFig_(figTFL,   opts, exportBaseTFL,   files, 'tfl');
    end

    % ------------------------------------------------------------
    % return
    % ------------------------------------------------------------
    r = struct();
    r.fig_force = figForce;
    r.fig_tfl   = figTFL;
    r.files     = files;

    r.meta = struct();
    r.meta.keep      = keep;
    r.meta.labelsSub = labelsSub;
    r.meta.uSub      = uSub;
    r.meta.WsubPlot  = WsubPlot;
    r.meta.xForce    = xForce;
    r.meta.yForce    = yForce;
end

% =====================================================================
function files = exportFig_(fig, opts, exportBase, files, prefix)

    if opts.ExportPDF
        pdfPath = fullfile(opts.ExportDir, sprintf('%s.pdf', char(exportBase)));
        if opts.PDFVector
            exportgraphics(fig, pdfPath, ...
                'ContentType', 'vector', ...
                'BackgroundColor', 'white');
        else
            exportgraphics(fig, pdfPath, ...
                'ContentType', 'image', ...
                'Resolution', opts.PNGResolution, ...
                'BackgroundColor', 'white');
        end
        files.(sprintf('%s_pdf', prefix)) = string(pdfPath);
    end

    if opts.ExportPNG
        pngPath = fullfile(opts.ExportDir, sprintf('%s.png', char(exportBase)));
        exportgraphics(fig, pngPath, ...
            'Resolution', opts.PNGResolution, ...
            'BackgroundColor', 'white');
        files.(sprintf('%s_png', prefix)) = string(pngPath);
    end
end

%{
function files = exportFig_(fig, opts, exportBase, files, prefix)
    if opts.ExportPDF
        pdfPath = fullfile(opts.ExportDir, exportBase + ".pdf");
        if opts.PDFVector
            exportgraphics(fig, pdfPath, 'ContentType', 'vector');
        else
            exportgraphics(fig, pdfPath);
        end
        files.(sprintf('%s_pdf', prefix)) = string(pdfPath);
    end

    if opts.ExportPNG
        pngPath = fullfile(opts.ExportDir, exportBase + ".png");
        exportgraphics(fig, pngPath, 'Resolution', opts.PNGResolution);
        files.(sprintf('%s_png', prefix)) = string(pngPath);
    end
end
%}

% =====================================================================
function s = mergeStructsLocal_(a, b)
% Merge scalar structs recursively: fields in b overwrite fields in a.
    s = a;
    if isempty(b)
        return;
    end
    fn = fieldnames(b);
    for i = 1:numel(fn)
        k = fn{i};
        if isfield(s, k) && isstruct(s.(k)) && isscalar(s.(k)) ...
                         && isstruct(b.(k)) && isscalar(b.(k))
            s.(k) = mergeStructsLocal_(s.(k), b.(k));
        else
            s.(k) = b.(k);
        end
    end
end