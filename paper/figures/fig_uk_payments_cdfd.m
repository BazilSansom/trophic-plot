function r = fig_uk_payments_cdfd(opts)
%FIG_UK_PAYMENTS_CDFD
% Build UK regional payments CDFD figures for paper.
%
% Produces two separate figures:
%   (1) main:    TFL with directional-flow overlay
%   (2) support: split CDFD figure (directional part + circular part)
%
% Designed to work with make_all_figs / callFig conventions.
%
% Required external functions:
%   cdfd_bff
%   plotCDFD
%
% Example:
%   r = fig_uk_payments_cdfd(struct( ...
%       'Export', true, ...
%       'ExportDir', fullfile(pwd,'outputs','figs')));

    if nargin < 1 || isempty(opts)
        opts = struct();
    end

    % ------------------------------------------------------------
    % defaults
    % ------------------------------------------------------------
    def = struct();

    % export controls
    def.Export        = true;
    def.ExportDir     = '';
    def.ExportBase    = 'fig_uk_payments_cdfd';
    def.ExportBaseOverlay = '';
    def.ExportBaseSplit   = '';
    def.ExportPDF     = true;
    def.ExportPNG     = true;
    def.PNGResolution = 300;
    def.PDFVector     = false;

    % data
    def.XlsxPath      = '';
    def.Sheet         = 'in';
    def.RemoveDiagonal = true;

    % decomposition
    def.DecompositionArgs = struct( ...
        'ToleranceZero', 1e-10, ...
        'Validate', true);

    % labels
    def.LabelColor     = [1 0 0];
    def.LabelFontSize  = 8;
    def.LabelHalo      = true;
    def.LabelHaloAlpha = 0.25;
    def.LabelHaloPadFrac = 0.004;
    def.AdjustLabels   = true;

    % titles
    def.TitleOverlay = 'UK regional payments: TFL with directional-flow overlay';
    def.TitleSplit   = 'UK regional payments: Circular--directional flow decomposition';

    % figure sizes
    def.PositionOverlay = [100 100 1400 1100];
    def.PositionSplit   = [120 120 1700 900];

    % common sizing / widths
    def.NodeSizeRange = [70 260];
    def.EdgeWidthRange = [0.05 6.0];
    def.EdgeWidthQuantiles = [0.20 0.99];
    def.XScale = 1.6;

    % optional user overrides passed through to wrappers
    def.OverlayTFLPlotOpts   = struct();
    def.OverlayTFLLayoutOpts = struct();
    def.SplitTFLPlotOpts     = struct();
    def.SplitTFLLayoutOpts   = struct();
    def.SplitRingPlotOpts    = struct();

    opts = mergeStructsLocal_(def, opts);


    if strlength(string(opts.ExportDir)) == 0
        thisFile = mfilename('fullpath');
        thisDir  = fileparts(thisFile);   % .../trophic-plot/paper/figures
        opts.ExportDir = fullfile(thisDir, '..', 'outputs', 'figs');
    end

    % ------------------------------------------------------------
    % locate workbook if not supplied
    % ------------------------------------------------------------
    opts.XlsxPath = resolveONSPath_(opts.XlsxPath);

    % ------------------------------------------------------------
    % load matrix
    % ------------------------------------------------------------
    [W, labels] = readRegionalPaymentsMatrix_(opts.XlsxPath, opts.Sheet);

    if opts.RemoveDiagonal
        W(1:size(W,1)+1:end) = 0;
    end
    W = sparse(double(W));

    % ------------------------------------------------------------
    % decompose once, then pass through to both figures
    % ------------------------------------------------------------
    cdfdArgs = structToNameValueLocal_(opts.DecompositionArgs);
    [C, D, info] = cdfd_bff(W, cdfdArgs{:});

    uW = full(sum(W,2) + sum(W,1)');
    uC = full(sum(C,2) + sum(C,1)');
    uD = full(sum(D,2) + sum(D,1)');

    labelsForPlot = cellstr(labels);

    % ------------------------------------------------------------
    % export names
    % ------------------------------------------------------------
    if strlength(string(opts.ExportBaseOverlay)) == 0
        exportBaseOverlay = sprintf('%s_overlay', opts.ExportBase);
    else
        exportBaseOverlay = char(opts.ExportBaseOverlay);
    end

    if strlength(string(opts.ExportBaseSplit)) == 0
        exportBaseSplit = sprintf('%s_split', opts.ExportBase);
    else
        exportBaseSplit = char(opts.ExportBaseSplit);
    end

    % ------------------------------------------------------------
    % overlay figure defaults
    % ------------------------------------------------------------
    overlayPlotOpts = struct( ...
        'Labels', {labelsForPlot}, ...
        'ShowLabels', true, ...
        'AdjustLabels', opts.AdjustLabels, ...
        'LabelHalo', opts.LabelHalo, ...
        'LabelHaloAlpha', opts.LabelHaloAlpha, ...
        'LabelHaloPadFrac', opts.LabelHaloPadFrac, ...
        'LabelColor', opts.LabelColor, ...
        'LabelFontSize', opts.LabelFontSize, ...
        'XScale', opts.XScale, ...
        'NodeSizeData', uW, ...
        'NodeSizeScaleMode', 'area', ...
        'NodeSizeRange', opts.NodeSizeRange, ...
        'EdgeWidthMode', 'weight', ...
        'EdgeWidthScale', 'log', ...
        'EdgeWidthRange', opts.EdgeWidthRange, ...
        'EdgeWidthQuantiles', opts.EdgeWidthQuantiles, ...
        'ArrowLenBasis', 'ltyp', ...
        'ArrowSizeMode', 'fixed', ...
        'ArrowSize', 0.06, ...
        'ArrowMinFrac', 0, ...
        'ArrowOverNodes', true, ...
        'ArrowNodeGapPts', 1.0, ...
        'PadFracX', 0.05, ...
        'PadFracY', 0.03, ...
        'CurvPadScaleX', 1.0, ...
        'CurvPadScaleY', 0.12, ...
        'RenderMode', 'overlay', ...
        'OverlayComponent', 'd', ...
        'OverlayBaseColorMode', 'plain', ...
        'ArrowColorMode', 'blend', ...
        'ArrowLayer', 'fullflow', ...
        'LevelSemantics', 'trophic');

    overlayLayoutOpts = struct( ...
        'InitMode', 'spectral', ...
        'InitYNoiseScale', 0.2, ...
        'FreePhaseFrac', 0.2, ...
        'Seed', 1, ...
        'NodeSizeData', uW, ...
        'SizeAwareDeoverlap', true, ...
        'NodeRadiusRange', [0.05 0.14], ...
        'OverlapBandwidthY', 0.6, ...
        'NodeGap', 0.02, ...
        'DeoverlapPasses', 5);

    overlayPlotOpts   = mergeStructsLocal_(overlayPlotOpts,   opts.OverlayTFLPlotOpts);
    overlayLayoutOpts = mergeStructsLocal_(overlayLayoutOpts, opts.OverlayTFLLayoutOpts);

    overlayFigOpts = struct( ...
        'Position', opts.PositionOverlay, ...
        'Color', 'w', ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact', ...
        'PanelLineWidth', 0.75, ...
        'MainTitle', opts.TitleOverlay, ...
        'ShowMainTitle', true, ...
        'ShowPanelTitles', false);

    % ------------------------------------------------------------
    % split figure defaults
    % ------------------------------------------------------------
    splitTFLPlotOpts = struct( ...
        'LevelSemantics', 'trophic', ...
        'Labels', {labelsForPlot}, ...
        'ShowLabels', true, ...
        'AdjustLabels', opts.AdjustLabels, ...
        'LabelHalo', opts.LabelHalo, ...
        'LabelHaloAlpha', opts.LabelHaloAlpha, ...
        'LabelHaloPadFrac', opts.LabelHaloPadFrac, ...
        'LabelColor', opts.LabelColor, ...
        'LabelFontSize', opts.LabelFontSize, ...
        'XScale', opts.XScale, ... %'NodeSizeData', uD, ...
        'NodeSizeScaleMode', 'area', ...
        'NodeSizeRange', opts.NodeSizeRange, ...
        'EdgeWidthMode', 'weight', ...
        'EdgeWidthScale', 'log', ...
        'EdgeWidthRange', opts.EdgeWidthRange, ...
        'EdgeWidthQuantiles', opts.EdgeWidthQuantiles, ...
        'ArrowLenBasis', 'ltyp', ...
        'ArrowSizeMode', 'fixed', ...
        'ArrowSize', 0.10, ...
        'ArrowMinFrac', 0, ...
        'ArrowOverNodes', true, ...
        'ArrowNodeGapPts', 1.0, ...
        'PadFracX', 0.05, ...
        'PadFracY', 0.03, ...
        'CurvPadScaleX', 1.0, ...
        'CurvPadScaleY', 0.12);

    splitTFLLayoutOpts = struct( ...
        'InitMode', 'spectral', ...
        'InitYNoiseScale', 0.2, ...
        'FreePhaseFrac', 0.2, ...
        'Seed', 1, ...
        'NodeSizeData', uW, ...   % full-network sizing signal for layout
        'SizeAwareDeoverlap', true, ...
        'NodeRadiusRange', [0.05 0.14], ...
        'OverlapBandwidthY', 0.6, ...
        'NodeGap', 0.02, ...
        'DeoverlapPasses', 5);

    splitRingPlotOpts = struct( ...
        'Labels', {labelsForPlot}, ...
        'ShowLabels', true, ...
        'LabelColor', opts.LabelColor, ...
        'LabelFontSize', opts.LabelFontSize, ...
        'LabelRadiusFactor', 1.12, ... %'NodeSizeData', uC, ...
        'NodeSizeScaleMode', 'area', ...
        'NodeSizeRange', opts.NodeSizeRange, ...
        'ShowRing', true, ...
        'RingWidth', 0.75, ...
        'EdgeAlpha', 0.75, ...
        'EdgeWidthMode', 'weight', ...
        'EdgeWidthScale', 'log', ...
        'EdgeWidthRange', opts.EdgeWidthRange, ...
        'EdgeWidthQuantiles', opts.EdgeWidthQuantiles, ...
        'EdgeColorMode', 'ringdirection', ...
        'MutualOverlay', true, ...
        'Curvature', 0.78, ...
        'ArrowAngle', pi/7, ...
        'ArrowNodeGapPts', 1.0, ...
        'MutualLaneOffset', 0.5, ...
        'ArrowSize', 0.06, ...
        'ArrowBacktrackFactor', 1.8, ...
        'ArrowPlacement', 'late', ...
        'ArrowPositionFrac', 0.92, ...
        'CurveNpts', 80, ...
        'StartAngle', pi/2, ...
        'Clockwise', false, ...
        'Radius', 1.0, ...
        'Title', '');

    splitTFLPlotOpts   = mergeStructsLocal_(splitTFLPlotOpts,   opts.SplitTFLPlotOpts);
    splitTFLLayoutOpts = mergeStructsLocal_(splitTFLLayoutOpts, opts.SplitTFLLayoutOpts);
    splitRingPlotOpts  = mergeStructsLocal_(splitRingPlotOpts,  opts.SplitRingPlotOpts);

    splitFigOpts = struct( ...
        'Position', opts.PositionSplit, ...
        'Color', 'w', ...
        'TileSpacing', 'compact', ...
        'Padding', 'loose', ...
        'PanelLineWidth', 0.75, ...
        'MainTitle', opts.TitleSplit, ...
        'ShowMainTitle', true, ...
        'ShowPanelTitles', true);

    % ------------------------------------------------------------
    % build overlay figure
    % ------------------------------------------------------------
    outOverlay = plotCDFD(W, ...
        'Mode', 'overlay', ...
        'C', C, 'D', D, 'Info', info, ...
        'Labels', labelsForPlot, ...
        'TFLPlotOpts', overlayPlotOpts, ...
        'TFLLayoutOpts', overlayLayoutOpts, ...
        'FigureOpts', overlayFigOpts);

    % ------------------------------------------------------------
    % build split figure
    % ------------------------------------------------------------
    %{
    outSplit = plotCDFD(W, ...
        'Mode', 'split', ...
        'C', C, 'D', D, 'Info', info, ...
        'Labels', labelsForPlot, ...
        'TFLPlotOpts', splitTFLPlotOpts, ...
        'TFLLayoutOpts', splitTFLLayoutOpts, ...
        'RingPlotOpts', splitRingPlotOpts, ...
        'FigureOpts', splitFigOpts);
        %}
    
    outSplit = plotCDFD(W, ...
        'Mode', 'split', ...
        'SplitLayoutSource', 'w', ...
        'SplitNodeSizeMode', 'component_shared', ...
        'C', C, 'D', D, 'Info', info, ...
        'Labels', labelsForPlot, ...
        'TFLPlotOpts', splitTFLPlotOpts, ...
        'TFLLayoutOpts', splitTFLLayoutOpts, ...
        'RingPlotOpts', splitRingPlotOpts, ...
        'FigureOpts', splitFigOpts);

    % ------------------------------------------------------------
    % export
    % ------------------------------------------------------------
    files = struct();
    if opts.Export
        if strlength(string(opts.ExportDir)) == 0
            error('fig_uk_payments_cdfd:MissingExportDir', ...
                'opts.ExportDir must be supplied when Export=true.');
        end
        if ~exist(opts.ExportDir, 'dir')
            mkdir(opts.ExportDir);
        end

        drawnow;
        %files = exportFig_(outOverlay.fig, opts, exportBaseOverlay, files, 'overlay');
        %$files = exportFig_(outSplit.fig,   opts, exportBaseSplit,   files, 'split');

        files = exportFig_(outOverlay.tiled, opts, exportBaseOverlay, files, 'overlay');
        files = exportFig_(outSplit.tiled,   opts, exportBaseSplit,   files, 'split');

        %files = exportFig_(outOverlay.tiled, ...
        %    opts, exportBaseOverlay, files, 'overlay');

        %files = exportFig_(outSplit.tiled, ...
        %    opts, exportBaseSplit, files, 'split');

    end

    % ------------------------------------------------------------
    % return
    % ------------------------------------------------------------
    r = struct();
    r.fig_overlay = outOverlay.fig;
    r.fig_split   = outSplit.fig;
    r.out_overlay = outOverlay;
    r.out_split   = outSplit;
    r.files       = files;

    r.meta = struct();
    r.meta.XlsxPath = string(opts.XlsxPath);
    r.meta.Sheet    = string(opts.Sheet);
    r.meta.labels   = labels;
    r.meta.W        = W;
    r.meta.C        = C;
    r.meta.D        = D;
    r.meta.info     = info;
    r.meta.uW       = uW;
    r.meta.uC       = uC;
    r.meta.uD       = uD;
end

% =====================================================================
function [W, labels] = readRegionalPaymentsMatrix_(xlsxPath, sheetName)
% ONS workbook convention:
%   rows = destination
%   cols = source
% toolbox convention:
%   rows = source
%   cols = destination
% so transpose after reading.

    rowLabels = string(readcell(xlsxPath, 'Sheet', sheetName, 'Range', 'A9:A20'));
    colLabels = string(readcell(xlsxPath, 'Sheet', sheetName, 'Range', 'B8:M8'));
    rawBody   = readcell(xlsxPath, 'Sheet', sheetName, 'Range', 'B9:M20');

    if ~isequal(size(rawBody), [12 12])
        error('fig_uk_payments_cdfd:BadMatrixSize', ...
            'Expected a 12x12 body in B9:M20, but got %dx%d.', ...
            size(rawBody,1), size(rawBody,2));
    end

    Wraw = nan(size(rawBody));
    for i = 1:size(rawBody,1)
        for j = 1:size(rawBody,2)
            Wraw(i,j) = coerceNumericCell_(rawBody{i,j});
        end
    end

    if any(~isfinite(Wraw(:)))
        error('fig_uk_payments_cdfd:BadMatrixEntries', ...
            'Matrix body contains non-numeric entries after conversion.');
    end

    if ~all(rowLabels(:) == colLabels(:))
        warning('fig_uk_payments_cdfd:LabelMismatch', ...
            'Row and column labels differ; using column labels after transpose.');
    end

    W = Wraw.';
    labels = colLabels(:);
end

% =====================================================================
function x = coerceNumericCell_(v)

    if isnumeric(v)
        if isscalar(v) && isfinite(v)
            x = double(v);
            return;
        end
        x = NaN;
        return;
    end

    if ismissing(v) || isempty(v)
        x = NaN;
        return;
    end

    s = string(v);
    s = strtrim(s);

    if strlength(s) == 0
        x = NaN;
        return;
    end

    % remove commas / currency / spaces that sometimes creep in
    s = replace(s, ",", "");
    s = replace(s, "£", "");
    s = replace(s, " ", "");

    y = str2double(s);
    if isfinite(y)
        x = y;
    else
        x = NaN;
    end
end

% =====================================================================
function p = resolveONSPath_(pIn)

    if strlength(string(pIn)) > 0 && exist(pIn, 'file')
        p = char(pIn);
        return;
    end

    thisFile = mfilename('fullpath');
    thisDir  = fileparts(thisFile);

    cand = {
        fullfile(thisDir, 'data', 'ons_payments', 'ons_regional_payment_flows.xlsx')
        fullfile(thisDir, 'paper', 'data', 'ons_payments', 'ons_regional_payment_flows.xlsx')
        fullfile(fileparts(thisDir), 'paper', 'data', 'ons_payments', 'ons_regional_payment_flows.xlsx')
        fullfile(fileparts(fileparts(thisDir)), 'paper', 'data', 'ons_payments', 'ons_regional_payment_flows.xlsx')
        fullfile(pwd, 'paper', 'data', 'ons_payments', 'ons_regional_payment_flows.xlsx')
        '/Users/bazilsansom/dev/trophic-plot/paper/data/ons_payments/ons_regional_payment_flows.xlsx'
    };

    found = "";
    for ii = 1:numel(cand)
        if exist(cand{ii}, 'file')
            found = string(cand{ii});
            break;
        end
    end

    if strlength(found) == 0
        msg = sprintf(['Could not locate ons_regional_payment_flows.xlsx automatically.\n' ...
            'Tried:\n  %s'], strjoin(cand, newline + "  "));
        error('fig_uk_payments_cdfd:DataFileMissing', '%s', msg);
    end

    p = char(found);
end

% =====================================================================
function files = exportFig_(h, opts, exportBase, files, prefix)

    if opts.ExportPDF
        %pdfPath = fullfile(opts.ExportDir, exportBase + ".pdf");
        pdfPath = fullfile(opts.ExportDir, sprintf('%s.pdf', char(exportBase)));
        if isfield(opts, 'PDFVector') && opts.PDFVector
            exportgraphics(h, pdfPath, 'ContentType', 'vector');
        else
            exportgraphics(h, pdfPath, 'ContentType', 'image');
        end
        files.(sprintf('%s_pdf', prefix)) = string(pdfPath);
    end

    if opts.ExportPNG
        %pngPath = fullfile(opts.ExportDir, exportBase + ".png");
        pngPath = fullfile(opts.ExportDir, sprintf('%s.png', char(exportBase)));
        exportgraphics(h, pngPath, 'Resolution', opts.PNGResolution);
        files.(sprintf('%s_png', prefix)) = string(pngPath);
    end
end

%{
function files = exportFig_(h, opts, exportBase, files, prefix)

    if nargin < 4 || isempty(files)
        files = struct();
    end

    if ~isgraphics(h)
        error('fig_uk_payments_cdfd_dev:BadExportHandle', ...
            'Export target is not a valid graphics object. Class: %s', class(h));
    end

    if opts.ExportPDF
        pdfPath = fullfile(opts.ExportDir, exportBase + ".pdf");
        if opts.PDFVector
            exportgraphics(h, pdfPath, 'ContentType', 'vector');
        else
            exportgraphics(h, pdfPath);
        end
        files.(sprintf('%s_pdf', prefix)) = string(pdfPath);
    end

    if opts.ExportPNG
        pngPath = fullfile(opts.ExportDir, exportBase + ".png");
        exportgraphics(h, pngPath, 'Resolution', opts.PNGResolution);
        files.(sprintf('%s_png', prefix)) = string(pngPath);
    end
end
%}
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
function nv = structToNameValueLocal_(s)
    f = fieldnames(s);
    nv = cell(1, 2*numel(f));
    for ii = 1:numel(f)
        nv{2*ii-1} = f{ii};
        nv{2*ii}   = s.(f{ii});
    end
end

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