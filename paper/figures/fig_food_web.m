function out = fig_food_web(E, opts)
%FIG_FOOD_WEB  Build the empirical food-web comparison figures used in paper.
%
% Produces two separate 1x2 figures:
%
%   Figure A: TFL vs Sugiyama layered
%       (1) Trophic Flow Layout (TFL)
%       (2) Sugiyama layered
%
%   Figure B: TFL vs TFL snapped
%       (1) Trophic Flow Layout (TFL)
%       (2) Trophic Flow Layout (TFL), snapped heights
%
% Inputs:
%   E    : edge list [source target weight] or format accepted by
%          edgelistToAdjacency(E)
%   opts : struct of options
%
% Label options:
%   opts.ShowLabels            : true/false (default true)
%   opts.Labels                : cellstr / string array of node labels.
%                                If empty, plotting falls back to numeric labels.
%
% Export options:
%   opts.Export                : true/false (default false)
%   opts.ExportDir             : output folder (default pwd)
%   opts.ExportBaseLayered     : basename for TFL-vs-Sugiyama figure
%                                (default 'fig_foodweb_tfl_vs_sugiyama')
%   opts.ExportBaseSnapped     : basename for TFL-vs-snapped figure
%                                (default 'fig_foodweb_tfl_vs_snapped')
%   opts.ExportPDF             : true/false (default true)
%   opts.ExportPNG             : true/false (default true)
%   opts.PNGResolution         : dpi (default 300)
%   opts.PDFVector             : true/false (default true)
%
% Plot options:
%   opts.ArrowSize             : default 0.03
%   opts.PadFrac               : default 0.02
%   opts.TileSpacing           : default 'compact'
%   opts.Padding               : default 'compact'
%   opts.FigWidthCM            : default 16
%   opts.FigHeightCM           : default 8.5
%   opts.SnappedMinSepFrac     : default 0.05
%
% Outputs:
%   out.fig_layered
%   out.fig_snapped
%   out.axes_layered
%   out.axes_snapped
%   out.files.layered_pdf
%   out.files.layered_png
%   out.files.snapped_pdf
%   out.files.snapped_png
%   out.W
%   out.labels

    if nargin < 2 || isempty(opts), opts = struct(); end
    opts = applyDefaults_(opts);

    S = paperStyle();

    % ---------------------------------------------------------------------
    % Data
    % ---------------------------------------------------------------------
    W = edgelistToAdjacency(E);
    n = size(W,1);

    labels = prepareLabels_(opts.Labels, n);
        
    plotOptsBase = struct( ...
        'ShowLabels',     opts.ShowLabels, ...
        'Labels',         {labels}, ...
        'LabelFontSize',  opts.LabelFontSize, ...
        'LabelColor',     opts.LabelColor, ...
        'LabelHalo',      opts.LabelHalo, ...
        'LabelHaloPadFrac', opts.LabelHaloPadFrac, ...
        'LabelHaloAlpha', opts.LabelHaloAlpha, ...
        'LabelHaloColor', opts.LabelHaloColor, ...
        'AdjustLabels',    opts.AdjustLabels, ...
        'EdgeWidthMode',  'fixed', ...
        'EdgeWidthScale', 'linear', ...
        'ArrowSize',      opts.ArrowSize);
    
    % ---------------------------------------------------------------------
    % Figure 1: TFL vs Sugiyama layered
    % ---------------------------------------------------------------------
    fig1 = figure('Color','w');
    tl1 = tiledlayout(fig1, 1, 2, ...
        'TileSpacing', opts.TileSpacing, ...
        'Padding',     opts.Padding);

    set(fig1, 'Units', 'centimeters');
    fig1.Position(3) = opts.FigWidthCM;
    fig1.Position(4) = opts.FigHeightCM;

    axs1 = gobjects(2,1);

    % Panel 1: TFL
    axs1(1) = nexttile(tl1, 1);
    tflPlot(W, 'PlotOpts', plotOptsBase);

    setPanelTitle(axs1(1), "Trophic Flow Layout (TFL)", S);
    fixTiledAxesCommon(axs1(1));
    makeAxesDataSquare(axs1(1), opts.PadFrac);

    % Panel 2: Sugiyama layered
    axs1(2) = nexttile(tl1, 2);
    g = digraph(W);
    p = plot(axs1(2), g, 'Layout', 'layered', 'Direction', 'up');
    Xg = p.XData(:);
    Yg = p.YData(:);
    delete(p);

    plotTFL(W, Xg, Yg, Yg, ...
        'Parent',          axs1(2), ...
        'ClearAxes',       true, ...
        'ShowBands',       false, ...
        'ShowLabels',      opts.ShowLabels, ...
        'Labels',          labels, ...
        'LabelFontSize',   opts.LabelFontSize, ...
        'LabelColor',      opts.LabelColor, ...
        'LabelHalo',       opts.LabelHalo, ...
        'AdjustLabels',    opts.AdjustLabels, ...
        'LabelHaloPadFrac',opts.LabelHaloPadFrac, ...
        'LabelHaloAlpha',  opts.LabelHaloAlpha, ...
        'LabelHaloColor',  opts.LabelHaloColor, ...
        'ArrowSize',       opts.ArrowSize);

    %{
    plotTFL(W, Xg, Yg, Yg, ...
        'Parent',     axs1(2), ...
        'ClearAxes',  true, ...
        'ShowBands',  false, ...
        'ShowLabels', opts.ShowLabels, ...
        'Labels',     labels, ...
        'ArrowSize',  opts.ArrowSize);
%}
    setPanelTitle(axs1(2), "Sugiyama layered", S);
    fixTiledAxesCommon(axs1(2));
    makeAxesDataSquare(axs1(2), opts.PadFrac);

    drawnow;

    % ---------------------------------------------------------------------
    % Figure 2: TFL vs TFL snapped
    % ---------------------------------------------------------------------
    fig2 = figure('Color','w');
    tl2 = tiledlayout(fig2, 1, 2, ...
        'TileSpacing', opts.TileSpacing, ...
        'Padding',     opts.Padding);

    set(fig2, 'Units', 'centimeters');
    fig2.Position(3) = opts.FigWidthCM;
    fig2.Position(4) = opts.FigHeightCM;

    axs2 = gobjects(2,1);

    % Panel 1: TFL
    axs2(1) = nexttile(tl2, 1);
    tflPlot(W, 'PlotOpts', plotOptsBase);

    setPanelTitle(axs2(1), "Trophic Flow Layout (TFL)", S);
    fixTiledAxesCommon(axs2(1));
    makeAxesDataSquare(axs2(1), opts.PadFrac);

    % Panel 2: TFL snapped
    axs2(2) = nexttile(tl2, 2);
    tflPlot(W, ...
        'PlotOpts', plotOptsBase, ...
        'LayoutOpts', struct( ...
            'HeightMode', 'snapped', ...
            'MinSepFrac', opts.SnappedMinSepFrac));

    setPanelTitle(axs2(2), "TFL snapped", S);
    fixTiledAxesCommon(axs2(2));
    makeAxesDataSquare(axs2(2), opts.PadFrac);

    drawnow;

    % ---------------------------------------------------------------------
    % Export
    % ---------------------------------------------------------------------
    files = struct( ...
        'layered_pdf', "", ...
        'layered_png', "", ...
        'snapped_pdf', "", ...
        'snapped_png', "" );

    if opts.Export
        exportDir   = char(string(opts.ExportDir));
        baseLayered = char(string(opts.ExportBaseLayered));
        baseSnapped = char(string(opts.ExportBaseSnapped));

        if ~exist(exportDir, 'dir')
            mkdir(exportDir);
        end

        pdf1 = fullfile(exportDir, sprintf('%s.pdf', baseLayered));
        pdf2 = fullfile(exportDir, sprintf('%s.pdf', baseSnapped));
        png1 = fullfile(exportDir, sprintf('%s.png', baseLayered));
        png2 = fullfile(exportDir, sprintf('%s.png', baseSnapped));

        if opts.ExportPDF
            if opts.PDFVector
                exportgraphics(fig1, pdf1, 'ContentType', 'vector');
                exportgraphics(fig2, pdf2, 'ContentType', 'vector');
            else
                exportgraphics(fig1, pdf1, 'ContentType', 'image', ...
                    'Resolution', opts.PNGResolution);
                exportgraphics(fig2, pdf2, 'ContentType', 'image', ...
                    'Resolution', opts.PNGResolution);
            end

            files.layered_pdf = string(pdf1);
            files.snapped_pdf = string(pdf2);
        end

        if opts.ExportPNG
            exportgraphics(fig1, png1, 'Resolution', opts.PNGResolution);
            exportgraphics(fig2, png2, 'Resolution', opts.PNGResolution);

            files.layered_png = string(png1);
            files.snapped_png = string(png2);
        end
    end

    %{
    % ---------------------------------------------------------------------
    % Export
    % ---------------------------------------------------------------------
    files = struct( ...
        'layered_pdf', "", ...
        'layered_png', "", ...
        'snapped_pdf', "", ...
        'snapped_png', "" );

    if opts.Export
        if ~exist(opts.ExportDir, 'dir'), mkdir(opts.ExportDir); end

        base1 = fullfile(opts.ExportDir, opts.ExportBaseLayered);
        base2 = fullfile(opts.ExportDir, opts.ExportBaseSnapped);

        if opts.ExportPDF
            pdf1 = string(base1) + ".pdf";
            pdf2 = string(base2) + ".pdf";

            if opts.PDFVector
                exportgraphics(fig1, pdf1, 'ContentType', 'vector');
                exportgraphics(fig2, pdf2, 'ContentType', 'vector');
            else
                exportgraphics(fig1, pdf1, 'ContentType', 'image', ...
                    'Resolution', opts.PNGResolution);
                exportgraphics(fig2, pdf2, 'ContentType', 'image', ...
                    'Resolution', opts.PNGResolution);
            end

            files.layered_pdf = pdf1;
            files.snapped_pdf = pdf2;
        end

        if opts.ExportPNG
            png1 = string(base1) + ".png";
            png2 = string(base2) + ".png";

            exportgraphics(fig1, png1, 'Resolution', opts.PNGResolution);
            exportgraphics(fig2, png2, 'Resolution', opts.PNGResolution);

            files.layered_png = png1;
            files.snapped_png = png2;
        end
    end
    %}

    % ---------------------------------------------------------------------
    % Return
    % ---------------------------------------------------------------------
    out = struct();
    out.fig_layered  = fig1;
    out.fig_snapped  = fig2;
    out.axes_layered = axs1;
    out.axes_snapped = axs2;
    out.files        = files;
    out.W            = W;
    out.labels       = labels;
end

% =========================================================================
function opts = applyDefaults_(opts)

    opts = setDefault_(opts, 'ShowLabels',        true);
    opts = setDefault_(opts, 'Labels',            {});
    opts = setDefault_(opts, 'ArrowSize',         0.03);
    opts = setDefault_(opts, 'PadFrac',           0.02);
    opts = setDefault_(opts, 'TileSpacing',       'compact');
    opts = setDefault_(opts, 'Padding',           'compact');
    opts = setDefault_(opts, 'FigWidthCM',        16);
    opts = setDefault_(opts, 'FigHeightCM',       8.5);
    opts = setDefault_(opts, 'SnappedMinSepFrac', 0.05);

    opts = setDefault_(opts, 'Export',            false);
    opts = setDefault_(opts, 'ExportDir',         pwd);
    opts = setDefault_(opts, 'ExportBaseLayered', 'fig_foodweb_tfl_vs_sugiyama');
    opts = setDefault_(opts, 'ExportBaseSnapped', 'fig_foodweb_tfl_vs_snapped');
    opts = setDefault_(opts, 'ExportPDF',         true);
    opts = setDefault_(opts, 'ExportPNG',         true);
    opts = setDefault_(opts, 'PNGResolution',     300);
    opts = setDefault_(opts, 'PDFVector',         true);

    opts = setDefault_(opts, 'LabelFontSize',     8);
    opts = setDefault_(opts, 'LabelColor',        [0.15 0.15 0.15]);
    opts = setDefault_(opts, 'AdjustLabels',      false);
    opts = setDefault_(opts, 'LabelHalo',         false);
    opts = setDefault_(opts, 'LabelHaloPadFrac',  0.006);
    opts = setDefault_(opts, 'LabelHaloAlpha',    0.30);
    opts = setDefault_(opts, 'LabelHaloColor',    [1 1 1]);
end

function opts = setDefault_(opts, fieldName, value)
    if ~isfield(opts, fieldName) || isempty(opts.(fieldName))
        opts.(fieldName) = value;
    end
end

%{
function opts = applyDefaults_(opts)
    if ~isfield(opts, 'ShowLabels'),        opts.ShowLabels = true; end
    if ~isfield(opts, 'Labels'),            opts.Labels = {}; end
    if ~isfield(opts, 'ArrowSize'),         opts.ArrowSize = 0.03; end
    if ~isfield(opts, 'PadFrac'),           opts.PadFrac = 0.02; end
    if ~isfield(opts, 'TileSpacing'),       opts.TileSpacing = 'compact'; end
    if ~isfield(opts, 'Padding'),           opts.Padding = 'compact'; end
    if ~isfield(opts, 'FigWidthCM'),        opts.FigWidthCM = 16; end
    if ~isfield(opts, 'FigHeightCM'),       opts.FigHeightCM = 8.5; end
    if ~isfield(opts, 'SnappedMinSepFrac'), opts.SnappedMinSepFrac = 0.05; end

    if ~isfield(opts, 'Export'),            opts.Export = false; end
    if ~isfield(opts, 'ExportDir'),         opts.ExportDir = pwd; end
    if ~isfield(opts, 'ExportBaseLayered'), opts.ExportBaseLayered = 'fig_foodweb_tfl_vs_sugiyama'; end
    if ~isfield(opts, 'ExportBaseSnapped'), opts.ExportBaseSnapped = 'fig_foodweb_tfl_vs_snapped'; end
    if ~isfield(opts, 'ExportPDF'),         opts.ExportPDF = true; end
    if ~isfield(opts, 'ExportPNG'),         opts.ExportPNG = true; end
    if ~isfield(opts, 'PNGResolution'),     opts.PNGResolution = 300; end
    if ~isfield(opts, 'PDFVector'),         opts.PDFVector = true; end

    if ~isfield(opts, 'LabelFontSize'),    opts.LabelFontSize = 8; end
    if ~isfield(opts, 'LabelColor'),       opts.LabelColor = [0.15 0.15 0.15]; end
    if ~isfield(opts, 'AdjustLabels'),     opts.AdjustLabels = false; end
    if ~isfield(opts, 'LabelHalo'),        opts.LabelHalo = false; end
    if ~isfield(opts, 'LabelHaloPadFrac'), opts.LabelHaloPadFrac = 0.006; end
    if ~isfield(opts, 'LabelHaloAlpha'),   opts.LabelHaloAlpha = 0.30; end
    if ~isfield(opts, 'LabelHaloColor'),   opts.LabelHaloColor = [1 1 1]; end
end
%}

% =========================================================================
function labels = prepareLabels_(labelsIn, n)
% Convert labels to cellstr and validate length if nonempty.

    if isempty(labelsIn)
        labels = {};
        return;
    end

    if isstring(labelsIn)
        labels = cellstr(labelsIn(:));
    elseif iscell(labelsIn)
        labels = labelsIn(:);
    else
        error('fig_food_web:BadLabels', ...
            'Labels must be empty, a cell array, or a string array.');
    end

    if numel(labels) ~= n
        error('fig_food_web:BadLabelsLength', ...
            'Labels must have length equal to number of nodes (%d).', n);
    end
end