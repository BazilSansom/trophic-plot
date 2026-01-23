function out = fig_scatter_tfl_vs_layered(opts)
%FIG_SCATTER_TFL_VS_LAYERED Compare Sugiyama layered y vs TFL trophic levels.
%
% out = fig_scatter_tfl_vs_layered(opts)
%
% 1x3 figure (when opts.showLayouts=true):
%   (1) scatter: layered YData vs trophic level h
%   (2) MATLAB layered layout
%   (3) TFL layout via tflPlot
%
% If opts.W is provided (NxN adjacency/weight matrix), the function uses it
% directly and skips all network generation.

    if nargin < 1, opts = struct(); end
    opts = applyDefaults_(opts);

    rng(opts.seed);

    % 1) Obtain network W
    if ~isempty(opts.W)
        W = opts.W;
        meta = opts.meta;
        if isempty(meta), meta = struct('source','opts.W'); end
    else
        [W, meta] = generateW_(opts);
    end

    W = double(W);
    n = size(W,1);
    if size(W,2) ~= n
        error('fig_scatter_tfl_vs_layered:BadW', 'opts.W must be square (NxN).');
    end
    W(1:n+1:end) = 0; % no self-loops (layout stability)

    if opts.requireDAG
        Gcheck = digraph(W);
        if ~isdag(Gcheck)
            error('fig_scatter_tfl_vs_layered:NotDAG', ...
                  'opts.requireDAG=true but provided/generated W is not a DAG.');
        end
    end

    G = digraph(W);

    % 2) Layered layout (MATLAB drawing coordinates)
    [y_layered, layeredX, layeredY] = getLayeredXY_(G, opts.direction);

    % 3) Trophic levels (TFL heights)
    [h, info] = trophic_levels(W, opts.trophicArgs{:});
    h = h(:);

    % Optional TFL outputs
    tflX = []; tflY = []; tflCompIdx = []; tflInfo = [];

    % 4) Plot
    wasVisible = get(0,'DefaultFigureVisible');
    if ~opts.visible
        set(0,'DefaultFigureVisible','on');
    end

    fig = figure('Color','w');
    if opts.showLayouts
        tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
    else
        tiledlayout(1,1,'Padding','compact','TileSpacing','compact');
    end

    % ---- Panel 1: scatter (raw rendered y vs trophic) ----
    ax1 = nexttile;
    s = scatter(ax1, y_layered, h, opts.markerSize, 'filled');
    if opts.markerFaceAlpha < 1
        s.MarkerFaceAlpha = opts.markerFaceAlpha;
        s.MarkerEdgeAlpha = opts.markerFaceAlpha;
    end
    grid(ax1,'on');
    xlabel(ax1, 'Sugiyama layered vertical coordinate (rendered y)');
    ylabel(ax1, 'Trophic level (TFL height)');

    if opts.addFitLine
        hold(ax1,'on');
        ok = isfinite(y_layered) & isfinite(h);
        if any(ok)
            pp = polyfit(y_layered(ok), h(ok), 1);
            xf = linspace(min(y_layered(ok)), max(y_layered(ok)), 200);
            yf = polyval(pp, xf);
            plot(ax1, xf, yf, '--', 'LineWidth', 1.2);
        end
        hold(ax1,'off');
    end
    if opts.addTitle
        title(ax1, makeTitle_(opts.titlePrefix, 'Scatter'));
    end
    % Make the plot box square (NOT “equal data units”)
    pbaspect(ax1,[1 1 1]);

    if opts.showLayouts

        % -- Panel 2: Sugiyama layout, but rendered via plotTFL for consistency
        ax2 = nexttile;
        
        % Use layeredY as the “height” argument so up/down + bands match Sugiyama layers
        plotTFL(W, layeredX, layeredY, y_layered, ...
            'Parent', ax2, ...
            'ClearAxes', true, ...
            'FillSquare', true, ...   % IMPORTANT: keep coordinates identical to scatter x-axis “rendered y”
            'ShowLabels', true, ...
            'LabelFontSize', 5, ...
            'LabelOffset', [0 0.3], ...
            'Title', '');   
        
        addPanelLabel(ax2, "Sugiyama layered layout");


        % ---- Panel 3: TFL layout via tflPlot ----
        ax3 = nexttile;

        plotOpts = opts.tflPlotOpts;
        if isempty(plotOpts), plotOpts = struct(); end
        plotOpts.Parent    = ax3;
        plotOpts.ClearAxes = true;
        plotOpts.ShowLabels = true;
        plotOpts.LabelFontSize = 5;
        plotOpts.LabelOffset   = [0 0.3];
        plotOpts.Title     = '';

        layoutOpts = opts.tflLayoutOpts;
        if isempty(layoutOpts), layoutOpts = struct(); end

        try
            [Xt, Yt, compIdx, infoTFL] = tflPlot(W, ...
                'LayoutOpts', layoutOpts, ...
                'PlotOpts',   plotOpts);

            tflX = Xt(:);
            tflY = Yt(:);
            tflCompIdx = compIdx;
            tflInfo = infoTFL;
        catch
            tflPlot(W, 'LayoutOpts', layoutOpts, 'PlotOpts', plotOpts);
        end

        addPanelLabel(ax3, "TFL layout");

        if opts.addTitle
            title(ax3,'TFL layout (tflPlot)');
        end
    end

    if ~opts.visible
        set(0,'DefaultFigureVisible', wasVisible);
    end

    
    % 5) Export (new convention) + legacy fallback
    files = struct('pdf',"",'png',"");

    if isfield(opts,'Export') && opts.Export
        if ~exist(opts.ExportDir,'dir'), mkdir(opts.ExportDir); end

        if isfield(opts,'ExportPDF') && opts.ExportPDF
            pdfPath = fullfile(opts.ExportDir, [opts.ExportBase '.pdf']);
            if isfield(opts,'PDFVector') && opts.PDFVector
                exportgraphics(fig, pdfPath, 'ContentType','vector');
            else
                exportgraphics(fig, pdfPath, 'ContentType','image', ...
                    'Resolution', opts.PNGResolution);
            end
            files.pdf = string(pdfPath);
        end

        if isfield(opts,'ExportPNG') && opts.ExportPNG
            pngPath = fullfile(opts.ExportDir, [opts.ExportBase '.png']);
            exportgraphics(fig, pngPath, 'Resolution', opts.PNGResolution);
            files.png = string(pngPath);
        end

    elseif isfield(opts,'makePdf') && opts.makePdf
        % Legacy behaviour
        if ~exist(opts.outDir,'dir'), mkdir(opts.outDir); end
        pdfPath = fullfile(opts.outDir, opts.fileName);
        set(fig,'PaperPositionMode','auto');
        print(fig, pdfPath, '-dpdf');
        files.pdf = string(pdfPath);
    end

    % 6) Return
    out = struct();
    out.fig       = fig;
    %out.pathPdf   = pathPdf;
    out.W         = W;
    out.meta      = meta;
    out.G         = G;
    out.y_layered = y_layered(:);
    out.h         = h(:);
    out.info      = info;
    out.layeredX  = layeredX(:);
    out.layeredY  = layeredY(:);
    out.tflX      = tflX;
    out.tflY      = tflY;
    out.tflCompIdx = tflCompIdx;
    out.tflInfo   = tflInfo;
    out.files = files;

end

%% ---------------- helpers ----------------

function opts = applyDefaults_(opts)
    d = struct();

    % Optional direct input
    d.W            = [];
    d.meta         = [];
    d.requireDAG   = false;

    % Generator choice
    d.generator  = 'gppm'; % 'gppm' or 'dg'

    % Common
    d.n          = 40;
    d.p          = 0.08;
    d.seed       = 1;

    % randGPPM
    d.T          = 0.35;
    d.gppmArgs   = {};

    % randDG
    d.connected  = true;
    d.weighted   = true;
    d.dgArgs     = {};

    % trophic_levels
    d.trophicArgs = {};

    % layered layout
    d.direction  = 'up';

    % plot controls
    d.showLayouts = true;
    d.visible     = true;
    d.markerSize  = 28;

    % styling
    d.markerFaceAlpha = 0.6;
    d.addFitLine       = false;
    d.addTitle         = false; % caption does the work
    d.titlePrefix      = '';

    % tflPlot passthrough
    d.tflLayoutOpts = struct();
    d.tflPlotOpts   = struct();

    % output
    d.makePdf     = false;
    d.outDir      = '.';
    d.fileName    = 'fig_scatter_tfl_vs_layered.pdf';

    f = fieldnames(d);
    for k = 1:numel(f)
        fn = f{k};
        if ~isfield(opts, fn) || isempty(opts.(fn))
            opts.(fn) = d.(fn);
        end
    end

    % ---------------- Export convention (paper figure pipeline) ----------------
    d.Export        = false;
    d.ExportDir     = '.';
    d.ExportBase    = 'fig_scatter_tfl_vs_layered';
    d.ExportPDF     = true;
    d.ExportPNG     = true;
    d.PNGResolution = 300;
    d.PDFVector     = true;   % if false, rasterize pdf at PNGResolution

end

function [W, meta] = generateW_(opts)
    meta = struct();
    gen = lower(string(opts.generator));

    if gen == "gppm"
        gargs = opts.gppmArgs;
        if ~hasNameValue_(gargs, 'Seed')
            gargs = [gargs, {'Seed', opts.seed}];
        end
        [W, meta] = randGPPM(opts.n, opts.p, opts.T, gargs{:});

    elseif gen == "dg"
        W = randDG(opts.n, opts.p, opts.connected, 'Weighted', opts.weighted, opts.dgArgs{:});
        meta = struct('generator','randDG','connected',opts.connected,'weighted',opts.weighted);

    else
        error('fig_scatter_tfl_vs_layered:BadGenerator', ...
              'opts.generator must be ''gppm'' or ''dg'' (got "%s").', gen);
    end
end

function tf = hasNameValue_(args, name)
    tf = false;
    if isempty(args) || ~iscell(args), return; end
    for k = 1:2:numel(args)
        if ischar(args{k}) || isstring(args{k})
            if strcmpi(string(args{k}), string(name))
                tf = true; return;
            end
        end
    end
end

function [y, X, Y] = getLayeredXY_(G, direction)
    figTmp = figure('Visible','off');
    hp = plot(G, 'Layout','layered', 'Direction', direction);
    X = hp.XData(:);
    Y = hp.YData(:);
    y = Y;
    close(figTmp);
end

function t = makeTitle_(prefix, main)
    if ~isempty(prefix)
        t = sprintf('%s - %s', prefix, main);
    else
        t = main;
    end
end
