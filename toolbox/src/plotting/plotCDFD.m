function out = plotCDFD(W, varargin)
%PLOTCDFD  Visualise a circular-directional flow decomposition.
%
%   out = plotCDFD(W)
%   out = plotCDFD(W, 'Mode', 'split', ...)
%
% Modes
% -----
%   'split'
%       Two-panel display:
%         left  = directional component D shown with tflPlot
%         right = circular component C shown with ringFlowPlot
%
%   'overlay'
%       Single-panel display using tflPlot(...,'RenderMode','overlay')
%       on the full matrix W, with C and D passed through.
%
%   'cdfdcolor'  (aliases: 'colorcdfd', 'color')
%       Single-panel display using tflPlot(...,'RenderMode','cdfdcolor')
%       on the full matrix W, with C and D passed through.
%
% Notes
% -----
% - In overlay and cdfdcolor modes, this function is intentionally thin:
%   it resolves the decomposition, prepares defaults, and then delegates
%   rendering to tflPlot.
% - In split mode, the default node-size policy is 'component_shared':
%   each panel uses node sizes based on the component being displayed, but
%   with a common size scale across the two panels.
%
% INPUT
%   W : square weighted adjacency matrix (rows = source, cols = destination)
%
% NAME-VALUE OPTIONS
%   'Mode'                : 'split' (default), 'overlay', or 'cdfdcolor'
%                           aliases: 'colorcdfd', 'color'
%   'SplitLayoutSource'   : split mode only; matrix used to compute the
%                           directional-panel TFL layout:
%                           'w' (default) or 'd'
%   'SplitNodeSizeMode'   : split mode only; how node sizes are handled
%                           across the two panels:
%                           'component_shared' (default)
%                               node size in each panel reflects incident
%                               flow in that component (uD on left, uC on
%                               right), using a common size scale across
%                               both panels
%                           'component_local'
%                               node size in each panel reflects incident
%                               flow in that component, with each panel
%                               scaled independently
%                           'full'
%                               both panels use full-network node sizes uW
%   'C'                   : optional precomputed circular component
%   'D'                   : optional precomputed directional component
%   'Info'                : optional decomposition info struct
%   'ComputeIfMissing'    : true (default)
%   'DecompositionMethod' : 'bff' (default)
%   'DecompositionArgs'   : struct passed through to cdfd_bff
%
%   'Labels'              : node labels (string/cellstr)
%   'RingOrder'           : optional explicit ring order for split mode
%
%   'TFLPlotOpts'         : struct passed to tflPlot as PlotOpts
%   'TFLLayoutOpts'       : struct passed to tflPlot as LayoutOpts
%   'RingPlotOpts'        : struct passed to ringFlowPlot (split mode only)
%   'FigureOpts'          : struct for figure/panel formatting
%   'OverlayOpts'         : deprecated; accepted but ignored
%
% OUTPUT
%   out : struct with fields including
%       .W, .C, .D, .Wplot, .Cplot, .Dplot, .info
%       .fig, .tiled
%       .uW, .uC, .uD
%   and mode-specific fields
%       split     : .axD, .axC, .ringOrder, .X, .Y, .h, .layoutInfo, .sceneD, .geomD
%       overlay/* : .ax, .X, .Y, .compIdx, .layoutInfo, .scene, .geom
%
% See also tflPlot, ringFlowPlot, cdfd_bff

    % -------------------- checks --------------------
    if ~ismatrix(W) || size(W,1) ~= size(W,2)
        error('plotCDFD:WNotSquare', 'W must be a square adjacency matrix.');
    end
    if ~isnumeric(W) || any(~isfinite(W(:)))
        error('plotCDFD:BadW', 'W must be a finite numeric matrix.');
    end
    W = sparse(double(W));
    n = size(W,1);

    % -------------------- parse --------------------
    opts = parsePlotCDFDOpts_(varargin{:});
    warnDeprecatedPlotCDFDOpts_(opts);

    labels = normaliseLabels_(opts.Labels, n);

    % -------------------- decomposition --------------------
    [C, D, info] = resolveCDFDPlot_(W, opts);
    info = fillCDFDInfoPlot_(W, C, D, info);

    % -------------------- raw node-size summaries --------------------
    uW = full(sum(W,2) + sum(W,1)');
    uC = full(sum(C,2) + sum(C,1)');
    uD = full(sum(D,2) + sum(D,1)');

    % Shared size-scaling reference for split mode
    %uCD = [uD(:); uC(:)];

    % For compatibility with older wrapper expectations
    Wplot = W;
    Cplot = C;
    Dplot = D;

    % -------------------- figure defaults --------------------
    figOpts = opts.FigureOpts;
    figOpts = setDefaultFieldPlot_(figOpts, 'Position',        [100 100 1500 780]);
    figOpts = setDefaultFieldPlot_(figOpts, 'Color',           'w');
    figOpts = setDefaultFieldPlot_(figOpts, 'TileSpacing',     'compact');
    figOpts = setDefaultFieldPlot_(figOpts, 'Padding',         'loose');
    figOpts = setDefaultFieldPlot_(figOpts, 'PanelLineWidth',  0.75);
    figOpts = setDefaultFieldPlot_(figOpts, 'MainTitle',       'CDFD (BFF)');
    figOpts = setDefaultFieldPlot_(figOpts, 'ShowMainTitle',   true);
    figOpts = setDefaultFieldPlot_(figOpts, 'ShowPanelTitles', true);

    mode = lower(strtrim(opts.Mode));

    switch mode

        % ==============================================================
        case 'split'

            % ---------- ring ordering ----------
            if ~isempty(opts.RingOrder)
                ringOrder = opts.RingOrder(:);
            elseif isfield(opts.RingPlotOpts, 'LayoutOrder') && ~isempty(opts.RingPlotOpts.LayoutOrder)
                ringOrder = opts.RingPlotOpts.LayoutOrder(:);
            else
                Sorder = C + C.';
                Sorder(1:size(Sorder,1)+1:end) = 0;
                ringOrder = circularSpectralOrderPlot_(Sorder);
            end

            %{
            if numel(ringOrder) ~= n || ...
               numel(unique(ringOrder)) ~= n || ...
               ~isequal(sort(ringOrder), (1:n)')
                error('plotCDFD:BadRingOrder', ...
                    'RingOrder must be a permutation of 1:n.');
            end
            %}

            ringOrder = sanitiseRingOrderPlot_(ringOrder, n);


            % ---------- split-mode node sizing policy ----------
            [nodeSizeDataD, nodeSizeDataC, nodeSizeScaleDataD, nodeSizeScaleDataC] = ...
                resolveSplitNodeSizingPlot_(uW, uD, uC, opts.SplitNodeSizeMode);


            % ---------- defaults: directional TFL panel ----------
            tflPlotOpts   = opts.TFLPlotOpts;
            tflLayoutOpts = opts.TFLLayoutOpts;

            tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'Labels', labels);
            tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'ShowLabels', true);
            tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'LabelHalo', true);
            tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'LabelHaloAlpha', 0.25);
            tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'LabelHaloPadFrac', 0.004);
            tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'LabelColor', [1 0 0]);
            tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'LabelFontSize', 8);
            tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'AdjustLabels', true);
            tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'XScale', 1.6);
            %tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'NodeSizeData', uD);
            tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'NodeSizeData', nodeSizeDataD);
            tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'NodeSizeScaleMode', 'area');
            tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'NodeSizeScaleData', nodeSizeScaleDataD);
            tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'NodeSizeRange', [70 260]);
            tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'EdgeWidthMode', 'weight');
            tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'EdgeWidthScale', 'log');
            tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'EdgeWidthRange', [0.05 6.0]);
            tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'EdgeWidthQuantiles', [0.20 0.99]);
            tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'ArrowLenBasis', 'ltyp');
            tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'ArrowSizeMode', 'fixed');
            tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'ArrowSize', 0.06);
            tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'ArrowMinFrac', 0);
            tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'ArrowOverNodes', true);
            tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'ArrowNodeGapPts', 1.0);
            tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'PadFracX', 0.05);
            tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'PadFracY', 0.03);
            tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'CurvPadScaleX', 1.0);
            tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'CurvPadScaleY', 0.12);



            % ---------- defaults: circular ring panel ----------
            ringPlotOpts = opts.RingPlotOpts;
            ringPlotOpts = setDefaultFieldPlot_(ringPlotOpts, 'LayoutOrder', ringOrder);
            ringPlotOpts = setDefaultFieldPlot_(ringPlotOpts, 'Labels', labels);
            ringPlotOpts = setDefaultFieldPlot_(ringPlotOpts, 'ShowLabels', true);
            ringPlotOpts = setDefaultFieldPlot_(ringPlotOpts, 'LabelColor', [1 0 0]);
            ringPlotOpts = setDefaultFieldPlot_(ringPlotOpts, 'LabelFontSize', 8);
            ringPlotOpts = setDefaultFieldPlot_(ringPlotOpts, 'LabelRadiusFactor', 1.12);
            %ringPlotOpts = setDefaultFieldPlot_(ringPlotOpts, 'NodeSizeData', uC);
            ringPlotOpts = setDefaultFieldPlot_(ringPlotOpts, 'NodeSizeData', nodeSizeDataC);
            ringPlotOpts = setDefaultFieldPlot_(ringPlotOpts, 'NodeSizeScaleMode', 'area');
            ringPlotOpts = setDefaultFieldPlot_(ringPlotOpts, 'NodeSizeScaleData', nodeSizeScaleDataC);
            ringPlotOpts = setDefaultFieldPlot_(ringPlotOpts, 'NodeSizeRange', [70 260]);
            ringPlotOpts = setDefaultFieldPlot_(ringPlotOpts, 'ShowRing', true);
            ringPlotOpts = setDefaultFieldPlot_(ringPlotOpts, 'RingWidth', 0.75);
            ringPlotOpts = setDefaultFieldPlot_(ringPlotOpts, 'EdgeAlpha', 0.75);
            ringPlotOpts = setDefaultFieldPlot_(ringPlotOpts, 'EdgeWidthMode', 'weight');
            ringPlotOpts = setDefaultFieldPlot_(ringPlotOpts, 'EdgeWidthScale', 'log');
            ringPlotOpts = setDefaultFieldPlot_(ringPlotOpts, 'EdgeWidthRange', [0.05 6.0]);
            ringPlotOpts = setDefaultFieldPlot_(ringPlotOpts, 'EdgeWidthQuantiles', [0.20 0.99]);
            ringPlotOpts = setDefaultFieldPlot_(ringPlotOpts, 'EdgeColorMode', 'ringdirection');
            ringPlotOpts = setDefaultFieldPlot_(ringPlotOpts, 'MutualOverlay', true);
            ringPlotOpts = setDefaultFieldPlot_(ringPlotOpts, 'Curvature', 0.78);
            ringPlotOpts = setDefaultFieldPlot_(ringPlotOpts, 'ArrowAngle', pi/7);
            ringPlotOpts = setDefaultFieldPlot_(ringPlotOpts, 'ArrowNodeGapPts', 1.0);
            ringPlotOpts = setDefaultFieldPlot_(ringPlotOpts, 'MutualLaneOffset', 0.5);
            ringPlotOpts = setDefaultFieldPlot_(ringPlotOpts, 'ArrowSize', 0.06);
            ringPlotOpts = setDefaultFieldPlot_(ringPlotOpts, 'ArrowBacktrackFactor', 1.8);
            ringPlotOpts = setDefaultFieldPlot_(ringPlotOpts, 'ArrowPlacement', 'late');
            ringPlotOpts = setDefaultFieldPlot_(ringPlotOpts, 'ArrowPositionFrac', 0.92);
            ringPlotOpts = setDefaultFieldPlot_(ringPlotOpts, 'CurveNpts', 80);
            ringPlotOpts = setDefaultFieldPlot_(ringPlotOpts, 'StartAngle', pi/2);
            ringPlotOpts = setDefaultFieldPlot_(ringPlotOpts, 'Clockwise', false);
            ringPlotOpts = setDefaultFieldPlot_(ringPlotOpts, 'Radius', 1.0);
            ringPlotOpts = setDefaultFieldPlot_(ringPlotOpts, 'Title', '');

            % ---------- figure ----------
            fig = figure('Color', figOpts.Color, 'Position', figOpts.Position);
            tl = tiledlayout(fig, 1, 2, ...
                'TileSpacing', figOpts.TileSpacing, ...
                'Padding', figOpts.Padding);

            if figOpts.ShowMainTitle
                title(tl, figOpts.MainTitle, 'Interpreter', 'none');
            end

            % ---------- choose layout source for split directional panel ----------
            splitLayoutSource = lower(strtrim(opts.SplitLayoutSource));

            switch splitLayoutSource
                case 'w'
                    Wlayout = W;
                    uLayout = uW;
                case 'd'
                    Wlayout = D;
                    uLayout = uD;
                otherwise
                    error('plotCDFD:BadSplitLayoutSource', ...
                        'SplitLayoutSource must be ''w'' or ''d''.');
            end

            % Layout defaults should follow the chosen layout source
            tflLayoutOpts = setDefaultFieldPlot_(tflLayoutOpts, 'InitMode', 'spectral');
            tflLayoutOpts = setDefaultFieldPlot_(tflLayoutOpts, 'InitYNoiseScale', 0.2);
            tflLayoutOpts = setDefaultFieldPlot_(tflLayoutOpts, 'FreePhaseFrac', 0.2);
            tflLayoutOpts = setDefaultFieldPlot_(tflLayoutOpts, 'Seed', 1);
            tflLayoutOpts = setDefaultFieldPlot_(tflLayoutOpts, 'NodeSizeData', uLayout);
            tflLayoutOpts = setDefaultFieldPlot_(tflLayoutOpts, 'SizeAwareDeoverlap', true);
            tflLayoutOpts = setDefaultFieldPlot_(tflLayoutOpts, 'NodeRadiusRange', [0.05 0.14]);
            tflLayoutOpts = setDefaultFieldPlot_(tflLayoutOpts, 'OverlapBandwidthY', 0.6);
            tflLayoutOpts = setDefaultFieldPlot_(tflLayoutOpts, 'NodeGap', 0.02);
            tflLayoutOpts = setDefaultFieldPlot_(tflLayoutOpts, 'DeoverlapPasses', 5);

            layoutArgs = structToNameValuePlot_(tflLayoutOpts);
            [Xlay, Ylay, hlay] = trophicLayout(Wlayout, layoutArgs{:});

            layoutInfoD = struct();
            layoutInfoD.LayoutSource = splitLayoutSource;
            layoutInfoD.LayoutMatrixName = upper(splitLayoutSource);

            
            % ---------- panel 1: D rendered on shared coordinates ----------
            axD = nexttile(tl, 1);
            set(axD, 'Color', 'w');

            tflPlotOpts.Parent = axD;
            tflPlotOpts.ClearAxes = true;
            tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'LevelSemantics', 'trophic');

            plotArgs = structToNameValuePlot_(tflPlotOpts);
            [XD, YD, sceneD, geomD] = plotTFL(D, Xlay, Ylay, hlay, plotArgs{:});

            if figOpts.ShowPanelTitles
                title(axD, sprintf('Directional part D  (share = %.3f)', info.directionality), ...
                    'Interpreter', 'none');
            end
            formatPanelAxesPlot_(axD, figOpts);

            %{
            % ---------- panel 1: D via tflPlot ----------
            axD = nexttile(tl, 1);
            set(axD, 'Color', 'w');

            tflPlotOpts.Parent = axD;
            tflPlotOpts.ClearAxes = true;

            [XD, YD, compIdxD, layoutInfoD, sceneD, geomD] = tflPlot(D, ...
                'LayoutOpts', tflLayoutOpts, ...
                'PlotOpts',   tflPlotOpts);

            if figOpts.ShowPanelTitles
                title(axD, sprintf('Directional part D  (share = %.3f)', info.directionality), ...
                    'Interpreter', 'none');
            end
            formatPanelAxesPlot_(axD, figOpts);

            %}

            % ---------- panel 2: C via ringFlowPlot ----------
            axC = nexttile(tl, 2);
            set(axC, 'Color', 'w');

            ringPlotOpts.Parent = axC;
            ringArgs = structToNameValuePlot_(ringPlotOpts);
            %ringFlowPlot(C, ringArgs{:});
            outRing = ringFlowPlot(C, ringArgs{:});

            if figOpts.ShowPanelTitles
                title(axC, sprintf('Circular part C  (share = %.3f)', info.circularity), ...
                    'Interpreter', 'none');
            end
            formatPanelAxesPlot_(axC, figOpts);

            drawnow;
            equaliseSplitPanelBoxesPlot_(axD, axC);
            drawnow;

            out = struct();
            out.mode = 'split';
            out.W = W;
            out.C = C;
            out.D = D;
            out.Wplot = Wplot;
            out.Cplot = Cplot;
            out.Dplot = Dplot;
            out.info = info;
            out.fig = fig;
            out.tiled = tl;
            out.axD = axD;
            out.axC = axC;
            out.ringOrder = ringOrder;
            out.X = XD;
            out.Y = YD;
            out.h = hlay;
            out.layoutInfo = layoutInfoD;
            out.sceneD = sceneD;
            out.geomD = geomD;
            %out.X = XD;
            %out.Y = YD;
            %out.compIdx = compIdxD;
            %out.layoutInfo = layoutInfoD;
            %out.sceneD = sceneD;
            %out.geomD = geomD;
            out.uW = uW;
            out.uC = uC;
            out.uD = uD;
            
            out.splitNodeSizeMode = opts.SplitNodeSizeMode;
            out.nodeSizeDataD = nodeSizeDataD;
            out.nodeSizeDataC = nodeSizeDataC;
            %out.nodeSizeScaleData = nodeSizeScaleDataD;
            out.nodeSizeScaleDataD = nodeSizeScaleDataD;
            out.nodeSizeScaleDataC = nodeSizeScaleDataC;

            % Temporary outputs for debugging / inspection
            out.ringOut = outRing;
            out.nodeSizeD = out.sceneD.nodeSize;
            out.nodeSizeC = out.ringOut.nodeSize;

        % ==============================================================
        case 'overlay'

            [out, fig, tl, ax] = makeSinglePanelCDFD_( ...
                W, C, D, info, labels, uW, opts, figOpts, 'overlay');

            out.mode = 'overlay';
            out.fig = fig;
            out.tiled = tl;
            out.ax = ax;
            out.axMain = ax;

        % ==============================================================
        case {'cdfdcolor','colorcdfd','color'}

            [out, fig, tl, ax] = makeSinglePanelCDFD_( ...
                W, C, D, info, labels, uW, opts, figOpts, 'cdfdcolor');

            out.mode = 'cdfdcolor';
            out.fig = fig;
            out.tiled = tl;
            out.ax = ax;
            out.axMain = ax;

        otherwise
            error('plotCDFD:UnknownMode', 'Unknown Mode: %s', opts.Mode);
    end
end

% =====================================================================
function [out, fig, tl, ax] = makeSinglePanelCDFD_(W, C, D, info, labels, uW, opts, figOpts, renderMode)

    tflPlotOpts   = opts.TFLPlotOpts;
    tflLayoutOpts = opts.TFLLayoutOpts;

    % Sensible wrapper defaults, but rendering is delegated entirely to tflPlot
    tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'Labels', labels);
    tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'ShowLabels', true);
    tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'LabelHalo', true);
    tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'LabelHaloAlpha', 0.25);
    tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'LabelHaloPadFrac', 0.004);
    tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'LabelColor', [1 0 0]);
    tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'LabelFontSize', 8);
    tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'AdjustLabels', true);
    tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'XScale', 1.6);
    tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'NodeSizeData', uW);
    tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'NodeSizeScaleMode', 'area');
    tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'NodeSizeRange', [70 260]);
    tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'EdgeWidthMode', 'weight');
    tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'EdgeWidthScale', 'log');
    tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'EdgeWidthRange', [0.05 6.0]);
    tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'EdgeWidthQuantiles', [0.20 0.99]);
    tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'ArrowLenBasis', 'ltyp');
    tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'ArrowSizeMode', 'fixed');
    tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'ArrowSize', 0.12);
    tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'ArrowMinFrac', 0);
    tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'ArrowOverNodes', true);
    tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'ArrowNodeGapPts', 1.0);
    tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'PadFracX', 0.05);
    tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'PadFracY', 0.03);
    tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'CurvPadScaleX', 1.0);
    tflPlotOpts = setDefaultFieldPlot_(tflPlotOpts, 'CurvPadScaleY', 0.12);

    tflLayoutOpts = setDefaultFieldPlot_(tflLayoutOpts, 'InitMode', 'spectral');
    tflLayoutOpts = setDefaultFieldPlot_(tflLayoutOpts, 'InitYNoiseScale', 0.2);
    tflLayoutOpts = setDefaultFieldPlot_(tflLayoutOpts, 'FreePhaseFrac', 0.2);
    tflLayoutOpts = setDefaultFieldPlot_(tflLayoutOpts, 'Seed', 1);
    tflLayoutOpts = setDefaultFieldPlot_(tflLayoutOpts, 'NodeSizeData', uW);
    tflLayoutOpts = setDefaultFieldPlot_(tflLayoutOpts, 'SizeAwareDeoverlap', true);
    tflLayoutOpts = setDefaultFieldPlot_(tflLayoutOpts, 'NodeRadiusRange', [0.05 0.14]);
    tflLayoutOpts = setDefaultFieldPlot_(tflLayoutOpts, 'OverlapBandwidthY', 0.6);
    tflLayoutOpts = setDefaultFieldPlot_(tflLayoutOpts, 'NodeGap', 0.02);
    tflLayoutOpts = setDefaultFieldPlot_(tflLayoutOpts, 'DeoverlapPasses', 5);

    % Force the rendering route
    tflPlotOpts.RenderMode = renderMode;
    tflPlotOpts.C = C;
    tflPlotOpts.D = D;
    tflPlotOpts.Info = info;

    fig = figure('Color', figOpts.Color, 'Position', figOpts.Position);
    tl = tiledlayout(fig, 1, 1, ...
        'TileSpacing', figOpts.TileSpacing, ...
        'Padding', figOpts.Padding);

    ax = nexttile(tl, 1);
    set(ax, 'Color', 'w');

    tflPlotOpts.Parent = ax;
    tflPlotOpts.ClearAxes = true;

    if figOpts.ShowMainTitle
        title(tl, figOpts.MainTitle, 'Interpreter', 'none');
    end

    [X, Y, compIdx, layoutInfo, scene, geom] = tflPlot(W, ...
        'LayoutOpts', tflLayoutOpts, ...
        'PlotOpts',   tflPlotOpts);

    formatPanelAxesPlot_(ax, figOpts);

    out = struct();
    out.W = W;
    out.C = C;
    out.D = D;
    out.Wplot = W;
    out.Cplot = C;
    out.Dplot = D;
    out.info = info;
    out.X = X;
    out.Y = Y;
    out.compIdx = compIdx;
    out.layoutInfo = layoutInfo;
    out.scene = scene;
    out.geom = geom;
    out.uW = uW;
    out.uC = full(sum(C,2) + sum(C,1)');
    out.uD = full(sum(D,2) + sum(D,1)');
end

% =====================================================================
function opts = parsePlotCDFDOpts_(varargin)

    opts = struct();

    opts.Mode = 'split';

    opts.SplitNodeSizeMode = 'component_shared';
    
    opts.SplitLayoutSource = 'w';

    opts.C = [];
    opts.D = [];
    opts.Info = [];

    opts.ComputeIfMissing = true;
    opts.DecompositionMethod = 'bff';
    opts.DecompositionArgs = struct( ...
        'ToleranceZero', 1e-10, ...
        'Validate', true);

    opts.Labels = {};

    opts.RingOrder = [];

    opts.TFLPlotOpts   = struct();
    opts.TFLLayoutOpts = struct();
    opts.RingPlotOpts  = struct();
    opts.OverlayOpts   = struct();    % deprecated, ignored

    opts.FigureOpts = struct();

    if isempty(varargin)
        return;
    end

    if mod(numel(varargin), 2) ~= 0
        error('plotCDFD:BadArgs', ...
            'Optional arguments must be name/value pairs.');
    end

    for k = 1:2:numel(varargin)
        name = varargin{k};
        value = varargin{k+1};

        if ~ischar(name) && ~isstring(name)
            error('plotCDFD:BadParamName', ...
                'Parameter names must be strings.');
        end

        switch lower(char(name))
            case 'mode'
                opts.Mode = char(value);

            case 'splitnodesizemode'
                opts.SplitNodeSizeMode = char(value);
            
            case 'splitlayoutsource'
                opts.SplitLayoutSource = char(value);

            case 'c'
                opts.C = value;
            case 'd'
                opts.D = value;
            case 'info'
                opts.Info = value;

            case 'computeifmissing'
                opts.ComputeIfMissing = logical(value);
            case 'decompositionmethod'
                opts.DecompositionMethod = char(value);
            case 'decompositionargs'
                opts.DecompositionArgs = value;

            case 'labels'
                opts.Labels = value;

            case 'ringorder'
                opts.RingOrder = value;

            case 'tflplotopts'
                opts.TFLPlotOpts = value;
            case 'tfllayoutopts'
                opts.TFLLayoutOpts = value;
            case 'ringplotopts'
                opts.RingPlotOpts = value;
            case 'overlayopts'
                opts.OverlayOpts = value;

            case 'figureopts'
                opts.FigureOpts = value;

            otherwise
                error('plotCDFD:UnknownOption', ...
                    'Unknown option: %s', char(name));
        end
    end
end

% =====================================================================
function warnDeprecatedPlotCDFDOpts_(opts)

    if isstruct(opts.OverlayOpts) && ~isempty(fieldnames(opts.OverlayOpts))
        warning('plotCDFD:OverlayOptsDeprecated', ...
            ['OverlayOpts is deprecated and ignored. ', ...
             'Configure overlay/cdfdcolor behaviour through TFLPlotOpts instead.']);
    end
end

% =====================================================================
function [C, D, info] = resolveCDFDPlot_(W, opts)

    n = size(W,1);

    if isempty(opts.C) || isempty(opts.D)
        if ~opts.ComputeIfMissing
            error('plotCDFD:MissingCD', ...
                'C and D were not supplied and ComputeIfMissing is false.');
        end

        switch lower(strtrim(opts.DecompositionMethod))
            case 'bff'
                cdfdArgs = structToNameValuePlot_(opts.DecompositionArgs);
                [C, D, info] = cdfd_bff(W, cdfdArgs{:});
            otherwise
                error('plotCDFD:UnknownDecompositionMethod', ...
                    'Unknown DecompositionMethod: %s', opts.DecompositionMethod);
        end
    else
        C = sparse(double(opts.C));
        D = sparse(double(opts.D));

        if ~isequal(size(C), [n n]) || ~isequal(size(D), [n n])
            error('plotCDFD:BadCDSize', ...
                'C and D must be the same size as W.');
        end

        if isempty(opts.Info)
            info = struct();
        else
            info = opts.Info;
        end
    end
end

% =====================================================================
function labels = normaliseLabels_(labelsIn, n)

    if isempty(labelsIn)
        labels = arrayfun(@(k) sprintf('%d', k), 1:n, 'UniformOutput', false);
        return;
    end

    if isstring(labelsIn)
        labels = cellstr(labelsIn(:));
    elseif iscell(labelsIn)
        labels = labelsIn(:);
    else
        error('plotCDFD:BadLabels', ...
            'Labels must be empty, a string array, or a cell array.');
    end

    if numel(labels) ~= n
        error('plotCDFD:BadLabelsLength', ...
            'Labels must have length n.');
    end
end

% =====================================================================
function info = fillCDFDInfoPlot_(W, C, D, info)

    if nargin < 4 || isempty(info)
        info = struct();
    end

    denom = full(sum(W(:)));
    if denom > 0
        circ = full(sum(C(:))) / denom;
        dire = full(sum(D(:))) / denom;
    else
        circ = 0;
        dire = 0;
    end

    if ~isfield(info, 'circularity') || isempty(info.circularity)
        info.circularity = circ;
    end
    if ~isfield(info, 'directionality') || isempty(info.directionality)
        info.directionality = dire;
    end
end

% =====================================================================
function s = setDefaultFieldPlot_(s, fieldName, value)
    if ~isfield(s, fieldName) || isempty(s.(fieldName))
        s.(fieldName) = value;
    end
end

% =====================================================================
function nv = structToNameValuePlot_(s)
    f = fieldnames(s);
    nv = cell(1, 2*numel(f));
    for ii = 1:numel(f)
        nv{2*ii-1} = f{ii};
        nv{2*ii}   = s.(f{ii});
    end
end

% =====================================================================
function formatPanelAxesPlot_(ax, figOpts)
    set(ax, ...
        'Color', 'w', ...
        'LineWidth', figOpts.PanelLineWidth, ...
        'Box', 'on');
end

% =====================================================================
function order = circularSpectralOrderPlot_(S)
% Compute a cyclic order from a symmetric affinity matrix S.

    S = sparse(double(S));
    n = size(S,1);

    if size(S,2) ~= n
        error('circularSpectralOrderPlot_:BadSize', 'S must be square.');
    end
    if n == 0
        order = [];
        return;
    elseif n == 1
        order = 1;
        order = order(:);
        return;
        %{
    elseif n == 1
        order = 1;
        return;
        %}
    end

    S = 0.5 * (S + S.');
    S(1:n+1:end) = 0;
    S = max(S, 0);

    G = graph(S, 'upper');
    comp = conncomp(G);
    nComp = max(comp);

    if nComp > 1
        compStrength = zeros(nComp,1);
        compOrderCell = cell(nComp,1);

        for c = 1:nComp
            idx = find(comp == c);
            Sc  = S(idx, idx);
            compStrength(c) = full(sum(Sc(:)));

            if numel(idx) == 1
                compOrderCell{c} = idx(:);
            else
                ordLocal = circularSpectralOrderPlot_(Sc);
                compOrderCell{c} = idx(ordLocal(:));
                compOrderCell{c} = compOrderCell{c}(:);
            end
            %{
            if numel(idx) == 1
                compOrderCell{c} = idx(:);
            else
                ordLocal = circularSpectralOrderPlot_(Sc);
                compOrderCell{c} = idx(ordLocal(:));
            end
            %}
        end

        [~, ordComp] = sort(compStrength, 'descend');

        blocks = compOrderCell(ordComp);
        blocks = cellfun(@(v) full(double(v(:))), blocks, 'UniformOutput', false);

        order = vertcat(blocks{:});
        order = order(:);
        return;

        %{
        [~, ordComp] = sort(compStrength, 'descend');
        order = vertcat(compOrderCell{ordComp});
        order = order(:);
        return;
        %}
    end

    %{
    if nComp > 1
        compStrength = zeros(nComp,1);
        compOrderCell = cell(nComp,1);

        for c = 1:nComp
            idx = find(comp == c);
            Sc  = S(idx, idx);
            compStrength(c) = full(sum(Sc(:)));

            if numel(idx) == 1
                compOrderCell{c} = idx;
            else
                ordLocal = circularSpectralOrderPlot_(Sc);
                compOrderCell{c} = idx(ordLocal);
            end
        end

        [~, ordComp] = sort(compStrength, 'descend');
        order = vertcat(compOrderCell{ordComp});
        return;
    end
    %}

    if n == 2
        d = full(sum(S,2));
        [~, order] = sort(d, 'descend');
        order = order(:);
        return;
    end

    %{
    if n == 2
        d = full(sum(S,2));
        [~, order] = sort(d, 'descend');
        return;
    end
    %}

    d = full(sum(S,2));
    L = spdiags(d, 0, n, n) - S;

    useFullEig = (n <= 12);
    if useFullEig
        [V, D] = eig(full(L));
        lam = diag(D);
        [~, ordLam] = sort(lam, 'ascend');
        V = V(:, ordLam);
    else
        try
            [V, D] = eigs(L, 3, 'smallestabs');
            lam = diag(D);
            [~, ordLam] = sort(lam, 'ascend');
            V = V(:, ordLam);
        catch
            [V, D] = eig(full(L));
            lam = diag(D);
            [~, ordLam] = sort(lam, 'ascend');
            V = V(:, ordLam);
        end
    end

    if size(V,2) < 3
        if size(V,2) >= 2
            [~, order] = sort(V(:,2), 'ascend');
        else
            order = (1:n).';
        end
        return;
    end

    x = V(:,2);
    y = V(:,3);

    if norm(x) <= eps || norm(y) <= eps
        [~, order] = sort(V(:,2), 'ascend');
        return;
    end

    ang = atan2(y, x);
    [~, order0] = sort(ang, 'ascend');

    affAdj = zeros(n,1);
    for k = 1:n
        a = order0(k);
        b = order0(mod(k, n) + 1);
        affAdj(k) = full(S(a,b));
    end

    [~, kBreak] = min(affAdj);
    %order = order0([kBreak+1:n, 1:kBreak]).';
    order = order0([kBreak+1:n, 1:kBreak]).';
    order = order(:);
end

% =====================================================================
function equaliseSplitPanelBoxesPlot_(ax1, ax2)
    axs = [ax1 ax2];

    set(axs, ...
        'DataAspectRatioMode', 'auto', ...
        'PlotBoxAspectRatioMode', 'manual', ...
        'PositionConstraint', 'innerposition');

    pbaspect(ax1, [1 1 1]);
    pbaspect(ax2, [1 1 1]);
end

function order = sanitiseRingOrderPlot_(order, n)

    order = full(double(order(:)));

    isValid = ...
        numel(order) == n && ...
        all(isfinite(order)) && ...
        all(abs(order - round(order)) < 1e-9);

    if isValid
        order = round(order);
        isValid = ...
            numel(unique(order)) == n && ...
            isequal(sort(order), (1:n)');
    end

    if ~isValid
        warning('plotCDFD:BadRingOrderFallback', ...
            'Computed RingOrder was invalid; falling back to 1:n.');
        order = (1:n).';
    end
end

function [nodeSizeDataD, nodeSizeDataC, nodeSizeScaleDataD, nodeSizeScaleDataC] = ...
    resolveSplitNodeSizingPlot_(uW, uD, uC, splitNodeSizeMode)

    mode = lower(strtrim(splitNodeSizeMode));

    switch mode
        case {'component_shared','componentshared','shared'}
            nodeSizeDataD = uD;
            nodeSizeDataC = uC;
            nodeSizeScaleDataD = [uD(:); uC(:)];
            nodeSizeScaleDataC = nodeSizeScaleDataD;

        case {'component_local','componentlocal','local'}
            nodeSizeDataD = uD;
            nodeSizeDataC = uC;
            nodeSizeScaleDataD = [];
            nodeSizeScaleDataC = [];

        case 'full'
            nodeSizeDataD = uW;
            nodeSizeDataC = uW;
            nodeSizeScaleDataD = [];
            nodeSizeScaleDataC = [];

        otherwise
            error('plotCDFD:BadSplitNodeSizeMode', ...
                ['SplitNodeSizeMode must be ''component_shared'', ', ...
                 '''component_local'', or ''full''.']);
    end
end

