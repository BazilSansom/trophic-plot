function [Xplot, Yplot, scene, geom] = plotTFL(W, X, Y, h, varargin)
%PLOTTFL Render a trophic-style layout with optional circular/directional overlays.
%
%   plotTFL(W, X, Y, h)
%   plotTFL(W, X, Y, h, Name, Value, ...)
%   [Xplot, Yplot, scene, geom] = plotTFL(...)
%
% Description
% -----------
% plotTFL renders a weighted directed network at supplied node coordinates
% X and Y. By default the vertical coordinate is interpreted as having
% trophic-step semantics, but the renderer can also be used with other
% supplied level vectors. Optional horizontal reference bands are shown only
% when 'LevelSemantics' is set to 'trophic'.
%
% It supports three edge rendering modes:
%
%   'plain'     : standard TFL rendering using full-flow up/down colours
%
%   'overlay'   : full-flow edge drawn at full width, with one selected
%                 component ('c' or 'd') highlighted on top. The base
%                 full-flow shaft can either use the ordinary plain-mode
%                 up/down colouring ('OverlayBaseColorMode' = 'plain') or
%                 the opposite component colour ('OverlayBaseColorMode' =
%                 'component').
%
%   'cdfdcolor' : full-flow edge drawn at full width, coloured by the
%                 edgewise circular/directional mixture
%
% In the decomposition-based modes, the function works with a circular–
% directional flow decomposition W = C + D. If C and D are not supplied,
% they are computed internally (by default using the BFF decomposition).
%
% The displayed full-flow edge width is always based on the total edge
% weight in W. In overlay mode, the highlighted component width is allocated
% proportionally using the raw edgewise shares C(i,j)/W(i,j) or D(i,j)/W(i,j).
%
%
% Inputs
% ------
% W : NxN numeric matrix
%     Weighted adjacency matrix of the directed network.
%
% X, Y : N-vector
%     Node coordinates to render.
%
% h : N-vector
%     Reference level vector associated with the nodes. By default this is
%     typically a trophic-level vector, but other supplied level vectors may
%     also be used. It is passed through to scene construction and related
%     annotations.
%
%
% Name-Value Options
% ------------------
% Rendering and decomposition
%   'RenderMode'           : 'plain' (default), 'overlay', or 'cdfdcolor'
%                            Accepted aliases for 'cdfdcolor' include
%                            'colorcdfd' and 'color'.
%   'LevelSemantics'       : 'trophic' or 'generic'
%   'C'                    : circular component matrix
%   'D'                    : directional component matrix
%   'Info'                 : decomposition metadata struct
%   'ComputeIfMissing'     : logical, whether to compute C and D if omitted
%                            (default true)
%   'DecompositionMethod'  : currently 'bff'
%   'DecompositionArgs'    : struct of arguments passed to the decomposition
%
% Overlay / colour controls
%   'OverlayComponent'     : component highlighted on top in overlay mode:
%                            'c' or 'd' (default 'c')
%   'OverlayBaseColorMode' : colouring rule for the full-flow base shaft in
%                            overlay mode:
%                              'plain' (default) = use the ordinary plain-mode full-
%                                  flow up/down colours underneath the highlight
%                              'component' = use the opposite
%                                  component colour underneath the highlight
%   'CircularColor'        : RGB colour for circular-flow highlighting
%   'CircularAlpha'        : alpha blend for circular colour
%   'DirectionalColor'     : RGB colour for directional-flow highlighting
%   'DirectionalAlpha'     : alpha blend for directional colour
%   'ArrowLayer' : which shaft carries arrowheads in RenderMode='overlay':
%                  'fullflow' (default), 'overlay', or 'both'
%
%                  'fullflow' = arrowheads annotate the total displayed edge W
%                  'overlay'  = arrowheads annotate only the overlaid component shaft
%                  'both'     = show both full-flow and overlay arrowheads
%   'ArrowColorMode'       : arrowhead colouring rule:
%                              'blend'    = blend circular/directional
%                                           component colours using raw
%                                           edgewise shares
%                              'shaft' = use the colour of the shaft the arrow is being drawn on / with
%                              'overlay'  = use the overlay highlight colour
%
% Edge appearance
%   'EdgeWidthMode'        : 'auto', 'fixed', or 'weight'
%   'EdgeWidthRange'       : [min max] linewidth range for weighted edges
%   'EdgeWidthScale'       : 'log' or 'linear'
%   'EdgeWidthQuantiles'   : quantile range used for clipping before width mapping
%   'UpEdgeColor'          : RGB colour for upward edges in plain mode
%   'DownEdgeColor'        : RGB colour for downward edges in plain mode
%   'EdgeAlpha'            : alpha blend for plain-mode edges
%
% Nodes and labels
%   'ShowNodes'            : logical
%   'NodeColor'            : RGB colour
%   'NodeSize'             : scalar marker area
%   'NodeSizeMode'         : node sizing mode
%   'NodeSizeData'         : optional data whose values are displayed by node size
%   'NodeSizeScaleData'    : optional reference data used to determine the
%                            node-size mapping; if omitted, NodeSizeData is used
%   'ShowLabels'           : logical
%   'Labels'               : cell array of node labels
%   'AdjustLabels'         : logical
%   'LabelHalo'            : logical
%
% Bands / axes / annotation
%   'ShowBands'            : logical
%   'BandTau'              : band spacing parameter
%   'BandColor'            : RGB colour
%   'BandStyle'            : line style
%   'BandWidth'            : line width
%   'Parent'               : target axes
%   'ClearAxes'            : logical
%   'Title'                : plot title
%   'ShowBFF'              : logical
%   'BFF'                  : numeric value to annotate
%
% Additional appearance and geometry controls are also supported; see
% parsePlotOpts_ for the full list.
%
%
% Outputs
% -------
% Xplot, Yplot : N-vectors
%     Coordinates actually rendered after scene construction.
%
% scene : struct
%     Scene data returned by computeTFLScene, including plotted coordinates,
%     node sizes, band positions, and suggested axis limits.
%
% geom : struct
%     Edge geometry returned by computeTFLEdgeGeometry, augmented in the
%     decomposition-based modes with per-edge circular/directional shares
%     and rendered linewidths.
%
%
% Notes
% -----
% - With 'EdgeWidthMode' = 'auto', plain rendering defaults to fixed edge
%   widths, while decomposition-based rendering defaults to weight-based
%   widths.
% - In overlay mode, the base shaft is always the displayed total-flow edge.
%   The selected component is then drawn on top using the corresponding raw
%   share of that displayed width.
% - With 'OverlayBaseColorMode' = 'plain', overlay mode can be used as a
%   component highlight on top of the ordinary TFL full-flow colouring.
% - 'ArrowColorMode' = 'blend' uses component-share blending based on
%   CircularColor and DirectionalColor.
% - If C and D are supplied, they must be the same size as W.
%
% See also computeTFLScene, computeTFLEdgeGeometry, cdfd_bff

    % ---- checks ----
    N = size(W,1);
    if size(W,2) ~= N
        error('plotTFL:WNotSquare', 'W must be a square adjacency matrix.');
    end

    if ~isvector(X) || ~isvector(Y) || ~isvector(h) || ...
            numel(X) ~= N || numel(Y) ~= N || numel(h) ~= N
        error('plotTFL:BadXYH', ...
            'X, Y, h must be vectors of length N = size(W,1).');
    end

    X = X(:); Y = Y(:); h = h(:);

    if any(~isfinite(X)) || any(~isfinite(Y)) || any(~isfinite(h))
        error('plotTFL:NonFiniteXYH', ...
            'X, Y, and h must be finite.');
    end

    % ---- parse opts ----
    opts = parsePlotOpts_(varargin{:});

    opts.ArrowLayer = validateArrowLayer_(opts.ArrowLayer);

    opts.ArrowColorMode = validateArrowColorMode_(opts.ArrowColorMode);

    if strcmpi(opts.EdgeWidthMode, 'auto')
        if any(strcmpi(opts.RenderMode, {'overlay','cdfdcolor','colorcdfd','color'}))
            opts.EdgeWidthMode = 'weight';
        else
            opts.EdgeWidthMode = 'fixed';
        end
    end

    % ---- axes setup ----
    ax = opts.Parent;
    if isempty(ax) || ~isgraphics(ax, 'axes')
        ax = gca;
    end

    holdState = ishold(ax);
    if opts.ClearAxes
        cla(ax);
    end
    hold(ax, 'on');
    ax.SortMethod = 'childorder';

    % ---- shared scene + shared geometry ----
    scene = computeTFLScene(X, Y, h, opts);
    geom  = computeTFLEdgeGeometry(W, scene, opts);

    % Expand suggested limits to include the actual realised edge geometry
    %scene = expandSceneLimitsToGeometry_(scene, geom, opts);
    scene = expandSceneLimitsToGeometry_(scene, geom);

    Xplot = scene.Xplot;
    Yplot = scene.Yplot;

    % ---- bands ----
    if opts.ShowBands && strcmpi(opts.LevelSemantics,'trophic')
        for k = 1:numel(scene.bandsY)
            yk = scene.bandsY(k);
            line(ax, scene.bandXLimits, [yk yk], ...
                'Color', opts.BandColor, ...
                'LineStyle', opts.BandStyle, ...
                'LineWidth', opts.BandWidth);
        end
    end

    
    %renderItems = struct([]);

    switch lower(strtrim(opts.RenderMode))
        case 'plain'
            renderItems = buildPlainRenderItems_(scene, geom, opts);

        case 'overlay'
            [geom, scene, renderItems] = ...
                buildOverlayRenderItems_(W, scene, geom, opts);

        case {'cdfdcolor','colorcdfd','color'}
            [geom, scene, renderItems] = ...
                buildCDFDColorRenderItems_(W, scene, geom, opts);

        otherwise
            error('plotTFL:BadRenderMode', ...
                'RenderMode must be ''plain'', ''overlay'', or ''cdfdcolor''.');
    end


    % ---- axes limits / formatting ----
    xlim(ax, scene.xlimSuggested);
    ylim(ax, scene.ylimSuggested);

    set(ax, 'XTick', [], 'YTick', []);
    set(ax, 'XTickLabel', [], 'YTickLabel', []);
    box(ax, 'on');

    nodeSize = scene.nodeSize;

    drawnow;

    if opts.ShowNodes && opts.ArrowOverNodes
        nodeRadiusData = markerRadiusData_(ax, nodeSize, opts.ArrowNodeGapPts);
    else
        nodeRadiusData = [];
    end

    % ---- draw shafts + arrows together ----
    drawRenderItems_(ax, renderItems, nodeRadiusData, opts);

    % ---- nodes last ----
    if opts.ShowNodes
        scatter(ax, Xplot, Yplot, nodeSize, ...
            'MarkerFaceColor', opts.NodeColor, ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 0.5);
    end

    % ---- labels ----
    labelHs = gobjects(0,1);

    if opts.ShowLabels
        if isempty(opts.Labels)
            lbls = arrayfun(@(k) sprintf('%d',k), 1:N, 'UniformOutput', false);
        else
            lbls = opts.Labels;
            if numel(lbls) ~= N
                error('plotTFL:BadLabels', ...
                    'Labels must be a cell array of N strings.');
            end
        end

        dx = opts.LabelOffset(1);
        dy = opts.LabelOffset(2);

        if strcmpi(opts.LabelOffsetMode, 'frac')
            dx = dx * rangeNonzero_(Xplot);
            dy = dy * rangeNonzero_(Yplot);
        end

        anchorXY = zeros(0,2);
        for ii = 1:N
            s = lbls{ii};
            if isempty(s) || (isstring(s) && strlength(s)==0)
                continue;
            end

            t = text(ax, Xplot(ii)+dx, Yplot(ii)+dy, s, ...
                'Units', 'data', ...
                'FontSize', opts.LabelFontSize, ...
                'Color', opts.LabelColor, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom');

            labelHs(end+1,1) = t; %#ok<AGROW>
            anchorXY(end+1,:) = [Xplot(ii)+dx, Yplot(ii)+dy]; %#ok<AGROW>
        end

        if opts.AdjustLabels && ~isempty(labelHs)
            adjustTextLabels_(ax, labelHs, anchorXY, opts);
        end
    end

    if opts.LabelHalo && ~isempty(labelHs)
        applyTextHalos_(ax, labelHs, opts.LabelHaloPadFrac, ...
            opts.LabelHaloAlpha, opts.LabelHaloColor);
    end

    % ---- title / BFF label ----
    if ~isempty(opts.Title)
        title(ax, opts.Title, 'Interpreter', 'none');
    end

    if opts.ShowBFF && ~isnan(opts.BFF)
        xl = xlim(ax);
        yl = ylim(ax);
        xText = xl(1) + 0.02 * (xl(2) - xl(1));
        yText = yl(2) - 0.05 * (yl(2) - yl(1));
        text(ax, xText, yText, sprintf('BFF = %.3f', opts.BFF), ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment',   'top', ...
            'FontSize',            10, ...
            'BackgroundColor',     'w', ...
            'Margin',              2);
    end

    if ~holdState
        hold(ax, 'off');
    end
end

% =====================================================================
function opts = parsePlotOpts_(varargin)

    opts = struct();

    opts.Parent = [];
    opts.ClearAxes = true;

    opts.RenderMode = 'plain';   % 'plain' | 'overlay' | 'cdfdcolor' (aliases: 'colorcdfd', 'color')
    opts.LevelSemantics = 'generic';   % 'trophic' | 'generic'

    opts.C = [];
    opts.D = [];
    opts.Info = [];

    opts.ComputeIfMissing = true;
    opts.DecompositionMethod = 'bff';
    opts.DecompositionArgs = struct( ...
        'ToleranceZero', 1e-10, ...
        'Validate', true);

    opts.OverlayComponent = 'c';      % 'c' or 'd'
    opts.OverlayBaseColorMode = 'plain';   % 'component' | 'plain'

    opts.CircularColor    = [1 0 0];
    opts.CircularAlpha    = 0.90;

    opts.DirectionalColor = [0 0.35 0.9]; %[0.45 0.45 0.45]; % [0 0.35 0.9]
    opts.DirectionalAlpha = 0.78;

    opts.ArrowLayer = 'fullflow';    % 'fullflow' | 'overlay' | 'both'
    opts.ArrowColorMode = 'blend';   % 'shaft' | 'overlay' | 'blend'

    opts.FillSquare = true;
    opts.YScale     = 1;
    opts.XScale     = 1;

    opts.PadFracX      = 0.05;
    opts.PadFracY      = 0.04;
    opts.CurvPadScaleX = 1.00;
    opts.CurvPadScaleY = 0.25;

    opts.ShowBands = true;
    opts.BandTau   = 1;
    opts.BandColor = [0.75 0.75 0.75];
    opts.BandStyle = ':';
    opts.BandWidth = 0.75;

    opts.ShowEdges = true;
    opts.ShowNodes = true;

    opts.NodeColor     = [0.1 0.3 0.7];
    opts.NodeSize      = 30;
    opts.NodeSizeMode  = 'auto';
    opts.NodeSizeRefN  = 20;
    opts.NodeSizeGamma = 0.3;
    opts.NodeSizeMin   = 10;
    opts.NodeSizeMax   = 80;
    opts.NodeSizeData      = [];
    opts.NodeSizeScaleData = [];
    opts.NodeSizeScaleMode = 'area';
    opts.NodeSizeRange     = [20 180];


    opts.UpEdgeColor   = [0.35 0.35 0.35];
    opts.DownEdgeColor = [0.68 0.68 0.68];
    opts.EdgeAlpha     = 0.78;
    opts.EdgeWidth     = 0.75;

    opts.ArrowLenBasis = 'span';
    opts.ArrowSize     = 0.03;
    opts.ArrowMinFrac  = 0.02;
    opts.ArrowMaxFrac  = 0.6;
    opts.ArrowAngle    = pi/7;
    opts.ArrowOverNodes = true;
    opts.ArrowNodeGapPts = 1.5;

    opts.ArrowSizeMode      = 'withNode';
    opts.ArrowSizeRefN      = 40;
    opts.ArrowSizeGamma     = 0.5;
    opts.ArrowSizeMinScale  = 0.35;
    opts.ArrowSizeMaxScale  = 1.2;

    opts.EdgeWidthScaleMode = 'withNode';
    opts.EdgeWidthRefN      = 40;
    opts.EdgeWidthGamma     = 0.35;
    opts.EdgeWidthMinScale  = 0.5;
    opts.EdgeWidthMaxScale  = 1.2;
    opts.EdgeWidthMin       = 0.25;
    opts.EdgeWidthMax       = 3.0;

    opts.ArrowScaleWithEdgeWidth = true;
    opts.ArrowEdgeGamma    = 0.7;
    opts.ArrowEdgeMinScale = 0.75;
    opts.ArrowEdgeMaxScale = 1.8;

    opts.ShowLabels    = false;
    opts.Labels        = {};
    opts.LabelFontSize = 9;
    opts.LabelColor    = [0 0 0];
    opts.LabelOffset   = [0, 0.03];
    opts.LabelOffsetMode = 'frac';

    opts.AdjustLabels       = false;
    opts.LabelRepelMode     = 'vertical';
    opts.LabelRepelIters    = 10;
    opts.LabelRepelStepFrac = 0.012;
    opts.LabelRepelPadFrac  = 0.003;
    opts.LabelMaxShiftFrac  = 0.06;

    opts.LabelHalo        = false;
    opts.LabelHaloPadFrac = 0.006;
    opts.LabelHaloAlpha   = 0.30;
    opts.LabelHaloColor   = [1 1 1];

    opts.Title   = '';
    opts.ShowBFF = false;
    opts.BFF     = NaN;

    opts.HitNode            = true;
    opts.HitNodeRadiusFrac  = 0.030;
    opts.HitNodeMarginFrac  = 0.005;
    opts.HitNodeCurvFrac    = 0.18;
    opts.HitNodeNpts        = 25;
    opts.HitNodeMinLenFrac  = 0.15;

    opts.LaneSplit        = true;
    opts.LaneSplitScaleWithEdgeLength = true;
    opts.LaneSplitTheta   = deg2rad(6);
    opts.LaneSplitSpacing = 0.15;
    opts.LaneSplitAlpha   = 0.25;
    opts.LaneSplitNpts    = 25;
    opts.LaneSplitMinLen  = 0;
    opts.LaneSplitCoherent = true;
    opts.LaneSplitCoherentTol = 0.25;
    opts.LaneSplitPreferDownSide = true;
    opts.LaneSplitPrimaryStraight = true;
    opts.LaneSplitPrimaryRule     = 'upweight';
    opts.LaneSplitOneSided        = true;
    opts.LaneSplitChainStraight   = true;

    opts.EdgeWidthMode      = 'auto';    % 'auto' | 'fixed' | 'weight'
    opts.EdgeWidthRange     = [0.6 2.5];
    opts.EdgeWidthScale     = 'log';
    opts.EdgeWidthQuantiles = [0.05 0.95];
    %opts.EdgeWidthMode      = 'fixed';
 

    if isempty(varargin)
        return;
    end
    if mod(numel(varargin),2) ~= 0
        error('plotTFL:BadArgs', ...
            'Optional arguments must be name/value pairs.');
    end

    for k = 1:2:numel(varargin)
        name  = varargin{k};
        value = varargin{k+1};

        if ~ischar(name) && ~isstring(name)
            error('plotTFL:BadParamName', ...
                'Parameter names must be strings.');
        end

        switch lower(char(name))
            case 'parent', opts.Parent = value;
            case 'clearaxes', opts.ClearAxes = logical(value);

            case 'rendermode'
                opts.RenderMode = char(value);

            case 'levelsemantics'
                opts.LevelSemantics = char(value);

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

            case 'overlaycomponent'
                opts.OverlayComponent = char(value);
            case 'overlaybasecolormode'
                opts.OverlayBaseColorMode = char(value);

            case 'circularcolor'
                opts.CircularColor = value;
            case 'circularalpha'
                opts.CircularAlpha = value;

            case 'directionalcolor'
                opts.DirectionalColor = value;
            case 'directionalalpha'
                opts.DirectionalAlpha = value;

            case 'arrowlayer'
                opts.ArrowLayer = char(value);
            case 'arrowcolormode'
                opts.ArrowColorMode = char(value);

            case 'fillsquare', opts.FillSquare = logical(value);
            case 'yscale', opts.YScale = value;
            case 'xscale', opts.XScale = value;
            case 'xstretch', opts.XScale = value;

            case 'padfracx', opts.PadFracX = value;
            case 'padfracy', opts.PadFracY = value;
            case 'curvpadscalex', opts.CurvPadScaleX = value;
            case 'curvpadscaley', opts.CurvPadScaleY = value;

            case 'showbands', opts.ShowBands = logical(value);
            case 'bandtau', opts.BandTau = value;
            case 'bandcolor', opts.BandColor = value;
            case 'bandstyle', opts.BandStyle = value;
            case 'bandwidth', opts.BandWidth = value;

            case 'showedges', opts.ShowEdges = logical(value);
            case 'shownodes', opts.ShowNodes = logical(value);

            case 'nodecolor', opts.NodeColor = value;
            case 'nodesize', opts.NodeSize = value;
            case 'nodesizemode', opts.NodeSizeMode = char(value);
            case 'nodesizerefn', opts.NodeSizeRefN = value;
            case 'nodesizegamma', opts.NodeSizeGamma = value;
            case 'nodesizemin', opts.NodeSizeMin = value;
            case 'nodesizemax', opts.NodeSizeMax = value;
            case 'nodesizedata', opts.NodeSizeData = value;
            case 'nodesizescaledata', opts.NodeSizeScaleData = value;
            case 'nodesizescalemode', opts.NodeSizeScaleMode = char(value);
            case 'nodesizerange', opts.NodeSizeRange = value;

            case 'upedgecolor', opts.UpEdgeColor = value;
            case 'downedgecolor', opts.DownEdgeColor = value;
            case 'edgewidth', opts.EdgeWidth = value;
            case 'edgealpha', opts.EdgeAlpha = value;

            case 'arrowlenbasis', opts.ArrowLenBasis = char(value);
            case 'arrowsize', opts.ArrowSize = value;
            case 'arrowangle', opts.ArrowAngle = value;
            case 'arrowminfrac', opts.ArrowMinFrac = value;
            case 'arrowmaxfrac', opts.ArrowMaxFrac = value;
            case 'arrowovernodes', opts.ArrowOverNodes = logical(value);
            case 'arrownodegappts', opts.ArrowNodeGapPts = value;

            case 'arrowsizemode', opts.ArrowSizeMode = char(value);
            case 'arrowsizerefn', opts.ArrowSizeRefN = value;
            case 'arrowsizegamma', opts.ArrowSizeGamma = value;
            case 'arrowsizeminscale', opts.ArrowSizeMinScale = value;
            case 'arrowsizemaxscale', opts.ArrowSizeMaxScale = value;

            case 'edgewidthscalemode', opts.EdgeWidthScaleMode = char(value);
            case 'edgewidthrefn', opts.EdgeWidthRefN = value;
            case 'edgewidthgamma', opts.EdgeWidthGamma = value;
            case 'edgewidthminscale', opts.EdgeWidthMinScale = value;
            case 'edgewidthmaxscale', opts.EdgeWidthMaxScale = value;
            case 'edgewidthmin', opts.EdgeWidthMin = value;
            case 'edgewidthmax', opts.EdgeWidthMax = value;

            case 'arrowscalewithedgewidth', opts.ArrowScaleWithEdgeWidth = logical(value);
            case 'arrowedgegamma', opts.ArrowEdgeGamma = value;
            case 'arrowedgeminscale', opts.ArrowEdgeMinScale = value;
            case 'arrowedgemaxscale', opts.ArrowEdgeMaxScale = value;

            case 'showlabels', opts.ShowLabels = logical(value);
            case 'labels', opts.Labels = value;
            case 'labelfontsize', opts.LabelFontSize = value;
            case 'labelcolor', opts.LabelColor = value;
            case 'labeloffset', opts.LabelOffset = value;
            case 'labeloffsetmode', opts.LabelOffsetMode = char(value);

            case 'adjustlabels', opts.AdjustLabels = logical(value);
            case 'labelrepelmode', opts.LabelRepelMode = char(value);
            case 'labelrepeliters', opts.LabelRepelIters = value;
            case 'labelrepelstepfrac', opts.LabelRepelStepFrac = value;
            case 'labelrepelpadfrac', opts.LabelRepelPadFrac = value;
            case 'labelmaxshiftfrac', opts.LabelMaxShiftFrac = value;

            case 'labelhalo', opts.LabelHalo = logical(value);
            case 'labelhalopadfrac', opts.LabelHaloPadFrac = value;
            case 'labelhaloalpha', opts.LabelHaloAlpha = value;
            case 'labelhalocolor', opts.LabelHaloColor = value;

            case 'title', opts.Title = value;
            case 'showbff', opts.ShowBFF = logical(value);
            case 'bff', opts.BFF = value;

            case 'hitnode', opts.HitNode = logical(value);
            case 'hitnoderadiusfrac', opts.HitNodeRadiusFrac = value;
            case 'hitnodemarginfrac', opts.HitNodeMarginFrac = value;
            case 'hitnodecurvfrac', opts.HitNodeCurvFrac = value;
            case 'hitnodenpts', opts.HitNodeNpts = value;
            case 'hitnodeminlenfrac', opts.HitNodeMinLenFrac = value;

            case 'lanesplit', opts.LaneSplit = logical(value);
            case 'lanesplitscalewithedgelength', opts.LaneSplitScaleWithEdgeLength = logical(value);
            case 'lanesplittheta', opts.LaneSplitTheta = value;
            case 'lanesplitspacing', opts.LaneSplitSpacing = value;
            case 'lanesplitalpha', opts.LaneSplitAlpha = value;
            case 'lanesplitnpts', opts.LaneSplitNpts = value;
            case 'lanesplitminlen', opts.LaneSplitMinLen = value;
            case 'lanesplitcoherent', opts.LaneSplitCoherent = logical(value);
            case 'lanesplitcoherenttol', opts.LaneSplitCoherentTol = value;
            case 'lanesplitpreferdownside', opts.LaneSplitPreferDownSide = logical(value);
            case 'lanesplitprimarystraight', opts.LaneSplitPrimaryStraight = logical(value);
            case 'lanesplitprimaryrule', opts.LaneSplitPrimaryRule = char(value);
            case 'lanesplitonesided', opts.LaneSplitOneSided = logical(value);
            case 'lanesplitchainstraight', opts.LaneSplitChainStraight = logical(value);

            case 'edgewidthmode', opts.EdgeWidthMode = char(value);
            case 'edgewidthrange', opts.EdgeWidthRange = value;
            case 'edgewidthscale', opts.EdgeWidthScale = char(value);
            case 'edgewidthquantiles', opts.EdgeWidthQuantiles = value;

            otherwise
                error('plotTFL:UnknownOption', ...
                    'Unknown option: %s', char(name));
        end
    end
end

function layer = validateArrowLayer_(layer)
%VALIDATEARROWLAYER_  Validate ArrowLayer for RenderMode='overlay'.

    layer = lower(strtrim(char(layer)));

    switch layer
        case {'fullflow','overlay','both'}
            % ok
        otherwise
            error('plotTFL:BadArrowLayer', ...
                'ArrowLayer must be ''fullflow'', ''overlay'', or ''both''.');
    end
end

function mode = validateArrowColorMode_(mode)
%VALIDATEARROWCOLORMODE_  Validate ArrowColorMode.

    mode = lower(strtrim(char(mode)));

    switch mode
        case {'shaft','overlay','blend'}
            % ok
        otherwise
            error('plotTFL:BadArrowColorMode', ...
                'ArrowColorMode must be ''shaft'', ''overlay'', or ''blend''.');
    end
end

% =====================================================================
function lw = mapWeightsToLineWidth_(wRaw, lwRange, scaleMode, qRange, fallbackLW)

    w = double(wRaw(:));

    lw = fallbackLW * ones(size(w));

    if isempty(w) || all(~isfinite(w)) || all(w <= 0)
        return;
    end

    wpos = w(isfinite(w) & w > 0);
    if isempty(wpos), return; end

    qLo = max(0, min(1, qRange(1)));
    qHi = max(0, min(1, qRange(2)));
    if qHi <= qLo
        qLo = 0.05;
        qHi = 0.95;
    end

    wSorted = sort(wpos);
    n = numel(wSorted);
    wLo = wSorted(max(1, round(qLo*(n-1))+1));
    wHi = wSorted(max(1, round(qHi*(n-1))+1));

    wClip = min(max(w, wLo), wHi);

    switch lower(strtrim(scaleMode))
        case 'log'
            x = log10(max(wClip, eps));
            xLo = log10(max(wLo, eps));
            xHi = log10(max(wHi, eps));
        otherwise
            x = wClip;
            xLo = wLo;
            xHi = wHi;
    end

    if ~isfinite(xLo) || ~isfinite(xHi) || xHi <= xLo
        lw = mean(lwRange) * ones(size(w));
        return;
    end

    t = (x - xLo) ./ (xHi - xLo);
    t = min(max(t,0),1);

    lwMin = lwRange(1);
    lwMax = lwRange(2);
    lw = lwMin + t .* (lwMax - lwMin);
end


% =====================================================================
function m = localRobustMedianPositive_(v, tol)
    v = v(:);
    v = v(isfinite(v) & v > tol);
    if isempty(v)
        m = NaN;
    else
        m = median(v);
    end
end

% =====================================================================
function r = markerRadiusData_(ax, markerAreaPts2, gapPts)

    s = double(markerAreaPts2(:));
    s = max(s, 0);

    rPts = sqrt(s / pi) + gapPts;

    ppi = get(groot, 'ScreenPixelsPerInch');
    rPx = rPts * ppi / 72;

    axPix = getpixelposition(ax, true);
    xl = xlim(ax);
    yl = ylim(ax);

    rx = rPx * diff(xl) / max(axPix(3), eps);
    ry = rPx * diff(yl) / max(axPix(4), eps);

    r = 0.5 * (rx + ry);
end

% =====================================================================
function drawArrowheadToNode_(ax, x1, y1, x2, y2, col, L, halfAngle, nodeRadius)

    dx = x2 - x1;
    dy = y2 - y1;
    len = hypot(dx, dy);
    if len <= 0
        return;
    end

    ux = dx / len;
    uy = dy / len;

    tipBack = max(nodeRadius, 0);
    avail = len - tipBack;
    if avail <= 0
        return;
    end

    Luse = min(L, 0.8 * avail);
    if Luse <= 0
        return;
    end

    tx = x2 - tipBack * ux;
    ty = y2 - tipBack * uy;

    px = tx - Luse * ux;
    py = ty - Luse * uy;

    px_perp = -uy;
    py_perp =  ux;

    w = tan(halfAngle) * Luse;

    bx1 = px + w * px_perp;
    by1 = py + w * py_perp;
    bx2 = px - w * px_perp;
    by2 = py - w * py_perp;

    edgeCol = darkenColor_(col, 0.28);

    patch(ax, [tx bx1 bx2], [ty by1 by2], col, ...
        'EdgeColor', edgeCol, ...
        'LineWidth', 0.45, ...
        'FaceColor', col, ...
        'HitTest', 'off', ...
        'PickableParts', 'none');
end

% =====================================================================
function adjustTextLabels_(ax, textHs, anchorXY, opts)

    if isempty(textHs)
        return;
    end

    drawnow;

    spanX = rangeNonzero_(xlim(ax));
    spanY = rangeNonzero_(ylim(ax));

    padX  = opts.LabelRepelPadFrac  * spanX;
    padY  = opts.LabelRepelPadFrac  * spanY;
    stepX = opts.LabelRepelStepFrac * spanX;
    stepY = opts.LabelRepelStepFrac * spanY;
    maxDx = opts.LabelMaxShiftFrac  * spanX;
    maxDy = opts.LabelMaxShiftFrac  * spanY;

    M = numel(textHs);

    for it = 1:opts.LabelRepelIters
        pos = getTextPositions_(textHs);
        ext = getTextExtentsData_(textHs);

        dpos = zeros(M,2);
        nOverlap = 0;

        for a = 1:M-1
            for b = a+1:M
                [tf, ovx, ovy] = boxesOverlap_(ext(a,:), ext(b,:), padX, padY);
                if ~tf
                    continue;
                end

                nOverlap = nOverlap + 1;

                switch lower(strtrim(opts.LabelRepelMode))
                    case 'xy'
                        if ovx < ovy
                            sgn = sign(pos(b,1) - pos(a,1));
                            if sgn == 0
                                sgn = sign(anchorXY(b,1) - anchorXY(a,1));
                            end
                            if sgn == 0, sgn = 1; end

                            mag = min(stepX, 0.55*ovx + padX);
                            dpos(a,1) = dpos(a,1) - 0.5 * sgn * mag;
                            dpos(b,1) = dpos(b,1) + 0.5 * sgn * mag;
                        else
                            sgn = sign(pos(b,2) - pos(a,2));
                            if sgn == 0
                                sgn = sign(anchorXY(b,2) - anchorXY(a,2));
                            end
                            if sgn == 0
                                sgn = sign(pos(b,1) - pos(a,1));
                            end
                            if sgn == 0, sgn = 1; end

                            mag = min(stepY, 0.55*ovy + padY);
                            dpos(a,2) = dpos(a,2) - 0.5 * sgn * mag;
                            dpos(b,2) = dpos(b,2) + 0.5 * sgn * mag;
                        end

                    otherwise
                        sgn = sign(pos(b,2) - pos(a,2));
                        if sgn == 0
                            sgn = sign(anchorXY(b,2) - anchorXY(a,2));
                        end
                        if sgn == 0
                            sgn = sign(pos(b,1) - pos(a,1));
                        end
                        if sgn == 0, sgn = 1; end

                        mag = min(stepY, 0.55*ovy + padY);
                        dpos(a,2) = dpos(a,2) - 0.5 * sgn * mag;
                        dpos(b,2) = dpos(b,2) + 0.5 * sgn * mag;
                end
            end
        end

        if nOverlap == 0
            break;
        end

        for m = 1:M
            newPos = pos(m,:) + dpos(m,:);
            delta  = newPos - anchorXY(m,:);

            delta(1) = min(max(delta(1), -maxDx), maxDx);
            delta(2) = min(max(delta(2), -maxDy), maxDy);

            p = textHs(m).Position;
            p(1:2) = anchorXY(m,:) + delta;
            textHs(m).Position = p;
        end

        drawnow limitrate nocallbacks;
    end
end

% =====================================================================
function pos = getTextPositions_(textHs)
    M = numel(textHs);
    pos = zeros(M,2);
    for m = 1:M
        p = textHs(m).Position;
        pos(m,:) = p(1:2);
    end
end

% =====================================================================
function ext = getTextExtentsData_(textHs)
    M = numel(textHs);
    ext = zeros(M,4);
    for m = 1:M
        oldUnits = textHs(m).Units;
        textHs(m).Units = 'data';
        ext(m,:) = textHs(m).Extent;
        textHs(m).Units = oldUnits;
    end
end

% =====================================================================
function [tf, ovx, ovy] = boxesOverlap_(A, B, padX, padY)

    ax1 = A(1) - padX;
    ax2 = A(1) + A(3) + padX;
    ay1 = A(2) - padY;
    ay2 = A(2) + A(4) + padY;

    bx1 = B(1) - padX;
    bx2 = B(1) + B(3) + padX;
    by1 = B(2) - padY;
    by2 = B(2) + B(4) + padY;

    ovx = min(ax2, bx2) - max(ax1, bx1);
    ovy = min(ay2, by2) - max(ay1, by1);

    tf = (ovx > 0) && (ovy > 0);
end

% =====================================================================
function applyTextHalos_(ax, labelHs, padFrac, faceAlpha, faceColor)

    if nargin < 3 || isempty(padFrac),   padFrac = 0.006; end
    if nargin < 4 || isempty(faceAlpha), faceAlpha = 0.30; end
    if nargin < 5 || isempty(faceColor), faceColor = [1 1 1]; end

    labelHs = labelHs(isgraphics(labelHs, 'text'));
    if isempty(labelHs)
        return;
    end

    drawnow;
    ax.SortMethod = 'childorder';

    xl = xlim(ax);
    yl = ylim(ax);
    padX = padFrac * (xl(2) - xl(1));
    padY = padFrac * (yl(2) - yl(1));

    haloHs = gobjects(0,1);

    for k = 1:numel(labelHs)
        t = labelHs(k);

        s = t.String;
        if (ischar(s) && isempty(s)) || (isstring(s) && all(strlength(s)==0))
            continue;
        end

        oldUnits = t.Units;
        t.Units = 'data';
        ext = t.Extent;
        t.Units = oldUnits;

        if numel(ext) ~= 4 || ext(3) <= 0 || ext(4) <= 0
            continue;
        end

        x0 = ext(1) - padX;
        x1 = ext(1) + ext(3) + padX;
        y0 = ext(2) - padY;
        y1 = ext(2) + ext(4) + padY;

        h = patch(ax, [x0 x1 x1 x0], [y0 y0 y1 y1], faceColor, ...
            'EdgeColor', 'none', ...
            'FaceAlpha', faceAlpha, ...
            'HitTest', 'off', ...
            'PickableParts', 'none', ...
            'HandleVisibility', 'off', ...
            'Tag', 'LabelHalo');

        haloHs(end+1,1) = h; %#ok<AGROW>
    end

    if isempty(haloHs)
        return;
    end

    allChildren = ax.Children;
    isLabel = ismember(allChildren, labelHs);
    isHalo  = ismember(allChildren, haloHs);
    others  = allChildren(~(isLabel | isHalo));

    ax.Children = [labelHs(:); haloHs(:); others(:)];
end

% =====================================================================
function c2 = darkenColor_(c, frac)
    c = c(:)';
    c2 = max(0, c * (1 - frac));
end

% =====================================================================
function cBlend = blendColor_(c, alpha, bg)
    c  = c(:)';
    bg = bg(:)';
    cBlend = alpha * c + (1-alpha) * bg;
end

% =====================================================================
function r = rangeNonzero_(x)
    r = max(x) - min(x);
    if r == 0
        r = 1;
    end
end


function [C, D, info] = resolveCDFD_(W, opts)

    if isempty(opts.C) || isempty(opts.D)
        if ~opts.ComputeIfMissing
            error('plotTFL:MissingCD', ...
                'C and D were not supplied and ComputeIfMissing is false.');
        end

        switch lower(strtrim(opts.DecompositionMethod))
            case 'bff'
                cdfdArgs = structToNameValue_(opts.DecompositionArgs);
                [C, D, info] = cdfd_bff(W, cdfdArgs{:});
            otherwise
                error('plotTFL:BadDecompositionMethod', ...
                    'Unknown DecompositionMethod: %s', opts.DecompositionMethod);
        end
    else
        C = sparse(double(opts.C));
        D = sparse(double(opts.D));

        if ~isequal(size(C), size(W)) || ~isequal(size(D), size(W))
            error('plotTFL:BadCDSize', ...
                'C and D must have the same size as W.');
        end

        if isempty(opts.Info)
            info = struct();
        else
            info = opts.Info;
        end
    end

    info = fillCDFDInfo_(W, C, D, info);
end

function geom = attachEdgeVisualData_(geom, Wraw, Craw, Draw)

    edges = geom.edges;
    if isempty(edges)
        return;
    end

    lin = vertcat(edges.lin);

    wRawW = full(Wraw(lin));
    wRawC = full(Craw(lin));
    wRawD = full(Draw(lin));

    rhoC = safeRatio_(wRawC, wRawW);
    rhoD = safeRatio_(wRawD, wRawW);

    for e = 1:numel(edges)
        edges(e).wRawW = wRawW(e);
        edges(e).wRawC = wRawC(e);
        edges(e).wRawD = wRawD(e);

        edges(e).rhoC = rhoC(e);
        edges(e).rhoD = rhoD(e);

        % Displayed widths are assigned later from the total-flow edge width
        edges(e).lwW = NaN;
        edges(e).lwC = NaN;
        edges(e).lwD = NaN;
    end

    geom.edges = edges;
end

%{
function [arrowX1, arrowY1, arrowX2, arrowY2, arrowTo, arrowLenRec, arrowCol] = ...
    appendArrowRecordFromEdge_(arrowX1, arrowY1, arrowX2, arrowY2, ...
    arrowTo, arrowLenRec, arrowCol, edge, L, col)

    if ~isfinite(edge.arrowX1) || ~isfinite(edge.arrowY1) || ...
       ~isfinite(edge.arrowX2) || ~isfinite(edge.arrowY2)
        return;
    end

    arrowX1(end+1,1) = edge.arrowX1;
    arrowY1(end+1,1) = edge.arrowY1;
    arrowX2(end+1,1) = edge.arrowX2;
    arrowY2(end+1,1) = edge.arrowY2;
    arrowTo(end+1,1) = edge.toNode;
    arrowLenRec(end+1,1) = L;
    arrowCol(end+1,:) = col(:).';
end

function tf = edgeFieldPositive_(edges, fieldName)
    if isempty(edges)
        tf = false(0,1);
        return;
    end
    tf = vertcat(edges.(fieldName)) > 0;
end

%}

%{
function vals = edgeFieldValues_(edges, fieldName, idx)
    if nargin < 3 || isempty(idx)
        idx = 1:numel(edges);
    end
    if isempty(idx)
        vals = zeros(0,1);
        return;
    end
    vals = vertcat(edges(idx).(fieldName));
end
%}

function col = resolveArrowColor_(edge, opts, colShaft, colOverlay)

    switch lower(strtrim(opts.ArrowColorMode))
        case 'shaft'
            col = colShaft;

        case 'overlay'
            col = colOverlay;

        case 'blend'
            switch lower(strtrim(opts.OverlayComponent))
                case 'c'
                    a = edge.rhoC;
                case 'd'
                    a = edge.rhoD;
                otherwise
                    a = 0.5;
            end

            a = min(max(a, 0), 1);
            col = (1 - a) * colShaft + a * colOverlay;

        otherwise
            error('plotTFL:BadArrowColorMode', ...
                'ArrowColorMode must be ''blend'', ''shaft'', or ''overlay''.');
    end
end


function r = safeRatio_(num, den)
    num = double(num(:));
    den = double(den(:));

    r = zeros(size(num));
    mask = (den > 0) & isfinite(den) & isfinite(num);
    r(mask) = num(mask) ./ den(mask);

    r = min(max(r, 0), 1);
end

function info = fillCDFDInfo_(W, C, D, info)
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

function nv = structToNameValue_(s)
    f = fieldnames(s);
    nv = cell(1, 2*numel(f));
    for ii = 1:numel(f)
        nv{2*ii-1} = f{ii};
        nv{2*ii}   = s.(f{ii});
    end
end

function renderItems = buildPlainRenderItems_(scene, geom, opts)
%BUILDPLAINRENDERITEMS_  Build sorted shaft+arrow render items for plain mode.

    renderItems = struct([]);

    E = numel(geom.edges);
    if ~opts.ShowEdges || E == 0
        return;
    end

    % ---- line widths ----
    edgeLW = computePlainEdgeLineWidths_(geom, scene, opts);

    % ---- arrow lengths ----
    arrowLenE = computeArrowLengths_(scene, geom, edgeLW, opts);

    % ---- colours ----
    bg = [1 1 1];
    colUp   = blendColor_(opts.UpEdgeColor,   opts.EdgeAlpha, bg);
    colDown = blendColor_(opts.DownEdgeColor, opts.EdgeAlpha, bg);

    % ---- build items ----
    for e = 1:E
        edge = geom.edges(e);

        if edge.isUp
            col = colUp;
        else
            col = colDown;
        end

        item = makeRenderItem_( ...
            edge, ...
            edgeLW(e), ...
            col, ...
            true, ...
            arrowLenE(e), ...
            col);

        renderItems = appendRenderItem_(renderItems, item);
    end

    % ---- largest first ----
    renderItems = sortRenderItems_(renderItems);
end


function [geom, scene, renderItems] = buildCDFDColorRenderItems_(W, scene, geom, opts)
%BUILDCDFDCOLORRENDERITEMS_  Build sorted shaft+arrow render items for cdfdcolor mode.

    renderItems = struct([]);

    if ~opts.ShowEdges || isempty(geom.edges)
        return;
    end

    % ---- resolve decomposition ----
    [C, D, info] = resolveCDFD_(W, opts);
    scene.C = C;
    scene.D = D;
    scene.CDFDInfo = info;

    % ---- attach raw per-edge circular/directional shares ----
    geom = attachEdgeVisualData_(geom, W, C, D);

    % ---- full-width shaft width should match plain rendering ----
    plainLW = computePlainEdgeLineWidths_(geom, scene, opts);
    for e = 1:numel(geom.edges)
        geom.edges(e).lwW = plainLW(e);
    end

    edges = geom.edges;

    % ---- arrow lengths from full displayed widths ----
    arrowLenW = computeArrowLengths_(scene, geom, vertcat(edges.lwW), opts);

    % ---- component colours ----
    colCirc = blendColor_(opts.CircularColor,    opts.CircularAlpha,    [1 1 1]);
    colDir  = blendColor_(opts.DirectionalColor, opts.DirectionalAlpha, [1 1 1]);

    % ---- build items ----
    for e = 1:numel(edges)
        edge = edges(e);

        if edge.lwW <= 0
            continue;
        end

        % Full-width colour blend from raw edgewise shares
        col = edge.rhoC * colCirc + edge.rhoD * colDir;

        item = makeRenderItem_( ...
            edge, ...
            edge.lwW, ...
            col, ...
            true, ...
            arrowLenW(e), ...
            col);

        renderItems = appendRenderItem_(renderItems, item);
    end

    % ---- largest first ----
    renderItems = sortRenderItems_(renderItems);
end

function [geom, scene, renderItems] = buildOverlayRenderItems_(W, scene, geom, opts)
%BUILDOVERLAYRENDERITEMS_  Build sorted shaft+arrow render items for overlay mode.

    renderItems = struct([]);

    if ~opts.ShowEdges || isempty(geom.edges)
        return;
    end

    % ---- resolve decomposition ----
    [C, D, info] = resolveCDFD_(W, opts);
    scene.C = C;
    scene.D = D;
    scene.CDFDInfo = info;

    % ---- attach per-edge visual data by edge identity ----
    geom = attachEdgeVisualData_(geom, W, C, D);

    plainLW = computePlainEdgeLineWidths_(geom, scene, opts);

    rhoC = vertcat(geom.edges.rhoC);
    rhoD = vertcat(geom.edges.rhoD);

    lwC = rhoC .* plainLW;
    lwD = rhoD .* plainLW;

    for e = 1:numel(geom.edges)
        geom.edges(e).lwW = plainLW(e);
        geom.edges(e).lwC = lwC(e);
        geom.edges(e).lwD = lwD(e);
    end

    edges = geom.edges;

    % ---- colours ----
    colCirc = blendColor_(opts.CircularColor,    opts.CircularAlpha,    [1 1 1]);
    colDir  = blendColor_(opts.DirectionalColor, opts.DirectionalAlpha, [1 1 1]);
    [colUp, colDown] = plainEdgeColors_(opts);

    switch lower(strtrim(opts.OverlayComponent))
        case 'c'
            overlayLWField   = 'lwC';
            colBaseComponent = colDir;
            colOverlay       = colCirc;

        case 'd'
            overlayLWField   = 'lwD';
            colBaseComponent = colCirc;
            colOverlay       = colDir;

        otherwise
            error('plotTFL:BadOverlayComponent', ...
                'OverlayComponent must be ''c'' or ''d''.');
    end

    % ---- arrow lengths ----
    arrowLenW = computeArrowLengths_(scene, geom, vertcat(edges.lwW), opts);
    arrowLenOverlay = computeArrowLengths_(scene, geom, ...
        vertcat(edges.(overlayLWField)), opts);

    % ---- build items ----
    for e = 1:numel(edges)
        edge = edges(e);

        colBaseEdge = resolveOverlayBaseColor_( ...
            edge, opts, colBaseComponent, colUp, colDown);

        % ----- base shaft item -----
        if edge.lwW > 0
            hasArrowBase = false;
            arrowColBase = colBaseEdge;
            arrowLenBase = arrowLenW(e);

            switch opts.ArrowLayer
                case 'fullflow'
                    hasArrowBase = true;
                    arrowColBase = resolveArrowColor_(edge, opts, colBaseEdge, colOverlay);

                case 'both'
                    hasArrowBase = true;
                    arrowColBase = colBaseEdge;
            end

            %{
            switch lower(strtrim(opts.ArrowLayer))
                case 'base'
                    hasArrowBase = true;
                    arrowColBase = resolveArrowColor_(edge, opts, colBaseEdge, colOverlay);

                case 'both'
                    hasArrowBase = true;
                    arrowColBase = colBaseEdge;
            end
            %}

            itemBase = makeRenderItem_( ...
                edge, ...
                edge.lwW, ...
                colBaseEdge, ...
                hasArrowBase, ...
                arrowLenBase, ...
                arrowColBase);

            renderItems = appendRenderItem_(renderItems, itemBase);
        end

        % ----- overlay shaft item -----
        lwOverlay = edge.(overlayLWField);
        if lwOverlay > 0
            hasArrowOverlay = false;
            arrowColOverlay = colOverlay;

            switch opts.ArrowLayer
                case 'both'
                    % In the two-arrow case, let the overlay arrow show the
                    % component share visually.
                    arrowLenOv = arrowLenOverlay(e);

                otherwise
                    % When only the overlay arrow is shown, keep its size
                    % comparable to plain / cdfdcolor by using total weight.
                    arrowLenOv = arrowLenW(e);
            end

            %{
            switch lower(strtrim(opts.ArrowLayer))
                case 'both'
                    % In the two-arrow case, let the overlay arrow show the
                    % component share visually.
                    arrowLenOv = arrowLenOverlay(e);

                otherwise
                    % When only the overlay arrow is shown, keep its size
                    % comparable to plain / cdfdcolor by using total weight.
                    arrowLenOv = arrowLenW(e);
            end
            %}

            switch opts.ArrowLayer
                case 'overlay'
                    hasArrowOverlay = true;
                    arrowColOverlay = resolveArrowColor_(edge, opts, colBaseEdge, colOverlay);

                case 'both'
                    hasArrowOverlay = true;
                    arrowColOverlay = colOverlay;
            end

            %{
            switch lower(strtrim(opts.ArrowLayer))
                case {'overlay','circular'}
                    hasArrowOverlay = true;
                    arrowColOverlay = resolveArrowColor_(edge, opts, colBaseEdge, colOverlay);

                case 'both'
                    hasArrowOverlay = true;
                    arrowColOverlay = colOverlay;
            end
            %}

            itemOverlay = makeRenderItem_( ...
                edge, ...
                lwOverlay, ...
                colOverlay, ...
                hasArrowOverlay, ...
                arrowLenOv, ...
                arrowColOverlay);

            renderItems = appendRenderItem_(renderItems, itemOverlay);
        end
    end

    % largest first, stable within ties
    renderItems = sortRenderItems_(renderItems);
end

function edgeLW = computePlainEdgeLineWidths_(geom, scene, opts)

    E = numel(geom.edges);
    if E == 0
        edgeLW = zeros(0,1);
        return;
    end

    % Raw per-edge total-flow weights
    edgeW = abs(double(geom.w(:)));
    N = numel(scene.Xplot);

    % Base linewidth
    edgeLW = opts.EdgeWidth * ones(E,1);

    % Weight-based displayed widths
    if strcmpi(opts.EdgeWidthMode, 'weight')
        edgeLW = mapWeightsToLineWidth_( ...
            edgeW, ...
            opts.EdgeWidthRange, ...
            opts.EdgeWidthScale, ...
            opts.EdgeWidthQuantiles, ...
            opts.EdgeWidth);
    end

    % Optional size-aware global rescaling
    nodeAutoOn = strcmpi(opts.NodeSizeMode, 'auto');
    edgeScale = 1;
    if strcmpi(opts.EdgeWidthScaleMode, 'auto') || ...
       (strcmpi(opts.EdgeWidthScaleMode, 'withnode') && nodeAutoOn)

        Nref = max(1, opts.EdgeWidthRefN);
        edgeScale = (Nref / N) ^ opts.EdgeWidthGamma;
        edgeScale = min(max(edgeScale, opts.EdgeWidthMinScale), ...
                        opts.EdgeWidthMaxScale);
    end

    edgeLW = edgeLW * edgeScale;
    edgeLW = min(max(edgeLW, opts.EdgeWidthMin), opts.EdgeWidthMax);
end

function [colUp, colDown] = plainEdgeColors_(opts)
    bg = [1 1 1];
    colUp   = blendColor_(opts.UpEdgeColor,   opts.EdgeAlpha, bg);
    colDown = blendColor_(opts.DownEdgeColor, opts.EdgeAlpha, bg);
end

function col = resolveOverlayBaseColor_(edge, opts, colBaseComponent, colUp, colDown)

    switch lower(strtrim(opts.OverlayBaseColorMode))
        case {'component','oppositecomponent','legacy'}
            col = colBaseComponent;

        case {'plain','fullflow'}
            if edge.isUp
                col = colUp;
            else
                col = colDown;
            end

        otherwise
            error('plotTFL:BadOverlayBaseColorMode', ...
                'OverlayBaseColorMode must be ''component'' or ''plain''.');
    end
end

%{
function [geom, scene, renderItems] = ...
    buildCDFDColorRenderItems_(W, scene, geom, opts)

    renderItems = struct([]);

    if ~opts.ShowEdges || isempty(geom.edges)
        return;
    end

    [C, D, info] = resolveCDFD_(W, opts);
    scene.C = C;
    scene.D = D;
    scene.CDFDInfo = info;

    geom = attachEdgeVisualData_(geom, W, C, D);

    plainLW = computePlainEdgeLineWidths_(geom, scene, opts);
    for e = 1:numel(geom.edges)
        geom.edges(e).lwW = plainLW(e);
    end

    edges = geom.edges;
    arrowLenW = computeArrowLengths_(scene, geom, vertcat(edges.lwW), opts);

    colCirc = blendColor_(opts.CircularColor,    opts.CircularAlpha,    [1 1 1]);
    colDir  = blendColor_(opts.DirectionalColor, opts.DirectionalAlpha, [1 1 1]);

    for e = 1:numel(edges)
        edge = edges(e);
        if edge.lwW <= 0
            continue;
        end

        col = edge.rhoC * colCirc + edge.rhoD * colDir;

        item = makeRenderItem_( ...
            edge, ...
            edge.lwW, ...
            col, ...
            true, ...
            arrowLenW(e), ...
            col);

        renderItems = appendRenderItem_(renderItems, item);
    end

    renderItems = sortRenderItems_(renderItems);
end

%}

function arrowLenE = computeArrowLengths_(scene, geom, lwVec, opts)

    E = numel(geom.edges);
    if E == 0
        arrowLenE = zeros(0,1);
        return;
    end

    lwVec = double(lwVec(:));
    if numel(lwVec) ~= E
        error('plotTFL:BadArrowLWVec', ...
            'Arrow linewidth vector must have one entry per edge.');
    end

    Xplot = scene.Xplot;
    Yplot = scene.Yplot;
    N = numel(Xplot);

    tolLen = 1e-12;
    spanMax = max(rangeNonzero_(Xplot), rangeNonzero_(Yplot));

    edgeLens = double(geom.len(:));
    Ltyp = localRobustMedianPositive_(edgeLens, tolLen);
    if ~isfinite(Ltyp) || Ltyp < tolLen
        Ltyp = spanMax;
    end
    if ~isfinite(spanMax) || spanMax < tolLen
        spanMax = 1;
    end

    nodeAutoOn = strcmpi(opts.NodeSizeMode, 'auto');

    arrowScale = 1;
    if strcmpi(opts.ArrowSizeMode, 'auto') || ...
       (strcmpi(opts.ArrowSizeMode, 'withnode') && nodeAutoOn)

        Nref = max(1, opts.ArrowSizeRefN);
        arrowScale = (Nref / N) ^ opts.ArrowSizeGamma;
        arrowScale = min(max(arrowScale, opts.ArrowSizeMinScale), ...
                         opts.ArrowSizeMaxScale);
    end

    switch lower(strtrim(opts.ArrowLenBasis))
        case 'ltyp'
            Lbase = Ltyp;
        otherwise
            Lbase = spanMax;
    end

    arrowLen = (opts.ArrowSize * arrowScale) * Lbase;

    if opts.ArrowMinFrac > 0
        arrowLen = max(arrowLen, opts.ArrowMinFrac * Lbase);
    end
    if opts.ArrowMaxFrac > 0
        arrowLen = min(arrowLen, opts.ArrowMaxFrac * Lbase);
    end

    refLW = median(lwVec(isfinite(lwVec) & lwVec > 0));
    if isempty(refLW) || ~isfinite(refLW) || refLW <= 0
        refLW = 1;
    end

    arrowLenE = arrowLen * ones(E,1);
    if opts.ArrowScaleWithEdgeWidth
        s = (lwVec / refLW) .^ opts.ArrowEdgeGamma;
        s = min(max(s, opts.ArrowEdgeMinScale), opts.ArrowEdgeMaxScale);
        arrowLenE = arrowLen .* s;
    end
end

%function scene = expandSceneLimitsToGeometry_(scene, geom, opts)
function scene = expandSceneLimitsToGeometry_(scene, geom)
%EXPANDSCENELIMITSTOGEOMETRY_  Ensure suggested limits include realised geometry.
%
% Limits are expanded to include:
%   - plotted node coordinates
%   - all realised edge polylines in geom.edges
%
% Horizontal reference bands do NOT affect limits; they are annotations.

    % ---- start from plotted node coordinates ----
    xMin = min(scene.Xplot);
    xMax = max(scene.Xplot);
    yMin = min(scene.Yplot);
    yMax = max(scene.Yplot);

    % ---- expand to include realised edge geometry ----
    for e = 1:numel(geom.edges)
        xe = geom.edges(e).x(:);
        ye = geom.edges(e).y(:);

        xe = xe(isfinite(xe));
        ye = ye(isfinite(ye));

        if ~isempty(xe)
            xMin = min(xMin, min(xe));
            xMax = max(xMax, max(xe));
        end
        if ~isempty(ye)
            yMin = min(yMin, min(ye));
            yMax = max(yMax, max(ye));
        end
    end

    % ---- reuse scene padding ----
    scene.xlimSuggested = [xMin - scene.padX, xMax + scene.padX];
    scene.ylimSuggested = [yMin - scene.padY, yMax + scene.padY];
end

function item = makeRenderItem_(edge, lw, shaftCol, hasArrow, arrowLen, arrowCol)
%MAKERENDERITEM_  Bundle one drawable shaft+arrow object.

    item = struct();

    item.x = edge.x;
    item.y = edge.y;
    item.lw = lw;
    item.shaftColor = shaftCol;

    item.hasArrow = logical(hasArrow) && ...
        isfinite(edge.arrowX1) && isfinite(edge.arrowY1) && ...
        isfinite(edge.arrowX2) && isfinite(edge.arrowY2);

    item.arrowX1 = edge.arrowX1;
    item.arrowY1 = edge.arrowY1;
    item.arrowX2 = edge.arrowX2;
    item.arrowY2 = edge.arrowY2;
    item.toNode  = edge.toNode;

    item.arrowLen   = arrowLen;
    item.arrowColor = arrowCol;
end


function items = appendRenderItem_(items, item)
%APPENDRENDERITEM_  Append one render item to an array of items.

    if isempty(items)
        items = item;
    else
        items(end+1,1) = item;
    end
end


function items = sortRenderItems_(items)
%SORTRENDERITEMS_  Draw largest first so smaller objects remain visible.

    if isempty(items)
        return;
    end

    lw = vertcat(items.lw);
    idx = (1:numel(items)).';

    % Largest linewidth first; stable insertion-order tie break
    [~, ord] = sortrows([-lw, idx], [1 2]);
    items = items(ord);
end


function drawRenderItems_(ax, items, nodeRadiusData, opts)
%DRAWRENDERITEMS_  Draw each shaft and its arrow together.

    if isempty(items)
        return;
    end

    useNodeTrim = ~isempty(nodeRadiusData);

    for k = 1:numel(items)
        it = items(k);

        % ---- shaft ----
        plot(ax, it.x, it.y, '-', ...
            'Color', it.shaftColor, ...
            'LineWidth', it.lw, ...
            'HitTest', 'off', ...
            'PickableParts', 'none');

        % ---- arrow ----
        if it.hasArrow
            nodeRadius = 0;
            if useNodeTrim
                nodeRadius = nodeRadiusData(it.toNode);
            end

            drawArrowheadToNode_(ax, ...
                it.arrowX1, it.arrowY1, ...
                it.arrowX2, it.arrowY2, ...
                it.arrowColor, it.arrowLen, opts.ArrowAngle, ...
                nodeRadius);
        end
    end
end