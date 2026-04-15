function out = ringFlowPlot(W, varargin)
%RINGFLOWPLOT  Plot a weighted directed network on a node ring with interior arcs.
%
%   out = RINGFLOWPLOT(W)
%   out = RINGFLOWPLOT(W, 'Name', Value, ...)
%
% Inputs
% ------
%   W : square weighted adjacency matrix, rows = source, cols = destination
%
% Name-value options
% ------------------
% Layout
%   'LayoutOrder'      : explicit node order around ring (default = [])
%   'SortBy'           : vector used to sort nodes around ring if LayoutOrder empty
%   'SortDirection'    : 'ascend' (default) or 'descend'
%   'StartAngle'       : starting angle in radians (default = pi/2)
%   'Clockwise'        : true/false (default = false)
%   'Radius'           : ring radius (default = 1)
%
% Labels
%   'Labels'           : cellstr/string labels (default = {'1','2',...})
%   'ShowLabels'       : true/false (default = true)
%   'LabelColor'       : label colour (default = [1 0 0])
%   'LabelFontSize'    : label font size (default = 9)
%   'LabelRadiusFactor': label radius relative to ring (default = 1.10)
%
% Nodes
%   'NodeSize'         : fixed marker area if NodeSizeData empty (default = 60)
%   'NodeSizeData'     : vector whose values are displayed by node size
%                        (default = [])
%   'NodeSizeScaleData': optional reference data used to determine the
%                        node-size mapping; if empty, NodeSizeData is used
%                        (default = [])
%   'NodeSizeRange'    : marker area range if NodeSizeData used (default = [40 180])
%   'NodeSizeScaleMode': 'area' (default), 'radius', or 'log'
%   'NodeColor'        : single node colour (default = [0.1 0.3 0.7])
%   'NodeColors'       : n-by-3 custom node colours (default = [])
%
% Ring
%   'ShowRing'         : true/false (default = true)
%   'RingColor'        : ring colour (default = [0.75 0.75 0.75])
%   'RingWidth'        : ring line width (default = 0.75)
%
% Edges
%   'EdgeColorMode'    : 'single' (default), 'source', 'target', or 'ringdirection'
%   'EdgeColor'        : single edge colour if EdgeColorMode='single'
%                        (default = [0 0 0])
%   'RingDirColorPos'  : colour for edges with positive wrapToPi(theta_j-theta_i)
%                        when EdgeColorMode='ringdirection' (default blue-ish)
%   'RingDirColorNeg'  : colour for edges with negative wrapToPi(theta_j-theta_i)
%                        when EdgeColorMode='ringdirection' (default orange-ish)
%   'RingDirZeroColor' : fallback colour if exactly tied (default grey)
%   'EdgeAlpha'        : blend with white (default = 0.70)
%   'EdgeWidthMode'    : 'fixed' (default) or 'weight'
%   'EdgeWidth'        : fixed width if EdgeWidthMode='fixed' (default = 0.75)
%   'EdgeWidthRange'   : width range if EdgeWidthMode='weight' (default = [0.6 2.5])
%   'EdgeWidthScale'   : 'log' (default) or 'linear'
%   'EdgeWidthQuantiles': quantile clip range (default = [0.05 0.95])
%
% Arc geometry
%   'Curvature'        : inward pull of arcs, in [0,1] (default = 0.70)
%   'CurveNpts'        : points per arc (default = 60)
%   'MutualOverlay'    : true/false (default = true). If true, mutual pairs
%                        use nested bilateral overlay on a shared centerline.
%   'MutualLaneOffset' : lateral separation for opposite directions between
%                        the same node pair when MutualOverlay=false,
%                        as fraction of ring radius (default = 0.06)
%
% Arrow geometry
%   'ArrowSize'        : arrow length as fraction of radius (default = 0.08)
%   'ArrowAngle'       : arrow half-angle in radians (default = pi/7)
%   'ArrowNodeGapPts'  : extra gap beyond node boundary, in points (default = 1)
%   'ArrowBacktrackFactor' : anchor arrowhead this many arrow-lengths back
%                        along the curve (default = 1.6)
%   'ArrowScaleWithEdgeWidth' : true/false (default = true)
%   'ArrowEdgeGamma'   : width-to-arrow scaling exponent (default = 0.7)
%   'ArrowEdgeMinScale': min arrow scale clamp (default = 0.75)
%   'ArrowEdgeMaxScale': max arrow scale clamp (default = 1.8)
%
% Lane scaling (used only when MutualOverlay=false)
%   'MutualLaneScaleWithEdgeWidth' : true/false (default = true)
%   'MutualLaneGamma'    : width-to-lane scaling exponent (default = 0.7)
%   'MutualLaneMinScale' : min lane scale clamp (default = 0.8)
%   'MutualLaneMaxScale' : max lane scale clamp (default = 3.0)
%
% Axes
%   'Parent'           : target axes (default = gca)
%   'ClearAxes'        : true/false (default = true)
%   'Title'            : title string (default = '')
%
% Output
% ------
%   out.X, out.Y, out.theta, out.order, out.ax
%
% Notes
% -----
%   - Self-loops are ignored.
%   - By default, mutual pairs use nested bilateral overlay:
%       larger direction first, smaller direction overlaid on the same path.
%   - Arrowheads are drawn after nodes, using an anchor point selected by
%     arc length rather than just the final sampled segment.
%
% -------------------------------------------------------------------------

    % ---- checks ----
    if ~ismatrix(W) || size(W,1) ~= size(W,2)
        error('ringFlowPlot:WNotSquare', 'W must be a square adjacency matrix.');
    end
    if any(~isfinite(W(:)))
        error('ringFlowPlot:NonFinite', 'W must contain only finite values.');
    end
    W = sparse(double(W));
    n = size(W,1);

    % ---- defaults ----
    opts = struct();

    % Layout
    opts.LayoutOrder       = [];
    opts.SortBy            = [];
    opts.SortDirection     = 'ascend';
    opts.StartAngle        = pi/2;
    opts.Clockwise         = false;
    opts.Radius            = 1;

    % Labels
    opts.Labels            = arrayfun(@(k) sprintf('%d', k), 1:n, 'UniformOutput', false);
    opts.ShowLabels        = true;
    opts.LabelColor        = [1 0 0];
    opts.LabelFontSize     = 9;
    opts.LabelRadiusFactor = 1.10;
    opts.LabelWrap          = true;
    opts.LabelWrapMaxChars  = 13;
    opts.FitAxesToLabels    = true;
    opts.LabelExtentPadFrac = 0.01;

    % Nodes
    opts.NodeSize          = 60;
    opts.NodeSizeData      = [];
    opts.NodeSizeScaleData = [];
    opts.NodeSizeRange     = [40 180];
    opts.NodeSizeScaleMode = 'area';
    opts.NodeColor         = [0.1 0.3 0.7];
    opts.NodeColors        = [];

    % Ring
    opts.ShowRing          = true;
    opts.RingColor         = [0.75 0.75 0.75];
    opts.RingWidth         = 0.75;

    % Edges / colours
    opts.EdgeColorMode      = 'single';
    opts.EdgeColor          = [0 0 0];
    opts.RingDirColorPos    = [0.35 0.35 0.35]; % changed to match tflPlot [0.10 0.55 0.85];
    opts.RingDirColorNeg    = [0.45 0.60 0.78]; %[0.85 0.45 0.15];
    opts.RingDirZeroColor   = [0.60 0.60 0.60]; %[0.40 0.40 0.40];
    opts.EdgeAlpha          = 0.70;
    opts.EdgeWidthMode      = 'fixed';
    opts.EdgeWidth          = 0.75;
    opts.EdgeWidthRange     = [0.6 2.5];
    opts.EdgeWidthScale     = 'log';
    opts.EdgeWidthQuantiles = [0.05 0.95];
    opts.MutualOverlayShoulder = 1.5;   % linewidth units

    % Arc geometry
    opts.Curvature        = 0.70;
    opts.CurveNpts        = 60;
    opts.MutualOverlay    = true;
    opts.MutualLaneOffset = 0.06;

    % Arrow geometry
    opts.ArrowSize              = 0.08;
    opts.ArrowAngle             = pi/7;
    opts.ArrowNodeGapPts        = 1;
    opts.ArrowBacktrackFactor   = 1.6;
    opts.ArrowScaleWithEdgeWidth = true;
    opts.ArrowEdgeGamma         = 0.7;
    opts.ArrowEdgeMinScale      = 0.75;
    opts.ArrowEdgeMaxScale      = 1.8;
    opts.ArrowPlacement    = 'late';   % 'late' or 'node'
    opts.ArrowPositionFrac = 0.82;     % used when ArrowPlacement='late'

    % Lane scaling (used only when MutualOverlay=false)
    opts.MutualLaneScaleWithEdgeWidth = true;
    opts.MutualLaneGamma    = 0.7;
    opts.MutualLaneMinScale = 0.8;
    opts.MutualLaneMaxScale = 3.0;

    % Axes
    opts.Parent            = [];
    opts.ClearAxes         = true;
    opts.Title             = '';

    % ---- parse ----
    if mod(numel(varargin),2) ~= 0
        error('ringFlowPlot:BadArgs', 'Optional arguments must be name/value pairs.');
    end

    for k = 1:2:numel(varargin)
        name = varargin{k};
        value = varargin{k+1};
        if ~ischar(name) && ~isstring(name)
            error('ringFlowPlot:BadParamName', 'Parameter names must be strings.');
        end
        name = lower(char(name));

        switch name
            case 'layoutorder'
                opts.LayoutOrder = value;
            case 'sortby'
                opts.SortBy = value;
            case 'sortdirection'
                opts.SortDirection = char(value);
            case 'startangle'
                opts.StartAngle = value;
            case 'clockwise'
                opts.Clockwise = logical(value);
            case 'radius'
                opts.Radius = value;

            case 'labels'
                opts.Labels = value;
            case 'showlabels'
                opts.ShowLabels = logical(value);
            case 'labelcolor'
                opts.LabelColor = value;
            case 'labelfontsize'
                opts.LabelFontSize = value;
            case 'labelradiusfactor'
                opts.LabelRadiusFactor = value;
            case 'labelwrap'
                opts.LabelWrap = logical(value);
            case 'labelwrapmaxchars'
                opts.LabelWrapMaxChars = value;
            case 'fitaxestolabels'
                opts.FitAxesToLabels = logical(value);
            case 'labelextentpadfrac'
                opts.LabelExtentPadFrac = value;

            case 'nodesize'
                opts.NodeSize = value;
            case 'nodesizedata'
                opts.NodeSizeData = value;
            case 'nodesizescaledata'
                opts.NodeSizeScaleData = value;
            case 'nodesizerange'
                opts.NodeSizeRange = value;
            case 'nodesizescalemode'
                opts.NodeSizeScaleMode = char(value);
            case 'nodecolor'
                opts.NodeColor = value;
            case 'nodecolors'
                opts.NodeColors = value;

            case 'showring'
                opts.ShowRing = logical(value);
            case 'ringcolor'
                opts.RingColor = value;
            case 'ringwidth'
                opts.RingWidth = value;

            case 'edgecolormode'
                opts.EdgeColorMode = char(value);
            case 'edgecolor'
                opts.EdgeColor = value;
            case 'ringdircolorpos'
                opts.RingDirColorPos = value;
            case 'ringdircolorneg'
                opts.RingDirColorNeg = value;
            case 'ringdirzerocolor'
                opts.RingDirZeroColor = value;
            case 'edgealpha'
                opts.EdgeAlpha = value;
            case 'edgewidthmode'
                opts.EdgeWidthMode = char(value);
            case 'edgewidth'
                opts.EdgeWidth = value;
            case 'edgewidthrange'
                opts.EdgeWidthRange = value;
            case 'edgewidthscale'
                opts.EdgeWidthScale = char(value);
            case 'edgewidthquantiles'
                opts.EdgeWidthQuantiles = value;
            case 'mutualoverlayshoulder'
                opts.MutualOverlayShoulder = value;

            case 'curvature'
                opts.Curvature = value;
            case 'curvenpts'
                opts.CurveNpts = value;
            case 'mutualoverlay'
                opts.MutualOverlay = logical(value);
            case 'mutuallaneoffset'
                opts.MutualLaneOffset = value;

            case 'arrowsize'
                opts.ArrowSize = value;
            case 'arrowangle'
                opts.ArrowAngle = value;
            case 'arrownodegappts'
                opts.ArrowNodeGapPts = value;
            case 'arrowbacktrackfactor'
                opts.ArrowBacktrackFactor = value;
            case 'arrowscalewithedgewidth'
                opts.ArrowScaleWithEdgeWidth = logical(value);
            case 'arrowedgegamma'
                opts.ArrowEdgeGamma = value;
            case 'arrowedgeminscale'
                opts.ArrowEdgeMinScale = value;
            case 'arrowedgemaxscale'
                opts.ArrowEdgeMaxScale = value;
            case 'arrowplacement'
                opts.ArrowPlacement = char(value);
            case 'arrowpositionfrac'
                opts.ArrowPositionFrac = value;

            case 'mutuallanescalewithedgewidth'
                opts.MutualLaneScaleWithEdgeWidth = logical(value);
            case 'mutuallanegamma'
                opts.MutualLaneGamma = value;
            case 'mutuallaneminscale'
                opts.MutualLaneMinScale = value;
            case 'mutuallanemaxscale'
                opts.MutualLaneMaxScale = value;

            case 'parent'
                opts.Parent = value;
            case 'clearaxes'
                opts.ClearAxes = logical(value);
            case 'title'
                opts.Title = value;

            otherwise
                error('ringFlowPlot:UnknownOption', 'Unknown option: %s', name);
        end
    end

    %{
    % ---- labels ----
    if isempty(opts.Labels)
        labels = arrayfun(@(k) sprintf('%d', k), 1:n, 'UniformOutput', false);
    else
        labels = opts.Labels;
        if isstring(labels)
            labels = cellstr(labels(:));
        end
        if numel(labels) ~= n
            error('ringFlowPlot:BadLabels', 'Labels must have length n.');
        end
    end
    %}


    % ---- order / ring layout ----
    if ~isempty(opts.LayoutOrder)
        order = opts.LayoutOrder(:);
        if numel(order) ~= n || numel(unique(order)) ~= n || any(sort(order) ~= (1:n)')
            error('ringFlowPlot:BadLayoutOrder', ...
                'LayoutOrder must be a permutation of 1:n.');
        end
    elseif ~isempty(opts.SortBy)
        s = opts.SortBy(:);
        if numel(s) ~= n
            error('ringFlowPlot:BadSortBy', 'SortBy must be length n.');
        end
        if strcmpi(opts.SortDirection, 'descend')
            [~, order] = sort(s, 'descend');
        else
            [~, order] = sort(s, 'ascend');
        end
    else
        order = (1:n).';
    end

    m = (0:n-1).';
    sgn = 1;
    if opts.Clockwise
        sgn = -1;
    end
    angOrdered = opts.StartAngle + sgn * 2*pi*m/n;

    theta = nan(n,1);
    theta(order) = angOrdered;

    R = opts.Radius;
    X = R * cos(theta);
    Y = R * sin(theta);

    
    % ---- labels ----

    if isempty(opts.Labels)
        labels = arrayfun(@(k) sprintf('%d', k), 1:n, 'UniformOutput', false);
    else
        labels = opts.Labels;
        if isstring(labels)
            labels = cellstr(labels(:));
        end
        if numel(labels) ~= n
            error('ringFlowPlot:BadLabels', 'Labels must have length n.');
        end
    end
    
    % ---- axes ----
    ax = opts.Parent;
    if isempty(ax) || ~isgraphics(ax, 'axes')
        ax = gca;
    end

    holdState = ishold(ax);
    if opts.ClearAxes
        cla(ax);
    end
    hold(ax, 'on');
    set(ax, 'Color', 'w');

    % ---- node sizes ----
    %nodeSize = computeNodeSizes_(opts.NodeSizeData, opts.NodeSize, ...
    %    opts.NodeSizeRange, opts.NodeSizeScaleMode, n);
    nodeSize = computeNodeSizes_( ...
        opts.NodeSizeData, ...
        opts.NodeSizeScaleData, ...
        opts.NodeSize, ...
        opts.NodeSizeRange, ...
        opts.NodeSizeScaleMode, ...
        n);

    % ---- node colours ----
    if ~isempty(opts.NodeColors)
        nodeColors = opts.NodeColors;
        if size(nodeColors,1) ~= n || size(nodeColors,2) ~= 3
            error('ringFlowPlot:BadNodeColors', ...
                'NodeColors must be n-by-3.');
        end
    else
        nodeColors = repmat(opts.NodeColor, n, 1);
    end

    % ---- limits before arrow radius conversion ----
    limR = 1.20 * opts.LabelRadiusFactor * R;
    xlim(ax, [-limR limR]);
    ylim(ax, [-limR limR]);
    daspect(ax, [1 1 1]);
    pbaspect(ax, [1 1 1]);
    set(ax, 'XTick', [], 'YTick', [], 'XTickLabel', [], 'YTickLabel', []);
    box(ax, 'on');

    % ---- ring ----
    if opts.ShowRing
        tt = linspace(0, 2*pi, 400);
        plot(ax, R*cos(tt), R*sin(tt), '-', ...
            'Color', opts.RingColor, ...
            'LineWidth', opts.RingWidth);
    end

    % ---- directed edge data ----
    [iDir, jDir, wDir] = find(W);
    maskNoLoops = (iDir ~= jDir);
    iDir = iDir(maskNoLoops);
    jDir = jDir(maskNoLoops);
    wDir = wDir(maskNoLoops);

    E = numel(wDir);
    if E == 0
        % plot nodes/labels only
        scatter(ax, X, Y, nodeSize, ...
            'MarkerFaceColor', 'flat', ...
            'CData', nodeColors, ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 0.6);

        
        labelHs = gobjects(0,1);
        if opts.ShowLabels
            labelHs = drawLabels_(ax, labels, theta, R, opts);

            if opts.FitAxesToLabels && ~isempty(labelHs)
                expandAxesToLabelExtents_(ax, labelHs, opts.LabelExtentPadFrac);
                daspect(ax, [1 1 1]);
                pbaspect(ax, [1 1 1]);
            end
        end

        %{
        if opts.ShowLabels
            drawLabels_(ax, labels, theta, R, opts);
        end
        %}

        if ~isempty(opts.Title)
            title(ax, opts.Title, 'Interpreter', 'none');
        end

        if ~holdState
            hold(ax, 'off');
        end

        out = struct('X',X,'Y',Y,'theta',theta,'order',order,'ax',ax);
        return;
    end

    edgeLW = opts.EdgeWidth * ones(E,1);
    if strcmpi(opts.EdgeWidthMode, 'weight')
        edgeLW = mapWeightsToLineWidth_(wDir, ...
            opts.EdgeWidthRange, opts.EdgeWidthScale, opts.EdgeWidthQuantiles, opts.EdgeWidth);
    end

    refLW = median(edgeLW(edgeLW > 0));
    if isempty(refLW) || ~isfinite(refLW) || refLW <= 0
        refLW = 1;
    end

    baseLaneOffset = opts.MutualLaneOffset * R;
    baseArrowLen   = opts.ArrowSize * R;

    % Per-directed-edge arrow scaling
    arrowLenDir = zeros(E,1);
    for e = 1:E
        if opts.ArrowScaleWithEdgeWidth
            s = (edgeLW(e) / refLW) ^ opts.ArrowEdgeGamma;
            s = min(max(s, opts.ArrowEdgeMinScale), opts.ArrowEdgeMaxScale);
        else
            s = 1;
        end
        arrowLenDir(e) = baseArrowLen * s;
    end

    % Map (i,j) -> directed-edge index
    eIdx = sparse(iDir, jDir, 1:E, n, n);

    % ---- build drawing jobs (unordered pairs) ----
    jobs = struct('kind', {}, 'i', {}, 'j', {}, 'e1', {}, 'e2', {}, 'repLW', {});
    for i = 1:n
        for j = i+1:n
            eij = full(eIdx(i,j));
            eji = full(eIdx(j,i));

            if eij > 0 && eji > 0
                jobs(end+1) = struct( ... %#ok<AGROW>
                    'kind', 'mutual', ...
                    'i', i, 'j', j, ...
                    'e1', eij, 'e2', eji, ...
                    'repLW', max(edgeLW([eij, eji])) );
            elseif eij > 0
                jobs(end+1) = struct( ... %#ok<AGROW>
                    'kind', 'single', ...
                    'i', i, 'j', j, ...
                    'e1', eij, 'e2', 0, ...
                    'repLW', edgeLW(eij) );
            elseif eji > 0
                jobs(end+1) = struct( ... %#ok<AGROW>
                    'kind', 'single', ...
                    'i', j, 'j', i, ...
                    'e1', eji, 'e2', 0, ...
                    'repLW', edgeLW(eji) );
            end
        end
    end

    if ~isempty(jobs)
        [~, ordJobs] = sort([jobs.repLW], 'ascend');
        jobs = jobs(ordJobs);
    end

    % ---- arrow storage: full curves, per directed edge ----
    curveX   = cell(0,1);
    curveY   = cell(0,1);
    arrowTo  = zeros(0,1);
    arrowCol = zeros(0,3);
    arrowLen = zeros(0,1);

   % ---- draw jobs ----
    for q = 1:numel(jobs)
        jb = jobs(q);

        if strcmp(jb.kind, 'single')
            i = jb.i;
            j = jb.j;
            e = jb.e1;

            col = edgeColorForEdge_(i, j, opts, nodeColors, theta);
            col = blendColor_(col, opts.EdgeAlpha, [1 1 1]);

            [xCurve, yCurve] = ringArcCurve_( ...
                X(i), Y(i), X(j), Y(j), ...
                theta(i), theta(j), R, opts.Curvature, opts.CurveNpts, 0);

            plot(ax, xCurve, yCurve, '-', 'Color', col, 'LineWidth', edgeLW(e));

            curveX{end+1,1} = xCurve; %#ok<AGROW>
            curveY{end+1,1} = yCurve; %#ok<AGROW>
            arrowTo(end+1,1)  = j; %#ok<AGROW>
            arrowCol(end+1,:) = col; %#ok<AGROW>
            arrowLen(end+1,1) = arrowLenDir(e); %#ok<AGROW>

        else
            i   = jb.i;
            j   = jb.j;
            eij = jb.e1;
            eji = jb.e2;

            colIJ = edgeColorForEdge_(i, j, opts, nodeColors, theta);
            colJI = edgeColorForEdge_(j, i, opts, nodeColors, theta);

            colIJ = blendColor_(colIJ, opts.EdgeAlpha, [1 1 1]);
            colJI = blendColor_(colJI, opts.EdgeAlpha, [1 1 1]);

            if opts.MutualOverlay
                % Shared centerline, larger direction first, smaller nested on top
                [xCurve, yCurve] = ringArcCurve_( ...
                    X(i), Y(i), X(j), Y(j), ...
                    theta(i), theta(j), R, opts.Curvature, opts.CurveNpts, 0);

                if edgeLW(eij) >= edgeLW(eji)
                    lwBase = edgeLW(eij);
                    lwTop  = min(edgeLW(eji), max(0, lwBase - opts.MutualOverlayShoulder));

                    plot(ax, xCurve, yCurve, '-', 'Color', colIJ, 'LineWidth', lwBase);
                    if lwTop > 0
                        plot(ax, xCurve, yCurve, '-', 'Color', colJI, 'LineWidth', lwTop);
                    end
                else
                    lwBase = edgeLW(eji);
                    lwTop  = min(edgeLW(eij), max(0, lwBase - opts.MutualOverlayShoulder));

                    plot(ax, xCurve, yCurve, '-', 'Color', colJI, 'LineWidth', lwBase);
                    if lwTop > 0
                        plot(ax, xCurve, yCurve, '-', 'Color', colIJ, 'LineWidth', lwTop);
                    end
                end

                % Arrow for i -> j
                curveX{end+1,1} = xCurve; %#ok<AGROW>
                curveY{end+1,1} = yCurve; %#ok<AGROW>
                arrowTo(end+1,1)  = j; %#ok<AGROW>
                arrowCol(end+1,:) = colIJ; %#ok<AGROW>
                arrowLen(end+1,1) = arrowLenDir(eij); %#ok<AGROW>

                % Arrow for j -> i (reverse same centerline)
                curveX{end+1,1} = fliplr(xCurve); %#ok<AGROW>
                curveY{end+1,1} = fliplr(yCurve); %#ok<AGROW>
                arrowTo(end+1,1)  = i; %#ok<AGROW>
                arrowCol(end+1,:) = colJI; %#ok<AGROW>
                arrowLen(end+1,1) = arrowLenDir(eji); %#ok<AGROW>

            else
                % Fallback: separate tracks for the two directions
                if opts.MutualLaneScaleWithEdgeWidth
                    s1 = (edgeLW(eij) / refLW) ^ opts.MutualLaneGamma;
                    s1 = min(max(s1, opts.MutualLaneMinScale), opts.MutualLaneMaxScale);

                    s2 = (edgeLW(eji) / refLW) ^ opts.MutualLaneGamma;
                    s2 = min(max(s2, opts.MutualLaneMinScale), opts.MutualLaneMaxScale);
                else
                    s1 = 1;
                    s2 = 1;
                end

                lane1 = mutualLaneOffset_(i, j, W, theta, baseLaneOffset * s1);
                lane2 = mutualLaneOffset_(j, i, W, theta, baseLaneOffset * s2);

                [x1, y1] = ringArcCurve_( ...
                    X(i), Y(i), X(j), Y(j), ...
                    theta(i), theta(j), R, opts.Curvature, opts.CurveNpts, lane1);

                [x2, y2] = ringArcCurve_( ...
                    X(j), Y(j), X(i), Y(i), ...
                    theta(j), theta(i), R, opts.Curvature, opts.CurveNpts, lane2);

                plot(ax, x1, y1, '-', 'Color', colIJ, 'LineWidth', edgeLW(eij));
                plot(ax, x2, y2, '-', 'Color', colJI, 'LineWidth', edgeLW(eji));

                curveX{end+1,1} = x1; %#ok<AGROW>
                curveY{end+1,1} = y1; %#ok<AGROW>
                arrowTo(end+1,1)  = j; %#ok<AGROW>
                arrowCol(end+1,:) = colIJ; %#ok<AGROW>
                arrowLen(end+1,1) = arrowLenDir(eij); %#ok<AGROW>

                curveX{end+1,1} = x2; %#ok<AGROW>
                curveY{end+1,1} = y2; %#ok<AGROW>
                arrowTo(end+1,1)  = i; %#ok<AGROW>
                arrowCol(end+1,:) = colJI; %#ok<AGROW>
                arrowLen(end+1,1) = arrowLenDir(eji); %#ok<AGROW>
            end
        end
    end

    % ---- nodes ----
    scatter(ax, X, Y, nodeSize, ...
        'MarkerFaceColor', 'flat', ...
        'CData', nodeColors, ...
        'MarkerEdgeColor', 'k', ...
        'LineWidth', 0.6);

    
    % ---- arrowheads ----
    drawnow;
    nodeRadiusData = markerRadiusData_(ax, nodeSize, opts.ArrowNodeGapPts);

    for k = 1:numel(arrowTo)
        j = arrowTo(k);

        if isempty(curveX{k}) || numel(curveX{k}) < 2
            continue;
        end

        switch lower(strtrim(opts.ArrowPlacement))
            case 'late'
                backDist = opts.ArrowBacktrackFactor * arrowLen(k);

                [x1, y1, x2, y2] = arrowSegmentAtFraction_( ...
                    curveX{k}, curveY{k}, opts.ArrowPositionFrac, backDist);

                % draw with nodeRadius = 0 because this is not node-terminal
                drawArrowheadToNode_(ax, ...
                    x1, y1, x2, y2, ...
                    arrowCol(k,:), arrowLen(k), opts.ArrowAngle, 0);

            otherwise   % 'node'
                backDist = max(opts.ArrowBacktrackFactor * arrowLen(k), ...
                           nodeRadiusData(j) + 0.25 * arrowLen(k));

                [x1, y1, x2, y2] = arrowAnchorFromCurve_( ...
                    curveX{k}, curveY{k}, backDist);

                drawArrowheadToNode_(ax, ...
                    x1, y1, x2, y2, ...
                    arrowCol(k,:), arrowLen(k), opts.ArrowAngle, ...
                    nodeRadiusData(j));
        end
    end
    
    %{
    
    % ---- arrowheads on top of nodes ----
    drawnow;
    nodeRadiusData = markerRadiusData_(ax, nodeSize, opts.ArrowNodeGapPts);

    for k = 1:numel(arrowTo)
        j = arrowTo(k);

        if isempty(curveX{k}) || numel(curveX{k}) < 2
            continue;
        end

        backDist = max(opts.ArrowBacktrackFactor * arrowLen(k), ...
                       nodeRadiusData(j) + 0.25 * arrowLen(k));

        [x1, y1, x2, y2] = arrowAnchorFromCurve_(curveX{k}, curveY{k}, backDist);

        drawArrowheadToNode_(ax, ...
            x1, y1, x2, y2, ...
            arrowCol(k,:), arrowLen(k), opts.ArrowAngle, ...
            nodeRadiusData(j));
    end

    %}  

    %{
    % ---- labels ----
    if opts.ShowLabels
        drawLabels_(ax, labels, theta, R, opts);
    end
    %}

    % ---- labels ----
    labelHs = gobjects(0,1);
    if opts.ShowLabels
        labelHs = drawLabels_(ax, labels, theta, R, opts);

        if opts.FitAxesToLabels && ~isempty(labelHs)
            expandAxesToLabelExtents_(ax, labelHs, opts.LabelExtentPadFrac);
            daspect(ax, [1 1 1]);
            pbaspect(ax, [1 1 1]);
        end
    end

    % ---- title ----
    if ~isempty(opts.Title)
        title(ax, opts.Title, 'Interpreter', 'none');
    end

    if ~holdState
        hold(ax, 'off');
    end

    out = struct();
    out.X = X;
    out.Y = Y;
    out.theta = theta;
    out.order = order;
    out.ax = ax;

    % Temporary outputs for debugging / inspection
    out.nodeSize = nodeSize;
    out.opts = opts;

    
end

% =====================================================================
function nodeSize = computeNodeSizes_(sData, sScale, baseSize, sizeRange, scaleMode, n)

    if isempty(sData)
        nodeSize = baseSize * ones(n,1);
        return;
    end

    sData = double(sData(:));
    if numel(sData) ~= n
        error('ringFlowPlot:BadNodeSizeData', ...
            'NodeSizeData must be empty or length n.');
    end
    if any(~isfinite(sData))
        error('ringFlowPlot:BadNodeSizeData', ...
            'NodeSizeData must contain only finite values.');
    end

    if isempty(sScale)
        sScale = sData;
    else
        sScale = double(sScale(:));
        if any(~isfinite(sScale))
            error('ringFlowPlot:BadNodeSizeScaleData', ...
                'NodeSizeScaleData must contain only finite values.');
        end
    end

    sData  = max(sData,  0);
    sScale = max(sScale, 0);

    if all(sData == 0) || isempty(sScale) || all(sScale == 0)
        nodeSize = baseSize * ones(n,1);
        return;
    end

    switch lower(strtrim(scaleMode))
        case 'log'
            zData  = log1p(sData);
            zScale = log1p(sScale);

        case 'radius'
            zData  = sqrt(sData);
            zScale = sqrt(sScale);

        otherwise   % 'area'
            % Keep current visual convention
            zData  = sqrt(sData);
            zScale = sqrt(sScale);
    end

    zLo = min(zScale);
    zHi = max(zScale);

    if ~isfinite(zLo) || ~isfinite(zHi) || zHi <= zLo
        nodeSize = mean(sizeRange) * ones(n,1);
        return;
    end

    z = (zData - zLo) ./ (zHi - zLo);
    z = min(max(z, 0), 1);

    nodeSize = sizeRange(1) + z * (sizeRange(2) - sizeRange(1));
end

%{
function nodeSize = computeNodeSizes_(s, baseSize, sizeRange, scaleMode, n)
    if isempty(s)
        nodeSize = baseSize * ones(n,1);
        return;
    end

    s = double(s(:));
    if numel(s) ~= n
        error('ringFlowPlot:BadNodeSizeData', ...
            'NodeSizeData must be empty or length n.');
    end

    s = max(s, 0);

    if all(s == 0) || all(~isfinite(s))
        nodeSize = baseSize * ones(n,1);
        return;
    end

    switch lower(strtrim(scaleMode))
        case 'log'
            z = log1p(s);
        case 'radius'
            z = sqrt(s);
        otherwise   % 'area'
            z = sqrt(s);
    end

    z = z - min(z);
    if max(z) > 0
        z = z / max(z);
    end

    nodeSize = sizeRange(1) + z * (sizeRange(2) - sizeRange(1));
end
%}

% =====================================================================
function labelHs = drawLabels_(ax, labels, theta, R, opts)
    rLab = opts.LabelRadiusFactor * R;
    n = numel(labels);

    labelHs = gobjects(0,1);

    for i = 1:n
        xL = rLab * cos(theta(i));
        yL = rLab * sin(theta(i));

        ha = 'center';
        if xL > 0.10 * R
            ha = 'left';
        elseif xL < -0.10 * R
            ha = 'right';
        end

        va = 'middle';
        if yL > 0.25 * R
            va = 'bottom';
        elseif yL < -0.25 * R
            va = 'top';
        end

        sLab = labels{i};
        if opts.LabelWrap
            sLab = wrapRingLabel_(sLab, opts.LabelWrapMaxChars);
        end

        labelHs(end+1,1) = text(ax, xL, yL, sLab, ...
            'Color', opts.LabelColor, ...
            'FontSize', opts.LabelFontSize, ...
            'HorizontalAlignment', ha, ...
            'VerticalAlignment', va); %#ok<AGROW>
    end
end

%{
function drawLabels_(ax, labels, theta, R, opts)
    rLab = opts.LabelRadiusFactor * R;
    n = numel(labels);

    for i = 1:n
        xL = rLab * cos(theta(i));
        yL = rLab * sin(theta(i));

        ha = 'center';
        if xL > 0.10 * R
            ha = 'left';
        elseif xL < -0.10 * R
            ha = 'right';
        end

        va = 'middle';
        if yL > 0.25 * R
            va = 'bottom';
        elseif yL < -0.25 * R
            va = 'top';
        end

        text(ax, xL, yL, labels{i}, ...
            'Color', opts.LabelColor, ...
            'FontSize', opts.LabelFontSize, ...
            'HorizontalAlignment', ha, ...
            'VerticalAlignment', va);
    end
end
%}

% =====================================================================
%{
function col = edgeColorForEdge_(i, j, opts, nodeColors, theta)
    switch lower(strtrim(opts.EdgeColorMode))
        case 'source'
            col = nodeColors(i,:);

        case 'target'
            col = nodeColors(j,:);

        case 'ringdirection'
            dth = wrapToPi_(theta(j) - theta(i));
            if dth > 0
                col = opts.RingDirColorPos;
            elseif dth < 0
                col = opts.RingDirColorNeg;
            else
                col = opts.RingDirZeroColor;
            end

        otherwise
            col = opts.EdgeColor;
    end
end
%}


function col = edgeColorForEdge_(i, j, opts, nodeColors, theta)
    switch lower(strtrim(opts.EdgeColorMode))
        case 'source'
            col = nodeColors(i,:);

        case 'target'
            col = nodeColors(j,:);

        case 'ringdirection'
            s = ringDirectionSign_(i, j, theta);
            if s > 0
                col = opts.RingDirColorPos;
            elseif s < 0
                col = opts.RingDirColorNeg;
            else
                col = opts.RingDirZeroColor;
            end

        otherwise
            col = opts.EdgeColor;
    end
end

% =====================================================================
function laneOffset = mutualLaneOffset_(i, j, W, theta, baseOffset)
    if i == j || W(j,i) == 0
        laneOffset = 0;
        return;
    end

    s = ringDirectionSign_(i, j, theta);
    if s == 0
        laneOffset = 0;
    else
        laneOffset = s * baseOffset;
    end
end

%{

function laneOffset = mutualLaneOffset_(i, j, W, theta, baseOffset)
% Small signed offset to separate i->j and j->i onto different tracks.

    if i == j || W(j,i) == 0
        laneOffset = 0;
        return;
    end

    dth = wrapToPi_(theta(j) - theta(i));
    if dth == 0
        if i < j
            s = 1;
        else
            s = -1;
        end
    else
        s = sign(dth);
    end

    laneOffset = s * baseOffset;
end

%}

% =====================================================================
function [xCurve, yCurve] = ringArcCurve_(x1, y1, x2, y2, th1, th2, R, curvature, npts, laneOffset)
% Cubic Bezier arc through the interior of the ring, with optional lateral
% lane offset.

    p0 = [x1; y1];
    p3 = [x2; y2];

    sep = abs(wrapToPi_(th2 - th1));

    u0 = p0 / max(norm(p0), eps);
    u3 = p3 / max(norm(p3), eps);
    um = u0 + u3;

    if norm(um) < 1e-12
        chord = p3 - p0;
        if norm(chord) < 1e-12
            um = [0; 0];
        else
            um = [-chord(2); chord(1)];
            um = um / norm(um);
        end
        rCtrl = 0;
    else
        um = um / norm(um);
        rCtrl = R * (1 - curvature * sin(sep/2));
        rCtrl = max(0, min(R, rCtrl));
    end

    pMid = rCtrl * um;

    if norm(um) > 1e-12
        nperp = [-um(2); um(1)];
    else
        chord = p3 - p0;
        if norm(chord) < 1e-12
            nperp = [0; 0];
        else
            chord = chord / norm(chord);
            nperp = [-chord(2); chord(1)];
        end
    end

    alpha = 0.85;
    p1 = (1-alpha)*p0 + alpha*pMid + laneOffset*nperp;
    p2 = (1-alpha)*p3 + alpha*pMid + laneOffset*nperp;

    [xCurve, yCurve] = cubicBezier_(p0, p1, p2, p3, npts);
end

% =====================================================================
function [xCurve, yCurve] = cubicBezier_(p0, p1, p2, p3, npts)
    t = linspace(0,1,npts).';
    b0 = (1-t).^3;
    b1 = 3*(1-t).^2.*t;
    b2 = 3*(1-t).*t.^2;
    b3 = t.^3;
    B  = b0.*p0.' + b1.*p1.' + b2.*p2.' + b3.*p3.';
    xCurve = B(:,1).';
    yCurve = B(:,2).';
end

% =====================================================================
function [x1, y1, x2, y2] = arrowAnchorFromCurve_(xCurve, yCurve, backDist)
% Return a point (x1,y1) approximately backDist from the end of the curve,
% together with the end point (x2,y2).

    x2 = xCurve(end);
    y2 = yCurve(end);

    ds = hypot(diff(xCurve), diff(yCurve));
    if isempty(ds) || sum(ds) <= 0
        x1 = NaN; y1 = NaN;
        return;
    end

    acc = 0;
    for k = numel(ds):-1:1
        if acc + ds(k) >= backDist
            t = (backDist - acc) / max(ds(k), eps);
            x1 = xCurve(k+1) - t * (xCurve(k+1) - xCurve(k));
            y1 = yCurve(k+1) - t * (yCurve(k+1) - yCurve(k));
            return;
        end
        acc = acc + ds(k);
    end

    x1 = xCurve(1);
    y1 = yCurve(1);
end

% =====================================================================
function rData = markerRadiusData_(ax, markerAreaPts2, gapPts)
% Approximate scatter marker radius in data units.

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

    rData = 0.5 * (rx + ry);
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

    edgeCol = max(0, 0.65 * col);   % darker shade of same fill

    patch(ax, [tx bx1 bx2], [ty by1 by2], col, ...
        'EdgeColor', edgeCol, ...
        'LineWidth', 0.25, ...
        'FaceColor', col, ...
        'HitTest', 'off', ...
        'PickableParts', 'none');
%{
    patch(ax, [tx bx1 bx2], [ty by1 by2], col, ...
        'EdgeColor', 'none', ...
        'FaceColor', col, ...
        'HitTest', 'off', ...
        'PickableParts', 'none');
        %}
end

% =====================================================================
function cBlend = blendColor_(c, alpha, bg)
    c  = c(:)';
    bg = bg(:)';
    cBlend = alpha * c + (1-alpha) * bg;
end

% =====================================================================
function lw = mapWeightsToLineWidth_(w, lwRange, scaleMode, qRange, fallbackLW)
    w = double(w(:));
    lw = fallbackLW * ones(size(w));

    if isempty(w) || all(~isfinite(w)) || all(w <= 0)
        return;
    end

    wpos = w(isfinite(w) & w > 0);
    if isempty(wpos)
        return;
    end

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
            x = log10(wClip);
            xLo = log10(wLo);
            xHi = log10(wHi);
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
    t = min(max(t, 0), 1);

    lw = lwRange(1) + t .* (lwRange(2) - lwRange(1));
end

% =====================================================================
function a = wrapToPi_(a)
    a = mod(a + pi, 2*pi) - pi;
end

% =====================================================================
function [x1, y1, x2, y2] = arrowSegmentAtFraction_(xCurve, yCurve, frac, backDist)
% Return an arrow segment whose tip lies at a chosen fraction of curve length,
% and whose base lies backDist earlier along the same curve.

    frac = max(0, min(1, frac));

    ds = hypot(diff(xCurve), diff(yCurve));
    if isempty(ds)
        x1 = NaN; y1 = NaN; x2 = NaN; y2 = NaN;
        return;
    end

    sCum = [0; cumsum(ds(:))];
    totalLen = sCum(end);

    if totalLen <= 0
        x1 = NaN; y1 = NaN; x2 = NaN; y2 = NaN;
        return;
    end

    sTip  = frac * totalLen;
    sBase = max(0, sTip - backDist);

    [x2, y2] = pointOnPolylineAtArc_(xCurve, yCurve, sCum, sTip);
    [x1, y1] = pointOnPolylineAtArc_(xCurve, yCurve, sCum, sBase);
end

% =====================================================================
function [x, y] = pointOnPolylineAtArc_(xCurve, yCurve, sCum, sTarget)
% Interpolate the point lying at arc-length sTarget along a polyline.

    n = numel(xCurve);

    if sTarget <= 0
        x = xCurve(1);
        y = yCurve(1);
        return;
    elseif sTarget >= sCum(end)
        x = xCurve(end);
        y = yCurve(end);
        return;
    end

    idx = find(sCum <= sTarget, 1, 'last');
    idx = min(max(idx, 1), n-1);

    segLen = sCum(idx+1) - sCum(idx);
    if segLen <= 0
        x = xCurve(idx);
        y = yCurve(idx);
        return;
    end

    t = (sTarget - sCum(idx)) / segLen;
    x = xCurve(idx) + t * (xCurve(idx+1) - xCurve(idx));
    y = yCurve(idx) + t * (yCurve(idx+1) - yCurve(idx));
end

function s = ringDirectionSign_(i, j, theta)
% Return +1 / -1 for direction around the ring, with a stable tie-break
% for antipodal pairs.

    tol = 1e-12;
    d = mod(theta(j) - theta(i), 2*pi);   % in [0, 2pi)

    if d < tol || abs(d - 2*pi) < tol
        s = 0;
    elseif d < pi - tol
        s = +1;
    elseif d > pi + tol
        s = -1;
    else
        % exactly opposite: break tie by node index so reverse edge flips sign
        if i < j
            s = +1;
        else
            s = -1;
        end
    end
end

function sOut = wrapRingLabel_(sIn, maxChars)
% Wrap label at word boundaries to roughly maxChars per line.

    if isstring(sIn)
        sIn = char(sIn);
    end
    sIn = strtrim(sIn);

    if numel(sIn) <= maxChars || ~contains(sIn, ' ')
        sOut = sIn;
        return;
    end

    words = strsplit(sIn, ' ');
    lines = {};
    cur = words{1};

    for k = 2:numel(words)
        candidate = [cur ' ' words{k}];
        if strlength(string(candidate)) <= maxChars
            cur = candidate;
        else
            lines{end+1} = cur; %#ok<AGROW>
            cur = words{k};
        end
    end
    lines{end+1} = cur; %#ok<AGROW>

    sOut = strjoin(lines, newline);
end

function expandAxesToLabelExtents_(ax, labelHs, padFrac)
% Expand x/y limits so all label extents fit inside the plot box.

    labelHs = labelHs(isgraphics(labelHs, 'text'));
    if isempty(labelHs)
        return;
    end

    drawnow;

    xl = xlim(ax);
    yl = ylim(ax);

    xmin = xl(1);
    xmax = xl(2);
    ymin = yl(1);
    ymax = yl(2);

    for k = 1:numel(labelHs)
        oldUnits = labelHs(k).Units;
        labelHs(k).Units = 'data';
        ext = labelHs(k).Extent;   % [x y w h]
        labelHs(k).Units = oldUnits;

        xmin = min(xmin, ext(1));
        xmax = max(xmax, ext(1) + ext(3));
        ymin = min(ymin, ext(2));
        ymax = max(ymax, ext(2) + ext(4));
    end

    dx = xmax - xmin;
    dy = ymax - ymin;
    padX = padFrac * max(dx, eps);
    padY = padFrac * max(dy, eps);

    xlim(ax, [xmin - padX, xmax + padX]);
    ylim(ax, [ymin - padY, ymax + padY]);
end