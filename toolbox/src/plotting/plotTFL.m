function plotTFL(W, X, Y, h, varargin)
%PLOTTFL  Plot a trophic layout with bands, arrows, labels, BFF, lane-splitting,
%         and optional edge-weight-based line widths.
%
%   plotTFL(W, X, Y, h) plots a directed network with adjacency matrix W,
%   node coordinates X, Y (N x 1) and trophic levels h (N x 1), and:
%       - draws faint horizontal bands for integer trophic levels,
%       - plots directed edges (with arrowheads),
%       - styles upward vs downward edges differently,
%       - uses curved edges for bidirectional links between a pair of nodes,
%       - uses curved edges for unidirectional links between nodes at the
%         same trophic level,
%       - optionally "lane-splits" near-collinear shared-endpoint edges,
%       - optionally maps edge weights to edge widths,
%       - plots nodes on top, optional labels.
%
% Author:
%   Bazil Sansom (Warwick Business School, University of Warwick)
%   Contact: bazil.sansom@wbs.ac.uk

% ---- Basic checks ----
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
    error('plotTFL:NonFiniteXYH', 'X, Y, and h must be finite.');
end

% ---- Parse options ----
opts = parsePlotOpts(varargin{:});


% ---- Setup figure/axes ----
ax = opts.Parent;
if isempty(ax) || ~isgraphics(ax,'axes')
    ax = gca;
end

holdState = ishold(ax);

if opts.ClearAxes
    cla(ax);
end
hold(ax,'on');


% ---- Optional: scale Y to fill a square box (without distorting X) ----
Xplot = X;
Yplot = Y;

spanX = rangeNonzero(Xplot);
spanY = rangeNonzero(Yplot);

fillScale = 1;
if opts.FillSquare && spanY > 1e-12
    fillScale = spanX / spanY;   % makes Y span ~= X span
end

yScaleTotal = fillScale * opts.YScale;

if abs(yScaleTotal - 1) > 1e-12
    y0 = mean(Yplot);
    Yplot = y0 + (Yplot - y0) * yScaleTotal;
end

Y=Yplot; X=Xplot;  % use scaled versions from here on


% ---- Compute band positions (linear mapping Y ≈ a + b*h) ----
rangeH = max(h) - min(h);
if rangeH < 1e-12
    bands_y = mean(Y);
else
    H = [h, ones(N,1)];
    p = H \ Y;
    b = p(1);
    a = p(2);

    h_min_raw = floor(min(h));
    h_max_raw = ceil(max(h));
    ks        = h_min_raw:h_max_raw;
    bands_y   = a + b * ks;
end

% ---- Draw bands (behind everything else) ----
if opts.ShowBands
    xlimsBands = getBandXLimits(X);
    for k = 1:numel(bands_y)
        yk = bands_y(k);
        line(ax,xlimsBands, [yk yk], ...
            'Color', opts.BandColor, ...
            'LineStyle', opts.BandStyle, ...
            'LineWidth', opts.BandWidth);
    end
end

% ---- Extract edges ----
[iList, jList] = find(W ~= 0);
E = numel(iList);

% Edge index map: (i,j) -> edge index e
eIdx = sparse(iList, jList, 1:E, N, N);

% Reverse-edge lookup for each edge (0 if none)
if E > 0
    eRev = full(eIdx(jList, iList));    % E x 1
else
    eRev = zeros(0,1);
end
isMutualEdge = (eRev > 0);

% ---- Optional: map edge weights to line widths (compute once) ----
if E > 0
    edgeW  = abs(full(W(sub2ind([N N], iList, jList))));
else
    edgeW  = zeros(0,1);
end

edgeLW = opts.EdgeWidth * ones(E,1); % default fixed
if strcmpi(opts.EdgeWidthMode, 'weight') && E > 0
    edgeLW = mapWeightsToLineWidth(edgeW, opts.EdgeWidthRange, ...
        opts.EdgeWidthScale, opts.EdgeWidthQuantiles, opts.EdgeWidth);
end

%---- Decide node size scaling for edges if needed ----

nodeAutoOn = strcmpi(opts.NodeSizeMode,'auto');

edgeScale = 1;
if strcmpi(opts.EdgeWidthScaleMode,'auto') || ...
   (strcmpi(opts.EdgeWidthScaleMode,'withnode') && nodeAutoOn)

    Nref = max(1, opts.EdgeWidthRefN);
    edgeScale = (Nref / N) ^ opts.EdgeWidthGamma;
    edgeScale = min(max(edgeScale, opts.EdgeWidthMinScale), opts.EdgeWidthMaxScale);
end

edgeLW = edgeLW * edgeScale;
edgeLW = min(max(edgeLW, opts.EdgeWidthMin), opts.EdgeWidthMax);


% ---- Typical geometric scales ----
tolLen = 1e-12;  % geometric tolerance in data units
spanMax = max(rangeNonzero(X), rangeNonzero(Y));

if E > 0
    edgeLens = hypot(X(jList)-X(iList), Y(jList)-Y(iList));
    Ltyp = localRobustMedianPositive(edgeLens, tolLen);
else
    Ltyp = NaN;
end
if ~isfinite(Ltyp) || Ltyp < tolLen, Ltyp = spanMax; end
if ~isfinite(spanMax) || spanMax < tolLen, spanMax = 1; end

% ---- Arrow scaling (complexity-aware) ----
arrowScale = 1;
if strcmpi(opts.ArrowSizeMode,'auto') || ...
   (strcmpi(opts.ArrowSizeMode,'withnode') && nodeAutoOn)

    Nref = max(1, opts.ArrowSizeRefN);
    arrowScale = (Nref / N) ^ opts.ArrowSizeGamma;
    arrowScale = min(max(arrowScale, opts.ArrowSizeMinScale), opts.ArrowSizeMaxScale);
end

% ---- Choose arrow basis length ----
switch lower(strtrim(opts.ArrowLenBasis))
    case 'ltyp'
        Lbase = Ltyp;
    otherwise % 'span'
        Lbase = spanMax;
end

% ---- Arrow length: ALL knobs are fractions of Lbase ----
arrowLen = (opts.ArrowSize * arrowScale) * Lbase;

if opts.ArrowMinFrac > 0
    arrowLen = max(arrowLen, opts.ArrowMinFrac * Lbase);
end
if opts.ArrowMaxFrac > 0
    arrowLen = min(arrowLen, opts.ArrowMaxFrac * Lbase);
end


% ---- Hit-node guardrail thresholds (computed once) ----
%spanMax = max(rangeNonzero(X), rangeNonzero(Y));  % we already compute this in arrow section
hitR    = opts.HitNodeRadiusFrac * spanMax + opts.HitNodeMarginFrac * spanMax;

% Robust typical edge length for min-length gating
hitMinLen = opts.HitNodeMinLenFrac * Ltyp;


% Categorise edges: upward vs downward in trophic level
isUp = false(E,1);
if E > 0
    isUp = (h(jList) >= h(iList));
end

% Pre-blend colours with white to simulate alpha
bg      = [1 1 1];
colUp   = blendColor(opts.UpEdgeColor,   opts.EdgeAlpha, bg);
colDown = blendColor(opts.DownEdgeColor, opts.EdgeAlpha, bg);

% ---- Parameters for curvature ----
drawnE      = false(E,1);        % track edges already drawn by edge index
curvFracMut = 0.15;              % curvature as fraction of distance (mutual)
curvFracSame = 0.12;             % curvature for same-level edges
tolSame     = 1e-8;

% ---- Lane-splitting preprocessing (shared-endpoint collinearity) ----
offS = zeros(E,1);
offT = zeros(E,1);

if opts.LaneSplit && E > 0
    isSameEdge = abs(h(iList) - h(jList)) < tolSame;
    isExcluded = isMutualEdge | isSameEdge;

    P = [X Y];

    % Determine a robust typical edge length for MinLen defaulting
    edgeLens = hypot(X(jList)-X(iList), Y(jList)-Y(iList));
    medL = median(edgeLens(edgeLens > tolLen));
    if isempty(medL) || ~isfinite(medL), medL = median(edgeLens); end
    if isempty(medL) || ~isfinite(medL) || medL < tolLen, medL = 1; end

    minLenEff = opts.LaneSplitMinLen;
    if minLenEff <= 0
        minLenEff = 1e-6 * medL;
    else
        minLenEff = max(minLenEff, 1e-6 * medL);
    end

    [offS, offT] = computeLaneOffsets(P, iList, jList, isExcluded, ...
        'Theta', opts.LaneSplitTheta, ...
        'LaneSpacing', opts.LaneSplitSpacing, ...
        'MinLen', minLenEff, ...
        'h', h, ...
        'W', W, ...
        'PrimaryStraight', opts.LaneSplitPrimaryStraight, ...
        'PrimaryRule', opts.LaneSplitPrimaryRule, ...
        'OneSided', opts.LaneSplitOneSided);

    if opts.LaneSplitCoherent
        [offS, offT] = enforceCoherentOffsets(offS, offT, opts.LaneSplitCoherentTol);
    end

    if opts.LaneSplitPreferDownSide
        isDown = (h(jList) < h(iList));
        sgn = -1; % choose -1 or +1 as the "return-flow" side
        offS(isDown) = sgn * abs(offS(isDown));
        offT(isDown) = sgn * abs(offT(isDown));
    end
end

% ---- Handle bidirectional edges with curves (mutual from topology, robust) ----
% Uses mutual = (W~=0) & (W.'~=0), then uses eIdx to map back to edge indices.
% Convention: the DOWN edge curves "downwards" (bulge has negative y).

Wdir   = (W ~= 0);
mutual = Wdir & Wdir.';     % mutual(i,j)=true iff i->j and j->i both exist

for i = 1:N
    for j = i+1:N
        if ~mutual(i,j)
            continue;
        end

        e_ij = full(eIdx(i,j));
        e_ji = full(eIdx(j,i));

        % Safety: should exist, but guard anyway
        if e_ij == 0 || e_ji == 0
            continue;
        end
        if drawnE(e_ij) || drawnE(e_ji)
            continue;
        end

        % Decide which direction is "up" vs "down"
        % Default: up if h increases (ties: treat as up by using Y as tie-break)
        if (h(j) > h(i)) || (abs(h(j)-h(i)) < tolSame && Y(j) >= Y(i))
            low  = i; high = j;
        else
            low  = j; high = i;
        end

        % Edge indices for up/down (low->high is "up")
        e_up   = full(eIdx(low, high));
        e_down = full(eIdx(high, low));

        % If something weird, skip
        if e_up == 0 || e_down == 0
            continue;
        end

        x1 = X(low);  y1 = Y(low);
        x2 = X(high); y2 = Y(high);

        dx  = x2 - x1;
        dy  = y2 - y1;
        len = hypot(dx, dy);

        if len < tolLen
            drawnE([e_up, e_down]) = true;
            continue;
        end

        % Unit normal to low->high segment
        nx = -dy / len;
        ny =  dx / len;

        % We draw UP on +n and DOWN on -n.
        % Want DOWN bulge to be "downwards": (-n)_y < 0  =>  n_y > 0
        if ny <= 0
            nx = -nx;
            ny = -ny;
        end

        c = curvFracMut * len;
        s = linspace(0,1,25);

        % ---- UP edge: low -> high on +n ----
        xCurveU = x1 + dx.*s + c*sin(pi*s).*nx;
        yCurveU = y1 + dy.*s + c*sin(pi*s).*ny;
        drawCurvedArrow(ax,xCurveU, yCurveU, colUp, edgeLW(e_up), arrowLen, opts.ArrowAngle);

        % ---- DOWN edge: high -> low on -n ----
        dx2 = x1 - x2;
        dy2 = y1 - y2;
        xCurveD = x2 + dx2.*s - c*sin(pi*s).*nx;
        yCurveD = y2 + dy2.*s - c*sin(pi*s).*ny;
        drawCurvedArrow(ax,xCurveD, yCurveD, colDown, edgeLW(e_down), arrowLen, opts.ArrowAngle);

        % Mark both as drawn so they won't be redrawn later as straight/lane-split
        drawnE([e_up, e_down]) = true;
    end
end

%{
% ---- Handle bidirectional edges with curves (robust + consistent) ----
for e = 1:E
    eR = eRev(e);
    if eR == 0, continue; end
    if eR == e, continue; end
    if e > eR, continue; end                 % each mutual pair once
    if drawnE(e) || drawnE(eR), continue; end

    i  = iList(e);    j  = jList(e);
    iR = iList(eR);   jR = jList(eR);

    % Sanity: e should be i->j and eR should be j->i
    % If not, just skip to avoid weirdness (or swap, but this should hold).
    if ~(iR == j && jR == i)
        continue;
    end

    % Identify which directed edge is "up" vs "down" using trophic levels
    % (ties treated as up, so down only if strictly decreasing).
    if h(j) >= h(i)
        eUp = e;   eDn = eR;
    else
        eUp = eR;  eDn = e;
    end

    % Geometry: base segment for normal is the UP direction (tail->head of up edge)
    a = iList(eUp);  b = jList(eUp);
    x1 = X(a); y1 = Y(a);
    x2 = X(b); y2 = Y(b);

    dx  = x2 - x1;
    dy  = y2 - y1;
    len = hypot(dx, dy);
    if len < tolLen
        drawnE(e)  = true;
        drawnE(eR) = true;
        continue;
    end

    % Unit normal to the UP segment
    nx = -dy / len;
    ny =  dx / len;

    % Choose sign so DOWN curves downward: ny > 0
    if ny <= 0
        nx = -nx; ny = -ny;
    end
    
    % If ny is tiny, curvature is mostly horizontal and can look straight.
    % % Bias the normal so there is always a visible vertical bulge.
    minNy = 0.5;  % tune 0.15–0.35
    if abs(ny) < minNy
        ny = minNy;
        nx = sign(nx + 1e-12) * sqrt(max(0, 1 - ny^2));
    end

    c = curvFracMut * len;
    s = linspace(0,1,25);

    % ---- Draw UP edge as +n curve ----
    xCurveU = x1 + dx.*s + c*sin(pi*s).*nx;
    yCurveU = y1 + dy.*s + c*sin(pi*s).*ny;

    colU = colUp;
    lwU  = edgeLW(eUp);
    drawCurvedArrow(xCurveU, yCurveU, colU, lwU, arrowLen, opts.ArrowAngle);

    % ---- Draw DOWN edge as -n curve (reverse direction of same geometric segment) ----
    % Down is b->a
    dx2 = x1 - x2;
    dy2 = y1 - y2;

    xCurveD = x2 + dx2.*s - c*sin(pi*s).*nx;
    yCurveD = y2 + dy2.*s - c*sin(pi*s).*ny;

    colD = colDown;
    lwD  = edgeLW(eDn);
    drawCurvedArrow(xCurveD, yCurveD, colD, lwD, arrowLen, opts.ArrowAngle);

    %drawnE(e)  = true;
    %drawnE(eR) = true;
    drawnE([e, eR, eUp, eDn]) = true;


end
%}

% ---- Plot remaining edges (non-bidirectional) ----
for e = 1:E
    if drawnE(e)
        continue;
    end

    i = iList(e);
    j = jList(e);

    x1 = X(i); y1 = Y(i);
    x2 = X(j); y2 = Y(j);

    % Choose colour based on direction
    col = colUp;
    if ~isUp(e), col = colDown; end

    % Choose line width (precomputed)
    lw = edgeLW(e);

    % Same-level edges: curved
    if abs(h(i) - h(j)) < tolSame
        dx  = x2 - x1;
        dy  = y2 - y1;
        len = hypot(dx, dy);
        if len < tolLen
            plot(ax,[x1 x2], [y1 y2], '-', 'Color', col, 'LineWidth', lw);
            drawArrowhead(ax,x1, y1, x2, y2, col, arrowLen, opts.ArrowAngle);
            continue;
        end

        % Normal
        ty = dy / len;
        nx = -ty;
        ny =  dx / len;

        c = curvFracSame * len;
        s = linspace(0, 1, 25);

        xCurve = x1 + dx.*s + c * sin(pi*s) .* nx;
        yCurve = y1 + dy.*s + c * sin(pi*s) .* ny;

        drawCurvedArrow(ax,xCurve, yCurve, col, lw, arrowLen, opts.ArrowAngle);
        continue;
    end

    % ---- Hit-node guardrail: curve long edges that pass near other nodes ----
    if opts.HitNode
        dx0 = x2 - x1;
        dy0 = y2 - y1;
        len0 = hypot(dx0, dy0);

        if len0 >= max(hitMinLen, tolLen)
            [hits, kHit, pClosest] = edgeHitsAnyNode(i, j, [x1 y1], [x2 y2], X, Y, hitR);

            if hits
                % Choose curvature direction to bulge AWAY from offending node kHit
                % Compute unit normal to edge
                tx = dx0 / len0;
                ty = dy0 / len0;
                nx = -ty;
                ny =  tx;

                % Vector from closest point on segment to the node
                vkx = X(kHit) - pClosest(1);
                vky = Y(kHit) - pClosest(2);

                % If node lies on +normal side, curve to -normal, else +normal
                sgn = sign(vkx*nx + vky*ny);
                if sgn == 0, sgn = 1; end

                c = opts.HitNodeCurvFrac * len0;
                s = linspace(0,1,opts.HitNodeNpts);

                xCurve = x1 + dx0.*s - sgn * c*sin(pi*s).*nx;
                yCurve = y1 + dy0.*s - sgn * c*sin(pi*s).*ny;

                drawCurvedArrow(ax, xCurve, yCurve, col, lw, arrowLen, opts.ArrowAngle);
                continue;
            end
        end
    end

    
    % Lane-splitting for near-collinear shared-endpoint edges (otherwise straight)
    if opts.LaneSplit && (offS(e) ~= 0 || offT(e) ~= 0)
        dx  = x2 - x1;
        dy  = y2 - y1;
        len = hypot(dx, dy);
        if len < tolLen
            plot(ax,[x1 x2], [y1 y2], '-', 'Color', col, 'LineWidth', lw);
            drawArrowhead(ax, x1, y1, x2, y2, col, arrowLen, opts.ArrowAngle);
            continue;
        end

        % Unit normal
        tx = dx / len;
        ty = dy / len;
        nx = -ty;
        ny =  tx;

        % Cubic Bezier control points with endpoint-local offsets
        alpha = opts.LaneSplitAlpha;
        p0 = [x1 y1];
        p3 = [x2 y2];
        d  = [dx dy];

        offSe = offS(e);
        offTe = offT(e);

        if opts.LaneSplitScaleWithEdgeLength
            offSe = offSe * len;
            offTe = offTe * len;
        end

        p1 = p0 + alpha * d + [nx ny] * offSe;
        p2 = p0 + (1-alpha) * d + [nx ny] * offTe;

        [xCurve, yCurve] = cubicBezier(p0, p1, p2, p3, opts.LaneSplitNpts);
        drawCurvedArrow(ax, xCurve, yCurve, col, lw, arrowLen, opts.ArrowAngle);
        continue;
    end

    % Default: straight edge with arrow
    plot(ax,[x1 x2], [y1 y2], '-', 'Color', col, 'LineWidth', lw);
    drawArrowhead(ax, x1, y1, x2, y2, col, arrowLen, opts.ArrowAngle);
end

% ---- Plot nodes ----
% Adjust node size if in 'auto' mode
nodeSize = opts.NodeSize;

if strcmpi(opts.NodeSizeMode,'auto')
    Nref = max(1, opts.NodeSizeRefN);
    scale = (Nref / N) ^ opts.NodeSizeGamma;
    nodeSize = nodeSize * scale;
    nodeSize = min(max(nodeSize, opts.NodeSizeMin), opts.NodeSizeMax);
end

scatter(ax, X, Y, nodeSize, ...
    'MarkerFaceColor', opts.NodeColor, ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 0.5);


% ---- Plot node labels (optional) ----
if opts.ShowLabels
    if isempty(opts.Labels)
        lbls = arrayfun(@(k) sprintf('%d',k), 1:N, 'UniformOutput', false);
    else
        lbls = opts.Labels;
        if numel(lbls) ~= N
            error('plotTFL:BadLabels', 'Labels must be a cell array of N strings.');
        end
    end

    dx = opts.LabelOffset(1);
    dy = opts.LabelOffset(2);

    for ii = 1:N
        text(ax,X(ii) + dx, Y(ii) + dy, lbls{ii}, ...
            'FontSize',  opts.LabelFontSize, ...
            'Color',     opts.LabelColor, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment',   'bottom');
    end
end

% ---- Axes scaling and limits ----
rangeX   = rangeNonzero(X);
rangeY   = rangeNonzero(Y);

basePadX = 0.05 * rangeX;
basePadY = 0.05 * rangeY;

curvMax = max(curvFracMut, curvFracSame);
padX = basePadX + curvMax * rangeY;
padY = basePadY + curvMax * rangeX;

xlim(ax, [min(X) - padX, max(X) + padX]);
ylim(ax, [min(Y) - padY, max(Y) + padY]);

daspect(ax, [1 1 1]);   % equal data units
pbaspect(ax, [1 1 1]);  % square plotting box


% Clean diagram
set(ax, 'XTick', [], 'YTick', []);
set(ax, 'XTickLabel', [], 'YTickLabel', []);
box(ax,'on');

% ---- Add title if requested ----
if ~isempty(opts.Title)
    title(ax, opts.Title, 'Interpreter', 'none');
end

% ---- Add BFF label if requested ----
if opts.ShowBFF && ~isnan(opts.BFF)
    xl = xlim(ax);
    yl = ylim(ax);
    xText = xl(1) + 0.02 * (xl(2) - xl(1));
    yText = yl(2) - 0.05 * (yl(2) - yl(1));
    text(ax,xText, yText, sprintf('BFF = %.3f', opts.BFF), ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment',   'top', ...
        'FontSize',            10, ...
        'BackgroundColor',     'w', ...
        'Margin',              2);
end

if ~holdState
    hold(ax,'off');
end

end % main function

% =====================================================================
function opts = parsePlotOpts(varargin)
%PARSEPLOTOPTS  Parse name/value options for plotTFL.

% Default options
opts.Parent = [];      % axes handle to draw into (default: gca)
opts.ClearAxes = true; % keep existing behaviour unless overridden

opts.FillSquare = true;   % make the drawing fill a square window
opts.YScale     = 1;      % additional manual multiplier (1 = none)


opts.ShowBands     = true;
opts.BandColor     = [0.75 0.75 0.75];
opts.BandStyle     = ':';
opts.BandWidth     = 0.75;

opts.NodeColor     = [0.1 0.3 0.7];

opts.NodeSize      = 30;      % base size at NodeSizeRefN
opts.NodeSizeMode  = 'auto'; % 'fixed' or 'auto'
opts.NodeSizeRefN  = 20;      % reference N for NodeSize
opts.NodeSizeGamma = 0.3;     % mild scaling
opts.NodeSizeMin   = 10;
opts.NodeSizeMax   = 80;


opts.UpEdgeColor   = [0 0 0];
opts.DownEdgeColor = [0.6 0.6 0.6];
opts.EdgeWidth     = 0.75;
opts.EdgeAlpha     = 0.7;

opts.ArrowLenBasis = 'span'; % length scale Lbase for arrows: 'span'->spanMax, 'ltyp'->median edge length
opts.ArrowSize     = 0.04;   % base arrow length as fraction of Lbase
opts.ArrowMinFrac  = 0.02;   % min arrow length as fraction of Lbase
opts.ArrowMaxFrac  = 0.06;   % max arrow length as fraction of Lbase (0 disables)
opts.ArrowAngle    = pi/7;


% ---- Arrow scaling (complexity-aware) ----
% ArrowSizeMode controls an extra N-based multiplier arrowScale:
%   'fixed'    : arrowScale = 1
%   'auto'     : arrowScale = clip((ArrowSizeRefN/N)^ArrowSizeGamma, MinScale, MaxScale)
%   'withNode' : same as 'auto' but only when NodeSizeMode='auto' (else arrowScale=1)
opts.ArrowSizeMode      = 'withNode';  % 'fixed'|'auto'|'withNode'
opts.ArrowSizeRefN      = 40;
opts.ArrowSizeGamma     = 0.5;
opts.ArrowSizeMinScale  = 0.35;
opts.ArrowSizeMaxScale  = 1.2;

% ---- Edge width scaling (complexity-aware) ----
opts.EdgeWidthScaleMode     = 'withNode'; % 'fixed'|'auto'|'withNode'
opts.EdgeWidthRefN          = 40;
opts.EdgeWidthGamma         = 0.35;
opts.EdgeWidthMinScale      = 0.5;
opts.EdgeWidthMaxScale      = 1.2;
opts.EdgeWidthMin           = 0.25;       % post-scale clamp
opts.EdgeWidthMax           = 3.0;

opts.ShowLabels    = false;
opts.Labels        = {};
opts.LabelFontSize = 9;
opts.LabelColor    = [0 0 0];
opts.LabelOffset   = [0, 0.02];

opts.Title         = '';
opts.ShowBFF       = false;
opts.BFF           = NaN;

% ---- Hit-node guardrail (curve edges that pass near other nodes) ----
opts.HitNode            = true;   % enable/disable
opts.HitNodeRadiusFrac  = 0.030;  % node "keep-out" radius as fraction of plot span
opts.HitNodeMarginFrac  = 0.005;  % extra margin as fraction of plot span
opts.HitNodeCurvFrac    = 0.18;   % curvature strength as fraction of edge length
opts.HitNodeNpts        = 25;     % samples along curve
opts.HitNodeMinLenFrac  = 0.15;   % only apply to edges longer than this * median edge length

% ---- Lane splitting (shared-endpoint collinearity) ----
opts.LaneSplit        = true;
opts.LaneSplitScaleWithEdgeLength = true;
opts.LaneSplitTheta   = deg2rad(6);
opts.LaneSplitSpacing = 0.15;
opts.LaneSplitAlpha   = 0.25;
opts.LaneSplitNpts    = 25;
opts.LaneSplitMinLen  = 0;     % if 0, derived from median edge length
opts.LaneSplitCoherent = true;
opts.LaneSplitCoherentTol = 0.25;
opts.LaneSplitPreferDownSide = true;

opts.LaneSplitPrimaryStraight = true;
opts.LaneSplitPrimaryRule     = 'upweight'; % 'up'|'weight'|'length'|'upweight'
opts.LaneSplitOneSided        = true;

% ---- Edge weighting options ----
opts.EdgeWidthMode      = 'fixed';     % 'fixed' or 'weight'
opts.EdgeWidthRange     = [0.6 2.5];   % [min max] linewidth in points
opts.EdgeWidthScale     = 'log';       % 'log' or 'linear'
opts.EdgeWidthQuantiles = [0.05 0.95]; % robust clipping for weight mapping

if isempty(varargin), return; end
if mod(numel(varargin),2) ~= 0
    error('plotTFL:BadArgs', 'Optional arguments must be name/value pairs.');
end

for k = 1:2:numel(varargin)
    name  = varargin{k};
    value = varargin{k+1};
    if ~ischar(name) && ~isstring(name)
        error('plotTFL:BadParamName', 'Parameter names must be strings.');
    end

    switch lower(char(name))
        case 'parent'
            opts.Parent = value;
        case 'clearaxes'
            opts.ClearAxes = logical(value);

        case 'fillsquare'
            opts.FillSquare = logical(value);
        case 'yscale'
            opts.YScale = value;
        
        case 'showbands'
            opts.ShowBands = logical(value);
        case 'bandcolor'
            opts.BandColor = value;
        case 'bandstyle'
            opts.BandStyle = value;
        case 'bandwidth'
            opts.BandWidth = value;

        case 'nodesize'
            opts.NodeSize = value;
        case 'nodesizemode'
            opts.NodeSizeMode = char(value);
        case 'nodesizerefn'
            opts.NodeSizeRefN = value;
        case 'nodesizegamma'
            opts.NodeSizeGamma = value;
        case 'nodesizemin'
            opts.NodeSizeMin = value;
        case 'nodesizemax'
            opts.NodeSizeMax = value;
        case 'nodecolor'
            opts.NodeColor = value;

        case 'upedgecolor'
            opts.UpEdgeColor = value;
        case 'downedgecolor'
            opts.DownEdgeColor = value;
        case 'edgewidth'
            opts.EdgeWidth = value;
        case 'edgealpha'
            opts.EdgeAlpha = value;

        case 'arrowlenbasis'
            opts.ArrowLenBasis = char(value);
        case 'arrowsize'
            opts.ArrowSize = value;
        case 'arrowangle'
            opts.ArrowAngle = value;
        case 'arrowminfrac'
            opts.ArrowMinFrac = value;

        case 'arrowsizemode',     opts.ArrowSizeMode = char(value);
        case 'arrowsizerefn',     opts.ArrowSizeRefN = value;
        case 'arrowsizegamma',    opts.ArrowSizeGamma = value;
        case 'arrowsizeminscale', opts.ArrowSizeMinScale = value;
        case 'arrowsizemaxscale', opts.ArrowSizeMaxScale = value;
        case 'arrowmaxfrac',      opts.ArrowMaxFrac = value;

        case 'edgewidthscalemode', opts.EdgeWidthScaleMode = char(value);
        case 'edgewidthrefn',      opts.EdgeWidthRefN = value;
        case 'edgewidthgamma',     opts.EdgeWidthGamma = value;
        case 'edgewidthminscale',  opts.EdgeWidthMinScale = value;
        case 'edgewidthmaxscale',  opts.EdgeWidthMaxScale = value;
        case 'edgewidthmin',       opts.EdgeWidthMin = value;
        case 'edgewidthmax',       opts.EdgeWidthMax = value;

        case 'showlabels'
            opts.ShowLabels = logical(value);
        case 'labels'
            opts.Labels = value;
        case 'labelfontsize'
            opts.LabelFontSize = value;
        case 'labelcolor'
            opts.LabelColor = value;
        case 'labeloffset'
            opts.LabelOffset = value;

        case 'title'
            opts.Title = value;
        case 'showbff'
            opts.ShowBFF = logical(value);
        case 'bff'
            opts.BFF = value;

        % ---- Hit-node guardrail options ----
        case 'hitnode'
            opts.HitNode = logical(value);
        case 'hitnoderadiusfrac'
            opts.HitNodeRadiusFrac = value;
        case 'hitnodemarginfrac'
            opts.HitNodeMarginFrac = value;
        case 'hitnodecurvfrac'
            opts.HitNodeCurvFrac = value;
        case 'hitnodenpts'
            opts.HitNodeNpts = value;
        case 'hitnodeminlenfrac'
            opts.HitNodeMinLenFrac = value;

        % ---- Lane splitting options ----
        case 'lanesplit'
            opts.LaneSplit = logical(value);
        case 'lanesplitscalewithedgelength'
            opts.LaneSplitScaleWithEdgeLength = logical(value);
        case 'lanesplittheta'
            opts.LaneSplitTheta = value;
        case 'lanesplitspacing'
            opts.LaneSplitSpacing = value;
        case 'lanesplitalpha'
            opts.LaneSplitAlpha = value;
        case 'lanesplitnpts'
            opts.LaneSplitNpts = value;
        case 'lanesplitminlen'
            opts.LaneSplitMinLen = value;
        case 'lanesplitcoherent'
            opts.LaneSplitCoherent = logical(value);
        case 'lanesplitcoherenttol'
            opts.LaneSplitCoherentTol = value;
        case 'lanesplitpreferdownside'
            opts.LaneSplitPreferDownSide = logical(value);
        case 'lanesplitprimarystraight'
            opts.LaneSplitPrimaryStraight = logical(value);
        case 'lanesplitprimaryrule'
            opts.LaneSplitPrimaryRule = char(value);
        case 'lanesplitonesided'
            opts.LaneSplitOneSided = logical(value);

        % ---- Edge weighting options ----
        case 'edgewidthmode'
            opts.EdgeWidthMode = char(value);
        case 'edgewidthrange'
            opts.EdgeWidthRange = value;
        case 'edgewidthscale'
            opts.EdgeWidthScale = char(value);
        case 'edgewidthquantiles'
            opts.EdgeWidthQuantiles = value;

        otherwise
            error('plotTFL:UnknownOption', 'Unknown option: %s', name);
    end
end

end

% =====================================================================
function xlims = getBandXLimits(X)
pad   = 0.05 * rangeNonzero(X);
xlims = [min(X)-pad, max(X)+pad];
end

function r = rangeNonzero(x)
r = max(x) - min(x);
if r == 0
    r = 1;
end
end

function cBlend = blendColor(c, alpha, bg)
c  = c(:)';
bg = bg(:)';
cBlend = alpha * c + (1-alpha) * bg;
end

function drawArrowhead(ax, x1, y1, x2, y2, col, L, halfAngle)
dx  = x2 - x1;
dy  = y2 - y1;
len = hypot(dx, dy);
if len == 0
    return;
end

ux = dx / len;
uy = dy / len;

px = x2 - L * ux;
py = y2 - L * uy;

px_perp = -uy;
py_perp =  ux;

w = tan(halfAngle) * L;

bx1 = px + w * px_perp;
by1 = py + w * py_perp;
bx2 = px - w * px_perp;
by2 = py - w * py_perp;

patch(ax, [x2 bx1 bx2], [y2 by1 by2], col, ...
    'EdgeColor', 'none', ...
    'FaceColor', col);
end

function drawCurvedArrow(ax, xCurve, yCurve, col, edgeWidth, arrowSize, arrowAngle)
plot(ax, xCurve, yCurve, '-', 'Color', col, 'LineWidth', edgeWidth);

n = numel(xCurve);
if n < 2
    return;
end
x1 = xCurve(end-1); y1 = yCurve(end-1);
x2 = xCurve(end);   y2 = yCurve(end);

drawArrowhead(ax, x1, y1, x2, y2, col, arrowSize, arrowAngle);
end

% 
function [hits, kHit, pClosest] = edgeHitsAnyNode(i, j, p1, p2, X, Y, hitR)
%EDGEHITSANYNODE  True if segment p1->p2 passes within hitR of any node k != i,j.
%
% Returns:
%   hits     : logical
%   kHit     : index of closest offending node
%   pClosest : closest point on segment to that node

N = numel(X);
hits = false;
kHit = 0;
pClosest = [NaN NaN];

bestD = inf;

for k = 1:N
    if k == i || k == j
        continue;
    end
    pk = [X(k) Y(k)];
    [d, pc] = pointToSegmentDist(pk, p1, p2);
    if d < bestD
        bestD = d;
        kHit = k;
        pClosest = pc;
    end
end

hits = (bestD < hitR);
end


function [d, pClosest] = pointToSegmentDist(p, a, b)
%POINTTOSEGMENTDIST  Distance from point p to segment a-b in R^2.
% Returns d and closest point pClosest on the segment.

ax = a(1); ay = a(2);
bx = b(1); by = b(2);
px = p(1); py = p(2);

abx = bx - ax; aby = by - ay;
apx = px - ax; apy = py - ay;

ab2 = abx*abx + aby*aby;
if ab2 <= 0
    pClosest = a;
    d = hypot(px-ax, py-ay);
    return;
end

t = (apx*abx + apy*aby) / ab2;
t = max(0, min(1, t));

cx = ax + t*abx;
cy = ay + t*aby;

pClosest = [cx cy];
d = hypot(px-cx, py-cy);
end


% =====================================================================
% Lane-splitting helper: compute per-edge offsets at start/end nodes
function [offS, offT] = computeLaneOffsets(P, ei, ej, isExcluded, varargin)
% computeLaneOffsets
% ------------------
% Compute per-edge lateral offsets (start/end) to separate near-collinear
% edges that share an endpoint (ignoring direction).
%
% - Groups incident edges at each node by *folded* direction (mod pi),
%   so incoming/outgoing collinear edges are treated together.
% - Assigns small signed "lane" offsets within each group.
% - Optionally keeps a "primary" edge straight (offset 0) in each group.
% - Optionally keeps a degree-2 straight-through UP-chain node *both* lanes 0.

% ---- defaults ----
opts.Theta           = deg2rad(6);
opts.LaneSpacing     = 0.10;        % dimensionless fraction (scale by edge length when drawing)
opts.MinLen          = 1e-6;

opts.PrimaryStraight = true;
opts.PrimaryRule     = 'up';        % 'up' | 'weight' | 'length' | 'upweight'
opts.OneSided        = true;
opts.ChainStraight   = true;

opts.h = [];
opts.W = [];

opts = parseLocalOpts(opts, varargin{:});

E = numel(ei);
N = size(P,1);

offS = zeros(E,1);
offT = zeros(E,1);

% ---- build incident lists (excluding special-case edges) ----
inc = cell(N,1);
for k = 1:E
    if isExcluded(k), continue; end
    inc{ei(k)}(end+1) = k; %#ok<AGROW>
    inc{ej(k)}(end+1) = k; %#ok<AGROW>
end

theta = opts.Theta;

for v = 1:N
    edges = inc{v};
    if numel(edges) < 2
        continue;
    end

    % Folded angles in [0, pi)
    a    = zeros(numel(edges),1);
    keep = true(numel(edges),1);

    for idx = 1:numel(edges)
        k = edges(idx);
        if ei(k) == v
            u = ej(k);
        else
            u = ei(k);
        end

        d   = P(u,:) - P(v,:);
        len = hypot(d(1), d(2));
        if len < opts.MinLen
            keep(idx) = false;
            continue;
        end

        tv = d / len;
        ai = atan2(tv(2), tv(1));
        a(idx) = mod(ai, pi);   % fold direction so opposite vectors match
    end

    edges = edges(keep);
    a     = a(keep);

    if numel(edges) < 2
        continue;
    end

    % Sort by folded angle and build contiguous groups
    [aS, ord] = sort(a);
    edgesS = edges(ord);

    groups = {};
    g = 1;
    groups{g} = 1; %#ok<AGROW>
    for ii = 2:numel(edgesS)
        if (aS(ii) - aS(ii-1)) <= theta
            groups{g}(end+1) = ii;
        else
            g = g+1;
            groups{g} = ii; %#ok<AGROW>
        end
    end

    % Wraparound merge near 0/pi
    if numel(groups) >= 2 && (aS(1) + pi - aS(end)) <= theta
        groups{1} = [groups{end}, groups{1}];
        groups(end) = [];
    end

    % ---- assign lanes per group ----
    for gi = 1:numel(groups)
        idxs = groups{gi};
        mG = numel(idxs);
        if mG < 2
            continue;
        end

        eG = edgesS(idxs);

        % Deterministic ordering within group
        [~, perm] = sort(eG);
        eG = eG(perm);

        lanes = zeros(mG,1);

        % Special case: degree-2 straight-through chain at node v
        % (one edge incoming to v, one edge outgoing from v)
        if opts.ChainStraight && mG == 2 && ~isempty(opts.h)
            k1 = eG(1); k2 = eG(2);

            in1  = (ej(k1) == v); out1 = (ei(k1) == v);
            in2  = (ej(k2) == v); out2 = (ei(k2) == v);
            isThrough = (in1 && out2) || (out1 && in2);

            if isThrough
                % Only keep both straight if BOTH edges are UP edges
                up1 = (opts.h(ej(k1)) >= opts.h(ei(k1)));
                up2 = (opts.h(ej(k2)) >= opts.h(ei(k2)));

                if up1 && up2
                    lanes(:) = 0;
                    writeOffsetsForGroup(v, eG, lanes);
                    continue;
                end
            end
        end

        if opts.PrimaryStraight && mG >= 2
            % Choose primary edge index p0 (1..mG) to keep at lane 0
            p0 = choosePrimaryEdge(v, eG);

            if opts.OneSided
                % Put all non-primary edges on same side: +1, +2, ...
                q = 1;
                for p = 1:mG
                    if p == p0, continue; end
                    lanes(p) = q * opts.LaneSpacing;
                    q = q + 1;
                end
            else
                % Symmetric around 0, but keep primary at 0
                others = setdiff(1:mG, p0, 'stable');
                mO = numel(others);
                base = (-(mO-1)/2 : (mO-1)/2) * opts.LaneSpacing;
                lanes(others) = base(:);
                lanes(p0) = 0;
            end
        else
            % Fallback: classic symmetric lanes
            lanes = (-(mG-1)/2 : (mG-1)/2) * opts.LaneSpacing;
            lanes = lanes(:);
        end

        writeOffsetsForGroup(v, eG, lanes);
    end
end

% ---------- nested helpers ----------
    function writeOffsetsForGroup(vv, eGrp, laneVals)
        for pp = 1:numel(eGrp)
            kk = eGrp(pp);
            if ei(kk) == vv
                offS(kk) = offS(kk) + laneVals(pp);
            else
                offT(kk) = offT(kk) + laneVals(pp);
            end
        end
    end

    function p0 = choosePrimaryEdge(~, eGrp)
        mLoc = numel(eGrp);

        switch lower(strtrim(opts.PrimaryRule))
            case 'up'
                if ~isempty(opts.h)
                    isUp = false(mLoc,1);
                    for pp = 1:mLoc
                        kk = eGrp(pp);
                        isUp(pp) = (opts.h(ej(kk)) >= opts.h(ei(kk)));
                    end
                    cand = find(isUp);
                    if ~isempty(cand)
                        % Prefer the most "axial" among candidates
                        score = -inf(numel(cand),1);
                        for ii2 = 1:numel(cand)
                            pp = cand(ii2);
                            kk = eGrp(pp);
                            d = P(ej(kk),:) - P(ei(kk),:);
                            L = hypot(d(1), d(2));
                            if L > 0
                                score(ii2) = abs(d(2)) / L;
                            end
                        end
                        [~, jj] = max(score);
                        p0 = cand(jj);
                        return;
                    end
                end
                p0 = 1;

            case 'weight'
                if ~isempty(opts.W)
                    wgt = -inf(mLoc,1);
                    for pp = 1:mLoc
                        kk = eGrp(pp);
                        wgt(pp) = abs(opts.W(ei(kk), ej(kk)));
                    end
                    [~, p0] = max(wgt);
                else
                    p0 = 1;
                end

            case 'length'
                lenv = -inf(mLoc,1);
                for pp = 1:mLoc
                    kk = eGrp(pp);
                    d = P(ej(kk),:) - P(ei(kk),:);
                    lenv(pp) = hypot(d(1), d(2));
                end
                [~, p0] = max(lenv);

            case 'upweight'
                if ~isempty(opts.h)
                    isUp = false(mLoc,1);
                    wgt  = zeros(mLoc,1);

                    for pp = 1:mLoc
                        kk = eGrp(pp);
                        isUp(pp) = (opts.h(ej(kk)) >= opts.h(ei(kk)));
                        if ~isempty(opts.W)
                            wgt(pp) = abs(opts.W(ei(kk), ej(kk)));
                        else
                            wgt(pp) = 0;
                        end
                    end

                    cand = find(isUp);
                    if ~isempty(cand)
                        [~, jj] = max(wgt(cand));
                        p0 = cand(jj);
                        return;
                    end
                end

                % If no up edge exists, fall back to max weight (if available)
                if ~isempty(opts.W)
                    wgt = zeros(mLoc,1);
                    for pp = 1:mLoc
                        kk = eGrp(pp);
                        wgt(pp) = abs(opts.W(ei(kk), ej(kk)));
                    end
                    [~, p0] = max(wgt);
                else
                    p0 = 1;
                end

            otherwise
                p0 = 1;
        end
    end
end

function opts = parseLocalOpts(opts, varargin)
for i = 1:2:numel(varargin)
    opts.(varargin{i}) = varargin{i+1};
end
end

% =====================================================================
function [xCurve, yCurve] = cubicBezier(p0, p1, p2, p3, npts)
t = linspace(0,1,npts).';
b0 = (1-t).^3;
b1 = 3*(1-t).^2.*t;
b2 = 3*(1-t).*t.^2;
b3 = t.^3;
B  = b0.*p0 + b1.*p1 + b2.*p2 + b3.*p3;
xCurve = B(:,1).';
yCurve = B(:,2).';
end

% =====================================================================
function [offS2, offT2] = enforceCoherentOffsets(offS, offT, tol)
offS2 = offS; offT2 = offT;
m = numel(offS2);

for k = 1:m
    a = offS2(k);
    b = offT2(k);
    if a == 0 || b == 0
        continue;
    end

    % Opposite signs -> S-curve tendency
    if a*b < 0
        if abs(a + b) <= tol * max(abs(a), abs(b))
            if abs(a) < abs(b)
                offS2(k) = 0;
            else
                offT2(k) = 0;
            end
        else
            s = sign(a + b);
            if s == 0, s = sign(a); end
            mag = max(abs(a), abs(b));
            offS2(k) = s * mag;
            offT2(k) = s * mag;
        end
    end
end
end

% =====================================================================
function lw = mapWeightsToLineWidth(w, lwRange, scaleMode, qRange, fallbackLW)
%MAPWEIGHTSTOLINEWIDTH  Robust mapping from weights to line widths.

w = double(w(:));
lw = fallbackLW * ones(size(w));

if isempty(w) || all(~isfinite(w)) || all(w <= 0)
    return;
end

wpos = w(isfinite(w) & w > 0);
if isempty(wpos), return; end

qLo = qRange(1); qHi = qRange(2);
qLo = max(0, min(1, qLo));
qHi = max(0, min(1, qHi));
if qHi <= qLo
    qLo = 0.05; qHi = 0.95;
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
t = min(max(t,0),1);

lwMin = lwRange(1);
lwMax = lwRange(2);
lw = lwMin + t .* (lwMax - lwMin);
end


function m = localRobustMedianPositive(v, tol)
v = v(:);
v = v(isfinite(v) & v > tol);
if isempty(v)
    m = NaN;
else
    m = median(v);
end
end
