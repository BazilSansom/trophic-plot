function plotTFL(W, X, Y, h, varargin)
%PLOTTFL  Plot a trophic layout with bands, arrows, labels, and BFF.
%
%   plotTFL(W, X, Y, h) plots a directed network with adjacency matrix W,
%   node coordinates X, Y (N x 1) and trophic levels h (N x 1), and:
%       - draws faint horizontal bands for integer trophic levels,
%       - plots directed edges (with arrowheads),
%       - styles upward vs downward edges differently,
%       - uses curved edges for bidirectional links between a pair of nodes,
%       - uses curved edges for unidirectional links between nodes at the
%         same trophic level,
%       - plots nodes on top.
%
%   plotTFL(..., 'Name', Value, ...) supports:
%
%   Bands:
%     'ShowBands'     - logical, draw unit trophic bands (default: true)
%     'BandColor'     - 1x3 RGB for band lines (default: [0.75 0.75 0.75])
%     'BandStyle'     - line style for bands (default: ':')
%     'BandWidth'     - line width for bands (default: 0.75)
%
%   Nodes:
%     'NodeSize'      - marker size for nodes (default: 50)
%     'NodeColor'     - RGB for node face (default: [0.1 0.3 0.7])
%
%   Edges:
%     'UpEdgeColor'   - RGB for edges with h(j) >= h(i) (default: [0 0 0])
%     'DownEdgeColor' - RGB for edges with h(j) < h(i) (default: [0.6 0.6 0.6])
%     'EdgeWidth'     - base line width for edges (default: 0.75)
%     'EdgeAlpha'     - "transparency" in [0,1]; implemented by blending
%                       with white (default: 0.7)
%
%   Arrows:
%     'ArrowSize'     - dimensionless factor setting arrowhead length
%                       relative to the smaller of the X/Y spans
%                       (default: 0.06).
%     'ArrowAngle'    - half-angle of arrowhead in radians (default: pi/7)
%
%   Labels:
%     'ShowLabels'    - logical, show node labels (default: false)
%     'Labels'        - cell array of N strings; if empty and ShowLabels
%                       is true, uses '1','2',...,'N'
%     'LabelFontSize' - font size for labels (default: 9)
%     'LabelColor'    - RGB for label text (default: [0 0 0])
%     'LabelOffset'   - [dx dy] offset for labels in data units
%                       (default: [0, 0.02])
%
%   Title / BFF:
%     'Title'         - title string (default: '')
%     'ShowBFF'       - logical, annotate BFF value (default: false)
%     'BFF'           - numeric BFF value to display (default: NaN)
%
%   Notes:
%   - Mapping between h and Y is assumed linear: Y ≈ a + b*h.
%     This function estimates a, b by least squares and draws bands at
%     Y = a + b*k for integer k from floor(min(h)) to ceil(max(h)).
%   - BFF is passed in by the caller; this function does not compute it.
%   - For node pairs with edges in both directions (i->j and j->i),
%     two curved edges are drawn on opposite sides of the straight line,
%     each with its own colour and arrowhead.
%   - For single-direction edges between nodes at the same trophic level,
%     a single curved edge is drawn to avoid overlap/ambiguity.
%
% Author:
%   Bazil Sansom (Warwick Business School, University of Warwick)
%   Contact: bazil.sansom@wbs.ac.uk
%
% -------------------------------------------------------------------------

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

% ---- Parse options ----
opts = parsePlotOpts(varargin{:});

% ---- Compute band positions (linear mapping Y ≈ a + b*h) ----
rangeH = max(h) - min(h);

if rangeH < 1e-12
    % Degenerate case: essentially one trophic level
    bands_h = mean(h);
    bands_y = mean(Y);
else
    H = [h, ones(N,1)];
    p = H \ Y;
    b = p(1);    % step in Y per unit trophic level
    a = p(2);    % offset

    h_min_raw = floor(min(h));
    h_max_raw = ceil(max(h));
    ks        = h_min_raw:h_max_raw;    % integer trophic levels

    bands_h = ks;             %#ok<NASGU>  % labels in raw trophic units
    bands_y = a + b * ks;     % actual Y positions of bands
end

% ---- Setup figure/axes ----
holdState = ishold;
cla; hold on;

% ---- Draw bands (behind everything else) ----
if opts.ShowBands
    xlimsBands = getBandXLimits(X);
    for k = 1:numel(bands_y)
        yk = bands_y(k);
        line(xlimsBands, [yk yk], ...
            'Color', opts.BandColor, ...
            'LineStyle', opts.BandStyle, ...
            'LineWidth', opts.BandWidth);
    end
end

% ---- Extract edges ----
[iList, jList] = find(W ~= 0);
E = numel(iList);

% ---- Decide arrow length based on plot span ----
rangeX   = rangeNonzero(X);
rangeY   = rangeNonzero(Y);
baseSpan = min(rangeX, rangeY);    % use the smaller dimension for scale
arrowLen = opts.ArrowSize * baseSpan;

% Categorise edges: upward vs downward in trophic level
isUp = h(jList) >= h(iList);

% Pre-blend colours with white to simulate alpha
bg      = [1 1 1];
colUp   = blendColor(opts.UpEdgeColor,   opts.EdgeAlpha, bg);
colDown = blendColor(opts.DownEdgeColor, opts.EdgeAlpha, bg);

% ---- Parameters for curvature ----
Wdir         = W ~= 0;
mutual       = Wdir & Wdir.';   % mutual(i,j) = true if i->j and j->i exist
drawn        = false(N);        % track edges already drawn (i,j)
curvFracMut  = 0.15;            % curvature as fraction of distance (mutual)
curvFracSame = 0.12;            % curvature for same-level edges
tolLen       = 1e-12;
tolSame      = 1e-8;

% ---- Handle bidirectional edges with curves ----
for i = 1:N
    for j = i+1:N
        if mutual(i,j)
            x1 = X(i); y1 = Y(i);
            x2 = X(j); y2 = Y(j);

            dx  = x2 - x1;
            dy  = y2 - y1;
            len = hypot(dx, dy);
            if len < tolLen
                % Degenerate case: nodes at (almost) same position.
                % For now, fall back to straight edges later.
                continue;
            end

            % Unit tangent and normal
            tx = dx / len;
            ty = dy / len;
            nx = -ty;
            ny =  tx;

            c = curvFracMut * len;
            s = linspace(0, 1, 25);

            % Colours for i->j and j->i directions
            if h(j) >= h(i)
                col_ij = colUp;
            else
                col_ij = colDown;
            end
            if h(i) >= h(j)
                col_ji = colUp;
            else
                col_ji = colDown;
            end

            % Curve for i -> j (offset +n)
            xCurve_ij = x1 + dx.*s + c * sin(pi*s) .* nx;
            yCurve_ij = y1 + dy.*s + c * sin(pi*s) .* ny;
            drawCurvedArrow(xCurve_ij, yCurve_ij, col_ij, ...
                opts.EdgeWidth, arrowLen, opts.ArrowAngle);

            % Curve for j -> i (offset -n)
            dx2 = x1 - x2;
            dy2 = y1 - y2;
            xCurve_ji = x2 + dx2.*s - c * sin(pi*s) .* nx;
            yCurve_ji = y2 + dy2.*s - c * sin(pi*s) .* ny;
            drawCurvedArrow(xCurve_ji, yCurve_ji, col_ji, ...
                opts.EdgeWidth, arrowLen, opts.ArrowAngle);

            drawn(i,j) = true;
            drawn(j,i) = true;
        end
    end
end

% ---- Plot remaining edges (non-bidirectional) ----
for e = 1:E
    i = iList(e); j = jList(e);
    if drawn(i,j)
        continue;   % already drawn as part of a mutual pair
    end

    x1 = X(i); y1 = Y(i);
    x2 = X(j); y2 = Y(j);

    % Choose colour based on direction
    if isUp(e)
        col = colUp;
    else
        col = colDown;
    end

    % Same-level edges: draw as curved, not straight
    if abs(h(i) - h(j)) < tolSame
        dx  = x2 - x1;
        dy  = y2 - y1;
        len = hypot(dx, dy);
        if len < tolLen
            % Degenerate: fall back to straight line & arrow
            plot([x1 x2], [y1 y2], '-', ...
                'Color', col, ...
                'LineWidth', opts.EdgeWidth);
            drawArrowhead(x1, y1, x2, y2, col, arrowLen, opts.ArrowAngle);
            continue;
        end

        % Tangent and normal
        tx = dx / len; %#ok<NASGU>
        ty = dy / len;
        nx = -ty;
        ny =  tx;

        c = curvFracSame * len;
        s = linspace(0, 1, 25);

        xCurve = x1 + dx.*s + c * sin(pi*s) .* nx;
        yCurve = y1 + dy.*s + c * sin(pi*s) .* ny;

        drawCurvedArrow(xCurve, yCurve, col, ...
            opts.EdgeWidth, arrowLen, opts.ArrowAngle);
        continue;
    end

    % Default: straight edge with arrow
    plot([x1 x2], [y1 y2], '-', ...
        'Color', col, ...
        'LineWidth', opts.EdgeWidth);

    drawArrowhead(x1, y1, x2, y2, col, arrowLen, opts.ArrowAngle);
end

% ---- Plot nodes on top ----
scatter(X, Y, opts.NodeSize, ...
    'MarkerFaceColor', opts.NodeColor, ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 0.5);

% ---- Plot node labels (optional) ----
if opts.ShowLabels
    % Determine labels: user-specified or default to node indices
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

    for ii = 1:N
        text(X(ii) + dx, Y(ii) + dy, lbls{ii}, ...
            'FontSize',  opts.LabelFontSize, ...
            'Color',     opts.LabelColor, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment',   'bottom');
    end
end

% ---- Axes scaling and limits ----
axis equal;

% Base padding from node spread
rangeX   = rangeNonzero(X);
rangeY   = rangeNonzero(Y);
basePadX = 0.05 * rangeX;
basePadY = 0.05 * rangeY;

% Extra padding to account for curved edges
curvMax = max(curvFracMut, curvFracSame);
extraX  = curvMax * rangeY;   % vertical span → horizontal bulge
extraY  = curvMax * rangeX;   % horizontal span → vertical bulge

padX = basePadX + extraX;
padY = basePadY + extraY;

% Initial limits based on node positions + padding
xLim0 = [min(X) - padX, max(X) + padX];
yLim0 = [min(Y) - padY, max(Y) + padY];

spanX = diff(xLim0);
spanY = diff(yLim0);

if spanY < 1e-12
    spanY = 1e-12;
end

aspect       = spanX / spanY;
aspectThresh = 2.5;

if aspect > aspectThresh || aspect < 1/aspectThresh
    % Extreme aspect ratio → square the window
    spanMax = max(spanX, spanY);
    xCenter = mean(xLim0);
    yCenter = mean(yLim0);
    xLim    = [xCenter - spanMax/2, xCenter + spanMax/2];
    yLim    = [yCenter - spanMax/2, yCenter + spanMax/2];
else
    % Mild aspect ratio → use padded limits as-is
    xLim = xLim0;
    yLim = yLim0;
end

xlim(xLim);
ylim(yLim);

% Remove ticks and axis labels for a clean network diagram
set(gca, 'XTick', [], 'YTick', []);
set(gca, 'XTickLabel', [], 'YTickLabel', []);
box on;

% ---- Add title if requested ----
if ~isempty(opts.Title)
    title(opts.Title, 'Interpreter', 'none');
end

% ---- Add BFF label if requested ----
if opts.ShowBFF && ~isnan(opts.BFF)
    xl = xlim;
    yl = ylim;
    xText = xl(1) + 0.02 * (xl(2) - xl(1));
    yText = yl(2) - 0.05 * (yl(2) - yl(1));
    text(xText, yText, sprintf('BFF = %.3f', opts.BFF), ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment',   'top', ...
        'FontSize',            10, ...
        'BackgroundColor',     'w', ...
        'Margin',              2);
end

if ~holdState
    hold off;
end

end % main function


% =====================================================================
function opts = parsePlotOpts(varargin)
%PARSEPLOTOPTS  Parse name/value options for plotTFL.

% Default options
opts.ShowBands     = true;
opts.BandColor     = [0.75 0.75 0.75];
opts.BandStyle     = ':';
opts.BandWidth     = 0.75;

opts.NodeSize      = 50;
opts.NodeColor     = [0.1 0.3 0.7];   % muted blue

opts.UpEdgeColor   = [0 0 0];
opts.DownEdgeColor = [0.6 0.6 0.6];
opts.EdgeWidth     = 0.75;
opts.EdgeAlpha     = 0.7;             % 0=white, 1=full colour (blend)

opts.ArrowSize     = 0.06;            % relative to min(X-span, Y-span)
opts.ArrowAngle    = pi/7;            % ~26 degrees

opts.ShowLabels    = false;
opts.Labels        = {};
opts.LabelFontSize = 9;
opts.LabelColor    = [0 0 0];
opts.LabelOffset   = [0, 0.02];       % (dx, dy) in data units

opts.Title         = '';
opts.ShowBFF       = false;
opts.BFF           = NaN;

if isempty(varargin), return; end
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

        case 'arrowsize'
            opts.ArrowSize = value;
        case 'arrowangle'
            opts.ArrowAngle = value;

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

        otherwise
            error('plotTFL:UnknownOption', 'Unknown option: %s', name);
    end
end

end


% =====================================================================
function xlims = getBandXLimits(X)
%GETBANDXLIMITS  Reasonable X-range for band lines.
pad   = 0.05 * rangeNonzero(X);
xlims = [min(X)-pad, max(X)+pad];
end

function r = rangeNonzero(x)
%RANGENONZERO  Like range(x), but avoid zero for constant vectors.
r = max(x) - min(x);
if r == 0
    r = 1;
end
end

function cBlend = blendColor(c, alpha, bg)
%BLENDCOLOR  Blend colour c with background bg using "alpha" in [0,1].
% alpha=1 -> c; alpha=0 -> bg.
c  = c(:)'; 
bg = bg(:)';
cBlend = alpha * c + (1-alpha) * bg;
end

function drawArrowhead(x1, y1, x2, y2, col, L, halfAngle)
%DRAWARROWHEAD  Draw a filled triangular arrowhead at (x2,y2).

dx  = x2 - x1;
dy  = y2 - y1;
len = hypot(dx, dy);
if len == 0
    return;
end

ux = dx / len;
uy = dy / len;

% Point where arrow base starts (a bit back from the tip)
px = x2 - L * ux;
py = y2 - L * uy;

% Perpendicular unit vectors
px_perp = -uy;
py_perp =  ux;

% Half-width of the base
w = tan(halfAngle) * L;

% Two base corners
bx1 = px + w * px_perp;
by1 = py + w * py_perp;
bx2 = px - w * px_perp;
by2 = py - w * py_perp;

patch([x2 bx1 bx2], [y2 by1 by2], col, ...
    'EdgeColor', 'none', ...
    'FaceColor', col);
end

function drawCurvedArrow(xCurve, yCurve, col, edgeWidth, arrowSize, arrowAngle)
%DRAWCURVEDARROW  Draw a curved edge with an arrow at the end.

% Draw the curve
plot(xCurve, yCurve, '-', 'Color', col, 'LineWidth', edgeWidth);

% Arrowhead at the end using last segment as tangent
n = numel(xCurve);
if n < 2
    return;
end
x1 = xCurve(end-1); y1 = yCurve(end-1);
x2 = xCurve(end);   y2 = yCurve(end);

drawArrowhead(x1, y1, x2, y2, col, arrowSize, arrowAngle);
end
