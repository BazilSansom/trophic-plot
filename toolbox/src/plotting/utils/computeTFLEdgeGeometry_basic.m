function geom = computeTFLEdgeGeometry_basic(W, X, Y, h, varargin)
%COMPUTETFLEDGEGEOMETRY_BASIC  Basic shared TFL edge geometry.
%
%   geom = computeTFLEdgeGeometry_basic(W, X, Y, h)
%   geom = computeTFLEdgeGeometry_basic(W, X, Y, h, 'Name', Value, ...)
%
% Minimal geometry builder for coincident overlay development.
%
% Behaviour
% ---------
% - straight edges by default
% - mutual pairs are put onto opposite curves
% - same-level unidirectional edges are curved
% - no hit-node rerouting
% - no lane splitting
%
% Inputs
% ------
% W : adjacency matrix, rows = source, cols = destination
% X : rendered x coordinates
% Y : rendered y coordinates
% h : trophic levels (used for up/down semantics)
%
% Options
% -------
% 'CurvFracMutual' : default 0.15
% 'CurvFracSame'   : default 0.12
% 'CurveNpts'      : default 25
% 'TolSame'        : default 1e-8
%
% Output
% ------
% geom.edges(e).i
% geom.edges(e).j
% geom.edges(e).lin
% geom.edges(e).w
% geom.edges(e).x
% geom.edges(e).y
% geom.edges(e).kind
% geom.edges(e).isUp
% geom.edges(e).len
%
% geom.allX, geom.allY
%
% Compatibility mirrors:
% geom.i, geom.j, geom.lin, geom.w, geom.x, geom.y, geom.kind, geom.isUp, geom.len
%
% -------------------------------------------------------------------------

    if ~ismatrix(W) || size(W,1) ~= size(W,2)
        error('computeTFLEdgeGeometry_basic:WNotSquare', 'W must be square.');
    end

    n = size(W,1);

    X = X(:);
    Y = Y(:);
    h = h(:);

    if numel(X) ~= n || numel(Y) ~= n || numel(h) ~= n
        error('computeTFLEdgeGeometry_basic:BadXYH', ...
            'X, Y, h must all have length size(W,1).');
    end

    opts = parseOpts_(varargin{:});

    [iList, jList] = find(W ~= 0);
    E = numel(iList);

    edgeTemplate = struct( ...
        'i', [], ...
        'j', [], ...
        'lin', [], ...
        'w', [], ...
        'x', [], ...
        'y', [], ...
        'kind', '', ...
        'isUp', false, ...
        'len', 0);

    if E > 0
        edges = repmat(edgeTemplate, E, 1);
    else
        edges = repmat(edgeTemplate, 0, 1);
    end

    for e = 1:E
        i = iList(e);
        j = jList(e);
        edges(e).i   = i;
        edges(e).j   = j;
        edges(e).lin = sub2ind([n n], i, j);
        edges(e).w   = full(W(i,j));
    end

    geom = struct();
    geom.edges = edges;
    geom.allX  = X(:);
    geom.allY  = Y(:);

    if E == 0
        geom = addCompatibilityMirrors_(geom);
        return;
    end

    eIdx = sparse(iList, jList, 1:E, n, n);
    drawn = false(E,1);

    Wdir = (W ~= 0);
    mutual = Wdir & Wdir.';

    % -------------------- mutual pairs --------------------
    for i = 1:n
        for j = i+1:n
            if ~mutual(i,j)
                continue;
            end

            eij = full(eIdx(i,j));
            eji = full(eIdx(j,i));

            if eij == 0 || eji == 0
                continue;
            end
            if drawn(eij) || drawn(eji)
                continue;
            end

            if (h(j) > h(i)) || (abs(h(j)-h(i)) < opts.TolSame && Y(j) >= Y(i))
                low  = i;
                high = j;
            else
                low  = j;
                high = i;
            end

            eUp   = full(eIdx(low,  high));
            eDown = full(eIdx(high, low));

            x1 = X(low);  y1 = Y(low);
            x2 = X(high); y2 = Y(high);

            dx = x2 - x1;
            dy = y2 - y1;
            len = hypot(dx, dy);

            if len <= 0
                edges(eUp).x = [x1 x2];
                edges(eUp).y = [y1 y2];
                edges(eUp).kind = 'straight';
                edges(eUp).isUp = true;
                edges(eUp).len = 0;

                edges(eDown).x = [x2 x1];
                edges(eDown).y = [y2 y1];
                edges(eDown).kind = 'straight';
                edges(eDown).isUp = false;
                edges(eDown).len = 0;

                drawn([eUp, eDown]) = true;
                continue;
            end

            nx = -dy / len;
            ny =  dx / len;

            % orient the pair so the "down" arc is visually the lower bulge
            if ny <= 0
                nx = -nx;
                ny = -ny;
            end

            c = opts.CurvFracMutual * len;
            s = linspace(0,1,opts.CurveNpts);

            % low -> high
            xUp = x1 + dx.*s + c*sin(pi*s).*nx;
            yUp = y1 + dy.*s + c*sin(pi*s).*ny;

            % high -> low on opposite side
            dx2 = x1 - x2;
            dy2 = y1 - y2;
            xDn = x2 + dx2.*s - c*sin(pi*s).*nx;
            yDn = y2 + dy2.*s - c*sin(pi*s).*ny;

            edges(eUp).x = xUp;
            edges(eUp).y = yUp;
            edges(eUp).kind = 'mutual';
            edges(eUp).isUp = true;
            edges(eUp).len = len;

            edges(eDown).x = xDn;
            edges(eDown).y = yDn;
            edges(eDown).kind = 'mutual';
            edges(eDown).isUp = false;
            edges(eDown).len = len;

            drawn([eUp, eDown]) = true;
        end
    end

    % -------------------- remaining edges --------------------
    for e = 1:E
        if drawn(e)
            continue;
        end

        i = edges(e).i;
        j = edges(e).j;

        x1 = X(i); y1 = Y(i);
        x2 = X(j); y2 = Y(j);

        dx = x2 - x1;
        dy = y2 - y1;
        len = hypot(dx, dy);

        edges(e).isUp = (h(j) >= h(i));
        edges(e).len = len;

        if len <= 0
            edges(e).x = [x1 x2];
            edges(e).y = [y1 y2];
            edges(e).kind = 'straight';
            continue;
        end

        if abs(h(i) - h(j)) < opts.TolSame
            tx = dx / len;
            ty = dy / len;
            nx = -ty;
            ny =  tx;

            c = opts.CurvFracSame * len;
            s = linspace(0,1,opts.CurveNpts);

            xCurve = x1 + dx.*s + c*sin(pi*s).*nx;
            yCurve = y1 + dy.*s + c*sin(pi*s).*ny;

            edges(e).x = xCurve;
            edges(e).y = yCurve;
            edges(e).kind = 'same';
        else
            edges(e).x = [x1 x2];
            edges(e).y = [y1 y2];
            edges(e).kind = 'straight';
        end
    end

    geom.edges = edges;
    geom = addCompatibilityMirrors_(geom);
end

% =====================================================================
function geom = addCompatibilityMirrors_(geom)

    E = numel(geom.edges);

    if E == 0
        geom.i = zeros(0,1);
        geom.j = zeros(0,1);
        geom.lin = zeros(0,1);
        geom.w = zeros(0,1);
        geom.x = cell(0,1);
        geom.y = cell(0,1);
        geom.kind = cell(0,1);
        geom.isUp = false(0,1);
        geom.len = zeros(0,1);
        return;
    end

    geom.i    = vertcat(geom.edges.i);
    geom.j    = vertcat(geom.edges.j);
    geom.lin  = vertcat(geom.edges.lin);
    geom.w    = vertcat(geom.edges.w);
    geom.x    = {geom.edges.x}.';
    geom.y    = {geom.edges.y}.';
    geom.kind = {geom.edges.kind}.';
    geom.isUp = vertcat(geom.edges.isUp);
    geom.len  = vertcat(geom.edges.len);
end

% =====================================================================
function opts = parseOpts_(varargin)

    opts = struct();
    opts.CurvFracMutual = 0.15;
    opts.CurvFracSame   = 0.12;
    opts.CurveNpts      = 25;
    opts.TolSame        = 1e-8;

    if isempty(varargin), return; end
    if mod(numel(varargin),2) ~= 0
        error('computeTFLEdgeGeometry_basic:BadArgs', ...
            'Optional arguments must be name/value pairs.');
    end

    for k = 1:2:numel(varargin)
        name = varargin{k};
        value = varargin{k+1};

        if ~ischar(name) && ~isstring(name)
            error('computeTFLEdgeGeometry_basic:BadParamName', ...
                'Parameter names must be strings.');
        end

        switch lower(char(name))
            case 'curvfracmutual'
                opts.CurvFracMutual = value;
            case 'curvfracsame'
                opts.CurvFracSame = value;
            case 'curvenpts'
                opts.CurveNpts = value;
            case 'tolsame'
                opts.TolSame = value;
            otherwise
                error('computeTFLEdgeGeometry_basic:UnknownOption', ...
                    'Unknown option: %s', char(name));
        end
    end
end