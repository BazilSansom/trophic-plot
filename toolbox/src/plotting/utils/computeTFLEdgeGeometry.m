function geom = computeTFLEdgeGeometry(W, scene, varargin)
%COMPUTETFLEDGEGEOMETRY  Build shared edge geometry for TFL-family renderers.
%
%   geom = computeTFLEdgeGeometry(W, scene)
%   geom = computeTFLEdgeGeometry(W, scene, optsStruct)
%   geom = computeTFLEdgeGeometry(W, scene, 'Name', Value, ...)
%
% Purpose
% -------
% Extract edge-centerline decisions from plotTFL into a reusable geometry
% object. This includes:
%   - mutual-pair curvature
%   - same-level curvature
%   - hit-node rerouting
%   - lane splitting
%   - arrow anchor segments near the destination
%
% This helper does NOT draw anything.
%
% Expected scene fields
% ---------------------
% scene.Xplot
% scene.Yplot
% scene.lev
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
% geom.edges(e).arrowX1
% geom.edges(e).arrowY1
% geom.edges(e).arrowX2
% geom.edges(e).arrowY2
% geom.edges(e).toNode
% geom.edges(e).isMutualTopological
%
% plus convenience mirrors:
% geom.i, geom.j, geom.lin, geom.w, geom.x, geom.y, geom.kind, geom.isUp,
% geom.len, geom.arrowX1, geom.arrowY1, geom.arrowX2, geom.arrowY2, geom.toNode
% geom.allX, geom.allY, geom.lev
%
% Notes
% -----
% - Geometry uses scene.Xplot / scene.Yplot.
% - Up/down classification uses scene.lev, matching current plotTFL.
% - The helper preserves current plotTFL logic as closely as possible.
%
% -------------------------------------------------------------------------

    % ---- checks ----
    if ~ismatrix(W) || size(W,1) ~= size(W,2)
        error('computeTFLEdgeGeometry:WNotSquare', ...
            'W must be a square adjacency matrix.');
    end

    N = size(W,1);

    req = {'Xplot','Yplot','lev'};
    for r = 1:numel(req)
        if ~isfield(scene, req{r})
            error('computeTFLEdgeGeometry:BadScene', ...
                'scene must contain field %s.', req{r});
        end
    end

    X = scene.Xplot(:);
    Y = scene.Yplot(:);
    lev = scene.lev(:);

    if numel(X) ~= N || numel(Y) ~= N || numel(lev) ~= N
        error('computeTFLEdgeGeometry:BadSceneSize', ...
            'scene.Xplot, scene.Yplot, and scene.lev must have length size(W,1).');
    end

    opts = parseGeomOpts_(varargin{:});
    opts.h = lev;
    opts.W = W;

    % ---- edge extraction ----
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
        'len', 0, ...
        'arrowX1', NaN, ...
        'arrowY1', NaN, ...
        'arrowX2', NaN, ...
        'arrowY2', NaN, ...
        'toNode', [], ...
        'isMutualTopological', false);

    if E > 0
        edges = repmat(edgeTemplate, E, 1);
    else
        edges = repmat(edgeTemplate, 0, 1);
    end

    for e = 1:E
        edges(e).i    = iList(e);
        edges(e).j    = jList(e);
        edges(e).lin  = sub2ind([N N], iList(e), jList(e));
        edges(e).w    = full(W(iList(e), jList(e)));
        edges(e).toNode = jList(e);
    end

    geom = struct();
    geom.edges = edges;
    geom.allX = X;
    geom.allY = Y;
    geom.lev  = lev;

    geom.iList = iList(:);
    geom.jList = jList(:);

    if E == 0
        geom.eIdx = sparse(N,N);
        geom.eRev = zeros(0,1);
        geom.offS = zeros(0,1);
        geom.offT = zeros(0,1);
        geom = addCompatibilityMirrors_(geom);
        return;
    end

    % ---- topology helpers ----
    eIdx = sparse(iList, jList, 1:E, N, N);
    eRev = full(eIdx(jList, iList));
    isMutualEdge = (eRev > 0);

    for e = 1:E
        edges(e).isMutualTopological = isMutualEdge(e);
    end

    geom.eIdx = eIdx;
    geom.eRev = eRev;

    % ---- geometric scales ----
    tolLen  = 1e-12;
    spanMax = max(rangeNonzero_(X), rangeNonzero_(Y));

    edgeLens = hypot(X(jList)-X(iList), Y(jList)-Y(iList));
    Ltyp = localRobustMedianPositive_(edgeLens, tolLen);
    if ~isfinite(Ltyp) || Ltyp < tolLen
        Ltyp = spanMax;
    end
    if ~isfinite(spanMax) || spanMax < tolLen
        spanMax = 1;
    end

    hitR      = opts.HitNodeRadiusFrac * spanMax + opts.HitNodeMarginFrac * spanMax;
    hitMinLen = opts.HitNodeMinLenFrac * Ltyp;

    % ---- lane splitting preprocessing ----
    offS = zeros(E,1);
    offT = zeros(E,1);

    if opts.LaneSplit
        isSameEdge = abs(lev(iList) - lev(jList)) < opts.TolSame;
        isExcluded = isMutualEdge | isSameEdge;

        P = [X Y];

        medL = median(edgeLens(edgeLens > tolLen));
        if isempty(medL) || ~isfinite(medL)
            medL = median(edgeLens);
        end
        if isempty(medL) || ~isfinite(medL) || medL < tolLen
            medL = 1;
        end

        minLenEff = opts.LaneSplitMinLen;
        if minLenEff <= 0
            minLenEff = 1e-6 * medL;
        else
            minLenEff = max(minLenEff, 1e-6 * medL);
        end

        [offS, offT] = computeLaneOffsets_(P, iList, jList, isExcluded, ...
            'Theta', opts.LaneSplitTheta, ...
            'LaneSpacing', opts.LaneSplitSpacing, ...
            'MinLen', minLenEff, ...
            'h', lev, ...
            'W', W, ...
            'PrimaryStraight', opts.LaneSplitPrimaryStraight, ...
            'PrimaryRule', opts.LaneSplitPrimaryRule, ...
            'OneSided', opts.LaneSplitOneSided, ...
            'ChainStraight', opts.LaneSplitChainStraight);

        if opts.LaneSplitCoherent
            [offS, offT] = enforceCoherentOffsets_(offS, offT, opts.LaneSplitCoherentTol);
        end

        if opts.LaneSplitPreferDownSide
            isDown = (lev(jList) < lev(iList));
            sgn = -1;
            offS(isDown) = sgn * abs(offS(isDown));
            offT(isDown) = sgn * abs(offT(isDown));
        end
    end

    geom.offS = offS;
    geom.offT = offT;

    % ---- classify up/down ----
    isUp = (lev(jList) >= lev(iList));
    for e = 1:E
        edges(e).isUp = isUp(e);
    end

    % ---- build paths ----
    drawnE = false(E,1);

    Wdir   = (W ~= 0);
    mutual = Wdir & Wdir.';

    % ---- mutual pairs ----
    for i = 1:N
        for j = i+1:N
            if ~mutual(i,j)
                continue;
            end

            e_ij = full(eIdx(i,j));
            e_ji = full(eIdx(j,i));

            if e_ij == 0 || e_ji == 0
                continue;
            end
            if drawnE(e_ij) || drawnE(e_ji)
                continue;
            end

            if (lev(j) > lev(i)) || (abs(lev(j)-lev(i)) < opts.TolSame && Y(j) >= Y(i))
                low  = i;
                high = j;
            else
                low  = j;
                high = i;
            end

            e_up   = full(eIdx(low,  high));
            e_down = full(eIdx(high, low));

            x1 = X(low);  y1 = Y(low);
            x2 = X(high); y2 = Y(high);

            dx  = x2 - x1;
            dy  = y2 - y1;
            len = hypot(dx, dy);

            if len < tolLen
                edges(e_up)   = assignStraightDegenerate_(edges(e_up),   x1, y1, x2, y2, 'mutual');
                edges(e_down) = assignStraightDegenerate_(edges(e_down), x2, y2, x1, y1, 'mutual');
                drawnE([e_up, e_down]) = true;
                continue;
            end

            nx = -dy / len;
            ny =  dx / len;

            if ny <= 0
                nx = -nx;
                ny = -ny;
            end

            c = opts.CurvFracMutual * len;
            s = linspace(0,1,opts.CurveNpts);

            % up edge: low -> high on +n
            xCurveU = x1 + dx.*s + c*sin(pi*s).*nx;
            yCurveU = y1 + dy.*s + c*sin(pi*s).*ny;

            % down edge: high -> low on -n
            dx2 = x1 - x2;
            dy2 = y1 - y2;
            xCurveD = x2 + dx2.*s - c*sin(pi*s).*nx;
            yCurveD = y2 + dy2.*s - c*sin(pi*s).*ny;

            edges(e_up)   = assignCurvedEdge_(edges(e_up),   xCurveU, yCurveU, len, 'mutual');
            edges(e_down) = assignCurvedEdge_(edges(e_down), xCurveD, yCurveD, len, 'mutual');

            drawnE([e_up, e_down]) = true;
        end
    end

    % ---- remaining edges ----
    for e = 1:E
        if drawnE(e)
            continue;
        end

        i = iList(e);
        j = jList(e);

        x1 = X(i); y1 = Y(i);
        x2 = X(j); y2 = Y(j);

        dx  = x2 - x1;
        dy  = y2 - y1;
        len = hypot(dx, dy);

        edges(e).len = len;

        % same-level edges
        if abs(lev(i) - lev(j)) < opts.TolSame
            if len < tolLen
                edges(e) = assignStraightDegenerate_(edges(e), x1, y1, x2, y2, 'same');
                continue;
            end

            ty = dy / len;
            nx = -ty;
            ny =  dx / len;

            c = opts.CurvFracSame * len;
            s = linspace(0,1,opts.CurveNpts);

            xCurve = x1 + dx.*s + c*sin(pi*s).*nx;
            yCurve = y1 + dy.*s + c*sin(pi*s).*ny;

            edges(e) = assignCurvedEdge_(edges(e), xCurve, yCurve, len, 'same');
            continue;
        end

        % hit-node rerouting
        if opts.HitNode && len >= max(hitMinLen, tolLen)
            [hits, kHit, pClosest] = edgeHitsAnyNode_(i, j, [x1 y1], [x2 y2], X, Y, hitR);

            if hits
                tx = dx / len;
                ty = dy / len;
                nx = -ty;
                ny =  tx;

                vkx = X(kHit) - pClosest(1);
                vky = Y(kHit) - pClosest(2);

                sgn = sign(vkx*nx + vky*ny);
                if sgn == 0
                    sgn = 1;
                end

                c = opts.HitNodeCurvFrac * len;
                s = linspace(0,1,opts.HitNodeNpts);

                xCurve = x1 + dx.*s - sgn * c*sin(pi*s).*nx;
                yCurve = y1 + dy.*s - sgn * c*sin(pi*s).*ny;

                edges(e) = assignCurvedEdge_(edges(e), xCurve, yCurve, len, 'hitnode');
                continue;
            end
        end

        % lane splitting
        if opts.LaneSplit && (offS(e) ~= 0 || offT(e) ~= 0)
            if len < tolLen
                edges(e) = assignStraightDegenerate_(edges(e), x1, y1, x2, y2, 'lanesplit');
                continue;
            end

            tx = dx / len;
            ty = dy / len;
            nx = -ty;
            ny =  tx;

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

            [xCurve, yCurve] = cubicBezier_(p0, p1, p2, p3, opts.LaneSplitNpts);
            edges(e) = assignCurvedEdge_(edges(e), xCurve, yCurve, len, 'lanesplit');
            continue;
        end

        % default straight edge
        if len < tolLen
            edges(e) = assignStraightDegenerate_(edges(e), x1, y1, x2, y2, 'straight');
        else
            edges(e) = assignStraightEdge_(edges(e), x1, y1, x2, y2, len, 'straight');
        end
    end

    geom.edges = edges;
    geom = addCompatibilityMirrors_(geom);
end

% =====================================================================
function opts = parseGeomOpts_(varargin)

    opts = defaultGeomOpts_();

    if isempty(varargin)
        return;
    end

    if numel(varargin) == 1 && isstruct(varargin{1})
        opts = mergeKnownFields_(opts, varargin{1});
        return;
    end

    if mod(numel(varargin),2) ~= 0
        error('computeTFLEdgeGeometry:BadArgs', ...
            'Optional arguments must be name/value pairs.');
    end

    for k = 1:2:numel(varargin)
        name = varargin{k};
        value = varargin{k+1};

        if ~ischar(name) && ~isstring(name)
            error('computeTFLEdgeGeometry:BadParamName', ...
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
            case 'lanesplitchainstraight'
                opts.LaneSplitChainStraight = logical(value);

            otherwise
                % ignore unknowns so a fuller plot opts struct can be passed
        end
    end
end

% =====================================================================
function opts = defaultGeomOpts_()

    opts = struct();

    opts.CurvFracMutual = 0.15;
    opts.CurvFracSame   = 0.12;
    opts.CurveNpts      = 25;
    opts.TolSame        = 1e-8;

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

    opts.h = [];
    opts.W = [];
end

% =====================================================================
function edge = assignStraightEdge_(edge, x1, y1, x2, y2, len, kind)
    edge.x = [x1 x2];
    edge.y = [y1 y2];
    edge.kind = kind;
    edge.len = len;
    edge.arrowX1 = x1;
    edge.arrowY1 = y1;
    edge.arrowX2 = x2;
    edge.arrowY2 = y2;
end

% =====================================================================
function edge = assignStraightDegenerate_(edge, x1, y1, x2, y2, kind)
    edge.x = [x1 x2];
    edge.y = [y1 y2];
    edge.kind = kind;
    edge.len = hypot(x2-x1, y2-y1);
    edge.arrowX1 = NaN;
    edge.arrowY1 = NaN;
    edge.arrowX2 = NaN;
    edge.arrowY2 = NaN;
end

% =====================================================================
function edge = assignCurvedEdge_(edge, xCurve, yCurve, len, kind)
    [xa1, ya1, xa2, ya2] = curveTerminalSegment_(xCurve, yCurve);

    edge.x = xCurve(:).';
    edge.y = yCurve(:).';
    edge.kind = kind;
    edge.len = len;
    edge.arrowX1 = xa1;
    edge.arrowY1 = ya1;
    edge.arrowX2 = xa2;
    edge.arrowY2 = ya2;
end

% =====================================================================
function [x1, y1, x2, y2] = curveTerminalSegment_(xCurve, yCurve)
    n = numel(xCurve);
    if n < 2
        x1 = NaN; y1 = NaN; x2 = NaN; y2 = NaN;
        return;
    end

    kBack = min(5, n-1);
    x1 = xCurve(end-kBack);
    y1 = yCurve(end-kBack);
    x2 = xCurve(end);
    y2 = yCurve(end);
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
        geom.arrowX1 = zeros(0,1);
        geom.arrowY1 = zeros(0,1);
        geom.arrowX2 = zeros(0,1);
        geom.arrowY2 = zeros(0,1);
        geom.toNode  = zeros(0,1);
        return;
    end

    geom.i = vertcat(geom.edges.i);
    geom.j = vertcat(geom.edges.j);
    geom.lin = vertcat(geom.edges.lin);
    geom.w = vertcat(geom.edges.w);
    geom.x = {geom.edges.x}.';
    geom.y = {geom.edges.y}.';
    geom.kind = {geom.edges.kind}.';
    geom.isUp = vertcat(geom.edges.isUp);
    geom.len = vertcat(geom.edges.len);
    geom.arrowX1 = vertcat(geom.edges.arrowX1);
    geom.arrowY1 = vertcat(geom.edges.arrowY1);
    geom.arrowX2 = vertcat(geom.edges.arrowX2);
    geom.arrowY2 = vertcat(geom.edges.arrowY2);
    geom.toNode  = vertcat(geom.edges.toNode);
end

% =====================================================================
function opts = mergeKnownFields_(opts, s)

    if isempty(s) || ~isstruct(s), return; end
    fn = fieldnames(s);
    for i = 1:numel(fn)
        if isfield(opts, fn{i})
            opts.(fn{i}) = s.(fn{i});
        end
    end
end

% =====================================================================
function [hits, kHit, pClosest] = edgeHitsAnyNode_(i, j, p1, p2, X, Y, hitR)

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
        [d, pc] = pointToSegmentDist_(pk, p1, p2);
        if d < bestD
            bestD = d;
            kHit = k;
            pClosest = pc;
        end
    end

    hits = (bestD < hitR);
end

% =====================================================================
function [d, pClosest] = pointToSegmentDist_(p, a, b)

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
function [offS, offT] = computeLaneOffsets_(P, ei, ej, isExcluded, varargin)

    opts.Theta           = deg2rad(6);
    opts.LaneSpacing     = 0.10;
    opts.MinLen          = 1e-6;

    opts.PrimaryStraight = true;
    opts.PrimaryRule     = 'up';
    opts.OneSided        = true;
    opts.ChainStraight   = true;

    opts.h = [];
    opts.W = [];

    opts = parseLocalOpts_(opts, varargin{:});

    E = numel(ei);
    N = size(P,1);

    offS = zeros(E,1);
    offT = zeros(E,1);

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

        a    = zeros(numel(edges),1);
        keep = true(numel(edges),1);

        for idx = 1:numel(edges)
            k = edges(idx);
            if ei(k) == v
                u = ej(k);
            else
                u = ei(k);
            end

            d = P(u,:) - P(v,:);
            len = hypot(d(1), d(2));
            if len < opts.MinLen
                keep(idx) = false;
                continue;
            end

            tv = d / len;
            ai = atan2(tv(2), tv(1));
            a(idx) = mod(ai, pi);
        end

        edges = edges(keep);
        a     = a(keep);

        if numel(edges) < 2
            continue;
        end

        [aS, ord] = sort(a);
        edgesS = edges(ord);

        groups = {};
        g = 1;
        groups{g} = 1; %#ok<AGROW>
        for ii = 2:numel(edgesS)
            if (aS(ii) - aS(ii-1)) <= theta
                groups{g}(end+1) = ii;
            else
                g = g + 1;
                groups{g} = ii; %#ok<AGROW>
            end
        end

        if numel(groups) >= 2 && (aS(1) + pi - aS(end)) <= theta
            groups{1} = [groups{end}, groups{1}];
            groups(end) = [];
        end

        for gi = 1:numel(groups)
            idxs = groups{gi};
            mG = numel(idxs);
            if mG < 2
                continue;
            end

            eG = edgesS(idxs);
            [~, perm] = sort(eG);
            eG = eG(perm);

            lanes = zeros(mG,1);

            if opts.ChainStraight && mG == 2 && ~isempty(opts.h)
                k1 = eG(1); k2 = eG(2);

                in1  = (ej(k1) == v); out1 = (ei(k1) == v);
                in2  = (ej(k2) == v); out2 = (ei(k2) == v);
                isThrough = (in1 && out2) || (out1 && in2);

                if isThrough
                    up1 = (opts.h(ej(k1)) >= opts.h(ei(k1)));
                    up2 = (opts.h(ej(k2)) >= opts.h(ei(k2)));

                    if up1 && up2
                        lanes(:) = 0;
                        writeOffsetsForGroup_(v, eG, lanes);
                        continue;
                    end
                end
            end

            if opts.PrimaryStraight && mG >= 2
                p0 = choosePrimaryEdge_(v, eG, opts, P, ei, ej);

                if opts.OneSided
                    q = 1;
                    for p = 1:mG
                        if p == p0, continue; end
                        lanes(p) = q * opts.LaneSpacing;
                        q = q + 1;
                    end
                else
                    others = setdiff(1:mG, p0, 'stable');
                    mO = numel(others);
                    base = (-(mO-1)/2 : (mO-1)/2) * opts.LaneSpacing;
                    lanes(others) = base(:);
                    lanes(p0) = 0;
                end
            else
                lanes = (-(mG-1)/2 : (mG-1)/2) * opts.LaneSpacing;
                lanes = lanes(:);
            end

            writeOffsetsForGroup_(v, eG, lanes);
        end
    end

    function writeOffsetsForGroup_(vv, eGrp, laneVals)
        for pp = 1:numel(eGrp)
            kk = eGrp(pp);
            if ei(kk) == vv
                offS(kk) = offS(kk) + laneVals(pp);
            else
                offT(kk) = offT(kk) + laneVals(pp);
            end
        end
    end
end

% =====================================================================
function p0 = choosePrimaryEdge_(~, eGrp, opts, P, ei, ej)

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
                    score = -inf(numel(cand),1);
                    for ii = 1:numel(cand)
                        pp = cand(ii);
                        kk = eGrp(pp);
                        d = P(ej(kk),:) - P(ei(kk),:);
                        L = hypot(d(1), d(2));
                        if L > 0
                            score(ii) = abs(d(2)) / L;
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
                    end
                end

                cand = find(isUp);
                if ~isempty(cand)
                    [~, jj] = max(wgt(cand));
                    p0 = cand(jj);
                    return;
                end
            end

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

% =====================================================================
function opts = parseLocalOpts_(opts, varargin)
    for i = 1:2:numel(varargin)
        opts.(varargin{i}) = varargin{i+1};
    end
end

% =====================================================================
function [xCurve, yCurve] = cubicBezier_(p0, p1, p2, p3, npts)

    t = linspace(0,1,npts).';
    b0 = (1-t).^3;
    b1 = 3*(1-t).^2.*t;
    b2 = 3*(1-t).*t.^2;
    b3 = t.^3;

    B = b0.*p0 + b1.*p1 + b2.*p2 + b3.*p3;
    xCurve = B(:,1).';
    yCurve = B(:,2).';
end

% =====================================================================
function [offS2, offT2] = enforceCoherentOffsets_(offS, offT, tol)

    offS2 = offS;
    offT2 = offT;
    m = numel(offS2);

    for k = 1:m
        a = offS2(k);
        b = offT2(k);

        if a == 0 || b == 0
            continue;
        end

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
function r = rangeNonzero_(x)

    r = max(x) - min(x);
    if r == 0
        r = 1;
    end
end