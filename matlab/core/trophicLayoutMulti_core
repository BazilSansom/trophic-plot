function [X, Y, compIdx, info] = trophicLayoutMulti_core(W, hGlobal, varargin)
%TROPHICLAYOUTMULTI_CORE  Multi-component layout wrapper given a node hierarchy.
%
%   [X, Y, compIdx, info] = TROPHICLAYOUTMULTI_CORE(W, hGlobal, ...)
%
% INPUT
%   W        : adjacency matrix (n x n), directed, weighted or unweighted.
%   hGlobal  : n x 1 vector of node "heights" or hierarchy scores
%              (e.g. trophic levels, topological layers, etc.).
%
% OPTIONAL NAME–VALUE PAIRS:
%   'LayoutFun'            : handle to underlying layout function
%                            (default: @trophicLayout). The layout function
%                            must accept at least:
%                              LayoutFun(Wc, 'hProvided', hc, 'RescaleLevels', false, ...)
%
%   'ComponentGapFrac'     : fraction of main span used for each inter-component
%                            gap (default: 0.15). This controls how much
%                            horizontal space is left between non-isolated
%                            components.
%
%   'IsolateColumnGapFrac' : gap from main layout to isolate column as a
%                            fraction of the main span (default: 0.1).
%
%   'IsolateSpreadFrac'    : total horizontal spread of isolates as a
%                            fraction of the main span (default: 0.05).
%
%   'GlobalRescaleH'       : whether to globally centre/scale hGlobal before
%                            passing to components (default: true).
%
% Any other name–value pairs are passed directly through to LayoutFun
% (e.g. 'NumIters', 'Seed', 'AttractStrength', etc. for TROPHICLAYOUT).
%
% OUTPUT
%   X, Y     : layout coordinates for all nodes (n x 1).
%
%   compIdx  : component index for each node:
%              0 for isolates (degree == 0),
%              1..nComp for non-isolated weak components.
%
%   info     : struct with diagnostics:
%              .nComp            number of non-isolated components
%              .offsetX          horizontal offsets applied to components
%              .width            approximate width of each component
%              .isoIdx           indices of isolates
%              .activeIdx        indices of non-isolates
%              .h_raw            input hierarchy hGlobal (as provided)
%              .h_global         possibly rescaled h used internally
%              .gapAbs           absolute inter-component gap in X
%              .ComponentGapFrac the gap fraction used
%              .mainSpan         span of main (non-isolate) layout in X
%              .packOrder        component indices in packing order
%              .compSizes        size (node count) of each component
%
% This function does not plot; it only produces coordinates and diagnostics.
% It is designed to be used together with TROPHICLAYOUT and PLOTTFL.
%
% Author:
%   Bazil Sansom (Warwick Business School, University of Warwick)
%   Contact: bazil.sansom@wbs.ac.uk
%
% -------------------------------------------------------------------------

    % ---------- basic checks ----------
    n = size(W,1);
    if size(W,2) ~= n
        error('trophicLayoutMulti_core:WNotSquare', ...
              'W must be a square adjacency matrix.');
    end
    if numel(hGlobal) ~= n
        error('trophicLayoutMulti_core:BadHGlobal', ...
              'hGlobal must be a vector of length n.');
    end

    if issparse(W)
        W = full(W);
    end

    h_raw = hGlobal(:);

    % ---------- defaults ----------
    opts.LayoutFun            = @trophicLayout;
    opts.ComponentGapFrac     = 0.15;
    opts.IsolateColumnGapFrac = 0.1;
    opts.IsolateSpreadFrac    = 0.05;
    opts.GlobalRescaleH       = true;

    passThrough = {};

    % ---------- parse name–value pairs ----------
    if mod(numel(varargin),2) ~= 0
        error('trophicLayoutMulti_core:BadArgs', ...
              'Optional arguments must be name/value pairs.');
    end

    for k = 1:2:numel(varargin)
        name  = varargin{k};
        value = varargin{k+1};
        if ~ischar(name) && ~isstring(name)
            error('trophicLayoutMulti_core:BadParamName', ...
                  'Parameter names must be strings.');
        end
        switch lower(char(name))
            case 'layoutfun'
                opts.LayoutFun = value;
            case 'componentgapfrac'
                opts.ComponentGapFrac = value;
            case 'isolatecolumngapfrac'
                opts.IsolateColumnGapFrac = value;
            case 'isolatespreadfrac'
                opts.IsolateSpreadFrac = value;
            case 'globalrescaleh'
                opts.GlobalRescaleH = logical(value);
            otherwise
                % Unrecognised options are passed through to LayoutFun
                passThrough{end+1} = name;  %#ok<AGROW>
                passThrough{end+1} = value; %#ok<AGROW>
        end
    end

    % ---------- identify isolates vs active nodes ----------
    deg      = sum(W + W.', 2);   % total degree (in + out) for each node
    isoMask  = (deg == 0);
    isoIdx   = find(isoMask);
    activeIdx = find(~isoMask);

    X = nan(n,1);
    Y = nan(n,1);

    % ---------- optionally rescale hierarchy globally ----------
    if opts.GlobalRescaleH
        u  = ones(n,1);
        mu = sum(u .* h_raw) / sum(u);
        h_centered = h_raw - mu;
        s2 = sum(u .* (h_centered.^2)) / sum(u);
        s  = sqrt(max(s2, eps));
        h_global = h_centered / s;
    else
        h_global = h_raw;
    end

    % ---------- edge case: all isolates ----------
    if isempty(activeIdx)
        X(:) = 0;
        Y(:) = (0:n-1)';   % arbitrary vertical line
        compIdx                  = zeros(n,1);
        info.nComp               = 0;
        info.offsetX             = [];
        info.width               = [];
        info.isoIdx              = isoIdx;
        info.activeIdx           = activeIdx;
        info.h_raw               = h_raw;
        info.h_global            = h_global;
        info.gapAbs              = 0;
        info.ComponentGapFrac    = opts.ComponentGapFrac;
        info.mainSpan            = 0;
        info.packOrder           = [];
        info.compSizes           = [];
        return;
    end

    % ---------- weakly connected components on active subgraph ----------
    Wactive   = W(activeIdx, activeIdx);
    G         = graph(Wactive + Wactive.', 'omitselfloops');
    compLocal = conncomp(G);           % labels 1..nComp on activeIdx
    nComp     = max(compLocal);

    % Component sizes (number of nodes per component)
    compSizes = accumarray(compLocal(:), 1, [nComp, 1]);
    % Order components by size: largest first (for packing)
    [~, packOrder] = sort(compSizes, 'descend');

    compIdx = zeros(n,1);
    compIdx(activeIdx) = compLocal;    % isolates remain 0

    % Storage for each component
    xLoc  = cell(nComp,1);
    yLoc  = cell(nComp,1);
    nodes = cell(nComp,1);
    width = zeros(nComp,1);

    layoutFun  = opts.LayoutFun;
    layoutArgs = passThrough;

    % ---------- layout each non-isolated component ----------
    for c = 1:nComp
        localIdx = activeIdx(compLocal == c);
        nodes{c} = localIdx;

        Wc = W(localIdx, localIdx);
        hc = h_global(localIdx);

        [xc, yc] = layoutFun( ...
                        Wc, ...
                        'hProvided',     hc, ...
                        'RescaleLevels', false, ...
                        layoutArgs{:});

        xc = xc(:);
        yc = yc(:);

        % recentre in x to mean zero so we can pack nicely
        if ~isempty(xc)
            xc = xc - mean(xc);
        end

        xLoc{c} = xc;
        yLoc{c} = yc;

        if numel(xc) > 1
            width(c) = max(xc) - min(xc);
        else
            width(c) = 1;  % default width for single-node component
        end
    end

    % ---------- compute inter-component gap (absolute) ----------
    Wsum = sum(width);
    if nComp > 1
        g = opts.ComponentGapFrac;
        % Prevent gaps from exceeding available span
        maxG = 1 / (nComp - 1) - 1e-6;
        if g > maxG
            g = maxG;
        elseif g < 0
            g = 0;
        end
        mainSpan = Wsum / (1 - (nComp - 1)*g);
        gapAbs   = g * mainSpan;
    else
        mainSpan = width(1);
        gapAbs   = 0;
    end

    % ---------- pack components horizontally (largest → smallest) ----------
    offsetX = zeros(nComp,1);
    xCursor = 0;
    first   = true;

    for k = 1:nComp
        c  = packOrder(k);   % component index in size order
        xc = xLoc{c};
        yc = yLoc{c};
        if isempty(xc)
            continue;
        end

        left  = min(xc);
        right = max(xc);

        if first
            % First (largest) component: left edge at 0
            shift = -left;
            first = false;
        else
            % Subsequent components: left edge at (previous right + gapAbs)
            shift = (xCursor + gapAbs) - left;
        end

        offsetX(c)   = shift;
        X(nodes{c})  = xc + shift;
        Y(nodes{c})  = yc;
        xCursor      = right + shift;
    end

    % ---------- place isolates in their own band/column ----------
    if ~isempty(isoIdx)
        numIso = numel(isoIdx);

        mainX = X(activeIdx);
        mainSpanX = max(mainX) - min(mainX);
        if mainSpanX <= 0
            mainSpanX = 1;
        end

        gapFrac = opts.IsolateColumnGapFrac;
        baseX   = max(mainX) + gapFrac * mainSpanX;

        spreadFrac = opts.IsolateSpreadFrac;
        if numIso > 1
            totalIsoSpan = spreadFrac * mainSpanX;
            dx = totalIsoSpan / (numIso - 1);
            offsets = ((0:numIso-1)' - (numIso-1)/2) * dx;
        else
            offsets = 0;
        end

        X(isoIdx) = baseX + offsets;
        Y(isoIdx) = h_global(isoIdx);
    end

    % ---------- info struct ----------
    info.nComp            = nComp;
    info.offsetX          = offsetX;
    info.width            = width;
    info.isoIdx           = isoIdx;
    info.activeIdx        = activeIdx;
    info.h_raw            = h_raw;
    info.h_global         = h_global;
    info.gapAbs           = gapAbs;
    info.ComponentGapFrac = opts.ComponentGapFrac;
    info.mainSpan         = max(X(activeIdx)) - min(X(activeIdx));
    info.packOrder        = packOrder;
    info.compSizes        = compSizes;

end
