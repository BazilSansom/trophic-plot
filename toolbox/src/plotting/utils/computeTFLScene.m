function scene = computeTFLScene(X, Y, h, varargin)
%COMPUTETFLSCENE  Build shared TFL scene/state from node coordinates.
%
%   scene = computeTFLScene(X, Y, h)
%   scene = computeTFLScene(X, Y, h, optsStruct)
%   scene = computeTFLScene(X, Y, h, 'Name', Value, ...)
%
% Purpose
% -------
% Separates scene construction from rendering. This computes:
%   - plotted coordinates after FillSquare / YScale / XScale
%   - tau-level semantics used for up/down and same-level logic
%   - band positions
%   - node sizes
%   - suggested axis limits / padding
%
% Notes
% -----
% - Incoming Y is assumed to be in tau-units.
% - lev is taken to be Y_tau, matching current plotTFL behaviour.
% - This helper is intended to be shared across plain TFL, CDFD overlay,
%   and colour-based CDFD rendering.
%
% Output fields
% -------------
% scene.N
% scene.Xtau, scene.Ytau, scene.h
% scene.Xplot, scene.Yplot
% scene.lev
% scene.x0, scene.y0
% scene.fillScale, scene.xScaleTotal, scene.yScaleTotal
% scene.rangeX, scene.rangeY
% scene.bandsTau, scene.bandsY, scene.bandXLimits
% scene.nodeSize
% scene.padX, scene.padY
% scene.xlimSuggested, scene.ylimSuggested
% scene.opts

    % ---- checks ----
    X = X(:);
    Y = Y(:);
    h = h(:);

    N = numel(X);
    if numel(Y) ~= N || numel(h) ~= N
        error('computeTFLScene:BadXYH', ...
            'X, Y, h must be vectors of the same length.');
    end
    if any(~isfinite(X)) || any(~isfinite(Y)) || any(~isfinite(h))
        error('computeTFLScene:NonFiniteXYH', ...
            'X, Y, and h must be finite.');
    end

    opts = parseSceneOpts_(varargin{:});

    sem = lower(strtrim(opts.LevelSemantics));
    switch sem
        case {'trophic','generic'}
            % ok
        otherwise
            error('computeTFLScene:BadLevelSemantics', ...
                'LevelSemantics must be ''trophic'' or ''generic''.');
    end
    opts.LevelSemantics = sem;

    % ---- semantic tau coordinates ----
    X_tau = X;
    Y_tau = Y;

    % Current plotTFL convention: use rendered tau-heights for all layer logic
    lev = Y_tau;

    % ---- plotted coordinates ----
    Xplot = X_tau;
    Yplot = Y_tau;

    spanX = rangeNonzero_(Xplot);
    spanY = rangeNonzero_(Yplot);

    fillScale = 1;
    if opts.FillSquare && spanY > 1e-12
        fillScale = spanX / spanY;
    end

    yScaleTotal = fillScale * opts.YScale;
    y0 = mean(Y_tau);
    if abs(yScaleTotal - 1) > 1e-12
        Yplot = y0 + (Y_tau - y0) * yScaleTotal;
    end

    xScaleTotal = opts.XScale;
    x0 = mean(X_tau);
    if abs(xScaleTotal - 1) > 1e-12
        Xplot = x0 + (X_tau - x0) * xScaleTotal;
    end

    rangeX = rangeNonzero_(Xplot);
    rangeY = rangeNonzero_(Yplot);

    
    % ---- bands ----
    tau = opts.BandTau;
    if ~(isfinite(tau) && tau > 0)
        error('computeTFLScene:BadBandTau', ...
            'BandTau must be a positive finite scalar.');
    end

    showSemanticBands = strcmpi(opts.LevelSemantics, 'trophic');

    if showSemanticBands
        k0 = floor(min(Y_tau) / tau);
        k1 = ceil(max(Y_tau) / tau);

        bandsTau = (k0:k1) * tau;
        bandsY   = y0 + (bandsTau - y0) * yScaleTotal;
    else
        bandsTau = zeros(0,1);
        bandsY   = zeros(0,1);
    end

    bandXLimits = getBandXLimits_(Xplot);

    %{

    % ---- bands ----
    tau = opts.BandTau;
    if ~(isfinite(tau) && tau > 0)
        error('computeTFLScene:BadBandTau', ...
            'BandTau must be a positive finite scalar.');
    end

    k0 = floor(min(Y_tau) / tau);
    k1 = ceil(max(Y_tau) / tau);

    bandsTau = (k0:k1) * tau;
    bandsY   = y0 + (bandsTau - y0) * yScaleTotal;
    bandXLimits = getBandXLimits_(Xplot);

    %}

    % ---- node sizes ----
    if ~isempty(opts.NodeSizeData)
        sData = double(opts.NodeSizeData(:));
        if numel(sData) ~= N
            error('computeTFLScene:BadNodeSizeData', ...
                'NodeSizeData must be empty or length N.');
        end

        if ~isempty(opts.NodeSizeScaleData)
            sScale = double(opts.NodeSizeScaleData(:));
            if any(~isfinite(sScale))
                error('computeTFLScene:BadNodeSizeScaleData', ...
                    'NodeSizeScaleData must contain only finite values.');
            end
        else
            sScale = sData;
        end

        if any(~isfinite(sData))
            error('computeTFLScene:BadNodeSizeData', ...
                'NodeSizeData must contain only finite values.');
        end

        sData  = max(sData,  0);
        sScale = max(sScale, 0);

        if all(sData == 0) || isempty(sScale) || all(sScale == 0)
            nodeSize = opts.NodeSize * ones(N,1);
        else
            switch lower(strtrim(opts.NodeSizeScaleMode))
                case 'log'
                    zData  = log1p(sData);
                    zScale = log1p(sScale);

                case 'radius'
                    zData  = sqrt(sData);
                    zScale = sqrt(sScale);

                otherwise  % 'area'
                    % Keep current visual convention: map sqrt(data) to marker area
                    zData  = sqrt(sData);
                    zScale = sqrt(sScale);
            end

            zLo = min(zScale);
            zHi = max(zScale);

            if ~isfinite(zLo) || ~isfinite(zHi) || zHi <= zLo
                nodeSize = mean(opts.NodeSizeRange) * ones(N,1);
            else
                z = (zData - zLo) ./ (zHi - zLo);
                z = min(max(z, 0), 1);

                r = opts.NodeSizeRange;
                nodeSize = r(1) + z * (r(2) - r(1));
            end
        end
    else
        nodeSize = opts.NodeSize;
        if strcmpi(opts.NodeSizeMode, 'auto')
            Nref = max(1, opts.NodeSizeRefN);
            scale = (Nref / N) ^ opts.NodeSizeGamma;
            nodeSize = nodeSize * scale;
            nodeSize = min(max(nodeSize, opts.NodeSizeMin), opts.NodeSizeMax);
        end
        nodeSize = nodeSize * ones(N,1);
    end

    %{

    % ---- node sizes ----
    if ~isempty(opts.NodeSizeData)
        s = double(opts.NodeSizeData(:));
        if numel(s) ~= N
            error('computeTFLScene:BadNodeSizeData', ...
                'NodeSizeData must be empty or length N.');
        end

        s = max(s, 0);

        if all(s == 0) || all(~isfinite(s))
            nodeSize = opts.NodeSize * ones(N,1);
        else
            switch lower(strtrim(opts.NodeSizeScaleMode))
                case 'log'
                    z = log1p(s);
                case 'radius'
                    z = sqrt(s);
                otherwise  % 'area'
                    z = sqrt(s);
            end

            z = z - min(z);
            if max(z) > 0
                z = z / max(z);
            end

            r = opts.NodeSizeRange;
            nodeSize = r(1) + z * (r(2) - r(1));
        end
    else
        nodeSize = opts.NodeSize;
        if strcmpi(opts.NodeSizeMode, 'auto')
            Nref = max(1, opts.NodeSizeRefN);
            scale = (Nref / N) ^ opts.NodeSizeGamma;
            nodeSize = nodeSize * scale;
            nodeSize = min(max(nodeSize, opts.NodeSizeMin), opts.NodeSizeMax);
        end
        nodeSize = nodeSize * ones(N,1);
    end

    %}

    % ---- suggested limits / padding ----
    curvMax = max(opts.CurvFracMutual, opts.CurvFracSame);

    basePadX = opts.PadFracX * rangeX;
    basePadY = opts.PadFracY * rangeY;

    padX = basePadX + opts.CurvPadScaleX * curvMax * rangeY;
    padY = basePadY + opts.CurvPadScaleY * curvMax * rangeX;

    xlimSuggested = [min(Xplot) - padX, max(Xplot) + padX];
    ylimSuggested = [min(Yplot) - padY, max(Yplot) + padY];

    % ---- output ----
    scene = struct();
    scene.N = N;

    scene.Xtau = X_tau;
    scene.Ytau = Y_tau;
    scene.h    = h;

    scene.Xplot = Xplot;
    scene.Yplot = Yplot;
    scene.lev   = lev;

    scene.x0 = x0;
    scene.y0 = y0;

    scene.fillScale   = fillScale;
    scene.xScaleTotal = xScaleTotal;
    scene.yScaleTotal = yScaleTotal;

    scene.rangeX = rangeX;
    scene.rangeY = rangeY;

    scene.levelSemantics = opts.LevelSemantics;
    scene.hasSemanticBands = showSemanticBands;

    scene.bandsTau   = bandsTau(:);
    scene.bandsY     = bandsY(:);
    scene.bandXLimits = bandXLimits;

    scene.nodeSize = nodeSize;

    scene.padX = padX;
    scene.padY = padY;

    scene.xlimSuggested = xlimSuggested;
    scene.ylimSuggested = ylimSuggested;

    scene.opts = opts;
end

% =====================================================================
function opts = parseSceneOpts_(varargin)

    opts = defaultSceneOpts_();

    if isempty(varargin)
        return;
    end

    if numel(varargin) == 1 && isstruct(varargin{1})
        opts = mergeKnownFields_(opts, varargin{1});
        return;
    end

    if mod(numel(varargin),2) ~= 0
        error('computeTFLScene:BadArgs', ...
            'Optional arguments must be name/value pairs.');
    end

    for k = 1:2:numel(varargin)
        name = varargin{k};
        value = varargin{k+1};

        if ~ischar(name) && ~isstring(name)
            error('computeTFLScene:BadParamName', ...
                'Parameter names must be strings.');
        end

        switch lower(char(name))
            
            case 'levelsemantics'
                opts.LevelSemantics = char(value);
            
            case 'fillsquare'
                opts.FillSquare = logical(value);
            case 'yscale'
                opts.YScale = value;
            case 'xscale'
                opts.XScale = value;
            case 'xstretch'
                opts.XScale = value;

            case 'bandtau'
                opts.BandTau = value;

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
            case 'nodesizedata'
                opts.NodeSizeData = value;
            case 'nodesizescaledata'
                opts.NodeSizeScaleData = value;
            case 'nodesizescalemode'
                opts.NodeSizeScaleMode = char(value);
            case 'nodesizerange'
                opts.NodeSizeRange = value;

            case 'padfracx'
                opts.PadFracX = value;
            case 'padfracy'
                opts.PadFracY = value;
            case 'curvpadscalex'
                opts.CurvPadScaleX = value;
            case 'curvpadscaley'
                opts.CurvPadScaleY = value;
            case 'curvfracmutual'
                opts.CurvFracMutual = value;
            case 'curvfracsame'
                opts.CurvFracSame = value;

            otherwise
                % ignore unknown options silently so a fuller plot opts struct
                % can be passed in without trouble
        end
    end
end

% =====================================================================
function opts = defaultSceneOpts_()

    opts = struct();

    opts.LevelSemantics = 'generic';   % 'trophic' | 'generic'

    opts.FillSquare = true;
    opts.YScale     = 1;
    opts.XScale     = 1;

    opts.BandTau = 1;

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

    opts.PadFracX      = 0.05;
    opts.PadFracY      = 0.04;
    opts.CurvPadScaleX = 1.00;
    opts.CurvPadScaleY = 0.25;

    % Match current plotTFL defaults for padding estimation
    opts.CurvFracMutual = 0.15;
    opts.CurvFracSame   = 0.12;
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
function xlims = getBandXLimits_(X)
    pad = 0.05 * rangeNonzero_(X);
    xlims = [min(X)-pad, max(X)+pad];
end

% =====================================================================
function r = rangeNonzero_(x)
    r = max(x) - min(x);
    if r == 0
        r = 1;
    end
end