function [Xplot, Yplot, compIdx, info, scene, geom] = tflPlot(W, varargin)
%TFLPLOT  Convenience wrapper for trophicLayoutMulti + plotTFL.
%
%   [Xplot, Yplot, compIdx, info, scene, geom] = tflPlot(W, ...)
%
% Description
% -----------
% tflPlot computes a trophic layout using trophicLayoutMulti and then
% renders it using plotTFL.
%
% The first two outputs, Xplot and Yplot, are the coordinates actually
% rendered by plotTFL after scene construction and any plotting-time
% rescaling. These are therefore the most natural coordinates to use for
% downstream plotting or annotation.
%
% The raw layout coordinates returned by trophicLayoutMulti before plotting
% transforms are preserved in:
%
%   info.XLayoutRaw
%   info.YLayoutRaw
%
% Name-value pairs
% ----------------
%   'LayoutOpts' : scalar struct passed to trophicLayoutMulti
%   'PlotOpts'   : scalar struct passed to plotTFL
%                  (for example NodeSizeData, NodeSizeScaleData, Labels,
%                  RenderMode options, etc.)
%   'Parent'     : axes handle to plot into
%
% Any other name-value pairs are treated as plot options and passed through
% to plotTFL.
%
% Default semantics
% -----------------
% If LevelSemantics is not explicitly supplied in PlotOpts:
%
%   - internally computed trophic levels  -> 'trophic'
%   - user-supplied hProvided             -> 'generic'
%
% Outputs
% -------
%   Xplot, Yplot : rendered node coordinates returned by plotTFL
%   compIdx      : component assignment from trophicLayoutMulti
%   info         : layout metadata from trophicLayoutMulti, augmented with
%                  info.XLayoutRaw and info.YLayoutRaw
%   scene        : scene struct returned by plotTFL
%   geom         : edge geometry struct returned by plotTFL
%
% Example
% -------
%   [Xplot, Yplot, compIdx, info, scene, geom] = tflPlot(W, ...
%       'PlotOpts', struct( ...
%           'NodeSizeData', uThisPanel, ...
%           'NodeSizeScaleData', uShared, ...
%           'RenderMode', 'overlay'));

    n = size(W,1);
    if size(W,2) ~= n
        error('tflPlot:WNotSquare', 'W must be square.');
    end
    if ~isnumeric(W)
        error('tflPlot:WNotNumeric', 'W must be numeric.');
    end

    layoutOptsStruct = struct();
    layoutArgs = {};
    plotArgs   = {};
    parentAx   = [];

    if mod(numel(varargin),2) ~= 0
        error('tflPlot:BadArgs', ...
            'Optional arguments must be name/value pairs.');
    end

    k = 1;
    while k <= numel(varargin)
        name = varargin{k};
        val  = varargin{k+1};

        if ~ischar(name) && ~isstring(name)
            error('tflPlot:BadParamName', ...
                'Parameter names must be strings.');
        end

        switch lower(char(name))
            case 'layoutopts'
                layoutOptsStruct = validateScalarStruct_(val, 'LayoutOpts');
                layoutArgs = [layoutArgs, structToNameValue_(layoutOptsStruct)]; %#ok<AGROW>

            case 'plotopts'
                plotArgs = [plotArgs, structToNameValue_( ...
                    validateScalarStruct_(val, 'PlotOpts'))]; %#ok<AGROW>

            case 'parent'
                parentAx = val;

            otherwise
                % Default: treat as plot option
                plotArgs = [plotArgs, {name, val}]; %#ok<AGROW>
        end

        k = k + 2;
    end

    if ~isempty(parentAx) && ~nameValueHas_(plotArgs, 'Parent')
        plotArgs = [plotArgs, {'Parent', parentAx}]; %#ok<AGROW>
    end

    % Default semantics:
    % - internally computed trophic levels -> 'trophic'
    % - user-supplied hProvided           -> 'generic'
    % unless explicitly overridden in plot options.
    if ~nameValueHas_(plotArgs, 'LevelSemantics')
        if isfield(layoutOptsStruct, 'hProvided') && ~isempty(layoutOptsStruct.hProvided)
            plotArgs = [plotArgs, {'LevelSemantics', 'generic'}]; %#ok<AGROW>
        else
            plotArgs = [plotArgs, {'LevelSemantics', 'trophic'}]; %#ok<AGROW>
        end
    end

    [Xraw, Yraw, compIdx, info] = trophicLayoutMulti(W, layoutArgs{:});

    [Xplot, Yplot, scene, geom] = plotTFL(W, Xraw, Yraw, info.h, plotArgs{:});

    % Preserve raw layout coordinates for downstream use if needed.
    info.XLayoutRaw = Xraw;
    info.YLayoutRaw = Yraw;
end

function s = validateScalarStruct_(s, argName)
    if isempty(s)
        s = struct();
        return;
    end
    if ~isstruct(s) || numel(s) ~= 1
        error('tflPlot:OptsNotStruct', ...
            '%s must be a scalar struct.', argName);
    end
end

function nv = structToNameValue_(s)
    if isempty(s)
        nv = {};
        return;
    end
    if ~isstruct(s) || numel(s) ~= 1
        error('tflPlot:OptsNotStruct', ...
            'LayoutOpts and PlotOpts must be scalar structs.');
    end

    fn = fieldnames(s);
    nv = cell(1, 2*numel(fn));
    for i = 1:numel(fn)
        nv{2*i-1} = fn{i};
        nv{2*i}   = s.(fn{i});
    end
end

function tf = nameValueHas_(nv, key)
    tf = false;
    if isempty(nv), return; end
    for i = 1:2:numel(nv)
        if (ischar(nv{i}) || isstring(nv{i})) && strcmpi(char(nv{i}), key)
            tf = true;
            return;
        end
    end
end