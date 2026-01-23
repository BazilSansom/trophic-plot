function [X, Y, compIdx, info] = tflPlot(W, varargin)
%TFLPLOT  Convenience wrapper: full TFL layout + plotting in one call.
%
%   [X, Y, compIdx, info] = TFLPLOT(W, ...)
%
% Name–value pairs:
%   'LayoutOpts' : scalar struct passed to TROPHICLAYOUTMULTI
%   'PlotOpts'   : scalar struct passed to PLOTTFL
%   'Parent'     : axes handle to plot into (shorthand for PlotOpts.Parent)
%
% Any other name–value pairs are treated as layout options:

% INPUT
%   W : adjacency matrix (n x n), directed, weighted or unweighted.
%
% OPTIONAL NAME–VALUE PAIRS
%   'LayoutOpts' : scalar struct of options passed to TROPHICLAYOUTMULTI
%   'PlotOpts'   : scalar struct of options passed to PLOTTFL
%
%   Any other name–value pairs are treated as layout options and passed
%   directly to TROPHICLAYOUTMULTI. This is mainly for backwards
%   convenience; for cleaner code, prefer using LayoutOpts / PlotOpts.
%
% EXAMPLES
%
%   % Simplest usage
%   tflPlot(W);
%
%   % With layout options
%   tflPlot(W, 'LayoutOpts', struct( ...
%       'UseBarycentre',      true, ...
%       'CoherenceThreshold', 0.8, ...
%       'UseXSmooth',         true));
%
%   % With plot options
%   tflPlot(W, ...
%       'LayoutOpts', struct('UseBarycentre', false), ...
%       'PlotOpts',   struct('ShowLabels', true, 'NodeSize', 60));
%
% Dependencies (in this repository):
%   trophicLayoutMulti.m
%   plotTFL.m
%
% -------------------------------------------------------------------------

    % ---------- basic checks ----------
    n = size(W,1);
    if size(W,2) ~= n
        error('tflPlot:WNotSquare', 'W must be a square adjacency matrix.');
    end
    if ~isnumeric(W)
        error('tflPlot:WNotNumeric', 'W must be numeric.');
    end

    layoutArgs = {};
    plotArgs   = {};

    % convenience capture
    parentAx   = [];

    if mod(numel(varargin),2) ~= 0
        error('tflPlot:BadArgs', 'Optional arguments must be name/value pairs.');
    end

    k = 1;
    while k <= numel(varargin)
        name = varargin{k};
        val  = varargin{k+1};

        if ~ischar(name) && ~isstring(name)
            error('tflPlot:BadParamName', 'Parameter names must be strings.');
        end
        lname = lower(char(name));

        switch lname
            case 'layoutopts'
                layoutArgs = [layoutArgs, structToNameValue(val)]; %#ok<AGROW>

            case 'plotopts'
                plotArgs = [plotArgs, structToNameValue(val)]; %#ok<AGROW>

            case 'parent'
                parentAx = val;

            otherwise
                % By default, treat everything else as a layout option
                layoutArgs = [layoutArgs, {name, val}]; %#ok<AGROW>
        end

        k = k + 2;
    end

    % If user passed Parent, push it into plot args (unless PlotOpts already has Parent)
    if ~isempty(parentAx)
        % Only add if not already specified in plotArgs
        if ~nameValueHas_(plotArgs, 'Parent')
            plotArgs = [plotArgs, {'Parent', parentAx}]; %#ok<AGROW>
        end
    end

    % ---------- 1) Layout ----------
    [X, Y, compIdx, info] = trophicLayoutMulti(W, layoutArgs{:});

    % ---------- 2) Plot ----------
    plotTFL(W, X, Y, info.h_raw, plotArgs{:});
    %plotTFL(W, X, Y, Y, plotArgs{:});

end

% =====================================================================
function nv = structToNameValue(s)
    if isempty(s)
        nv = {};
        return;
    end
    if ~isstruct(s) || numel(s) ~= 1
        error('tflPlot:OptsNotStruct', 'LayoutOpts and PlotOpts must be scalar structs.');
    end

    fn = fieldnames(s);
    nv = cell(1, 2*numel(fn));
    for i = 1:numel(fn)
        nv{2*i-1} = fn{i};
        nv{2*i}   = s.(fn{i});
    end
end

function tf = nameValueHas_(nv, key)
%NAMEVALUEHAS_ true if nv contains parameter name key (case-insensitive)
    tf = false;
    if isempty(nv), return; end
    for i = 1:2:numel(nv)
        if (ischar(nv{i}) || isstring(nv{i})) && strcmpi(char(nv{i}), key)
            tf = true;
            return;
        end
    end
end
