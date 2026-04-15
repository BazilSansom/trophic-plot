function out = plotTFL_focal_corridor(W, X, Y, h, focalIdx, labels, corrOpts, plotOpts)
%PLOTTFL_FOCAL_CORRIDOR
% Two-pass plot:
%   1) full graph faded
%   2) corridor overlay highlighted
%
% Inputs:
%   X, Y : already-computed display coordinates
%   h    : trophic levels used for corridor selection
%
% Bands are aligned to the supplied plotted vertical coordinates Y.

    if nargin < 7 || isempty(corrOpts), corrOpts = struct(); end
    if nargin < 8 || isempty(plotOpts), plotOpts = struct(); end

    ax = [];
    if isfield(plotOpts,'Parent'), ax = plotOpts.Parent; end
    if isempty(ax) || ~isgraphics(ax,'axes'), ax = gca; end

    % Use true trophic levels for corridor selection
    corr = tfl_focal_corridor_mask(W, h, focalIdx, corrOpts);

    % Use plotted vertical coordinates for reference-band placement
    %hPlot = Y(:);

    % ---- Background full graph ----
    bgOpts = struct( ...
        'Parent', ax, ...
        'ClearAxes', true, ...
        'ShowBands', true, ...
        'LevelSemantics', 'trophic', ...
        'FillSquare', false, ...
        'ShowLabels', false, ...
        'NodeColor', [0.78 0.82 0.90], ...
        'UpEdgeColor', [0.50 0.50 0.50], ...
        'DownEdgeColor', [0.75 0.75 0.75], ...
        'EdgeAlpha', 0.18, ...
        'ArrowSize', 0.03, ...
        'ArrowSizeMode', 'fixed', ...
        'EdgeWidthMode', 'fixed');

    bgOpts = mergeStructs_(bgOpts, rmfieldIfExists_(plotOpts, {'Parent'}));

    bgNV = structToNV_(bgOpts);
    %[Xdisp, Ydisp] = plotTFL(W, X, Y, hPlot, bgNV{:});
    [Xdisp, Ydisp] = plotTFL(W, X, Y, h, bgNV{:});


    xl = xlim(ax);
    yl = ylim(ax);

    % ---- Corridor overlay subgraph ----
    idx = find(corr.nodeMask);
    if isempty(idx)
        holdState = ishold(ax);
        hold(ax, 'on');
        scatter(ax, Xdisp(focalIdx), Ydisp(focalIdx), 80, ...
            'MarkerFaceColor', [0.85 0.20 0.20], ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 1.0);
        if ~holdState
            hold(ax, 'off');
        end

        xlim(ax, xl);
        ylim(ax, yl);

        out = struct();
        out.corr = corr;
        out.idx = idx;
        out.Xdisp = Xdisp;
        out.Ydisp = Ydisp;
        return;
    end

    Wc = W(idx, idx);
    Ec = corr.edgeMask(idx, idx);
    Wc(~Ec) = 0;

    Xc = Xdisp(idx);
    Yc = Ydisp(idx);
    %hPlotC = Yc(:);

    labelsC = repmat({''}, numel(idx), 1);
    if nargin >= 6 && ~isempty(labels)
        for k = 1:numel(idx)
            labelsC{k} = labels{idx(k)};
        end
    end

    % ---- Foreground corridor ----
    fgOpts = struct( ...
        'Parent', ax, ...
        'ClearAxes', false, ...
        'ShowBands', false, ...
        'LevelSemantics', 'trophic', ...
        'FillSquare', false, ...
        'ShowLabels', true, ...
        'Labels', {labelsC}, ...
        'LabelFontSize', 7, ...
        'LabelOffset', [0 0.02], ...
        'LabelOffsetMode', 'frac', ...
        'LabelHalo', true, ...
        'LabelHaloPadFrac', 0.006, ...
        'LabelHaloAlpha', 0.30, ...
        'LabelHaloColor', [1 1 1], ...
        'NodeColor', [0.10 0.30 0.75], ...
        'UpEdgeColor', [0.15 0.15 0.15], ...
        'DownEdgeColor', [0.55 0.55 0.55], ...
        'EdgeAlpha', 0.95, ...
        'ArrowSize', 0.03, ...
        'ArrowSizeMode', 'fixed', ...
        'EdgeWidthMode', 'fixed');

    fgNV = structToNV_(fgOpts);
    %plotTFL(Wc, Xc, Yc, hPlotC, fgNV{:});
    hc = h(idx);
    plotTFL(Wc, Xc, Yc, hc, fgNV{:});

    xlim(ax, xl);
    ylim(ax, yl);

    % ---- Highlight focal node ----
    holdState = ishold(ax);
    hold(ax, 'on');

    scatter(ax, Xdisp(focalIdx), Ydisp(focalIdx), 80, ...
        'MarkerFaceColor', [0.85 0.20 0.20], ...
        'MarkerEdgeColor', 'k', ...
        'LineWidth', 1.0);

    if ~holdState
        hold(ax, 'off');
    end

    xlim(ax, xl);
    ylim(ax, yl);

    out = struct();
    out.corr = corr;
    out.idx = idx;
    out.Xdisp = Xdisp;
    out.Ydisp = Ydisp;
end

function s = mergeStructs_(a,b)
s = a;
fn = fieldnames(b);
for k = 1:numel(fn)
    s.(fn{k}) = b.(fn{k});
end
end

function s = rmfieldIfExists_(s, fns)
for k = 1:numel(fns)
    if isfield(s, fns{k})
        s = rmfield(s, fns{k});
    end
end
end

function nv = structToNV_(s)
nv = {};
fn = fieldnames(s);
for k = 1:numel(fn)
    nv(end+1:end+2) = {fn{k}, s.(fn{k})}; %#ok<AGROW>
end
end
