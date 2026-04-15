function out = fig_early_picture_foodweb_net(net, opts)
%FIG_EARLY_PICTURE_FOODWEB_NET
% 1x3 figure: TFL vs force-directed vs Sugiyama layered
% Input is a single network struct, e.g. from mangal_get_example_network.

if nargin < 2 || isempty(opts), opts = struct(); end

% ----------------------------- Defaults ----------------------------------
if ~isfield(opts,'ShowLabels');     opts.ShowLabels = true; end
if ~isfield(opts,'LabelFontSize');  opts.LabelFontSize = 7; end
if ~isfield(opts,'LabelColor');     opts.LabelColor = [0.6350, 0.0780, 0.1840]; end
if ~isfield(opts,'ArrowSize');      opts.ArrowSize = 0.03; end
if ~isfield(opts,'EdgeWidthMode');  opts.EdgeWidthMode = 'fixed'; end
if ~isfield(opts,'EdgeWidthRange'); opts.EdgeWidthRange = [0.6 2.2]; end

if ~isfield(opts,'PadFrac');        opts.PadFrac = 0.02; end
if ~isfield(opts,'TileSpacing');    opts.TileSpacing = 'compact'; end
if ~isfield(opts,'Padding');        opts.Padding = 'compact'; end

if ~isfield(opts,'ForceSeed');      opts.ForceSeed = 1; end
if ~isfield(opts,'LayerDirection'); opts.LayerDirection = 'up'; end

% Export defaults
if ~isfield(opts,'Export');         opts.Export = false; end
if ~isfield(opts,'ExportBase');     opts.ExportBase = 'fig_early_picture_foodweb'; end
if ~isfield(opts,'ExportDir');      opts.ExportDir = pwd; end
if ~isfield(opts,'ExportPDF');      opts.ExportPDF = true; end
if ~isfield(opts,'ExportPNG');      opts.ExportPNG = true; end
if ~isfield(opts,'PNGResolution');  opts.PNGResolution = 300; end
if ~isfield(opts,'PDFVector');      opts.PDFVector = true; end

% ------------------------------ Data -------------------------------------
if ~isfield(net,'W') || isempty(net.W)
    error('fig_early_picture_foodweb_net:NoW', 'net.W is missing or empty.');
end

W = net.W;
if issparse(W), W = full(W); end
n = size(W,1);

if isfield(net,'nodes') && ~isempty(net.nodes)
    labelsAll = cellstr(string(net.nodes));
else
    labelsAll = arrayfun(@(k) sprintf('%d',k), 1:n, 'UniformOutput', false);
end
labelsAll = labelsAll(:);

% Optional short-reader labels: only replace where you have a hand-made mapping
labelsReadable = applyShortLabelMap_(labelsAll);


% Sparse label mask
[h, infoH] = trophic_levels(W, 'ComputeCoherence', true);

topK = 4;
botK = 4;

[~, ordTopH] = sort(h, 'descend');
[~, ordBotH] = sort(h, 'ascend');

keepAB = unique([ordTopH(1:topK); ordBotH(1:botK)], 'stable');

labelsAB = labelsReadable;
labelsAB(setdiff((1:n)', keepAB)) = {''};

commonPlotNV_AB = { ...
    'ShowLabels',    opts.ShowLabels, ...
    'Labels',        labelsAB, ...
    'LabelFontSize', opts.LabelFontSize, ...
    'LabelColor',    opts.LabelColor, ...
    'LabelHalo', true, ...
    'LabelHaloPadFrac', 0.006, ...
    'LabelHaloAlpha', 0.30, ...
    'LabelHaloColor', [1 1 1], ...
    'ArrowSize',     opts.ArrowSize, ...
    'EdgeWidthMode', opts.EdgeWidthMode, ...
    'EdgeWidthRange',opts.EdgeWidthRange, ...
    'ShowBFF',       false ...
    };

% ------------------------------ Figure -----------------------------------
fig = figure('Color','w');
tl = tiledlayout(fig, 1, 3, ...
    'TileSpacing', opts.TileSpacing, ...
    'Padding', opts.Padding);

set(fig, 'Units', 'centimeters');
fig.Position(3) = 18;
fig.Position(4) = 6.5;

axs = gobjects(3,1);

% ------------------------------ Panel 1: TFL -----------------------------
axs(1) = nexttile(tl,1);

layoutOpts1 = struct( ...
    'HeightMode', 'continuous', ...
    'RepelStrength', 2, ...
    'AttractStrength', 0.01);

plotOpts1 = struct( ...
    'ShowBands', true, ...
    'FillSquare', true, ...
    'ShowLabels', opts.ShowLabels, ...
    'Labels', {labelsAB}, ...
    'LabelFontSize', opts.LabelFontSize, ...
    'LabelColor', opts.LabelColor, ...
    'LabelHalo', true, ...
    'LabelHaloPadFrac', 0.006, ...
    'LabelHaloAlpha', 0.30, ...
    'LabelHaloColor', [1 1 1], ...
    'ArrowSize', opts.ArrowSize, ...
    'EdgeWidthMode', opts.EdgeWidthMode, ...
    'EdgeWidthRange', opts.EdgeWidthRange);

[X1,Y1,compIdx1,info1] = tflPlot(W, ...
    'Parent', axs(1), ...
    'LayoutOpts', layoutOpts1, ...
    'PlotOpts', plotOpts1);

setPanelTitle_(axs(1), 'TFL (continuous, weight-aware)');
panelTag_(axs(1), 'A');
squareAxes_(axs(1), opts.PadFrac);

% ------------------------------ Panel 2: Force ---------------------------
axs(2) = nexttile(tl,2);

G = digraph(W);
rng(opts.ForceSeed);

p = plot(axs(2), G, 'Layout', 'force');
Xf = p.XData(:);
Yf = p.YData(:);
delete(p);

[Xf,Yf] = normalizeXY_(Xf,Yf);

plotTFL(W, Xf, Yf, Yf, ...
    'Parent', axs(2), ...
    'ClearAxes', true, ...
    'ShowBands', false, ...
    commonPlotNV_AB{:}, ...
    'FillSquare', true);

setPanelTitle_(axs(2), 'Force-directed');
panelTag_(axs(2), 'B');
squareAxes_(axs(2), opts.PadFrac);

% ------------------------------ Panel 3: Layered -------------------------
axs(3) = nexttile(tl,3);

p = plot(axs(3), G, 'Layout', 'layered', 'Direction', opts.LayerDirection);
Xl = p.XData(:);
Yl = p.YData(:);
delete(p);

% Optional: flip vertically if you want higher trophic nodes to appear higher
% Uncomment if needed.
% if safeCorr_(Yl, h) < 0
%     Yl = -Yl;
% end

% Convert raw Sugiyama coordinates into ordinal layer IDs
% (robust even if values are not exactly 1:6 in future examples)
[~, ~, layerOrd] = unique(Yl, 'sorted');

% Hide labels on selected Sugiyama layers (default: layer 2 only)
if ~isfield(opts, 'SugiyamaHideLayers')
    opts.SugiyamaHideLayers = 2;
end

keepC = ~ismember(layerOrd, opts.SugiyamaHideLayers);

labelsC = labelsReadable;
labelsC(~keepC) = {''};

commonPlotNV_C = commonPlotNV_AB;
idx = find(strcmpi(commonPlotNV_C, 'Labels'), 1, 'first');
if isempty(idx)
    error('Could not find ''Labels'' in commonPlotNV_C.');
end
commonPlotNV_C{idx+1} = labelsC;

% Preserve imposed vertical levels; rescale horizontal spread only
[Xl, Yl] = normalizeXOnly_(Xl, Yl);

plotTFL(W, Xl, Yl, Yl, ...
    'Parent', axs(3), ...
    'ClearAxes', true, ...
    'ShowBands', true, ...
    'levelsemantics', 'trophic', ...
     'AdjustLabels', true, ...
    'LabelRepelMode', 'vertical', ...
    'LabelRepelIters', 12, ...
    'LabelMaxShiftFrac', 0.08, ...
    commonPlotNV_C{:}, ...
    'FillSquare', true);

setPanelTitle_(axs(3), 'Sugiyama layered');
panelTag_(axs(3), 'C');
squareAxes_(axs(3), opts.PadFrac);

%{

% ------------------------------ Panel 3: Layered -------------------------
axs(3) = nexttile(tl,3);

p = plot(axs(3), G, 'Layout', 'layered', 'Direction', opts.LayerDirection);
Xl = p.XData(:);
Yl = p.YData(:);
delete(p);

% Rank positions: larger h / Y means "higher"
rankH = zeros(n,1);
rankYl = zeros(n,1);

[~, ordH] = sort(h,  'ascend');
[~, ordL] = sort(Yl, 'ascend');

rankH(ordH)  = 1:n;
rankYl(ordL) = 1:n;

% Positive value = node sits higher in Sugiyama than in TFL
promoteScore = rankYl - rankH;

% Exclude nodes already labelled in A/B
promoteScore(keepAB) = -Inf;

% Take top few promoted nodes
nExtra = 3;
[~, ordProm] = sort(promoteScore, 'descend');
extraC = ordProm(1:min(nExtra, sum(isfinite(promoteScore))));

% Optionally require actual promotion gap
extraC = extraC(promoteScore(extraC) >= 2);   % tune threshold

keepC = unique([keepAB; extraC], 'stable');

labelsC = labelsReadable;
labelsC(setdiff((1:n)', keepC)) = {''};


commonPlotNV_C = commonPlotNV_AB;

idx = find(strcmpi(commonPlotNV_C, 'Labels'), 1, 'first');
if isempty(idx)
    error('Could not find ''Labels'' in commonPlotNV_C.');
end

commonPlotNV_C{idx+1} = labelsC;

% Keep original Sugiyama levels before normalisation
yLevelsRaw = unique(Yl(:), 'sorted');

%[Xl,Yl] = normalizeXY_(Xl,Yl);
[Xl, Yl] = normalizeXOnly_(Xl, Yl);


plotTFL(W, Xl, Yl, Yl, ...
    'Parent', axs(3), ...
    'ClearAxes', true, ...
    'ShowBands', true, ...
    commonPlotNV_C{:}, ...
    'FillSquare', true);

setPanelTitle_(axs(3), 'Sugiyama layered');
panelTag_(axs(3), 'C');
squareAxes_(axs(3), opts.PadFrac);

drawnow;

%}

%{
% Halos on all panels if labels are shown
if opts.ShowLabels
    applyLabelHalos_(axs(1), 0.006, 0.3);
    applyLabelHalos_(axs(2), 0.006, 0.3);
    applyLabelHalos_(axs(3), 0.006, 0.3);
end
%}

drawnow;

% ------------------------------ Export -----------------------------------
files = struct('pdf',"",'png',"");

if opts.Export
    exportDir  = char(string(opts.ExportDir));
    exportBase = char(string(opts.ExportBase));

    if ~exist(exportDir, 'dir')
        mkdir(exportDir);
    end

    if opts.ExportPDF
        pdfFile = fullfile(exportDir, sprintf('%s.pdf', exportBase));
        if opts.PDFVector
            exportgraphics(fig, pdfFile, 'ContentType', 'vector');
        else
            exportgraphics(fig, pdfFile, 'ContentType', 'image');
        end
        files.pdf = string(pdfFile);
    end

    if opts.ExportPNG
        pngFile = fullfile(exportDir, sprintf('%s.png', exportBase));
        exportgraphics(fig, pngFile, 'Resolution', opts.PNGResolution);
        files.png = string(pngFile);
    end
end

%{

% ------------------------------ Export -----------------------------------
files = struct('pdf',"",'png',"");
if opts.Export
    if ~exist(opts.ExportDir,'dir')
        mkdir(opts.ExportDir);
    end
    base = fullfile(opts.ExportDir, opts.ExportBase);

    if opts.ExportPDF
        pdfFile = string(base) + ".pdf";
        if opts.PDFVector
            exportgraphics(fig, pdfFile, 'ContentType', 'vector');
        else
            exportgraphics(fig, pdfFile, 'ContentType', 'image');
        end
        files.pdf = pdfFile;
    end

    if opts.ExportPNG
        pngFile = string(base) + ".png";
        exportgraphics(fig, pngFile, 'Resolution', opts.PNGResolution);
        files.png = pngFile;
    end
end

%}

% ------------------------------ Outputs ----------------------------------
out = struct();
out.fig      = fig;
out.axes     = axs;
out.files    = files;
out.W        = W;
out.h        = h;
out.infoH    = infoH;
out.labelsAB = labelsAB;
out.labelsC  = labelsC;
out.keepAB   = keepAB;
out.keepC    = keepC;
out.labelsAll= labelsAll;
out.panel1   = struct('X',X1,'Y',Y1,'compIdx',compIdx1,'info',info1);
out.panel2   = struct('X',Xf,'Y',Yf);
out.panel3   = struct('X',Xl,'Y',Yl);

end

% ======================================================================
function [Xn,Yn] = normalizeXY_(X,Y)
X = X(:); Y = Y(:);
X = X - mean(X);
Y = Y - mean(Y);
sx = max(X) - min(X); if sx <= 0, sx = 1; end
sy = max(Y) - min(Y); if sy <= 0, sy = 1; end
s = max(sx, sy);
Xn = X / s;
Yn = Y / s;
end

function setPanelTitle_(ax, txt)
title(ax, txt, 'FontWeight','normal', 'Interpreter','none');
end

function panelTag_(ax, tag)
text(ax, 0.02, 0.98, tag, ...
    'Units','normalized', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top', ...
    'FontWeight','bold', ...
    'Clipping','off');
end

function squareAxes_(ax, padFrac)
xl = xlim(ax); yl = ylim(ax);
cx = mean(xl); cy = mean(yl);
rx = xl(2) - xl(1); ry = yl(2) - yl(1);
r  = max(rx, ry);
pad = padFrac * r;

xlim(ax, [cx - r/2 - pad, cx + r/2 + pad]);
ylim(ax, [cy - r/2 - pad, cy + r/2 + pad]);
daspect(ax, [1 1 1]);
pbaspect(ax, [1 1 1]);
set(ax,'XTick',[],'YTick',[]);
box(ax,'on');
end

%{
function applyLabelHalos_(ax, padFrac, faceAlpha)
% Transparent white halos behind sparse node labels.
% Robust version: create all halos first, then reorder children once.

if nargin < 2 || isempty(padFrac),   padFrac = 0.006; end
if nargin < 3 || isempty(faceAlpha), faceAlpha = 0.35; end

drawnow;  % make sure text extents are current

% Respect child order when rendering
ax.SortMethod = 'childorder';

xl = xlim(ax);
yl = ylim(ax);
padX = padFrac * (xl(2) - xl(1));
padY = padFrac * (yl(2) - yl(1));

ths = findall(ax, 'Type', 'text');

labelHs = gobjects(0,1);
haloHs  = gobjects(0,1);

for k = 1:numel(ths)
    t = ths(k);

    % Skip title / axis labels / panel tags
    if isequal(t, ax.Title) || isequal(t, ax.XLabel) || isequal(t, ax.YLabel) || isequal(t, ax.ZLabel)
        continue;
    end
    if strcmpi(t.Units, 'normalized')
        continue;
    end

    % Skip empty labels
    s = t.String;
    if (ischar(s) && isempty(s)) || (isstring(s) && all(strlength(s)==0))
        continue;
    end

    % Heuristic for node labels from plotTFL/tflPlot
    if ~strcmpi(t.HorizontalAlignment, 'center') || ~strcmpi(t.VerticalAlignment, 'bottom')
        continue;
    end

    % Measure in data units
    oldUnits = t.Units;
    t.Units = 'data';
    ext = t.Extent;   % [x y w h]
    t.Units = oldUnits;

    if numel(ext) ~= 4 || ext(3) <= 0 || ext(4) <= 0
        continue;
    end

    x0 = ext(1) - padX;
    x1 = ext(1) + ext(3) + padX;
    y0 = ext(2) - padY;
    y1 = ext(2) + ext(4) + padY;

    h = patch(ax, [x0 x1 x1 x0], [y0 y0 y1 y1], [1 1 1], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', faceAlpha, ...
        'HitTest', 'off', ...
        'PickableParts', 'none', ...
        'HandleVisibility', 'off', ...
        'Tag', 'LabelHalo');

    labelHs(end+1,1) = t; %#ok<AGROW>
    haloHs(end+1,1)  = h; %#ok<AGROW>
end

if isempty(haloHs)
    return;
end

% Reorder once:
% first children are drawn on top in axes child order
allChildren = ax.Children;
isLabel = ismember(allChildren, labelHs);
isHalo  = ismember(allChildren, haloHs);
others  = allChildren(~(isLabel | isHalo));

ax.Children = [labelHs(:); haloHs(:); others(:)];
end

%}

function labelsOut = applyShortLabelMap_(labelsIn)

labelsOut = cellstr(string(labelsIn));

from = { ...
    'Melita palmata';
    'Crangon crangon';
    'Cerastoderma edule';
    'Scrobicularia plana';
    'Carcinus maenas';
    'Peringia ulvae';
    'Cyathura carinata';
    'Alkmaria romijni';
    'Hediste diversicolor'
    };

to = { ...
    'Amphipod species';
    'Brown shrimp';
    'Cockle';
    'Furrow shell';
    'Shore crab';
    'Spire shell';
    'Isopod species';
    'Tentacled lagoon worm';
    'Ragworm'
    };


assert(numel(from) == numel(to), 'applyShortLabelMap_: map lengths differ.');

for i = 1:numel(from)
    hit = strcmpi(labelsOut, from{i});
    labelsOut(hit) = to(i);
end

end

function yLevelsN = normalizeYLevels_(yLevelsRaw)
    yLevelsRaw = yLevelsRaw(:);
    yc = mean(yLevelsRaw);

    % For layered plot, Y span is just max-min of the raw integer levels.
    sy = max(yLevelsRaw) - min(yLevelsRaw);
    if sy <= 0, sy = 1; end

    % But your normalizeXY_ uses s = max(sx,sy), not sy alone.
    % So this helper should really be passed sx too if you want exact agreement.
    % Simpler: use the same normalization function on a dummy X.
    error('Use normalizeYLevelsFromXY_ instead.');
end

function [Xn, Yn] = normalizeXOnly_(X, Y)
X = X(:);
Y = Y(:);

X = X - mean(X);
sx = max(X) - min(X);
if sx <= 0, sx = 1; end

Xn = X / sx;
Yn = Y;
end