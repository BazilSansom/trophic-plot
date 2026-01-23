function out = fig_tfl_vs_sugiyama(opts)
%FIG_TFL_VS_SUGIYAMA  Compare TFL vs Sugiyama layered layout for an input network.
%
% Produces a 1x2 figure consistent with your paper figure style:
%   Left:  Trophic Flow Layout (TFL) via tflPlot
%   Right: Sugiyama layered layout (MATLAB) but drawn with plotTFL aesthetics
%
% INPUT (opts fields)
%   Network input (choose ONE):
%     opts.W        : NxN adjacency (sparse or full)
%     opts.E        : edge list Nx2 or Nx3 (passed to edgelistToAdjacency)
%
%   Plot options:
%     opts.ShowLabels (true)
%     opts.PadFrac    (0.02)
%     opts.TileSpacing ('compact')
%     opts.Padding     ('compact')
%
%   Panel titles:
%     opts.LeftTitle  ("Trophic Flow Layout (TFL)")
%     opts.RightTitle ("Sugiyama layered (MATLAB)")
%     opts.TitleWeight ('bold') or ('normal') depending on your style
%
%   Sugiyama options:
%     opts.Direction ('up')     : 'up'/'down' etc supported by MATLAB
%     opts.FlipRightX (true)    : mirror right panel in x to align orientation
%     opts.ShowSugiyamaBands (false) : show integer-ish bands on Sugiyama panel
%
%   Export options:
%     opts.Export (false)
%     opts.ExportDir ("figures")
%     opts.ExportBase ("fig_tfl_vs_sugiyama")
%     opts.ExportPDF (true)
%     opts.ExportPNG (false)
%     opts.PNGResolution (300)
%     opts.PDFVector (true)
%
% OUTPUT
%   out.fig, out.axes, out.W, out.files, out.meta
%
% Dependencies (your repo):
%   paperFigure, paperStyle, fixTiledAxesCommon, makeAxesDataSquare
%   tflPlot, plotTFL, edgelistToAdjacency

% -------------------- defaults --------------------
if nargin < 1 || isempty(opts), opts = struct(); end

def.ShowLabels = true;
def.PadFrac = 0.02;
def.TileSpacing = 'compact';
def.Padding = 'compact';

def.LeftTitle  = "Trophic Flow Layout (TFL)";
def.RightTitle = "Sugiyama layered (MATLAB)";
def.TitleWeight = "bold";

def.Direction = "up";
def.FlipRightX = true;
def.ShowSugiyamaBands = false;

def.Export = false;
def.ExportDir = "figures";
def.ExportBase = "fig_tfl_vs_sugiyama";
def.ExportPDF = true;
def.ExportPNG = false;
def.PNGResolution = 300;
def.PDFVector = true;

opts = applyDefaults(opts, def);

% -------------------- build W --------------------
if isfield(opts,'W') && ~isempty(opts.W)
    W = opts.W;
elseif isfield(opts,'E') && ~isempty(opts.E)
    W = edgelistToAdjacency(opts.E);
else
    error('fig_tfl_vs_sugiyama:NoInput', 'Provide opts.W (adjacency) or opts.E (edge list).');
end

% Basic checks
N = size(W,1);
if size(W,2) ~= N
    error('fig_tfl_vs_sugiyama:WNotSquare','W must be square.');
end

% -------------------- figure + style --------------------
fig = paperFigure('double');
try
    S = paperStyle(); %#ok<NASGU>
catch
    S = struct();
end

tl = tiledlayout(fig, 1, 2, ...
    'TileSpacing', opts.TileSpacing, ...
    'Padding',     opts.Padding);

axs = gobjects(2,1);

% -------------------- (1) TFL --------------------
axs(1) = nexttile(tl, 1);
tflPlot(W, 'plotopts', struct('ShowLabels', opts.ShowLabels));
fixTiledAxesCommon(axs(1));
makeAxesDataSquare(axs(1), opts.PadFrac);

setPanelTitleCompat(axs(1), opts.LeftTitle, opts.TitleWeight);

% -------------------- (2) Sugiyama layered coords, re-plotted with plotTFL style --------------------
axs(2) = nexttile(tl, 2);

g = digraph(W);

% Compute layered coords
p = plot(axs(2), g, 'Layout','layered', 'Direction', char(opts.Direction));
Xg = p.XData(:);
Yg = p.YData(:);
delete(p);

% Optional mirror in x to align with TFL orientation (helps reduce cognitive load)
if opts.FlipRightX
    Xg = max(Xg) - Xg;
end

% Build "layer index" from discrete Y levels if requested (for bands)
if opts.ShowSugiyamaBands
    [~, ~, layerIdx] = unique(Yg, 'sorted');
    hSug = layerIdx - min(layerIdx);     % basal layer = 0
    Yplot = double(hSug);
    hplot = double(hSug);
    showBands = true;
else
    % No bands: use raw Y for both Y and "h" (bands off anyway)
    Yplot = Yg;
    hplot = Yg;
    showBands = false;
end

plotTFL(W, Xg, Yplot, hplot, ...
    'ShowBands',  showBands, ...
    'ShowLabels', opts.ShowLabels);

fixTiledAxesCommon(axs(2));
makeAxesDataSquare(axs(2), opts.PadFrac);

setPanelTitleCompat(axs(2), opts.RightTitle, opts.TitleWeight);

% -------------------- optional export --------------------
files = struct('pdf',"",'png',"");
if opts.Export
    if ~exist(opts.ExportDir,'dir')
        mkdir(opts.ExportDir);
    end
    base = fullfile(opts.ExportDir, opts.ExportBase);

    if opts.ExportPDF
        pdfFile = base + ".pdf";
        if opts.PDFVector
            exportgraphics(fig, pdfFile, 'ContentType','vector');
        else
            exportgraphics(fig, pdfFile, 'ContentType','image');
        end
        files.pdf = pdfFile;
    end

    if opts.ExportPNG
        pngFile = base + ".png";
        exportgraphics(fig, pngFile, 'Resolution', opts.PNGResolution);
        files.png = pngFile;
    end
end

% -------------------- outputs --------------------
out = struct();
out.fig  = fig;
out.axes = axs;
out.W    = W;
out.files = files;

out.meta = struct();
out.meta.N = N;
out.meta.direction = string(opts.Direction);
out.meta.flipRightX = logical(opts.FlipRightX);
out.meta.showSugiyamaBands = logical(opts.ShowSugiyamaBands);

end

% ======================================================================
function s = applyDefaults(s, def)
fn = fieldnames(def);
for k = 1:numel(fn)
    f = fn{k};
    if ~isfield(s,f) || isempty(s.(f))
        s.(f) = def.(f);
    end
end
end

% ======================================================================
function setPanelTitleCompat(ax, txt, weight)
% Uses your setPanelTitle(ax,txt,S) if it exists; otherwise falls back to title().
if nargin < 4, weight = "bold"; end
if exist('setPanelTitle','file') == 2
    try
        S = paperStyle();
    catch
        S = struct();
    end
    setPanelTitle(ax, txt, S);
else
    title(ax, txt, 'FontWeight', char(weight), 'Interpreter','none');
end
end
