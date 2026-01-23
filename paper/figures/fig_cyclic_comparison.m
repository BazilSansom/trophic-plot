function out = fig_cyclic_comparison(opts)
%FIG_CYCLIC_COMPARISON  Regime 3 comparison figure: cyclic digraph.
%
% Left:  Trophic Flow Layout (TFL)
% Right: Sugiyama layered layout (MATLAB layered), re-plotted with plotTFL aesthetics
%
% out.fig, out.axes, out.W, out.info, out.files

% -- Example usage --:
%
%fig_cyclic_comparison(struct( ...
%    'ShowLabels', true, ...
%    'Export', true, ...
%    'ExportDir', fullfile(pwd,'figures'), ...
%    'ExportBase','fig_cyclic_comparison', ...
%    'PNGResolution', 300));

% -------------------- defaults --------------------
if nargin < 1 || isempty(opts), opts = struct(); end

def.ShowLabels   = true;
def.TileSpacing  = 'compact';
def.Padding      = 'compact';
def.PadFrac      = 0.02;

% Export
def.Export        = false;
def.ExportDir     = fullfile("paper","outputs","figs");
def.ExportBase    = "fig_cyclic_comparison";
def.ExportPDF     = true;
def.ExportPNG     = false;
def.PNGResolution = 300;
def.PDFVector     = true;

% Example selection / inputs
def.W = [];              % if provided, overrides ExampleEdgeList
def.UseWeights = false;  % if true, treat 3rd col as weights if present
opts = applyDefaults(opts, def);

S = paperStyle();

% -------------------- build example cyclic network --------------------
% A small, readable cyclic example:
%   - core 3-cycle: 1 -> 2 -> 3 -> 1
%   - plus a feed-forward branch and a 2-cycle to show bidirectional pair rendering
%
% Columns: [src dst] (or [src dst w] if you want)
E = [
     1,3;
     2,3;
     5,3;
     1,4;
     2,5;
     3,5;
     4,6;
     5,6];

if ~isempty(opts.W)
    W = opts.W;
else
    W = edgelistToAdjacency(E);
end

% -------------------- trophic analysis (for info / optional annotations) --------------------
[h_raw, infoLevels] = trophic_levels(W, 'ComputeCoherence', true);

% -------------------- figure --------------------
fig = paperFigure('double');
tl = tiledlayout(fig, 1, 2, ...
    'TileSpacing', opts.TileSpacing, ...
    'Padding',     opts.Padding);

axs = gobjects(2,1);

% ============================================================
% (1) TFL
% ============================================================
axs(1) = nexttile(tl,1);
tflPlot(W, 'plotopts', struct('ShowLabels', opts.ShowLabels));
fixTiledAxesCommon(axs(1));
makeAxesDataSquare(axs(1), opts.PadFrac);
setPanelTitle(axs(1), "Trophic Flow Layout (TFL)", S);

% ============================================================
% (2) "Sugiyama layered (MATLAB)" but drawn with plotTFL aesthetics
% ============================================================
axs(2) = nexttile(tl,2);

g = digraph(W);

% Compute MATLAB layered coords (MATLAB will internally handle cycles)
p = plot(axs(2), g, 'Layout','layered', 'Direction','up');
Xg = p.XData(:);
Yg = p.YData(:);
delete(p);

% Optional: flip horizontally to align “viewing orientation” with TFL (like you did before)
Xg = max(Xg) - Xg;

% Layer reference lines: for cyclic graphs the Y returned may still be discrete-ish;
% if you prefer integer banding, discretise:
[yVals,~,layerIdx] = unique(Yg, 'sorted');
hSug   = double(layerIdx - min(layerIdx));
Yg_int = hSug;

plotTFL(W, Xg, Yg_int, Yg_int, ...
    'ShowBands',  true, ...
    'ShowLabels', opts.ShowLabels);

fixTiledAxesCommon(axs(2));
makeAxesDataSquare(axs(2), opts.PadFrac);
setPanelTitle(axs(2), "Sugiyama layered layout", S);

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
out.fig   = fig;
out.axes  = axs;
out.W     = W;
out.info  = struct('h_raw',h_raw,'levelsInfo',infoLevels);
out.files = files;

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
