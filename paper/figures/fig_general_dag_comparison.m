function out = fig_general_dag_comparison(opts)
%FIG_GENERAL_DAG_COMPARISON  TFL vs Sugiyama layered layout on a (chosen) DAG.
%
% out.fig, out.axes, out.W, out.info, out.files
%
% -- Example usage --:
%
%fig_general_dag_comparison(struct( ...
%    'ShowLabels', true, ...
%    'Export', true, ...
%    'ExportDir', fullfile(pwd,'figures'), ...
%    'ExportBase','fig_general_dag_comparison', ...
%    'PNGResolution', 300));

if nargin < 1 || isempty(opts), opts = struct(); end

% Layout/plot
def.ShowLabels  = true;
def.TileSpacing = 'compact';
def.Padding     = 'compact';
def.PadFrac     = 0.02;

% Export
def.Export        = false;
def.ExportDir     = "figures";
def.ExportBase    = "fig_general_dag_comparison";
def.ExportPDF     = true;
def.ExportPNG     = true;
def.PNGResolution = 300;
def.PDFVector     = true;

opts = applyDefaults(opts, def);

% -------------------- network --------------------
E = [1  3
     2  3
     3  4
     2  5
     3  5];
W = edgelistToAdjacency(E);

[h_raw, infoLevels] = trophic_levels(W, 'ComputeCoherence', true);

% -------------------- figure --------------------
fig = paperFigure('double');

tl = tiledlayout(fig, 1, 2, ...
    'TileSpacing', opts.TileSpacing, ...
    'Padding',     opts.Padding);

axs = gobjects(2,1);

% --- (1) TFL ---
axs(1) = nexttile(tl,1);
tflPlot(W, 'plotopts', struct('ShowLabels', opts.ShowLabels));
fixTiledAxesCommon(axs(1));
makeAxesDataSquare(axs(1), opts.PadFrac);
S = paperStyle();
setPanelTitle(axs(1), 'Trophic Flow Layout (TFL)', S);

% --- (2) Sugiyama layered (re-plotted with plotTFL aesthetics + bands) ---
axs(2) = nexttile(tl,2);
g = digraph(W);

p = plot(axs(2), g, 'Layout','layered', 'Direction','up');
Xg = p.XData(:);
Xg = max(Xg) - Xg;   % flip left-right to match TFL
Yg = p.YData(:);
delete(p);

% Build integer-ish layers from Y (discrete levels from layered layout)
[yVals, ~, layerIdx] = unique(Yg, 'sorted'); %#ok<ASGLU>
hSug = layerIdx - min(layerIdx);
Yg_int = double(hSug);

plotTFL(W, Xg, Yg_int, Yg_int, ...
    'ShowBands',  true, ...
    'ShowLabels', opts.ShowLabels);

fixTiledAxesCommon(axs(2));
makeAxesDataSquare(axs(2), opts.PadFrac);
S = paperStyle();
setPanelTitle(axs(2), 'Sugiyama layered layout', S);

% -------------------- optional export --------------------
files = struct('pdf',"",'png',"");
if opts.Export
    if ~exist(opts.ExportDir,'dir'), mkdir(opts.ExportDir); end
    base = fullfile(opts.ExportDir, opts.ExportBase);

    if opts.ExportPDF
        pdfFile = base + ".pdf";
        if opts.PDFVector
            exportgraphics(fig, pdfFile, 'ContentType','vector', 'BackgroundColor','white');
        else
            exportgraphics(fig, pdfFile, 'ContentType','image',  'BackgroundColor','white');
        end
        files.pdf = pdfFile;
    end

    if opts.ExportPNG
        pngFile = base + ".png";
        exportgraphics(fig, pngFile, 'Resolution', opts.PNGResolution, 'BackgroundColor','white');
        files.png = pngFile;
    end
end

% -------------------- outputs --------------------
out = struct();
out.fig   = fig;
out.axes  = axs;
out.W     = W;

out.info = struct();
out.info.h_raw = h_raw;
out.info.levelsInfo = infoLevels;
if isfield(infoLevels,'F0') && isfield(infoLevels,'C')
    out.info.F0 = infoLevels.F0;
    out.info.C  = infoLevels.C;
end

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
