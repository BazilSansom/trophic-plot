function out = fig_rand_cyclic_net_comparison(opts)
%FIG_REGIME2_GENERAL_DAG_COMPARISON  TFL vs Sugiyama layered layout on a GPPM DAG.
%
% Generates a random DAG using randGPPM (ForceDAG=true), then plots:
%   Left:  Trophic Flow Layout (TFL)
%   Right: Sugiyama-style layered layout (MATLAB layered) drawn with plotTFL aesthetics
%
% Returns out.fig, out.axes, out.W, out.info, out.files, out.meta (incl. seed).

% -------------------- defaults --------------------
if nargin < 1 || isempty(opts), opts = struct(); end

% Network draw
def.N = 8;
def.p = 0.18;
def.T = 1.0;
def.Seed = [];                 % if empty => random
def.ForceDAG = false;           % important for regime-2 DAG story
def.Weighted = false;          % start unweighted; you can switch later
def.NumBasal = 3;

% Layout/plot
def.ShowLabels = true;
def.TileSpacing = 'compact';
def.Padding = 'compact';
def.PadFrac = 0.02;

% Export
def.Export = false;
def.ExportDir = "figures";
def.ExportBase = "fig_regime2_general_dag";
def.ExportPDF = true;
def.ExportPNG = true;
def.PNGResolution = 300;
def.PDFVector = true;

opts = applyDefaults(opts, def);

% -------------------- generate network --------------------
genArgs = { ...
    'ForceDAG',       logical(opts.ForceDAG), ...
    'Weighted',       logical(opts.Weighted), ...
    'NumBasal',       opts.NumBasal ...
};

if ~isempty(opts.Seed)
    genArgs = [genArgs, {'Seed', opts.Seed}];
end

[W, metaGen] = randGPPM(opts.N, opts.p, opts.T, genArgs{:});
W = sparse(W+0);

% Basic diagnostics (optional but useful while hunting examples)
[h_raw, infoLevels] = trophic_levels(W, 'ComputeCoherence', true);

% -------------------- figure --------------------
fig = figure('Color','w');

tl = tiledlayout(fig, 1, 2, ...
    'TileSpacing', opts.TileSpacing, ...
    'Padding',     opts.Padding);

axs = gobjects(2,1);

% Make figure roughly wide enough for two square panels
set(fig,'Units','centimeters');
pos = fig.Position;
height = pos(4);
fig.Position = [pos(1) pos(2) 2.05*height height];  % ~2 squares + small gutter

% --- (1) TFL ---
axs(1) = nexttile(tl,1);
tflPlot(W, 'plotopts', struct('ShowLabels', opts.ShowLabels));
fixTiledAxesCommon(axs(1));
makeAxesDataSquare(axs(1), opts.PadFrac);
title(axs(1), 'Trophic Flow Layout (TFL)', 'Interpreter','none');

% --- (2) Sugiyama layered (but re-plotted with plotTFL aesthetics + bands) ---
axs(2) = nexttile(tl,2);
g = digraph(W);

% Compute layered coords
p = plot(axs(2), g, 'Layout','layered', 'Direction','up');
Xg = p.XData(:);
Xg = max(Xg) - Xg;  % flip left-right to match TFL
Yg = p.YData(:);
delete(p);

% ---- Build "integer-ish" layer index from Yg ----
% 1) Identify layers by unique Y values (MATLAB layered uses discrete Y)
[yVals, ~, layerIdx] = unique(Yg, 'sorted');

% 2) Make basal layer = 0 and step = 1
%    With Direction='up', basal is usually the minimum Y layer.
%    If you ever flip Direction, just swap min/max logic accordingly.
hSug = layerIdx - min(layerIdx);

% 3) Re-embed Y to match integer spacing, but keep *relative* ordering
%    (This keeps plot clean and makes bands meaningful.)
Yg_int = double(hSug);

% Plot with plotTFL aesthetics; show bands at integer "layers"
plotTFL(W, Xg, Yg_int, Yg_int, ...
    'ShowBands',  true, ...
    'ShowLabels', opts.ShowLabels);

fixTiledAxesCommon(axs(2));
makeAxesDataSquare(axs(2), opts.PadFrac);
title(axs(2), 'Sugiyama layered (MATLAB)', 'Interpreter','none');


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

out.meta = struct();
out.meta.randGPPM = metaGen;
out.meta.optsUsed = opts;

out.info = struct();
out.info.h_raw = h_raw;
out.info.levelsInfo = infoLevels;

% handy for screening candidate networks
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
