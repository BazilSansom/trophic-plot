function out = fig_perfectly_layered_comparison(opts)
%FIG_PERFECTLY_LAYERED_COMPARISON  TFL vs Sugiyama layered layout on a perfectly layered DAG.
%
% out.fig, out.axes, out.files

% -- Example usage --:
%
%fig_perfectly_layered_comparison(struct( ...
%    'ShowLabels', true, ...
%    'Export', true, ...
%    'ExportDir', fullfile(pwd,'figures'), ...
%    'ExportBase','fig_perfectly_layered_comparison', ...
%    'PNGResolution', 300));


if nargin < 1 || isempty(opts), opts = struct(); end

% plotting
if ~isfield(opts,'ShowLabels');  opts.ShowLabels  = true;  end
if ~isfield(opts,'PadFrac');     opts.PadFrac     = 0.02; end
if ~isfield(opts,'TileSpacing'); opts.TileSpacing = 'compact'; end
if ~isfield(opts,'Padding');     opts.Padding     = 'compact'; end

% export
if ~isfield(opts,'Export');         opts.Export         = false; end
if ~isfield(opts,'ExportDir');      opts.ExportDir      = pwd;   end
if ~isfield(opts,'ExportBase');     opts.ExportBase     = "fig_perfectly_layered"; end
if ~isfield(opts,'ExportPDF');      opts.ExportPDF      = true;  end
if ~isfield(opts,'ExportPNG');      opts.ExportPNG      = false; end
if ~isfield(opts,'PDFVector');      opts.PDFVector      = true;  end
if ~isfield(opts,'PNGResolution');  opts.PNGResolution  = 300;   end

% ================= data =================
E = [
    1 2
    3 2
    4 3
    5 4
    5 6
];
W = edgelistToAdjacency(E);

% ================= figure =================
fig = paperFigure('double');

tl = tiledlayout(fig,1,2, ...
    'TileSpacing', opts.TileSpacing, ...
    'Padding',     opts.Padding);

axs = gobjects(2,1);

% ---- (1) TFL ----
axs(1) = nexttile(tl,1);
tflPlot(W,'plotopts',struct('ShowLabels',opts.ShowLabels));
fixTiledAxesCommon(axs(1));
makeAxesDataSquare(axs(1), opts.PadFrac);
S = paperStyle();
setPanelTitle(axs(1), 'Trophic Flow Layout (TFL)', S);
%title(axs(1),'Trophic Flow Layout','FontWeight','normal');

% ---- (2) Sugiyama layered (drawn with plotTFL aesthetic) ----
axs(2) = nexttile(tl,2);

g = digraph(W);
p = plot(axs(2), g, 'Layout','layered', 'Direction','up');
Xg = p.XData(:);
Xg = max(Xg) - Xg;     % flip left-right to match TFL orientation
Yg = p.YData(:);
delete(p);

plotTFL(W, Xg, Yg, Yg, ...
    'ShowBands',  true, ...
    'ShowLabels', opts.ShowLabels);

fixTiledAxesCommon(axs(2));
makeAxesDataSquare(axs(2), opts.PadFrac);
S = paperStyle();
setPanelTitle(axs(2), 'Sugiyama layered layout', S);

% ================= optional export =================
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

out = struct();
out.fig   = fig;
out.axes  = axs;
out.files = files;

end
