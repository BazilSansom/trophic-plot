function out = fig_topology_vs_weighted_flow(opts)
%FIG_TOPOLOGY_VS_WEIGHTED_FLOW  2x2 figure: columns = method, rows = weights.
%
% Layout:
%   (1,1) TFL, weighted         (1,2) Sugiyama layered, weighted
%   (2,1) TFL, uniform          (2,2) Sugiyama layered, uniform
%
% Export options:
%   opts.Export        : true/false (default false)
%   opts.ExportBase    : base filename (default 'fig_topology_vs_weighted_flow')
%   opts.ExportDir     : output folder (default pwd)
%   opts.ExportPDF     : true/false (default true)
%   opts.ExportPNG     : true/false (default true)
%   opts.PNGResolution : dpi (default 300)
%   opts.PDFVector     : true/false (default true) -> ContentType 'vector'

% ----------------------------- Parse options ---------------------------------
if nargin < 1 || isempty(opts), opts = struct(); end

if ~isfield(opts,'ShowLabels');      opts.ShowLabels = true; end
%if ~isfield(opts,'ArrowSize');      opts.ArrowSize = 0.10; end
if ~isfield(opts,'PadFrac');         opts.PadFrac    = 0.02; end
if ~isfield(opts,'TileSpacing');     opts.TileSpacing = 'compact'; end
if ~isfield(opts,'Padding');         opts.Padding     = 'compact'; end
if ~isfield(opts,'ShowRowLabels');   opts.ShowRowLabels = true; end

% ---- export defaults ----
if ~isfield(opts,'Export');         opts.Export = false; end
if ~isfield(opts,'ExportBase');     opts.ExportBase = 'fig_topology_vs_weighted_flow'; end
if ~isfield(opts,'ExportDir');      opts.ExportDir  = pwd; end
if ~isfield(opts,'ExportPDF');      opts.ExportPDF  = true; end
if ~isfield(opts,'ExportPNG');      opts.ExportPNG  = true; end
if ~isfield(opts,'PNGResolution');  opts.PNGResolution = 300; end
if ~isfield(opts,'PDFVector');      opts.PDFVector = true; end

S = paperStyle();

% ------------------------------ Data -----------------------------------------
E = [3 2   0.1
     4 2   0.1
     2 3 200.0
     8 4 200.0
     1 5 200.0
     7 5 200.0
     3 6 200.0
     5 6 200.0
     6 7   0.1
     1 8 200.0
     4 8   0.1];

w_weighted = edgelistToAdjacency(E);
w_uniform  = spones(w_weighted);  % same topology, uniform weights

% ------------------------------ Figure ---------------------------------------
fig = figure('Color','w');

tl = tiledlayout(fig,2,2, ...
    'TileSpacing', opts.TileSpacing, ...
    'Padding',     opts.Padding);

axs = gobjects(4,1);

% Make figure roughly square (helps with square axes)
set(fig,'Units','centimeters');
pos  = fig.Position;
side = min(pos(3), pos(4));
fig.Position = [pos(1) pos(2) side side];

% ------------------------------ Panels ---------------------------------------
% (1) TFL weighted
axs(1) = nexttile(tl,1);
tflPlot(w_weighted, 'plotopts', struct('ShowLabels', opts.ShowLabels));
setPanelTitle(axs(1), "Trophic Flow Layout (TFL)",S);
fixTiledAxesCommon(axs(1));
makeAxesDataSquare(axs(1), opts.PadFrac);

% (2) Layered weighted (drawn using plotTFL for aesthetic consistency)
axs(2) = nexttile(tl,2);
g = digraph(w_weighted);
p = plot(axs(2), g, 'Layout','layered', 'Direction','up');  % compute positions
Xg = p.XData(:); Yg = p.YData(:);
delete(p);

plotTFL(w_weighted, Xg, Yg, Yg, ...
    'Parent',     axs(2), ...
    'ClearAxes',  true, ...
    'ShowBands',  false, ...
    'ShowLabels', opts.ShowLabels);
setPanelTitle(axs(2), "Sugiyama layered", S);
fixTiledAxesCommon(axs(2));
makeAxesDataSquare(axs(2), opts.PadFrac);

% (3) TFL uniform
axs(3) = nexttile(tl,3);
tflPlot(w_uniform, 'plotopts', struct('ShowLabels', opts.ShowLabels));
fixTiledAxesCommon(axs(3));
makeAxesDataSquare(axs(3), opts.PadFrac);

% (4) Layered uniform
axs(4) = nexttile(tl,4);
g = digraph(w_uniform);
p = plot(axs(4), g, 'Layout','layered', 'Direction','up');
Xg = p.XData(:); Yg = p.YData(:);
delete(p);

plotTFL(w_uniform, Xg, Yg, Yg, ...
    'Parent',     axs(4), ...
    'ClearAxes',  true, ...
    'ShowBands',  false, ...
    'ShowLabels', opts.ShowLabels);
fixTiledAxesCommon(axs(4));
makeAxesDataSquare(axs(4), opts.PadFrac);

% ------------------------- Row labels (optional) ------------------------------
if opts.ShowRowLabels
    addRowLabel_(axs(1), 'Weighted');
    addRowLabel_(axs(3), 'Uniform weights');
end

drawnow;

% ------------------------------ Export ---------------------------------------
files = struct('pdf',"",'png',"");
if opts.Export
    if ~exist(opts.ExportDir,'dir'), mkdir(opts.ExportDir); end
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

% ------------------------------ Outputs --------------------------------------
out = struct();
out.fig   = fig;
out.axes  = axs;
out.files = files;

end

% =====================================================================
function addRowLabel_(ax, txt)
% Put a row header just to the left of the tile, without changing layout.
text(ax, -0.10, 0.50, txt, ...
    'Units','normalized', ...
    'Rotation', 90, ...
    'FontWeight','bold', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'Clipping','off');
end
