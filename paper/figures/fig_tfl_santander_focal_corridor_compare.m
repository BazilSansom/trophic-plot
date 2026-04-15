function out = fig_tfl_santander_focal_corridor_compare(opts)
%FIG_TFL_SANTANDER_FOCAL_CORRIDOR_COMPARE
% Compare trophically aligned focal-node corridors for a TfL Santander Cycles
% station in morning and evening windows.
%
% Assumes helpers already exist on path:
%   - build_tfl_santander_week_window
%   - trophic_levels
%   - tflPlot
%   - plotTFL_focal_corridor
%   - tfl_focal_corridor_mask
%
% Usage:
%   out = fig_tfl_santander_focal_corridor_compare();
%
%   out = fig_tfl_santander_focal_corridor_compare(struct( ...
%       'TargetDate', datetime(2025,6,2), ...
%       'MorningWindow', [7 10], ...
%       'EveningWindow', [17 20], ...
%       'K', 30, ...
%       'FocalStationName', "Waterloo Station 3, Waterloo", ...
%       'MaxHops', 2, ...
%       'Export', true));
%
% Outputs:
%   out.fig
%   out.files
%   out.morning
%   out.evening
%   out.common
%   out.compare
%   out.focal
%   out.corridor_morning
%   out.corridor_evening

    if nargin < 1 || isempty(opts), opts = struct(); end
    opts = applyDefaults_(opts);


    % ---------------------------------------------------------------------
    % 1) Build morning and evening datasets
    % ---------------------------------------------------------------------
    buildOptsMorning = makeBuildOpts_(opts, opts.TargetDate, opts.MorningWindow, opts.OverwriteDownload);
    morning = build_tfl_santander_week_window(buildOptsMorning);

    % Use resolved morning.target_date so evening is guaranteed same day
    buildOptsEvening = makeBuildOpts_(opts, morning.target_date, opts.EveningWindow, false);
    evening = build_tfl_santander_week_window(buildOptsEvening);

    % ---------------------------------------------------------------------
    % 2) Resolve focal station on full data, then keep common top-K stations
    % ---------------------------------------------------------------------
    focal = resolveFocalStation_(morning, evening, opts.FocalStationName);

    common = keepCommonTopStations_(morning, evening, opts.K, focal.station_id);

    focalIdx = find(common.station_table.station_id == focal.station_id, 1, 'first');
    if isempty(focalIdx)
        error('Focal station was not retained in common station set.');
    end

    Wm = full(common.W_morning);
    We = full(common.W_evening);
    T  = common.station_table;

    % ---------------------------------------------------------------------
    % 3) Trophic levels and comparison table
    % ---------------------------------------------------------------------
    [hM, infoM] = trophic_levels(Wm, 'ComputeCoherence', true);
    [hE, infoE] = trophic_levels(We, 'ComputeCoherence', true);

    hM = hM(:);
    hE = hE(:);

    cmp = T;
    cmp.h_morning = hM;
    cmp.h_evening = hE;
    cmp.delta_h   = hE - hM;

    inM  = full(sum(Wm,1))';
    outM = full(sum(Wm,2));
    totM = inM + outM;
    rhoM = nan(size(totM));
    okM = totM > 0;
    rhoM(okM) = (inM(okM) - outM(okM)) ./ totM(okM);

    inE  = full(sum(We,1))';
    outE = full(sum(We,2));
    totE = inE + outE;
    rhoE = nan(size(totE));
    okE = totE > 0;
    rhoE(okE) = (inE(okE) - outE(okE)) ./ totE(okE);

    cmp.in_morning   = inM;
    cmp.out_morning  = outM;
    cmp.rho_morning  = rhoM;
    cmp.in_evening   = inE;
    cmp.out_evening  = outE;
    cmp.rho_evening  = rhoE;

    % ---------------------------------------------------------------------
    % 4) Labels
    % ---------------------------------------------------------------------
    namesFull  = cellstr(T.station_name);
    %namesShort = maybeShortenLabels_(namesFull);
    namesShort = shortenSantanderStationLabels(namesFull);

    labelsCommon = makeComparisonLabels_( ...
        namesShort, hM, hE, common.combined_flow, cmp.delta_h, opts);

    % Always label focal node
    labelsCommon{focalIdx} = namesShort{focalIdx};

    % ---------------------------------------------------------------------
    % 5) Compute TFL coordinates on hidden axes
    % ---------------------------------------------------------------------
    layoutOpts = struct('HeightMode', 'continuous');

    %[Xm, Ym] = computeTFLcoords_(Wm, layoutOpts);
    %[Xe, Ye] = computeTFLcoords_(We, layoutOpts);

    [Xm, Ym, hPlotM] = computeTFLcoords_(Wm, layoutOpts);
    [Xe, Ye, hPlotE] = computeTFLcoords_(We, layoutOpts);

    % ---------------------------------------------------------------------
    % 6) Plot corridor comparison
    % ---------------------------------------------------------------------
    fig = figure('Color','w');
    tl = tiledlayout(fig, 1, 2, 'Padding','compact', 'TileSpacing','compact');

    set(fig,'Units','centimeters');
    fig.Position(3) = opts.FigWidthCM;
    fig.Position(4) = opts.FigHeightCM;

    corrOpts = struct( ...
        'MaxHops',   opts.MaxHops, ...
        'MinDeltaH', opts.MinDeltaH, ...
        'MinWeight', opts.MinWeight);

    % Morning
    ax1 = nexttile(tl,1);
    %outCorrM = plotTFL_focal_corridor( ...
    %    Wm, Xm, Ym, hM, focalIdx, labelsCommon, corrOpts, struct('Parent', ax1));
    outCorrM = plotTFL_focal_corridor( ...
        Wm, Xm, Ym, hPlotM, focalIdx, labelsCommon, corrOpts, struct('Parent', ax1));
    %title(ax1, sprintf('%s: morning %02d:00--%02d:00  (F_0 = %.3f)', ...
    %    namesShort{focalIdx}, opts.MorningWindow(1), opts.MorningWindow(2), infoM.F0_global), ...
    %    'Interpreter','none');
    title(ax1, sprintf('%s: %s, %02d:00--%02d:00  (F_0 = %.3f)', ...
        namesShort{focalIdx}, ...
        datestr(morning.target_date, 'ddd dd-mmm-yyyy'), ...
        opts.MorningWindow(1), opts.MorningWindow(2), infoM.F0_global), ...
        'Interpreter','none');
    panelTag_(ax1, 'A');
    setCleanAxes_(ax1);

    % Evening
    ax2 = nexttile(tl,2);
    %outCorrE = plotTFL_focal_corridor( ...
    %    We, Xe, Ye, hE, focalIdx, labelsCommon, corrOpts, struct('Parent', ax2));
    outCorrE = plotTFL_focal_corridor( ...
        We, Xe, Ye, hPlotE, focalIdx, labelsCommon, corrOpts, struct('Parent', ax2));
    %title(ax2, sprintf('%s: evening %02d:00--%02d:00  (F_0 = %.3f)', ...
    %    namesShort{focalIdx}, opts.EveningWindow(1), opts.EveningWindow(2), infoE.F0_global), ...
    %    'Interpreter','none');
    title(ax2, sprintf('%s: %s, %02d:00--%02d:00  (F_0 = %.3f)', ...
        namesShort{focalIdx}, ...
        datestr(evening.target_date, 'ddd dd-mmm-yyyy'), ...
        opts.EveningWindow(1), opts.EveningWindow(2), infoE.F0_global), ...
        'Interpreter','none');
    panelTag_(ax2, 'B');
    setCleanAxes_(ax2);

    % ---------------------------------------------------------------------
    % 7) Export
    % ---------------------------------------------------------------------
    files = struct('pdf',"",'png',"");
    if opts.Export
        if ~exist(opts.ExportDir,'dir'), mkdir(opts.ExportDir); end
        base = fullfile(opts.ExportDir, opts.ExportBase);

        if opts.ExportPDF
            pdfFile = string(base) + ".pdf";
            if opts.PDFVector
                exportgraphics(fig, pdfFile, 'ContentType','vector');
            else
                exportgraphics(fig, pdfFile, 'ContentType','image', ...
                    'Resolution', opts.PNGResolution);
            end
            files.pdf = pdfFile;
        end

        if opts.ExportPNG
            pngFile = string(base) + ".png";
            exportgraphics(fig, pngFile, 'Resolution', opts.PNGResolution);
            files.png = pngFile;
        end
    end

    % ---------------------------------------------------------------------
    % 8) Return
    % ---------------------------------------------------------------------
    out = struct();
    out.fig   = fig;
    out.files = files;

    out.morning = morning;
    out.evening = evening;
    out.common  = common;
    out.compare = cmp;

    out.focal = struct( ...
        'station_id', focal.station_id, ...
        'name_full',  string(T.station_name(focalIdx)), ...
        'name_short', string(namesShort{focalIdx}), ...
        'idx_common', focalIdx);

    out.h_morning    = hM;
    out.h_evening    = hE;
    out.info_morning = infoM;
    out.info_evening = infoE;

    out.corridor_morning = outCorrM;
    out.corridor_evening = outCorrE;
end

% =========================================================================
function opts = applyDefaults_(opts)
    d = struct();

    %d.URL = "https://cycling.data.tfl.gov.uk/usage-stats/421JourneyDataExtract01Jun2025-15Jun2025.csv";
    %d.DataDir = fullfile(pwd, 'data', 'tfl_santander');
    %d.Filename = "421JourneyDataExtract01Jun2025-15Jun2025.csv";

    thisFile = mfilename('fullpath');
    thisDir  = fileparts(thisFile);
    repoRoot = getRepoRoot_(thisDir);

    d = struct();
    d.URL = "";
    %d.DataDir = fullfile(repoRoot, 'data', 'tfl_santander');
    d.DataDir = fullfile(repoRoot, 'paper', 'data', 'tfl_santander');
    %d.DataDir = fullfile(pwd, 'data', 'tfl_santander');
    d.Filename = "";
    d.OverwriteDownload = false;

    d.TargetDate = [];
    d.MorningWindow = [7 10];
    d.EveningWindow = [17 20];
    d.RemoveSelfLoops = true;

    d.K = 30;
    d.FocalStationName = "Waterloo Station 3, Waterloo";

    d.MaxHops   = 2;
    d.MinDeltaH = 1e-6;
    d.MinWeight = 0;

    d.NExtremes = 3;
    d.NMovers   = 4;
    d.NFlow     = 5;
    d.NMiddle   = 3;

    d.FigWidthCM  = 22;
    d.FigHeightCM = 9.5;

    d.Export = false;
    d.ExportDir = '.';
    d.ExportBase = 'fig_tfl_santander_focal_corridor_compare';
    d.ExportPDF = true;
    d.ExportPNG = true;
    d.PNGResolution = 300;
    d.PDFVector = true;

    fn = fieldnames(d);
    for k = 1:numel(fn)
        f = fn{k};
        if ~isfield(opts,f) || isempty(opts.(f))
            opts.(f) = d.(f);
        end
    end
end

% =========================================================================
function focal = resolveFocalStation_(morning, evening, focalName)
% Resolve focal station using exact name first, then contains-match.
% Must exist in both windows.

    Tm = morning.station_table;
    Te = evening.station_table;

    nmM = string(Tm.station_name);
    nmE = string(Te.station_name);

    exactM = find(strcmpi(nmM, string(focalName)), 1, 'first');
    exactE = find(strcmpi(nmE, string(focalName)), 1, 'first');

    if ~isempty(exactM) && ~isempty(exactE)
        if Tm.station_id(exactM) ~= Te.station_id(exactE)
            error('Exact focal-name match found in both windows but IDs differ.');
        end
        focal = struct('station_id', Tm.station_id(exactM), 'name', nmM(exactM));
        return;
    end

    hitM = find(contains(lower(nmM), lower(string(focalName))), 1, 'first');
    hitE = find(contains(lower(nmE), lower(string(focalName))), 1, 'first');

    if ~isempty(hitM) && ~isempty(hitE)
        if Tm.station_id(hitM) ~= Te.station_id(hitE)
            error(['Contains-match for focal station is ambiguous across windows. ' ...
                   'Please use exact full station name.']);
        end
        focal = struct('station_id', Tm.station_id(hitM), 'name', nmM(hitM));
        return;
    end

    error('Could not resolve focal station "%s" in both morning and evening datasets.', string(focalName));
end

% =========================================================================
function common = keepCommonTopStations_(morning, evening, K, forcedIDs)
% Keep same top-K stations in both windows, based on combined flow and
% restricted to stations active in both windows. Forced station IDs are
% always retained.

    if nargin < 4, forcedIDs = []; end
    forcedIDs = unique(double(forcedIDs(:)));

    Tm = morning.station_table;
    Te = evening.station_table;

    flowM = Tm.in_trips + Tm.out_trips;
    flowE = Te.in_trips + Te.out_trips;

    idsM = Tm.station_id;
    idsE = Te.station_id;
    commonIDs = intersect(idsM, idsE);

    n = numel(commonIDs);
    fM = zeros(n,1);
    fE = zeros(n,1);

    mapM = containers.Map(double(idsM), double(1:height(Tm)));
    mapE = containers.Map(double(idsE), double(1:height(Te)));

    for k = 1:n
        sid = commonIDs(k);
        iM = mapM(sid);
        iE = mapE(sid);
        fM(k) = flowM(iM);
        fE(k) = flowE(iE);
    end

    combined = fM + fE;
    [~, ord] = sort(combined, 'descend');
    keepIDs = commonIDs(ord(1:min(K, numel(ord))));

    % Force focal station(s) in if they are common to both windows
    forcedIDs = intersect(forcedIDs, commonIDs);
    keepIDs = unique([keepIDs(:); forcedIDs(:)], 'stable');

    [Wm, Tsub] = induceByIDs_(morning, keepIDs);
    [We, ~]    = induceByIDs_(evening, keepIDs);

    % Recompute combined flow in retained order
    combinedFlow = zeros(numel(keepIDs),1);
    for k = 1:numel(keepIDs)
        sid = keepIDs(k);
        i = find(commonIDs == sid, 1, 'first');
        combinedFlow(k) = combined(i);
    end

    common = struct();
    common.station_table  = Tsub;
    common.W_morning      = Wm;
    common.W_evening      = We;
    common.combined_flow  = combinedFlow;
    common.station_ids    = keepIDs;
end

% =========================================================================
function [Wsub, Tsub] = induceByIDs_(out, keepIDs)
    T = out.station_table;
    W = out.W;

    ids = T.station_id;
    map = containers.Map(double(ids), double(1:height(T)));

    idx = zeros(numel(keepIDs),1);
    for k = 1:numel(keepIDs)
        idx(k) = map(keepIDs(k));
    end

    Wsub = W(idx, idx);
    Tsub = T(idx, :);
    Tsub.idx = (1:height(Tsub))';
end

% =========================================================================
function labels = makeComparisonLabels_(names, hM, hE, flow, deltaH, opts)
% Common informative label set for both panels.

    n = numel(names);
    labels = repmat({''}, n, 1);
    keep = [];

    [~, ordMAsc]  = sort(hM, 'ascend');
    [~, ordMDesc] = sort(hM, 'descend');
    keep = [keep; ordMAsc(1:min(opts.NExtremes,n)); ordMDesc(1:min(opts.NExtremes,n))];

    [~, ordEAsc]  = sort(hE, 'ascend');
    [~, ordEDesc] = sort(hE, 'descend');
    keep = [keep; ordEAsc(1:min(opts.NExtremes,n)); ordEDesc(1:min(opts.NExtremes,n))];

    [~, ordD] = sort(abs(deltaH), 'descend');
    keep = [keep; ordD(1:min(opts.NMovers,n))];

    [~, ordF] = sort(flow, 'descend');
    keep = [keep; ordF(1:min(opts.NFlow,n))];

    hMid = 0.5 * (hM + hE);
    [~, ordMid] = sort(hMid, 'ascend');
    if opts.NMiddle > 0
        midPos = round(linspace(max(2, floor(n*0.30)), min(n-1, ceil(n*0.70)), opts.NMiddle));
        keep = [keep; ordMid(midPos(:))];
    end

    keep = unique(keep, 'stable');

    for ii = 1:numel(keep)
        labels{keep(ii)} = names{keep(ii)};
    end
end

% =========================================================================
function [X, Y, h, info] = computeTFLcoords_(W, layoutOpts)
% Compute RAW TFL layout coordinates directly from trophicLayoutMulti.
% Do not use tflPlot here, because tflPlot now returns rendered coordinates.

    if nargin < 2 || isempty(layoutOpts)
        layoutOpts = struct();
    end

    nv = structToNVLocal_(layoutOpts);
    [X, Y, ~, info] = trophicLayoutMulti(W, nv{:});

    if isfield(info, 'h') && ~isempty(info.h)
        h = info.h(:);
    else
        h = trophic_levels(W);
        h = h(:);
    end
end

function nv = structToNVLocal_(s)
    nv = {};
    fn = fieldnames(s);
    for k = 1:numel(fn)
        nv(end+1:end+2) = {fn{k}, s.(fn{k})}; %#ok<AGROW>
    end
end

%{
function [X, Y] = computeTFLcoords_(W, layoutOpts)
% Compute TFL coordinates on hidden axes using tflPlot outputs.

    tmpFig = figure('Visible','off', 'Color','w');
    tmpAx = axes('Parent', tmpFig); %#ok<LAXES>

    plotOpts = struct( ...
        'Parent', tmpAx, ...
        'ClearAxes', true, ...
        'ShowLabels', false, ...
        'ShowBands', false, ...
        'FillSquare', true);

    try
        [X, Y] = tflPlot(W, ...
            'LayoutOpts', layoutOpts, ...
            'PlotOpts',   plotOpts);
    catch ME
        close(tmpFig);
        rethrow(ME);
    end

    if isvalid(tmpFig)
        close(tmpFig);
    end
end
%}

% =========================================================================
function setCleanAxes_(ax)
    set(ax,'XTick',[],'YTick',[]);
    box(ax,'on');
end

% =========================================================================
function panelTag_(ax, tag)
    text(ax, 0.02, 0.98, tag, ...
        'Units','normalized', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','top', ...
        'FontWeight','bold', ...
        'Clipping','off');
end

function b = makeBuildOpts_(opts, targetDate, timeWindow, overwriteDownload)
%MAKEBUILDOPTS_  Build clean opts struct for build_tfl_santander_week_window.
%
% Pass URL/Filename only if explicitly supplied, so the builder can
% resolve the source file from TargetDate when those overrides are empty.

    b = struct();
    b.DataDir = opts.DataDir;
    b.OverwriteDownload = overwriteDownload;
    b.TargetDate = targetDate;
    b.TimeWindow = timeWindow;
    b.RemoveSelfLoops = opts.RemoveSelfLoops;

    if isfield(opts, 'URL') && strlength(string(opts.URL)) > 0
        b.URL = opts.URL;
    end

    if isfield(opts, 'Filename') && strlength(string(opts.Filename)) > 0
        b.Filename = opts.Filename;
    end
end

function repoRoot = getRepoRoot_(startDir)
    repoRoot = startDir;
    while true
        if exist(fullfile(repoRoot, '.git'), 'dir') || ...
           exist(fullfile(repoRoot, 'data'), 'dir') || ...
           exist(fullfile(repoRoot, 'src'), 'dir')
            return;
        end

        parent = fileparts(repoRoot);
        if strcmp(parent, repoRoot)
            error('Could not locate repo root starting from: %s', startDir);
        end
        repoRoot = parent;
    end
end