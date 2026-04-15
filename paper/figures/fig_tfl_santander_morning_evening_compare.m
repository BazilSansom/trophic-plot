function out = fig_tfl_santander_morning_evening_compare(opts)
%FIG_TFL_SANTANDER_MORNING_EVENING_COMPARE
% Build weekday morning and evening TfL Santander Cycles OD networks,
% keep a common top-K station set, and export two 2-panel figures:
%
%   1) Network figure:
%      (A) Morning TFL layout
%      (B) Evening TFL layout
%
%   2) Diagnostic figure:
%      (A) Trophic level scatter: morning vs evening
%      (B) Largest station-level changes in trophic level

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
    % 2) Keep common top-K stations
    % ---------------------------------------------------------------------
    common = keepCommonTopStations_(morning, evening, opts.K);

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
    cmp.station_name_short = shortenSantanderStationLabels(cellstr(cmp.station_name));

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
    % 4) Sparse labels for network panels
    % ---------------------------------------------------------------------
    namesShort = shortenSantanderStationLabels(cellstr(T.station_name));

    labelsCommon = makeComparisonLabels_( ...
        namesShort, hM, hE, common.combined_flow, cmp.delta_h, opts);

    labelsM = labelsCommon;
    labelsE = labelsCommon;

    % ---------------------------------------------------------------------
    % 5) Figure 1: Morning vs evening networks
    % ---------------------------------------------------------------------
    figNet = figure('Color','w');
    tl1 = tiledlayout(figNet, 1, 2, 'Padding','compact', 'TileSpacing','compact');
    set(figNet,'Units','centimeters');
    figNet.Position(3) = opts.FigWidthCM;
    figNet.Position(4) = opts.FigHeightCM_Networks;

    layoutOpts = struct('HeightMode','continuous',...
        'RepelStrength', 0.6);

    %'LabelColor',       [0.6350, 0.0780, 0.1840], ... [ 0.5843 0.8157 0.9882] [0 0.7 0.7]

    plotOptsBase = struct( ...
        'ShowBands',        true, ...
        'FillSquare',       true, ...
        'ShowLabels',       opts.ShowLabels, ...
        'LabelFontSize',    opts.LabelFontSize, ...
        'LabelColor',       [0, 0, 1], ...  
        'LabelOffset',      [0 0.02], ...
        'LabelOffsetMode',  'frac', ...
        'LabelHalo',        true, ...
        'LabelHaloPadFrac', 0.006, ...
        'LabelHaloAlpha',   0.30, ...
        'LabelHaloColor',   [1 1 1], ...
        'AdjustLabels',     true, ...
        'ArrowSize',        opts.ArrowSize, ...
        'ArrowSizeMode',    'fixed', ...
        'EdgeWidthMode',    'fixed', ...
        'DownEdgeColor',    [0.8500, 0.3250, 0.0980]);

    % Morning
    ax1 = nexttile(tl1,1);
    plotOptsM = plotOptsBase;
    plotOptsM.Labels = labelsM;

    tflPlot(Wm, ...
        'Parent',     ax1, ...
        'LayoutOpts', layoutOpts, ...
        'PlotOpts',   plotOptsM);

    title(ax1, sprintf('Morning: %s, %02d:00--%02d:00  (F_0 = %.3f)', ...
        datestr(morning.target_date, 'ddd dd-mmm-yyyy'), ...
        opts.MorningWindow(1), opts.MorningWindow(2), infoM.F0_global), ...
        'Interpreter','none');
    panelTag_(ax1, 'A');
    setCleanAxes_(ax1);

    % Evening
    ax2 = nexttile(tl1,2);
    plotOptsE = plotOptsBase;
    plotOptsE.Labels = labelsE;

    tflPlot(We, ...
        'Parent',     ax2, ...
        'LayoutOpts', layoutOpts, ...
        'PlotOpts',   plotOptsE);

    title(ax2, sprintf('Evening: %s, %02d:00--%02d:00  (F_0 = %.3f)', ...
        datestr(evening.target_date, 'ddd dd-mmm-yyyy'), ...
        opts.EveningWindow(1), opts.EveningWindow(2), infoE.F0_global), ...
        'Interpreter','none');
    panelTag_(ax2, 'B');
    setCleanAxes_(ax2);

    % ---------------------------------------------------------------------
    % 6) Figure 2: Diagnostics
    % ---------------------------------------------------------------------
    figDiag = figure('Color','w');
    tl2 = tiledlayout(figDiag, 1, 2, 'Padding','compact', 'TileSpacing','compact');
    set(figDiag,'Units','centimeters');
    figDiag.Position(3) = opts.FigWidthCM;
    figDiag.Position(4) = opts.FigHeightCM_Diagnostics;

    % Scatter
    ax3 = nexttile(tl2,1);
    scatter(ax3, hM, hE, opts.ScatterSize, common.combined_flow, 'filled');
    hold(ax3,'on');
    xl = [min([hM; hE]) max([hM; hE])];
    plot(ax3, xl, xl, '--', 'LineWidth', 1.0);
    hold(ax3,'off');
    grid(ax3,'on');
    xlabel(ax3, 'Morning trophic level');
    ylabel(ax3, 'Evening trophic level');
    title(ax3, sprintf('Station trophic levels: %s', datestr(morning.target_date, 'ddd dd-mmm-yyyy')));
    colorbar(ax3);
    panelTag_(ax3, 'A');

    [~, ordMove] = sort(abs(cmp.delta_h), 'descend');
    nLab = min(opts.LabelMoversCount, height(cmp));
    for k = 1:nLab
        ii = ordMove(k);
        %text(ax3, hM(ii), hE(ii), ['  ' char(cmp.station_name(ii))], ...
        %    'FontSize', opts.MoverLabelFontSize, ...
        %    'Interpreter','none');
        text(ax3, hM(ii), hE(ii), ['  ' cmp.station_name_short{ii}], ...
            'FontSize', opts.MoverLabelFontSize, ...
            'Interpreter','none');
    end

    % Largest changes
    ax4 = nexttile(tl2,2);
    [~, ordAbs] = sort(abs(cmp.delta_h), 'descend');
    kShow = min(opts.ShowChangeCount, height(cmp));

    D = cmp(ordAbs(1:kShow), :);
    D = sortrows(D, 'delta_h', 'ascend');

    barh(ax4, D.delta_h);
    xline(ax4, 0, ':');
    yticks(ax4, 1:height(D));
    yticklabels(ax4, D.station_name_short);

    %D = cmp(ordAbs(1:kShow), :);
    %D = sortrows(D, 'delta_h', 'ascend');

    %barh(ax4, D.delta_h);
    %xline(ax4, 0, ':');
    %yticks(ax4, 1:height(D));
    %yticklabels(ax4, D.station_name);

    xlabel(ax4, '\Delta trophic level  (evening - morning)');
    title(ax4, sprintf('Largest changes (top %d by |\\Delta h|)', kShow));
    grid(ax4,'on');
    panelTag_(ax4, 'B');

    % ---------------------------------------------------------------------
    % 7) Export if requested
    % ---------------------------------------------------------------------
    files = struct( ...
        'networks_pdf',    "", ...
        'networks_png',    "", ...
        'diagnostics_pdf', "", ...
        'diagnostics_png', "" );

    if opts.Export
        if ~exist(opts.ExportDir,'dir'), mkdir(opts.ExportDir); end

        baseNet  = fullfile(opts.ExportDir, opts.ExportBaseNetworks);
        baseDiag = fullfile(opts.ExportDir, opts.ExportBaseDiagnostics);

        if opts.ExportPDF
            pdfNet = string(baseNet) + ".pdf";
            pdfDiag = string(baseDiag) + ".pdf";

            if opts.PDFVector
                exportgraphics(figNet,  pdfNet,  'ContentType','vector');
                exportgraphics(figDiag, pdfDiag, 'ContentType','vector');
            else
                exportgraphics(figNet,  pdfNet,  'ContentType','image', 'Resolution', opts.PNGResolution);
                exportgraphics(figDiag, pdfDiag, 'ContentType','image', 'Resolution', opts.PNGResolution);
            end

            files.networks_pdf = pdfNet;
            files.diagnostics_pdf = pdfDiag;
        end

        if opts.ExportPNG
            pngNet = string(baseNet) + ".png";
            pngDiag = string(baseDiag) + ".png";

            exportgraphics(figNet,  pngNet,  'Resolution', opts.PNGResolution);
            exportgraphics(figDiag, pngDiag, 'Resolution', opts.PNGResolution);

            files.networks_png = pngNet;
            files.diagnostics_png = pngDiag;
        end
    end

    % ---------------------------------------------------------------------
    % 8) Return
    % ---------------------------------------------------------------------
    out = struct();
    out.fig_networks    = figNet;
    out.fig_diagnostics = figDiag;
    out.files = files;

    out.morning = morning;
    out.evening = evening;
    out.source_morning = morning.source;
    out.source_evening = evening.source;

    out.common  = common;
    out.compare = cmp;

    out.h_morning    = hM;
    out.h_evening    = hE;
    out.info_morning = infoM;
    out.info_evening = infoE;
end

% =========================================================================
function opts = applyDefaults_(opts)
    thisDir  = fileparts(mfilename('fullpath'));
    repoRoot = getRepoRoot_(thisDir);

    d = struct();

    d.URL = "";
    %d.DataDir = fullfile(repoRoot, 'data', 'tfl_santander');
    d.DataDir = fullfile(repoRoot, 'paper', 'data', 'tfl_santander');
    d.Filename = "";
    d.OverwriteDownload = false;

    d.TargetDate = [];
    d.MorningWindow = [7 10];
    d.EveningWindow = [17 20];
    d.RemoveSelfLoops = true;

    d.K = 30;

    d.ShowLabels = true;
    d.LabelCount = 4;
    d.LabelFontSize = 6;
    d.NExtremes = 3;
    d.NMovers   = 4;
    d.NFlow     = 4;
    d.NMiddle   = 4;

    d.ArrowSize = 0.03;

    d.ScatterSize = 40;
    d.LabelMoversCount = 6;
    d.MoverLabelFontSize = 7;
    d.ShowChangeCount = 10;

    d.FigWidthCM = 22;
    d.FigHeightCM_Networks = 9;
    d.FigHeightCM_Diagnostics = 11;

    d.Export = false;
    d.ExportDir = '.';
    d.ExportBaseNetworks = 'fig_tfl_santander_morning_evening_networks';
    d.ExportBaseDiagnostics = 'fig_tfl_santander_morning_evening_diagnostics';
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
function common = keepCommonTopStations_(morning, evening, K)
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

    [Wm, Tsub] = induceByIDs_(morning, keepIDs);
    [We, ~]    = induceByIDs_(evening, keepIDs);

    common = struct();
    common.station_table = Tsub;
    common.W_morning = Wm;
    common.W_evening = We;
    common.combined_flow = combined(ord(1:min(K, numel(ord))));
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

    if opts.NMiddle > 0
        hMid = 0.5 * (hM + hE);
        [~, ordMid] = sort(hMid, 'ascend');
        midPos = round(linspace(max(2, floor(n*0.30)), min(n-1, ceil(n*0.70)), opts.NMiddle));
        keep = [keep; ordMid(midPos(:))];
    end

    keep = unique(keep, 'stable');

    for i = 1:numel(keep)
        labels{keep(i)} = names{keep(i)};
    end
end

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

% =========================================================================
function b = makeBuildOpts_(opts, targetDate, timeWindow, overwriteDownload)
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

% =========================================================================
function repoRoot = getRepoRoot_(startDir)
    repoRoot = startDir;
    while true
        if exist(fullfile(repoRoot, '.git'), 'dir') || ...
           (exist(fullfile(repoRoot, 'src'), 'dir') && exist(fullfile(repoRoot, 'data'), 'dir'))
            return;
        end

        parent = fileparts(repoRoot);
        if strcmp(parent, repoRoot)
            error('Could not locate repo root starting from: %s', startDir);
        end
        repoRoot = parent;
    end
end