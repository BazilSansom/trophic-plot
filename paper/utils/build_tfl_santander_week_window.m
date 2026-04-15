function out = build_tfl_santander_week_window(opts)
%BUILD_TFL_SANTANDER_WEEK_WINDOW
% Load one raw TfL Santander Cycles weekly CSV, standardise fields,
% filter to a single weekday time window, and build a station-level OD matrix.
%
% Preferred usage:
%   out = build_tfl_santander_week_window(struct( ...
%       'TargetDate', datetime(2025,6,2), ...
%       'TimeWindow', [7 10], ...
%       'DataDir', fullfile(repoRoot,'data','tfl_santander')));
%
% Optional overrides:
%   'Filename'  : force a specific local CSV filename in DataDir
%   'URL'       : download this URL if the chosen file is missing
%
% Outputs:
%   out.filepath
%   out.source
%   out.raw
%   out.trips
%   out.trips_window
%   out.target_date
%   out.station_table
%   out.W
%   out.W_duration
%   out.summary

    if nargin < 1 || isempty(opts), opts = struct(); end
    opts = applyDefaults_(opts);

    if ~exist(opts.DataDir, 'dir')
        mkdir(opts.DataDir);
    end

    % ---------------------------------------------------------------------
    % 1) Resolve source file from TargetDate / Filename / URL
    % ---------------------------------------------------------------------
    src = resolveTflSourceForDate_(opts);
    localFile = src.localFile;

    % Download if needed and possible
    if ~exist(localFile, 'file')
    if strlength(string(src.url)) > 0
        fprintf('Downloading TfL file (%s)...\n', src.mode);
        websave(localFile, src.url);
    else
        error(['Resolved source file does not exist locally and no URL was available.\n' ...
               'Missing file: %s'], localFile);
    end
    else
        fprintf('Using existing local file (%s): %s\n', src.mode, localFile);
    end

    %{
    if ~exist(localFile, 'file')
        if strlength(string(src.url)) > 0
            fprintf('Downloading TfL file...\n');
            websave(localFile, src.url);
        else
            error(['Resolved source file does not exist locally and no URL was provided.\n' ...
                   'Missing file: %s'], localFile);
        end
    else
        fprintf('Using existing local file: %s\n', localFile);
    end
    %}

    % ---------------------------------------------------------------------
    % 2) Read raw CSV
    % ---------------------------------------------------------------------
    Traw = readTflCsv_(localFile);

    % ---------------------------------------------------------------------
    % 3) Standardise fields
    % ---------------------------------------------------------------------
    T = standardiseTflTrips_(Traw);

    % ---------------------------------------------------------------------
    % 4) Choose / validate target weekday
    % ---------------------------------------------------------------------
    allDates = unique(dateshift(T.start_date, 'start', 'day'));
    wd = weekday(allDates); % 1=Sun, 2=Mon, ..., 7=Sat
    isWeekday = ismember(wd, 2:6);

    if isempty(opts.TargetDate)
        if ~any(isWeekday)
            error('No weekday dates found in file: %s', localFile);
        end
        targetDate = allDates(find(isWeekday, 1, 'first'));
    else
        targetDate = dateshift(opts.TargetDate, 'start', 'day');

        if ~ismember(targetDate, allDates)
            minD = min(allDates);
            maxD = max(allDates);
            error(['Requested TargetDate %s is not present in source file.\n' ...
                   'File: %s\n' ...
                   'Observed date span in file: %s to %s'], ...
                   datestr(targetDate), localFile, datestr(minD), datestr(maxD));
        end

        if ~ismember(weekday(targetDate), 2:6)
            warning('Requested TargetDate %s is not a weekday.', datestr(targetDate));
        end
    end

    % ---------------------------------------------------------------------
    % 5) Filter to single weekday time window
    % ---------------------------------------------------------------------
    h0 = opts.TimeWindow(1);
    h1 = opts.TimeWindow(2);

    tt = timeofday(T.start_date);
    dayMask = dateshift(T.start_date, 'start', 'day') == targetDate;
    timeMask = tt >= hours(h0) & tt < hours(h1);

    Tw = T(dayMask & timeMask, :);

    if isempty(Tw)
        warning('No trips found for %s between %02d:00 and %02d:00.', ...
            datestr(targetDate), h0, h1);
    end

    % ---------------------------------------------------------------------
    % 6) Remove rows with missing station IDs, optionally remove self-loops
    % ---------------------------------------------------------------------
    ok = isfinite(Tw.start_station_id) & isfinite(Tw.end_station_id);
    Tw = Tw(ok, :);

    if opts.RemoveSelfLoops
        Tw = Tw(Tw.start_station_id ~= Tw.end_station_id, :);
    end

    % ---------------------------------------------------------------------
    % 7) Build station lookup and OD matrices
    % ---------------------------------------------------------------------
    [stationTable, W, Wdur] = buildStationOd_(Tw);

    % ---------------------------------------------------------------------
    % 8) Summary
    % ---------------------------------------------------------------------
    summary = struct();
    summary.target_date = targetDate;
    summary.time_window = opts.TimeWindow;
    summary.n_trips_total = height(T);
    summary.n_trips_window = height(Tw);
    summary.n_stations = height(stationTable);
    summary.n_edges = nnz(W);
    summary.total_weight = full(sum(W(:)));
    summary.total_duration_sec = full(sum(Wdur(:)));

    % ---------------------------------------------------------------------
    % 9) Return
    % ---------------------------------------------------------------------
    out = struct();
    out.filepath      = localFile;
    out.source        = src;
    out.raw           = Traw;
    out.trips         = T;
    out.trips_window  = Tw;
    out.target_date   = targetDate;
    out.station_table = stationTable;
    out.W             = W;
    out.W_duration    = Wdur;
    out.summary       = summary;
end

% =========================================================================
function opts = applyDefaults_(opts)
    thisDir  = fileparts(mfilename('fullpath'));
    repoRoot = getRepoRoot_(thisDir);

    d = struct();

    % Optional URL. If omitted, function only searches DataDir locally.
    d.URL = "";

    % Directory holding TfL weekly CSVs
    d.DataDir = fullfile(repoRoot, 'paper', 'data', 'tfl_santander');

    % Optional explicit filename override
    d.Filename = "";

    d.OverwriteDownload = false;

    % If empty, choose first weekday present in chosen file
    d.TargetDate = [];

    % [start_hour end_hour), e.g. [7 10] = 07:00 to 09:59:59
    d.TimeWindow = [7 10];

    d.RemoveSelfLoops = true;

    fn = fieldnames(d);
    for k = 1:numel(fn)
        f = fn{k};
        if ~isfield(opts, f) || isempty(opts.(f))
            opts.(f) = d.(f);
        end
    end
end

%{
function opts = applyDefaults_(opts)
    d = struct();

    % Optional URL. If omitted, function only searches DataDir locally.
    d.URL = "";

    % Directory holding TfL weekly CSVs
    %d.DataDir = fullfile(pwd, 'data', 'tfl_santander');
    d.DataDir = fullfile(repoRoot, 'paper', 'data', 'tfl_santander');

    % Optional explicit filename override
    d.Filename = "";

    d.OverwriteDownload = false;

    % If empty, choose first weekday present in chosen file
    d.TargetDate = [];

    % [start_hour end_hour), e.g. [7 10] = 07:00 to 09:59:59
    d.TimeWindow = [7 10];

    d.RemoveSelfLoops = true;

    fn = fieldnames(d);
    for k = 1:numel(fn)
        f = fn{k};
        if ~isfield(opts, f) || isempty(opts.(f))
            opts.(f) = d.(f);
        end
    end
end
%}


% =========================================================================
function src = resolveTflSourceForDate_(opts)
% Resolve local/remote TfL source file from:
%   1) explicit Filename, else
%   2) local DataDir scan using TargetDate, else
%   3) explicit URL, else
%   4) remote discovery from TargetDate via discoverRemoteTflJourneyFile_

    src = struct();
    src.mode = "";
    src.filename = "";
    src.localFile = "";
    src.url = string(opts.URL);
    src.range_start = NaT;
    src.range_end   = NaT;

    % ------------------------------------------------------------------
    % Case 1: explicit filename override
    % ------------------------------------------------------------------
    if strlength(string(opts.Filename)) > 0
        src.mode = "explicit_filename";
        src.filename = string(opts.Filename);
        src.localFile = fullfile(opts.DataDir, char(src.filename));

        [d0, d1, ok] = parseTflDateRangeFromFilename_(src.filename);
        if ok
            src.range_start = d0;
            src.range_end   = d1;
        end
        return;
    end

    % ------------------------------------------------------------------
    % Case 2: no TargetDate
    % ------------------------------------------------------------------
    if isempty(opts.TargetDate)
        % Prefer URL-derived filename if URL supplied
        if strlength(string(opts.URL)) > 0
            fname = filepartsFromUrl_(opts.URL);
            src.mode = "url_default";
            src.filename = string(fname);
            src.localFile = fullfile(opts.DataDir, fname);

            [d0, d1, ok] = parseTflDateRangeFromFilename_(src.filename);
            if ok
                src.range_start = d0;
                src.range_end   = d1;
            end
            return;
        end

        % Otherwise use first matching local CSV
        candidates = listLocalTflCsvs_(opts.DataDir);
        if isempty(candidates)
            error(['No local TfL CSVs found in %s and no URL/Filename supplied.'], opts.DataDir);
        end

        src.mode = "first_local";
        src.filename = string(candidates(1).name);
        src.localFile = fullfile(opts.DataDir, candidates(1).name);
        src.range_start = candidates(1).range_start;
        src.range_end   = candidates(1).range_end;
        return;
    end

    % ------------------------------------------------------------------
    % Case 3: TargetDate provided -> try local files first
    % ------------------------------------------------------------------
    targetDate = dateshift(opts.TargetDate, 'start', 'day');
    candidates = listLocalTflCsvs_(opts.DataDir);

    hits = [];
    for i = 1:numel(candidates)
        if ~isnat(candidates(i).range_start) && ~isnat(candidates(i).range_end)
            if targetDate >= candidates(i).range_start && targetDate <= candidates(i).range_end
                hits = [hits; i]; %#ok<AGROW>
            end
        end
    end

    if ~isempty(hits)
        widths = days([candidates(hits).range_end] - [candidates(hits).range_start]);
        [~, j] = min(widths);
        c = candidates(hits(j));

        src.mode = "local_by_date";
        src.filename = string(c.name);
        src.localFile = fullfile(opts.DataDir, c.name);
        src.range_start = c.range_start;
        src.range_end   = c.range_end;
        return;
    end

    % ------------------------------------------------------------------
    % Case 4: no local hit, but explicit URL supplied
    % ------------------------------------------------------------------
    if strlength(string(opts.URL)) > 0
        fname = filepartsFromUrl_(opts.URL);
        src.mode = "url_fallback";
        src.filename = string(fname);
        src.localFile = fullfile(opts.DataDir, fname);
        src.url = string(opts.URL);

        [d0, d1, ok] = parseTflDateRangeFromFilename_(src.filename);
        if ok
            src.range_start = d0;
            src.range_end   = d1;
        end
        return;
    end

    % ------------------------------------------------------------------
    % Case 5: no local hit and no explicit URL -> discover remote file
    % ------------------------------------------------------------------
    [url, filename, meta] = discoverRemoteTflJourneyFile_(targetDate);

    src.mode = "remote_discovery";
    src.filename = filename;
    src.localFile = fullfile(opts.DataDir, char(filename));
    src.url = url;
    src.range_start = meta.match_row.range_start;
    src.range_end   = meta.match_row.range_end;
    src.remote_meta = meta;
end
%{
function src = resolveTflSourceForDate_(opts)
% Resolve local TfL source file from:
%   1) explicit Filename, else
%   2) local DataDir scan using TargetDate, else
%   3) URL-derived filename if URL supplied

    src = struct();
    src.mode = "";
    src.filename = "";
    src.localFile = "";
    src.url = string(opts.URL);
    src.range_start = NaT;
    src.range_end   = NaT;

    % ------------------------------------------------------------------
    % Case 1: explicit filename override
    % ------------------------------------------------------------------
    if strlength(string(opts.Filename)) > 0
        src.mode = "explicit_filename";
        src.filename = string(opts.Filename);
        src.localFile = fullfile(opts.DataDir, char(src.filename));

        [d0, d1, ok] = parseTflDateRangeFromFilename_(src.filename);
        if ok
            src.range_start = d0;
            src.range_end   = d1;
        end
        return;
    end

    % ------------------------------------------------------------------
    % Case 2: if no TargetDate, use URL-derived filename if possible
    % ------------------------------------------------------------------
    if isempty(opts.TargetDate)
        if strlength(string(opts.URL)) > 0
            fname = filepartsFromUrl_(opts.URL);
            src.mode = "url_default";
            src.filename = string(fname);
            src.localFile = fullfile(opts.DataDir, fname);

            [d0, d1, ok] = parseTflDateRangeFromFilename_(src.filename);
            if ok
                src.range_start = d0;
                src.range_end   = d1;
            end
            return;
        end

        % Fallback: use the first matching TfL CSV found locally
        candidates = listLocalTflCsvs_(opts.DataDir);
        if isempty(candidates)
            error(['No local TfL CSVs found in %s and no URL/Filename supplied.'], opts.DataDir);
        end

        src.mode = "first_local";
        src.filename = string(candidates(1).name);
        src.localFile = fullfile(opts.DataDir, candidates(1).name);
        src.range_start = candidates(1).range_start;
        src.range_end   = candidates(1).range_end;
        return;
    end

    % ------------------------------------------------------------------
    % Case 3: TargetDate provided -> scan local files by filename date span
    % ------------------------------------------------------------------
    targetDate = dateshift(opts.TargetDate, 'start', 'day');
    candidates = listLocalTflCsvs_(opts.DataDir);

    hits = [];
    for i = 1:numel(candidates)
        if ~isnat(candidates(i).range_start) && ~isnat(candidates(i).range_end)
            if targetDate >= candidates(i).range_start && targetDate <= candidates(i).range_end
                hits = [hits; i]; %#ok<AGROW>
            end
        end
    end

    if ~isempty(hits)
        % Choose narrowest covering range
        widths = days(candidates(hits).range_end - candidates(hits).range_start);
        [~, j] = min(widths);
        c = candidates(hits(j));

        src.mode = "local_by_date";
        src.filename = string(c.name);
        src.localFile = fullfile(opts.DataDir, c.name);
        src.range_start = c.range_start;
        src.range_end   = c.range_end;
        return;
    end

    % ------------------------------------------------------------------
    % Case 4: no local hit, but URL supplied
    % ------------------------------------------------------------------
    if strlength(string(opts.URL)) > 0
        fname = filepartsFromUrl_(opts.URL);
        src.mode = "url_fallback";
        src.filename = string(fname);
        src.localFile = fullfile(opts.DataDir, fname);

        [d0, d1, ok] = parseTflDateRangeFromFilename_(src.filename);
        if ok
            src.range_start = d0;
            src.range_end   = d1;
        end
        return;
    end

    error(['No local TfL CSV in %s appears to cover TargetDate %s, and no URL/Filename was supplied.'], ...
        opts.DataDir, datestr(targetDate));
end
%}

% =========================================================================
function candidates = listLocalTflCsvs_(dataDir)
% List local CSVs and parse filename date ranges where possible.

    dd = dir(fullfile(dataDir, '*.csv'));
    candidates = struct('name', {}, 'range_start', {}, 'range_end', {});

    for i = 1:numel(dd)
        nm = string(dd(i).name);
        [d0, d1, ok] = parseTflDateRangeFromFilename_(nm);
        if ok
            candidates(end+1).name = dd(i).name; %#ok<AGROW>
            candidates(end).range_start = d0;
            candidates(end).range_end   = d1;
        end
    end
end

% =========================================================================
function fname = filepartsFromUrl_(url)
    [~, nm, ext] = fileparts(char(string(url)));
    fname = [nm ext];
end

% =========================================================================
function T = readTflCsv_(localFile)
    opts = detectImportOptions(localFile, ...
        'VariableNamingRule', 'preserve');

    try
        opts = setvaropts(opts, opts.VariableNames, 'WhitespaceRule', 'preserve');
        opts = setvaropts(opts, opts.VariableNames, 'EmptyFieldRule', 'auto');
    catch
    end

    T = readtable(localFile, opts);
end

% =========================================================================
function T = standardiseTflTrips_(Traw)
    vn = string(Traw.Properties.VariableNames);

    cNum       = findVar_(vn, "Number");
    cStartDate = findVar_(vn, "Start date");
    cStartID   = findVar_(vn, "Start station number");
    cStartName = findVar_(vn, "Start station");
    cEndDate   = findVar_(vn, "End date");
    cEndID     = findVar_(vn, "End station number");
    cEndName   = findVar_(vn, "End station");
    cBikeNum   = findVar_(vn, "Bike number");
    cBikeModel = findVar_(vn, "Bike model");
    cDur       = findVar_(vn, "Total duration");
    cDurMs     = findVar_(vn, "Total duration (ms)");

    T = table();

    T.journey_id          = toDoubleCol_(Traw.(cNum));
    T.start_date          = parseTfLDate_(Traw.(cStartDate));
    T.start_station_id    = toDoubleCol_(Traw.(cStartID));
    T.start_station_name  = strip(string(Traw.(cStartName)));

    T.end_date            = parseTfLDate_(Traw.(cEndDate));
    T.end_station_id      = toDoubleCol_(Traw.(cEndID));
    T.end_station_name    = strip(string(Traw.(cEndName)));

    T.bike_number         = strip(string(Traw.(cBikeNum)));
    T.bike_model          = strip(string(Traw.(cBikeModel)));

    T.total_duration_sec  = toDoubleCol_(Traw.(cDur));
    T.total_duration_ms   = toDoubleCol_(Traw.(cDurMs));

    missSec = ~isfinite(T.total_duration_sec) & isfinite(T.total_duration_ms);
    T.total_duration_sec(missSec) = T.total_duration_ms(missSec) ./ 1000;
end

% =========================================================================
function v = findVar_(vn, target)
    idx = find(strcmpi(vn, target), 1, 'first');
    if isempty(idx)
        error('Could not find required column: %s', target);
    end
    v = vn(idx);
end

% =========================================================================
function x = toDoubleCol_(v)
    if isnumeric(v)
        x = double(v);
        return;
    end
    s = strip(string(v));
    s(s == "") = missing;
    x = double(str2double(s));
end

% =========================================================================
function dt = parseTfLDate_(v)
    s = strip(string(v));
    s(s == "") = missing;

    fmts = { ...
        'dd/MM/yyyy HH:mm' ...
        'dd/MM/yyyy HH:mm:ss' ...
        'dd-MMM-yyyy HH:mm:ss' ...
        'yyyy-MM-dd HH:mm:ss' ...
        'yyyy-MM-dd''T''HH:mm:ss' ...
        };

    dt = NaT(size(s));
    for k = 1:numel(fmts)
        try
            tmp = datetime(s, 'InputFormat', fmts{k}, 'Format', 'yyyy-MM-dd HH:mm:ss');
            fill = isnat(dt) & ~isnat(tmp);
            dt(fill) = tmp(fill);
            if all(~isnat(dt) | ismissing(s))
                return;
            end
        catch
        end
    end

    try
        tmp = datetime(s, 'Format', 'yyyy-MM-dd HH:mm:ss');
        fill = isnat(dt) & ~isnat(tmp);
        dt(fill) = tmp(fill);
    catch
    end
end

% =========================================================================
function [stationTable, W, Wdur] = buildStationOd_(Tw)
    startIDs = Tw.start_station_id;
    endIDs   = Tw.end_station_id;

    allIDs = unique([startIDs; endIDs]);
    n = numel(allIDs);

    idxMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
    for k = 1:n
        idxMap(allIDs(k)) = k;
    end

    I = zeros(height(Tw),1);
    J = zeros(height(Tw),1);
    for r = 1:height(Tw)
        I(r) = idxMap(Tw.start_station_id(r));
        J(r) = idxMap(Tw.end_station_id(r));
    end

    W = sparse(I, J, 1, n, n);

    dur = Tw.total_duration_sec;
    dur(~isfinite(dur)) = 0;
    Wdur = sparse(I, J, dur, n, n);

    stationName = strings(n,1);
    startNameMap = firstNameMap_(Tw.start_station_id, Tw.start_station_name);
    endNameMap   = firstNameMap_(Tw.end_station_id,   Tw.end_station_name);

    for k = 1:n
        sid = allIDs(k);
        if isKey(startNameMap, sid)
            stationName(k) = startNameMap(sid);
        elseif isKey(endNameMap, sid)
            stationName(k) = endNameMap(sid);
        else
            stationName(k) = "";
        end
    end

    stationTable = table((1:n)', allIDs, stationName, ...
        full(sum(W,2)), full(sum(W,1))', ...
        'VariableNames', {'idx','station_id','station_name','out_trips','in_trips'});
end

% =========================================================================
function M = firstNameMap_(ids, names)
    M = containers.Map('KeyType', 'double', 'ValueType', 'char');
    for r = 1:numel(ids)
        sid = ids(r);
        nm  = strtrim(char(string(names(r))));
        if isfinite(sid) && ~isempty(nm) && ~isKey(M, sid)
            M(sid) = nm;
        end
    end
end

function [url, filename, meta] = discoverRemoteTflJourneyFile_(targetDate, varargin)
%DISCOVERREMOTETFLJOURNEYFILE_  Find the TfL Santander Cycles weekly CSV
% covering a requested date by querying the public S3 XML listing.
%
% Example:
%   [url, filename, meta] = discoverRemoteTflJourneyFile_(datetime(2024,6,2));

    p = inputParser;
    p.addParameter('ListingURL', ...
        "https://s3-eu-west-1.amazonaws.com/cycling.data.tfl.gov.uk/?list-type=2&max-keys=5000", ...
        @(x)ischar(x) || isstring(x));
    p.addParameter('PreferNarrowestMatch', true, @(x)islogical(x) && isscalar(x));
    p.addParameter('Verbose', false, @(x)islogical(x) && isscalar(x));
    p.parse(varargin{:});
    opt = p.Results;

    targetDate = datetime(targetDate);
    targetDate = dateshift(targetDate, 'start', 'day');

    listingURL = string(opt.ListingURL);
    xmlText = webread(char(listingURL));

    keys = extractJourneyKeysFromXml_(xmlText);
    if isempty(keys)
        error('Could not find any usage-stats JourneyDataExtract CSV keys in XML listing: %s', listingURL);
    end

    rows = struct('key', {}, 'filename', {}, 'range_start', {}, 'range_end', {}, ...
                  'url', {}, 'span_days', {});

    for i = 1:numel(keys)
        key = string(keys{i});  % e.g. usage-stats/148JourneyDataExtract06Feb2019-12Feb2019.csv
        [~, nm, ext] = fileparts(char(key));
        fn = string(nm) + string(ext);

        [d0, d1, ok] = parseTflDateRangeFromFilename_(fn);
        if ~ok
            continue;
        end

        rows(end+1).key = char(key); %#ok<AGROW>
        rows(end).filename    = char(fn);
        rows(end).range_start = d0;
        rows(end).range_end   = d1;
        rows(end).url         = char("https://cycling.data.tfl.gov.uk/" + key);
        rows(end).span_days   = days(d1 - d0);
    end

    if isempty(rows)
        error('Found candidate XML keys, but none could be parsed as JourneyDataExtract date ranges.');
    end

    T = struct2table(rows);
    hit = (targetDate >= T.range_start) & (targetDate <= T.range_end);

    if ~any(hit)
        [~, ord] = sort(abs(days(T.range_start - targetDate)));
        showN = min(5, height(T));
        nearby = T(ord(1:showN), {'filename','range_start','range_end'});
        error(['No remote TfL journey file appears to cover %s.\n' ...
               'Nearest parsed candidates:\n%s'], ...
            datestr(targetDate), evalc('disp(nearby)'));
    end

    Thit = T(hit, :);
    if opt.PreferNarrowestMatch && height(Thit) > 1
        [~, j] = min(Thit.span_days);
        pick = Thit(j, :);
    else
        pick = Thit(1, :);
    end

    url = string(pick.url{1});
    filename = string(pick.filename{1});

    meta = struct();
    meta.target_date  = targetDate;
    meta.listing_url  = listingURL;
    meta.n_candidates = height(T);
    meta.table        = T;
    meta.match_row    = pick;

    if opt.Verbose
        fprintf('Matched remote TfL file: %s\n', filename);
        fprintf('URL: %s\n', url);
        fprintf('Range: %s to %s\n', datestr(pick.range_start), datestr(pick.range_end));
    end
end

% =========================================================================
function keys = extractJourneyKeysFromXml_(xmlText)
% Extract <Key>...</Key> entries for usage-stats JourneyDataExtract CSVs.

    if isstring(xmlText), xmlText = char(xmlText); end

    toks = regexp(xmlText, '<Key>([^<]*usage-stats/[^<]*JourneyDataExtract[^<]*\.csv)</Key>', ...
        'tokens');

    keys = cell(size(toks));
    for i = 1:numel(toks)
        keys{i} = toks{i}{1};
    end

    % unique, stable
    if ~isempty(keys)
        keys = unique(string(keys), 'stable');
        keys = cellstr(keys);
    end
end

%{

function [url, filename, meta] = discoverRemoteTflJourneyFile_(targetDate, varargin)
%DISCOVERREMOTETFLJOURNEYFILE_  Find the TfL Santander Cycles weekly CSV
% covering a requested date by parsing the remote bucket listing.

    p = inputParser;
    p.addParameter('ListingURL', "https://cycling.data.tfl.gov.uk/", ...
        @(x)ischar(x) || isstring(x));
    p.addParameter('PreferNarrowestMatch', true, @(x)islogical(x) && isscalar(x));
    p.addParameter('Verbose', false, @(x)islogical(x) && isscalar(x));
    p.parse(varargin{:});
    opt = p.Results;

    targetDate = datetime(targetDate);
    targetDate = dateshift(targetDate, 'start', 'day');

    listingURL = string(opt.ListingURL);
    if ~endsWith(listingURL, "/")
        listingURL = listingURL + "/";
    end

    html = webread(char(listingURL));

    entries = extractJourneyEntries_(html);
    if isempty(entries)
        error('Could not find any JourneyDataExtract CSV links in remote listing: %s', listingURL);
    end

    rows = struct('filename', {}, 'href_rel', {}, 'range_start', {}, 'range_end', {}, ...
                  'url', {}, 'span_days', {});

    for i = 1:numel(entries)
        fn = string(entries(i).filename);
        hrefRel = string(entries(i).href_rel);

        [d0, d1, ok] = parseTflDateRangeFromFilename_(fn);
        if ~ok
            continue;
        end

        rows(end+1).filename = char(fn); %#ok<AGROW>
        rows(end).href_rel = char(hrefRel);
        rows(end).range_start = d0;
        rows(end).range_end   = d1;
        rows(end).url         = char(listingURL + hrefRel);
        rows(end).span_days   = days(d1 - d0);
    end

    if isempty(rows)
        error('Found candidate links, but none could be parsed as JourneyDataExtract date ranges.');
    end

    T = struct2table(rows);
    hit = (targetDate >= T.range_start) & (targetDate <= T.range_end);

    if ~any(hit)
        [~, ord] = sort(abs(days(T.range_start - targetDate)));
        showN = min(5, height(T));
        nearby = T(ord(1:showN), {'filename','range_start','range_end'});
        error(['No remote TfL journey file appears to cover %s.\n' ...
               'Nearest parsed candidates:\n%s'], ...
            datestr(targetDate), evalc('disp(nearby)'));
    end

    Thit = T(hit, :);
    if opt.PreferNarrowestMatch && height(Thit) > 1
        [~, j] = min(Thit.span_days);
        pick = Thit(j, :);
    else
        pick = Thit(1, :);
    end

    url = string(pick.url{1});
    filename = string(pick.filename{1});

    meta = struct();
    meta.target_date  = targetDate;
    meta.listing_url  = listingURL;
    meta.n_candidates = height(T);
    meta.table        = T;
    meta.match_row    = pick;

    if opt.Verbose
        fprintf('Matched remote TfL file: %s\n', filename);
        fprintf('URL: %s\n', url);
        fprintf('Range: %s to %s\n', datestr(pick.range_start), datestr(pick.range_end));
    end
end

%}

function entries = extractJourneyEntries_(html)
% Return struct array with fields:
%   .href_rel
%   .filename

    if isstring(html), html = char(html); end

    toks = regexp(html, 'href\s*=\s*"([^"]*JourneyDataExtract[^"]+\.csv)"', 'tokens');

    entries = struct('href_rel', {}, 'filename', {});
    seen = strings(0,1);

    for i = 1:numel(toks)
        hrefRel = string(toks{i}{1});

        % ignore absolute URLs to elsewhere
        if startsWith(hrefRel, "http://") || startsWith(hrefRel, "https://")
            continue;
        end

        [~, nm, ext] = fileparts(char(hrefRel));
        fn = string(nm) + string(ext);

        key = hrefRel + "||" + fn;
        if any(seen == key)
            continue;
        end
        seen(end+1,1) = key; %#ok<AGROW>

        entries(end+1).href_rel = char(hrefRel); %#ok<AGROW>
        entries(end).filename = char(fn);
    end
end

%{

function [d0, d1, ok] = parseTflDateRangeFromFilename_(fname)
    d0 = NaT; d1 = NaT; ok = false;

    fname = char(string(fname));
    tok = regexp(fname, ...
        'JourneyDataExtract(\d{2}[A-Za-z]{3}\d{2,4})-(\d{2}[A-Za-z]{3}\d{2,4})\.csv$', ...
        'tokens', 'once');

    if isempty(tok)
        return;
    end

    try
        d0 = parseCompactDateToken_(tok{1});
        d1 = parseCompactDateToken_(tok{2});
        ok = ~(isnat(d0) || isnat(d1));
    catch
        ok = false;
    end
end



function d = parseCompactDateToken_(s)
    s = char(string(s));
    if numel(s) == 9
        d = datetime(s, 'InputFormat', 'ddMMMyy');
    elseif numel(s) == 11
        d = datetime(s, 'InputFormat', 'ddMMMyyyy');
    else
        d = NaT;
    end
end

%}

%{

function [url, filename, meta] = discoverRemoteTflJourneyFile_(targetDate, varargin)
%DISCOVERREMOTETFLJOURNEYFILE_  Find the TfL Santander Cycles weekly CSV
% covering a requested date by parsing the remote usage-stats listing.
%
% [url, filename, meta] = discoverRemoteTflJourneyFile_(targetDate)
%
% Inputs:
%   targetDate : datetime or date-like value
%
% Name/value options:
%   'ListingURL' : default "https://cycling.data.tfl.gov.uk/usage-stats/"
%   'PreferNarrowestMatch' : true (default)
%   'Verbose' : false (default)
%
% Outputs:
%   url      : full remote URL to the matching CSV
%   filename : remote filename only
%   meta     : struct with fields
%                .target_date
%                .listing_url
%                .n_candidates
%                .table          (all parsed candidates as a table)
%                .match_row      (1-row table for selected file)
%
% Notes:
% - Parses filenames containing JourneyDataExtract<start>-<end>.csv
% - Supports both 2-digit and 4-digit years in the embedded date tokens,
%   e.g. 21Feb16 or 15Jun2025
% - Ignores the opaque prefix before JourneyDataExtract, e.g. 104 or 02b

    p = inputParser;
    p.addParameter('ListingURL', "https://cycling.data.tfl.gov.uk/#!usage-stats%2F", ...
        @(x)ischar(x) || isstring(x));
    p.addParameter('PreferNarrowestMatch', true, @(x)islogical(x) && isscalar(x));
    p.addParameter('Verbose', false, @(x)islogical(x) && isscalar(x));
    p.parse(varargin{:});
    opt = p.Results;

    targetDate = datetime(targetDate);
    targetDate = dateshift(targetDate, 'start', 'day');

    listingURL = string(opt.ListingURL);
    if ~endsWith(listingURL, "/")
        listingURL = listingURL + "/";
    end

    % ---------------------------------------------------------------------
    % 1) Read remote listing HTML
    % ---------------------------------------------------------------------
    html = webread(char(listingURL));

    % ---------------------------------------------------------------------
    % 2) Extract candidate filenames from hrefs and/or visible text
    % ---------------------------------------------------------------------
    filenames = extractJourneyFilenames_(html);

    if isempty(filenames)
        error(['Could not find any JourneyDataExtract CSV filenames in remote listing: %s'], ...
            listingURL);
    end

    filenames = unique(string(filenames), 'stable');

    % ---------------------------------------------------------------------
    % 3) Parse date ranges from filenames
    % ---------------------------------------------------------------------
    rows = struct('filename', {}, 'range_start', {}, 'range_end', {}, ...
                  'url', {}, 'span_days', {});

    for i = 1:numel(filenames)
        fn = filenames(i);
        [d0, d1, ok] = parseTflDateRangeFromFilename_(fn);
        if ~ok
            continue;
        end

        rows(end+1).filename = char(fn); %#ok<AGROW>
        rows(end).range_start = d0;
        rows(end).range_end   = d1;
        rows(end).url         = char(listingURL + fn);
        rows(end).span_days   = days(d1 - d0);
    end

    if isempty(rows)
        error(['Found candidate filenames in remote listing, but none could be parsed ' ...
               'as JourneyDataExtract date ranges.']);
    end

    T = struct2table(rows);

    % ---------------------------------------------------------------------
    % 4) Select row(s) covering targetDate
    % ---------------------------------------------------------------------
    hit = (targetDate >= T.range_start) & (targetDate <= T.range_end);

    if ~any(hit)
        % Helpful message with nearby spans
        [~, ord] = sort(abs(days(T.range_start - targetDate)));
        showN = min(5, height(T));
        nearby = T(ord(1:showN), {'filename','range_start','range_end'});
        error(['No remote TfL journey file in the listing appears to cover %s.\n' ...
               'Nearest parsed candidates:\n%s'], ...
            datestr(targetDate), evalc('disp(nearby)'));
    end

    Thit = T(hit, :);

    if opt.PreferNarrowestMatch && height(Thit) > 1
        [~, j] = min(Thit.span_days);
        pick = Thit(j, :);
    else
        pick = Thit(1, :);
    end

    url = string(pick.url{1});
    filename = string(pick.filename{1});

    meta = struct();
    meta.target_date  = targetDate;
    meta.listing_url  = listingURL;
    meta.n_candidates = height(T);
    meta.table        = T;
    meta.match_row    = pick;

    if opt.Verbose
        fprintf('Matched remote TfL file: %s\n', filename);
        fprintf('URL: %s\n', url);
        fprintf('Range: %s to %s\n', datestr(pick.range_start), datestr(pick.range_end));
    end
end

%}

% =========================================================================
function filenames = extractJourneyFilenames_(html)
% Extract candidate JourneyDataExtract CSV filenames from directory listing HTML.
%
% We look in:
%   1) href="...csv"
%   2) visible text containing JourneyDataExtract...csv
%
% Return as cellstr.

    if isstring(html), html = char(html); end

    % href targets
    hrefs = regexp(html, 'href\s*=\s*"([^"]+JourneyDataExtract[^"]+\.csv)"', ...
        'tokens');
    hrefNames = strings(0,1);
    for i = 1:numel(hrefs)
        href = string(hrefs{i}{1});
        [~, nm, ext] = fileparts(char(href));
        hrefNames(end+1,1) = nm + ext; %#ok<AGROW>
    end

    % visible text fallback
    txt = regexp(html, '([A-Za-z0-9_-]*JourneyDataExtract[0-9A-Za-z-]+\.csv)', ...
        'tokens');
    txtNames = strings(0,1);
    for i = 1:numel(txt)
        txtNames(end+1,1) = string(txt{i}{1}); %#ok<AGROW>
    end

    allNames = [hrefNames; txtNames];
    allNames = unique(allNames, 'stable');

    filenames = cellstr(allNames);
end

% =========================================================================
function [d0, d1, ok] = parseTflDateRangeFromFilename_(fname)
% Parse filenames like:
%   104JourneyDataExtract04Apr2018-10Apr2018.csv
%   02bJourneyDataExtract21Feb16-05Mar2016.csv

    d0 = NaT; d1 = NaT; ok = false;

    fname = char(string(fname));

    tok = regexp(fname, ...
        'JourneyDataExtract(\d{2}[A-Za-z]{3}\d{2,4})-(\d{2}[A-Za-z]{3}\d{2,4})\.csv$', ...
        'tokens', 'once');

    if isempty(tok)
        return;
    end

    try
        d0 = parseCompactDateToken_(tok{1});
        d1 = parseCompactDateToken_(tok{2});
        ok = ~(isnat(d0) || isnat(d1));
    catch
        ok = false;
    end
end


% =========================================================================
function d = parseCompactDateToken_(s)
% Support:
%   04Apr2018
%   21Feb16

    s = char(string(s));

    if numel(s) == 9
        % ddMMMyy
        d = datetime(s, 'InputFormat', 'ddMMMyy');
    elseif numel(s) == 11
        % ddMMMyyyy
        d = datetime(s, 'InputFormat', 'ddMMMyyyy');
    else
        d = NaT;
    end
end

function repoRoot = getRepoRoot_(startDir)
%GETREPOROOT_  Climb upward until the actual repository root is found.
%
% We prefer:
%   1) a .git directory, or
%   2) a directory containing both paper/ and toolbox/, or
%   3) a directory containing both paper/ and src/
%
% This avoids falsely stopping at nested directories such as paper/paper/.

    repoRoot = startDir;

    while true
        hasGit     = exist(fullfile(repoRoot, '.git'), 'dir') == 7;
        hasPaper   = exist(fullfile(repoRoot, 'paper'), 'dir') == 7;
        hasToolbox = exist(fullfile(repoRoot, 'toolbox'), 'dir') == 7;
        hasSrc     = exist(fullfile(repoRoot, 'src'), 'dir') == 7;

        if hasGit || (hasPaper && hasToolbox) || (hasPaper && hasSrc)
            return;
        end

        parent = fileparts(repoRoot);
        if strcmp(parent, repoRoot)
            error('Could not locate repo root starting from: %s', startDir);
        end
        repoRoot = parent;
    end
end