function [X, Y, compIdx, info] = trophicLayoutMulti(W, varargin)
%TROPHICLAYOUTMULTI  Full trophic layout pipeline with multi-component handling.
%
%   [X, Y, compIdx, info] = TROPHICLAYOUTMULTI(W, ...)
%
% Key change vs earlier versions:
%   After the core layout (which may use a globally rescaled height h_global),
%   we *always* convert the returned coordinates back into "tau-units"
%   (i.e., 1 unit on the Y-axis corresponds to one full trophic step tau=1).
%   This makes reference bands unambiguous in continuous mode.
%
% INPUT
%   W          : adjacency matrix (n x n), directed, weighted or unweighted.
%
% WRAPPER OPTIONS (handled here):
%   'UseBarycentre'      : true (default) | false
%   'CoherenceThreshold' : scalar in [0,1], default 0.6
%   'BaryNumSweeps'      : integer, default 4
%   'BaryDirection'      : 'both' (default) | 'topdown' | 'bottomup'
%   'BaryAnchor'         : 'bottom' (default) | 'top' | 'both'
%   'UseXSmooth'         : false (default) | true
%   'XSmoothNumIter'     : integer, default 10
%   'XSmoothStepSize'    : scalar in (0,1], default 0.3
%
% RENDERING OPTIONS:
%   'HeightMode'      : 'continuous' (default) | 'snapped'
%   'SnapDelta'       : default 1
%   'SnapMethod'      : 'round' (default) | 'floor' | 'ceil'
%   'SnapRefineIters' : default 200
%   'MinSepFrac'      : default 0.06
%
% Snapped-mode warnings:
%   'SnapWarn'       : true (default)
%   'SnapWarnF0Max'  : 0.1  (default)
%   'SnapWarnScope'  : 'component' (default) | 'global' | 'both'
%   'SnapWarnAction' : 'warn' (default) | 'none'
%
% Refinement force params (snapped only):
%   'SnapRefineEta','SnapRefineDamp','SnapRefineKa','SnapRefineKr','SnapRefineVClip'
%
% TROPHIC-ANALYSIS OPTIONS (passed to trophic_levels):
%   'h0', 'symtest', 'targets'
%
% LAYOUT OPTIONS (forwarded to trophicLayoutMulti_core / its LayoutFun):
%   Any other name–value pairs (NumIters, Seed, AttractStrength, etc.).
%
% OUTPUT
%   X, Y     : final coordinates (n x 1).
%              In 'continuous' mode, Y is in tau-units (so bands at integers
%              correspond to full trophic steps).
%   compIdx  : component index per node (0 isolates, 1..nComp components)
%   info     : diagnostics + mapping metadata
%
% -------------------------------------------------------------------------

    % ---------- basic checks ----------
    n = size(W,1);
    if size(W,2) ~= n
        error('trophicLayoutMulti:WNotSquare', 'W must be a square adjacency matrix.');
    end
    if ~isnumeric(W)
        error('trophicLayoutMulti:WNotNumeric', 'W must be numeric.');
    end

    % ---------- default wrapper options ----------
    optsW.UseBarycentre      = true;
    optsW.CoherenceThreshold = 0.6;
    optsW.BaryNumSweeps      = 4;
    optsW.BaryDirection      = 'both';
    optsW.BaryAnchor         = 'bottom';   % 'bottom' | 'top' | 'both'

    optsW.UseXSmooth         = false;
    optsW.XSmoothNumIter     = 10;
    optsW.XSmoothStepSize    = 0.3;

    % ---------- snapping / rendering options ----------
    optsW.HeightMode      = 'continuous';   % 'continuous' | 'snapped'
    optsW.SnapDelta       = 1;              % layer spacing in tau-units
    optsW.SnapMethod      = 'round';        % 'round' | 'floor' | 'ceil'
    optsW.SnapRefineIters = 200;            % x-only refinement iters (snapped only)
    optsW.MinSepFrac      = 0.06;           % snapped-only, fraction of component span (0 disables)

    % snapped-mode warnings / gating
    optsW.SnapWarn       = true;
    optsW.SnapWarnF0Max  = 0.1;
    optsW.SnapWarnScope  = 'component';   % 'component' | 'global' | 'both'
    optsW.SnapWarnAction = 'warn';        % 'warn' | 'none'

    % refinement force defaults (snapped only)
    optsW.SnapRefineEta   = 0.02;
    optsW.SnapRefineDamp  = 0.9;
    optsW.SnapRefineKa    = 1.0;
    optsW.SnapRefineKr    = 0.2;
    optsW.SnapRefineVClip = 0.5;

    % Buckets for where options go
    argsLevels = {};   % to trophic_levels
    argsCore   = {};   % to trophicLayoutMulti_core / LayoutFun

    % ---------- parse varargin into wrapper / levels / core ----------
    if mod(numel(varargin),2) ~= 0
        error('trophicLayoutMulti:BadArgs', 'Optional arguments must be name/value pairs.');
    end

    k = 1;
    while k <= numel(varargin)
        name = varargin{k};
        val  = varargin{k+1};

        if ~ischar(name) && ~isstring(name)
            error('trophicLayoutMulti:BadParamName', 'Parameter names must be strings.');
        end
        lname = lower(char(name));

        switch lname
            % wrapper-specific
            case 'usebarycentre'
                optsW.UseBarycentre = logical(val);
            case 'coherencethreshold'
                optsW.CoherenceThreshold = val;
            case 'barynumsweeps'
                optsW.BaryNumSweeps = val;
            case 'barydirection'
                optsW.BaryDirection = val;
            case 'baryanchor'
                optsW.BaryAnchor = lower(char(val));

            case 'usexsmooth'
                optsW.UseXSmooth = logical(val);
            case 'xsmoothnumiter'
                optsW.XSmoothNumIter = val;
            case 'xsmoothstepsize'
                optsW.XSmoothStepSize = val;

            % rendering
            case 'heightmode'
                optsW.HeightMode = lower(char(val));
            case 'snapdelta'
                optsW.SnapDelta = val;
            case 'snapmethod'
                optsW.SnapMethod = lower(char(val));
            case 'snaprefineiters'
                optsW.SnapRefineIters = val;
            case 'minsepfrac'
                optsW.MinSepFrac = val;

            % snapped warnings
            case 'snapwarn'
                optsW.SnapWarn = logical(val);
            case 'snapwarnf0max'
                optsW.SnapWarnF0Max = val;
            case 'snapwarnscope'
                optsW.SnapWarnScope = lower(char(val));
            case 'snapwarnaction'
                optsW.SnapWarnAction = lower(char(val));

            % snapped refinement forces
            case 'snaprefineeta'
                optsW.SnapRefineEta = val;
            case 'snaprefinedamp'
                optsW.SnapRefineDamp = val;
            case 'snaprefineka'
                optsW.SnapRefineKa = val;
            case 'snaprefinekr'
                optsW.SnapRefineKr = val;
            case 'snaprefinevclip'
                optsW.SnapRefineVClip = val;

            % trophic_levels options
            case 'h0'
                argsLevels = [argsLevels, {name, val}]; %#ok<AGROW>
            case 'symtest'
                argsLevels = [argsLevels, {name, val}]; %#ok<AGROW>
            case 'targets'
                argsLevels = [argsLevels, {name, val}]; %#ok<AGROW>

            % everything else to core/layout
            otherwise
                argsCore = [argsCore, {name, val}]; %#ok<AGROW>
        end

        k = k + 2;
    end

    % ---------- validate wrapper string options ----------
    optsW.BaryDirection  = tfl.validateEnum(optsW.BaryDirection, {'both','topdown','bottomup'}, 'BaryDirection');
    optsW.BaryAnchor     = tfl.validateEnum(optsW.BaryAnchor, {'bottom','top','both'}, 'BaryAnchor');

    optsW.HeightMode     = tfl.validateEnum(optsW.HeightMode, {'continuous','snapped'}, 'HeightMode');
    optsW.SnapMethod     = tfl.validateEnum(optsW.SnapMethod, {'round','floor','ceil'}, 'SnapMethod');

    optsW.SnapWarnScope  = tfl.validateEnum(optsW.SnapWarnScope, {'component','global','both'}, 'SnapWarnScope');
    optsW.SnapWarnAction = tfl.validateEnum(optsW.SnapWarnAction, {'warn','none'}, 'SnapWarnAction');

    % ===================== main pipeline =====================

    % ---------- 1) Trophic analysis: compute canonical h (tau-units) ----------
    [h, infoLevels] = trophic_levels(W, argsLevels{:}, 'ComputeCoherence', false);
    h = h(:);

    % ---------- 2) Core multi-component layout using h (core may internally rescale) ----------
    [X, Y, compIdx, infoCore] = trophicLayoutMulti_core(W, h, argsCore{:});
    nComp    = infoCore.nComp;
    h_global = infoCore.h_global(:); %#ok<NASGU> % exposed in info (diagnostic)

    % ---------- 3) Compute global and per-component F0, C in canonical h-units ----------
    [i_all, j_all, w_all] = find(W);

    if isempty(i_all)
        F0_global = 0;
    else
        dh_all    = h(j_all) - h(i_all) - 1;
        F0_global = sum(w_all .* (dh_all.^2)) / sum(w_all);
        F0_global = max(0, min(1, F0_global));
        F0_global = round(F0_global, 10);
    end
    C_global = 1 - F0_global;

    F0_comp = nan(nComp,1);
    for c = 1:nComp
        S = find(compIdx == c);
        if numel(S) <= 1
            F0_comp(c) = 0;
            continue;
        end

        inComp = ismember(i_all, S) & ismember(j_all, S);
        ic = i_all(inComp);
        jc = j_all(inComp);
        wc = w_all(inComp);

        if isempty(ic)
            F0_comp(c) = 0;
            continue;
        end

        dh_c   = h(jc) - h(ic) - 1;
        num_c  = sum(wc .* (dh_c.^2));
        den_c  = sum(wc);

        if den_c > 0
            F0_c       = num_c / den_c;
            F0_c       = max(0, min(1, F0_c));
            F0_comp(c) = round(F0_c, 10);
        else
            F0_comp(c) = 0;
        end
    end
    C_comp = 1 - F0_comp;

    % ---------- 4) Optional barycentre post-process (use canonical h for layering) ----------
    if optsW.UseBarycentre && nComp > 0
        C_thresh = optsW.CoherenceThreshold;

        switch optsW.BaryAnchor
            case 'bottom'
                anchorsToTry = {'bottom'};
            case 'top'
                anchorsToTry = {'top'};
            case 'both'
                anchorsToTry = {'bottom','top'};
            otherwise
                error('trophicLayoutMulti:BadBaryAnchor', 'BaryAnchor must be bottom/top/both.');
        end

        relTol = 0.05;

        for c = 1:nComp
            if C_comp(c) < C_thresh
                continue;
            end

            S = find(compIdx == c);
            if numel(S) <= 2
                continue;
            end

            Wc  = W(S,S);
            hc  = h(S);      % IMPORTANT: use canonical h, not h_global
            Xc0 = X(S);

            candidates = {};
            scores     = [];

            candidates{end+1} = Xc0;
            scores(end+1)     = layoutScore(Wc, hc, Xc0);

            for aa = 1:numel(anchorsToTry)
                Xa = tflBarycentreSweep( ...
                        Wc, hc, Xc0, ...
                        'NumSweeps', optsW.BaryNumSweeps, ...
                        'Direction', optsW.BaryDirection, ...
                        'Anchor',    anchorsToTry{aa});

                if optsW.UseXSmooth
                    Xa = tflXSmooth( ...
                            Wc, hc, Xa, ...
                            'NumIter',  optsW.XSmoothNumIter, ...
                            'StepSize', optsW.XSmoothStepSize);
                end

                candidates{end+1} = Xa; %#ok<AGROW>
                scores(end+1)     = layoutScore(Wc, hc, Xa); %#ok<AGROW>
            end

            Jbase     = scores(1);
            isSweep   = false(size(scores));
            isSweep(2:end) = true;
            JsweepMin = min(scores(isSweep));

            if JsweepMin < Jbase * (1 - relTol)
                [~, bestIdx] = min(scores);
            elseif Jbase < JsweepMin * (1 - relTol)
                bestIdx = 1;
            else
                sweepIdx = find(isSweep & scores == JsweepMin, 1, 'first');
                bestIdx  = sweepIdx;
            end

            X(S) = candidates{bestIdx};
        end
    end

    % ---------- 4b) Convert the core Y into tau-units (so 1 = one full trophic step) ----------
    Y_core = Y(:);
    H = [h, ones(numel(h),1)];
    p = H \ Y_core;                 % Y_core ≈ slope*h + intercept
    slopeY     = p(1);
    interceptY = p(2);

    if ~isfinite(slopeY) || abs(slopeY) < 1e-12
        warning('trophicLayoutMulti:BadYScale', ...
            'Degenerate Y~h fit (slope=%.3g). Skipping tau-unit rescale.', slopeY);
        X_tau = X(:);
        Y_tau = Y_core;
    else
        Y_tau = (Y_core - interceptY) / slopeY;
        x0    = mean(X);
        X_tau = x0 + (X - x0) / slopeY;    % isotropic rescale to preserve geometry
    end

    % ---------- 4c) HeightMode: continuous vs snapped (now in tau-units) ----------
    switch optsW.HeightMode
        case 'continuous'
            X = X_tau;
            Y = Y_tau;

        case 'snapped'
            % warn if snapping is used on poorly layered material
            if optsW.SnapWarn
                switch optsW.SnapWarnScope
                    case 'global'
                        bad = (F0_global > optsW.SnapWarnF0Max);
                    case 'both'
                        bad = (F0_global > optsW.SnapWarnF0Max) || any(F0_comp > optsW.SnapWarnF0Max);
                    otherwise % 'component'
                        bad = any(F0_comp > optsW.SnapWarnF0Max);
                end

                if bad && strcmpi(optsW.SnapWarnAction,'warn')
                    worstF0 = max(F0_comp);
                    nBad    = nnz(F0_comp > optsW.SnapWarnF0Max);
                    warning('trophicLayoutMulti:SnapNotRecommended', ...
                        ['HeightMode=''snapped'' requested, but %d/%d components have incoherence F0 > %.3g ' ...
                         '(worst F0 = %.3g). Snapped rendering may be visually misleading; ' ...
                         'consider HeightMode=''continuous'' for those cases.'], ...
                        nBad, max(nComp,1), optsW.SnapWarnF0Max, worstF0);
                end
            end

            X = X_tau;
            Y = snapHeights_(Y_tau, optsW.SnapDelta, optsW.SnapMethod);

            % Refine X conditional on snapped Y (per component)
            if nComp > 0 && optsW.SnapRefineIters > 0
                for c = 1:nComp
                    S = find(compIdx == c);
                    if numel(S) <= 2, continue; end
                    X(S) = refineX_givenY_(W(S,S), X(S), Y(S), optsW);
                end
            end

            % Smooth within snapped layers
            if nComp > 0
                for c = 1:nComp
                    S = find(compIdx == c);
                    if numel(S) <= 2, continue; end
                    X(S) = tflXSmooth(W(S,S), Y(S), X(S), 'NumIter', 10, 'StepSize', 0.4);
                end
            end

            % Hard min-separation within each snapped layer
            if optsW.MinSepFrac > 0 && nComp > 0
                for c = 1:nComp
                    S = find(compIdx == c);
                    if numel(S) <= 2, continue; end
                    Xc = X(S);
                    span = max(Xc) - min(Xc);
                    minSep = optsW.MinSepFrac * max(span, 1);
                    X(S) = enforceMinLayerSpacing_(X(S), Y(S), minSep);
                end
            end

        otherwise
            error('trophicLayoutMulti:BadHeightMode', ...
                  'HeightMode must be ''continuous'' or ''snapped''.');
    end

    % ---------- 5) Assemble info ----------
    info            = infoCore;
    info.h          = h;                % canonical trophic levels (tau-units)
    info.compIdx    = compIdx;

    info.F0_global  = F0_global;
    info.C_global   = C_global;
    info.F0_comp    = F0_comp;
    info.C_comp     = C_comp;

    info.HeightMode = optsW.HeightMode;
    info.SnapDelta  = optsW.SnapDelta;
    info.SnapMethod = optsW.SnapMethod;

    % mapping metadata
    info.Y_core            = Y_core;       % pre tau-unit conversion (from core)
    info.X_core            = infoCore;     %#ok<NASGU> % (kept in infoCore if needed)
    info.y_from_h_slope    = slopeY;
    info.y_from_h_intercept= interceptY;

    info.X_tau      = X_tau;
    info.Y_tau      = Y_tau;
    info.Y_render   = Y;

    info.levelsInfo = infoLevels;
end

% =====================================================================
function X = tflBarycentreSweep(W, h, X, varargin)
%TFLBARYCENTRESWEEP  Sugiyama-style barycentre sweeps on horizontal ordering.

    p = inputParser;
    p.addParameter('NumSweeps', 4, @(z) isnumeric(z) && isscalar(z) && z >= 0);
    p.addParameter('Direction', 'both', @(s) ischar(s) || isstring(s));
    p.addParameter('Anchor',    'both', @(s) ischar(s) || isstring(s));
    p.parse(varargin{:});
    opts = p.Results;

    direction = lower(string(opts.Direction));
    if ~ismember(direction, ["both","topdown","bottomup"])
        error('tflBarycentreSweep:BadDirection', ...
              'Direction must be ''both'', ''topdown'', or ''bottomup''.');
    end

    anchor = lower(string(opts.Anchor));
    if ~ismember(anchor, ["both","bottom","top","none"])
        error('tflBarycentreSweep:BadAnchor', ...
              'Anchor must be ''both'', ''bottom'', ''top'', or ''none''.');
    end

    [nRows, nCols] = size(W);
    if nRows ~= nCols
        error('tflBarycentreSweep:BadSize', 'W must be square.');
    end
    n = nRows;

    h = h(:);
    X = X(:);
    if numel(h) ~= n || numel(X) ~= n
        error('tflBarycentreSweep:SizeMismatch', ...
              'Length of h and X must match size of W.');
    end

    if ~issparse(W)
        W = sparse(W);
    end

    hDisc = round(h);
    [~, ~, layerIdx] = unique(hDisc, 'sorted');
    nL = max(layerIdx);

    layerNodes = cell(nL,1);
    for a = 1:nL
        layerNodes{a} = find(layerIdx == a);
    end

    [ei, ej] = find(W);
    if ~isempty(ei)
        dLayer      = layerIdx(ej) - layerIdx(ei);
        fracPerfect = nnz(dLayer == 1) / numel(dLayer);
        if fracPerfect < 0.5
            warning('tflBarycentreSweep:NotLayered', ...
                ['Graph does not look layered (%.1f%%%% of edges go ', ...
                 'to the next layer). Returning X unchanged.'], ...
                100*fracPerfect);
            return;
        end
    end

    switch anchor
        case "both"
            firstMovable = 2;
            lastMovable  = nL - 1;
        case "bottom"
            firstMovable = 2;
            lastMovable  = nL;
        case "top"
            firstMovable = 1;
            lastMovable  = nL - 1;
        case "none"
            firstMovable = 1;
            lastMovable  = nL;
    end

    nSweeps = round(opts.NumSweeps);
    if nSweeps == 0 || firstMovable > lastMovable || nL <= 1
        return;
    end

    for s = 1:nSweeps %#ok<NASGU>

        if direction == "both" || direction == "topdown"
            startTD = max(firstMovable, 2);
            endTD   = min(lastMovable, nL);
            for a = startTD:endTD
                prevLayer = a - 1;
                S = layerNodes{a};
                if numel(S) <= 1, continue; end

                Xi_sorted = sort(X(S));
                b = computeBarycentres(W, X, S, layerIdx, prevLayer, 'parents');
                order = stableOrder(b, X(S));

                S_new         = S(order);
                layerNodes{a} = S_new;
                X(S_new)      = Xi_sorted;
            end
        end

        if direction == "both" || direction == "bottomup"
            startBU = min(lastMovable, nL - 1);
            endBU   = max(firstMovable, 1);
            for a = startBU:-1:endBU
                nextLayer = a + 1;
                S = layerNodes{a};
                if numel(S) <= 1, continue; end

                Xi_sorted = sort(X(S));
                b = computeBarycentres(W, X, S, layerIdx, nextLayer, 'children');
                order = stableOrder(b, X(S));

                S_new         = S(order);
                layerNodes{a} = S_new;
                X(S_new)      = Xi_sorted;
            end
        end
    end
end

% =====================================================================
function b = computeBarycentres(W, X, S, layerIdx, targetLayer, mode)
    nS = numel(S);
    b  = zeros(nS,1);
    targetMask = (layerIdx == targetLayer);

    switch mode
        case 'parents'
            for k = 1:nS
                i = S(k);
                js = find(W(:,i) ~= 0 & targetMask);
                if isempty(js)
                    b(k) = X(i);
                else
                    wji = full(W(js,i));
                    b(k) = sum(wji .* X(js)) / sum(wji);
                end
            end
        case 'children'
            for k = 1:nS
                i = S(k);
                js = find(W(i,:) ~= 0 & targetMask.');
                if isempty(js)
                    b(k) = X(i);
                else
                    wij = full(W(i,js));
                    b(k) = sum(wij .* X(js).') / sum(wij);
                end
            end
        otherwise
            error('computeBarycentres:BadMode', 'mode must be ''parents'' or ''children''.');
    end
end

% =====================================================================
function order = stableOrder(b, Xi)
    M = [b(:), Xi(:)];
    [~, order] = sortrows(M, [1 2]);
end

% =====================================================================
function X = tflXSmooth(W, h, X, varargin)
    p = inputParser;
    p.addParameter('NumIter', 10, @(z) isnumeric(z) && isscalar(z) && z >= 0);
    p.addParameter('StepSize', 0.3, @(z) isnumeric(z) && isscalar(z) && z > 0 && z <= 1);
    p.parse(varargin{:});
    opts = p.Results;

    nIter = round(opts.NumIter);
    alpha = opts.StepSize;

    [nRows, nCols] = size(W);
    if nRows ~= nCols, error('tflXSmooth:BadSize', 'W must be square.'); end
    n = nRows;

    h = h(:);
    X = X(:);
    if numel(h) ~= n || numel(X) ~= n
        error('tflXSmooth:SizeMismatch', 'Length of h and X must match size of W.');
    end
    if nIter == 0, return; end
    if ~issparse(W), W = sparse(W); end

    X0 = X;

    Wsym = W + W.';
    Wsym = spones(Wsym) .* max(Wsym, Wsym.');

    hDisc = round(h);
    [~, ~, layerIdx] = unique(hDisc, 'sorted');
    nL = max(layerIdx);

    layerOrder = cell(nL,1);
    for a = 1:nL
        S = find(layerIdx == a);
        if numel(S) <= 1
            layerOrder{a} = S;
        else
            [~, ord] = sort(X(S));
            layerOrder{a} = S(ord);
        end
    end

    gamma = 0.2;

    for it = 1:nIter %#ok<NASGU>
        X_new = X;

        for i = 1:n
            js = find(Wsym(i,:) ~= 0);
            if isempty(js), continue; end
            w_ij = full(Wsym(i,js));
            if all(w_ij == 0)
                x_bar = mean(X(js));
            else
                x_bar = sum(w_ij .* X(js).') / sum(w_ij);
            end
            X_new(i) = (1 - alpha)*X(i) + alpha*x_bar;
        end

        X = X_new;

        for a = 1:nL
            S = layerOrder{a};
            if numel(S) <= 1, continue; end

            Xi_monotone = sort(X(S));

            X0_layer = X0(S);
            min0 = min(X0_layer); max0 = max(X0_layer);
            min1 = min(Xi_monotone); max1 = max(Xi_monotone);

            if max1 > min1 && max0 > min0
                Xi_rescaled = (Xi_monotone - min1) / (max1 - min1);
                Xi_rescaled = Xi_rescaled * (max0 - min0) + min0;
            else
                Xi_rescaled = Xi_monotone;
            end

            m      = numel(S);
            x_grid = linspace(min0, max0, m).';
            X(S)   = (1 - gamma)*Xi_rescaled + gamma*x_grid;
        end
    end
end

% =====================================================================
function J = layoutScore(Wc, hc, Xc) %#ok<INUSD>
    Xc = Xc(:);
    [nRows, nCols] = size(Wc);
    if nRows ~= nCols, error('layoutScore:BadSize', 'Wc must be square.'); end
    if numel(Xc) ~= nRows, error('layoutScore:SizeMismatch', 'Xc length must match Wc.'); end

    [i, j, w] = find(Wc);
    if isempty(i), J = 0; return; end

    width = max(Xc) - min(Xc);
    if width <= 0, J = 0; return; end

    dx = (Xc(j) - Xc(i)) / width;
    J  = sum(w .* (dx.^2));
end

% =====================================================================
function Ysnap = snapHeights_(Y, delta, method)
    if nargin < 2 || isempty(delta), delta = 1; end
    if nargin < 3 || isempty(method), method = 'round'; end
    if delta <= 0
        error('snapHeights_:BadDelta','SnapDelta must be > 0.');
    end

    z = Y(:) / delta;
    switch lower(method)
        case 'round', zq = round(z);
        case 'floor', zq = floor(z);
        case 'ceil',  zq = ceil(z);
        otherwise
            error('snapHeights_:BadMethod','SnapMethod must be round/floor/ceil.');
    end
    Ysnap = delta * zq;
end

% =====================================================================
function X = refineX_givenY_(W, X, Y, optsW)
    if ~issparse(W), W = sparse(W); end
    X = X(:); Y = Y(:);
    n = size(W,1);
    if numel(X) ~= n || numel(Y) ~= n
        error('refineX_givenY_:SizeMismatch','X,Y must match W.');
    end

    T     = max(0, round(optsW.SnapRefineIters));
    eta   = optsW.SnapRefineEta;
    damp  = optsW.SnapRefineDamp;
    ka    = optsW.SnapRefineKa;
    kr    = optsW.SnapRefineKr;
    vclip = optsW.SnapRefineVClip;

    if T == 0 || n <= 2
        return;
    end

    S = W + W.';
    S = triu(S,1);
    [si,sj,sw] = find(S);

    x0 = X;
    r0 = max(x0) - min(x0);
    m0 = mean(x0);

    vx = zeros(n,1);
    eps2 = 1e-12;

    for t = 1:T %#ok<NASGU>
        Fx = zeros(n,1);

        for e = 1:numel(sw)
            i = si(e); j = sj(e); w = sw(e);
            dx = X(j) - X(i);
            f  = ka * w * dx;
            Fx(i) = Fx(i) + f;
            Fx(j) = Fx(j) - f;
        end

        for i = 1:n-1
            xi = X(i); yi = Y(i);
            for j = i+1:n
                dx = xi - X(j);
                dy = yi - Y(j);
                r2 = dx*dx + dy*dy + eps2;
                invr3 = 1 / (r2 * sqrt(r2));
                f = kr * dx * invr3;
                Fx(i) = Fx(i) + f;
                Fx(j) = Fx(j) - f;
            end
        end

        vx = damp*vx + eta*Fx;
        if vclip > 0
            vx = max(min(vx, vclip), -vclip);
        end
        X = X + vx;
    end

    r1 = max(X) - min(X);
    if r0 > 0 && r1 > 0
        X = (X - mean(X)) * (r0 / r1) + m0;
    else
        X = X - mean(X) + m0;
    end
end

% =====================================================================
function X = enforceMinLayerSpacing_(X, Y, minSep)
    X = X(:); Y = Y(:);
    levs = unique(Y);
    for k = 1:numel(levs)
        idx = find(Y == levs(k));
        if numel(idx) <= 1, continue; end

        [~, ord] = sort(X(idx));
        idx = idx(ord);

        for t = 2:numel(idx)
            gap = X(idx(t)) - X(idx(t-1));
            if gap < minSep
                X(idx(t)) = X(idx(t-1)) + minSep;
            end
        end
    end
end