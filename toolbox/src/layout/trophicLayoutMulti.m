function [X, Y, compIdx, info] = trophicLayoutMulti(W, varargin)
%TROPHICLAYOUTMULTI  Full trophic layout pipeline with multi-component handling.
%
%   [X, Y, compIdx, info] = TROPHICLAYOUTMULTI(W, ...)
%
% INPUT
%   W          : adjacency matrix (n x n), directed, weighted or unweighted.
%
% WRAPPER OPTIONS (handled here):
%   'UseBarycentre'      : true (default) | false
%   'CoherenceThreshold' : scalar in [0,1], default 0.6
%                          Components with coherence C_comp(c) below this
%                          threshold are skipped by the barycentre sweep.
%
%   'BaryNumSweeps'      : integer, default 4
%   'BaryDirection'      : 'both' (default) | 'topdown' | 'bottomup'
%
%   'BaryAnchor'         : which levels to treat as fixed during barycentre:
%                          'bottom' (default) - anchor bottom level only
%                          'top'              - anchor top level only
%                          'both'             - try both and pick best
%
%   'UseXSmooth'         : false (default) | true
%   'XSmoothNumIter'     : integer, default 10
%   'XSmoothStepSize'    : scalar in (0,1], default 0.3
%
% TROPHIC-ANALYSIS OPTIONS (passed to trophic_levels):
%   'h0'                 : 'min' (default), 'wm', or 'sm'
%   'symtest'            : 1 (default) or 0
%   'targets'            : tau-matrix, same size as W
%
% LAYOUT OPTIONS (passed to trophicLayoutMulti_core and its LayoutFun):
%   Any other name–value pairs are forwarded, e.g. for TROPHICLAYOUT:
%   'NumIters', 'Seed', 'AttractStrength', etc.
%
% OUTPUT
%   X, Y     : layout coordinates for all nodes (n x 1).
%
%   compIdx  : component index for each node:
%              0 for isolates, 1..nComp for non-isolated weak components.
%
%   info     : struct combining:
%              - core layout info (from TROPHICLAYOUTMULTI_CORE)
%              - global and per-component incoherence/coherence:
%                   .F0_global, .C_global
%                   .F0_comp,   .C_comp
%              - metadata from trophic_levels (in info.levelsInfo)
%
% External dependencies (other files in this repository):
%   trophic_levels.m          – compute trophic levels for W
%   trophicLayout.m           – single-component trophic layout
%   trophicLayoutMulti_core.m – multi-component packing wrapper
%
% Internal subfunctions (defined below in this file):
%   tflBarycentreSweep  – Sugiyama-style barycentre reordering
%   tflXSmooth          – x-only smoothing / polishing of a layered layout
%   layoutScore         – scalar score used to compare candidate layouts
%   computeBarycentres  – helper for barycentre computation
%   stableOrder         – helper for stable sorting by barycentre and X
%
% Author:
%   Bazil Sansom (Warwick Business School, University of Warwick)
%   Contact: bazil.sansom@wbs.ac.uk
%
% -------------------------------------------------------------------------

    % ---------- basic checks ----------
    n = size(W,1);
    if size(W,2) ~= n
        error('trophicLayoutMulti:WNotSquare', ...
              'W must be a square adjacency matrix.');
    end
    if ~isnumeric(W)
        error('trophicLayoutMulti:WNotNumeric', ...
              'W must be numeric.');
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
    optsW.SnapDelta       = 1;              % layer spacing in y-units
    optsW.SnapMethod      = 'round';        % 'round' | 'floor' | 'ceil'
    optsW.SnapRefineIters = 200;            % x-only refinement iters (snapped only)
    optsW.MinSepFrac      = 0.06;           % snapped-only, fraction of component span (0 disables)
    % ---------- snapped-mode warnings / gating ----------
    optsW.SnapWarn       = true;          % warn if snapping on not-very-layered comps
    optsW.SnapWarnF0Max  = 0.1;           % threshold on incoherence F0 (lower = more layered)
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
        error('trophicLayoutMulti:BadArgs', ...
              'Optional arguments must be name/value pairs.');
    end

    k = 1;
    while k <= numel(varargin)
        name = varargin{k};
        val  = varargin{k+1};

        if ~ischar(name) && ~isstring(name)
            error('trophicLayoutMulti:BadParamName', ...
                  'Parameter names must be strings.');
        end
        lname = lower(char(name));

        switch lname
            % --- wrapper-specific options ---
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

            % --- Snapped rendering options ---
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
            % --- snapped-mode warnings / gating ---
            case 'snapwarn'
                 optsW.SnapWarn = logical(val);
            case 'snapwarnf0max'
                optsW.SnapWarnF0Max = val;
            case 'snapwarnscope'
                optsW.SnapWarnScope = lower(char(val));
            case 'snapwarnaction'
                optsW.SnapWarnAction = lower(char(val));
            
                % --- refinement force options (snapped only) ---
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


            % --- trophic_levels options ---
            case 'h0'
                argsLevels = [argsLevels, {name, val}];
            case 'symtest'
                argsLevels = [argsLevels, {name, val}];
            case 'targets'
                argsLevels = [argsLevels, {name, val}];

            % --- everything else goes to core/layout ---
            otherwise
                argsCore = [argsCore, {name, val}];
        end

        k = k + 2;
    end

    
    % ---------- validate wrapper string options ----------
    optsW.BaryDirection   = tfl.validateEnum(optsW.BaryDirection,   {'both','topdown','bottomup'}, 'BaryDirection');
    optsW.BaryAnchor      = tfl.validateEnum(optsW.BaryAnchor,      {'bottom','top','both'},      'BaryAnchor');

    optsW.HeightMode      = tfl.validateEnum(optsW.HeightMode,      {'continuous','snapped'},     'HeightMode');
    optsW.SnapMethod      = tfl.validateEnum(optsW.SnapMethod,      {'round','floor','ceil'},     'SnapMethod');

    optsW.SnapWarnScope   = tfl.validateEnum(optsW.SnapWarnScope,   {'component','global','both'}, 'SnapWarnScope');
    optsW.SnapWarnAction  = tfl.validateEnum(optsW.SnapWarnAction,  {'warn','none'},               'SnapWarnAction');
    
    
    % ---------- main pipeline ----------
    
    % ---------- 1. Trophic analysis: compute h (no coherence here) ----------
    % We deliberately set ComputeCoherence = false to avoid duplicate work.
    [h, infoLevels] = trophic_levels(W, argsLevels{:}, 'ComputeCoherence', false);
    h = h(:);

    % ---------- 2. Core multi-component layout using h ----------
    [X, Y, compIdx, infoCore] = trophicLayoutMulti_core(W, h, argsCore{:});

    nComp    = infoCore.nComp;
    h_global = infoCore.h_global(:);   % rescaled heights actually used

    % ---------- 3. Compute global and per-component F0, C ----------
    % We do this here so the components used for gating match the layout.
    [i_all, j_all, w_all] = find(W);

    if isempty(i_all)
        F0_global = 0;
    else
        % Here we assume a default target of tau_ij = 1 on each existing edge.
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
            F0_c      = num_c / den_c;
            F0_c      = max(0, min(1, F0_c));
            F0_comp(c) = round(F0_c, 10);
        else
            F0_comp(c) = 0;
        end
    end
    C_comp = 1 - F0_comp;

    % ---------- 4. Component-wise barycentre search with coherence gate ----------
    if optsW.UseBarycentre && nComp > 0
        C_thresh = optsW.CoherenceThreshold;

        % Decide which anchors to try based on BaryAnchor option
        switch optsW.BaryAnchor
            case 'bottom'
                anchorsToTry = {'bottom'};
            case 'top'
                anchorsToTry = {'top'};
            case 'both'
                anchorsToTry = {'bottom', 'top'};
            otherwise
                error('trophicLayoutMulti:BadBaryAnchor', ...
                      'BaryAnchor must be ''bottom'', ''top'', or ''both''.');
        end

        % Relative tolerance for "indistinguishable" scores (e.g. 5%)
        relTol = 0.05;

        for c = 1:nComp
            if C_comp(c) < C_thresh
                continue;   % skip less-layered / highly circulatory components
            end

            S = find(compIdx == c);
            if numel(S) <= 2
                continue;   % nothing interesting to reorder
            end

            Wc  = W(S,S);
            hc  = h_global(S);
            Xc0 = X(S);   % baseline from multi-component trophic layout

            candidates = {};
            scores     = [];

            % (a) baseline: no barycentre, no x-smoothing
            candidates{end+1} = Xc0;
            scores(end+1)     = layoutScore(Wc, hc, Xc0);

            % (b) barycentre candidates for different anchors
            for a = 1:numel(anchorsToTry)
                Xa = tflBarycentreSweep( ...
                        Wc, hc, Xc0, ...
                        'NumSweeps', optsW.BaryNumSweeps, ...
                        'Direction', optsW.BaryDirection, ...
                        'Anchor',    anchorsToTry{a});

                if optsW.UseXSmooth
                    Xa = tflXSmooth( ...
                            Wc, hc, Xa, ...
                            'NumIter',  optsW.XSmoothNumIter, ...
                            'StepSize', optsW.XSmoothStepSize);
                end

                candidates{end+1} = Xa;
                scores(end+1)     = layoutScore(Wc, hc, Xa);
            end

            % (c) choose best candidate with tolerance band favouring sweeps
            % candidates{1} = baseline, candidates{2..K} = barycentre sweeps
            Jbase     = scores(1);
            isSweep   = false(size(scores));
            isSweep(2:end) = true;
            JsweepMin = min(scores(isSweep));
            Jmin      = min(scores);

            if JsweepMin < Jbase * (1 - relTol)
                % A sweep is clearly better: just take the overall best
                [~, bestIdx] = min(scores);

            elseif Jbase < JsweepMin * (1 - relTol)
                % Baseline clearly better by more than tolerance
                bestIdx = 1;
            else
                % Within tolerance band: prefer best sweep
                sweepIdx = find(isSweep & scores == JsweepMin, 1, 'first');
                bestIdx  = sweepIdx;
            end

            X_best = candidates{bestIdx};

            % Write back into global X
            X(S) = X_best;
        end
    end

    % ---------- 4b. HeightMode rendering: continuous vs snapped ----------
    % Use the canonical heights actually used by the core layout.
    Y_cont = h_global(:); % Continuous plotting height (what core used)

    % Fit affine map: Y_cont ≈ a*h + b
    [a,b] = affineMap_(h, Y_cont);    % <-- returns scalars
    
    switch lower(optsW.HeightMode)
        case 'continuous'
            % Keep whatever Y core gave you (usually already ~ h_global),
            % but return the canonical continuous heights for safety.
            Y = Y_cont;

        case 'snapped'
            
        
            % --- warn if snapping is being used on not-very-layered components ---
            if optsW.SnapWarn
                switch lower(optsW.SnapWarnScope)
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

            % Snap heights for rendering
            Y_snap = snapHeights_(h, optsW.SnapDelta, optsW.SnapMethod);
            % Snap in raw-h units, then map to plot scale
            %h_snap_raw = snapHeights_(h, optsW.SnapDelta, optsW.SnapMethod); % usually SnapDelta=1
            %Y_snap     = a*h_snap_raw + b;


            % Refine X conditional on snapped Y (per component)
            if nComp > 0 && optsW.SnapRefineIters > 0
                for c = 1:nComp
                    S = find(compIdx == c);
                    if numel(S) <= 2, continue; end

                    Wc = W(S,S);
                    Xc = X(S);
                    Yc = Y_snap(S);

                    Xc = refineX_givenY_(Wc, Xc, Yc, optsW);
                    X(S) = Xc;
                end
            end

            % 2) smooth within snapped layers (prevents nasty clumps)
            %X = tflXSmooth(W, Y_snap, X, 'NumIter', 10, 'StepSize', 0.4);
            if nComp > 0
                for c = 1:nComp
                    S = find(compIdx == c);
                    if numel(S) <= 2, continue; end
                    X(S) = tflXSmooth(W(S,S), Y_snap(S), X(S), 'NumIter', 10, 'StepSize', 0.4);
                end
            end
            
            % (3) --- hard anti-overlap guarantee (snapped only) ---
            if optsW.MinSepFrac > 0 && nComp > 0
                for c = 1:nComp
                    S = find(compIdx == c);
                    if numel(S) <= 2, continue; end
                    
                    Xc = X(S);
                    Yc = Y_snap(S);
                        
                    span = max(Xc) - min(Xc);
                    minSep = optsW.MinSepFrac * max(span, 1);
                    
                    Xc = enforceMinLayerSpacing_(Xc, Yc, minSep);
                    X(S) = Xc;
                end
            end

            % 3) hard anti-overlap guarantee
            %span = max(X) - min(X);
            %minSep = 0.05 * max(span, 1);
            %X = enforceMinLayerSpacing_(X, Y_snap, minSep);

            % Final rendered Y is snapped
            Y = Y_snap;

        otherwise
            error('trophicLayoutMulti:BadHeightMode', ...
                  'HeightMode must be ''continuous'' or ''snapped''.');
    end

    
    % ---------- 5. Assemble info ----------
    info           = infoCore;      % start from core layout info
    info.h_global  = h_global;      % ensure we expose the actual heights used
    info.h_raw     = infoCore.h_raw;
    info.compIdx   = compIdx;
    info.F0_global = F0_global;
    info.C_global  = C_global;
    info.F0_comp   = F0_comp;
    info.C_comp    = C_comp;

    info.HeightMode = optsW.HeightMode;
    info.Y_cont     = Y_cont;
    info.Y_render   = Y;                 % equals Y_cont in continuous mode
    info.SnapDelta  = optsW.SnapDelta;
    info.SnapMethod = optsW.SnapMethod;


    % Expose levels-info (without coherence) for completeness
    info.levelsInfo = infoLevels;
end

% =====================================================================

function X = tflBarycentreSweep(W, h, X, varargin)
%TFLBARYCENTRESWEEP  Sugiyama-style barycentre sweeps on horizontal ordering.
%
%   X = TFLBARYCENTRESWEEP(W, h, X) post-processes a trophic layout by
%   reordering nodes horizontally within each layer to reduce edge crossings,
%   keeping vertical coordinates (h) fixed.
%
%   X = TFLBARYCENTRESWEEP(W, h, X, 'Name', Value, ...) supports:
%
%   'NumSweeps' : number of sweeps (default: 4).
%                 A sweep consists of a top-down pass and/or bottom-up pass.
%
%   'Direction' : 'both' (default), 'topdown', or 'bottomup'.
%
%   'Anchor'    : which layers to keep fixed:
%                 'both'   (default)  - freeze bottom & top layers
%                 'bottom'            - freeze bottom only
%                 'top'               - freeze top only
%                 'none'              - all layers movable
%
% INPUTS
%   W  : n x n adjacency matrix (directed). Can be weighted or unweighted.
%        Intended for well-layered graphs (edges mostly from layer k to k+1).
%   h  : n x 1 trophic levels / layer indices (numeric). Layers are defined
%        by discretising h with round(h).
%   X  : n x 1 horizontal coordinates from a baseline layout (e.g. TFL).
%
% OUTPUT
%   X  : n x 1 updated horizontal coordinates, same set of values per layer
%        as the input X, but permuted within each (movable) layer.
%
% NOTES
%   - This only permutes the mapping node -> X within each layer; the set of
%     X-positions on a given layer is preserved exactly.
%   - Designed initially for reasonably layered, single-component DAGs; the
%     function emits a warning and returns X unchanged if the graph does
%     not look layered.

    % ---------- Parse options ----------
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

    % ---------- Basic checks ----------
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

    % ---------- Build layers (discretised from h) ----------
    % For well-layered graphs, h should already be integers; round() is
    % safe and robust to tiny numerical noise.
    hDisc = round(h);  % integer "layer" labels

    [layerVals, ~, layerIdx] = unique(hDisc, 'sorted');  %#ok<ASGLU>
    nL = numel(layerVals);

    % Nodes in each layer
    layerNodes = cell(nL,1);
    for a = 1:nL
        layerNodes{a} = find(layerIdx == a);
    end

    % Optional sanity check for layering
    [ei, ej] = find(W);  % edges i -> j
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

    % ---------- Determine movable layer range from Anchor ----------
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
        return;  % nothing to do
    end

    % ---------- Main sweep loop ----------
    for s = 1:nSweeps %#ok<NASGU>

        % ----- Top-down pass -----
        if direction == "both" || direction == "topdown"
            % Need a previous layer, so a >= 2
            startTD = max(firstMovable, 2);
            endTD   = min(lastMovable, nL);
            for a = startTD:endTD
                prevLayer = a - 1;
                S         = layerNodes{a};      % nodes in this layer
                if numel(S) <= 1
                    continue;
                end

                Xi        = X(S);
                Xi_sorted = sort(Xi);           % preserve positions

                % Compute barycentres based on parents in previous layer
                b = computeBarycentres(W, X, S, layerIdx, prevLayer, 'parents');

                % Sort by barycentre with X as tiebreaker (stable-ish)
                order = stableOrder(b, Xi);

                S_new         = S(order);
                layerNodes{a} = S_new;
                X(S_new)      = Xi_sorted;
            end
        end

        % ----- Bottom-up pass -----
        if direction == "both" || direction == "bottomup"
            % Need a next layer, so a <= nL-1
            startBU = min(lastMovable, nL - 1);
            endBU   = max(firstMovable, 1);
            for a = startBU:-1:endBU
                nextLayer = a + 1;
                S         = layerNodes{a};
                if numel(S) <= 1
                    continue;
                end

                Xi        = X(S);
                Xi_sorted = sort(Xi);

                % Compute barycentres based on children in next layer
                b = computeBarycentres(W, X, S, layerIdx, nextLayer, 'children');

                order = stableOrder(b, Xi);

                S_new         = S(order);
                layerNodes{a} = S_new;
                X(S_new)      = Xi_sorted;
            end
        end

    end

end


% =====================================================================
function b = computeBarycentres(W, X, S, layerIdx, targetLayer, mode)
%COMPUTEBARYCENTRES  Barycentres for nodes in S, based on edges to a layer.
%
%   b = COMPUTEBARYCENTRES(W, X, S, layerIdx, targetLayer, mode)
%
%   mode = 'parents'  : use incoming edges from target layer j -> i.
%   mode = 'children' : use outgoing edges to target layer i -> j.
%
%   If a node has no neighbours in the target layer, its barycentre falls
%   back to its current X(i).

    nS = numel(S);
    b  = zeros(nS,1);

    % Boolean mask for target layer
    targetMask = (layerIdx == targetLayer);

    switch mode
        case 'parents'
            % Use incoming edges from target layer: j -> i
            for k = 1:nS
                i = S(k);
                % Parents j in target layer with W(j,i) ~= 0
                js = find(W(:,i) ~= 0 & targetMask);
                if isempty(js)
                    b(k) = X(i);
                else
                    wji = full(W(js,i));
                    b(k) = sum(wji .* X(js)) / sum(wji);  % weighted mean
                end
            end

        case 'children'
            % Use outgoing edges to target layer: i -> j
            for k = 1:nS
                i = S(k);
                % Children j in target layer with W(i,j) ~= 0
                js = find(W(i,:) ~= 0 & targetMask.');
                if isempty(js)
                    b(k) = X(i);
                else
                    wij = full(W(i,js));
                    b(k) = sum(wij .* X(js).') / sum(wij);  % weighted mean
                end
            end

        otherwise
            error('computeBarycentres:BadMode', ...
                  'mode must be ''parents'' or ''children''.');
    end
end


% =====================================================================
function order = stableOrder(b, Xi)
%STABLEORDER  Sort indices by barycentre b, with X as tiebreaker.
%
%   order = STABLEORDER(b, Xi) returns a permutation of indices such that
%   nodes are ordered by increasing barycentre; when barycentres are equal
%   or very close, their current X positions are used as a secondary key.
%   This tends to preserve relative order for near-ties.

    M = [b(:), Xi(:)];
    [~, order] = sortrows(M, [1 2]);  % sort by b, then by Xi
end


% =====================================================================

function X = tflXSmooth(W, h, X, varargin)
%TFLXSMOOTH  X-only smoothing / polishing step for a layered trophic layout.
%
%   X = TFLXSMOOTH(W, h, X) performs an x-only smoothing step on a layout,
%   keeping vertical coordinates (h) fixed. Nodes are smoothed towards the
%   mean x-position of their neighbours, while preserving within-layer
%   left–right ordering and the total horizontal span of each layer.
%
%   X = TFLXSMOOTH(W, h, X, 'Name', Value, ...) supports:
%
%   'NumIter'  : number of smoothing iterations (default: 10).
%   'StepSize' : smoothing step size alpha in (0,1], default: 0.3.
%                Update rule:
%                   X_new = (1 - alpha)*X + alpha*X_neighbour_mean
%
% INPUTS
%   W  : n x n adjacency matrix (directed or undirected, weighted or not).
%        Used only to define "neighbours" for smoothing.
%   h  : n x 1 trophic levels / layer indices (numeric).
%        Used to group nodes into layers (via round(h)).
%   X  : n x 1 horizontal coordinates (e.g. after a barycentre sweep).
%
% OUTPUT
%   X  : n x 1 polished horizontal coordinates.
%
% NOTES
%   - Vertical positions (h) are NOT changed by this function.
%   - Within each layer, the initial left–right ordering is preserved:
%     nodes may spread/compress but never overtake each other.
%   - The total horizontal span of each layer is kept equal to its
%     original span (from the input X), to avoid global clumping.
%   - A small internal "de-clumping" weight nudges positions towards an
%     evenly spaced grid within each layer, while preserving the original
%     span and order.

    % ---------- Parse options ----------
    p = inputParser;
    p.addParameter('NumIter', 10, @(z) isnumeric(z) && isscalar(z) && z >= 0);
    p.addParameter('StepSize', 0.3, @(z) isnumeric(z) && isscalar(z) && z > 0 && z <= 1);
    p.parse(varargin{:});
    opts = p.Results;

    nIter = round(opts.NumIter);
    alpha = opts.StepSize;

    % ---------- Basic checks ----------
    [nRows, nCols] = size(W);
    if nRows ~= nCols
        error('tflXSmooth:BadSize', 'W must be square.');
    end
    n = nRows;

    h = h(:);
    X = X(:);
    if numel(h) ~= n || numel(X) ~= n
        error('tflXSmooth:SizeMismatch', ...
              'Length of h and X must match size of W.');
    end

    if nIter == 0
        return;
    end

    if ~issparse(W)
        W = sparse(W);
    end

    % Store original positions for each layer to preserve spread
    X0 = X;

    % Symmetrise: neighbours treated as undirected for smoothing
    Wsym = W + W.';
    Wsym = spones(Wsym) .* max(Wsym, Wsym.');  % keep max weight if asymmetric

    % ---------- Layers and fixed within-layer order ----------
    hDisc = round(h);
    [layerVals, ~, layerIdx] = unique(hDisc, 'sorted'); %#ok<ASGLU>
    nL = numel(layerVals);

    layerNodes = cell(nL,1);
    layerOrder = cell(nL,1);  % fixed ordering to preserve
    for a = 1:nL
        S = find(layerIdx == a);
        layerNodes{a} = S;
        if numel(S) <= 1
            layerOrder{a} = S;
        else
            % Initial left–right order within this layer
            [~, ord] = sort(X(S));
            layerOrder{a} = S(ord);
        end
    end

    % Small de-clumping weight for within-layer spacing
    gamma = 0.2;  % gamma = 0 => no grid nudging; gamma = 1 => perfect grid

    % ---------- Iterative smoothing ----------
    for it = 1:nIter %#ok<NASGU>
        X_new = X;

        % Smooth each node towards the mean x of its neighbours
        for i = 1:n
            js = find(Wsym(i,:) ~= 0);
            if isempty(js)
                continue;  % isolated: leave as is
            end
            w_ij = full(Wsym(i,js));
            if all(w_ij == 0)
                x_bar = mean(X(js));
            else
                % Weighted mean of neighbour x-positions
                x_bar = sum(w_ij .* X(js).') / sum(w_ij);
            end
            X_new(i) = (1 - alpha)*X(i) + alpha*x_bar;
        end

        % Re-impose fixed ordering and original spread within each layer
        X = X_new;
        for a = 1:nL
            S = layerOrder{a};  % fixed order from initial layout
            if numel(S) <= 1
                continue;
            end

            % Current (smoothed) positions in this layer, in fixed order
            Xi = X(S);

            % Enforce monotone x in this fixed order by sorting Xi
            Xi_monotone = sort(Xi);

            % Rescale to match original layer span
            X0_layer = X0(S);
            min0 = min(X0_layer);
            max0 = max(X0_layer);
            min1 = min(Xi_monotone);
            max1 = max(Xi_monotone);

            if max1 > min1 && max0 > min0
                % Linear rescale from [min1,max1] -> [min0,max0]
                Xi_rescaled = (Xi_monotone - min1) / (max1 - min1);
                Xi_rescaled = Xi_rescaled * (max0 - min0) + min0;
            else
                % Degenerate case: keep as is
                Xi_rescaled = Xi_monotone;
            end

            % Gently de-clump within the layer towards an even grid
            m      = numel(S);
            x_grid = linspace(min0, max0, m).';   % perfectly even positions
            Xi_final = (1 - gamma)*Xi_rescaled + gamma*x_grid;

            X(S) = Xi_final;
        end
    end

end


% =====================================================================
function J = layoutScore(Wc, hc, Xc)
%LAYOUTSCORE  Scale-invariant scalar "goodness" measure for a 1D layout.
%
%   J = LAYOUTSCORE(Wc, hc, Xc) computes a simple, scale-invariant score
%   for a 1D horizontal layout on a single component. It is used to compare
%   candidate x-arrangements (e.g. before vs. after a barycentre sweep).
%
% INPUTS
%   Wc : adjacency/weight matrix for this component (n_c x n_c).
%   hc : trophic (or other) heights for this component (n_c x 1).
%        (Currently unused, but included for future alternative scores.)
%   Xc : x-coordinates for this component (n_c x 1).
%
% OUTPUT
%   J  : scalar score; lower is better.
%
% CURRENT DEFINITION:
%   Let width = max(Xc) - min(Xc).
%   If width > 0,
%       J = sum_{i->j} w_ij * ((x_j - x_i) / width)^2
%   else
%       J = 0.
%
%   This is invariant under global rescalings Xc -> a * Xc + b.
%   Intuitively, it penalises long edges in the horizontal direction.

    %#ok<INUSD>  % hc is currently unused but kept for extensibility

    Xc = Xc(:);   % ensure column
    [nRows, nCols] = size(Wc);
    if nRows ~= nCols
        error('layoutScore:BadSize', 'Wc must be square.');
    end
    if numel(Xc) ~= nRows
        error('layoutScore:SizeMismatch', ...
              'Length of Xc must match size of Wc.');
    end

    [i, j, w] = find(Wc);
    if isempty(i)
        J = 0;
        return;
    end

    width = max(Xc) - min(Xc);
    if width <= 0
        J = 0;
        return;
    end

    dx = (Xc(j) - Xc(i)) / width;
    J  = sum(w .* (dx.^2));
end


% =====================================================================
function Ysnap = snapHeights_(Y, delta, method)
%SNAPHEIGHTS_  Quantise heights to multiples of delta.
    if nargin < 2 || isempty(delta), delta = 1; end
    if nargin < 3 || isempty(method), method = 'round'; end
    if delta <= 0
        error('snapHeights_:BadDelta','SnapDelta must be > 0.');
    end

    z = Y(:) / delta;
    switch lower(method)
        case 'round'
            zq = round(z);
        case 'floor'
            zq = floor(z);
        case 'ceil'
            zq = ceil(z);
        otherwise
            error('snapHeights_:BadMethod','SnapMethod must be round/floor/ceil.');
    end
    Ysnap = delta * zq;
end


% =====================================================================
function X = refineX_givenY_(W, X, Y, optsW)
%REFINEX_GIVENY_  X-only force refinement with Y fixed.
%
% Uses:
%   - attraction on undirected support S = W+W'
%   - repulsion on all pairs using full (dx,dy) distance
% Updates only X (y fixed).

    if ~issparse(W), W = sparse(W); end
    X = X(:); Y = Y(:);
    n = size(W,1);
    if numel(X) ~= n || numel(Y) ~= n
        error('refineX_givenY_:SizeMismatch','X,Y must match W.');
    end

    T    = max(0, round(optsW.SnapRefineIters));
    eta  = optsW.SnapRefineEta;
    a    = optsW.SnapRefineDamp;
    ka   = optsW.SnapRefineKa;
    kr   = optsW.SnapRefineKr;
    vclip = optsW.SnapRefineVClip;

    if T == 0 || n <= 2
        return;
    end


    % Undirected spring support with weights
    %S = W + W.';
    %[si, sj, sw] = find(S);
    S = W + W.';
    S = triu(S,1);
    [si,sj,sw] = find(S);

    % Preserve original range to avoid component drift/scale explosion
    x0 = X;
    r0 = max(x0) - min(x0);
    m0 = mean(x0);

    vx = zeros(n,1);
    eps2 = 1e-12;  % numerical safeguard

    for t = 1:T
        Fx = zeros(n,1);

        % --- attraction: sum_j S_ij (xj-xi) ---
        % accumulate via edge list for speed
        %for e = 1:numel(sw)
         %   i = si(e); j = sj(e); w = sw(e);
         %   dx = X(j) - X(i);
        %    f  = ka * w * dx;
        %    Fx(i) = Fx(i) + f;
            % note: (j,i) also appears in (si,sj) because S is symmetric sparse,
            % so we don't double-add an explicit opposite force here.
        %end

        for e=1:numel(sw)
            i=si(e); j=sj(e); w=sw(e);
            dx = X(j)-X(i);
            f  = ka * w * dx;
            Fx(i) = Fx(i) + f;
            Fx(j) = Fx(j) - f;
        end

        % --- repulsion: all pairs (exact, O(n^2)) ---
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

        % --- velocity + position update (x-only) ---
        vx = a*vx + eta*Fx;

        if vclip > 0
            vx = max(min(vx, vclip), -vclip);
        end

        X = X + vx;
    end

    % Re-normalise to preserve original range and mean (if range is nonzero)
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

        % push rightwards to enforce spacing
        for t = 2:numel(idx)
            gap = X(idx(t)) - X(idx(t-1));
            if gap < minSep
                X(idx(t)) = X(idx(t-1)) + minSep;
            end
        end

        % recentre this layer (optional, keeps things tidy)
        X(idx) = X(idx) - mean(X(idx)) + mean(X(idx));
    end
end

% =====================================================================
function [a,b] = affineMap_(x, y)
%AFFINEMAP_  Scalar least-squares affine fit y ≈ a*x + b.
    x = x(:); y = y(:);
    xc = x - mean(x);
    yc = y - mean(y);
    den = sum(xc.^2);

    if den < eps
        a = 1;
    else
        a = sum(xc .* yc) / den;   % scalar
    end
    b = mean(y) - a*mean(x);
end
