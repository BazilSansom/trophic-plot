function out = fig_foodweb_FW24(opts)
%FIG_FOODWEB_FW24  Build the FW24 empirical food-web comparison figures.
%
% Uses load_foodweb_FW24() and passes improved short ecological-group labels
% through to fig_food_web.
%
% Example:
%   out = fig_foodweb_FW24();
%   out = fig_foodweb_FW24(struct('ShowLabels', false, 'Export', true));

    if nargin < 1 || isempty(opts), opts = struct(); end

    data = load_foodweb_FW24();

    % Default to improved short labels unless caller overrides
    if ~isfield(opts, 'Labels') || isempty(opts.Labels)
        opts.Labels = abbreviateFW24Labels_(data.nodes);
    end

    if ~isfield(opts, 'ExportBaseLayered') || isempty(opts.ExportBaseLayered)
        opts.ExportBaseLayered = 'fig_foodweb_FW24_tfl_vs_sugiyama';
    end

    if ~isfield(opts, 'ExportBaseSnapped') || isempty(opts.ExportBaseSnapped)
        opts.ExportBaseSnapped = 'fig_foodweb_FW24_tfl_vs_snapped';
    end

    out = fig_food_web(data.E, opts);
    out.data = data;
    out.labels_short = opts.Labels;
    out.labels_raw   = data.nodes;
end

function labelsOut = abbreviateFW24Labels_(labelsIn)
%ABBREVIATEFW24LABELS_  Map raw FW24 group labels to cleaner short labels.

    if isstring(labelsIn)
        labels = cellstr(labelsIn(:));
    elseif iscell(labelsIn)
        labels = labelsIn(:);
    else
        error('fig_foodweb_FW24:BadLabels', ...
            'labelsIn must be a cell array or string array.');
    end

    labelsOut = labels;

    for k = 1:numel(labels)
        raw = char(labels{k});
        key = regexprep(lower(raw), '[^a-z0-9]', '');

        if startsWith(key, 'sha')
            labelsOut{k} = 'Shark';

        elseif contains(key, 'pelampred') || contains(key, 'pelmedpred') || strcmp(key, 'pelmpred')
            labelsOut{k} = 'PelMPred';

        elseif contains(key, 'demmpred') || contains(key, 'demmedpred')
            labelsOut{k} = 'DemMPred';

        elseif contains(key, 'pelaspred') || contains(key, 'pelspred')
            labelsOut{k} = 'PelSPred';

        elseif contains(key, 'demspred')
            labelsOut{k} = 'DemSPred';

        elseif contains(key, 'demsomn') || contains(key, 'demsmomn') || contains(key, 'demsmallomn')
            labelsOut{k} = 'DemSOmn';

        elseif contains(key, 'pony')
            labelsOut{k} = 'PonyFish';

        elseif contains(key, 'squid')
            labelsOut{k} = 'Squid';

        elseif contains(key, 'pelapl') || contains(key, 'pelplankt') || contains(key, 'pelplank') ...
                || contains(key, 'pelpl') || contains(key, 'pelplanktiv')
            labelsOut{k} = 'PelPlankt';

        elseif contains(key, 'macroepi')
            labelsOut{k} = 'MacroEpi';

        elseif contains(key, 'benthinfa') || contains(key, 'benthinfau') || contains(key, 'benthinv')
            labelsOut{k} = 'BenthInfa';

        elseif contains(key, 'zoopl')
            labelsOut{k} = 'Zoopl';

        elseif contains(key, 'demherb')
            labelsOut{k} = 'DemHerb';

        elseif contains(key, 'phyto')
            labelsOut{k} = 'Phyto';

        elseif contains(key, 'macroph') || contains(key, 'macrophy')
            labelsOut{k} = 'Macrophy';

        elseif contains(key, 'detrit')
            labelsOut{k} = 'Detrit';

        else
            labelsOut{k} = raw;
        end
    end
end
