% ======================================================================
function flat = flatten_mangal(allDs)
% Flatten common allDs shapes into a network-level struct array with .W.
flat = struct([]);

if isempty(allDs)
    return;
end

% If cell array: concatenate flatten of each cell
if iscell(allDs)
    tmp = cell(size(allDs));
    for i = 1:numel(allDs)
        tmp{i} = flatten_mangal_(allDs{i});
    end
    flat = vertcat(tmp{:});
    return;
end

if ~isstruct(allDs)
    return;
end

% Case A: already network-level
if isfield(allDs, 'W')
    flat = allDs;
    return;
end

% Case B: dataset-level, with .ds inside
if isfield(allDs, 'ds')
    parts = {};
    for d = 1:numel(allDs)
        if isempty(allDs(d).ds), continue; end
        ds = allDs(d).ds;
        % propagate dataset metadata if present
        for k = 1:numel(ds)
            if ~isfield(ds,'dataset_id') && isfield(allDs(d),'dataset_id')
                ds(k).dataset_id = allDs(d).dataset_id;
            end
            if ~isfield(ds,'dataset_name') && isfield(allDs(d),'name')
                ds(k).dataset_name = allDs(d).name;
            end
            if ~isfield(ds,'dataset_name') && isfield(allDs(d),'dataset_name')
                ds(k).dataset_name = allDs(d).dataset_name;
            end
            if ~isfield(ds,'network_id') && isfield(ds,'network_id')
                % nothing
            end
        end
        parts{end+1,1} = ds; %#ok<AGROW>
    end
    if ~isempty(parts)
        flat = vertcat(parts{:});
    else
        flat = struct([]);
    end
    return;
end

% Fallback: try common field 'networks'
if isfield(allDs,'networks') && isstruct(allDs.networks) && isfield(allDs.networks,'W')
    flat = allDs.networks;
    return;
end

end