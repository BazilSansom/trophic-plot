function data = mangal_api_get(endpoint, params)
%MANGAL_API_GET  GET JSON from Mangal API v2.
%
% data = mangal_api_get('network', struct('public',true,'count',1000,'page',0));
%
% Docs: https://mangal.io/api/v2/ ... :contentReference[oaicite:1]{index=1}

    if nargin < 2 || isempty(params), params = struct(); end
    base = "https://mangal.io/api/v2/";
    url  = base + string(endpoint);

    % MATLAB webread likes name/value pairs; convert struct -> pairs
    nv = {};
    f = fieldnames(params);
    for k = 1:numel(f)
        v = params.(f{k});
        if islogical(v), v = double(v); end
        nv(end+1:end+2) = {f{k}, v}; %#ok<AGROW>
    end

    opts = weboptions('ContentType','json','Timeout',60);
    data = webread(url, nv{:}, opts);
end